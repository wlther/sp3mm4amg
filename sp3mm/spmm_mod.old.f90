module spmm_mod
    use omp_lib
    use psb_base_mod
    use sp3mm_const_mod
    use sp3mm_config_mod
    use sp3mm_acc_mod
    use idx_tree_mod
    ! use simple_tree_mod
    implicit none

    contains
        subroutine spmm_size_uperbound(a,b,rows_sizes, info)
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            integer, allocatable, intent(out) :: rows_sizes(:)
            integer, intent(out) :: info

            integer :: full_mat_bound
            integer :: a_m
            integer :: r, jj, j, rlen
            real(8) :: start, end

            start = omp_get_wtime()
            a_m = a%get_nrows()
            call psb_realloc(a_m + 1, rows_sizes, info)

            rows_sizes = 0
            full_mat_bound = 0
            do r = 1, a_m
                do jj = a%irp(r), a%irp(r + 1) - 1
                    j = a%ja(jj)
                    rlen = b%irp(j+1) - b%irp(j)
                    rows_sizes(r) = rows_sizes(r) + rlen
                    full_mat_bound = full_mat_bound + rlen
                end do
            end do

            rows_sizes(a_m + 1) = full_mat_bound

            end = omp_get_wtime()
            write (*, '("spMMSizeUpperbound:", I0, A, ES12.6, " s")')&
            full_mat_bound, achar(9), end - start
        end subroutine spmm_size_uperbound

        ! subroutine spmm_size_uperbound_col_parts(a, b, cols, &
        !                     b_col_offsets, rows_parts_sizes, info)
        !     type(psb_d_csr_sparse_mat), intent(in) :: a,b
        !     integer(psb_ipk_), intent(in) :: cols
        !     integer(psb_ipk_), allocatable, intent(in) :: b_col_offsets(:)
        !     integer(psb_ipk_), allocatable, intent(out):: rows_parts_sizes(:)
        !     integer(psb_ipk_), intent(out) :: info

        ! end subroutine spmm_size_uperbound_col_parts
        
        subroutine spmm_upper_bound_serial(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info
            
            integer(psb_ipk_)   :: a_m, b_n
            integer(psb_ipk_), allocatable :: rows_sizes(:)
            integer(psb_ipk_)   :: row, col
            type(sparse_acc), allocatable :: sp_accs(:)

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)
            
            call spmm_size_uperbound(a, b, rows_sizes, info)
            
            allocate(sp_accs(a_m))

            do row = 1, a_m
                call sp_accs(row)%init(rows_sizes(row))
                do col = a%irp(row), a%irp(row + 1) - 1
                    call upper_bound_scalar_sparse_mul_row(sp_accs(row), a%val(col), b, a%ja(col))
                end do
                call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%ja, info)
                call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%nnz_idxs, info)
                call psb_msort(sp_accs(row)%ja, sp_accs(row)%nnz_idxs)
            end do
            call merge_rows(sp_accs, c)

        end subroutine spmm_upper_bound_serial

        subroutine spmm_upper_bound_row_by_row(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info
            
            integer(psb_ipk_)   :: a_m, b_n
            integer(psb_ipk_), allocatable :: rows_sizes(:)
            integer(psb_ipk_)   :: row, col
            type(sparse_acc), allocatable :: sp_accs(:)
            real(8) :: tic, toc
            ! , avg_tic, avg_toc, avg

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            tic = omp_get_wtime()
            call c%allocate(a_m, b_n)
            toc = omp_get_wtime()

            call spmm_size_uperbound(a, b, rows_sizes, info)
            
            allocate(sp_accs(a_m))

            tic = omp_get_wtime()
            !$omp parallel do schedule(runtime)
            do row = 1, a_m
                call sp_accs(row)%init(rows_sizes(row))
                do col = a%irp(row), a%irp(row + 1) - 1
                    call upper_bound_scalar_sparse_mul_row(sp_accs(row), a%val(col), b, a%ja(col))
                end do
                call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%ja, info)
                call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%nnz_idxs, info)
                call psb_qsort(sp_accs(row)%ja, sp_accs(row)%nnz_idxs)
            end do
            !$omp end parallel do
            toc = omp_get_wtime()

            write (*,*)
            write (*,'("- computing rows: ", ES12.3)') toc - tic
            tic = omp_get_wtime()
            call merge_rows(sp_accs, c)
            toc = omp_get_wtime()
            write (*,*)
            write (*,'("- merging rows: ", ES12.3)') toc - tic

        end subroutine spmm_upper_bound_row_by_row

        subroutine spmm_upper_bound_row_by_row_1D_blocks(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info
            
            integer(psb_ipk_)   :: a_m, b_n, row_block, row_block_rem
            integer(psb_ipk_)   :: grid_rows, block_size, row_start
            integer(psb_ipk_), allocatable :: rows_sizes(:)
            integer(psb_ipk_)   :: block_num, row, col
            type(sparse_acc), allocatable :: sp_accs(:)
            real(8) :: tic, toc

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            tic = omp_get_wtime()
            call c%allocate(a_m, b_n)
            toc = omp_get_wtime()

            call spmm_size_uperbound(a, b, rows_sizes, info)
            allocate(sp_accs(a_m))

            grid_rows = omp_get_max_threads() * spmm_1D_block_iteration_factor
            row_block = a_m / grid_rows
            row_block_rem = mod(a_m, grid_rows)
            tic = omp_get_wtime()
            !$omp parallel do schedule(runtime) private(block_size, row_start)
            do block_num = 0, grid_rows - 1
                if (block_num < row_block_rem) then
                    block_size = row_block + 1
                else
                    block_size = row_block
                end if
                row_start = (block_num) * row_block + min(block_num, row_block_rem) + 1
                do row = row_start, row_start + block_size - 1
                    call sp_accs(row)%init(rows_sizes(row))
                    do col = a%irp(row), a%irp(row + 1) - 1
                        call upper_bound_scalar_sparse_mul_row(sp_accs(row), &
                                                    a%val(col), b, a%ja(col))
                    end do
                    call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%ja, info)
                    call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%nnz_idxs, info)
                    call psb_msort(sp_accs(row)%ja, sp_accs(row)%nnz_idxs)
                end do
            end do
            !$omp end parallel do
            toc = omp_get_wtime()
            write (*,*)
            write (*,'("- computing rows: ", ES12.3)') toc - tic
            tic = omp_get_wtime()
            call merge_rows(sp_accs, c)
            toc = omp_get_wtime()
            write (*,*)
            write (*,'("- merging rows: ", ES12.3)') toc - tic

        end subroutine spmm_upper_bound_row_by_row_1D_blocks

        subroutine spmm_rb_tree_serial(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n
            integer(psb_ipk_)   :: row, col
            type(idx_tree), allocatable :: row_accs(:)

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            allocate(row_accs(a_m))
            call c%allocate(a_m, b_n)

            do row = 1, a_m
                row_accs(row)%nnz = 0
                nullify(row_accs(row)%root)
                do col = a%irp(row), a%irp(row + 1) - 1
                    call rb_tree_scalar_sparse_row_mul(row_accs(row), a%val(col), b, a%ja(col))
                end do 
            end do
            call merge_trees(row_accs, c)

            deallocate(row_accs)

            info = 0
        end subroutine spmm_rb_tree_serial

        subroutine spmm_rb_tree_row_by_row(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n
            integer(psb_ipk_)   :: row, col
            type(idx_tree), allocatable :: row_accs(:)
            real(8) :: tic, toc

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n
            call c%allocate(a_m, b_n)

            allocate(row_accs(a_m))
            call c%allocate(a_m, b_n)

            tic = omp_get_wtime()
            !$omp parallel do schedule(runtime)
            do row = 1, a_m
                row_accs(row)%nnz = 0
                nullify(row_accs(row)%root)
                do col = a%irp(row), a%irp(row + 1) - 1
                    call rb_tree_scalar_sparse_row_mul(row_accs(row), a%val(col), b, a%ja(col))
                end do 
            end do
            !$omp end parallel do
            toc = omp_get_wtime()

            write (*,*)
            write (*,'("- computing rows: ", ES12.3)') toc - tic

            tic = omp_get_wtime()
            call merge_trees_distrib(row_accs, c)
            toc = omp_get_wtime()

            write (*,*)
            write (*,'("- merging rows: ", ES12.3)') toc - tic

            deallocate(row_accs)
            info = 0
        end subroutine spmm_rb_tree_row_by_row

        subroutine spmm_rb_tree_row_by_row_1D_blocks(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info
            
            integer(psb_ipk_)   :: a_m, b_n

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            c%irp(1) = 1

            info = 0
        end subroutine spmm_rb_tree_row_by_row_1D_blocks

        subroutine spmm_rb_tree_row_by_row_2D_blocks(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            info = 0
        end subroutine spmm_rb_tree_row_by_row_2D_blocks

        subroutine compute_indices(a, b, c, info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info

            integer :: full_mat_bound
            integer :: row, col, i, j, k, nnz

            full_mat_bound = 0
            !omp parallel do schedule(static) reduction(+:full_mat_bound)
            do row = 1, a%get_nrows()
                do col = a%irp(row), a%irp(row + 1) - 1
                    j = a%ja(col)
                    full_mat_bound = full_mat_bound + b%irp(j+1) - b%irp(j)
                end do
            end do
            !omp end parallel do

            call psb_realloc(a%get_nrows() + 1, c%irp, info)
            call psb_realloc(full_mat_bound, c%ja, info)
            c%ja = 0
            c%irp(1) = 1

            nnz = 0
            
            do row = 1, a%get_nrows()
                do col = a%irp(row), a%irp(row + 1) - 1
                    do i = b%irp(a%ja(col)), b%irp(a%ja(col) + 1) - 1
                        k = 0
                        do while(c%ja(c%irp(row) + k) /= 0 .and. c%ja(c%irp(row) + k) /= b%ja(i))
                            k = k + 1
                        end do
                        if (c%ja(c%irp(row) + k) == 0) then
                            c%ja(c%irp(row)+k) = b%ja(i)
                            nnz = nnz + 1
                        end if
                    end do
                end do
                c%irp(row + 1) = nnz + 1
                call psb_qsort(c%ja(c%irp(row):c%irp(row + 1)-1))
            end do


            call psb_realloc(nnz, c%ja, info)
            call psb_realloc(nnz, c%val, info)

            c%val = 0
        end subroutine compute_indices

        subroutine direct_scalar_sparse_row_mul(out_mat, out_row_num, scalar, mat, trgt_row_num)
            type(psb_d_csr_sparse_mat), intent(inout) :: out_mat
            integer(psb_ipk_), intent(in) :: out_row_num
            real(psb_dpk_), intent(in) :: scalar
            type(psb_d_csr_sparse_mat), intent(in) :: mat
            integer(psb_ipk_), intent(in) :: trgt_row_num

            integer(psb_ipk_) :: i, k, row_start, row_end

            row_start = out_mat%irp(out_row_num)
            row_end = out_mat%irp(out_row_num + 1) - 1

            do i = mat%irp(trgt_row_num), mat%irp(trgt_row_num + 1) - 1
                do k = out_mat%irp(out_row_num), out_mat%irp(out_row_num + 1) - 1
                    if (out_mat%ja(k) == mat%ja(i)) then
                        out_mat%val(k) = out_mat%val(k) + scalar * mat%val(i)
                        exit
                    end if
                end do
            end do

        end subroutine direct_scalar_sparse_row_mul

        subroutine spmm_idx_map_row_by_row(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n, row, col
            real(8) :: tic, toc

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)
            
            tic = omp_get_wtime()
            call compute_indices(a, b, c, info)
            toc = omp_get_wtime()
            write (*,*)
            write (*,'("- computing indices: ", ES12.3)') toc - tic

            tic = omp_get_wtime()
            !$omp parallel do schedule(runtime)
            do row = 1, a_m
                do col = a%irp(row), a%irp(row + 1) - 1
                    call direct_scalar_sparse_row_mul(c, row, a%val(col), b, a%ja(col))
                end do 
            end do
            !$omp end parallel do
            toc = omp_get_wtime()
            write (*,*)
            write (*,'("- computing rows: ", ES12.3)') toc - tic

        end subroutine spmm_idx_map_row_by_row

        subroutine spmm_omp(a,b,c,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_) :: a_m, b_n, row
            integer(psb_ipk_) :: i,j,k, nnz
            integer(psb_ipk_), allocatable :: rows_sizes(:), sorted_ja(:), ix(:)
            real(8) :: tic, toc

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            call psb_realloc(a_m + 1, c%irp, info)

            call spmm_size_uperbound(a, b, rows_sizes, info)
            
            call psb_realloc(rows_sizes(a_m + 1), c%ja, info)
            call psb_realloc(rows_sizes(a_m + 1), c%val, info)
            
            c%ja(1) = 0
            c%val(1) = 0

            ! call psb_realloc(b%get_nzeros(), sorted_ja, info)
            call psb_realloc(b%get_nzeros(), ix, info)
            sorted_ja = b%ja
            call psb_msort(sorted_ja, ix)

            i = 1
            nnz = 1
            do row = 1, a_m
                c%irp(i) = nnz
                do j = 1, b_n
                    k = a%irp(row)
                    do while(sorted_ja(i) == j .and. k < a%irp(row + 1))
                        ! if the row of b is lesser than the column of a
                        ! go to the next row of b
                        if (ix(i) < b%irp(a%ja(k))) then
                            i = i + 1
                        ! if the column of a is lesser than the row of b
                        ! go to the next column of a
                        else if (ix(i) >= b%irp(a%ja(k) + 1)) then
                            k = k + 1
                        ! right row and right column
                        else
                            c%ja(nnz) = j
                            c%val(nnz) = c%val(nnz) + a%val(k) * b%val(ix(i))
                            i = i + 1
                            k = k + 1
                        end if
                    end do
                    if (c%ja(nnz) /= 0) then
                        nnz = nnz + 1
                        c%ja(nnz) = 0
                        c%val(nnz) = 0
                    end if
                end do
            end do

            c%irp(a_m + 1) = nnz
            call psb_realloc(c%get_nzeros(), c%ja, info)
            call psb_realloc(c%get_nzeros(), c%val, info)
        end subroutine spmm_omp
end module spmm_mod