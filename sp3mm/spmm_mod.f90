module spmm_mod
    use omp_lib
    use psb_base_mod
    use sp3mm_config_mod
    use sp3mm_acc_mod
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
            ! TODO fix inconsitencies in omp do
            !omp parallel do schedule(static) reduction(+:full_mat_bound)
            do r = 1, a_m
                do jj = a%irp(r), a%irp(r + 1) - 1
                    j = a%ja(jj)
                    rlen = b%irp(j+1) - b%irp(j)
                    rows_sizes(r) = rows_sizes(r) + rlen
                    full_mat_bound = full_mat_bound + rlen
                end do
            end do
            !omp end parallel do

            rows_sizes(a_m + 1) = full_mat_bound

            end = omp_get_wtime()
            write (*, '("spMMSizeUpperbound:", I0, A, ES12.6, " s")')&
            full_mat_bound, achar(9), end - start
        end subroutine spmm_size_uperbound

        subroutine spmm_size_uperbound_col_parts(a, b, cols, &
                            b_col_offsets, rows_parts_sizes, info)
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            integer(psb_ipk_), intent(in) :: cols
            integer(psb_ipk_), allocatable, intent(in) :: b_col_offsets(:)
            integer(psb_ipk_), allocatable, intent(out):: rows_parts_sizes(:)
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_) :: full_mat_bound
            integer(psb_ipk_) :: row, jj, j, rlen, col_group, b_part_id

            call psb_realloc(a%get_nrows() * cols + 1, rows_parts_sizes, info)

            rows_parts_sizes = 0

            !omp parallel do schedule(static) reduction(+:full_mat_bound)
            do row = 0, a%get_nrows() - 1
                do jj = a%irp(row), a%irp(row + 1) - 1
                    j = a%ja(jj)
                    b_part_id = col_group + j * cols + 1
                    do col_group = 0, cols - 1
                        rlen = b_col_offsets(b_part_id + 1) - b_col_offsets(b_part_id)
                        rows_parts_sizes(col_group + row * cols + 1) = rlen + &
                                    rows_parts_sizes(col_group + row * cols + 1)
                        full_mat_bound = full_mat_bound + rlen
                    end do
                end do
            end do
            !omp end parallel do

            rows_parts_sizes(a%get_nrows() * cols + 1) = full_mat_bound

        end subroutine spmm_size_uperbound_col_parts
        
        subroutine spmm_serial(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n
            type(dense_acc)    :: acc
            integer(psb_ipk_)   :: row, col

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate_mnnz(a_m, b_n)

            c%irp(1) = 1

            call acc%init(b_n)

            do row=1, a_m
                do col = a%irp(row), a%irp(row + 1) - 1
                    call scalar_sparse_row_mul(acc, a%val(col), b, a%ja(col))
                end do
                call sparsify_direct(acc, c, row)
                call acc%reset()
            end do

        info = 0
        end subroutine spmm_serial

        subroutine spmm_row_by_row(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n
            type(dense_acc), allocatable :: accs(:)
            type(dense_acc) :: acc
            type(sparse_acc), allocatable:: sp_accs(:)
            integer(psb_ipk_)   :: i, row, col
            integer(psb_ipk_), allocatable :: rows_sizes(:)
            real(8) :: tic, toc
            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            c%irp(1) = 1

            ! call spmm_size_uperbound(a, b, rows_sizes, info)

            ! allocating the dense accumulators
            allocate(accs(cfg%thread_num))

            do i = 1, cfg%thread_num
                call accs(i)%init(b_n)
            end do

            ! allocating the sparse accumulators (rows)
            allocate(sp_accs(a_m))

            call cfg%chunk_distrib_func(a_m, info)

            tic = omp_get_wtime()
            !$omp parallel do schedule(runtime), private(acc)
            do row=1, a_m
                acc = accs(omp_get_thread_num() + 1)
                do col = a%irp(row), a%irp(row + 1) - 1
                    call scalar_sparse_row_mul(acc, a%val(col), b, a%ja(col))
                end do
                call sparsify_ub(acc, sp_accs(row), 0)
                call acc%reset()
            end do
            !$omp end parallel do
            
            toc = omp_get_wtime()
            print *, "calculated rows in ", toc - tic
            tic = omp_get_wtime()
            call merge_rows(sp_accs, c)
            toc = omp_get_wtime()
            print *, "merge rows in ", toc - tic
        end subroutine spmm_row_by_row

        subroutine spmm_row_by_row_1D_blocks(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info
            
            integer(psb_ipk_)   :: a_m, b_n
            type(dense_acc), allocatable :: accs(:)
            type(dense_acc) :: acc
            type(sparse_acc), allocatable:: sp_accs(:)
            integer(psb_ipk_)   :: i, row, col
            integer(psb_ipk_) :: row_block, row_block_rem, block_num, blck, start_row

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            c%irp(1) = 1

            ! call spmm_size_uperbound(a, b, rows_sizes, info)

            ! allocating the dense accumulators
            allocate(accs(cfg%thread_num))

            do i = 1, cfg%thread_num
                call accs(i)%init(b_n)
            end do

            allocate(sp_accs(a_m))

            row_block = a_m / cfg%rows
            row_block_rem = mod(a_m, cfg%rows)

            call cfg%chunk_distrib_func(cfg%rows, info)

            !$omp parallel do schedule(runtime) private(acc, start_row, blck)
            do block_num = 0, cfg%rows - 1
                if (block_num < row_block_rem) then
                    blck = row_block + 1
                else
                    blck = row_block
                end if

                start_row = block_num * row_block + min(block_num, row_block_rem)
                
                acc = accs(omp_get_thread_num() + 1) 
                
                do row = start_row + 1, start_row + blck
                    do col = a%irp(row), a%irp(row + 1) - 1
                        call scalar_sparse_row_mul(acc, a%val(col), b, a%ja(col))
                    end do
                    call sparsify_ub(acc, sp_accs(row), 0)
                    call acc%reset()
                end do
            end do
            !$omp end parallel do

            call merge_rows(sp_accs, c)
        end subroutine spmm_row_by_row_1D_blocks

        subroutine spmm_row_by_row_2D_blocks(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n
            type(dense_acc), allocatable :: accs(:)
            type(dense_acc) :: acc
            type(sparse_acc), allocatable:: sp_accs(:)
            integer(psb_ipk_) :: i
            integer(psb_ipk_), allocatable :: b_col_offsets(:)
            integer(psb_ipk_) :: grid_size, a_sub_row_num, rows_parts_sizes_num
            integer(psb_ipk_), allocatable :: rows_parts_sizes(:)
            integer(psb_ipk_) :: row_block, row_block_rem
            integer(psb_ipk_) :: col_block, col_block_rem
            ! integer(psb_ipk_) :: row_block, col_block
            integer(psb_ipk_) :: row_start, col_start


            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            grid_size = cfg%rows * cfg%cols
            a_sub_row_num = a_m * cfg%cols

            row_block = a_m / cfg%rows
            row_block_rem = mod(a_m, cfg%rows)

            col_block = b_n / cfg%cols
            col_block_rem = mod(b_n, cfg%cols)

            call cols_offset_partitioning_unif_ranges(b, cfg%cols, b_col_offsets)

            rows_parts_sizes_num = a_sub_row_num

            call spmm_size_uperbound_col_parts(a, b, cfg%cols, b_col_offsets, &
                                                rows_parts_sizes, info)

            allocate(accs(grid_size))

            do i = 1, grid_size
                if (col_block_rem == 0) then
                    call accs(i)%init(col_block)
                else
                    call accs(i)%init(col_block + 1)
                end if
            end do

            allocate(sp_accs(a_m))

        end subroutine spmm_row_by_row_2D_blocks

        subroutine cols_offset_partitioning_unif_ranges(mat, grid_cols, offsets)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: mat
            integer(psb_ipk_), intent(in) :: grid_cols
            integer(psb_ipk_), allocatable, intent(out):: offsets(:)

            integer(psb_ipk_) :: sub_rows_num, col_block, col_block_rem
            integer(psb_ipk_) :: info
            integer(psb_ipk_) :: row, j, col_group, col_group_col_start

            sub_rows_num = mat%get_nrows() * grid_cols
            col_block = mat%get_ncols() / grid_cols
            col_block_rem = mod(mat%get_ncols(), grid_cols)

            call psb_realloc(sub_rows_num + 1, offsets, info)

            j = 1
            do row = 0, mat%get_nrows() - 1
                j = mat%irp(row + 1)
                offsets(1 + row * grid_cols) = j
                do col_group = 1, grid_cols - 1
                    col_group_col_start = col_group * col_block + min(col_group, col_block_rem)
                    do while (j < mat%irp(row + 2) .and. mat%ja(j) < col_group_col_start)
                        j = j + 1
                    end do
                    offsets(col_group + row * grid_cols + 1) = j
                end do
            end do

            offsets(sub_rows_num + 1) = mat%get_nzeros() + 1
        end subroutine
end module spmm_mod