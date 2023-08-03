module spmm_mod
    use omp_lib
    use psb_base_mod
    use sp3mm_config_mod
    use sp3mm_acc_mod
    use idx_tree_mod
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

        end subroutine spmm_size_uperbound_col_parts
        
        subroutine spmm_upper_bound_serial(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
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

        subroutine spmm_upper_bound_row_by_row(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
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

            !$omp parallel do schedule(runtime)
            do row = 1, a_m
                call sp_accs(row)%init(rows_sizes(row))
                do col = a%irp(row), a%irp(row + 1) - 1
                    call upper_bound_scalar_sparse_mul_row(sp_accs(row), a%val(col), b, a%ja(col))
                end do
                call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%ja, info)
                call psb_realloc(sp_accs(row)%nnz, sp_accs(row)%nnz_idxs, info)
                call psb_msort(sp_accs(row)%ja, sp_accs(row)%nnz_idxs)
            end do
            !$omp end parallel do
            call merge_rows(sp_accs, c)


        end subroutine spmm_upper_bound_row_by_row

        subroutine spmm_rb_tree_serial(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
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

        subroutine spmm_rb_tree_row_by_row(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n
            integer(psb_ipk_)   :: row, col
            type(idx_tree), allocatable :: row_accs(:)

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            allocate(row_accs(a_m))
            call c%allocate(a_m, b_n)

            !$omp parallel do schedule(runtime)
            do row = 1, a_m
                row_accs(row)%nnz = 0
                nullify(row_accs(row)%root)
                do col = a%irp(row), a%irp(row + 1) - 1
                    call rb_tree_scalar_sparse_row_mul(row_accs(row), a%val(col), b, a%ja(col))
                end do 
            end do
            !$omp end parallel do
            call merge_trees_distrib(row_accs, c)

            deallocate(row_accs)
            info = 0
        end subroutine spmm_rb_tree_row_by_row

        subroutine spmm_rb_tree_row_by_row_1D_blocks(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info
            
            integer(psb_ipk_)   :: a_m, b_n

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

            c%irp(1) = 1

        end subroutine spmm_rb_tree_row_by_row_1D_blocks

        subroutine spmm_rb_tree_row_by_row_2D_blocks(a,b,c,cfg,info)
            implicit none
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate(a_m, b_n)

        end subroutine spmm_rb_tree_row_by_row_2D_blocks

end module spmm_mod