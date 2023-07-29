module spmm_mod
    implicit none

    contains
        subroutine spmm_size_uperbound(a,b,rows_sizes, info)
            use psb_base_mod
            use omp_lib
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


        
        subroutine spmm_serial(a,b,c,cfg,info)
            use psb_base_mod
            use sp3mm_config_mod
            use sp3mm_acc_mod
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
        end subroutine spmm_serial

        subroutine spmm_row_by_row(a,b,c,cfg,info)
            use psb_base_mod
            use sp3mm_config_mod
            use sp3mm_acc_mod
            type(psb_d_csr_sparse_mat), intent(in) :: a,b
            type(psb_d_csr_sparse_mat), intent(out):: c
            type(sp3mm_config), intent(in) :: cfg
            integer(psb_ipk_), intent(out) :: info

            integer(psb_ipk_)   :: a_m, b_n
            type(dense_acc), allocatable :: accs(:)
            integer(psb_ipk_)   :: i, row, col

            a_m = a%get_nrows()
            b_n = b%get_ncols()

            write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
            achar(9), achar(9), achar(9), a_m, b_n

            call c%allocate_mnnz(a_m, b_n)

            c%irp(1) = 1

            ! allocating the dense accumulators
            allocate(accs(cfg%thread_num))

            do i = 1, cfg%thread_num
                call accs(i)%init(b_n)
            end do

            

        end subroutine spmm_row_by_row

end module spmm_mod