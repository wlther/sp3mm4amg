module spmm_mod
    use omp_lib
    use psb_base_mod
    implicit none
    
contains
    ! gustavson's algorithm using perfect hashing 
    ! and OpenMP parallelisation
    subroutine spmm_omp(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info

        integer(psb_ipk_)   :: ma, nb
        real(psb_dpk_), allocatable :: vals(:), acc(:)
        integer(psb_ipk_), allocatable :: col_inds(:)
        integer(psb_ipk_)   :: irw, jj, j, k, nnz, rwnz, thread_upperbound, start_idx
        real(8) :: tic, toc

        ma = a%get_nrows()
        nb = b%get_ncols()

        write (*, '("spmm", A, "rows of A,", A, "full B", A, "M=", I0, " x N=", I0)')&
        achar(9), achar(9), achar(9), ma, nb

        call c%allocate(ma, nb)
        c%irp(1) = 1
        !$omp parallel private(vals,col_inds,nnz,rwnz,thread_upperbound, acc, start_idx, tic, toc) shared(a,b,c)
        thread_upperbound = 0
        start_idx = 0
        !$omp do schedule(static) private(irw, jj, j)
        do irw = 1, ma
            if (start_idx == 0) start_idx = irw
            do jj = a%irp(irw), a%irp(irw + 1) - 1
                j = a%ja(jj)
                thread_upperbound = thread_upperbound + b%irp(j+1) - b%irp(j)
            end do
        end do
        !$omp end do
        
        ! dense accumulator
        call psb_realloc(nb, acc, info)
        call psb_realloc(thread_upperbound, vals, info)
        call psb_realloc(thread_upperbound, col_inds, info)
        
        ! possible bottleneck
        acc = 0
        
        nnz = 0
        !$omp do schedule(static) private(irw, jj)
        do irw = 1, ma
            rwnz = 0
            do jj = a%irp(irw), a%irp(irw + 1) - 1
                j = a%ja(jj)
                do k = b%irp(j), b%irp(j + 1) - 1
                    if (acc(b%ja(k)) == 0) then
                        nnz = nnz + 1
                        rwnz = rwnz + 1
                        col_inds(nnz) =  b%ja(k)
                    end if
                    acc(b%ja(k)) = acc(b%ja(k)) + a%val(jj) * b%val(k) 
                end do
            end do
            call psb_qsort(col_inds(nnz - rwnz + 1:nnz))
            
            do k = nnz - rwnz + 1, nnz
                vals(k) = acc(col_inds(k))
                acc(col_inds(k)) = 0
            end do
            c%irp(irw + 1) = rwnz
        end do
        !$omp end do
        
        !$omp barrier
        
        ! possible bottleneck
        !$omp single
        do irw = 1, ma
            c%irp(irw + 1) = c%irp(irw + 1) + c%irp(irw)
        end do
        call psb_realloc(c%irp(ma + 1), c%val, info)
        call psb_realloc(c%irp(ma + 1), c%ja, info)
        !$omp end single
        !$omp barrier

        c%val(c%irp(start_idx):c%irp(start_idx) + nnz) = vals(1:nnz)
        c%ja(c%irp(start_idx):c%irp(start_idx) + nnz) = col_inds(1:nnz)
        !$omp end parallel
    end subroutine
end module spmm_mod