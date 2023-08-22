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
        
        real(psb_dpk_), allocatable :: vals(:), acc(:)
        real(8) :: tic, toc
        integer(psb_ipk_)   :: ma, nb
        integer(psb_ipk_), allocatable :: col_inds(:), offsets(:)
        integer(psb_ipk_)   :: irw, jj, j, k, nnz, rwnz, thread_upperbound, start_idx, end_idx

        ma = a%get_nrows()
        nb = b%get_ncols()

        call c%allocate(ma, nb)
        c%irp(1) = 1

        ! dense accumulator
        ! https://sc18.supercomputing.org/proceedings/workshops/workshop_files/ws_lasalss115s2-file1.pdf
        call psb_realloc(nb, acc, info)
        
        allocate(offsets(omp_get_max_threads()))
        !$omp parallel private(vals,col_inds,nnz,rwnz,thread_upperbound,acc,start_idx,end_idx) shared(a,b,c,offsets)
        thread_upperbound = 0
        start_idx = 0
        !$omp do schedule(static) private(irw, jj, j)
        do irw = 1, ma
            if (start_idx == 0) then
                start_idx = irw
            end if
            end_idx = irw
            do jj = a%irp(irw), a%irp(irw + 1) - 1
                j = a%ja(jj)
                thread_upperbound = thread_upperbound + b%irp(j+1) - b%irp(j)
            end do
        end do
        !$omp end do
        
        call psb_realloc(thread_upperbound, vals, info)
        call psb_realloc(thread_upperbound, col_inds, info)
        
        ! possible bottleneck
        acc = 0
        
        nnz = 0
        !$omp do schedule(static) private(irw, jj, j, k)
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
        
        offsets(omp_get_thread_num() + 1) = nnz
        !$omp barrier
        
        ! possible bottleneck
        !$omp single
        do k = 1, omp_get_num_threads() - 1
            offsets(k + 1) = offsets(k + 1) + offsets(k)
        end do
        !$omp end single

        !$omp barrier

        if (omp_get_thread_num() /= 0) then
            c%irp(start_idx) = offsets(omp_get_thread_num()) + 1
        end if

        do irw = start_idx, end_idx - 1
            c%irp(irw + 1) = c%irp(irw + 1) + c%irp(irw)
        end do

        !$omp barrier

        !$omp single
        c%irp(ma + 1) = c%irp(ma + 1) + c%irp(ma)
        call psb_realloc(c%irp(ma + 1), c%val, info)        
        call psb_realloc(c%irp(ma + 1), c%ja, info)
        !$omp end single
        
        c%val(c%irp(start_idx):c%irp(start_idx) + nnz) = vals(1:nnz)
        c%ja(c%irp(start_idx):c%irp(start_idx) + nnz) = col_inds(1:nnz)
        !$omp end parallel
    end subroutine spmm_omp

    subroutine spmm_omp_1d(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info
        
        real(psb_dpk_), allocatable :: vals(:), acc(:)
        real(8) :: tic, toc
        integer(psb_ipk_)   :: ma, nb
        integer(psb_ipk_), allocatable :: col_inds(:), offsets(:)
        integer(psb_ipk_)   :: irw, jj, j, k, nnz, rwnz, thread_upperbound, &
                                start_idx, end_idx , blk, blk_size, rwstart,&
                                rwblk, rwblkrem, nblks

        ma = a%get_nrows()
        nb = b%get_ncols()

        call c%allocate(ma, nb)
        c%irp(1) = 1

        ! dense accumulator
        ! https://sc18.supercomputing.org/proceedings/workshops/workshop_files/ws_lasalss115s2-file1.pdf
        call psb_realloc(nb, acc, info)
        allocate(offsets(omp_get_max_threads()))

        nblks = 4 * omp_get_max_threads()
        rwblk = (ma / nblks)
        rwblkrem = modulo(ma, nblks)
        !$omp parallel private(vals,col_inds,nnz,thread_upperbound,acc,start_idx,end_idx) shared(a,b,c,offsets)
        thread_upperbound = 0
        start_idx = 0
        !$omp do schedule(static) private(irw, jj, j)
        do irw = 1, ma
            do jj = a%irp(irw), a%irp(irw + 1) - 1
                j = a%ja(jj)
                thread_upperbound = thread_upperbound + b%irp(j+1) - b%irp(j)
            end do
        end do
        !$omp end do
        
        call psb_realloc(thread_upperbound, vals, info)
        call psb_realloc(thread_upperbound, col_inds, info)
        
        ! possible bottleneck
        acc = 0
        
        nnz = 0
        !$omp do schedule(static) private(irw,jj,j,k,rwnz,blk,blk_size,rwstart)
        do blk = 0, nblks - 1
            if (blk < rwblkrem) then
                blk_size = rwblk + 1
                rwstart = blk * rwblk + blk + 1
            else
                blk_size = rwblk
                rwstart = blk * rwblk &
                + rwblkrem + 1
            end if
            do irw = rwstart, rwstart + blk_size - 1
                if (start_idx == 0) then
                    start_idx = irw
                end if
                end_idx = irw
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
        end do
        !$omp end do
        
        offsets(omp_get_thread_num() + 1) = nnz
        !$omp barrier
        
        ! possible bottleneck
        !$omp single
        do k = 1, omp_get_num_threads() - 1
            offsets(k + 1) = offsets(k + 1) + offsets(k)
        end do
        !$omp end single

        !$omp barrier

        if (omp_get_thread_num() /= 0) then
            c%irp(start_idx) = offsets(omp_get_thread_num()) + 1
        end if

        do irw = start_idx, end_idx - 1
            c%irp(irw + 1) = c%irp(irw + 1) + c%irp(irw)
        end do

        !$omp barrier

        !$omp single
        c%irp(ma + 1) = c%irp(ma + 1) + c%irp(ma)
        call psb_realloc(c%irp(ma + 1), c%val, info)        
        call psb_realloc(c%irp(ma + 1), c%ja, info)
        !$omp end single
        
        c%val(c%irp(start_idx):c%irp(start_idx) + nnz) = vals(1:nnz)
        c%ja(c%irp(start_idx):c%irp(start_idx) + nnz) = col_inds(1:nnz)
        !$omp end parallel
    end subroutine spmm_omp_1d
end module spmm_mod