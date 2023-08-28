module spmm_mod
    use omp_lib
    use psb_base_mod
    use idx_tree_mod
    implicit none
    
contains

    subroutine spmm_serial(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info

        call psb_dcsrspspmm(a,b,c,info)
    end subroutine spmm_serial

    ! gustavson's algorithm using perfect hashing 
    ! and OpenMP parallelisation
    subroutine spmm_omp_gustavson(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info
        
        real(psb_dpk_), allocatable :: vals(:), acc(:)
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
        !$omp parallel private(vals,col_inds,nnz,rwnz,thread_upperbound,acc,start_idx,end_idx) &
        !$omp shared(a,b,c,offsets) 
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
        
        c%val(c%irp(start_idx):c%irp(end_idx + 1) - 1) = vals(1:nnz)
        c%ja(c%irp(start_idx):c%irp(end_idx + 1) - 1) = col_inds(1:nnz)
        deallocate(acc)
        !$omp end parallel
        deallocate(offsets)
    end subroutine spmm_omp_gustavson

    subroutine spmm_omp_gustavson_1d(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info
        
        real(psb_dpk_), allocatable :: vals(:), acc(:)
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
        
        c%val(c%irp(start_idx):c%irp(end_idx + 1) - 1) = vals(1:nnz)
        c%ja(c%irp(start_idx):c%irp(end_idx + 1) - 1) = col_inds(1:nnz)
        !$omp end parallel
        deallocate(offsets)
    end subroutine spmm_omp_gustavson_1d

    subroutine spmm_serial_rb_tree(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info

        integer(psb_ipk_)   :: a_m, b_n
        integer(psb_ipk_)   :: row, col
        type(idx_tree), allocatable :: row_accs(:)

        a_m = a%get_nrows()
        b_n = b%get_ncols()

        call c%allocate(a_m, b_n)
        allocate(row_accs(a_m))

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
    end subroutine spmm_serial_rb_tree

    subroutine spmm_omp_rb_tree(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info

        integer(psb_ipk_)   :: a_m, b_n
        integer(psb_ipk_)   :: row, col
        type(idx_tree), allocatable :: row_accs(:)

        a_m = a%get_nrows()
        b_n = b%get_ncols()

        call c%allocate(a_m, b_n)
        allocate(row_accs(a_m))

        !$omp parallel do schedule(static)
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
    end subroutine spmm_omp_rb_tree

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

    subroutine spmm_omp_two_pass(a,b,c,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out):: c
        integer(psb_ipk_), intent(out) :: info

        integer(psb_ipk_)   :: a_m, b_n, row, col

        a_m = a%get_nrows()
        b_n = b%get_ncols()

        call c%allocate(a_m, b_n)
        
        call compute_indices(a, b, c, info)

        !$omp parallel do schedule(static)
        do row = 1, a_m
            do col = a%irp(row), a%irp(row + 1) - 1
                call direct_scalar_sparse_row_mul(c, row, a%val(col), b, a%ja(col))
            end do 
        end do
        !$omp end parallel do
    end subroutine spmm_omp_two_pass

end module spmm_mod