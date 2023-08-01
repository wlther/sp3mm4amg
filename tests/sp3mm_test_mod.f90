module sp3mm_test_mod
    use psb_base_mod
    use sp3mm_base_mod
    use omp_lib
    implicit none
    
contains
    function spmm_is_eq(a,b,eps) result(is_eq)
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        real(psb_dpk_), optional, intent(in) :: eps
        logical :: is_eq

        integer(psb_ipk_) :: i
        real(psb_dpk_) :: eps_

        if (present(eps)) then
            eps_ = eps
        else
            eps_ = 1e-6
        end if

        is_eq = .false.

        if (a%get_ncols() /= b%get_ncols()) then
            print *, 'incompatible number or columns', a%get_ncols(), b%get_ncols()
            return
        end if

        if (a%get_nrows() /= b%get_nrows()) then
            print *, 'incompatible number of rows', a%get_nrows(), b%get_nrows()
            return
        end if

        if (a%get_nzeros() /= b%get_nzeros()) then
            print *, 'incompatible number of non zeros', a%get_nzeros(), b%get_nzeros()
            return
        end if

        do i=1, a%get_nrows() + 1
            if (a%irp(i) /= b%irp(i)) then
                print *, 'incompatible irp values', a%irp(i), b%irp(i)
                return
            end if
        end do

        do i=1, a%get_nzeros()
            if (a%ja(i) /= b%ja(i)) then
                print *, 'incompatible ja values', a%ja(i), b%ja(i)
                return
            end if
            if (abs(a%val(i) - b%val(i)) >= eps_) then
                print *, 'incompatible data values', a%val(i), b%val(i)
                return
            end if
        end do

        is_eq = .true.
        return
    end function spmm_is_eq

    subroutine test_sp3mm_pair_serial(r,ac,p,oracle,cfg,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: r,ac,p,oracle
        type(sp3mm_config), intent(in) :: cfg
        integer(psb_ipk_), intent(out) :: info

        real(8) :: end,start
        type(psb_d_csr_sparse_mat) :: rac, out_to_check

        start = omp_get_wtime()
        call spmm_serial(r, ac, rac, cfg, info)
        if (info /= 0) return
        call spmm_serial(rac, p, out_to_check, cfg, info)
        if (info /= 0) return
        end = omp_get_wtime()
        print *, 'sp3mm as pair of spmm_serial :', end - start

        if (.not. spmm_is_eq(out_to_check, oracle)) info = 1

        return
    end subroutine test_sp3mm_pair_serial

    subroutine test_sp3mm_pair_ub(r,ac,p,oracle,cfg,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: r,ac,p,oracle
        type(sp3mm_config), intent(in) :: cfg
        integer(psb_ipk_), intent(out) :: info

        real(8) :: end,start
        type(psb_d_csr_sparse_mat) :: rac, out_to_check

        start = omp_get_wtime()
        call spmm_row_by_row(r, ac, rac, cfg, info)
        if (info /= 0) return
        call spmm_row_by_row(rac, p, out_to_check, cfg, info)
        if (info /= 0) return
        end = omp_get_wtime()
        print *, 'sp3mm as pair of spmm_row_by_row_ub :', end - start


        if (.not. spmm_is_eq(out_to_check, oracle)) info = 1

        return
    end subroutine test_sp3mm_pair_ub
    

    subroutine test_sp3mm_pair_1D_block_ub(r,ac,p,oracle,cfg,info)
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: r,ac,p,oracle
        type(sp3mm_config), intent(in) :: cfg
        integer(psb_ipk_), intent(out) :: info

        real(8) :: end,start
        type(psb_d_csr_sparse_mat) :: rac, out_to_check

        start = omp_get_wtime()
        call spmm_row_by_row_1D_blocks(r, ac, rac, cfg, info)
        if (info /= 0) return
        call spmm_row_by_row_1D_blocks(rac, p, out_to_check, cfg, info)
        if (info /= 0) return
        end = omp_get_wtime()
        print *, 'sp3mm as pair of spmm_row_by_row_1D_block :', end - start


        if (.not. spmm_is_eq(out_to_check, oracle)) info = 1

        return
    end subroutine test_sp3mm_pair_1D_block_ub
    
end module sp3mm_test_mod