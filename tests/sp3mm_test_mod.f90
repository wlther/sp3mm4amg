module sp3mm_test_mod
    use psb_base_mod
    use psb_util_mod
    use omp_lib
    use spmm_mod
    implicit none
    integer(psb_ipk_), parameter :: iunit=12
    character(len=3), parameter :: out_fmt = 'CSR'

contains
    subroutine get_test_matrices(r, ac, p, oracle_out)
        type(psb_d_csr_sparse_mat), intent(out):: r, ac, p, oracle_out
        
        type(psb_d_csr_sparse_mat):: rac
        integer(psb_ipk_) :: info
        character(len = 128) :: input_argument_1, input_argument_2, &
        input_argument_3, input_argument_4
        type(psb_dspmat_type) :: tmp

        !$omp parallel private(tmp)
        !$omp single
        call get_command_argument(1, input_argument_1)
        call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_1)
        call tmp%cscnv(info, type=out_fmt)
        call r%mv_from_fmt(tmp%a, info)
        if (info /= 0) then
            error stop 'Error during conversion MM -> CSR of R'
        end if
        !$omp end single

        !$omp single 
        call get_command_argument(2, input_argument_2)
        call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_2)
        call tmp%cscnv(info, type=out_fmt)
        call ac%mv_from_fmt(tmp%a, info)
        if (info /= 0) then
            error stop 'Error during conversion MM -> CSR of AC_{i}'
        end if 
        !$omp end single

        !$omp single 
        call get_command_argument(3, input_argument_3)
        call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_3)
        call tmp%cscnv(info, type=out_fmt)
        call p%mv_from_fmt(tmp%a, info)
        if (info /= 0) then
            error stop 'Error during conversion MM -> CSR of P'
        end if
        !$omp end single

        !$omp single
        if (command_argument_count() > 3) then
            call get_command_argument(4, input_argument_4)
            if (input_argument_4 == '++') then
                input_argument_4 = ''
                write (*, '(A)') 'No oracle was given to compare results of Sp3mm : generating with spmm_serial'
                call psb_dcsrspspmm(r, ac, rac, info)
                call psb_dcsrspspmm(rac, p, oracle_out, info)
            else
                call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_4)
                call tmp%cscnv(info, type=out_fmt)
                call oracle_out%mv_from_fmt(tmp%a, info)
                if (info /= 0) then
                    error stop 'Error during conversion MM -> CSR of AC_{i+1}'
                end if
            end if
        else
            input_argument_4 = ''
            write (*, '(A)') 'No oracle was given to compare results of Sp3mm : generating with spmm_serial'
    
        end if
        !$omp end single
        !$omp end parallel

        ! consistency check (checking dimensions)
        if (r%get_ncols() /= ac%get_nrows()) then
            error stop 'Error: incompatible dimensions between r and ac'
        end if
        if (ac%get_ncols() /= p%get_nrows()) then
            error stop 'Error: incompatible dimensions between ac and p'
        end if
        if (command_argument_count() > 3) then
            if (r%get_nrows() /= oracle_out%get_nrows()) then
                error stop 'Error: oracle has an incorrect number of rows'
            end if
            if (p%get_ncols() /= oracle_out%get_ncols()) then
                error stop 'Error: oracle has an incorrect number of columns'
            end if
        end if
    end subroutine get_test_matrices

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
                print *, 'incompatible data values at index', i, a%val(i), b%val(i)
                return
            end if
        end do

        is_eq = .true.
        return
    end function spmm_is_eq

end module sp3mm_test_mod
