program sp3mm_test
    use omp_lib
    use psb_base_mod
    use psb_util_mod
    use sp3mm_base_mod
    use sp3mm_test_mod
    implicit none

    real(8) :: end,start,elapsed
    type(psb_dspmat_type) :: tmp
    type(psb_d_csr_sparse_mat) :: r, ac, p, oracle_out, out_to_check
    integer(psb_ipk_) :: info
    character(len = 128) :: input_argument_1, input_argument_2, &
    input_argument_3, input_argument_4
    integer(psb_ipk_), parameter :: iunit=12
    character(len=3), parameter :: out_fmt = 'CSR'
    type(sp3mm_config) :: cfg
    integer, dimension(3) :: kind_chunk_monotonic
    
    cfg%rows = 20
    cfg%cols = 2
    cfg%symb_mm_row_impl_id = 0

    if (command_argument_count() < 3) then
        print *, 'USAGE: ./sp3mm_test R_{i+1}, AC_{i}, P_{i+1}, AC_{i+1}'
        goto 9999
    end if

    start = omp_get_wtime()
    call get_command_argument(1, input_argument_1)
    call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_1)
    call tmp%cscnv(info, type=out_fmt)
    call r%mv_from_fmt(tmp%a, info)

    call get_command_argument(2, input_argument_2)
    call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_2)
    call tmp%cscnv(info, type=out_fmt)
    call ac%mv_from_fmt(tmp%a, info)

    call get_command_argument(3, input_argument_3)
    call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_3)
    call tmp%cscnv(info, type=out_fmt)
    call p%mv_from_fmt(tmp%a, info)

    if (command_argument_count() > 3) then
        call get_command_argument(4, input_argument_4)
        call mm_mat_read(tmp,info,iunit=iunit,filename=input_argument_4)
        call tmp%cscnv(info, type=out_fmt)
        call oracle_out%mv_from_fmt(tmp%a, info)
        if (info /= 0) then
            print *, 'Error during conversion MM -> CSR of AC_{i+1}'
            goto 9999
        end if
    else
        input_argument_4 = ''
        print *, 'No oracle was given to compare results of Sp3mm : generating with sp3mmRowByRowPair'
    end if


    ! consistency check (checking dimensions)
    if (r%get_ncols() /= ac%get_nrows()) then
        print *, 'Error: incompatible dimensions between r and ac'
        goto 9999
    end if
    if (ac%get_ncols() /= p%get_nrows()) then
        print *, 'Error: incompatible dimensions between ac and p'
        goto 9999
    end if
    if (command_argument_count() > 3) then
        if (r%get_nrows() /= oracle_out%get_nrows()) then
            print *, 'Error: oracle has an incorrect number of rows'
            goto 9999
        end if
        if (p%get_ncols() /= oracle_out%get_ncols()) then
            print *, 'Error: oracle has an incorrect number of columns'
            goto 9999
        end if
    end if


    call get_config(cfg, info)
    if (info /= 0) then
        print *, 'configuration changed from environnement'
        info = 0
    end if

    cfg%thread_num = omp_get_max_threads()

    write (*, '("omp_get_max_threads:", A, I0)') achar(9), cfg%thread_num

    write (*, '("#", A, A, A, A)') trim(input_argument_1), &
    trim(input_argument_2), trim(input_argument_3), trim(input_argument_4)

    call omp_get_runtime_schedule(kind_chunk_monotonic)
    
    call set_chunk_distrib_func_impl(cfg, 2)

    call print_sp3mm_core(r, ac, p, cfg)

    end = omp_get_wtime()
    elapsed = end - start

    write (*, '("preparing time: ", ES12.6)') elapsed

    ! test SP3MM as pair of SPMM: RAC = R * AC; RACP = RAC * P
    write (*, '(A, "CHECKING UPPERBOUND IMPLEMENTATIONS", A)')&
    achar(27) // "[1;32m", achar(27) // "[0m"

    write (*, '(A, "@computing Sp3MM as pair of serial SpMM", A)')&
    achar(27) // "[1;32m", achar(27) // "[0m"
    start = omp_get_wtime()
    call test_sp3mm_pair_serial(r, ac, p, oracle_out, cfg, info)
    if (info /= 0) then
        print *, 'sp3mm as pair of spmm_serial failed'
        goto 9999
    end if
    end = omp_get_wtime()

    write (*, '(A, "@computing Sp3MM as pair of SpMM UpperBounded", A)')&
    achar(27) // "[1;32m", achar(27) // "[0m"
    start = omp_get_wtime()
    call test_sp3mm_pair_UB(r, ac, p, oracle_out, cfg, info)
    if (info /= 0) then
        print *, 'sp3mm as pair of spmm_row_by_row_ub failed'
        goto 9999
    end if
    end = omp_get_wtime()

    print *, 'sp3mm as pair of spmm_row_by_row_ub :', end - start
9999 continue
end program sp3mm_test