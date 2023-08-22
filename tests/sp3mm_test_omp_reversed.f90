program sp3mm_test_omp
    use omp_lib
    use psb_base_mod
    use sp3mm_test_mod
    implicit none
    
    character(len = 3) :: level
    character(len = 32) :: program_name, collection, size, smoothing
    character(len = 128) :: input_argument_1, input_argument_2, &
                            input_argument_3, input_argument_4, &
                            collection_size_smoothing
    real(8) :: preparing_time, spmm1_time, spmm2_time
    type(psb_d_csr_sparse_mat) :: r, ac, p, oracle_out, acp, out_to_check
    integer(psb_ipk_) :: info, i
    
    call get_command_argument(0, program_name)
    if (command_argument_count() < 3) then
        write (*, '("USAGE: ", A, "R_{i+1}, AC_{i}, P_{i+1}, AC_{i+1}")') trim(program_name)
        error stop "incorrect arguments count"
    end if
    
    preparing_time = omp_get_wtime()
    call get_test_matrices(r, ac, p, oracle_out)
    preparing_time = omp_get_wtime() - preparing_time
    
    call get_command_argument(1, input_argument_1)
    call get_command_argument(2, input_argument_2)
    call get_command_argument(3, input_argument_3)
    if (command_argument_count() < 4) then
        input_argument_4 = ''
    else
        call get_command_argument(4, input_argument_4)
    end if

    spmm1_time = omp_get_wtime()
    call spmm_omp_1d(ac, p, acp, info)
    spmm1_time = omp_get_wtime() - spmm1_time
    if (info /= 0) error stop
    spmm2_time = omp_get_wtime()
    call spmm_omp_1d(r, acp, out_to_check,info)
    spmm2_time = omp_get_wtime() - spmm2_time
    if (info /= 0) error stop

    if (.not. spmm_is_eq(out_to_check, oracle_out)) error stop

    if (command_argument_count() >= 5) then
        if ('++' == trim(input_argument_4)) then
            if (command_argument_count() < 6) then
                write (*, '("USAGE: ", A, "R_{i+1}, AC_{i}, P_{i+1}, AC_{i+1} # Collection/Size/Smoothing level")') &
                        trim(program_name)
                error stop "incorrect arguments count"
            end if
            call get_command_argument(5, collection_size_smoothing)
            call get_command_argument(6, level)
        else if (command_argument_count() >= 7) then
            call get_command_argument(6, collection_size_smoothing)
            call get_command_argument(7, level)
        end if

        do i = 1, len_trim(collection_size_smoothing)
            if (collection_size_smoothing(i:i) == '/') collection_size_smoothing(i:i) = ' '
        end do

        read (collection_size_smoothing, *) collection, size, smoothing

        write (*, '(A,",",A,",",A,",",A,",",A,",",ES12.6,",",I0,",",I0,",",I0,",",I0,",",ES12.6,",",I0,",",I0,",",I0,",",ES12.6)') &
                                    "gustavson_reversed", trim(collection), &
                                    trim(size), trim(smoothing),&
                                    trim(level), preparing_time,&
                                    omp_get_max_threads(),&
                                    acp%get_nrows(), acp%get_ncols(),&
                                    acp%get_nzeros(), spmm1_time,&
                                    out_to_check%get_nrows(),&
                                    out_to_check%get_ncols(),&
                                    out_to_check%get_nzeros(),&
                                    spmm2_time
    else
        write (*, '("#### ", A, " ####")') trim(program_name)
        write (*, '("# R: ", A, A, A)') achar(9), achar(9), trim(input_argument_1)
        write (*, '("# AC: ", A, A, A)') achar(9), achar(9), trim(input_argument_2)
        write (*, '("# P: ", A, A, A)') achar(9), achar(9), trim(input_argument_3)
        write (*, '("# AC_{i+1}: ", A, A)') achar(9), trim(input_argument_4)

        write (*,*)
        write (*, '("omp_get_max_threads:", A, I0)') achar(9), omp_get_max_threads()
        
        write (*,*)
        write (*,'("preparing time: ", ES12.6)') preparing_time

        write (*,*)
        write (*, '("AC * P: ", ES12.6)') spmm1_time

        write (*,*)
        write (*, '("R * ACP: ", ES12.6)') spmm2_time
    end if
end program sp3mm_test_omp