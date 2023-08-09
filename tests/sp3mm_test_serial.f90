program sp3mm_test_serial
    use omp_lib
    use psb_base_mod
    use sp3mm_test_mod
    implicit none
    
    character(len = 32) :: program_name
    character(len = 128) :: input_argument_1, input_argument_2, &
    input_argument_3, input_argument_4
    real(8) :: tic, toc
    type(psb_d_csr_sparse_mat) :: r, ac, p, oracle_out, rac, out_to_check
    type(psb_d_csc_sparse_mat) :: csc
    integer(psb_ipk_) :: info
    
    call get_command_argument(0, program_name)
    if (command_argument_count() < 3) then
        write (*, '("USAGE: ", A, "R_{i+1}, AC_{i}, P_{i+1}, AC_{i+1}")') trim(program_name)
        error stop "incorrect arguments count"
    end if
    
    tic = omp_get_wtime()
    call get_test_matrices(r, ac, p, oracle_out)
    toc = omp_get_wtime()
    
    call get_command_argument(1, input_argument_1)
    call get_command_argument(2, input_argument_2)
    call get_command_argument(3, input_argument_3)
    write (*, '("#### ", A, " ####")') trim(program_name)
    write (*, '("# R: ", A, A, A)') achar(9), achar(9), trim(input_argument_1)
    write (*, '("# AC: ", A, A, A)') achar(9), achar(9), trim(input_argument_2)
    write (*, '("# P: ", A, A, A)') achar(9), achar(9), trim(input_argument_3)
    if (command_argument_count() < 4) then
        write (*, '("# AC_{i+1} automatically generated")')
    else
    call get_command_argument(4, input_argument_4)
        write (*, '("# AC_{i+1}: ", A, A)') achar(9), trim(input_argument_4)
    end if

    write (*,*)
    write (*, '("omp_get_max_threads:", A, I0)') achar(9), omp_get_max_threads()
    
    write (*,*)
    write (*,'("preparing time: ", ES12.6)') toc - tic

    tic = omp_get_wtime()
    call psb_dcsrspspmm(r, ac, rac, info)
    if (info /= 0) error stop
    call psb_dcsrspspmm(rac, p, out_to_check,info)
    if (info /= 0) error stop
    toc = omp_get_wtime()

    write (*, '(A, ": ", ES12.6)') trim(program_name), toc - tic

    if (.not. spmm_is_eq(out_to_check, oracle_out)) error stop

end program sp3mm_test_serial