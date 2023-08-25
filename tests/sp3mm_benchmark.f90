program sp3mm_benchmark
    use omp_lib
    use psb_base_mod
    use sp3mm_test_mod
    implicit none

    interface
    subroutine spmm_func(a,b,c,info_)
        import psb_d_csr_sparse_mat
        import psb_ipk_
        type(psb_d_csr_sparse_mat), intent(in) :: a,b
        type(psb_d_csr_sparse_mat), intent(out) :: c
        integer(psb_ipk_), intent(out) :: info_
    end subroutine spmm_func
    end interface

    integer :: num_iterations
    logical :: thread_iterations
    character(32) :: collection, &
                     size, &
                     smoothing, &
                     level
    character(128):: collection_size_smoothing, input_argument_4
    integer :: num_implementations
    character(32), allocatable :: implementations(:)
    integer :: i, j, max_threads, num_threads
    integer(psb_ipk_) :: info, nrows, ncols, nnz
    logical :: reversed
    real(8) :: preparing_time, spmm1_time, spmm2_time
    type(psb_d_csr_sparse_mat) :: r, ac, p, oracle_out, rac, acp, out_to_check
    procedure(spmm_func), pointer :: spmm

    open(unit=10, file='test_settings.conf', status='old')
    read (10, *) num_iterations
    read (10, *) thread_iterations
    read (10, *) num_implementations

    allocate(implementations(num_implementations))

    do i = 1, num_implementations
        read (10, *) implementations(i)
    end do

    close(10)

    preparing_time = omp_get_wtime()
    call get_test_matrices(r, ac, p, oracle_out)
    preparing_time = omp_get_wtime() - preparing_time

    if (command_argument_count() >= 5) then
        call get_command_argument(4, input_argument_4)
        if ('++' == trim(input_argument_4)) then
            if (command_argument_count() < 6) then
                write (*, *) "USAGE: ./sp3mm_benchmark R_{i+1}, AC_{i}, P_{i+1}, AC_{i+1} # Collection/Size/Smoothing level"
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
    else
        collection = ""
        size = ""
        smoothing = ""
        level = ""
    end if

    max_threads = omp_get_max_threads()

    do i = 1, num_implementations
        select case(trim(implementations(i)))
        case("serial")
            reversed=.false.
            spmm => spmm_serial
        case("serial_reversed")
            reversed=.true.
            spmm => spmm_serial
        case("omp_gustavson")
            reversed = .false.
            spmm => spmm_omp_gustavson
        case("omp_gustavson_reversed")
            reversed = .true.
            spmm => spmm_omp_gustavson
        case("omp_gustavson_1d")
            reversed = .false.
            spmm => spmm_omp_gustavson_1d
        case("omp_gustavson_1d_reversed")
            reversed = .true.
            spmm => spmm_omp_gustavson_1d
        case("serial_rb_tree")
            reversed = .false.
            spmm => spmm_serial_rb_tree
        case("serial_rb_tree_reversed")
            reversed = .true.
            spmm => spmm_serial_rb_tree
        case("omp_rb_tree")
            reversed = .false.
            spmm => spmm_omp_rb_tree
        case("omp_rb_tree_reversed")
            reversed = .true.
            spmm => spmm_omp_rb_tree
        case("omp_two_pass")
            reversed = .false.
            spmm => spmm_omp_rb_tree
        case("omp_two_pass_reversed")
            reversed = .true.
            spmm => spmm_omp_rb_tree
        case default
            cycle
        end select

        num_threads = 1
        do while (num_threads <= max_threads)
            if (thread_iterations) then
                call omp_set_num_threads(num_threads)
                num_threads = 2 * num_threads
            else
                num_threads = max_threads + 1
            end if

            do j = 1, num_iterations
                if (reversed) then
                    spmm1_time = omp_get_wtime()
                    call spmm(ac, p, acp, info)
                    spmm1_time = omp_get_wtime() - spmm1_time
                    if (info /= 0) error stop
                    spmm2_time = omp_get_wtime()
                    call spmm(r, acp, out_to_check,info)
                    spmm2_time = omp_get_wtime() - spmm2_time
                    if (info /= 0) error stop
                    nrows = acp%get_nrows()
                    ncols = acp%get_ncols()
                    nnz = acp%get_nzeros()
                else 
                    spmm1_time = omp_get_wtime()
                    call spmm(r, ac, rac, info)
                    spmm1_time = omp_get_wtime() - spmm1_time
                    if (info /= 0) error stop
                    spmm2_time = omp_get_wtime()
                    call spmm(rac, p, out_to_check,info)
                    spmm2_time = omp_get_wtime() - spmm2_time
                    if (info /= 0) error stop
                    nrows = rac%get_nrows()
                    ncols = rac%get_ncols()
                    nnz = rac%get_nzeros()
                end if
                write (*, '(A,",",A,",",A,",",A,",",A,",",ES12.6,",",I0,",",I0,",",I0,",",&
                        I0,",",ES12.6,",",I0,",",I0,",",I0,",",ES12.6)') &
                                    trim(implementations(i)), trim(collection), &
                                    trim(size), trim(smoothing),&
                                    trim(level), preparing_time,&
                                    omp_get_max_threads(),&
                                    nrows, ncols,&
                                    nnz, spmm1_time,&
                                    out_to_check%get_nrows(),&
                                    out_to_check%get_ncols(),&
                                    out_to_check%get_nzeros(),&
                                    spmm2_time
            end do
        end do
    end do
end program sp3mm_benchmark