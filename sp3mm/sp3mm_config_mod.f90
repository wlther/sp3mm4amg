module sp3mm_config_mod
    use psb_base_mod
    implicit none

    type :: sp3mm_config

    integer(psb_ipk_) :: rows, cols
    integer(psb_ipk_) :: symb_mm_row_impl_id
    integer(psb_ipk_) :: thread_num
    
    procedure(chunk_distrib), pointer :: chunk_distrib_func

    contains
    
    procedure :: set_chunk_distrib_func => set_chunk_distrib_func_impl

    end type sp3mm_config

    abstract interface
    subroutine chunk_distrib(cfg, r, info)
        use psb_base_mod
        import sp3mm_config
        class(sp3mm_config), intent(in) ::cfg
        integer(psb_ipk_), intent(in) :: r
        integer(psb_ipk_), intent(out) :: info
    end subroutine chunk_distrib
    end interface

    contains 

    subroutine chunks_noop(cfg, r, info)
        ! does nothing
        class(sp3mm_config), intent(in) :: cfg
        integer(psb_ipk_), intent(in) :: r
        integer(psb_ipk_), intent(out) :: info
        info = 0
    end subroutine

    subroutine chunks_fair(cfg, r, info)
        use omp_lib
        use omp_lib_kinds
        class(sp3mm_config), intent(in) :: cfg
        integer(psb_ipk_), intent(in) :: r
        integer(psb_ipk_), intent(out) :: info

        integer(omp_sched_kind) :: k, kind
        integer(psb_ipk_) :: chunk_size, chunk_size_new, monotonic
        integer(psb_ipk_), dimension(3) :: kind_chunk_monotonic

        chunk_size_new = 0

        if (cfg%thread_num <= 0) then
            print *, 'Error: thread number should be a positive integer greater than zero'
            info = 1
            return
        end if

        call omp_get_schedule(kind, chunk_size)

        k = kind

        select case(k)
        case (omp_sched_static)
            ! static is always fair
            info = 0
            return
        case (omp_sched_dynamic)
            chunk_size_new = max(1, r/(cfg%thread_num))
            info = 0
        ! missing omp_sched_auto and omp_sched_guided
        ! (not present in the c code either)
        end select

        if (chunk_size == chunk_size_new) return

        call omp_set_schedule(kind, chunk_size_new)
        
        write (*, '("chunksFair:", A, "chunk adapted to ", I0)')&
        achar(9), chunk_size_new
        call omp_get_runtime_schedule(kind_chunk_monotonic)

    end subroutine

    subroutine chunks_fair_folded(cfg, r, info)
        use omp_lib
        use omp_lib_kinds
        class(sp3mm_config), intent(in) :: cfg
        integer(psb_ipk_), intent(in) :: r
        integer(psb_ipk_), intent(out) :: info

        integer(omp_sched_kind) :: k, kind
        integer(psb_ipk_) :: chunk_size, chunk_size_new, monotonic
        integer(psb_ipk_), dimension(3) :: kind_chunk_monotonic

        chunk_size_new = 0

        if (cfg%thread_num <= 0) then
            print *, 'Error: thread number should be a positive integer(psb_ipk_) greater than zero'
            info = 1
            return
        end if

        call omp_get_schedule(kind, chunk_size)

        k = kind

        select case(k)
        case (omp_sched_static)
            ! static is always fair
            info = 0
            return
        case (omp_sched_dynamic)
            chunk_size_new = max(1, r/(cfg%thread_num * 4))
            info = 0
        ! missing omp_sched_auto and omp_sched_guided
        ! (not present in the c code either)
        end select

        if (chunk_size == chunk_size_new) return

        call omp_set_schedule(kind, chunk_size_new)
        
        write (*, '("chunksFairFolded:", A, "chunk adapted to ", I0)')&
        achar(9), chunk_size_new
        call omp_get_runtime_schedule(kind_chunk_monotonic)

    end subroutine chunks_fair_folded

    subroutine omp_get_runtime_schedule(kind_chunk_monotonic)
        use omp_lib
        use omp_lib_kinds
        use sp3mm_const_mod
        integer(psb_ipk_), dimension(3), intent(out) :: kind_chunk_monotonic
        
        integer(omp_sched_kind) :: kind
        integer(psb_ipk_) :: chunk_size, monotonic
        character(len=17) :: schedule
        character(len=1) :: monotonic_status

        call omp_get_schedule(kind, chunk_size)

        monotonic = 0

        select case (kind)
        case (1)
            schedule = 'OMP_SCHED_STATIC'
        case (2)
            schedule = 'OMP_SCHED_DYNAMIC'
        case (3)
            schedule = 'OMP_SCHED_GUIDED'
        case (4)
            schedule = 'OMP_SCHED_AUTO'
        end select

        if (monotonic /= 0)  then
            monotonic_status = 'Y'
        else
            monotonic_status = 'N'
        end if

        write (*, '("omp sched gather:", A, "kind:", A, A, "omp chunkSize:", I0, &
         A, "monotonic:", A, A, "fairChunkFolding:", I0)') &
        achar(9), schedule, achar(9), chunk_size, achar(9), monotonic_status, &
        achar(9), fair_chunk_folding

        kind_chunk_monotonic(1) = kind
        kind_chunk_monotonic(2) = chunk_size
        kind_chunk_monotonic(3) = monotonic
    end subroutine omp_get_runtime_schedule

    subroutine set_chunk_distrib_func_impl(cfg, implementation_type)
        class(sp3mm_config), intent(inout) :: cfg
        integer(psb_ipk_), intent(in) :: implementation_type

        select case(implementation_type)
        case (1)
            cfg%chunk_distrib_func => chunks_fair
        case (2)
            cfg%chunk_distrib_func => chunks_fair_folded
        case default
            cfg%chunk_distrib_func => chunks_noop
        end select
    end subroutine set_chunk_distrib_func_impl

    subroutine get_config(cfg, info)
        class(sp3mm_config), intent(inout) :: cfg
        integer(psb_ipk_), intent(out) :: info

        character(len=64) :: str_value
        integer(psb_ipk_) :: len

        info = 0
        
        call get_environment_variable('GRID_COLS', value=str_value, length=len)
        if (len /= 0) then
            read(str_value, *, iostat = info) cfg%cols
            info = 1
        end if

        call get_environment_variable('GRID_ROWS',value=str_value, length=len)
        if (len /= 0) then
            read(str_value, *, iostat = info) cfg%rows
            info = 1
        end if
    end subroutine get_config
end module sp3mm_config_mod
