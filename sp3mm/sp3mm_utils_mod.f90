module sp3mm_utils_mod
    implicit none
    
contains
    subroutine init_urndfd(info)
        implicit none
        integer, intent(out) :: info
        open(unit=9, file='/dev/urandom', iostat=info)
    end subroutine init_urndfd

    subroutine extract_in_tmp_fs(path, tmp_fs_decompress_path, info)
        implicit none
        character(len = *), intent(in) :: path,tmp_fs_decompress_path
        integer, intent(out) :: info

        ! file extension
        character(len=:), allocatable :: ext
        integer :: dot_idx
        character(len=:), allocatable :: command
        dot_idx = index(path, '.')

        if (dot_idx < 1) then
            info = 1
            return
        end if

        ext = path(dot_idx:)

        select case(ext)
        case(".gz")
            command = "gzip -d -c "//trim(path)//" > "//tmp_fs_decompress_path 
        case(".xz")
            command = "xz -d -c "//trim(path)//" > "//tmp_fs_decompress_path 
        case(".bz2")
            command = "bzip2 -d -c "//trim(path)//" > "//tmp_fs_decompress_path 
        case(".zip")
            command = "unzip -c "//trim(path)//" > "//tmp_fs_decompress_path 
        case default
            info = 1
            return
        end select

        call system(command, info)

        if (info /= 0) then
            print *, "Error : decompression failed"
        else
            print *, "Decompression successful"
        end if
    end subroutine extract_in_tmp_fs

    subroutine print_sp3mm_core(r, ac, p, cfg)
        use psb_base_mod
        use sp3mm_config_mod
        use sp3mm_const_mod
        use iso_c_binding
        implicit none
        type(psb_d_csr_sparse_mat), intent(in) :: r,ac,p
        type(sp3mm_config), intent(in) :: cfg

        integer(psb_ipk_) :: limb

        write (*, '("COARSENING AC: ", I0, "x", I0, "  --> ", I0, "x", I0, &
        A, "conf grid: ", I0, "x", I0, ", ", A, "NNZ:", I0, "-", I0, "-", I0, A,&
        "AVG_TIMES_ITERATION:", I0, A, "symbUBAssignType:", A, A, "bitmapLimbSize:", I0)') &
        ac%get_nrows(), ac%get_ncols(), r%get_nrows(), p%get_ncols(), achar(9), &
        cfg%rows, cfg%cols, achar(9), r%get_nzeros(), ac%get_nzeros(), p%get_nzeros(), &
        achar(9), avg_time_iteration, achar(9), 'STATIC_ASSIGN', achar(9), bit_size(limb)

    end subroutine print_sp3mm_core

end module sp3mm_utils_mod