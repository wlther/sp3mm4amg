module sp3mm_acc_mod
    use psb_base_mod
    use sp3mm_const_mod
    implicit none

    type :: dense_acc
        ! values
        real(psb_dpk_), allocatable :: dense_val(:)
        ! col indices
        integer(psb_ipk_), allocatable :: nz_idxs(:)
        ! number of non zero elements
        integer :: nz
        ! bitmap
        integer(psb_ipk_), allocatable :: bitset(:)

        contains

        procedure :: init   => init_dense_acc
        procedure :: reset  => reset_dense_acc
    end type dense_acc

    contains

    subroutine init_dense_acc(sp_acc, size)
        implicit none
        class(dense_acc), intent(inout)  :: sp_acc
        integer(psb_ipk_), intent(in)   :: size
        
        integer(psb_ipk_) :: info

        ! call psb_realloc(size, sp_acc%dense_val, info)
        allocate(sp_acc%dense_val(size))
        allocate(sp_acc%nz_idxs(size))
        allocate(sp_acc%bitset(size))

        sp_acc%dense_val = 0
        sp_acc%bitset = 0
        sp_acc%nz = 0

    end subroutine init_dense_acc

    subroutine reset_dense_acc(sp_acc)
        implicit none
        class(dense_acc), intent(inout) :: sp_acc

        sp_acc%dense_val = 0
        sp_acc%bitset = 0
        sp_acc%nz = 0
    end subroutine reset_dense_acc

    subroutine scalar_sparse_row_mul(sp_acc, scalar, mat, row_num)
        implicit none
        type(dense_acc), intent(inout)  :: sp_acc
        real(psb_dpk_) :: scalar
        type(psb_d_csr_sparse_mat), intent(in) :: mat
        integer(psb_ipk_) :: row_num

        integer(psb_ipk_) :: i, j

        do i = mat%irp(row_num), mat%irp(row_num + 1) - 1
            j = mat%ja(i)

            sp_acc%dense_val(j) = sp_acc%dense_val(j) + scalar * mat%val(i)

            if (.not. idx_in_bitset(j, sp_acc)) then
                sp_acc%nz_idxs(sp_acc%nz) = j
            end if
        end do
    end subroutine scalar_sparse_row_mul

    function idx_in_bitset(idx, acc) result(is_in_bitset)
        integer(psb_ipk_), intent(in):: idx
        type(dense_acc), intent(inout) :: acc
        logical :: is_in_bitset

        integer(psb_ipk_) :: limb_idx, bit_idx, mask

        limb_idx = (idx - 1) / bit_size(acc%bitset(1)) + 1

        bit_idx = mod(idx - 1, bit_size(acc%bitset(1))) + 1

        mask = 2 ** (bit_idx - 1)

        is_in_bitset = iand(acc%bitset(limb_idx), mask) /= 0

        if (.not. is_in_bitset) then
            acc%bitset(limb_idx) = ior(acc%bitset(limb_idx), mask)
            acc%nz = acc%nz + 1
        end if
    end function idx_in_bitset

    subroutine sparsify_direct(acc, mat, row_num)
        implicit none
        type(dense_acc), target, intent(inout) :: acc
        type(psb_d_csr_sparse_mat), intent(inout) :: mat
        integer(psb_ipk_) :: row_num

        integer(psb_ipk_) :: i, j, info

        mat%irp(row_num + 1) = mat%irp(row_num) + acc%nz

        call psb_qsort(acc%nz_idxs(1:acc%nz))
        
        do i = 1, acc%nz
            j = acc%nz_idxs(i)
            mat%ja(mat%irp(row_num) + i - 1) = j
            mat%val(mat%irp(row_num) + i - 1) = acc%dense_val(j)
        end do
    end subroutine sparsify_direct
end module sp3mm_acc_mod