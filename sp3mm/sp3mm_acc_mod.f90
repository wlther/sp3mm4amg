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
        integer(psb_ipk_) :: nz
        ! bitmap
        integer(psb_ipk_), allocatable :: bitset(:)

        contains

        procedure :: init   => init_dense_acc
        procedure :: reset  => reset_dense_acc
    end type dense_acc

    type :: sparse_acc
        
        real(psb_dpk_), allocatable :: as(:)
        integer(psb_ipk_), allocatable :: ja(:)
        integer(psb_ipk_) :: nz

        contains

        procedure :: init => init_sparse_acc
    end type sparse_acc

    contains

    subroutine init_dense_acc(acc, size)
        implicit none
        class(dense_acc), intent(inout)  :: acc
        integer(psb_ipk_), intent(in)   :: size

        ! call psb_realloc(size, sp_acc%dense_val, info)
        allocate(acc%dense_val(size))
        allocate(acc%nz_idxs(size))
        allocate(acc%bitset(size))

        acc%dense_val = 0
        acc%bitset = 0
        acc%nz = 0
    end subroutine init_dense_acc

    subroutine reset_dense_acc(acc)
        implicit none
        class(dense_acc), intent(inout) :: acc

        acc%dense_val = 0
        acc%bitset = 0
        acc%nz = 0
    end subroutine reset_dense_acc

    subroutine init_sparse_acc(sp_acc)
        implicit none
        class(sparse_acc), intent(out) :: sp_acc

        sp_acc%nz = 0
    end subroutine init_sparse_acc

    subroutine scalar_sparse_row_mul(acc, scalar, mat, row_num)
        implicit none
        type(dense_acc), intent(inout)  :: acc
        real(psb_dpk_) :: scalar
        type(psb_d_csr_sparse_mat), intent(in) :: mat
        integer(psb_ipk_) :: row_num

        integer(psb_ipk_) :: i, j

        do i = mat%irp(row_num), mat%irp(row_num + 1) - 1
            j = mat%ja(i)

            acc%dense_val(j) = acc%dense_val(j) + scalar * mat%val(i)

            if (.not. idx_in_bitset(j, acc)) then
                acc%nz_idxs(acc%nz) = j
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

        integer(psb_ipk_) :: i, j

        mat%irp(row_num + 1) = mat%irp(row_num) + acc%nz

        call psb_qsort(acc%nz_idxs(1:acc%nz))
        
        do i = 1, acc%nz
            j = acc%nz_idxs(i)
            mat%ja(mat%irp(row_num) + i - 1) = j
            mat%val(mat%irp(row_num) + i - 1) = acc%dense_val(j)
        end do
    end subroutine sparsify_direct

    subroutine sparsify_ub(acc, sp_acc, start_col_acc)
        implicit none
        type(dense_acc), intent(inout) :: acc
        type(sparse_acc), intent(inout) :: sp_acc
        integer(psb_ipk_), intent(in) :: start_col_acc
        
        integer(psb_ipk_) :: i, j, info

        call psb_qsort(acc%nz_idxs(1:acc%nz))
        
        sp_acc%nz = acc%nz

        if (sp_acc%nz > 0) then
            call psb_realloc(sp_acc%nz, sp_acc%ja, info)
            call psb_realloc(sp_acc%nz, sp_acc%as, info)
        end if

        do i = 1, acc%nz
            j = acc%nz_idxs(i)
            sp_acc%ja(i) = j + start_col_acc
            sp_acc%as(i) = acc%dense_val(j)
        end do
    end subroutine sparsify_ub

    subroutine merge_rows(sp_accs, mat)
        type(sparse_acc), allocatable, intent(in) :: sp_accs(:)
        type(psb_d_csr_sparse_mat), intent(inout) :: mat

        integer(psb_ipk_) :: nnz, row

        nnz = 0
        mat%irp(1) = 1

        do row = 1, mat%get_nrows()
            nnz = nnz + sp_accs(row)%nz
            mat%irp(row + 1) = nnz + 1
        end do
        !$omp parallel do schedule(runtime)
        do row = 1, mat%get_nrows()
            mat%val(mat%irp(row):mat%irp(row+1)-1) = sp_accs(row)%as(1:sp_accs(row)%nz)
            mat%ja(mat%irp(row):mat%irp(row+1)-1) = sp_accs(row)%ja(1:sp_accs(row)%nz)
        end do
        !$omp end parallel do

    end subroutine merge_rows
end module sp3mm_acc_mod