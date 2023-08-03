module sp3mm_acc_mod
    use psb_base_mod
    use sp3mm_const_mod
    implicit none

    type :: sparse_acc
        
        real(psb_dpk_), allocatable :: as(:)
        integer(psb_ipk_), allocatable :: ja(:)
        integer(psb_ipk_), allocatable :: nnz_idxs(:)
        integer(psb_ipk_) :: nnz

        contains

        procedure :: init => init_sparse_acc
    end type sparse_acc

    contains

    subroutine init_sparse_acc(sp_acc, size)
        class(sparse_acc), intent(out) :: sp_acc
        integer(psb_ipk_), intent(in) :: size
        integer(psb_ipk_) :: info

        call psb_realloc(size, sp_acc%as, info)
        call psb_realloc(size, sp_acc%ja, info)
        call psb_realloc(size, sp_acc%nnz_idxs, info)

    end subroutine init_sparse_acc

    subroutine upper_bound_scalar_sparse_mul_row(sp_acc, scalar, mat, row_num)
        type(sparse_acc), intent(inout) :: sp_acc
        real(psb_dpk_), intent(in) :: scalar
        type(psb_d_csr_sparse_mat), intent(in) :: mat
        integer(psb_ipk_), intent(in) :: row_num

        integer(psb_ipk_) :: i, k
        sp_acc%nnz = 0

        do i = mat%irp(row_num), mat%irp(row_num + 1) - 1
            do k = 1, sp_acc%nnz
                if (sp_acc%ja(k) == mat%ja(i)) then
                    sp_acc%as(k) = sp_acc%as(k) + scalar * mat%val(i)
                    exit
                end if
            end do
            sp_acc%nnz = sp_acc%nnz + 1
            sp_acc%ja(sp_acc%nnz) = mat%ja(i)
            sp_acc%as(sp_acc%nnz) = scalar * mat%val(i)
        end do
    end subroutine upper_bound_scalar_sparse_mul_row

    subroutine merge_rows(sp_accs, mat)
        use psb_util_mod
        type(sparse_acc), allocatable, intent(inout) :: sp_accs(:)
        type(psb_d_csr_sparse_mat), intent(inout) :: mat

        integer(psb_ipk_) :: nnz, row, col
        integer(psb_ipk_) :: info

        nnz = 0
        mat%irp(1) = 1

        do row = 1, mat%get_nrows()
            nnz = nnz + sp_accs(row)%nnz
            mat%irp(row + 1) = nnz + 1
        end do

        call psb_realloc(nnz, mat%ja, info)
        call psb_realloc(nnz, mat%val, info)

        !omp parallel do schedule(runtime)
        print *,  mat%get_nrows()
        do row = 1, mat%get_nrows()
            call psb_msort(sp_accs(row)%ja(1:sp_accs(row)%nnz), sp_accs(row)%nnz_idxs(1:sp_accs(row)%nnz))
            do col = 1, sp_accs(row)%nnz
                ! print *, sp_accs(row)%nnz_idxs(1:sp_accs(row)%nnz)
                mat%val(mat%irp(row) + col - 1) = sp_accs(row)%as(sp_accs(row)%nnz_idxs(col))
                mat%ja(mat%irp(row) + col - 1) = sp_accs(row)%ja(col)
            end do
        end do
        !omp end parallel do

    end subroutine merge_rows
end module sp3mm_acc_mod