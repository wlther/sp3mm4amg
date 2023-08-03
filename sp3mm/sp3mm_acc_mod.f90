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
        sp_acc%nnz = 0
        sp_acc%ja = 0

    end subroutine init_sparse_acc

    subroutine upper_bound_scalar_sparse_mul_row(sp_acc, scalar, mat, row_num)
        type(sparse_acc), intent(inout) :: sp_acc
        real(psb_dpk_), intent(in) :: scalar
        type(psb_d_csr_sparse_mat), intent(in) :: mat
        integer(psb_ipk_), intent(in) :: row_num

        integer(psb_ipk_) :: i, k

        do i = mat%irp(row_num), mat%irp(row_num + 1) - 1
            k = 1
            do while (k <= sp_acc%nnz)
                if (sp_acc%ja(k) == mat%ja(i)) then
                    sp_acc%as(k) = sp_acc%as(k) + scalar * mat%val(i)
                    exit
                end if
                k = k + 1
            end do
            if (k > sp_acc%nnz) then 
                sp_acc%nnz = sp_acc%nnz + 1
                sp_acc%ja(sp_acc%nnz) = mat%ja(i)
                sp_acc%as(sp_acc%nnz) = scalar * mat%val(i)
            end if
        end do
    end subroutine upper_bound_scalar_sparse_mul_row

    subroutine merge_rows(sp_accs, mat)
        type(sparse_acc), allocatable, intent(inout) :: sp_accs(:)
        type(psb_d_csr_sparse_mat), intent(inout) :: mat

        integer(psb_ipk_) :: nnz, i, j
        integer(psb_ipk_) :: info

        nnz = 0
        mat%irp(1) = 1

        do i = 1, mat%get_nrows()
            nnz = nnz + sp_accs(i)%nnz
            mat%irp(i + 1) = nnz + 1
        end do

        call psb_realloc(nnz, mat%val, info)
        call psb_realloc(nnz, mat%ja, info)

        !$omp parallel do schedule(runtime)
        do i = 1, mat%get_nrows()
            do j = 1, sp_accs(i)%nnz
                mat%val(mat%irp(i) + j - 1) = sp_accs(i)%as(sp_accs(i)%nnz_idxs(j))
                mat%ja(mat%irp(i) + j - 1) = sp_accs(i)%ja(j)
            end do
        end do
        !$omp end parallel do

    end subroutine merge_rows
end module sp3mm_acc_mod