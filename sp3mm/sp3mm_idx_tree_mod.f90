module idx_tree_mod
    use psb_base_mod

    implicit none
    type :: idx_node
        integer(psb_ipk_) :: idx
        real(psb_dpk_) :: val
        type(idx_node), pointer :: left, right, parent
        logical :: is_red
    end type idx_node

    type :: idx_tree
        type(idx_node), pointer :: root
        integer(psb_ipk_) :: nnz

        contains

        procedure :: insert => insert_tree
    end type idx_tree

    contains

    subroutine insert_tree(this, idx, val)
        implicit none
        class(idx_tree), intent(inout) :: this
        integer(psb_ipk_), intent(in) :: idx
        real(psb_dpk_), intent(in) :: val

        type(idx_node), pointer :: new_node
        type(idx_node), pointer :: current, previous
        
        allocate(new_node)
        new_node%idx = idx
        new_node%val = val
        nullify(new_node%left)
        nullify(new_node%right)
        nullify(new_node%parent)
        new_node%is_red = .true.


        if (.not. associated(this%root)) then
            this%root => new_node
            this%nnz = 1
            new_node%is_red = .false.
            return
        end if

        current => this%root

        do while (associated(current))
            previous => current

            if (idx == current%idx) then
                current%val = current%val + val
                deallocate(new_node)
                return
            else if (idx < current%idx) then
                current => current%left
            else

                current => current%right
            end if
        end do

        if (idx < previous%idx) then
            new_node%parent => previous
            previous%left => new_node
        else
            new_node%parent => previous
            previous%right => new_node
        end if

        call fix_insertion(this, new_node)

        this%nnz = this%nnz + 1
    end subroutine insert_tree

    subroutine fix_insertion(this, node)
        class(idx_tree), intent(inout) :: this
        type(idx_node), pointer, intent(inout) :: node

        type(idx_node), pointer :: current, parent, grand_parent, uncle

        current => node
        parent => current%parent
        do while(associated(parent) .and. parent%is_red)
            ! grand parent exist because root can't be red
            grand_parent => parent%parent
            if (parent%idx < grand_parent%idx) then
                uncle => grand_parent%right
            else
                uncle => grand_parent%left
            end if

            if (associated(uncle) .and. uncle%is_red) then
                parent%is_red = .false.
                uncle%is_red = .false.
                grand_parent%is_red = .true.
                current => grand_parent
                parent => current%parent

            ! Left-Left case
            else if (current%idx < parent%idx .and. &
                parent%idx < grand_parent%idx) then
                call rotate_right(grand_parent)
                call swap_colors(parent, grand_parent)

                if (this%root%idx == grand_parent%idx) this%root => parent

                return
            ! Left-Right case
            else if (current%idx > parent%idx .and. &
                parent%idx < grand_parent%idx) then
                call rotate_left(parent)
                call rotate_right(grand_parent)
                call swap_colors(current, grand_parent)

                if (this%root%idx == grand_parent%idx) this%root => current

                return
            ! Right-Right case
            else if (current%idx > parent%idx .and. &
                parent%idx > grand_parent%idx) then
                call rotate_left(grand_parent)
                call swap_colors(parent, grand_parent)

                if (this%root%idx == grand_parent%idx) this%root => parent

                return
            ! Right-Left case
            else
                call rotate_right(parent)
                call rotate_left(grand_parent)
                call swap_colors(current, grand_parent)

                if (this%root%idx == grand_parent%idx) this%root => current
                
                return
            end if
        end do

        this%root%is_red = .false.
    end subroutine fix_insertion

    subroutine swap_colors(n1, n2)
        type(idx_node), pointer, intent(inout) :: n1, n2

        logical :: tmp

        tmp = n1%is_red
        n1%is_red = n2%is_red
        n2%is_red = tmp
    end subroutine swap_colors
    
    subroutine rotate_right(node)
        type(idx_node), pointer, intent(inout) :: node

        type(idx_node), pointer :: l, lr

        if (.not. associated(node%left)) return

        l => node%left
        lr => l%right
        node%left => lr

        if (associated(lr)) lr%parent => node

        if (associated(node%parent)) then
            if (node%idx < node%parent%idx) then
                node%parent%left => l
            else
                node%parent%right => l
            end if
        end if

        l%parent => node%parent
        node%parent => l

        l%right => node
    end subroutine rotate_right

    subroutine rotate_left(node)
        type(idx_node), pointer, intent(inout) :: node

        type(idx_node), pointer :: r, rl

        if (.not. associated(node%right)) return

        r => node%right
        rl => r%left
        node%right => rl

        if (associated(rl)) rl%parent => node

        if (associated(node%parent)) then
            if (node%idx < node%parent%idx) then
                node%parent%left => r
            else
                node%parent%right => r
            end if
        end if

        r%parent => node%parent
        node%parent => r

        r%left => node
    end subroutine rotate_left
    
    subroutine rb_tree_scalar_sparse_row_mul(tree, scalar, mat, row_num)
        implicit none
        type(idx_tree), intent(inout) :: tree
        real(psb_dpk_), intent(in) :: scalar
        type(psb_d_csr_sparse_mat), intent(in) :: mat
        integer(psb_ipk_), intent(in) :: row_num

        integer(psb_ipk_) :: i

        do i = mat%irp(row_num), mat%irp(row_num + 1) - 1
            call tree%insert(mat%ja(i),scalar * mat%val(i))
        end do

    end subroutine rb_tree_scalar_sparse_row_mul

    subroutine merge_trees(trees, mat)
        use omp_lib
        type(idx_tree), allocatable, intent(inout) :: trees(:)
        type(psb_d_csr_sparse_mat), intent(inout) :: mat

        integer(psb_ipk_) :: i, j, rows, info
        type(idx_node), pointer :: current, previous

        rows = size(trees)

        mat%irp(1) = 1
        
        do i=1, rows
            mat%irp(i + 1) = mat%irp(i) + trees(i)%nnz
        end do

        call psb_realloc(mat%irp(rows + 1), mat%val, info)
        call psb_realloc(mat%irp(rows + 1), mat%ja, info)

        do i = 1, size(trees)
            j = 0
            current => trees(i)%root
            do while(associated(current))
                ! go to the left-most node
                do while(associated(current%left))
                    current => current%left
                end do
                mat%val(j + mat%irp(i)) = current%val
                mat%ja(j + mat%irp(i)) = current%idx
                j = j + 1

                previous => current
                if (associated(current%right)) then
                    if (associated(current%parent)) then
                        current%parent%left => current%right
                    end if
                    current%right%parent => current%parent
                    current => current%right
                else
                    current => current%parent
                    if (associated(current)) nullify(current%left)
                end if
                deallocate(previous)
            end do
        end do
    end subroutine merge_trees

    subroutine merge_trees_distrib(trees, mat)
        use omp_lib
        type(idx_tree), allocatable, intent(inout) :: trees(:)
        type(psb_d_csr_sparse_mat), intent(inout) :: mat

        integer(psb_ipk_) :: i, j, rows, info
        type(idx_node), pointer :: current, previous

        rows = size(trees)

        mat%irp(1) = 1
        
        do i=1, rows
            mat%irp(i + 1) = mat%irp(i) + trees(i)%nnz
        end do

        call psb_realloc(mat%irp(rows + 1), mat%val, info)
        call psb_realloc(mat%irp(rows + 1), mat%ja, info)

        !$omp parallel do schedule(static), private(current, previous, j)
        do i = 1, size(trees)
            j = 0
            current => trees(i)%root
            do while(associated(current))
                ! go to the left-most node
                do while(associated(current%left))
                    current => current%left
                end do
                mat%val(j + mat%irp(i)) = current%val
                mat%ja(j + mat%irp(i)) = current%idx
                j = j + 1

                previous => current
                if (associated(current%right)) then
                    if (associated(current%parent)) then
                        current%parent%left => current%right
                    end if
                    current%right%parent => current%parent
                    current => current%right
                else
                    current => current%parent
                    if (associated(current)) nullify(current%left)
                end if
                deallocate(previous)
            end do
        end do
        !$omp end parallel do
    end subroutine merge_trees_distrib
end module idx_tree_mod