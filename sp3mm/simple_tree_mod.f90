module simple_tree_mod
    use psb_base_mod
    implicit none

    type :: simple_node
        integer(psb_ipk_) :: idx
        real(psb_dpk_) :: val
        type(simple_node), pointer :: left, right
    end type simple_node

    type :: simple_tree
        integer(psb_ipk_) :: nnz
        type(simple_node), pointer :: root
    end type simple_tree

    contains
        subroutine insert_in_tree(this, idx, val)
            class(simple_tree), intent(inout) :: this
            integer(psb_ipk_), intent(in) :: idx
            real(psb_dpk_), intent(in) :: val

            type(simple_node), pointer :: new_node, current, previous

            if (.not. associated(this%root)) then
                allocate(new_node)
                new_node%idx = idx
                new_node%val = val
                nullify(new_node%left)
                nullify(new_node%right)
                this%root => new_node
                this%nnz = 1
                return
            end if

            current => this%root

            do while (associated(current))
                if (idx == current%idx) then
                    current%val = current%val + val
                    return
                else if (idx < current%idx) then
                    previous => current
                    current => current%left
                else
                    previous => current
                    current => current%right
                end if
            end do

            allocate(new_node)
            new_node%idx = idx
            new_node%val = val
            nullify(new_node%left)
            nullify(new_node%right)

            if (idx < previous%idx) then
                previous%left => new_node
            else
                previous%right => new_node
            end if

            this%nnz = this%nnz + 1
        end subroutine insert_in_tree

        subroutine merge_tree_to_arrays(tree, idx_array, val_array)
            type(simple_tree), intent(in) :: tree
            integer(psb_ipk_), intent(out) :: idx_array(:)
            real(psb_dpk_), intent(out) :: val_array(:)
            integer(psb_ipk_) :: array_index, stack_size
            type(simple_node), allocatable, target :: stack(:)
            type(simple_node), pointer :: current_node
        
            allocate(stack(tree%nnz))
        
            array_index = 1
            stack_size = 0
            current_node => tree%root
        
            ! Traverse the tree in-order using an iterative approach
            do while (associated(current_node) .or. stack_size > 0)
                if (associated(current_node)) then
                    ! Push nodes onto the stack and move to the left child
                    stack_size = stack_size + 1
                    stack(stack_size) = current_node
                    current_node => current_node%left
                else if (stack_size > 0) then
                    ! Pop a node from the stack, process it, and move to the right child
                    current_node => stack(stack_size)
                    stack_size = stack_size - 1
                    idx_array(array_index) = current_node%idx
                    val_array(array_index) = current_node%val
                    array_index = array_index + 1
                    current_node => current_node%right
                end if
            end do
        
            deallocate(stack)
        end subroutine merge_tree_to_arrays

end module simple_tree_mod