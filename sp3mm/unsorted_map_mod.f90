module unordered_map_mod
    use psb_base_mod
    implicit none
    
    type, public :: unordered_map
    private
    integer(psb_ipk_) :: MAX_BUCKETS = 1000 ! Adjust this as needed
    type(bucket), allocatable :: buckets(:)
    integer(psb_ipk_) :: size
    contains 
        procedure :: create_map => unordered_map_create_map
        procedure :: insert => unordered_map_insert
        procedure :: insert_add => unordered_map_insert_add
        procedure :: find => unordered_map_find
        procedure :: get_size => unordered_map_get_size
        procedure :: clear_map => unordered_map_clear_map
        procedure :: move_to_arrays => unordered_map_move_to_arrays
        procedure :: destroy_map => unordered_map_destroy_map 
    end type unordered_map
  
    type :: key_val
      integer(psb_ipk_) :: key
      real(psb_dpk_) :: val
    end type key_val
  
    type :: bucket
      type(key_val), pointer :: entries(:)
    end type bucket
  
  contains
  
    subroutine unordered_map_create_map(map, max_buckets)
      class(unordered_map), intent(out) :: map
      integer(psb_ipk_), intent(in), optional :: max_buckets
      type(key_val), allocatable, target :: new_entry(:)
      integer(psb_ipk_) :: i

      map%size = 0
      
      if (present(max_buckets)) map%MAX_BUCKETS = max_buckets
      
      allocate(map%buckets(map%MAX_BUCKETS))

      do i = 1, map%MAX_BUCKETS
        allocate(new_entry(0))
        map%buckets(i)%entries => new_entry
      end do
    end subroutine unordered_map_create_map
  
    subroutine unordered_map_insert(map, key, val)
      class(unordered_map), intent(inout) :: map
      integer(psb_ipk_), intent(in) :: key
      real(psb_dpk_), intent(in) :: val
      type(key_val) :: new_entry
  
      integer(psb_ipk_) :: hash, i
  
      hash = modulo(key, map%MAX_BUCKETS) + 1
  
      do i = 1, size(map%buckets(hash)%entries)
        if (map%buckets(hash)%entries(i)%key == key) then
          map%buckets(hash)%entries(i)%val = val
          return
        end if
      end do
  
      new_entry%key = key
      new_entry%val = val
      map%buckets(hash)%entries = [map%buckets(hash)%entries, new_entry]
      map%size = map%size + 1
    end subroutine unordered_map_insert

    subroutine unordered_map_insert_add(map, key, val)
        class(unordered_map), intent(inout) :: map
        integer(psb_ipk_), intent(in) :: key
        real(psb_dpk_), intent(in) :: val
        type(key_val) :: new_entry
    
        integer(psb_ipk_) :: hash, i
    
        hash = modulo(key, map%MAX_BUCKETS) + 1
    
        do i = 1, size(map%buckets(hash)%entries)
          print *, map%buckets(hash)%entries(i)%key
          if (map%buckets(hash)%entries(i)%key == key) then
            map%buckets(hash)%entries(i)%val = map%buckets(hash)%entries(i)%val + val
            return
          end if
        end do
        
        new_entry%key = key
        new_entry%val = val
        map%buckets(hash)%entries = [map%buckets(hash)%entries, new_entry]
        map%size = map%size + 1
      end subroutine unordered_map_insert_add
  
    subroutine unordered_map_find(map, key, val, found)
      class(unordered_map), intent(in) :: map
      integer(psb_ipk_), intent(in) :: key
      real(psb_dpk_), intent(out) :: val
      logical, intent(out) :: found
  
      integer(psb_ipk_) :: hash, i
  
      hash = modulo(key, map%MAX_BUCKETS) + 1
  
      found = .false.
  
      do i = 1, size(map%buckets(hash)%entries)
        if (map%buckets(hash)%entries(i)%key == key) then
          val = map%buckets(hash)%entries(i)%val
          found = .true.
          return
        end if
      end do
    end subroutine unordered_map_find
  
    function unordered_map_get_size(map) result(n)
      class(unordered_map), intent(in) :: map
      integer(psb_ipk_) :: n

      n = map%size
    end function unordered_map_get_size

    subroutine unordered_map_clear_map(map)
      class(unordered_map), intent(inout) :: map
      integer(psb_ipk_) :: i
  
      do i = 1, map%MAX_BUCKETS
        deallocate(map%buckets(i)%entries)
      end do
  
      map%size = 0
    end subroutine unordered_map_clear_map
  
    subroutine unordered_map_move_to_arrays(map, idxs, vals)
      class(unordered_map), intent(inout) :: map
      integer(psb_ipk_), intent(inout) :: idxs(:)
      real(psb_dpk_), intent(inout) :: vals(:)

      integer(psb_ipk_), allocatable :: ix(:)
      integer(psb_ipk_) :: i, j, k

      k = 1
      do i = 1, map%MAX_BUCKETS
        do j = 1, size(map%buckets(i)%entries)
          idxs(k) = map%buckets(i)%entries(j)%key
          vals(k) = map%buckets(i)%entries(j)%val
          k = k + 1
        end do
      end do

      call psb_qsort(idxs, ix)

      vals(ix) = vals

      call map%clear_map()
    end subroutine unordered_map_move_to_arrays

    subroutine unordered_map_destroy_map(map)
      class(unordered_map), intent(inout) :: map
      call map%clear_map()
      deallocate(map%buckets)
    end subroutine unordered_map_destroy_map
  
  end module unordered_map_mod
  