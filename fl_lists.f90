module fl_lists
implicit none

type fl_string256_list_item
	character(256), pointer :: value
	type(fl_string256_list_item), pointer :: next
end type fl_string256_list_item

type fl_string256_list
	type(fl_string256_list_item), pointer :: first, last, iter
	integer length
end type fl_string256_list


type fl_string128_list_item
	character(128), pointer :: value
	type(fl_string128_list_item), pointer :: next
end type fl_string128_list_item

type fl_string128_list
	type(fl_string128_list_item), pointer :: first, last, iter
	integer length
end type fl_string128_list


type fl_string80_list_item
	character(80), pointer :: value
	type(fl_string80_list_item), pointer :: next
end type fl_string80_list_item

type fl_string80_list
	type(fl_string80_list_item), pointer :: first, last, iter
	integer length
end type fl_string80_list


type fl_string64_list_item
	character(64), pointer :: value
	type(fl_string64_list_item), pointer :: next
end type fl_string64_list_item

type fl_string64_list
	type(fl_string64_list_item), pointer :: first, last, iter
	integer length
end type fl_string64_list


type fl_string32_list_item
	character(32), pointer :: value
	type(fl_string32_list_item), pointer :: next
end type fl_string32_list_item

type fl_string32_list
	type(fl_string32_list_item), pointer :: first, last, iter
	integer length
end type fl_string32_list


type fl_string16_list_item
	character(16), pointer :: value
	type(fl_string16_list_item), pointer :: next
end type fl_string16_list_item

type fl_string16_list
	type(fl_string16_list_item), pointer :: first, last, iter
	integer length
end type fl_string16_list


type fl_string8_list_item
	character(8), pointer :: value
	type(fl_string8_list_item), pointer :: next
end type fl_string8_list_item

type fl_string8_list
	type(fl_string8_list_item), pointer :: first, last, iter
	integer length
end type fl_string8_list


type fl_real_list_item
	real(4), pointer :: value
	type(fl_real_list_item), pointer :: next
end type fl_real_list_item

type fl_real_list
	type(fl_real_list_item), pointer :: first, last, iter
	integer length
end type fl_real_list


type fl_double_list_item
	real(8), pointer :: value
	type(fl_double_list_item), pointer :: next
end type fl_double_list_item

type fl_double_list
	type(fl_double_list_item), pointer :: first, last, iter
	integer length
end type fl_double_list


type fl_int_list_item
	integer(4), pointer :: value
	type(fl_int_list_item), pointer :: next
end type fl_int_list_item

type fl_int_list
	type(fl_int_list_item), pointer :: first, last, iter
	integer length
end type fl_int_list


type fl_long_list_item
	integer(8), pointer :: value
	type(fl_long_list_item), pointer :: next
end type fl_long_list_item

type fl_long_list
	type(fl_long_list_item), pointer :: first, last, iter
	integer length
end type fl_long_list

contains

subroutine fl_create_string256_list(list)
	type(fl_string256_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_string256_list

subroutine fl_append_string256_list(list,value)
	type(fl_string256_list) :: list
	character(256) :: value
	
	type(fl_string256_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_string256_iterator_list(list)
end subroutine fl_append_string256_list

subroutine fl_destroy_string256_list(list)
	type(fl_string256_list) :: list
	type(fl_string256_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_string256_list

function fl_string256_list_to_array(list) result(arr)
	type(fl_string256_list) :: list
	character(256), pointer, dimension(:) :: arr
	type(fl_string256_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_string256_list_to_array

function fl_member_string256_list(list,value) result(B)
	logical :: B
	type(fl_string256_list) :: list
	type(fl_string256_list_item), pointer :: item
	character(256) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_string256_list

function fl_pop_last_string256_list(list) result(pop)
	type(fl_string256_list) :: list
	type(fl_string256_list_item), pointer :: item	
	character(256) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_string256_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_string256_iterator_list(list)	
end function fl_pop_last_string256_list

function fl_pop_first_string256_list(list) result(pop)
	type(fl_string256_list) :: list
	type(fl_string256_list_item), pointer :: item	
	character(256) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_string256_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_string256_iterator_list(list)	
end function fl_pop_first_string256_list

function fl_get_string256_list(list,n) result(x)
	type(fl_string256_list) :: list
	integer n,i
	type(fl_string256_list_item), pointer :: item	
	character(256) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_string256_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_string256_list

function fl_iterate_string256_list(list) result(x)
	type(fl_string256_list) :: list
	character(256) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_string256_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_string256_list

subroutine fl_insert_string256_list(list,value,after)
	type(fl_string256_list) :: list
	character(256) :: value
	type(fl_string256_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_string256_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_string256_iterator_list(list)	
end subroutine fl_insert_string256_list

subroutine fl_reset_string256_iterator_list(list)
	type(fl_string256_list) :: list
	nullify(list%iter)
end subroutine fl_reset_string256_iterator_list




subroutine fl_create_string128_list(list)
	type(fl_string128_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_string128_list

subroutine fl_append_string128_list(list,value)
	type(fl_string128_list) :: list
	character(128) :: value
	
	type(fl_string128_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_string128_iterator_list(list)
end subroutine fl_append_string128_list

subroutine fl_destroy_string128_list(list)
	type(fl_string128_list) :: list
	type(fl_string128_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_string128_list

function fl_string128_list_to_array(list) result(arr)
	type(fl_string128_list) :: list
	character(128), pointer, dimension(:) :: arr
	type(fl_string128_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_string128_list_to_array

function fl_member_string128_list(list,value) result(B)
	logical :: B
	type(fl_string128_list) :: list
	type(fl_string128_list_item), pointer :: item
	character(128) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_string128_list

function fl_pop_last_string128_list(list) result(pop)
	type(fl_string128_list) :: list
	type(fl_string128_list_item), pointer :: item	
	character(128) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_string128_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_string128_iterator_list(list)	
end function fl_pop_last_string128_list

function fl_pop_first_string128_list(list) result(pop)
	type(fl_string128_list) :: list
	type(fl_string128_list_item), pointer :: item	
	character(128) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_string128_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_string128_iterator_list(list)	
end function fl_pop_first_string128_list

function fl_get_string128_list(list,n) result(x)
	type(fl_string128_list) :: list
	integer n,i
	type(fl_string128_list_item), pointer :: item	
	character(128) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_string128_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_string128_list

function fl_iterate_string128_list(list) result(x)
	type(fl_string128_list) :: list
	character(128) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_string128_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_string128_list

subroutine fl_insert_string128_list(list,value,after)
	type(fl_string128_list) :: list
	character(128) :: value
	type(fl_string128_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_string128_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_string128_iterator_list(list)	
end subroutine fl_insert_string128_list

subroutine fl_reset_string128_iterator_list(list)
	type(fl_string128_list) :: list
	nullify(list%iter)
end subroutine fl_reset_string128_iterator_list




subroutine fl_create_string80_list(list)
	type(fl_string80_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_string80_list

subroutine fl_append_string80_list(list,value)
	type(fl_string80_list) :: list
	character(80) :: value
	
	type(fl_string80_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_string80_iterator_list(list)
end subroutine fl_append_string80_list

subroutine fl_destroy_string80_list(list)
	type(fl_string80_list) :: list
	type(fl_string80_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_string80_list

function fl_string80_list_to_array(list) result(arr)
	type(fl_string80_list) :: list
	character(80), pointer, dimension(:) :: arr
	type(fl_string80_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_string80_list_to_array

function fl_member_string80_list(list,value) result(B)
	logical :: B
	type(fl_string80_list) :: list
	type(fl_string80_list_item), pointer :: item
	character(80) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_string80_list

function fl_pop_last_string80_list(list) result(pop)
	type(fl_string80_list) :: list
	type(fl_string80_list_item), pointer :: item	
	character(80) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_string80_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_string80_iterator_list(list)	
end function fl_pop_last_string80_list

function fl_pop_first_string80_list(list) result(pop)
	type(fl_string80_list) :: list
	type(fl_string80_list_item), pointer :: item	
	character(80) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_string80_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_string80_iterator_list(list)	
end function fl_pop_first_string80_list

function fl_get_string80_list(list,n) result(x)
	type(fl_string80_list) :: list
	integer n,i
	type(fl_string80_list_item), pointer :: item	
	character(80) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_string80_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_string80_list

function fl_iterate_string80_list(list) result(x)
	type(fl_string80_list) :: list
	character(80) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_string80_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_string80_list

subroutine fl_insert_string80_list(list,value,after)
	type(fl_string80_list) :: list
	character(80) :: value
	type(fl_string80_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_string80_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_string80_iterator_list(list)	
end subroutine fl_insert_string80_list

subroutine fl_reset_string80_iterator_list(list)
	type(fl_string80_list) :: list
	nullify(list%iter)
end subroutine fl_reset_string80_iterator_list




subroutine fl_create_string64_list(list)
	type(fl_string64_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_string64_list

subroutine fl_append_string64_list(list,value)
	type(fl_string64_list) :: list
	character(64) :: value
	
	type(fl_string64_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_string64_iterator_list(list)
end subroutine fl_append_string64_list

subroutine fl_destroy_string64_list(list)
	type(fl_string64_list) :: list
	type(fl_string64_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_string64_list

function fl_string64_list_to_array(list) result(arr)
	type(fl_string64_list) :: list
	character(64), pointer, dimension(:) :: arr
	type(fl_string64_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_string64_list_to_array

function fl_member_string64_list(list,value) result(B)
	logical :: B
	type(fl_string64_list) :: list
	type(fl_string64_list_item), pointer :: item
	character(64) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_string64_list

function fl_pop_last_string64_list(list) result(pop)
	type(fl_string64_list) :: list
	type(fl_string64_list_item), pointer :: item	
	character(64) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_string64_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_string64_iterator_list(list)	
end function fl_pop_last_string64_list

function fl_pop_first_string64_list(list) result(pop)
	type(fl_string64_list) :: list
	type(fl_string64_list_item), pointer :: item	
	character(64) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_string64_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_string64_iterator_list(list)	
end function fl_pop_first_string64_list

function fl_get_string64_list(list,n) result(x)
	type(fl_string64_list) :: list
	integer n,i
	type(fl_string64_list_item), pointer :: item	
	character(64) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_string64_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_string64_list

function fl_iterate_string64_list(list) result(x)
	type(fl_string64_list) :: list
	character(64) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_string64_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_string64_list

subroutine fl_insert_string64_list(list,value,after)
	type(fl_string64_list) :: list
	character(64) :: value
	type(fl_string64_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_string64_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_string64_iterator_list(list)	
end subroutine fl_insert_string64_list

subroutine fl_reset_string64_iterator_list(list)
	type(fl_string64_list) :: list
	nullify(list%iter)
end subroutine fl_reset_string64_iterator_list




subroutine fl_create_string32_list(list)
	type(fl_string32_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_string32_list

subroutine fl_append_string32_list(list,value)
	type(fl_string32_list) :: list
	character(32) :: value
	
	type(fl_string32_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_string32_iterator_list(list)
end subroutine fl_append_string32_list

subroutine fl_destroy_string32_list(list)
	type(fl_string32_list) :: list
	type(fl_string32_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_string32_list

function fl_string32_list_to_array(list) result(arr)
	type(fl_string32_list) :: list
	character(32), pointer, dimension(:) :: arr
	type(fl_string32_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_string32_list_to_array

function fl_member_string32_list(list,value) result(B)
	logical :: B
	type(fl_string32_list) :: list
	type(fl_string32_list_item), pointer :: item
	character(32) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_string32_list

function fl_pop_last_string32_list(list) result(pop)
	type(fl_string32_list) :: list
	type(fl_string32_list_item), pointer :: item	
	character(32) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_string32_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_string32_iterator_list(list)	
end function fl_pop_last_string32_list

function fl_pop_first_string32_list(list) result(pop)
	type(fl_string32_list) :: list
	type(fl_string32_list_item), pointer :: item	
	character(32) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_string32_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_string32_iterator_list(list)	
end function fl_pop_first_string32_list

function fl_get_string32_list(list,n) result(x)
	type(fl_string32_list) :: list
	integer n,i
	type(fl_string32_list_item), pointer :: item	
	character(32) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_string32_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_string32_list

function fl_iterate_string32_list(list) result(x)
	type(fl_string32_list) :: list
	character(32) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_string32_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_string32_list

subroutine fl_insert_string32_list(list,value,after)
	type(fl_string32_list) :: list
	character(32) :: value
	type(fl_string32_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_string32_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_string32_iterator_list(list)	
end subroutine fl_insert_string32_list

subroutine fl_reset_string32_iterator_list(list)
	type(fl_string32_list) :: list
	nullify(list%iter)
end subroutine fl_reset_string32_iterator_list




subroutine fl_create_string16_list(list)
	type(fl_string16_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_string16_list

subroutine fl_append_string16_list(list,value)
	type(fl_string16_list) :: list
	character(16) :: value
	
	type(fl_string16_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_string16_iterator_list(list)
end subroutine fl_append_string16_list

subroutine fl_destroy_string16_list(list)
	type(fl_string16_list) :: list
	type(fl_string16_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_string16_list

function fl_string16_list_to_array(list) result(arr)
	type(fl_string16_list) :: list
	character(16), pointer, dimension(:) :: arr
	type(fl_string16_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_string16_list_to_array

function fl_member_string16_list(list,value) result(B)
	logical :: B
	type(fl_string16_list) :: list
	type(fl_string16_list_item), pointer :: item
	character(16) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_string16_list

function fl_pop_last_string16_list(list) result(pop)
	type(fl_string16_list) :: list
	type(fl_string16_list_item), pointer :: item	
	character(16) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_string16_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_string16_iterator_list(list)	
end function fl_pop_last_string16_list

function fl_pop_first_string16_list(list) result(pop)
	type(fl_string16_list) :: list
	type(fl_string16_list_item), pointer :: item	
	character(16) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_string16_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_string16_iterator_list(list)	
end function fl_pop_first_string16_list

function fl_get_string16_list(list,n) result(x)
	type(fl_string16_list) :: list
	integer n,i
	type(fl_string16_list_item), pointer :: item	
	character(16) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_string16_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_string16_list

function fl_iterate_string16_list(list) result(x)
	type(fl_string16_list) :: list
	character(16) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_string16_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_string16_list

subroutine fl_insert_string16_list(list,value,after)
	type(fl_string16_list) :: list
	character(16) :: value
	type(fl_string16_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_string16_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_string16_iterator_list(list)	
end subroutine fl_insert_string16_list

subroutine fl_reset_string16_iterator_list(list)
	type(fl_string16_list) :: list
	nullify(list%iter)
end subroutine fl_reset_string16_iterator_list




subroutine fl_create_string8_list(list)
	type(fl_string8_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_string8_list

subroutine fl_append_string8_list(list,value)
	type(fl_string8_list) :: list
	character(8) :: value
	
	type(fl_string8_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_string8_iterator_list(list)
end subroutine fl_append_string8_list

subroutine fl_destroy_string8_list(list)
	type(fl_string8_list) :: list
	type(fl_string8_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_string8_list

function fl_string8_list_to_array(list) result(arr)
	type(fl_string8_list) :: list
	character(8), pointer, dimension(:) :: arr
	type(fl_string8_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_string8_list_to_array

function fl_member_string8_list(list,value) result(B)
	logical :: B
	type(fl_string8_list) :: list
	type(fl_string8_list_item), pointer :: item
	character(8) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_string8_list

function fl_pop_last_string8_list(list) result(pop)
	type(fl_string8_list) :: list
	type(fl_string8_list_item), pointer :: item	
	character(8) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_string8_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_string8_iterator_list(list)	
end function fl_pop_last_string8_list

function fl_pop_first_string8_list(list) result(pop)
	type(fl_string8_list) :: list
	type(fl_string8_list_item), pointer :: item	
	character(8) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_string8_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_string8_iterator_list(list)	
end function fl_pop_first_string8_list

function fl_get_string8_list(list,n) result(x)
	type(fl_string8_list) :: list
	integer n,i
	type(fl_string8_list_item), pointer :: item	
	character(8) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_string8_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_string8_list

function fl_iterate_string8_list(list) result(x)
	type(fl_string8_list) :: list
	character(8) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_string8_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_string8_list

subroutine fl_insert_string8_list(list,value,after)
	type(fl_string8_list) :: list
	character(8) :: value
	type(fl_string8_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_string8_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_string8_iterator_list(list)	
end subroutine fl_insert_string8_list

subroutine fl_reset_string8_iterator_list(list)
	type(fl_string8_list) :: list
	nullify(list%iter)
end subroutine fl_reset_string8_iterator_list




subroutine fl_create_real_list(list)
	type(fl_real_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_real_list

subroutine fl_append_real_list(list,value)
	type(fl_real_list) :: list
	real(4) :: value
	
	type(fl_real_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_real_iterator_list(list)
end subroutine fl_append_real_list

subroutine fl_destroy_real_list(list)
	type(fl_real_list) :: list
	type(fl_real_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_real_list

function fl_real_list_to_array(list) result(arr)
	type(fl_real_list) :: list
	real(4), pointer, dimension(:) :: arr
	type(fl_real_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_real_list_to_array

function fl_member_real_list(list,value) result(B)
	logical :: B
	type(fl_real_list) :: list
	type(fl_real_list_item), pointer :: item
	real(4) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_real_list

function fl_pop_last_real_list(list) result(pop)
	type(fl_real_list) :: list
	type(fl_real_list_item), pointer :: item	
	real(4) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_real_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_real_iterator_list(list)	
end function fl_pop_last_real_list

function fl_pop_first_real_list(list) result(pop)
	type(fl_real_list) :: list
	type(fl_real_list_item), pointer :: item	
	real(4) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_real_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_real_iterator_list(list)	
end function fl_pop_first_real_list

function fl_get_real_list(list,n) result(x)
	type(fl_real_list) :: list
	integer n,i
	type(fl_real_list_item), pointer :: item	
	real(4) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_real_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_real_list

function fl_iterate_real_list(list) result(x)
	type(fl_real_list) :: list
	real(4) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_real_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_real_list

subroutine fl_insert_real_list(list,value,after)
	type(fl_real_list) :: list
	real(4) :: value
	type(fl_real_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_real_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_real_iterator_list(list)	
end subroutine fl_insert_real_list

subroutine fl_reset_real_iterator_list(list)
	type(fl_real_list) :: list
	nullify(list%iter)
end subroutine fl_reset_real_iterator_list




subroutine fl_create_double_list(list)
	type(fl_double_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_double_list

subroutine fl_append_double_list(list,value)
	type(fl_double_list) :: list
	real(8) :: value
	
	type(fl_double_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_double_iterator_list(list)
end subroutine fl_append_double_list

subroutine fl_destroy_double_list(list)
	type(fl_double_list) :: list
	type(fl_double_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_double_list

function fl_double_list_to_array(list) result(arr)
	type(fl_double_list) :: list
	real(8), pointer, dimension(:) :: arr
	type(fl_double_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_double_list_to_array

function fl_member_double_list(list,value) result(B)
	logical :: B
	type(fl_double_list) :: list
	type(fl_double_list_item), pointer :: item
	real(8) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_double_list

function fl_pop_last_double_list(list) result(pop)
	type(fl_double_list) :: list
	type(fl_double_list_item), pointer :: item	
	real(8) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_double_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_double_iterator_list(list)	
end function fl_pop_last_double_list

function fl_pop_first_double_list(list) result(pop)
	type(fl_double_list) :: list
	type(fl_double_list_item), pointer :: item	
	real(8) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_double_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_double_iterator_list(list)	
end function fl_pop_first_double_list

function fl_get_double_list(list,n) result(x)
	type(fl_double_list) :: list
	integer n,i
	type(fl_double_list_item), pointer :: item	
	real(8) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_double_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_double_list

function fl_iterate_double_list(list) result(x)
	type(fl_double_list) :: list
	real(8) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_double_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_double_list

subroutine fl_insert_double_list(list,value,after)
	type(fl_double_list) :: list
	real(8) :: value
	type(fl_double_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_double_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_double_iterator_list(list)	
end subroutine fl_insert_double_list

subroutine fl_reset_double_iterator_list(list)
	type(fl_double_list) :: list
	nullify(list%iter)
end subroutine fl_reset_double_iterator_list




subroutine fl_create_int_list(list)
	type(fl_int_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_int_list

subroutine fl_append_int_list(list,value)
	type(fl_int_list) :: list
	integer(4) :: value
	
	type(fl_int_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_int_iterator_list(list)
end subroutine fl_append_int_list

subroutine fl_destroy_int_list(list)
	type(fl_int_list) :: list
	type(fl_int_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_int_list

function fl_int_list_to_array(list) result(arr)
	type(fl_int_list) :: list
	integer(4), pointer, dimension(:) :: arr
	type(fl_int_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_int_list_to_array

function fl_member_int_list(list,value) result(B)
	logical :: B
	type(fl_int_list) :: list
	type(fl_int_list_item), pointer :: item
	integer(4) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_int_list

function fl_pop_last_int_list(list) result(pop)
	type(fl_int_list) :: list
	type(fl_int_list_item), pointer :: item	
	integer(4) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_int_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_int_iterator_list(list)	
end function fl_pop_last_int_list

function fl_pop_first_int_list(list) result(pop)
	type(fl_int_list) :: list
	type(fl_int_list_item), pointer :: item	
	integer(4) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_int_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_int_iterator_list(list)	
end function fl_pop_first_int_list

function fl_get_int_list(list,n) result(x)
	type(fl_int_list) :: list
	integer n,i
	type(fl_int_list_item), pointer :: item	
	integer(4) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_int_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_int_list

function fl_iterate_int_list(list) result(x)
	type(fl_int_list) :: list
	integer(4) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_int_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_int_list

subroutine fl_insert_int_list(list,value,after)
	type(fl_int_list) :: list
	integer(4) :: value
	type(fl_int_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_int_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_int_iterator_list(list)	
end subroutine fl_insert_int_list

subroutine fl_reset_int_iterator_list(list)
	type(fl_int_list) :: list
	nullify(list%iter)
end subroutine fl_reset_int_iterator_list




subroutine fl_create_long_list(list)
	type(fl_long_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_long_list

subroutine fl_append_long_list(list,value)
	type(fl_long_list) :: list
	integer(8) :: value
	
	type(fl_long_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_long_iterator_list(list)
end subroutine fl_append_long_list

subroutine fl_destroy_long_list(list)
	type(fl_long_list) :: list
	type(fl_long_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_long_list

function fl_long_list_to_array(list) result(arr)
	type(fl_long_list) :: list
	integer(8), pointer, dimension(:) :: arr
	type(fl_long_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_long_list_to_array

function fl_member_long_list(list,value) result(B)
	logical :: B
	type(fl_long_list) :: list
	type(fl_long_list_item), pointer :: item
	integer(8) :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_long_list

function fl_pop_last_long_list(list) result(pop)
	type(fl_long_list) :: list
	type(fl_long_list_item), pointer :: item	
	integer(8) :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_long_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_long_iterator_list(list)	
end function fl_pop_last_long_list

function fl_pop_first_long_list(list) result(pop)
	type(fl_long_list) :: list
	type(fl_long_list_item), pointer :: item	
	integer(8) :: pop
	if (list%length==0) stop 'Requested first item from empty fl_long_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_long_iterator_list(list)	
end function fl_pop_first_long_list

function fl_get_long_list(list,n) result(x)
	type(fl_long_list) :: list
	integer n,i
	type(fl_long_list_item), pointer :: item	
	integer(8) :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_long_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_long_list

function fl_iterate_long_list(list) result(x)
	type(fl_long_list) :: list
	integer(8) :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_long_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_long_list

subroutine fl_insert_long_list(list,value,after)
	type(fl_long_list) :: list
	integer(8) :: value
	type(fl_long_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_long_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_long_iterator_list(list)	
end subroutine fl_insert_long_list

subroutine fl_reset_long_iterator_list(list)
	type(fl_long_list) :: list
	nullify(list%iter)
end subroutine fl_reset_long_iterator_list



end module fl_lists

