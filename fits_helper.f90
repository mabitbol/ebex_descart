! A few FITSIO helper routines.

MODULE DS_FITSTOOLS

USE healpix_types
USE DS_TYPES
USE DS_UTILS

IMPLICIT NONE

PRIVATE
PUBLIC ::  ds_write_fits_map, ds_write_fits_covariance

TYPE fits_header
	INTEGER :: nrecords
	INTEGER,ALLOCATABLE,DIMENSION(:) :: record_int
	REAL(DP),ALLOCATABLE,DIMENSION(:) :: record_dp
	CHARACTER(len=80),ALLOCATABLE,DIMENSION(:) :: name, record_str 
	CHARACTER(len=1),ALLOCATABLE,DIMENSION(:) :: htype !see FITSIO guide (p.17) for conventions: J=integer, D=real(8), S=string
END TYPE fits_header


CONTAINS

! Cookbook routine to delete file
subroutine deletefile(filename,status)
    integer(4), intent(out) :: status
    character(len=80), intent(in) :: filename
    integer(4) :: unit,blocksize

    unit = ds_get_lun()

    open(file=filename, unit=unit, form='unformatted')
    close(unit=unit, status='delete')

    ! ! Simply return if status is greater than zero
    ! if (status .gt. 0)return
    ! 
    ! ! Get an unused Logical Unit Number to use to open the FITS file
    ! call ftgiou(unit,status)
    ! 
    ! ! Try to open the file, to see if it exists
    ! call ftopen(unit,filename,1,blocksize,status)
    ! 
    ! if (status .eq. 0)then!
    ! ! file was opened;  so now delete it 
    !    call ftdelt(unit,status)
    ! else if (status .eq. 103)then
    ! ! file doesn't exist, so just reset status to zero and clear errors
    !    status=0
    !    call ftcmsg
    ! else
    ! ! there was some other error opening the file; delete the file anyway
    !    status=0
    !    call ftcmsg
    !    call ftdelt(unit,status)
    ! end if
    ! 
    ! ! Free the unit number for later reuse
    ! call ftfiou(unit, status)

end subroutine


subroutine ds_write_fits_map(filename,maps,nside,units,isRing,coordinate_system, hitCount, covariance, extra_headers)
    use fitstools, only : output_map
    use head_fits, only : write_minimal_header, add_card
    type(ds_trimap) :: maps
    type(ds_covariance), optional :: covariance
    integer(8), dimension(1:), optional :: hitCount
    character(len=80), dimension(1:), optional :: extra_headers

    integer nside
    real(dp), dimension(:), pointer :: dummy_map
    character(*) :: filename
    integer card_pos
    logical write_card
    integer status, i
    character(1), optional :: coordinate_system
    logical, optional :: isRing
    integer(4) order_code
    integer(4), dimension(:), pointer :: indices
    character(8), optional :: units
    character(8) :: unitsCode
    integer :: nheal
    !integer,parameter :: map_type=2
    real(sp),allocatable,dimension(:,:) :: map
    character(len=80),dimension(1:600) :: header
    real(dp), parameter :: hpx_dbadval = -1.6375e30_dp
    character(len=20) :: dtype, coordsys
    character(len=8) :: ordering
    integer npix ,nmap
    integer t
    call ds_assert(maps%has_T .or. maps%has_P, "Maps with neither T nor P passed to write_fits_maps")

    if (maps%has_t) then
            npix = maps%T%npix
    else
            npix = maps%Q%npix	
    endif

    if (.not. (maps%has_T .and. maps%has_P)) then
            allocate(dummy_map(npix))
            dummy_map = 0
    endif

    if (present(coordinate_system)) then
            coordsys = coordinate_system
    else
            coordsys = 'C'
    endif



    status=0
    call deletefile(filename,status)

    if (.not. present(isRing)) then
       order_code = 2
            ordering='NESTED'
    elseif (isRing) then
       order_code = 1
    ordering='RING'
    else
            ordering='NESTED'
       order_code=2
    endif

    if (.not. present(units)) then
       unitsCode = 'unknown '
    else
       unitsCode= units
    endif

    nmap = 3
    if (present(covariance)) nmap = nmap + 6
    if (present(hitCount)) nmap = nmap + 1
    nheal= nside**2 *12

    allocate(map(0:nheal-1,1:nmap))

    map= hpx_dbadval
    if (maps%has_T) map(maps%T%indices,1)= maps%T%map
    if (maps%has_P) map(maps%Q%indices,2)= maps%Q%map
    if (maps%has_P) map(maps%U%indices,3)= maps%U%map
    if (present(covariance)) then
            if (covariance%has_t) map(maps%T%indices,4)= covariance%TT
            if (covariance%has_p) map(maps%Q%indices,5)= covariance%QQ
            if (covariance%has_p) map(maps%U%indices,6)= covariance%UU
            if (covariance%has_p .and. covariance%has_t) map(maps%U%indices,7)= covariance%TQ
            if (covariance%has_p .and. covariance%has_t) map(maps%U%indices,8)= covariance%TU
            if (covariance%has_p) map(maps%Q%indices,9)= covariance%QU
    endif
    if (present(hitCount) .and. present(covariance)) then
            map(:,10)= hitCount(1:nheal)
    elseif (present(hitCount))	then
            map(:,4)= hitCount(1:nheal)
    endif

    dtype= 'MAP'

    call write_minimal_header( header, dtype, nside=nside, &
            ordering=ordering, coordsys=coordsys,polar=.true.,units=unitsCode)

    if (present(extra_headers)) then
            
            card_pos = size(header)+1
            write_card=.true.
            do
                    card_pos=card_pos-1
                    if (header(card_pos) .ne. '') exit
                    if (card_pos==0) then
                            write(*,*) "FAILED TO WRITE ANY EXTRA HEADERS TO OUTPUT"
                            write_card=.false.
                    endif
            enddo
            header(card_pos) = "COMMENT  "
            card_pos = card_pos + 1
            header(card_pos) = "COMMENT  List of cuts applied follows"
            card_pos = card_pos + 1
            header(card_pos) = "COMMENT  "
            card_pos = card_pos + 1
            if (write_card) then
                    do i=1,size(extra_headers)
                            if (card_pos>size(header)) then
                                    write(*,*) "FAILED TO WRITE SOME EXTRA HEADERS TO OUTPUT"
                                    exit
                            endif
                            header(card_pos) = extra_headers(i)
                            card_pos = card_pos + 1
                    enddo
            endif
    endif

    if (present(covariance)) then
            call add_card(header,"TTYPE4","COV_TT")
            call add_card(header,"TTYPE5","COV_QQ")
            call add_card(header,"TTYPE6","COV_UU")
            call add_card(header,"TTYPE7","COV_TQ")
            call add_card(header,"TTYPE8","COV_TU")
            call add_card(header,"TTYPE9","COV_QU")
            
            call add_card(header,"TUNIT4",unitsCode//"^2")
            call add_card(header,"TUNIT5",unitsCode//"^2")
            call add_card(header,"TUNIT6",unitsCode//"^2")
            call add_card(header,"TUNIT7",unitsCode//"^2")
            call add_card(header,"TUNIT8",unitsCode//"^2")
            call add_card(header,"TUNIT9",unitsCode//"^2")


    endif

    if (present(hitCount) .and. present(covariance)) then
            call add_card(header,"TTYPE10","N_OBS")
    elseif (present(hitCount)) then
            call add_card(header,"TTYPE4","N_OBS")
    endif

            
    call output_map(map,header,filename)

    deallocate(map)
    if (.not. (maps%has_T .and. maps%has_P)) deallocate(dummy_map)	

end subroutine ds_write_fits_map


subroutine ds_write_fits_covariance(filename,covariance,indices,nside,units,isRing,coordinate_system,extra_headers)
    use fitstools, only : output_map
    use head_fits, only : write_minimal_header, add_card
    type(ds_covariance) :: covariance
    integer(8), dimension(:) :: indices
    character(len=80), dimension(1:), optional :: extra_headers
    real(8), dimension(:,:), allocatable :: cov_data
    integer nside
    real(dp), dimension(:), pointer :: dummy_map
    character(*) :: filename
    integer card_pos
    logical write_card
    integer status, i
    character(1), optional :: coordinate_system
    logical, optional :: isRing
    integer(4) order_code
    character(8), optional :: units
    character(8) :: unitsCode
    integer :: nheal
    !integer,parameter :: map_type=2
    real(sp),allocatable,dimension(:,:) :: map
    character(len=80),dimension(1:600) :: header
    real(dp), parameter :: hpx_dbadval = -1.6375e30_dp
    character(len=20) :: dtype, coordsys
    character(len=8) :: ordering
    integer npix , ncov
    integer t
    call ds_assert(covariance%has_T .or. covariance%has_P, "covariance with neither T nor P passed to write_fits_covariance")

    npix = covariance%npix

    if (covariance%has_t .and. covariance%has_p) then
            ncov = 6
    else if (covariance%has_t) then
            ncov = 1
    else 
            ncov = 3	
    endif

    if (.not. (covariance%has_T .and. covariance%has_P)) then
            allocate(dummy_map(npix))
            dummy_map = 0
    endif

    if (present(coordinate_system)) then
            coordsys = coordinate_system
    else
            coordsys = 'C'
    endif



    status=0
    call deletefile(filename,status)

    if (.not. present(isRing)) then
       order_code = 2
            ordering='NESTED'
    else if (isRing) then
       order_code = 1
    ordering='RING'
    else
            ordering='NESTED'
       order_code=2
    endif

    if (.not. present(units)) then
       unitsCode = 'unknown '
    else
       unitsCode= units
    endif
        
    nheal= nside**2 *12
    allocate(cov_data(0:nheal-1,1:ncov))


    cov_data = hpx_dbadval

    if (covariance%has_t .and. covariance%has_p) then
            cov_data(indices,1) = covariance%TT
            cov_data(indices,2) = covariance%QQ
            cov_data(indices,3) = covariance%QU
            cov_data(indices,4) = covariance%TQ
            cov_data(indices,5) = covariance%TU
            cov_data(indices,6) = covariance%QU
    else if (covariance%has_t) then
            cov_data(indices,1) = covariance%TT
    else 
            cov_data(indices,1) = covariance%QQ
            cov_data(indices,2) = covariance%QU
            cov_data(indices,3) = covariance%QU
    endif	

    dtype= 'MAP'

    call write_minimal_header( header, dtype, nside=nside, &
            ordering=ordering, coordsys=coordsys,polar=covariance%has_p)

    if (present(extra_headers)) then
            
            card_pos = size(header)+1
            write_card=.true.
            do
                    card_pos=card_pos-1
                    if (header(card_pos) .ne. '') exit
                    if (card_pos==0) then
                            write(*,*) "FAILED TO WRITE ANY EXTRA HEADERS TO OUTPUT"
                            write_card=.false.
                    endif
            enddo
            header(card_pos) = "COMMENT  "
            card_pos = card_pos + 1
            header(card_pos) = "COMMENT  List of cuts applied follows"
            card_pos = card_pos + 1
            header(card_pos) = "COMMENT  "
            card_pos = card_pos + 1
            if (write_card) then
                    do i=1,size(extra_headers)
                            if (card_pos>size(header)) then
                                    write(*,*) "FAILED TO WRITE SOME EXTRA HEADERS TO OUTPUT"
                                    exit
                            endif
                            header(card_pos) = extra_headers(i)
                            card_pos = card_pos + 1
                    enddo
            endif
    endif


    if (covariance%has_t .and. covariance%has_p) then
            call add_card(header,"TTYPE1","COV_TT", update=.true.)
            call add_card(header,"TTYPE2","COV_QQ", update=.true.)
            call add_card(header,"TTYPE3","COV_UU", update=.true.)
            call add_card(header,"TTYPE4","COV_TQ", update=.true.)
            call add_card(header,"TTYPE5","COV_TU", update=.true.)
            call add_card(header,"TTYPE6","COV_QU", update=.true.)
            
            call add_card(header,"TUNIT1",unitsCode//"^2", update=.true.)
            call add_card(header,"TUNIT2",unitsCode//"^2", update=.true.)
            call add_card(header,"TUNIT3",unitsCode//"^2", update=.true.)
            call add_card(header,"TUNIT4",unitsCode//"^2", update=.true.)
            call add_card(header,"TUNIT5",unitsCode//"^2", update=.true.)
            call add_card(header,"TUNIT6",unitsCode//"^2", update=.true.)
            
    elseif (covariance%has_t) then
            call add_card(header,"TTYPE1","COV_TT", update=.true.)
            call add_card(header,"TUNIT1",unitsCode//"^2", update=.true.)
            
    else
            call add_card(header,"TTYPE1","COV_QQ", update=.true.)
            call add_card(header,"TTYPE2","COV_UU", update=.true.)
            call add_card(header,"TTYPE3","COV_QU", update=.true.)
            
            call add_card(header,"TUNIT1",unitsCode//"^2", update=.true.)
            call add_card(header,"TUNIT2",unitsCode//"^2", update=.true.)
            call add_card(header,"TUNIT3",unitsCode//"^2", update=.true.)
            
    endif

            
    call output_map(cov_data,header,filename)

    deallocate(cov_data)

end subroutine ds_write_fits_covariance



END MODULE DS_FITSTOOLS

module fits_helper
implicit none
interface get_column_by_name
	module procedure get_column_by_name_J_allocatable,get_column_by_name_E_allocatable,get_column_by_name_L_allocatable,get_column_by_name_I_allocatable,get_column_by_name_D_allocatable,get_column_by_name_B_allocatable
end interface
interface get_column_by_name_preallocated
	module procedure get_column_by_name_J,get_column_by_name_E,get_column_by_name_L,get_column_by_name_I,get_column_by_name_D,get_column_by_name_B
end interface
contains

subroutine fatal_fits_error(status,message)
	integer :: status
	character(*) :: message

	call report_fits_error(status,message)
	if (status .ne. 0) stop

end subroutine fatal_fits_error

subroutine report_fits_error(status,message)
	integer :: status
	character(*) :: message
	
	if (status .ne. 0) then
		write(*,*) "FITS Error No:", status
		write(*,*) "from operation: ", trim(message)
	endif
end subroutine report_fits_error


subroutine get_column_by_name_J_allocatable(unit,name,data,status,case_sensitive,preallocated)
	integer(4) :: unit
	character(*) :: name
	integer(4), allocatable, dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	logical(4), optional :: preallocated
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	
do_alloc = .true.
if (present(preallocated)) then
	do_alloc = .not. preallocated
endif
if (do_alloc) then
	allocate(data(nelements))
else
	if (.not. allocated(data)) stop "Preallocate flag set but array not allocated in get_column_by_name"
	if (size(data)<nelements) stop "Preallocated array too small in get_column_by_name"
endif


	call ftgcvJ(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_J(unit,name,data,status,case_sensitive)
	integer(4) :: unit
	character(*) :: name
	integer(4),  dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	

	call ftgcvJ(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_E_allocatable(unit,name,data,status,case_sensitive,preallocated)
	integer(4) :: unit
	character(*) :: name
	real(4), allocatable, dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	logical(4), optional :: preallocated
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	
do_alloc = .true.
if (present(preallocated)) then
	do_alloc = .not. preallocated
endif
if (do_alloc) then
	allocate(data(nelements))
else
	if (.not. allocated(data)) stop "Preallocate flag set but array not allocated in get_column_by_name"
	if (size(data)<nelements) stop "Preallocated array too small in get_column_by_name"
endif


	call ftgcvE(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_E(unit,name,data,status,case_sensitive)
	integer(4) :: unit
	character(*) :: name
	real(4),  dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	

	call ftgcvE(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_L_allocatable(unit,name,data,status,case_sensitive,preallocated)
	integer(4) :: unit
	character(*) :: name
	logical, allocatable, dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	logical(4), optional :: preallocated
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	
do_alloc = .true.
if (present(preallocated)) then
	do_alloc = .not. preallocated
endif
if (do_alloc) then
	allocate(data(nelements))
else
	if (.not. allocated(data)) stop "Preallocate flag set but array not allocated in get_column_by_name"
	if (size(data)<nelements) stop "Preallocated array too small in get_column_by_name"
endif


	call ftgcvL(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_L(unit,name,data,status,case_sensitive)
	integer(4) :: unit
	character(*) :: name
	logical,  dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	

	call ftgcvL(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_I_allocatable(unit,name,data,status,case_sensitive,preallocated)
	integer(4) :: unit
	character(*) :: name
	integer(2), allocatable, dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	logical(4), optional :: preallocated
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	
do_alloc = .true.
if (present(preallocated)) then
	do_alloc = .not. preallocated
endif
if (do_alloc) then
	allocate(data(nelements))
else
	if (.not. allocated(data)) stop "Preallocate flag set but array not allocated in get_column_by_name"
	if (size(data)<nelements) stop "Preallocated array too small in get_column_by_name"
endif


	call ftgcvI(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_I(unit,name,data,status,case_sensitive)
	integer(4) :: unit
	character(*) :: name
	integer(2),  dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	

	call ftgcvI(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_D_allocatable(unit,name,data,status,case_sensitive,preallocated)
	integer(4) :: unit
	character(*) :: name
	real(8), allocatable, dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	logical(4), optional :: preallocated
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	
do_alloc = .true.
if (present(preallocated)) then
	do_alloc = .not. preallocated
endif
if (do_alloc) then
	allocate(data(nelements))
else
	if (.not. allocated(data)) stop "Preallocate flag set but array not allocated in get_column_by_name"
	if (size(data)<nelements) stop "Preallocated array too small in get_column_by_name"
endif


	call ftgcvD(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_D(unit,name,data,status,case_sensitive)
	integer(4) :: unit
	character(*) :: name
	real(8),  dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	

	call ftgcvD(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_B_allocatable(unit,name,data,status,case_sensitive,preallocated)
	integer(4) :: unit
	character(*) :: name
	integer(1), allocatable, dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	logical(4), optional :: preallocated
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	
do_alloc = .true.
if (present(preallocated)) then
	do_alloc = .not. preallocated
endif
if (do_alloc) then
	allocate(data(nelements))
else
	if (.not. allocated(data)) stop "Preallocate flag set but array not allocated in get_column_by_name"
	if (size(data)<nelements) stop "Preallocated array too small in get_column_by_name"
endif


	call ftgcvB(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

subroutine get_column_by_name_B(unit,name,data,status,case_sensitive)
	integer(4) :: unit
	character(*) :: name
	integer(1),  dimension(:) :: data
	logical :: do_alloc
	character(8), parameter :: naxis2_keyword = "NAXIS2"
	integer nelements
	character(80) :: comment
	
	logical(4), optional :: case_sensitive
	
	
	integer(4) :: case_sen
	integer(4) colnum,status
	real(8) :: nullval = -1.6375E30
	logical :: any_undefined_values
	
	status=0
	case_sen = 0
	if (present(case_sensitive)) then
		if (case_sensitive) case_sen = 1
	endif
	call ftgcno(unit,case_sen,name,colnum,status)
	call report_fits_error(status, "Could not find column: "//trim(name))
	if (status .ne. 0) return
	nelements=0
	comment=""
	call ftgkyj(unit,naxis2_keyword,nelements,comment,status)
	call report_fits_error(status, "Could not read naxis2 keyword")
	if (status .ne. 0) return

	

	call ftgcvB(unit,colnum,1,1,nelements,nullval,data(1:nelements),any_undefined_values,status)
	call report_fits_error(status, "Could not load column data "//trim(name))
	if (status .ne. 0) return
	
end subroutine 

end module fits_helper
