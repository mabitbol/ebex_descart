module ds_types
  use healpix_types
  use ds_utils
  implicit none
  

real(kind=DP), parameter, public :: HOUR2RAD = TWOPI / 24.0_DP
integer, parameter :: MAGIC_VALUE = 250582
real(dp), parameter :: IGNORE_THIS_VALUE = 2.0
  

type ds_map
    real(dp), allocatable, dimension(:) :: map
    integer, allocatable, dimension(:) :: indices
    integer :: npix
end type ds_map
  
type ds_trimap
    type(ds_map) :: T,Q,U
    logical :: has_p, has_t
end type ds_trimap

type ds_detpointing
    integer, allocatable, dimension(:) :: pixel
    real(dp), allocatable, dimension(:) :: theta
    integer :: npix
    integer :: nt
    type(ds_map) :: globalHitMap
end type ds_detpointing
  
type ds_modulescan
    integer :: ntod
    integer magic_check
    real(dp) :: inv_Cw 
    integer,pointer,dimension(:) :: pointing
    real(dp),pointer,dimension(:) :: theta !radians!
    logical :: has_p, has_t
    !type(ds_flagged),pointer, dimension(:) :: flagged
    !type(ds_timestream),pointer, dimension(:) :: timestreams !will have dimension ndiodes_max
    !type(ds_offsets),pointer, dimension(:) :: offsets  !Will rename later.
    type(ds_flagged),pointer :: flagged
    type(ds_timestream),pointer :: timestreams !will have dimension ndiodes_max
    type(ds_offsets),pointer :: offsets  !Will rename later.
    integer, dimension(3) :: id !JAZ - run, scan, seg
    real(dp) :: scanfreq
    logical :: has_leakage_matrix
    real(dp), dimension(:), pointer :: leakage
end type ds_modulescan

type ds_covariance
    real(dp), allocatable, dimension(:) :: TT
    real(dp), allocatable, dimension(:) :: QQ
    real(dp), allocatable, dimension(:) :: UU
    real(dp), allocatable, dimension(:) :: TQ
    real(dp), allocatable, dimension(:) :: TU
    real(dp), allocatable, dimension(:) :: QU
    logical :: has_p, has_t
    integer npix
end type ds_covariance

type ds_timestream
    integer :: nt, my_module
    real(dp), allocatable, dimension(:) :: timestream
end type ds_timestream

type ds_flagged
    integer :: nt, my_module
    integer, allocatable, dimension(:) :: flagged
end type ds_flagged
  
type ds_offsets
    integer :: na, my_module
    real(dp), allocatable, dimension(:) :: values
    integer :: length
    complex(dpc), pointer, dimension(:) :: Qinv
    complex(dpc), pointer, dimension(:) :: Pinv
    logical ownsCorrelators
end type ds_offsets
  
  
type ds_converter
    integer, allocatable, dimension(:) :: pixelTable
    integer length
    integer currentIndex_
    logical final_
end type ds_converter

type ds_noiseinfo
    real(dp),pointer :: sigma, fknee, alpha
end type ds_noiseinfo

interface prepareOffsets
    module procedure prepareOffsets_stan
end interface
 
interface saveOffsetsToFile
    module procedure saveOffsetsToFileName, saveOffsetsToFileNum
end interface
  
contains

subroutine setup_moduleScan(ms)
        type(ds_modulescan) :: ms
        real dummy

        dummy = -1
        ms%magic_check = MAGIC_VALUE

        !allocate(ms%inv_Cw)
        allocate(ms%flagged)
        allocate(ms%timestreams)
        allocate(ms%offsets)
        allocate(ms%leakage(3))
        ms%leakage = sqrt(dummy)
        ms%has_leakage_matrix=.false. !By default, make this false to preserve current behviour

end subroutine setup_moduleScan


subroutine prepare_covariance(cov, npix, has_t, has_p)
        type(ds_covariance) :: cov
        logical has_t, has_p
        integer npix
        
        cov%has_t = has_t
        cov%has_p = has_p
        cov%npix = npix
        if (has_t) then
            allocate(cov%TT(npix))
            cov%TT = 0
        endif
        if (has_p) then
            allocate(cov%QQ(npix))
            allocate(cov%UU(npix))
            allocate(cov%QU(npix))
            cov%QQ = 0
            cov%UU = 0
            cov%QU = 0
        endif
        if (has_p .and. has_t) then
            allocate(cov%TQ(npix))
            allocate(cov%TU(npix))
            cov%TQ = 0
            cov%TU = 0
        endif
end subroutine prepare_covariance


subroutine destroy_covariance(cov)
        type(ds_covariance) :: cov
        cov%npix = -1
        if (allocated(cov%TT)) deallocate(cov%TT)
        if (allocated(cov%TQ)) deallocate(cov%TQ)
        if (allocated(cov%TU)) deallocate(cov%TU)
        if (allocated(cov%QQ)) deallocate(cov%QQ)
        if (allocated(cov%QU)) deallocate(cov%QU)
        if (allocated(cov%UU)) deallocate(cov%UU)
        
        cov%has_t = .false.
        cov%has_p = .false.
end subroutine destroy_covariance


!JZ Notes
!Using the offset-owned correlators.
!In your code to build/load the correlators, allocate
!complex(dpc) vectors with the "target" attribute
!for Pinv and Qinv.
!Call the setOffsetCorrelators(offset,Pinv,Qinv) routine
!Do not deallocate your Pinv and Qinv.  Just let them be.
!They will be deallocated when the destroyOffsets is called.


!!                                            !!
!!    SUBROUTINES TO DO WITH THE DATA TYPES   !!
!!                                            !!


subroutine prepareOffsets_stan(offsets,na,offsetLength)
        type(ds_offsets) :: offsets
        integer, intent(in) :: na, offsetLength

        offsets%length = offsetLength
        offsets%na = na
        call ds_checkAllocate(offsets%values, offsets%na)
        nullify(offsets%Qinv)
        nullify(offsets%Pinv)
        offsets%ownsCorrelators=.false.
end subroutine prepareOffsets_stan


subroutine saveOffsetsToFileName(offsets,filename)
	type(ds_offsets) :: offsets
        character(*) :: filename
	integer unit

	unit=ds_get_lun()
	open(unit=unit,file=filename,form='unformatted',action='write')
	call saveOffsetsToFileNum(offsets,unit)
	close(unit)
end subroutine saveOffsetsToFileName


subroutine saveOffsetsToFileNum(offsets,unit)
	type(ds_offsets) :: offsets
	integer unit

	write(unit) offsets%length
	write(unit) offsets%na
	write(unit) offsets%values
end subroutine saveOffsetsToFileNum


subroutine loadOffsetsFromFile(offsets,filename)
	type(ds_offsets) :: offsets
        character(*) :: filename
	integer unit
	integer na,lc

	unit=ds_get_lun()
	open(unit=unit,file=filename,form='unformatted',action='read')
	read(unit) lc
	read(unit) na
	call prepareOffsets(offsets,na,lc)
	read(unit) offsets%values
	close(unit)
end subroutine loadOffsetsFromFile

  
subroutine copyOffsets(source, dest)
        type(ds_offsets) :: source, dest
        integer i

        call prepareOffsets(dest, source%na, source%length)
        do i= 1, source%na
            dest%values(i) = source%values(i)
        enddo
        if (associated(source%Qinv)) dest%Qinv=>source%Qinv
        if (associated(source%Pinv)) dest%Pinv=>source%Pinv
        dest%ownsCorrelators=.false.
end subroutine copyOffsets


subroutine setOffsetCorrelators(offset,pinv,qinv)
	complex(dpc), target, dimension(:) :: pinv,qinv
	type(ds_offsets) :: offset
	integer n
	
	n=size(pinv)
	call ds_assert(2*n + 1 == offset%na, "Wrong size correlators in setOffsetCorrelators")
	call ds_assert(size(qinv) == n, "Different size correlators in setOffsetCorrelators")
	offset%pinv => pinv
	offset%qinv => qinv
	offset%ownsCorrelators = .true.
end subroutine setOffsetCorrelators
  
  
subroutine destroyOffsets(offsets)
        type(ds_offsets),intent(inout) :: offsets

        if(allocated(offsets%values)) deallocate(offsets%values)
        offsets%na= -1
        offsets%length= -1
        if (offsets%ownsCorrelators) then
            deallocate(offsets%qinv)
            deallocate(offsets%pinv)
        endif
        if (associated(offsets%Qinv)) nullify(offsets%Qinv)
        if (associated(offsets%Pinv)) nullify(offsets%Pinv)
end subroutine destroyOffsets


subroutine prepareFlagged(flagged,nt)
        type(ds_flagged) :: flagged
        integer :: nt
        call ds_checkAllocate(flagged%flagged, nt)
        flagged%nt = nt
end subroutine prepareFlagged


subroutine prepareTimestream(timestream,nt)
        type(ds_timestream) :: timestream
        integer :: nt
        call ds_checkAllocate(timestream%timestream, nt)
        timestream%nt = nt
end subroutine prepareTimestream
  

subroutine copyTimestream(source, dest)
    !Copy a timestream into another.
    !If the destination is not allocated with the right size, (re-)allocate it.
    type(ds_timestream) :: source, dest
    integer i
    call prepareTimestream(dest,source%nt)
    do i= 1, source%nt
        dest%timestream(i) = source%timestream(i)
    enddo
end subroutine copyTimestream
 

subroutine copyFlagged(source, dest)
    type(ds_flagged) :: source, dest
    integer i
    call prepareFlagged(dest,source%nt)
    do i= 1, source%nt
        dest%flagged(i) = source%flagged(i)
    enddo
end subroutine copyFlagged


subroutine destroyTimestream(timestream)
    type(ds_timestream) :: timestream
    if (allocated(timestream%timestream)) deallocate(timestream%timestream)
    timestream%nt = -1
end subroutine destroyTimestream

subroutine destroyFlagged(flagged)
    type(ds_flagged) :: flagged
    if (allocated(flagged%flagged)) deallocate(flagged%flagged)
    flagged%nt = -1
end subroutine destroyFlagged

  
subroutine prepareMap(map, npix, zero_based)
    type(ds_map) :: map
    integer npix
    logical, optional :: zero_based
    integer :: start, finish
    
    start = 1
    if (present(zero_based)) then
            if (zero_based) start=0
    endif
    finish=start+npix-1

    allocate(map%map(start:finish))
    allocate(map%indices(start:finish))
    map%map = 0.0D0
    map%indices = 0
    map%npix = npix
end subroutine prepareMap


subroutine prepareTriMap(maps,npix,temperature, polarization, zero_based)
    type(ds_trimap) :: maps
    integer npix
    logical temperature, polarization
    logical, optional :: zero_based
    logical :: zero_based_
    maps%has_t = temperature
    maps%has_p = polarization

    zero_based_ = .false.
    if (present(zero_based)) zero_based_ = zero_based
    maps%T%npix=0
    maps%Q%npix=0
    maps%U%npix=0

    if (temperature)  call prepareMap(maps%T,npix,zero_based_)
    if (polarization) call prepareMap(maps%Q,npix,zero_based_)
    if (polarization) call prepareMap(maps%U,npix,zero_based_)
end subroutine prepareTriMap


subroutine destroyMap(map)
    type(ds_map) :: map
    if (allocated(map%map)) deallocate(map%map)
    if (allocated(map%indices)) deallocate(map%indices)
    map%npix=-1
end subroutine destroyMap
  

subroutine destroyTriMap(maps)
    type(ds_trimap) :: maps
    call destroyMap(maps%T)
    call destroyMap(maps%Q)
    call destroyMap(maps%U)

end subroutine destroyTriMap
	

subroutine copyTriMap(source,dest)
    type(ds_trimap) :: source, dest
    if (source%has_t) call copyMap(source%T,dest%T)
    if (source%has_p) call copyMap(source%Q,dest%Q)
    if (source%has_p) call copyMap(source%U,dest%U)
end subroutine copyTriMap


subroutine copyMap(source, dest)
    type(ds_map) :: source, dest
    integer i
    call prepareMap(dest, source%npix)
    do i= 1, source%npix
        dest%indices(i) = source%indices(i)
        dest%map(i) = source%map(i)
    enddo
end subroutine copyMap
  

subroutine preparePointing(pointing, nt, np,do_hitmap)
    type(ds_detpointing) :: pointing
    integer nt
    integer, optional :: np
    logical, optional :: do_hitmap
    integer npix
    
    if (present(np)) then
            npix = np
    else
            npix = 0
    endif
    
    call ds_checkAllocate(pointing%pixel, nt)
    call ds_checkAllocate(pointing%theta, nt)
    pointing%nt = nt
    pointing%npix = npix
    if (present(do_hitmap)) then 
        if (do_hitmap) then
            call prepareMap(pointing%globalHitMap, npix)
            pointing%globalHitMap%map = 0
        endif
    endif 
end subroutine preparePointing
  

subroutine destroyPointing(pointing)
    type(ds_detpointing) :: pointing
    if (allocated(pointing%pixel)) deallocate(pointing%pixel)
    if (allocated(pointing%theta)) deallocate(pointing%theta)
    pointing%npix=-1
    pointing%nt = -1
    call destroyMap(pointing%globalHitMap)
end subroutine destroyPointing


subroutine copyPointing(source, dest, do_hitmap)
    type(ds_detpointing) :: source,dest
    logical,optional :: do_hitmap
    integer i

    call ds_checkAllocate(dest%pixel,source%nt)
    call ds_checkAllocate(dest%theta,source%nt)
    
    do i= 1,source%nt
        dest%pixel(i) = source%pixel(i)
        dest%theta(i) = source%theta(i)
    enddo
    dest%npix = source%npix
    dest%nt = source%nt
end subroutine copyPointing
  

subroutine copy_modulescan(source,dest,duplicateoffsets,duplicatetimestreams,shareoffsets,sharetimestreams)
!copies modulescan object.  The pointing and theta pointers are shared.  By default, the offsets and
!timestream pointers are nullified. To share either, set shareoffsets and/or sharetimestreams to 1.
!To duplicate either (make own hard copy with allocated memory), set duplicateoffsets and/or 
!duplicatetimestreams to 1.

    type(ds_modulescan),intent(in) :: source
    type(ds_modulescan),intent(inout) :: dest
    logical,optional :: duplicateoffsets
    logical,optional :: duplicatetimestreams
    logical,optional :: shareoffsets
    logical,optional :: sharetimestreams
    integer i

    call check_modulescan_initialized(dest, "Copy modulescan")
    !eventually add an assertion that sharing and copying cannot both be true.
    if (present(sharetimestreams).and. present(duplicatetimestreams))  then
        call ds_assert(.not.(sharetimestreams .and. duplicatetimestreams),"Cannot both copy and share timestreams")
    endif
    if (present(shareoffsets).and. present(duplicateoffsets))  then
        call ds_assert(.not.(shareoffsets .and. duplicateoffsets),"Cannot both copy and share timestreams")
    endif

    dest%ntod    = source%ntod
    dest%id      = source%id
    dest%inv_Cw  = source%inv_Cw
    dest%pointing => source%pointing
    dest%theta => source%theta
    dest%has_T = source%has_T
    dest%has_P = source%has_P

    dest%has_leakage_matrix = source%has_leakage_matrix
    dest%leakage = source%leakage
	
    call copyFlagged(source%flagged,dest%flagged)

    if(present(sharetimestreams)) then
        if(sharetimestreams) then
            call ds_assert(.false., "Coding Error: sharing of timestreams is broken in copy_modulescan")
            dest%timestreams => source%timestreams
        endif
    endif

    if(present(duplicatetimestreams)) then
        if(duplicatetimestreams) then
            call copyTimestream(source%timestreams,dest%timestreams)
        endif
    endif
    
    if(present(shareoffsets)) then
        call ds_assert(.false., "Coding Error: sharing of offsets is broken in copy_modulescan")
        if(shareoffsets) dest%offsets => source%offsets
    endif

    if(present(duplicateoffsets)) then
        if(duplicateoffsets) then
            call copyOffsets(source%offsets,dest%offsets)
        endif
    endif

end subroutine copy_modulescan


subroutine destroy_modulescan(arg,deallocateoffsets)
!destroys pointer array of modulescans. Should work on an instance too.
!Default is to nullify offset pointer. Set deallocateoffsets=.true. to deallocate
    type(ds_modulescan),pointer,dimension(:) :: arg
    logical,optional :: deallocateoffsets
    integer :: m, nmodules
    logical :: deall
    type(ds_modulescan), pointer :: A

    deall= .false.
    if(present(deallocateoffsets)) deall= deallocateoffsets

    if(.not.associated(arg)) return

    nmodules= size(arg)
    do m= 0,nmodules-1
        A => arg(m)
        A%id=-1
        A%magic_check = 0
        !destroy instance of arg
        nullify(A%theta)
        nullify(A%pointing)
        if(deall) then
            if(associated(A%offsets)) deallocate(A%offsets)
        else
            nullify(A%offsets)
        endif
        if (associated(A%timestreams)) deallocate(A%timestreams)
        nullify(A%timestreams)
        if (associated(A%flagged)) deallocate(A%flagged)
        nullify(A%flagged)
        if (associated(A%leakage)) deallocate(A%leakage)
        nullify(A%leakage)
    enddo
    
    deallocate(arg)

end subroutine destroy_modulescan


subroutine check_modulescan_initialized(ms, message)
    type(ds_modulescan) :: ms
    character(*) :: message
    
    call ds_assert(ms%magic_check==MAGIC_VALUE, message)
end subroutine check_modulescan_initialized


subroutine saveAllOffsetsToFiles(moduleScans,dir)
    type(ds_modulescan), pointer, dimension(:) :: moduleScans
    character(*) dir
    type(ds_modulescan), pointer :: moduleScan
    integer i,d
    character(512) :: filename
    do i=0,size(moduleScans)-1
        moduleScan=>moduleScans(i)
        filename=filenameForSavedOffset(dir,moduleScan%id(1),moduleScan%id(2))
        call saveOffsetsToFile(moduleScan%offsets,filename)
    enddo
end subroutine saveAllOffsetsToFiles


function filenameForSavedOffset(directory,run,scan) result(f)
    character(512) :: f
    character(*) :: directory
    character(*), parameter :: fmt = '( A, "/", I5.5, "_",I5.5,"_", I5.5, ".off"  )'
    integer :: run, scan
    write(f,fmt) trim(directory),run,scan
end function filenameForSavedOffset


end module
