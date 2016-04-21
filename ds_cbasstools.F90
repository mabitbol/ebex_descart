! This file is the main one to change when creating a new executable for a new experiment.
!
! Its job is to read in the data and get it in the format used by the rest of the code.
! It also makes a naive map as it goes along.  Since there are a number of pieces of auxiliary information
! that need to be loaded in this code is a little complicated.

!Stolen from some nice person on comp.lang.fortran circa 2006


module ds_cbass_fitstools
    use ds_types
    use ds_utils
    use ds_simple_prior
    use healpix_types
    use pix_tools
    use ds_multidetector
    use ds_cbass_option_utils
    use fits_helper

#ifndef __GFORTRAN__
    use ieee_arithmetic
#endif

    use fl_lists
    implicit none

    real(dp), parameter :: DEGRA = 0.0174532925
    integer, parameter :: NCHANNEL = 2

    integer, parameter :: FITS_READ_ONLY = 0
    integer, parameter :: FITS_CASE_INSENSITIVE = 0
    character(8), parameter :: NAXIS2_KEYWORD = "NAXIS2  "
    type ds_moduleScanInfo
        integer id1,id2,id3
        integer owner
        character(len=256) filename
        integer first_index,last_index
        integer modScan !The local index into moduleScans
        type(ds_moduleScanInfo),pointer :: next
        integer n
        type(ds_noiseInfo) :: noise
    end type ds_moduleScanInfo
    type ds_file_data
        real(sp), dimension(:,:), allocatable :: channel_data
        integer, dimension(:,:), allocatable :: flagged
        real(dp), dimension(:), allocatable :: MJD,Ra,Dec,Az,El,theta,parangle,hwp_angle
        integer n
    end type ds_file_data
    type ds_moduleScanList
        type(ds_moduleScanInfo), pointer :: first,last
        integer length
    end type ds_moduleScanList
    type ds_detector_list
        integer nactive
        integer, dimension(NCHANNEL) :: active
        real(dp),dimension(NCHANNEL) :: delta_ra, delta_dec, delta_parangle
    end type ds_detector_list
    character(len=512), private :: message
    contains

subroutine readDataAssignWork(fileList,scans,noiseInfo,correlator,originalIndices,maps,covariance,opt,output_fits_cards,hitCountOutput)
	! fileList - input, array of strings of length 256
	! scans - output, array of (module,scan) pairs (ds_modulescan)
	! noiseInfo - output, array of noise info structures (ds_noiseinfo)
	! correlator - input/output.  object containing data on which processor owns what, and noise correlations.  On input, should have %comm set.  Other stuff set on output
	! originalIndices - output.  contains map from our pixel numbering to healpix.
	! maps - output.  (I,Q,U map.)  Contains naive map on output.
	! covariance - output.  the naive map covariance. (ds_covariance)
	! opt - input.  The options file.
	! output_fits_cards - output.  List of cards to propagate to the FITS header of the output map.
	!This is the only subroutine that you need to re-write to use 

	implicit none
	type(ds_moduleScan), pointer, dimension(:) :: scans
	type(ds_noiseInfo), pointer, dimension(:) :: noiseInfo
	type(ds_correlator) :: correlator
	integer(i8b), allocatable, dimension(:) :: originalIndices
	type(ds_trimap) :: maps, buffer_maps
	type(ds_covariance) :: covariance
	type(ds_cbass_options) :: opt
	integer(i8b) :: full_map_pixels
	character(256), dimension(:) :: fileList
	type(fl_string80_list) :: output_fits_cards
	integer progress
	!Information about the data to be loaded.
	integer nmodules
	integer(i4b) npix
	!Work assignment information
	integer rank,nproc
	integer nmin,nextra
	integer nScan,my_nScan
	integer, allocatable, dimension(:) :: nScan_proc
        !Things needed to load the data from file
	type(ds_file_data) :: full_file_data
	logical haveBadScans, thisScanBad
	integer status
	integer ordering
	integer scan_number_in_file
	character(256) :: current_filename
	type(ds_detector_list) :: detector_list
	!Things needed to process the data and put it in the moduleScan structure.
	integer ntod,na
	type(ds_moduleScan), pointer :: scan
	type(ds_moduleScanInfo), pointer :: scanInfo
	!Iteration variables and misc
	integer(i4b) :: unit,one, ierr,dummy
	integer ms,i,j,t
	real(dp) :: mjd
	integer my_ms
	integer(i8b), allocatable, dimension(:) :: hitCountOutput
	type(ds_moduleScanList) :: scanList

	!Set up the size and rank information and save it in the correlator too.
	!This requires that correlator%comm has already been set in the driver.
	unit = ds_get_lun()

#ifndef NO_MPI
	call MPI_Comm_size(correlator%comm, nproc,ierr)
	call MPI_Comm_rank(correlator%comm, rank, ierr)
#else
	nproc=1
	rank=0
	ierr=0
#endif

	correlator%nproc = nproc
	correlator%proc = rank

   	npix = nside2npix(opt%nside)
	dummy=1

#ifndef NO_MPI
	call MPI_Allreduce(MPI_IN_PLACE, dummy, 1, MPI_INTEGER,MPI_SUM,correlator%comm,ierr)
#endif

	!GET AND DIVIDE LIST OF SCANS
	!Turn the target_list data into a form slightly more useful for us.
	!The moduleScanList wraps a linked list of ds_moduleScanInfo types, each of which 
	!contains the module number, run, segment, etc.

	call buildModuleScanList(scanList, fileList, detector_list, opt, output_fits_cards)

	!This is the total number of moduleScans for all the processors
	nScan = scanList%length
	if(rank==0) then
	write(message,*) "Total number of (module,scans) = ", nScan
            call ds_log(message,ds_feedback_quiet)
	endif

	if (nScan==0) then
	    write(*,*) "NO SCANS!  CANNOT MAKE A MAP WITHOUT SCANS!"
#ifndef NO_MPI		    
            call MPI_Finalize(ierr)
#endif
	    stop
	endif		

	!Determine the number of moduleScans each process is responsible for.
	!Each process knows this number for each process.
	!nmin is the minimum number of moduleScans held by a process.
	!Some have more to make up the difference.
	allocate(nScan_proc(0:nproc-1))
	nmin=floor((1.0*nScan)/nproc)
	nextra = nScan - nmin*nproc
	nScan_proc = nmin
	call ds_assert(nextra>=0,"Broken nextra - this is Joe's fault.")
	if (nextra>0) nScan_proc(0:nextra-1) = nmin +1

	!Some useful checks and feedback about the allocation of moduleScans
	call ds_assert(sum(nScan_proc) == nScan, "Mistake in assignment of dets - not all assigned.")

	scanInfo=>scanList%first
	do i=0,nproc-1
            do j=0,nScan_proc(i)-1
                scanInfo%owner=i
                scanInfo=>scanInfo%next
            enddo
	enddo
	
	my_nScan = nScan_proc(rank)
	correlator%my_nmodules = my_nScan
	if(rank==0) then
            write(message,*) "Approx number of modules per proc:", my_nScan
            call ds_log(message,ds_feedback_quiet)
	endif

        ! Setup space for this proc
        !Now we have that info we check how many moduleScans this process has and allocate space
        !and initialize the pointers inside
        
        scanInfo=>scanList%first
        allocate(scans(0:my_nScan-1))
        allocate(noiseInfo(0:my_nscan-1))
        my_ms = 0
        do ms=0,nScan-1
            if (scanInfo%owner == rank) then
                scan=>scans(my_ms)
                call setup_moduleScan(scan)
                call buildLeakageInfo(scan)
                nullify(scan%theta)
                nullify(scan%pointing) !DWPS: always nullify pointers when you make them
                allocate(noiseInfo(my_ms)%sigma)
                allocate(noiseInfo(my_ms)%fknee)
                allocate(noiseInfo(my_ms)%alpha)
                my_ms = my_ms + 1
            endif
            scanInfo=>scanInfo%next
        enddo

        call prepareTriMap(maps,npix, opt%do_temperature, opt%do_polarization, zero_based = .true.)  !Healpix maps are zero-based.  For now this is a healpix map.
        call setupScans(scans,scanList,rank)

        if(rank==0 .and. opt%data_prior) call ds_log("Data-prior mode active: priors will be computed direct from the data. This may be slower.", ds_feedback_quiet)
        if (rank==0) call ds_log_milestone("DATA_SELECTION")

        scanInfo=>scanList%first
        haveBadScans=.false.
        progress=0
        current_filename=""

        !Loop through the moduleScanList.  Check if this proc is responsible for the corresponding moduleScan.
        !If so, read it in and process it.
        do 
            progress = progress +1
            ! End loop if finished
            ! Cycle loop if not owner of this scan
            if (.not. scanInfo%owner==rank) then
                if (.not. associated(scanInfo%next)) exit
                scanInfo=>scanInfo%next
                cycle
            endif

            ! Load new FITS file
            call read_detector_list_from_fits(scanInfo%filename,detector_list)

            !Load the new file data in if this is a new filename.
            if (.not. current_filename==scanInfo%filename  )     then
                call destroy_file_data(full_file_data)
                scan_number_in_file = 1
                call load_full_file(scanInfo%filename,detector_list,full_file_data,opt)
                write(message,'("Loaded file: ",A)') ,trim(scanInfo%filename)
                call ds_log(message,ds_feedback_noisy)
                if (rank==0) then
                    write(message,'(A,F7.1,A)') "Approximate load progress: ", (progress*100.0)/my_nScan,"%"
                    call ds_log(message,ds_feedback_noisy)
                endif
                thisScanBad=.false.
                current_filename = scanInfo%filename
            endif
            
            ms = scanInfo%modScan
            scan=>scans(ms)
            write(message, '(A,I3)') " - Loaded scan ", scan_number_in_file
            call ds_log(message,ds_feedback_debug)
            scan_number_in_file = scan_number_in_file + 1
            ntod = scanInfo%last_index - scanInfo%first_index + 1
            if (opt%offsetLength .ne. -1) then
                na = ntod/opt%offsetLength
                ntod = na * opt%offsetLength
            endif
            scan%ntod=ntod

            call buildNoiseInfo(scanInfo,NoiseInfo(ms),mjd,opt)
            call get_scan_from_fits(scanInfo,full_file_data,detector_list,noiseinfo(ms),scan,opt)
            if (opt%data_prior) call prepare_one_data_prior(correlator,opt%offsetLength,noiseinfo(ms),scan)
            call make_inv_Cw(noiseInfo(ms),scan)
            call invCw_mult_tod(scan)
            call add2rhs(scan,maps)

            call deprojectTimestreamOntoOffset(scan%timestreams,scan%offsets,scan%flagged)
            call destroyTimestream(scan%timestreams)
            if (.not. associated(scanInfo%next)) exit
            scanInfo=>scanInfo%next
        enddo

        if (allocated(nScan_proc)) deallocate(nScan_proc)
        write(message,*) "Rank ", rank," loaded all data."
        call ds_log(message,ds_feedback_debug)
        call ds_assert(.not. haveBadScans, "Bad scans reported")
        if (rank==0) call ds_log_milestone("DATA_LOADED")
        full_map_pixels = npix

        call repixelizeData(correlator,scans,full_map_pixels,originalIndices,maps,opt,hitCountOutput)

        if (rank==0) call ds_log_milestone("REPIXELIZED")
        if (opt%do_temperature) then
            npix = maps%T%npix
        else 
            npix = maps%Q%npix
        endif
        
        !Sum the accumulated maps accross all the processes, thereby summing for all modules for all scans
#ifndef NO_MPI		    
        if (opt%do_temperature)  call MPI_Allreduce(MPI_IN_PLACE, maps%T%map, maps%T%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
        if (opt%do_polarization) call MPI_Allreduce(MPI_IN_PLACE, maps%Q%map, maps%Q%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
        if (opt%do_polarization) call MPI_Allreduce(MPI_IN_PLACE, maps%U%map, maps%U%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
#endif

        call make_naivecov(npix,covariance,correlator,scans, opt%do_temperature, opt%do_polarization)
        call cov_mult(maps,covariance)

        if (.not. opt%naive_mode) then
            do ms=0,correlator%my_nmodules-1
                scan=>scans(ms)
                call prepareTimestream(scan%timestreams,scan%ntod )
                
                call map2tod(scan,maps)
                call invCw_mult_tod(scan)

                call subtractTimestreamFromOffset(scan%timestreams,scan%offsets)
                call destroyTimestream(scan%timestreams)
                call destroyFlagged(scan%flagged)
            enddo
        endif
        
        call destroyModuleScanList(scanList)

end subroutine readDataAssignWork


function filename_extension(filename) result(s)
	character(*) :: filename
	character(8) :: s
	!Get the final "." in the filename
	integer n
	integer l
	
	l=len_trim(filename)
	n = scan(filename,".",back=.true.)
	s(:) = filename(n+1:l)
end function filename_extension


function filename_extension_is_fits(filename)
	character(*) :: filename
	logical filename_extension_is_fits
	filename_extension_is_fits = .false.
	if (string_ends_with(filename,"fits") .or.  &
	    string_ends_with(filename,"FITS") .or.  &
	    string_ends_with(filename,"fit" ) .or.   &
	    string_ends_with(filename,"FIT" )) filename_extension_is_fits=.true.
end function


subroutine load_full_file(filename,detector_list,data,opt)
	type(ds_file_data) :: data
	type(ds_cbass_options) :: opt
	type(ds_detector_list) :: detector_list
	character(*) :: filename

	call ds_assert(opt%do_temperature .or. opt%do_polarization, "Elected to load neither temperature nor polarization in load_full_file")
	call load_full_file_fits(filename,detector_list,data,opt)
end subroutine


subroutine load_full_file_fits(filename,detector_list,data,opt)
	type(ds_file_data) :: data
	type(ds_cbass_options) :: opt
	type(ds_detector_list) :: detector_list
	character(*) :: filename
	logical :: do_T, do_P
	integer :: unit,status,hdutype
	integer n,i
	logical anyf
	character(80) :: comment
	integer nchannels, column_number, flag_num
	integer :: first_row = 1
	integer :: first_element = 1
	real(sp), parameter :: null_value = 0.0_4
	integer ndim, interval
	character(32) :: column_name, flag_column
	integer, dimension(2) :: first_pixels, last_pixels, dimensions	
	real(sp), allocatable, dimension(:) :: column
        integer, allocatable, dimension(:) :: buffer_int	
	status = 0
	!Open the file
	call FTGIOU(unit, status)
	call FTNOPN(unit,trim(filename),FITS_READ_ONLY,status)

	!Now move to the extension and get the data	
	call FTMAHD(unit,2,hdutype,status)
	call FTGKYJ(unit,NAXIS2_KEYWORD,n,comment,status)
	data%n = n
	call allocate_and_get(unit,n,"MJD",data%MJD,status)
	call allocate_and_get(unit,n,"RA",data%RA,status)
	call allocate_and_get(unit,n,"DEC",data%DEC,status)
	!might not need this - check this out.
	call allocate_and_get(unit,n,"AZIMUTH",data%AZ,status)
	call allocate_and_get(unit,n,"ELEVATIO",data%EL,status)
	call allocate_and_get(unit,n,"PARANGLE",data%parangle,status)
        call allocate_and_get(unit,n,"HWP",data%hwp_angle,status)

	!Covnert into radians from degrees.
	!Because Healpix expects this.
	data%RA = data%RA * DEGRA
	data%DEC = data%DEC * DEGRA
	data%AZ = data%AZ * DEGRA
	data%EL = data%EL * DEGRA
	data%parangle=data%parangle * DEGRA
        data%hwp_angle=data%hwp_angle * DEGRA

	nchannels = NCHANNEL
	!Load the complete channel data.  Have swapped around the order from the old cbass code
	!Because we load the data in as a complete array instead of as a series of columns.
	!Have altered the bit where the data is taken out accordingly.
	allocate(data%channel_data(nchannels,data%n))
        allocate(data%flagged(nchannels,data%n))
	data%channel_data = 0.0
        data%flagged = 1
	
	!Okay, so now we will be loading the channel data from the fits file.
	!We could do this by loading the entire data set in one, and then just
	!ignore the data we do not need, or we could loop through and only do the 
	!ones that we need.  There will be some threshold for which is faster.
	!We should keep that at the top somewhere and experiment with it if needed.
	!But for now just load detector-by-detector
	do i=1,NCHANNEL
            if (detector_list%active(i) .ne. 1) cycle
            write(column_name,'("CH",I0)') i
            write(flag_column,'("FLAG",I0)') i
            call get_column_by_name(unit,flag_column,buffer_int,status)
            data%flagged(i,:) = buffer_int
            !write(*,*) "flag buff", data%flagged(i,1)
            deallocate(buffer_int)
            call FTGCNO(unit,FITS_CASE_INSENSITIVE,column_name,column_number,status)
            !call FTGCNO(unit,FITS_CASE_INSENSITIVE,flag_column,flag_num,status)
            first_row=1
            first_element=1
            call FTGCVE(unit,column_number,first_row,first_element,data%n,null_value, &
                   data%channel_data(i,:),anyf,status)
            !call FTGCVI(unit,flag_num,first_row,first_element,data%n,null_value, &
                  ! data%flagged(i,:),anyf,status)
	enddo	
	!Check for errors and close the file.
	call fatal_fits_error(status, "Error reading from file:"//trim(filename))
	call FTCLOS(unit,status)
	call FTFIOU(unit,status)
	call report_fits_error(status, "Closing file or releasing unit.")
end subroutine load_full_file_fits


subroutine read_detector_list_from_fits(fits_file,list)
	character(*) :: fits_file
	type(ds_detector_list) :: list
	integer :: unit,status,hdutype,n,i
	real(dp),dimension(:),allocatable :: buffer_real
	integer, dimension(:),allocatable :: buffer_int
	
	status=0
	n=0
	
	!Open file to read from
	call FTGIOU(unit,status)
	call FTNOPN(unit,trim(fits_file),FITS_READ_ONLY,status)
	!Move to extension with detector info
	call FTMAHD(unit,4,hdutype,status)
	
	!Copy info into detector_list
	call get_column_by_name(unit,"ONOFF",buffer_int,status)
	list%active=buffer_int
	deallocate(buffer_int)
	
	call get_column_by_name(unit,"AZOFF",buffer_real,status)
	list%delta_ra=buffer_real
	deallocate(buffer_real)
	
        call get_column_by_name(unit,"ELOFF",buffer_real,status)
	list%delta_dec=buffer_real
	deallocate(buffer_real)
	
	call get_column_by_name(unit,"DELTA_PARANGLE",buffer_real,status)
	list%delta_parangle=buffer_real
	deallocate(buffer_real)
	
	!Close file
	call FTCLOS(unit,status)
	call FTFIOU(unit,status)
	
	do i=1,NCHANNEL
            if(list%active(i)==1) then
                n=n+1
            endif 
	enddo
	
	list%nactive=n
end subroutine read_detector_list_from_fits


subroutine allocate_and_get(unit,size,name,data,status)
	   integer unit,size
	   character(*) :: name
	   real(dp), dimension(:), allocatable :: data
		integer status
		if (status .ne. 0) return
	   call get_column_by_name(unit,name,data,status)
end subroutine 


subroutine destroy_file_data(F)
	type(ds_file_data) :: F
	if (allocated(F%channel_data)) deallocate(F%channel_data)
        if (allocated(F%flagged)) deallocate(F%flagged)
	if (allocated(F%MJD)) deallocate(F%MJD)
	if (allocated(F%RA)) deallocate(F%RA)
	if (allocated(F%Dec)) deallocate(F%DEC)
	if (allocated(F%Az)) deallocate(F%AZ)
	if (allocated(F%El)) deallocate(F%EL)
	if (allocated(F%theta)) deallocate(F%THETA)
	if (allocated(F%parangle)) deallocate(F%parangle)
        if (allocated(F%hwp_angle)) deallocate(F%hwp_angle)
	F%n = -1
end subroutine 


subroutine load_noise_columns(unit, nscan, noise_data, status, opt)
	integer unit,status
	integer nscan
	integer file_format
	integer, parameter :: number_noise_parameters = 3
!	character(len=5), dimension(number_noise_parameters), parameter :: noise_parameters = (/ "sigma", "fknee" ,"alpha" /)
	integer i,d
	logical anyf
	character(16) :: column_name
	real(dp), allocatable, dimension(:) :: noise_column
	real(dp), allocatable, dimension(:,:,:) :: noise_data
	type(ds_cbass_options) :: opt
	integer channel, nchannels
	character(*), parameter :: channel_data_name = "DATA"
	integer column_number
	integer, parameter :: first_row = 1
	integer, parameter :: first_element = 1
	real(sp), parameter :: null_value = 0.0_4
	
	nchannels = NCHANNEL
	allocate(noise_data(nchannels,nscan,number_noise_parameters))
	noise_data = -1.0

	call FTGCNO(unit,FITS_CASE_INSENSITIVE,"SIGMA",column_number,status)
	call FTGCVD(unit, column_number, first_row, first_element, nscan*nchannels, null_value, noise_data(:,:,1), anyf, status)
	call FTGCNO(unit,FITS_CASE_INSENSITIVE,"FKNEE",column_number,status)
	call FTGCVD(unit, column_number, first_row, first_element, nscan*nchannels, null_value, noise_data(:,:,2), anyf, status)
	call FTGCNO(unit,FITS_CASE_INSENSITIVE,"ALPHA",column_number,status)
	call FTGCVD(unit, column_number, first_row, first_element, nscan*nchannels, null_value, noise_data(:,:,3), anyf, status)
end subroutine


subroutine get_scan_from_fits(info,full_data,detector_list,noise,moduleScan,opt)
	type(ds_cbass_options) :: opt
	type(ds_moduleScanInfo) :: info
	type(ds_noiseInfo) :: noise
	type(ds_moduleScan) :: moduleScan
	type(ds_file_data) :: full_data
	type(ds_detector_list) :: detector_list
	real(dp) :: mu
	real(dp), dimension(2)::linearFit
	integer, parameter :: psi_field = 1
	logical, save :: firstTime=.true.
	real(dp) :: delta_psi
	integer start,finish
	integer pixel
	integer k
	real(dp) :: theta,phi
	integer ntod,na,i,t
	logical :: do_T, do_P
	integer offsetLength
	integer detector_index
	real(dp) :: delta_ra, delta_dec

        real(dp), allocatable, dimension(:) :: rot_roll
	
	do_T = opt%do_temperature
	do_P = opt%do_polarization
	
	!Set up the ID numbers in the moduleScans so we can keep track of them later
	moduleScan%id(1) = info%id1
	moduleScan%id(2) = info%id2
	moduleScan%id(3) = info%id3
	
        write(*,*) "Detector, ", info%id3
	detector_index = info%id3
	delta_ra = detector_list%delta_ra(detector_index) * DEGRA
	delta_dec = detector_list%delta_dec(detector_index) * DEGRA
	
	moduleScan%has_T = do_T
	moduleScan%has_P = do_P
	
	!Determine the size of this timestream
	
	start =  info%first_index
	finish = info%last_index
	ntod = finish - start + 1
	
	call ds_assert(ntod>0,"Zero-sized (or negative) TOD passed to get_scan_from_fits")
	! If offsetLength has special value -1 then
	! we just use a single offset
	if (opt%offsetLength==-1) then
            na = 1
            offsetLength = ntod
	else
            na = ntod/opt%offsetLength
            ntod = na * opt%offsetLength
            offsetLength = opt%offsetLength
	endif
	finish = start+ntod-1
	
	!Set up the offsets
        call prepareOffsets(moduleScan%offsets,na,opt%offsetLength)
	
	!Set up the pixel pointing
	allocate(moduleScan%pointing(ntod))
	call compute_pixel_number(modulescan%pointing,full_data,start,opt,delta_ra,delta_dec,rot_roll)
	moduleScan%inv_Cw=0.0

	if (do_P) then
            allocate(moduleScan%theta(ntod))
            !Get theta from FITS
            call get_parallactic_angle(ntod, rot_roll(start:finish), full_data%hwp_angle(start:finish), moduleScan%theta)
	endif

        call prepareTimestream(moduleScan%timestreams,ntod)
        call prepareFlagged(moduleScan%flagged,ntod)
        moduleScan%timestreams%timestream = full_data%channel_data(detector_index,start:finish)
        moduleScan%flagged%flagged = full_data%flagged(detector_index,start:finish)
        do t=1,ntod
            moduleScan%timestreams%timestream(t) = moduleScan%timestreams%timestream(t)
        enddo
end subroutine get_scan_from_fits


subroutine read_applied_cuts(unit, cards)
	integer unit
	type(fl_string80_list) :: cards
	character(80) :: card
	character(8) :: keyword
	integer i, status
	i=0
	call FTPMRK
	keyword=""
	status=0
	do
            write(keyword,'("DSCUT",I0.3)') i
            card=""
            call FTGCRD(unit,keyword,card,status)
            if (status .ne. 0) exit
            call fl_append_string80_list(cards,card)
            i=i+1
	enddo
	call FTCMRK
end subroutine read_applied_cuts


subroutine buildModuleScanList(modScanList,filename_list, detector_list, opt, cards)
! The job of this function is to build up a list of scans to be read.
! These are not yet read at this stage but will be later.
! The modScanList object is a linked-list of modScanInfo objects.


	type(ds_moduleScanList) :: modScanList
	type(ds_cbass_options) :: opt
	type(fl_string80_list), optional :: cards
	type(ds_detector_list) :: detector_list

	character(256), dimension(:) :: filename_list
	character(256) :: filename
	integer id1,id2,id3
	integer unit,hdutype
	integer nelements
	character(80) :: comment
	integer first_index, last_index, f, status
	character(8) :: keyword
	character(8) :: value
	logical :: do_T, do_P
	integer(4), allocatable, dimension(:) :: starts, ends
	integer nscan,s
	integer i
	integer number_of_hdu, cuts_hdu
	integer, parameter :: number_noise_parameters = 3
	real(dp), allocatable, dimension(:,:,:) :: noise_data
	integer det

	do_T = opt%do_temperature
	do_P = opt%do_polarization

	call init_moduleScanList(modScanList)

	!Loop through the file names in the list.
	do f=1,size(filename_list)
		filename = filename_list(f)
                call read_detector_list_from_fits(filename,detector_list) 	
		if (f==1) then
			call module_scans_in_file(filename,opt,starts,ends,noise_data,status,cards)
		else
			call module_scans_in_file(filename,opt,starts,ends,noise_data,status)
		endif

		if (status .ne. 0) cycle

		! Set up a few numbers.  id1 is just a reference number
		nscan = size(starts)
		id1 = f
                write(*,*) nscan	
		
		!Go through the chunks we found and record its stats in the modScanList.
		do s=1,nscan
                    do det=1,NCHANNEL
                        !Include only detectors listed as active in the detector list file
                        if (detector_list%active(det) .ne. 1) cycle
                        id2=s
                        id3=det
                        first_index = starts(s)
                        last_index = ends(s)
                        !Create and append the object.
                
                        call appendModuleScanInfo(modScanList,id1,id2,id3,first_index,last_index,filename)
                        allocate(modScanList%last%noise%sigma)
                        allocate(modScanList%last%noise%fknee)
                        allocate(modScanList%last%noise%alpha)
                
                        modScanList%last%noise%sigma = noise_data(det,s,1)
                        modScanList%last%noise%fknee = noise_data(det,s,2)
                        modScanList%last%noise%alpha = noise_data(det,s,3)
                    enddo
		enddo
		if (allocated(starts)) deallocate(starts)
		if (allocated(ends)) deallocate(ends)
		if (allocated(noise_data)) deallocate(noise_data)
	enddo
end subroutine buildModuleScanList


subroutine module_scans_in_file(filename, opt, starts, ends, noise_data, status,cards)
	character(*) :: filename
	integer status
	type(ds_cbass_options) :: opt
	integer(4), allocatable, dimension(:) :: starts, ends
	real(dp), allocatable, dimension(:,:,:) :: noise_data
	type(fl_string80_list), optional :: cards
	call module_scans_in_fits_file(filename,opt,starts,ends,noise_data,status,cards)
end subroutine


subroutine module_scans_in_fits_file(filename,opt,starts,ends,noise_data,status,cards)
	integer number_of_hdu, cuts_hdu
	integer, parameter :: number_noise_parameters = 3
	character(*) :: filename
	integer status
	type(ds_cbass_options) :: opt
	integer unit, hdutype
	integer nelements
	character(80) :: comment
	integer(4), allocatable, dimension(:) :: starts, ends
	real(dp), allocatable, dimension(:,:,:) :: noise_data
	type(fl_string80_list), optional :: cards
	nelements = 0
	status=0
	! Open the FITS file and count the number of headers.
	call FTGIOU(unit, status)
	call FTNOPN(unit,trim(filename),FITS_READ_ONLY,status)
	call FTTHDU(unit, number_of_hdu, status)

	if (status==0 .and. number_of_hdu .gt. 2) then
            !Choose which chunk extension to use based on the options supplied in the parameter file.
            cuts_hdu=3
            !Move to the right header
            call FTMAHD(unit,cuts_hdu,hdutype,status)
            !Check how many chunks there are here
            comment=""
            nelements=0
            call FTGKYJ(unit,NAXIS2_KEYWORD,nelements,comment,status)
            if (status .eq. 0 .and. nelements==0) then
                !If no chunks (all data cut) just skip this file and move on to the next
                call ds_log("File "//trim(filename)//" has no scans in - skipping", ds_feedback_quiet)
                call FTCLOS(unit,status)
                call FTFIOU(unit,status)
                status=1				
                return
            endif
            !Read the START and END columns, and the noise (FKNEE, etc) columns.
            if (status==0) call get_column_by_name(unit,"START",starts,status)
            if (status==0) call get_column_by_name(unit,"END",ends,status)
            call load_noise_columns(unit,nelements,noise_data,status,opt)
            call FTCLOS(unit,status)
            call FTFIOU(unit,status)
            if (status==0) starts = starts + 1
	endif

	if (status .ne. 0) then
		write(*,*) "WARNING!: Could not open FITS file or wrong format.  Skipping: ",trim(filename),status
		if (allocated(starts)) deallocate(starts)
		if (allocated(ends)) deallocate(ends)
		if (allocated(noise_data)) deallocate(noise_data)
	endif
end subroutine


subroutine init_moduleScanList(L)
	type(ds_moduleScanList) :: L
	nullify(L%first)
	nullify(L%last)
	L%length=0
end subroutine init_moduleScanList


function makeModuleScanInfo(id1,id2,id3,first_index, last_index,filename,nt,np) result(output)
	type(ds_moduleScanInfo), pointer :: output
	integer id1,id2,id3
	integer first_index, last_index
	character(len=256) filename
	integer nt, np

	allocate(output)
	output%first_index = first_index
	output%last_index = last_index
	output%id1=id1
	output%id2=id2
	output%id3=id3
	output%filename=filename
	nullify(output%next)
end function makeModuleScanInfo


subroutine appendModuleScanInfo(list,id1,id2,id3,first_index, last_index,filename)
	type(ds_moduleScanList) :: list
	integer id1,id2,id3
	integer first_index, last_index
	character(len=256) filename
	type(ds_moduleScanInfo), pointer :: info
	integer nt, np
	info => makeModuleScanInfo(id1,id2,id3,first_index, last_index,filename,nt, np)
	if (list%length==0) then
		list%first=>info
		list%last=>info
		list%length=1
	else
		list%last%next=>info
		list%last=>info
		list%length=list%length+1
	endif
end subroutine appendModuleScanInfo


subroutine destroyModuleScanList(list)
	type(ds_moduleScanList) :: list
	type(ds_moduleScanInfo), pointer :: info,next
	info => list%first
	if (.not. associated(info)) return
	do
		next=>info%next
		deallocate(info)
		if (.not. associated(next)) exit
		info=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	list%length=0
end subroutine destroyModuleScanList


subroutine repixelizeData(correlator,moduleScans,maxIndex,originalIndices,maps,opt,hitCountOutput)
	use udgrade_nr
	!This subroutine resets the pointing of a ds_pointing instance so that the 
	!pixel indices it uses run from 1..npix (and can therefore be used as array indices)
	!The maxIndex argument specifies the number of pixels in the input class,
	!for example if you are converting a set of healpix indices then maxIndex should be 12*Nside^2
	type(ds_correlator) :: correlator
	type(ds_cbass_options) :: opt
	type(ds_modulescan), dimension(0:correlator%my_nmodules-1) :: moduleScans
	integer(i8b) :: maxIndex
	integer min_hits
	integer(i8b), allocatable, dimension(:) :: originalIndices
	integer(i8b), allocatable, dimension(:) :: hitCountOutput
	type(ds_trimap) :: maps
	integer (i4b), allocatable, dimension(:) :: hitCount
	integer i,d,n,p,q,t,nbad
	integer ierror

	min_hits = opt%min_hits
	!Each process builds its own hit count map in healpix format
	allocate(hitCount(0:maxIndex-1))
	
	hitCount=0
	do i=0,correlator%my_nmodules-1
            do t=1,moduleScans(i)%ntod
                p=moduleScans(i)%pointing(t)
                !if (moduleScans(i)%flagged(d)%flagged(t) .eq. 0 ) cycle
                if (p>=0) hitCount(p)=hitCount(p)+1
            enddo
	enddo

	!The processes all sum their hit count map to get the global hit count map.
	ierror=0
#ifndef NO_MPI		    		
	call MPI_AllReduce(MPI_IN_PLACE,hitCount,maxIndex,MPI_INTEGER,MPI_SUM,correlator%comm,ierror)
	call ds_assert(ierror==0,"Error in MPI_AllReduce in repixeizeData")
#endif

	!If requested, the hit count map is saved to file
	if ((trim(opt%hits_filename).ne."") .and. correlator%proc==0) then
            open(file=trim(opt%hits_filename),unit=25)
            do p=0,maxIndex-1
                write(25,*) p,hitCount(p)
            enddo
	endif

	!Count the number of hit pixels to get the total npix
	n=0
	do p=0,maxIndex-1
            if (hitCount(p)>0) n=n+1
	enddo

	!added by DWPS - code to handle degenerate pixels
	!if the pixels have only a single hit, they are bad. So set pointing to bad_pixel flag for those
	!pixels (ds_mapping knows how to handle them).  Also remove from the map, by removing number of 
	!bad pixels from npix and setting the hitcount for those pixels to 0.
	nbad=0
	do p=0,maxIndex-1	
            if(hitCount(p)==1 .or. (hitCount(p)<min_hits .and. hitCount(p)>0)) nbad = nbad+1
	enddo
	do i=0,correlator%my_nmodules-1
            do t=1,moduleScans(i)%ntod
                p=moduleScans(i)%pointing(t)
                if (p==bad_pixel) cycle
                if (hitCount(p)==1 .or. (hitCount(p)<min_hits .and. hitCount(p)>0)) then
                    moduleScans(i)%pointing(t) = bad_pixel
                endif
            enddo
	enddo

	allocate(hitCountOutput(0:maxIndex-1))
	hitCountOutput = hitCount
	!remove bad pixels from the map
	forall(p=0:maxIndex-1 , hitCount(p)==1 .or. (hitCount(p)<min_hits .and. hitCount(p)>0)) hitCount(p)=0
	!remove bad pixels from npix	
	n=n-nbad

	write(message,*) "Built global hit map: npix = ",n, " pixels cut=", nbad
	if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)
	!The root process should preserve the mapping that takes the new indices back to the old ones.
	!The originalIndices array stores that data.
	if (correlator%proc==0) then
            allocate(originalIndices(n))
            originalIndices=-1
            q=1
            do p=0,maxIndex-1
                if (hitCount(p)>0) then
                    originalIndices(q)=p
                    q=q+1
                endif
            enddo
            write(*,*) minval(originalIndices)
	endif
	!Relabel the hit count map so that it contains the new index for each healpix pixel
	!instead of the hit count.  This will be useful for re-pixelizing the pointings and maps
	q=1
	do p=0,maxIndex-1
            if (hitCount(p)>0) then
                hitCount(p)=q
                q=q+1
            else
                hitCount(p)=-1
            endif
	enddo
	!Re-pixelize the pointings in the moduleScans
	do i=0,correlator%my_nmodules-1
            do t=1,moduleScans(i)%ntod
                p=moduleScans(i)%pointing(t)
                if (p>=0) then
                    moduleScans(i)%pointing(t) = hitCount(p)
                else 
                    moduleScans(i)%pointing(t) = bad_pixel
                endif
            enddo
	enddo

	!Repixelize the maps, reducing them in size from the full healpix sized ones to
	!ones containing only the pixels hit somewhere in one of the timestreams
	if (maps%has_T) call repixelizeMap(maps%T,hitCount,n)
	if (maps%has_P) call repixelizeMap(maps%Q,hitCount,n)
	if (maps%has_P) call repixelizeMap(maps%U,hitCount,n)
	
	deallocate(hitCount)
end subroutine repixelizeData


subroutine repixelizeMap(map,new_pix,npix)
	type(ds_map) :: map
	integer, dimension(0:map%npix-1) :: new_pix
	integer npix
	integer p,i
	integer, dimension(:), allocatable :: new_indices
	real(dp), dimension(:), allocatable :: new_map
	
	allocate(new_map(npix))
	allocate(new_indices(npix))
	do i=0,map%npix-1
            p = new_pix(i)
            if (p .ge. 1) then
                new_map(p) = map%map(i)
                new_indices(p) = i
            endif
	enddo
	call destroyMap(map)
	call prepareMap(map,npix)
	map%map = new_map
	map%indices = new_indices
	deallocate(new_map)
	deallocate(new_indices)
end subroutine
	

subroutine setupScans(scans,scanList,rank)
	!Perform misc setup
	type(ds_moduleScanList) :: scanList
	type(ds_moduleScan), pointer, dimension(:) :: scans
	integer :: rank
	type(ds_moduleScan), pointer :: scan
	type(ds_moduleScanInfo), pointer :: scanInfo
	integer ms,m,run,seg

	scanInfo=>scanList%first
	ms=0
	do 	
            if (scanInfo%owner .ne. rank) then
                if (.not. associated(scanInfo%next)) exit
                scanInfo=>scanInfo%next
                cycle		
            endif
            scan=>scans(ms)
            scanInfo%modscan=ms
            if (.not. associated(scanInfo%next)) exit
            scanInfo=>scanInfo%next
            ms=ms+1				
	enddo
end subroutine setupScans


subroutine buildNoiseInfo(info,noise,mjd,opt)
	type(ds_moduleScanInfo) :: info
	type(ds_noiseInfo) :: noise
	type(ds_cbass_options) :: opt
	integer, parameter :: number_of_noise_parameters = 3
	real(dp) :: mjd
	
	call ds_assert(associated(noise%sigma),"Noise array (sigma) not associated in buildNoiseInfo")
	call ds_assert(associated(noise%fknee),"Noise array (fknee) not associated in buildNoiseInfo")
	call ds_assert(associated(noise%alpha),"Noise array (alpha) not associated in buildNoiseInfo")

        noise%sigma = info%noise%sigma
        noise%alpha = info%noise%alpha
        noise%fknee = info%noise%fknee

        call testValueForEvil(noise%sigma,"sigma",zero_is_bad=.true.)
        call testValueForEvil(noise%fknee,"fknee",zero_is_bad=.false.)
        call testValueForEvil(noise%alpha,"alpha",zero_is_bad=.false.)
	
	deallocate(info%noise%sigma)
	deallocate(info%noise%alpha)
	deallocate(info%noise%fknee)
end subroutine buildNoiseInfo


subroutine buildLeakageInfo(scan)
	type(ds_modulescan) :: scan

        scan%leakage(1) = 1.0
        scan%leakage(2) = 0.5
        scan%has_leakage_matrix = .true.
end subroutine buildLeakageInfo


subroutine testValueForEvil(value,variable_name,zero_is_bad)
	real(dp) :: value
	character(*) :: variable_name
	character(256) ::  message
	logical, optional :: zero_is_bad
	
#ifdef USE_IEEE
	write(message,*) variable_name, " is NaN :", value
	call ds_assert(.not. ieee_is_nan(value),message)
	write(message,*) variable_name, " is infinite: ", value
	call ds_assert(ieee_is_finite(value),message)
#else
        write(message,*) variable_name, " is NaN/infinite :", value
	call ds_assert(ds_isfinite(value),message)
#endif

	if (present(zero_is_bad)) then
		if (zero_is_bad) then
			write(message,*) variable_name, " is zero"
			call ds_assert(value/=0,message)
		endif
	endif
end subroutine testValueForEvil


subroutine get_parallactic_angle(n,parangle,hwp,par)
	integer n,i
	real(dp), dimension(1:n) :: parangle,hwp,par

	do i=1,n
            par(i) = 2.0*hwp(i) + parangle(i)
	enddo
end subroutine


subroutine compute_pixel_number(pointing,full_data,start,opt, delta_ra, delta_dec, rot_roll)
	type(ds_file_data) :: full_data
	type(ds_cbass_options) :: opt
	real(dp) :: delta_ra, delta_dec
	real(dp) :: ra,dec
	real(dp) :: delta_ra_rot,delta_dec_rot
	integer :: status, pixel_x, pixel_y, pixel
	real(dp) :: el, theta, phi
	integer, dimension(:) :: pointing
	integer ntod, k, npix_x, t
	integer start
	integer ierr
	real(dp) :: el_min, el_max, el_delta
        real(dp), allocatable, dimension(:) :: rot_ra, rot_dec, rot_roll
	
	ntod = size(pointing)
        call offset_rotation(full_data%ra, full_data%dec, full_data%parangle, delta_ra, delta_dec, rot_ra, rot_dec, rot_roll)
	do t=1,ntod
            k = start+t-1
            dec = rot_dec(k)
            ra = rot_ra(k)
            theta = HALFPI - dec
            phi = ra
            call ang2pix_ring(opt%nside,theta,phi,pixel)
            pointing(t) = pixel
	enddo
end subroutine 


subroutine offset_rotation(ra, dec, roll, az_off, el_off, rot_ra, rot_dec, rot_roll)
    real(dp), dimension(:) :: ra, dec, roll
    real(dp), dimension(3,3) :: rotation
    real(dp), allocatable, dimension(:,:,:) :: attitude_mat, detector_mat
    real(dp), allocatable, dimension(:) :: rot_ra, rot_dec, rot_roll !output
    real(dp) :: az_off, el_off
    real(dp) :: alpha, beta
    real(dp) :: pi
    integer :: N, i

    pi = 3.1415926535897932384626
    alpha = -az_off
    beta = el_off
    N = size(ra)

    allocate( attitude_mat(N,3,3))
    allocate( detector_mat(N,3,3))
    allocate( rot_ra(N))
    allocate( rot_dec(N))
    allocate( rot_roll(N))

    rotation(1,1) = cos(alpha)
    rotation(1,2) = -sin(alpha)*cos(beta)
    rotation(1,3) = sin(alpha)*sin(beta)
    rotation(2,1) = sin(alpha)
    rotation(2,2) = cos(alpha)*cos(beta)
    rotation(2,3) = -cos(alpha)*sin(beta)
    rotation(3,1) = 0.
    rotation(3,2) = sin(beta)
    rotation(3,3) = cos(beta)

    do i=1,N
        attitude_mat(i,1,1) = cos(roll(i))*cos(ra(i)) - sin(roll(i))*sin(dec(i))*sin(ra(i))
        attitude_mat(i,2,1) = cos(roll(i))*sin(ra(i)) + sin(roll(i))*sin(dec(i))*cos(ra(i))
        attitude_mat(i,3,1) = -sin(roll(i))*cos(dec(i))
        attitude_mat(i,1,2) = -cos(dec(i))*sin(ra(i))
        attitude_mat(i,2,2) = cos(dec(i))*cos(ra(i))
        attitude_mat(i,3,2) = sin(dec(i))
        attitude_mat(i,1,3) = sin(roll(i))*cos(ra(i)) + cos(roll(i))*sin(dec(i))*sin(ra(i))
        attitude_mat(i,2,3) = sin(roll(i))*sin(ra(i)) - cos(roll(i))*sin(dec(i))*cos(ra(i))
        attitude_mat(i,3,3) = cos(roll(i))*cos(dec(i))
    end do

    do i=1,N
        detector_mat(i,:,:) = matmul(attitude_mat(i,:,:), rotation)
    end do

    do i=1,N
        rot_dec(i) = asin(detector_mat(i,3,2))
        rot_ra(i) = -atan2(detector_mat(i,1,2),detector_mat(i,2,2))
        rot_roll(i) = -atan2(detector_mat(i,3,1), detector_mat(i,3,3))
    end do
end subroutine


end module ds_cbass_fitstools
