program descart_cbass

!use mpi
use ds_types
use ds_fitstools
use ds_multidetector
use ds_solver
use inifile
use ds_utils
use ds_simple_prior
use ds_cbass_fitstools
use ds_cbass_option_utils
implicit none

character(256) :: arg
integer :: arg_count
integer(i4b)       :: i, j, k, l, m, n, q, p
integer(i4b)       :: myid, numprocs, ierr, root
integer(i4b)       :: nside, nmaps
character(len=512) :: message
type(ds_cbass_options) :: options
type(ds_modulescan), pointer, dimension(:) :: offsetTarget, offsetSolution
type(ds_trimap) :: data_maps, offset_maps
type(ds_correlator) :: correlator
type(ds_covariance) :: map_covariance
integer :: nd, d
real(dp), parameter :: sampling_rate = 190.73486328125
real(dp), parameter :: nyquist = sampling_rate / 2
character(256) :: offset_dir
integer ndet_total
integer(i8b), allocatable, dimension(:) :: originalIndices
integer :: totalNpix
character(128) :: cwd
logical feedback_set
integer xsize, ysize
integer(i8b), allocatable, dimension(:) :: hitCount
type(ds_noiseinfo), pointer, dimension(:) :: noiseInfo
character(256), dimension(:), allocatable :: file_list
type(fl_string80_list) :: fits_output_cards


! Initialize MPI environment
ierr=0
call ds_init_milestone()                                                                                    
#ifndef NO_MPI
call mpi_init(ierr)
if (ierr .ne. 0) then
    write(*,*) "Failed to initialize MPI!"
    write(*,*) "Something is wrong with your MPI setup"
endif
call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
#else
myid=0
numprocs=1
ierr=0
#endif
root = 0
! Read the command line.
arg_count=iargc()
write(*,*) "args:", arg_count
call flush(6)
if (arg_count==0) then
    write(*,*) 'Syntax : descart_cbass params.ini'
    write(*,*) 'The list of files will be given by the parameter file_list from the ini file'
    write(*,*) 'See default_params.ini for an explanation of the parameters'
    call flush(6)
#ifndef NO_MPI
    call MPI_Finalize(ierr)
#endif
    stop
endif
call flush(6)
#ifndef NO_MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
call flush(6)
if (myid==0) then
   call getarg(1,arg)
   write(*,*) "Root arg:", arg
endif
!call get_command_argument(1,arg)
#ifndef NO_MPI
call MPI_Bcast(arg,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
call ds_log("Reading options from: "//trim(arg),ds_feedback_quiet)
call read_options(arg,options)
!Read any command line options.
call command_line_override_options(options, 2)
call read_file_list(options%file_list, file_list)
!This parameter controls how many log messages are sent to the stdout.
ds_global_feedback = options%verbosity
!Setup the MPI communicator.
#ifndef NO_MPI
correlator%comm = MPI_COMM_WORLD
#else
correlator%comm = 0
#endif
correlator%traditionalMode = options%traditional_mode
if (myid==root) call ds_log_milestone("MODULE_SETUP")
call fl_create_string80_list(fits_output_cards)

call readDataAssignWork(file_list,offsetTarget,noiseInfo,correlator,originalIndices,data_maps, map_covariance,options, fits_output_cards, hitCount)

if (myid==root) call ds_log_milestone("DATA_LOAD_COMPLETE")
if(options%naive_mode) then
    if (myid==0) call ds_log("IN NAIVE MODE! NO PCG.",ds_feedback_quiet)
    if(correlator%proc==0) then
            call save_map_any_format(data_maps, hitCount, options%output_filename, options, fits_output_cards)
    endif
else 
    if(correlator%proc==0) then
        call save_map_any_format(data_maps, hitCount, "naive_" // options%output_filename, options, fits_output_cards)
    endif
    ! Set up the noise prior.  If traditional mode is on this does nothing.
    if (.not. options%data_prior) call prepare_prior(correlator,options%offsetLength,nyquist,noiseinfo,offsetTarget)
    if (myid==root) call ds_log_milestone("PRIOR_PREPARED")
    if (myid==root) call ds_log("Ready for PCG.", ds_feedback_quiet)

    ! This subroutine fills in the offsetSolution structure.
    call PCG(correlator,offsetTarget,offsetSolution,map_covariance,options%pcg_tolerance,options%pcg_iterations)
    if (myid==root) call ds_log_milestone("PCG_COMPLETE")

    ! Reclaim the memory used in the target data since we now have the solution.
    if (options%save_offsets) call saveAllOffsetsToFiles(offsetSolution,options%offset_dir)
    call destroy_modulescan(offsetTarget,deallocateoffsets=.true.)

    ! We need to wait for all the processes to finish destriping to make maps.
#ifndef NO_MPI
    call mpi_barrier(correlator%comm,ierr)
#endif
    if (myid==0) call ds_log("Completed PCG",ds_feedback_quiet)
    ! Make a naive map from all the calculated offsets.
    ! This map represents the contribution to the naive maps of correlated noise.
    ! Make navie map in multidetector

    call makeNaiveMap(correlator,offsetSolution,offset_maps,map_covariance,options%do_temperature,options%do_polarization)

    !The destriped map is the difference between the naive map and the offset map.
    !Subtract to get the destriped maps.

    call subtractTriMap(data_maps,offset_maps)
    
    if(correlator%proc==0) then
            call save_map_any_format(data_maps, hitCount, options%output_filename, options, fits_output_cards)
    endif
endif

if(correlator%proc==0 .and. options%save_covariance) then
    call save_cov_any_format(options%covariance_filename, originalIndices, map_covariance, options, fits_output_cards)
endif

#ifndef NO_MPI
call MPI_Finalize(ierr)
#endif
if (myid==root) call ds_log_milestone("END")


contains 

subroutine save_map_any_format(maps, hitCount, filename, opt, cards)
	type(ds_trimap) :: maps
	type(ds_covariance) :: covariance
	integer(8), dimension(1:) :: hitCount
	type(ds_cbass_options) :: opt
	type(fl_string80_list) :: cards
	character(1) :: coordinate_system
	character(*) :: filename
	character(len=80), dimension(:), pointer :: cards_array
	
        call ds_log("Saving map in NESTED ordering", ds_feedback_quiet)		
	!Add the options used to the FITS cards to be put in the output map.
	call options_to_fits_cards(opt, cards)
	cards_array => fl_string80_list_to_array(cards) 
        !Write the map as a healpix map
        call ds_write_fits_map(filename, maps, opt%nside, units='unknown ',isRing=opt%ring_ordering, &
            extra_headers=cards_array, hitCount=hitCount)
	deallocate(cards_array)
end subroutine save_map_any_format

subroutine save_cov_any_format(filename, originalIndices, covariance, opt, cards)
	character(*) :: filename
	integer(8), dimension(:) :: originalIndices 
	type(ds_covariance) :: covariance
	type(ds_cbass_options) :: opt
	type(fl_string80_list) :: cards
	character(len=80), dimension(:), pointer :: cards_array
	character(1) :: coordinate_system
	
	cards_array => fl_string80_list_to_array(cards) 	
        call ds_write_fits_covariance(filename,covariance,originalIndices,opt%nside, units="unknown ",isRing=opt%ring_ordering,&
                                        extra_headers=cards_array)
	deallocate(cards_array)
end subroutine

subroutine read_file_list(filename, file_list)
    !Read a list of strings from a file into an array.
    character(*) :: filename
    character(256), dimension(:), allocatable :: file_list
    character(512) :: line
    integer unit,io
    !Count the number of lines in the file by reading each in turn until end of file
    !Count only lines the are 
    unit=ds_get_lun()
    open(unit=unit,file=filename)
    n=0
    do
        read(unit,'(A)',iostat=io) line
        if (io<0) exit
        if (trim(line) .ne. "") then
            line = adjustl(line)
            if(line(1:1) .ne. '#') n=n+1
        endif
    enddo
    rewind(unit)
    call ds_assert(n>0,"No uncommented written lines found in file: "//trim(filename))
    allocate(file_list(n))
    n=1
    do
        read(unit,'(A)',iostat=io) line
        if (io<0) exit
        if (trim(line) .ne. "") then
            line = adjustl(line)
            if(line(1:1) .ne. '#') then
                file_list(n)=trim(line)
                n=n+1
            endif
        endif
    enddo
end subroutine read_file_list

end program
