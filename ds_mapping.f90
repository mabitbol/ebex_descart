!Once pointing information has been worked out and loaded into the modulescan type:
!  call make_naivecov to produce the naive covariance matrix
!When a naive map is required, initialise q and u map arrays and
!  call add2rhs in a loop to add each modulescan's data to the sum
!  call cov_mult on the maps to finish the map
!When you are finished making naive maps
!  call deallocate_cov

!NB: the inv_Cw matrix will have to be used in place of whiteNoise in simple_prior and driver too.
!Use invCW_mult to apply it, only after call to make_naivecov

!For data read-in, require pointing info (and angles) to be read first, to run make_naivecov.  Can then read-in
!the data and add2rhs as we go along, producing the first map. We must also convert each timestream into offset
!space at read-in. After these 2 ops are done on it, we can forget it. Z op is factorised thus:
! Zd = ( F^t - F^t P naive) d
!where naive includes Cw^-1 ops where appropriate.

!Support for bad diodes: set the diode_flag array of modulescan structure to zero AND set the sigma 
!of that diode to zero.  This will cause it to (i) be rejected from the mapping ops and (ii) have no
!effect on the accepted timestreams when inv_Cw is convolved with them. There might be a numerical 
!problem with using zero though.

module ds_mapping
use ds_types
use healpix_types
use ds_utils
implicit none

integer,parameter :: bad_pixel= -1

contains

subroutine make_inv_Cw(noiseinfo,modulescan,cross_correlated) !called per modulescan

    type(ds_noiseinfo) :: noiseinfo
    type(ds_modulescan),intent(inout) :: modulescan
    logical :: cross_correlated
    integer :: i, j, ierror
    logical :: corr
    real(dp),pointer :: sigma

    modulescan%inv_Cw= 0.0_8
    modulescan%has_white_correlations = cross_correlated
    sigma => noiseinfo%sigma
    modulescan%inv_Cw = sigma**2

    if(cross_correlated) then
        write(*,*) "shouldnt be in here cross"
        !make off diagonals
        modulescan%inv_Cw = modulescan%inv_Cw
        !inversion
        ierror= 0 
        stop 'JZ or DWPS sort this out'
        call dpotrf('U',modulescan%inv_Cw,ierror)
        call dpotri('U',modulescan%inv_Cw,ierror)
        modulescan%inv_Cw = modulescan%inv_Cw
    else
        !simple inversion
        modulescan%inv_Cw= 1.0_8 / modulescan%inv_Cw
    endif
end subroutine make_inv_Cw


subroutine invCw_mult_tod(modulescan)
    type(ds_modulescan),intent(inout) :: modulescan
    integer :: i, k
    real(dp) :: vec
    type(ds_timestream),pointer :: timestream
    integer length_check

    timestream => modulescan%timestreams

    if (modulescan%has_white_correlations) then
        write(*,*) "shouldnt be in here correlatios"
        do i=1,timestream%nt
            vec=0
            vec = timestream%timestream(i)
            vec = moduleScan%inv_Cw*vec
            timestream%timestream(i) = vec
        enddo
    else !No cross correlations and so this is simpler
        do i=1,timestream%nt
            timestream%timestream(i) = timestream%timestream(i) * moduleScan%inv_Cw
        enddo
    endif
end subroutine invCw_mult_tod


subroutine invCw_mult_offsets(matinv,modulescan)
    type(ds_modulescan),intent(inout) :: modulescan
    real(dp),intent(in) :: matinv
    integer :: i, k
    real(dp) :: vec
    type(ds_offsets),pointer :: offsets

    offsets => modulescan%offsets
    do i= 1, offsets%na
        vec = offsets%values(i)
        vec = matinv*vec
        offsets%values(i) = vec
    enddo
   
end subroutine invCw_mult_offsets


subroutine add2rhs(modulescan,rhs)   !called once per modulescan  !need an allreduce op after end of use
    !NB: do the C_w^-1 op separately
    type(ds_modulescan),intent(in) :: modulescan !one modulescan per module per scan
    type(ds_trimap) :: rhs
    integer :: t, pixel, i, i_p
    real(dp) :: theta, R_I, R_Q, R_U, S_t, cos_i, sin_i

    if (modulescan%has_p) call ds_assert(rhs%has_p, "no polarization map available in call to add2rhs")
    if (modulescan%has_t) call ds_assert(rhs%has_t, "no temperature map available in call to add2rhs")
    if (rhs%has_t) call ds_assert(allocated(rhs%T%map),"imap not allocated in add2rhs!")
    if (rhs%has_p) call ds_assert(allocated(rhs%Q%map),"qmap not allocated in add2rhs!")
    if (rhs%has_p) call ds_assert(allocated(rhs%U%map),"umap not allocated in add2rhs!")

    if (modulescan%has_leakage_matrix) then
	call ds_assert(modulescan%has_t .and. modulescan%has_p, "Can only use leakage matrix if using full T+P mapping")
        R_I = modulescan%leakage(1)
        R_Q = modulescan%leakage(2)
        R_U = modulescan%leakage(3)
        do t = 1, modulescan%ntod
            pixel = modulescan%pointing(t)
            if (pixel==bad_pixel) cycle
            if (modulescan%flagged%flagged(t) .eq. 0) then
                rhs%T%map(pixel) = rhs%T%map(pixel) + 0.0
                rhs%Q%map(pixel) = rhs%Q%map(pixel) + 0.0
                rhs%U%map(pixel) = rhs%U%map(pixel) + 0.0
            else
                theta = 2.0 * modulescan%theta(t)
                cos_i = cos(theta)
                sin_i = sin(theta)
                S_t = modulescan%timestreams%timestream(t)
                rhs%T%map(pixel) = rhs%T%map(pixel) + R_I * S_t
                rhs%Q%map(pixel) = rhs%Q%map(pixel) + (R_Q*cos_i )*S_t
                rhs%U%map(pixel) = rhs%U%map(pixel) + (R_Q*sin_i )*S_t
            end if 
        enddo
        return 
    endif

    !Temperature Modules
    if (modulescan%has_t) then
        do t = 1, modulescan%ntod
            pixel = modulescan%pointing(t)
            if (pixel==bad_pixel) cycle
            rhs%T%map(pixel) = rhs%T%map(pixel) + modulescan%timestreams%timestream(t)
        enddo
    endif

    !Polarization Modules
    if (modulescan%has_p) then
        do t = 1, modulescan%ntod
            pixel = modulescan%pointing(t)
            if (pixel==bad_pixel) cycle
            rhs%Q%map(pixel) = rhs%Q%map(pixel) + cos(2.0_8 * (modulescan%theta(t))) * modulescan%timestreams%timestream(t)
            rhs%U%map(pixel)= rhs%U%map(pixel) + sin(2.0_8 * (modulescan%theta(t))) * modulescan%timestreams%timestream(t)        
        enddo
    endif
end subroutine add2rhs


subroutine add2cov(modulescan, cov) !called once per module (per code)
    type(ds_modulescan),intent(in) :: modulescan
    type(ds_covariance) :: cov
    integer :: t, p, i, j
    real(dp) :: theta_i, theta_j
    real(dp) :: theta, sin_t, cos_t, cw
    real(dp) :: Ri_I, Ri_Q, Ri_U, Rj_I, Rj_Q, Rj_U
    real(dp) :: Fi_I, Fi_Q, Fi_U, Fj_I, Fj_Q, Fj_U

    if (modulescan%has_t) call ds_assert(cov%has_t, "Naive cov has no temperature in add2cov")
    if (modulescan%has_p) call ds_assert(cov%has_p, "Naive cov has no temperature in add2cov")
    if (modulescan%has_leakage_matrix) then
        call ds_assert(modulescan%has_t .and. modulescan%has_p, "Can only use leakage matrix if using full T+P mapping")
        Ri_I = modulescan%leakage(1)
        Ri_Q = modulescan%leakage(2)
        Ri_U = modulescan%leakage(3)
        !If no noise (between these)/(in this) channel(s) then short-cut
        Rj_I = modulescan%leakage(1)
        Rj_Q = modulescan%leakage(2)
        Rj_U = modulescan%leakage(3)
        do t=1,modulescan%ntod
            p=modulescan%pointing(t)
            if (p==bad_pixel) cycle 
            if (modulescan%flagged%flagged(t) .ne. 0) then
                theta = 2*modulescan%theta(t)			
                cos_t = cos(theta)
                sin_t = sin(theta)
                Fi_I = Ri_I
                Fi_Q = Ri_Q*cos_t 
                Fi_U = Ri_Q*sin_t 
                Fj_I = Rj_I
                Fj_Q = Rj_Q*cos_t 
                Fj_U = Rj_Q*sin_t 
                cw = modulescan%inv_Cw
                cov%TT(p) = cov%TT(p) + Fi_I*cw*Fj_I
                cov%TQ(p) = cov%TQ(p) + Fi_I*cw*Fj_Q
                cov%TU(p) = cov%TU(p) + Fi_I*cw*Fj_U
                cov%QQ(p) = cov%QQ(p) + Fi_Q*cw*Fj_Q
                cov%QU(p) = cov%QU(p) + Fi_Q*cw*Fj_U
                cov%UU(p) = cov%UU(p) + Fi_U*cw*Fj_U
            else 
                cov%TT(p) = cov%TT(p) + 0.0
                cov%TQ(p) = cov%TQ(p) + 0.0
                cov%TU(p) = cov%TU(p) + 0.0
                cov%QQ(p) = cov%QQ(p) + 0.0
                cov%QU(p) = cov%QU(p) + 0.0
                cov%UU(p) = cov%UU(p) + 0.0
           end if  
        enddo
        return
    endif

    if (modulescan%has_t) then
        do t=1,modulescan%ntod
            p=modulescan%pointing(t)
            if (p==bad_pixel) cycle
            cov%TT(p) = cov%TT(p) + modulescan%inv_Cw
        enddo
    endif

    if (modulescan%has_p) then
        do t=1,modulescan%ntod
            p=modulescan%pointing(t)
            if (p==bad_pixel) cycle
            theta_i = 2*(modulescan%theta(t)) 
            theta_j = 2*(modulescan%theta(t))
            
            cov%QQ(p) = cov%QQ(p) + cos(theta_i) * modulescan%inv_Cw * cos(theta_j)
            cov%QU(p) = cov%QU(p) + cos(theta_i) * modulescan%inv_Cw * sin(theta_j)
            cov%UU(p) = cov%UU(p) + sin(theta_i) * modulescan%inv_Cw * sin(theta_j)
        enddo
    endif

    if (modulescan%has_p .and. modulescan%has_t) then
        do t=1,modulescan%ntod
            p=modulescan%pointing(t)
            if (p==bad_pixel) cycle
            theta_j = 2*(modulescan%theta(t))
            cov%TQ(p) = cov%TQ(p) + modulescan%inv_Cw * cos(theta_j)
            cov%TU(p) = cov%TU(p) + modulescan%inv_Cw * sin(theta_j)
        enddo
    endif
end subroutine add2cov


subroutine map2tod(modulescan,maps) !called once per module per projection
    type(ds_modulescan),intent(inout) :: modulescan
	type(ds_trimap) :: maps
    integer :: t, pixel
    real(dp) :: theta_i, cos_i, sin_i, i_t, q_t, u_t
    real(dp) :: R_I, R_Q, R_U

    if (modulescan%has_p) call ds_assert(maps%has_p, "no polarization map available in call to map2tod")
    if (modulescan%has_t) call ds_assert(maps%has_t, "no temperature map available in call to map2tod")

    if (maps%has_t) call ds_assert(allocated(maps%T%map),"imap not allocated in map2modulescan!")
    if (maps%has_p) call ds_assert(allocated(maps%Q%map),"qmap not allocated in map2modulescan!")
    if (maps%has_p) call ds_assert(allocated(maps%U%map),"umap not allocated in map2modulescan!")

    if (modulescan%has_leakage_matrix) then
        call ds_assert(modulescan%has_t .and. modulescan%has_p, "Can only use leakage matrix if using full T+P mapping")
        R_I = modulescan%leakage(1)
        R_Q = modulescan%leakage(2)
        R_U = modulescan%leakage(3)
        do t=1, modulescan%ntod
            pixel = modulescan%pointing(t)
            if(pixel==bad_pixel) cycle
            !if (modulescan%flagged(i)%flagged(t).eq.0.0) cycle
            theta_i = 2*modulescan%theta(t)
            cos_i = cos(theta_i)
            sin_i = sin(theta_i)
            i_t = maps%T%map(pixel)
            q_t = maps%Q%map(pixel)
            u_t = maps%U%map(pixel)
            modulescan%timestreams%timestream(t) = R_I*i_t + R_Q*( cos_i*q_t + sin_i*u_t) 
        enddo
        return 
    endif

    if (modulescan%has_t) then
        do t=1, modulescan%ntod
            pixel = modulescan%pointing(t)
            if(pixel==bad_pixel) cycle	
            modulescan%timestreams%timestream(t) = maps%T%map(pixel)
        enddo
    endif

    if (modulescan%has_p) then
        do t= 1, modulescan%ntod
            pixel = modulescan%pointing(t)
            if(pixel==bad_pixel) cycle	
            modulescan%timestreams%timestream(t) = &
                cos(2.0_8 * (modulescan%theta(t))) * maps%Q%map(pixel) + &
                sin(2.0_8 * (modulescan%theta(t))) * maps%U%map(pixel)
        enddo
    endif
end subroutine map2tod


subroutine cov_mult(maps,cov)
    !cov is only the lower triangle!
    !map has indices and npix
    type(ds_trimap),intent(inout) :: maps
    type(ds_covariance) :: cov
    real(dp) :: t,q,u
    integer :: p, npix
    
    if (maps%has_p) call ds_assert(cov%has_p,"Inconsistent cov and maps in cov_mult (pol)!")
    if (maps%has_t) call ds_assert(cov%has_t,"Inconsistent cov and maps in cov_mult (temp)!")
    if (maps%has_t) call ds_assert(cov%npix==maps%T%npix,"Inconsistent cov and maps in cov_mult (I size)!")
    if (maps%has_p) call ds_assert(cov%npix==maps%Q%npix,"Inconsistent cov and maps in cov_mult (I size)!")
    if (maps%has_p) call ds_assert(cov%npix==maps%U%npix,"Inconsistent cov and maps in cov_mult (I size)!")
    
    if (maps%has_t) then
       npix = maps%T%npix
    else
       npix = maps%Q%npix
    endif
    if (maps%has_t .and. maps%has_p) then
        do p=1,npix
            t = maps%T%map(p)
            q = maps%Q%map(p)
            u = maps%U%map(p)
            maps%T%map(p) = cov%TT(p)*t + cov%TQ(p)*q + cov%TU(p)*u
            maps%Q%map(p) = cov%TQ(p)*t + cov%QQ(p)*q + cov%QU(p)*u
            maps%U%map(p) = cov%TU(p)*t + cov%QU(p)*q + cov%UU(p)*u
        enddo
    else if (maps%has_p) then
        do p=1,npix
            q = maps%Q%map(p)
            u = maps%U%map(p)
            maps%Q%map(p) = cov%QQ(p)*q + cov%QU(p)*u
            maps%U%map(p) = cov%QU(p)*q + cov%UU(p)*u
        enddo		
    else if (maps%has_t) then
        do p=1,npix
            t = maps%T%map(p)
            maps%T%map(p) = cov%TT(p)*t
        enddo
    else
        call ds_assert(.false.,"In cov_mult: neither temperature nor polarization found")
    endif

end subroutine cov_mult


subroutine invert_weight(cov)  !to be used once per code. cov= W^-1
    integer :: p
    type(ds_covariance) :: cov
    real(dp) :: det
    real(dp) :: a11,  a12,  a13
    real(dp) :: a21,  a22,  a23
    real(dp) :: a31,  a32,  a33
    real(dp) :: x1,x2,x3

    if (cov%has_t) call ds_assert(allocated(cov%TT),"TT not allocated in invert_cov!")
    if (cov%has_p) call ds_assert(allocated(cov%QQ),"QQ not allocated in invert_cov!")
    if (cov%has_p) call ds_assert(allocated(cov%QU),"QU not allocated in invert_cov!")
    if (cov%has_p) call ds_assert(allocated(cov%UU),"UU not allocated in invert_cov!")
    if (cov%has_t .and. cov%has_p) call ds_assert(allocated(cov%TQ),"TQ not allocated in invert_cov!")
    if (cov%has_t .and. cov%has_p) call ds_assert(allocated(cov%TU),"TU not allocated in invert_cov!")

    if (cov%has_t .and. cov%has_p) then  !The covariance is a 3x3 matrix.
        do p= 1, cov%npix
            a11 = cov%TT(p)
            a12 = cov%TQ(p)
            a13 = cov%TU(p)
            a21 = a12
            a22 = cov%QQ(p)
            a23 = cov%QU(p)
            a31 = a13
            a32 = a23
            a33 = cov%UU(p)
            x1 = a22*a33 - a23*a32
            x2 = a23*a31 - a21*a33
            x3 = a21*a32 - a22*a31
            det = a11*x1 + a22*x2 + a33*x3
            det = 1.0/det
            
            cov%TT(p) = x1*det
            cov%TQ(p) = x2*det
            cov%TU(p) = x3*det
            cov%QQ(p) = (a11*a33-a13*a31)*det
            cov%QU(p) = (a13*a21-a11*a23)*det
            cov%UU(p) = (a11*a22-a12*a21)*det
        enddo
    else if (cov%has_t) then  !The covariance is a 1x1 matrix
        do p=1,cov%npix
            cov%TT(p) = 1/cov%TT(p)
        enddo
    else if (cov%has_p) then  !The covariance is a 2x2 matrix
        do p=1,cov%npix
            a11 = cov%QQ(p)
            a12 = cov%QU(p)
            a21 = a12
            a22 = cov%UU(p)
            det = a11*a22-a12*a21
            det = 1.0/det
            cov%QQ(p) =  a22*det
            cov%QU(p) = -a12*det
            cov%UU(p) =  a11*det
        enddo
    else
        call ds_assert(.false.,"In invert_weight: neither temperature nor polarization found")	
    endif

end subroutine invert_weight

!I/O routines

subroutine writePixelsFile(cov,filename,nside,maps)
    !writes input file for the PPCL program, containing pixel list and rmsnoise per pixel
    character(len=*) :: filename
    integer :: nside
    type(ds_trimap),intent(in) :: maps
    type(ds_covariance) :: cov
    character(len=50) :: mess
    integer :: i, npix, err

    !open file
    err=0
    open(unit=5004,file=filename,status="replace",iostat=err)
    if(err.ne.0) then
        write(*,*) "Could not open file:",filename
        return
    endif
    !write file header
    write(mess,*) nside
    mess= "nside="//trim(adjustl(mess))
    write(5004,"(a)",iostat=err) trim(adjustl(mess))

    npix= maps%T%npix
    write(mess,*) npix
    mess= "npix="//trim(adjustl(mess))
    write(5004,"(a)") trim(adjustl(mess))
    if(err.ne.0) then
            write(*,*) "Error writing header in file:",filename
    endif
    !write file
    do i= 1,npix
        !Is the rms noise correct - like the magnitude of the complex intensity of polarised noise
        ! rms_noise = sqrt( sigma_qq^2 + sigma_uu^2 )
        write(5004,"('  ',I6,A)") maps%T%indices(i), "NOISE OUTPUT BROKEN" !sqrt(naive_cov(i,1) + naive_cov(i,3)) !naive_cov is a variance
    enddo
    !close file
    close(5004)

end subroutine writePixelsFile


subroutine addMap(source, dest)
    type(ds_map) :: source, dest
    integer i

    call ds_assert(source%npix==dest%npix, "Incompatible maps added")
    do i=lbound(source%map,1), ubound(source%map,1)
        dest%map(i) = dest%map(i) + source%map(i)
    enddo
end subroutine addMap


subroutine subtractMap(subtractee, subtractor)
    !subtractee= subtractee- subtractor
    type(ds_map),intent(inout) :: subtractee
    type(ds_map),intent(in) :: subtractor
    integer i
    
    call ds_assert(subtractee%npix==subtractor%npix, "Incompatible maps subtracted")
    do i=lbound(subtractee%map,1), ubound(subtractee%map,1)
            subtractee%map(i) = subtractee%map(i) - subtractor%map(i)
    enddo
end subroutine subtractMap


subroutine subtractTriMap(a, b)
    !  a -= b
    type(ds_trimap),intent(inout) :: a
    type(ds_trimap),intent(in) :: b
    
    if (a%has_T) then
        call ds_assert(b%has_T, "Tried to subtracted incompatible trimaps (no Temperature in one map).")
        call subtractMap(a%T,b%T)
    endif
    
    if (a%has_P) then
        call ds_assert(b%has_P, "Tried to subtracted incompatible trimaps (no Temperature in one map).")
        call subtractMap(a%Q,b%Q)
        call subtractMap(a%U,b%U)
    endif
end subroutine subtractTriMap


subroutine addOffsetToTimestream(offsets,timestream)
    type(ds_timestream), intent(inout) :: timestream
    type(ds_offsets), intent(inout) :: offsets
    integer a,basisStart, basisEnd

    do a = 1, offsets%na
        basisStart = (a-1)*offsets%length+1
        basisEnd = a*offsets%length
        timestream%timestream(basisStart:basisEnd) = timestream%timestream(basisStart:basisEnd) + offsets%values(a)
    enddo
end subroutine addOffsetToTimestream


subroutine subtractTimestreamFromOffset(timestream,offsets)
    type(ds_timestream), intent(in) :: timestream
    type(ds_offsets), intent(inout) :: offsets
    integer(dp) a,basisStart,basisEnd

    do a = 1,offsets%na
        basisStart = 1+offsets%length*(a-1)
        basisEnd = a*offsets%length
        offsets%values(a) = offsets%values(a) - sum( timestream%timestream(basisStart:basisEnd) )
    enddo
end subroutine


subroutine deprojectTimestreamOntoOffset(timestream,offsets,flagged)
!This is the operator F^T in the PCG.
    type(ds_timestream), intent(in) :: timestream
    type(ds_flagged), intent(in) :: flagged
    type(ds_offsets), intent(inout) :: offsets
    integer(dp) a,i,basisStart,basisEnd
    real(dp) bad_tod
    do a = 1,offsets%na
        basisStart = 1+offsets%length*(a-1)
        basisEnd = a*offsets%length
        bad_tod = 0
        offsets%values(a) = sum(timestream%timestream(basisStart:basisEnd)) - bad_tod 
        !do i=basisStart,basisEnd
        !   if (flagged%flagged(i) .eq. 1) cycle  !if good cycle if bad add. 
        !   bad_tod = bad_tod+timestream%timestream(i)
        !enddo
    enddo
end subroutine 	deprojectTimestreamOntoOffset

end module ds_mapping

