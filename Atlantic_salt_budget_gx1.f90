!  Compile on Cartesius with
!    ifort -O2 -convert big_endian -o Atlantic_salt_budget_gx1 Atlantic_salt_budget_gx1.f90 -lnetcdf -lnetcdff
!
   program Atlantic_salt_budget
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Calculates the salt budget for the Atlantic area determined by create_Atlantic_mask. It computes the salt transport across z=0, and across the
!  two bounding zonal sections, separating the latter into overturning (M_ov) and azonal (M_az) contributions according to
!  De Vries and Weber (GRL, 2005), Dijkstra (TellusA, 2007), Huisman et al. (JPO, 2010), Drijfhout et al. (Clim. Dyn, 2010)
!
!  Michael Kliphuis & Matthijs den Toom (IMAU, August 2011)
! 
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   use netcdf 
   implicit none



!=====================================================================================================================================================
!
!  variables
!

   ! user input 
   integer          :: imt,jmt,km
   integer          :: j_s,j_e,jump
   integer          :: nrec_VVEL,nrec_SALT,nrec_SSH,nrec_SFWF
   character*500    :: salt_file,vvel_file,ssh_file,sfwf_file,in_depths

   ! program output
   double precision :: area,volume,SSH_int,SALT_SSH_int,EPR,SALT_int
   double precision :: sS0,sT,sM_ov,sM_az,sM_tot
   double precision :: eS0,eT,eM_ov,eM_az,eM_tot
   real             :: missing_value_salt, missing_value_vvel, missing_value_ssh, missing_value_sfwf

   ! internal model variables 
   integer                                           :: rec_length,i,ip,j,k
   integer                                           :: ncid,maskid,varid
   integer,          dimension(:,:),     allocatable :: KMT, atl_mask
   
   real,             dimension(:,:),     allocatable :: VVEL_k,SALT_k ! read this as real !
   real,             dimension(:,:,:),   allocatable :: SALT3D
   real,             dimension(:,:,:),   allocatable :: VVEL3D
   real,             dimension(:,:),     allocatable :: SSH, SFWF
   
   double precision, dimension(:),       allocatable :: z1,z2,dz
   double precision, dimension(:,:),     allocatable :: HTN,HTE,WORK,DXT,DYT,TAREA,DXU,DYU
   double precision, dimension(:,:),     allocatable :: sVVEL,sSALT,eVVEL,eSALT,secDXDZU
   double precision, dimension(:,:,:),   allocatable :: DZT,SALT

   ! parameters
   double precision, parameter :: c0 = 0., p5 = 0.5, c1 = 1., S0ref = 0.035

   !netCDF
   include 'netcdf.inc'


!=====================================================================================================================================================
!
!  user input
!

   ! domain info
   write(*,*)' enter imt, jmt, km'
   read(*,*)imt,jmt,km

   ! section info (error checking is done in create_Atlantic_mask)
   write(*,*)' enter begin and end latitude indices'
   write(*,*)' j_s, j_e'
   read(*,*)j_s,j_e

   write(*,*)'type 1 if j_e is the Pacific, any other integer if it is in the Atlantic'
   read(*,*)jump

   ! names of files, etc
   write(*,*)' enter name of SALT file'
   read(*,'(a500)')salt_file
   write(*,*)' enter name of VVEL file'
   read(*,'(a500)')vvel_file
   write(*,*)' enter name of SSH file'
   read(*,'(a500)')ssh_file
   write(*,*)' enter name of SFWF file'
   read(*,'(a500)')sfwf_file
   write(*,*)' enter name of in_depths file'
   read(*,'(a500)')in_depths



!=====================================================================================================================================================
!
!  allocate arrays
!

   ! 2-D
   allocate( atl_mask(imt,jmt), SSH(imt,jmt),SFWF(imt,jmt), TAREA(imt,jmt),DXU(imt,jmt),DYU(imt,jmt) )

   ! 2-D (sections)
   allocate( sVVEL(imt,km),sSALT(imt,km),eVVEL(imt,km),eSALT(imt,km),secDXDZU(imt,km) )

   ! 3-D
   allocate( DZT(imt,jmt,km),SALT(imt,jmt,km) )

!=====================================================================================================================================================
!
!  read and create horizontal grid spacing, define TAREA
!  read SALT, VVEL, KMT, SSH and SFWF

   allocate( HTN(imt,jmt),HTE(imt,jmt),WORK(imt,jmt),DXT(imt,jmt),DYT(imt,jmt) )
   allocate( KMT(imt,jmt) )
   allocate( SALT3D(imt, jmt, km) )
   allocate( VVEL3D(imt, jmt, km)  )

   ! open salt_file 
   ! this file should contain variables SALT, HTN, HTE, KMT, z_t, z_w_bot and dz)
   call check( nf_open(salt_file, NF_NOWRITE, ncid) )
   call check( nf_inq_varid(ncid, "SALT", varid) )
   call check( nf_get_var(ncid, varid, SALT3D) )
   ! make sure you later set missing values to 0 (like in the original binary POP data)
   call check( nf_get_att(ncid, varid, "missing_value", missing_value_salt) )
   write(6,*) 'missing_value_salt is: ',missing_value_salt
   call check( nf_inq_varid(ncid, "HTN", varid) )
   call check( nf_get_var(ncid, varid, HTN) )
   call check( nf_inq_varid(ncid, "HTE", varid) )
   call check( nf_get_var(ncid, varid, HTE) )
   call check( nf_inq_varid(ncid, "KMT", varid) )
   call check( nf_get_var(ncid, varid, KMT) )
   call check(nf_close(ncid))
   
   write(*,*)' done reading file: ',salt_file
   
   ! open vvel_file 
   call check( nf_open(vvel_file, NF_NOWRITE, ncid) )
   call check( nf_inq_varid(ncid, "VVEL", varid) )
   call check( nf_get_var(ncid, varid, VVEL3D) )
   ! make sure you later set missing values to 0 (like in the original binary POP data)
   call check( nf_get_att(ncid, varid, "missing_value", missing_value_vvel) )
   write(6,*) 'missing_value_vvel is: ',missing_value_vvel
   call check(nf_close(ncid))


   write(*,*)' done reading file: ',vvel_file
   
   ! open ssh_file 
   call check( nf_open(ssh_file, NF_NOWRITE, ncid) )
   call check( nf_inq_varid(ncid, "SSH", varid) )
   call check( nf_get_var(ncid, varid, SSH) )
   ! make sure you later set missing values to 0 (like in the original binary POP data)
   call check( nf_get_att(ncid, varid, "missing_value", missing_value_ssh) )
   write(6,*) 'missing_value_ssh is: ',missing_value_ssh
   call check(nf_close(ncid))

   write(*,*)' done reading file: ',ssh_file

   ! open sfwf_file 
   call check( nf_open(sfwf_file, NF_NOWRITE, ncid) )
   call check( nf_inq_varid(ncid, "SFWF", varid) )
   call check( nf_get_var(ncid, varid, SFWF) )
   ! make sure you later set missing values to 0 (like in the original binary POP data)
   call check( nf_get_att(ncid, varid, "missing_value", missing_value_sfwf) )
   write(6,*) 'missing_value_sfwf is: ',missing_value_sfwf
   call check(nf_close(ncid))

   write(*,*)' done reading file: ',sfwf_file
  
   ! in binary case POP_gx1 all missing values are 0
   ! set them to 0 here too 
   where (SALT3D .eq. missing_value_salt ) 
     SALT3D = c0
   endwhere
   
   where (VVEL3D .eq. missing_value_vvel ) 
     VVEL3D = c0
   endwhere
   
   where (SSH .eq. missing_value_ssh) 
     SSH = c0
   endwhere
   
   where (SFWF .eq. missing_value_sfwf) 
     SFWF = c0
   endwhere
   
   ! Why do we do this?
   where (HTN <= c0) HTN = c1
   where (HTE <= c0) HTE = c1

   call s_rshift(WORK,HTN,imt,jmt)
   DXT   = p5*(HTN + WORK)
   call w_rshift(WORK,HTE,imt,jmt)
   DYT   = p5*(HTE + WORK)
   TAREA = DXT*DYT

   call e_rshift(WORK,HTN,imt,jmt)
   DXU = p5*(HTN + WORK)
   call n_rshift(WORK,HTE,imt,jmt)
   DYU = p5*(HTE + WORK)

   deallocate( HTN,HTE,WORK,DXT,DYT )


!=====================================================================================================================================================
!
!  read depth level variables, define DZT
!

   allocate( z1(km),z2(km),dz(km) )

   ! dz
   open(1,file=in_depths,status='old')
   do k = 1, km
     read(1,*) dz(k),z1(k),z2(k)
   enddo
   close(1)
   write(*,*)' done reading file: ',in_depths

   ! partial bottom cell depths

   !DZT
   do k=1,km
     where(KMT>=k)
       DZT(:,:,k) = dz(k)
     elsewhere
       DZT(:,:,k) = c0
     endwhere
   enddo

   deallocate( z1,z2,dz,KMT )



!=====================================================================================================================================================
!
!  read mask from NETCDF
!

   call check( nf_open('Atlantic_mask_gx1.nc', nf_nowrite, ncid) )
   call check( nf_inq_varid(ncid,'ATL_MASK',maskid) )
   call check( nf_get_var_int(ncid,maskid,atl_mask) )
   call check( nf_close(ncid) )



!=====================================================================================================================================================
!
!  read velocity, salinity, sea surface height, and virtual salt flux
!

   allocate( VVEL_k(imt,jmt),SALT_k(imt,jmt) )

   !then read 3-D fields, storing full salinity field, and velocity field for sections
   do k = 1, km
     SALT_k = SALT3D(:,:,k) / 1000 ! should be in order of 0.035
     VVEL_k = VVEL3D(:,:,k)

     ! sections of meridional velocity and salinity (interpolated to latitude of u points)
     sVVEL(:,k) = VVEL_k(:,j_s-1)
     sSALT(:,k) = p5*( SALT_k(:,j_s-1) + SALT_k(:,j_s  ) )
     if (jump.ne.1) then   ! j_e is in Atlantic
       eVVEL(:,k) = VVEL_k(:,j_e  )
       eSALT(:,k) = p5*( SALT_k(:,j_e  ) + SALT_k(:,j_e+1) )
     else                  ! j_e is in Pacific
       eVVEL(:,k) = VVEL_k(:,j_e-1)
       eSALT(:,k) = p5*( SALT_k(:,j_e-1) + SALT_k(:,j_e  ) )
     endif
     ! salinity field
     SALT(:,:,k) = SALT_k
   enddo

   deallocate( VVEL_k,SALT_k )


!=====================================================================================================================================================
!
!  calculate averaged quantities
!

   ! fixed in time
   area          = sum(TAREA,                 atl_mask.ne.0)/1.0E04          ! surface area
   volume        = sum(TAREA*sum(DZT,3),      atl_mask.ne.0)/1.0E06          ! volume between bottom and z=0

   ! variable, 2-D
   SSH_int       = sum(TAREA*SSH,             atl_mask.ne.0)/1.0E06          ! volume above z=0, due to SSH
   SALT_SSH_int  = sum(TAREA*SALT(:,:,1)*SSH, atl_mask.ne.0)/(S0ref*1.0E06)  ! salt content above z=0
   EPR           = sum(TAREA*SFWF,            atl_mask.ne.0)/(1.0E13)        ! integral of freshwater flux

   ! variable, 3-D
   SALT_int      = sum(TAREA*sum(DZT*SALT,3), atl_mask.ne.0)/(S0ref*1.0E06)  ! salt content below z=0 (this should be sufficiently accurate; if not
                                                                             !   subtract fixed volume before summations)


!=====================================================================================================================================================
!
!  calculate transports for section
!

   ! create integral operator for section with j_s
   secDXDZU = c0
   do k = 1, km
   do i = 1, imt
     ip = i+1
     if (ip>imt) ip = 1
     if ( (atl_mask(i,j_s)==1).AND.(atl_mask(ip,j_s)==1) ) then
       secDXDZU(i,k) = DXU(i,j_s-1)*min( DZT(i,j_s-1,k),DZT(ip,j_s-1,k),DZT(i,j_s,k),DZT(ip,j_s,k) )
     endif
   enddo
   enddo
   call section_transport(imt,km,secDXDZU,sVVEL,sSALT,sS0,sT,sM_ov,sM_az,sM_tot)

   ! create integral operator for section with j_e
   secDXDZU = c0
   do k = 1, km
   do i = 1, imt
     ip = i+1
     if (ip>imt) ip = 1
     if ( (jump.ne.1).AND.(atl_mask(i,j_e)== 1).AND.(atl_mask(ip,j_e)== 1) ) then       ! j_e in Atlantic
       secDXDZU(i,k) = DXU(i,j_e  )*min( DZT(i,j_e  ,k),DZT(ip,j_e  ,k),DZT(i,j_e+1,k),DZT(ip,j_e+1,k) )
     endif
     if ( (jump == 1).AND.(atl_mask(i,j_e)==-1).AND.(atl_mask(ip,j_e)==-1) ) then       ! j_e in Pacific
       secDXDZU(i,k) = DXU(i,j_e-1)*min( DZT(i,j_e-1,k),DZT(ip,j_e-1,k),DZT(i,j_e  ,k),DZT(ip,j_e  ,k) )
     endif
   enddo
   enddo
   if (jump.ne.1) then ! j_e in Atlantic
     call section_transport(imt,km,secDXDZU, eVVEL,eSALT,eS0,eT,eM_ov,eM_az,eM_tot)
   else                ! j_e in Pacific
     call section_transport(imt,km,secDXDZU,-eVVEL,eSALT,eS0,eT,eM_ov,eM_az,eM_tot)
   endif



!=====================================================================================================================================================
!
!  output
!

   open(2,file='Atlantic_salt_budget_value')
   write(2,'(16ES18.10)') area,volume,SSH_int,sT,eT,SALT_int,SALT_SSH_int,sS0*1000.,sM_ov,sM_az,sM_tot,eS0*1000,eM_ov,eM_az,eM_tot,EPR
!                         1    2      3       4  5  6        7            8         9     10    11     12       13    14    15     16
   close(2)



   contains
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



      subroutine check(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
        write(*,*) trim(nf90_strerror(status))
        stop
      endif
 
      end subroutine check


 
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


 
      subroutine section_transport(imt,km,DXDZU,VVEL,SALT,S0,T,M_ov,M_az,M_tot)
!
!     section tranports, separated into M_ov, M_az
!
      implicit none

      ! input/output variables
      integer,                             intent(in)  :: imt,km
      double precision, dimension(imt,km), intent(in)  :: DXDZU,VVEL,SALT
      double precision,                    intent(out) :: S0,T,M_ov,M_az,M_tot

      ! local variables
      double precision                    :: VVEL_tilde,check
      double precision, dimension(    km) :: DXDZU_int,VVEL_int,SALT_int,VVEL_x_SALT_int,VVEL_prime_x_SALT_prime_int
      double precision, dimension(    km) :: VVEL_brackets,SALT_brackets,VVEL_star
      double precision, dimension(imt,km) :: ONE,VVEL_prime,SALT_prime

      ! parameters
      double precision, parameter :: c0 = 0., toSv = 1.E-012, S0ref = 0.035

      !auxililary array
      ONE = 1.

      ! calculate transports
      call zonalintegral(imt,km,DXDZU, ONE, ONE,  DXDZU_int      )
      call zonalintegral(imt,km,DXDZU, VVEL,ONE,  VVEL_int       )
      call zonalintegral(imt,km,DXDZU, ONE, SALT, SALT_int       )
      call zonalintegral(imt,km,DXDZU, VVEL,SALT, VVEL_x_SALT_int)

      VVEL_brackets = c0
      SALT_brackets = c0
      where (DXDZU_int > 0.0 ) 
        VVEL_brackets = VVEL_int/DXDZU_int
        SALT_brackets = SALT_int/DXDZU_int
      endwhere
      do k = 1, km
        VVEL_prime(:,k) = VVEL(:,k) - VVEL_brackets(k)
        SALT_prime(:,k) = SALT(:,k) - SALT_brackets(k)
      enddo
      call zonalintegral(imt,km,DXDZU, VVEL_prime, SALT_prime, VVEL_prime_x_SALT_prime_int)

      VVEL_tilde  = sum(VVEL_int)/sum(DXDZU_int) 
      VVEL_star   = VVEL_brackets - VVEL_tilde
      S0          = sum(SALT_int)/sum(DXDZU_int)
      T     =  toSv*sum(VVEL_int)
      M_ov  = -toSv*sum(VVEL_star*SALT_brackets*DXDZU_int)/S0ref
      M_az  = -toSv*sum(VVEL_prime_x_SALT_prime_int)/S0ref
      M_tot = -toSv*sum(VVEL_x_SALT_int)/S0ref

      ! check
      check = -T+(M_ov+M_az-M_tot)*S0ref/S0
      write(6,*)'total correct?', check

      end subroutine section_transport



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



      subroutine zonalintegral(imt,km,A,inV,inS,output)
!
!     Zonal integral of inV times inS over section. Array inV is defined on velocity points, array inS on midpoints between tracer points.
!     Array A is the water filled-area in longitude depth plane of velocity cells.
!
      implicit none

      ! input/output variables
      integer,                             intent(in)  :: imt, km
      double precision, dimension(imt,km), intent(in)  :: A, inV
      double precision, dimension(imt,km), intent(in)  :: inS
      double precision, dimension(    km), intent(out) :: output

      ! local variable
      integer :: k

      ! computation of integral
      do k = 1, km
        output(k) =                  0.5*( A(imt    ,k)*inV(imt,    k) + A(1,    k)*inV(1,    k) ) * inS(1,    k)
        output(k) = output(k) + sum( 0.5*( A(1:imt-1,k)*inV(1:imt-1,k) + A(2:imt,k)*inV(2:imt,k) ) * inS(2:imt,k) )
      enddo

      end subroutine zonalintegral



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   end program
!
!  Below are some useful (model-native?) subroutines
!

!***********************************************************************

      subroutine e_rshift(XOUT,X,imt,jmt)

      implicit none

!-----------------------------------------------------------------------
!
!     shift from east (double precision arrays)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!

!     inputs
!
!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(in) ::  &
         X                  ! array to be shifted

!-----------------------------------------------------------------------
!
!     outputs
!
!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(out) ::  &
         XOUT               ! shifted result

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer ::  &
         i,j,n,ip &                     ! dummy counters
        ,imt,jmt

!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          ip = i + 1
          if(i == imt) ip = 1
          XOUT(i,j) = X(ip,j)
        end do
      end do

!-----------------------------------------------------------------------

      end subroutine e_rshift

!***********************************************************************

      subroutine s_rshift(XOUT,X,imt,jmt)

      implicit none

!-----------------------------------------------------------------------
!
!     shift from south (real arrays)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     inputs
!
!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(in) ::  &
         X                  ! array to be shifted

!-----------------------------------------------------------------------
!
!     outputs
!
!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(out) ::  &
         XOUT               ! shifted result

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer ::  &
         i,j,n,jm &                     ! dummy counters
        ,imt,jmt

!-----------------------------------------------------------------------

      do j=1,jmt
        jm = j - 1
        if(j == 1) jm = 1
        do i=1,imt
          XOUT(i,j) = X(i,jm)
        end do
      end do

!-----------------------------------------------------------------------

      end subroutine s_rshift

!***********************************************************************

      subroutine w_rshift(XOUT,X,imt,jmt)

      implicit none

!-----------------------------------------------------------------------
!
!     shift from west (real arrays)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     inputs
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(in) ::  &
         X                  ! array to be shifted

!-----------------------------------------------------------------------
!
!     outputs
!
!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(out) ::  &
         XOUT               ! shifted result

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer ::  &
         i,j,n,im &                     ! dummy counters
        ,imt,jmt

!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          im = i - 1
          if(i == 1) im = imt
          XOUT(i,j) = X(im,j)
        end do
      end do

!-----------------------------------------------------------------------

      end subroutine w_rshift

!***********************************************************************

      subroutine n_rshift(XOUT,X,imt,jmt)

      implicit none

!-----------------------------------------------------------------------
!
!     shift from north (double precision arrays)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     inputs
!
!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(in) ::  &
         X                  ! array to be shifted

!-----------------------------------------------------------------------
!
!     outputs
!
!-----------------------------------------------------------------------

      double precision, dimension(imt,jmt), intent(out) ::  &
         XOUT               ! shifted result

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer ::  &
         i,j,n,jp &                     ! dummy counters
        ,imt,jmt

!-----------------------------------------------------------------------

      do j=1,jmt
        jp = j+1
        if(j == jmt) jp = jmt
        do i=1,imt
          XOUT(i,j) = X(i,jp)
        end do
      end do

!-----------------------------------------------------------------------

      end subroutine n_rshift

!***********************************************************************
