!  Compile on Cartesius with
!    ifort -O2 -o create_in_depths_file_gx1 create_in_depths_file_gx1.f90 -lnetcdf -lnetcdff
!
   program create_in_depths_file_gx1
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Creates a file with depth values, 1st column contains thickness of each
!  layer dz (cm), 2nd column contains depth from surface to midpoint of layer z_t (m)
!  and the 3th column contains depth from surface to bottom of layer z_w_bot (m)
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   use netcdf 
   implicit none



!=====================================================================================================================================================
!
!  variables
!

   ! user input 
   integer          :: km
   character*500    :: depth_file   ! should be a file containing dz, z_t and z_w_bot of POP model

   ! internal model variables 
   integer                                         :: k
   integer                                         :: ncid,varid
   real,             dimension(:),     allocatable :: z_t,z_w_bot,dz

   !netCDF
   include 'netcdf.inc'


!=====================================================================================================================================================
!
!  user input
!

   ! domain info
   write(*,*)' enter km'
   read(*,*)km

   ! name of file containing depth values
   write(*,*)' enter name of file containing depth values'
   read(*,'(a500)')depth_file



!=====================================================================================================================================================
!
!  allocate and read depth values
!

   allocate( z_t(km),z_w_bot(km),dz(km) )
   call check( nf_open(depth_file, NF_NOWRITE, ncid) )
   call check( nf_inq_varid(ncid, "z_t", varid) )
   call check( nf_get_var(ncid, varid, z_t) )
   call check( nf_inq_varid(ncid, "z_w_bot", varid) )
   call check( nf_get_var(ncid, varid, z_w_bot) )
   call check( nf_inq_varid(ncid, "dz", varid) )
   call check( nf_get_var(ncid, varid, dz) )
   call check(nf_close(ncid))
   write(*,*)' done reading file: ',depth_file

   do k = 1, km
     open(3,file='in_depths.dat')
     write(3,'(3f20.3)') dz(k), z_t(k)/100, z_w_bot(k)/100
   enddo
   close(3)

!=====================================================================================================================================================
!
contains


      subroutine check(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
        write(*,*) trim(nf90_strerror(status))
        stop
      endif
 
      end subroutine check



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   end program
