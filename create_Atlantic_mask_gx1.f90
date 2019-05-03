!  Compile on Cartesius with
!    ifort -O2 -o create_Atlantic_mask create_Atlantic_mask.f90 -lnetcdf -lnetcdff
!
   program create_Atlantic_mask
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Creates a mask of the Atlantic based on the j-indices of the southern and northern zonal sections that bound the domain of interest.
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
   character*500    :: kmt_file

   ! internal model variables 
   integer                                         :: rec_length,i,j
   integer                                         :: ncid,idimid,jdimid,i_indexid,j_indexid,maskid,varid
   integer,          dimension(:),     allocatable :: i_index,j_index
   integer,          dimension(:,:),   allocatable :: KMT, atl_mask

   !netCDF
   include 'netcdf.inc'


!=====================================================================================================================================================
!
!  user input
!

   ! domain info
   write(*,*)' enter imt, jmt, km'
   read(*,*)imt,jmt,km

   ! section info
   write(*,*)' enter begin and end latitude indices'
   write(*,*)' j_s, j_e'
   read(*,*)j_s,j_e

   write(*,*)'type 1 if j_e is the Pacific, any other integer if it is in the Atlantic'
   read(*,*)jump

   ! some error checking
   if (j_s<85) then
     write(*,*)' j_s is not allowed to be in Southern Ocean'
     stop
   endif
   if (jump.ne.1) then                             !j_e is in the Atlantic
     if (j_e<=j_s) then
       write(*,*)'  j_e must be larger than j_s'
       stop
     endif
   else                                            !j_e is in the Pacific
     if (j_e<85) then
       write(*,*)' jump=1, j_e is not allowed to be in Southern Ocean'
       stop
     endif
   endif
   if ( (j_s>357).or.(j_e>357) ) then
     write(*,*)' j_s/j_e is not allowed to be in Arctic (ill defined section)'
     stop
   endif

   ! name of kmt file
   write(*,*)' enter name of KMT file'
   read(*,'(a500)')kmt_file



!=====================================================================================================================================================
!
!  allocate and read bathymetry
!

   allocate( KMT(imt,jmt),atl_mask(imt,jmt) )
   call check( nf_open(kmt_file, NF_NOWRITE, ncid) )
   call check( nf_inq_varid(ncid, "KMT", varid) )
   call check( nf_get_var(ncid, varid, KMT) )
   call check(nf_close(ncid))
   write(*,*)' done reading file: ',kmt_file



!=====================================================================================================================================================
!
!  create Atlantic mask
!

   atl_mask = 0
   if (jump.ne.1) then                             !j_e is in the Atlantic
     ! create the seed for the Atlantic
     i = 10
     j = j_s
     call flood(i,j,KMT,imt,jmt,j_s,j_e,atl_mask)
   else                                            !j_e is in the Pacific
     ! create the seed for the Pacific (may not be fully adequate)
     i = 200
     j = j_e
     call flood(i,j,KMT,imt,jmt,j_e,357,atl_mask)
     where (atl_mask==1) atl_mask = -1             !set Pacific mask to -1
     i = 10                                        !then find Atlantic proper
     j = 357
     call flood(i,j,KMT,imt,jmt,j_s,jmt,atl_mask)
   endif


!=====================================================================================================================================================
!
!  write mask to NETCDF
!

   ! create file
   call check( nf_create('Atlantic_mask_gx1.nc', nf_write, ncid) )

   ! define dimensions
   call check( nf_def_dim(ncid,'i_index',imt,idimid) )
   call check( nf_def_dim(ncid,'j_index',jmt,jdimid) )

   ! define variables
   call check( nf_def_var(ncid,'i_index' ,nf_int,  1,idimid,                  i_indexid) )
   call check( nf_def_var(ncid,'j_index' ,nf_int,  1,jdimid,                  j_indexid) )
   call check( nf_def_var(ncid,'ATL_MASK',nf_int,  2,(/idimid,jdimid/),       maskid   ) )
   call check( nf_enddef(ncid) )

   ! put variables
   allocate( i_index(imt),j_index(jmt) )
   i_index = (/(i,i=1,imt)/)
   j_index = (/(j,j=1,jmt)/)
   call check( nf_put_var_int( ncid,i_indexid,i_index ) )
   call check( nf_put_var_int( ncid,j_indexid,j_index ) )
   call check( nf_put_var_int( ncid,maskid   ,atl_mask) )
   call check( nf_close(ncid) )



   contains
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



      recursive subroutine flood(i,j,KMT,imt,jmt,j_s,j_e,atl_mask)
!
!     'flood' the area bounded by j_s and j_e
!

      ! input/output variables
      integer,                     intent(in)    :: i,j,imt,jmt,j_s,j_e
      integer, dimension(imt,jmt), intent(in)    :: KMT 
      integer, dimension(imt,jmt), intent(inout) :: atl_mask

      ! local variable
      integer :: xx 

      ! flooding
      if( (KMT(i,j)>0).AND.(atl_mask(i,j)==0) ) then
        atl_mask(i,j) = 1

        xx = i+1
        if (xx>imt) xx = 1
        call flood(xx,j,KMT,imt,jmt,j_s,j_e,atl_mask)

        xx = i-1
        if (xx<1  ) xx = imt
        call flood(xx,j,KMT,imt,jmt,j_s,j_e,atl_mask)

        xx = j+1
        if (xx>j_e) xx = j_e 
        call flood(i,xx,KMT,imt,jmt,j_s,j_e,atl_mask)

        xx = j-1
        if (xx<j_s) xx = j_s
        call flood(i,xx,KMT,imt,jmt,j_s,j_e,atl_mask)

      endif

      end subroutine flood



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



      subroutine check(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
        write(*,*) trim(nf90_strerror(status))
        stop
      endif
 
      end subroutine check



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   end program
