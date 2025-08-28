MODULE isfutils
   !!======================================================================
   !!                       ***  MODULE  isfutils  ***
   !! istutils module : miscelenious useful routines
   !!======================================================================
   !! History :  4.1  !  2019-09  (P. Mathiot) original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isfutils       : - read_2dcstdta to read a constant input file with iom_get
   !!                    - debug to print array sum, min, max in ocean.output
   !!----------------------------------------------------------------------

   USE iom           , ONLY: iom_open, iom_get, iom_close, jpdom_global      ! read input file
   USE lib_fortran
   USE par_oce                                                               ! domain size
   USE dom_oce       , ONLY: narea                                           ! local domain
   USE in_out_manager                                                        ! miscelenious
   USE lib_mpp
   USE domutl, ONLY : is_tile

   IMPLICIT NONE

   PRIVATE

   INTERFACE debug
      MODULE PROCEDURE debug2d, debug3d
   END INTERFACE debug

   PUBLIC read_2dcstdta, debug

   !! * Substitutions
#  include "do_loop_substitute.h90"
CONTAINS

   SUBROUTINE read_2dcstdta(cdfile, cdvar, pvar)
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE read_2dcstdta  ***
      !!
      !! ** Purpose : read input file
      !!
      !!--------------------------------------------------------------------
      CHARACTER(len=*)         , INTENT(in   ) ::   cdfile   ! input file name
      CHARACTER(len=*)         , INTENT(in   ) ::   cdvar    ! variable name
      REAL(wp), DIMENSION(:,:) , INTENT(inout) ::   pvar     ! output variable
      !!--------------------------------------------------------------------
      INTEGER :: inum
      !!--------------------------------------------------------------------

      CALL iom_open( TRIM(cdfile), inum )
      CALL iom_get( inum, jpdom_global, TRIM(cdvar), pvar )
      CALL iom_close(inum)

   END SUBROUTINE read_2dcstdta

   SUBROUTINE debug2d(cdtxt,pvar)
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE isf_debug2d  ***
      !!
      !! ** Purpose : add debug print for 2d variables
      !!
      !!--------------------------------------------------------------------
      CHARACTER(LEN=*)        , INTENT(in) ::   cdtxt
      REAL(wp), DIMENSION(:,:), INTENT(in) ::   pvar
      !!--------------------------------------------------------------------
      REAL(wp)    ::   zmin, zmax, zsum
      INTEGER(i8) ::   imodd, ip
      INTEGER     ::   imods
      INTEGER     ::   isums, idums
      INTEGER     ::   ji, jj, jk
      INTEGER, DIMENSION(jpnij) ::   itmps

      INTEGER    ::   ipi, ipj, ipk, ipl        ! dimensions
      INTEGER    ::   iilsht, ijlsht            ! loop shift indices
      INTEGER    ::   iiasht, ijasht            ! array shift indices
      !!--------------------------------------------------------------------
      !
      ! global min/max/sum to check data range and NaN
      zsum = glob_2Dsum( 'debug', pvar(:,:) )
      zmin = glob_2Dmin( 'debug', pvar(:,:) )
      zmax = glob_2Dmax( 'debug', pvar(:,:) )
      !
      ! basic check sum to check reproducibility
      ! TRANSFER function find out the integer corresponding to pvar(i,j) bit pattern
      ! MOD allow us to keep only the latest digits during the sum
      ! imod is not choosen to be very large as at the end there is a classic mpp_sum
      imodd=65521 ! highest prime number < 2**16 with i8 type
      imods=65521 ! highest prime number < 2**16 with default integer for mpp_sum subroutine
      isums=0 ; itmps(:)=0 ;
      !
      ! local MOD sum
      !!-----------------------------------------------------------------------
      !
      ipi = SIZE(pvar,1)   ! 1st dimension
      ipj = SIZE(pvar,2)   ! 2nd dimension
      !
      IF( .NOT. is_tile(SIZE(pvar,1), SIZE(pvar,2)) ) THEN
         iilsht = ( jpi - ipi ) / 2
         ijlsht = ( jpj - ipj ) / 2   ! should be the same as iisht...
      ELSE ! Tile sized array
         iilsht = ( ntei - ntsi + 1 - ipi ) / 2 + nn_hls
         ijlsht = ( ntej - ntsj + 1 - ipj ) / 2 + nn_hls
      END IF
      iiasht = iilsht + ntsi - 1 - nn_hls
      ijasht = ijlsht + ntsj - 1 - nn_hls
      !
      DO_2D( 0, 0, 0, 0 )
         idums = ABS(MOD(TRANSFER(pvar(ji-iiasht,jj-ijasht), ip),imodd))
         itmps(narea) = MOD(itmps(narea) + idums, imods)
      END_2D
      !
      ! global MOD sum
      CALL mpp_max('debug',itmps(:))
      DO jk = 1,jpnij
         isums = MOD(isums + itmps(jk),imods)
      END DO
      !
      ! print out
      IF (lwp) THEN
         WRITE(numout,*) TRIM(cdtxt),' (min, max, sum, tag) : ',zmin, zmax, zsum, isums
         CALL FLUSH(numout)
      END IF
      !
   END SUBROUTINE debug2d

   SUBROUTINE debug3d(cdtxt,pvar)
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE isf_debug3d  ***
      !!
      !! ** Purpose : add debug print for 3d variables
      !!
      !!--------------------------------------------------------------------
      CHARACTER(LEN=*)               , INTENT(in) ::   cdtxt
      REAL(wp), DIMENSION(A2D(0),jpk), INTENT(in) ::   pvar
      !!--------------------------------------------------------------------
      REAL(wp)    ::   zmin, zmax, zsum
      INTEGER(i8) ::   imodd, ip
      INTEGER     ::   imods
      INTEGER     ::   isums, idums
      INTEGER     ::   ji, jj, jk
      INTEGER, DIMENSION(jpnij) ::   itmps
      !!--------------------------------------------------------------------
      !
      ! global min/max/sum to check data range and NaN
      zsum = glob_3Dsum( 'debug', pvar(:,:,:) )
      zmin = glob_3Dmin( 'debug', pvar(:,:,:) )
      zmax = glob_3Dmax( 'debug', pvar(:,:,:) )
      !
      ! basic check sum to check reproducibility
      ! TRANSFER function find out the integer corresponding to pvar(i,j) bit pattern
      ! MOD allow us to keep only the latest digits during the sum
      ! imod is not choosen to be very large as at the end there is a classic mpp_sum
      imodd=65521 ! highest prime number < 2**16 with i8 type
      imods=65521 ! highest prime number < 2**16 with default integer for mpp_sum subroutine
      itmps=0; isums=0
      !
      ! local MOD sum
      DO jk=1,jpk
         DO jj=Njs0,Nje0
            DO ji=Nis0,Nie0
               idums = ABS(MOD(TRANSFER(pvar(ji,jj,jk), ip),imodd))
               itmps(narea) = MOD(itmps(narea) + idums, imods)
            END DO
         END DO
      END DO
      !
      ! global MOD sum
      CALL mpp_max('debug',itmps)
      DO jk = 1,jpnij
         isums = MOD(isums+itmps(jk),imods)
      END DO
      !
      ! print out
      IF (lwp) THEN
         WRITE(numout,*) TRIM(cdtxt),' (min, max, sum, tag) : ',zmin, zmax, zsum, isums
         CALL FLUSH(numout)
      END IF
      !
   END SUBROUTINE debug3d

END MODULE isfutils
