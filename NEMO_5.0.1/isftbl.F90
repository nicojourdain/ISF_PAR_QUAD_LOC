MODULE isftbl
   !!======================================================================
   !!                       ***  MODULE  isftbl  ***
   !! isftbl module :  compute properties of top boundary layer
   !!======================================================================
   !! History :  4.1  !  2019-09  (P. Mathiot) original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isftbl       : routine to compute : 
   !!                  - geometry of the ice shelf tbl (isf_tbl_lvl, isftbl_ktop, isftbl_kbot)
   !!                    (top and bottom level, thickness and fraction of deepest level affected)
   !!                  - tbl averaged properties (isf_tbl, isf_tbl_avg)
   !!----------------------------------------------------------------------

   USE isf_oce ! ice shelf variables

   USE par_oce ! ocean space and time domain
   USE dom_oce ! vertical scale factor and depth
   USE domutl, ONLY : lbnd_ij
   USE lib_fortran, ONLY: glob_2Dsum

   IMPLICIT NONE

   PRIVATE

   PUBLIC isf_tbl_avg, isf_tbl_lvl, isf_tbl_ktop
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
CONTAINS

!!$   SUBROUTINE isf_tbl( Kmm, pvarin, cd_ptin, ktop, phtbl, pvarout, kbot, pfrac )
!!$      !!--------------------------------------------------------------------
!!$      !!                  ***  SUBROUTINE isf_tbl  ***
!!$      !!
!!$      !! ** Purpose : compute mean T/S/U/V in the boundary layer at T- point
!!$      !!
!!$      !! ** Method : Average properties over a specific thickness
!!$      !!
!!$      !! ** Reference : inspired from : Losch, Modeling ice shelf cavities in a z coordinate ocean general circulation model
!!$      !!                https://doi.org/10.1029/2007JC004368 , 2008
!!$      !!
!!$      !!--------------------------------------------------------------------
!!$      INTEGER                              , INTENT(in ) ::   Kmm           ! ocean time level index
!!$      REAL(wp), DIMENSION(A2D(0),jpk)      , INTENT(in ) ::   pvarin        ! 3d variable to average over the tbl
!!$      CHARACTER(len=1)                     , INTENT(in ) ::   cd_ptin       ! point of variable in/out
!!$      INTEGER,  DIMENSION(A2D(0))          , INTENT(in ) ::   ktop          ! top level
!!$      REAL(wp), DIMENSION(A2D(0))          , INTENT(in ) ::   phtbl         ! tbl thickness
!!$      REAL(wp), DIMENSION(A2D(0))          , INTENT(out) ::   pvarout       ! 2d average of pvarin
!!$      INTEGER,  DIMENSION(A2D(0)), OPTIONAL, INTENT(in ) ::   kbot          ! bottom level
!!$      REAL(wp), DIMENSION(A2D(0)), OPTIONAL, INTENT(in ) ::   pfrac         ! fraction of bottom cell affected by tbl
!!$      !!--------------------------------------------------------------------
!!$      INTEGER ::   ji, jj, jk, ikt             ! loop index
!!$      REAL(wp), DIMENSION(A2D(0)) ::   zhtbl   ! temporary array for thickness
!!$      INTEGER , DIMENSION(A2D(0)) ::   ikbot   ! bottom level of the tbl
!!$      REAL(wp), DIMENSION(A2D(0)) ::   zfrac   ! thickness of the tbl
!!$      REAL(wp), DIMENSION(A2D(0),jpk) ::   ze3 ! e3
!!$      !!--------------------------------------------------------------------
!!$      !
!!$      SELECT CASE ( cd_ptin )
!!$      CASE ( 'U' )
!!$         !
!!$         DO_3D( 0 ,0, 0, 0, 1, jpk )
!!$            ze3(ji,jj,jk) = e3u(ji,jj,jk,Kmm)
!!$         END_3D
!!$         !
!!$         ! compute tbl lvl and thickness
!!$         DO_2D( 0 ,0, 0, 0 )
!!$            ikt = ktop(ji,jj)  ! tbl top indices
!!$            zhtbl(ji,jj) = MAX( MIN( phtbl(ji,jj), hu(ji,jj,Kmm) ), ze3(ji,jj,ikt) )
!!$         END_2D
!!$         CALL isf_tbl_lvl( ze3, ktop, zhtbl, & ! <<== in
!!$            &                   ikbot, zfrac ) ! ==>> out
!!$         !
!!$         ! compute tbl property at U point
!!$         CALL isf_tbl_avg( miku(A2D(0)), ikbot, zhtbl, zfrac, ze3, pvarin, & ! <<== in
!!$            &                                                      pvarout ) ! ==>> out
!!$         !
!!$      CASE ( 'V' )
!!$         !
!!$         DO_3D( 0 ,0, 0, 0, 1, jpk )
!!$            ze3(ji,jj,jk) = e3v(ji,jj,jk,Kmm)
!!$         END_3D
!!$         !
!!$         ! compute tbl lvl and thickness
!!$         DO_2D( 0 ,0, 0, 0 )
!!$            ikt = ktop(ji,jj)  ! tbl top indices
!!$            zhtbl(ji,jj) = MAX( MIN( phtbl(ji,jj), hv(ji,jj,Kmm) ), ze3(ji,jj,ikt) )
!!$         END_2D
!!$         CALL isf_tbl_lvl( ze3, ktop, zhtbl, & ! <<== in
!!$            &                   ikbot, zfrac ) ! ==>> out
!!$         !
!!$         ! compute tbl property at V point
!!$         CALL isf_tbl_avg( mikv(A2D(0)), ikbot, zhtbl, zfrac, ze3, pvarin, & ! <<== in
!!$            &                                                      pvarout ) ! ==>> out
!!$         !
!!$      CASE ( 'T' )
!!$         !
!!$         DO_3D( 0 ,0, 0, 0, 1, jpk )
!!$            ze3(ji,jj,jk) = e3t(ji,jj,jk,Kmm)
!!$         END_3D
!!$         !
!!$         ! compute tbl property at T point
!!$         CALL isf_tbl_avg( ktop, kbot, phtbl, pfrac, ze3, pvarin, & ! <<== in
!!$            &                                             pvarout ) ! ==>> out
!!$         !
!!$      END SELECT
!!$      !
!!$   END SUBROUTINE isf_tbl


   SUBROUTINE isf_tbl_avg( ktop, kbot, phtbl, pfrac, pe3, pvarin, pvarout )
      !!--------------------------------------------------------------------
      INTEGER,  DIMENSION(A2D(1))    , INTENT(in ) ::   ktop         ! top level of the top boundary layer
      INTEGER,  DIMENSION(A2D(1))    , INTENT(in ) ::   kbot         ! bottom level of the top boundary layer
      REAL(wp), DIMENSION(A2D(1))    , INTENT(in ) ::   phtbl, pfrac ! fraction of bottom level to be affected by the tbl
      REAL(wp), DIMENSION(A2D(0),jpk), INTENT(in ) ::   pe3          ! vertical scale factor
      REAL(wp), DIMENSION(:,:,:)     , INTENT(in ) ::   pvarin       ! tbl property to average between ktop, kbot over phtbl
      REAL(wp), DIMENSION(A2D(0))    , INTENT(out) ::   pvarout      ! tbl property averaged over phtbl between level ktop and kbot
      !!--------------------------------------------------------------------
      CALL isf_tbl_avg_t( ktop, kbot, phtbl, pfrac, pe3, pvarin, lbnd_ij(pvarin), pvarout )
   END SUBROUTINE isf_tbl_avg

   SUBROUTINE isf_tbl_avg_t( ktop, kbot, phtbl, pfrac, pe3, pvarin, ktvarin, pvarout )
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE isf_tbl_avg  ***
      !!
      !! ** Purpose : compute mean property in the boundary layer
      !!
      !! ** Method  : Depth average is made between the top level ktop and the bottom level kbot
      !!              over a thickness phtbl. The bottom level is partially counted (pfrac).
      !!
      !!--------------------------------------------------------------------
      INTEGER,  DIMENSION(2)                , INTENT(in ) ::   ktvarin
      INTEGER,  DIMENSION(A2D(1))           , INTENT(in ) ::   ktop         ! top level of the top boundary layer
      INTEGER,  DIMENSION(A2D(1))           , INTENT(in ) ::   kbot         ! bottom level of the top boundary layer
      REAL(wp), DIMENSION(A2D(1))           , INTENT(in ) ::   phtbl, pfrac ! fraction of bottom level to be affected by the tbl
      REAL(wp), DIMENSION(A2D(0),jpk)       , INTENT(in ) ::   pe3          ! vertical scale factor
      REAL(wp), DIMENSION(AB2D(ktvarin),JPK), INTENT(in ) ::   pvarin       ! tbl property to average between ktop, kbot over phtbl
      REAL(wp), DIMENSION(A2D(0))           , INTENT(out) ::   pvarout      ! tbl property averaged over phtbl between level ktop and kbot
      !!--------------------------------------------------------------------
      INTEGER  ::   ji, jj                        ! loop indices
      INTEGER  ::   ikt, ikb                      ! top and bottom levels
      !!--------------------------------------------------------------------
      !
      ! compute tbl top.bottom level and thickness
      DO_2D( 0, 0, 0, 0 )
         !
         ! tbl top/bottom indices initialisation
         ikt = ktop(ji,jj) ; ikb = kbot(ji,jj)
         !
         ! level fully include in the ice shelf boundary layer
         pvarout(ji,jj) = SUM( pvarin(ji,jj,ikt:ikb-1) * pe3(ji,jj,ikt:ikb-1) ) / phtbl(ji,jj)
         !
         ! level partially include in ice shelf boundary layer 
         pvarout(ji,jj) = pvarout(ji,jj) + pvarin(ji,jj,ikb) * pe3(ji,jj,ikb) / phtbl(ji,jj) * pfrac(ji,jj)
         !
      END_2D

   END SUBROUTINE isf_tbl_avg_t

   SUBROUTINE isf_tbl_lvl( pe3, ktop, phtbl, kbot, pfrac )
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE isf_tbl_lvl  ***
      !!
      !! ** Purpose : - compute bottom level off the top boundary layer
      !!              - thickness of the top boundary layer
      !!              - fraction of the bottom level affected by the tbl
      !!
      !!--------------------------------------------------------------------
      REAL(wp), DIMENSION(A2D(1),jpk), INTENT(in ) ::   pe3    ! vertical scale factor
      INTEGER,  DIMENSION(A2D(1))    , INTENT(in ) ::   ktop   ! top level of the top boundary layer
      REAL(wp), DIMENSION(A2D(1))    , INTENT(in ) ::   phtbl  ! top boundary layer thickness
      INTEGER,  DIMENSION(A2D(1))    , INTENT(out) ::   kbot   ! bottom level of the top boundary layer
      REAL(wp), DIMENSION(A2D(1))    , INTENT(out) ::   pfrac  ! fraction of bottom level in the tbl
      !!---------------------------------------------------------------------
      INTEGER ::   ji, jj
      INTEGER ::   ikt, ikb
      !!---------------------------------------------------------------------
      !
      DO_2D( 1, 1, 1, 1 )
         !
         ikt = ktop(ji,jj)  ! tbl top indices initialisation
         !
         ! --- get kbot --- !
         !                  ! determine the deepest level influenced by the boundary layer
         ikb = ikt+1
         DO WHILE( SUM( pe3(ji,jj,ikt:ikb-1) ) < phtbl(ji,jj ) ) ;  ikb = ikb + 1 ;  END DO
         kbot(ji,jj) = ikb - 1
         !
         ikb = kbot(ji,jj)   ! tbl bottom indices initialisation
         !
         ! --- get pfrac --- !
         !                   ! proportion of the bottom cell included in ice shelf boundary layer 
         pfrac(ji,jj) = ( phtbl(ji,jj) - SUM( pe3(ji,jj,ikt:ikb-1) ) ) / pe3(ji,jj,ikb)
         !
      END_2D
      !
   END SUBROUTINE isf_tbl_lvl
   !
   SUBROUTINE isf_tbl_ktop( pdep, ktop )
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE isf_tbl_top  ***
      !!
      !! ** Purpose : compute top level of the isf top boundary layer in case of an ice shelf parametrisation
      !!
      !!--------------------------------------------------------------------
      REAL(wp), DIMENSION(A2D(1)), INTENT(inout) :: pdep        ! top depth of the parametrisation influence
      INTEGER,  DIMENSION(A2D(1)), INTENT(  out) :: ktop        ! top level affected by the ice shelf parametrisation
      !!--------------------------------------------------------------------
      INTEGER :: ji, jj
      INTEGER :: ikt
      !!--------------------------------------------------------------------
      !
      ! if we need to recompute the top level at every time stepcompute top level (z*, z~) 
      ! in case of weak ht variation we can assume the top level of htbl to be constant
      ! => only done using gdepw_0
      ! be sure pdep is already correctly bounded
      ! test: this routine run on isfdraft should return mikt
      ! test: this routine run with pdep = 0 should return 1
      !
      DO_2D( 1, 1, 1, 1 )
         ! compute ktop
         ikt = 2
         DO WHILE ( gdepw_0(ji,jj,ikt) <= pdep(ji,jj) ) ;  ikt = ikt + 1 ;  END DO
         ktop(ji,jj) = ikt - 1
         !
         ! update pdep
         ikt=ktop(ji,jj)
         pdep(ji,jj) = gdepw_0(ji,jj,ikt)
      END_2D
      !
   END SUBROUTINE isf_tbl_ktop

END MODULE isftbl
