MODULE isfcav
   !!======================================================================
   !!                       ***  MODULE  isfcav  ***
   !! Ice shelf cavity module :  update ice shelf melting under ice
   !!                            shelf
   !!======================================================================
   !! History :  3.2  !  2011-02  (C.Harris  ) Original code isf cav
   !!            3.4  !  2013-03  (P. Mathiot) Merging + parametrization
   !!            4.1  !  2019-09  (P. Mathiot) Split ice shelf cavity and ice shelf parametrisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isf_cav       : update ice shelf melting under ice shelf
   !!----------------------------------------------------------------------
   USE isf_oce        ! ice shelf public variables
   !
   USE isfrst   , ONLY: isfrst_read                ! ice shelf restart read/write subroutine
   USE isfutils , ONLY: debug                      ! ice shelf debug subroutine
   USE isftbl                                      ! ice shelf top boundary layer properties subroutine
   USE isfcavmlt, ONLY: sf_isfcav_fwf, isfcav_mlt  ! ice shelf melt formulation subroutine
   USE isfcavgam, ONLY: isfcav_gammats             ! ice shelf melt exchange coeficient subroutine
   USE isfdiags , ONLY: isf_diags_flx              ! ice shelf diags subroutine
   USE zdfdrg   , ONLY: rCd0_top, r_ke0_top        ! vertical physics: top/bottom drag coef.
   !
   USE oce      , ONLY: ts, uu, vv, rn2                 ! ocean dynamics and tracers
   USE dom_oce                                          ! ocean space and time domain
   USE par_oce                                          ! ocean space and time domain
   USE phycst                                           ! physical constants
   USE eosbn2   , ONLY: eos_fzp, ln_teos10              ! equation of state & use ln_teos10 or not
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE fldread        ! read input field at current time step
   USE lbclnk         ! lbclnk
   USE lib_mpp        ! MPP library

   IMPLICIT NONE

   PRIVATE

   PUBLIC   isf_cav, isf_cav_init ! routine called in isfmlt
   !
   ! Cavity file information
   TYPE(FLD_N), PUBLIC :: sn_isfcav_fwf   !: information about the isf melting file to be read
   !
   !! * Substitutions   
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE isf_cav( kt, Kmm, ptsc, pfwf )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE isf_cav  ***
      !!
      !! ** Purpose :   handle surface boundary condition under ice shelf
      !!
      !! ** Method  :   based on Mathiot et al. (2017)
      !!
      !! ** Action  :   - compute geometry of the Losch top bournary layer (see Losch et al. 2008)
      !!                - depending on the chooses option
      !!                   - compute temperature/salt in the tbl
      !!                   - compute exchange coeficient
      !!                   - compute heat and fwf fluxes
      !!                   - output
      !! 
      !! ** Convention : all fluxes are from isf to oce
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      REAL(wp), DIMENSION(A2D(0),jpts), INTENT(inout) ::   ptsc   ! T & S ice shelf cavity contents
      REAL(wp), DIMENSION(A2D(0))     , INTENT(inout) ::   pfwf   ! ice shelf fwf
      !!
      INTEGER  ::   iconv, ji, jj, jk, ikt
      REAL(wp) ::   zerr
      REAL(wp) ::   zcoef, zdku, zdkv
      INTEGER , PARAMETER ::   iconv_max = 100
      LOGICAL , DIMENSION(A2D(0)) ::   l_converged              ! true when melting converges (per grid point)
      REAL(wp), DIMENSION(A2D(0)) ::   ztfrz                    ! tbl freezing temperature
      REAL(wp), DIMENSION(A2D(0)) ::   zqlat, zqoce, zqhc       ! heat fluxes
      REAL(wp), DIMENSION(A2D(0)) ::   zqh_b, zRc               !
      REAL(wp), DIMENSION(A2D(0)) ::   zgammat, zgammas         ! exchange coeficient
      REAL(wp), DIMENSION(A2D(0)) ::   zttbl, zstbl             ! temp. and sal. in top boundary layer
      REAL(wp), DIMENSION(A2D(0)) ::   zustar                   ! u*
      REAL(wp), DIMENSION(A2D(0)) ::   zutbl, zvtbl             ! top boundary layer velocity
      REAL(wp), DIMENSION(A2D(0),jpk) ::   zvel, ze3            ! u and v at T points and e3t
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztmp           ! temporary array
      !!---------------------------------------------------------------------
      !
      DO_3D( 0 ,0, 0, 0, 1, jpk )
         ze3(ji,jj,jk) = e3t(ji,jj,jk,Kmm)
      END_3D
      !
      !===============================
      ! 1.: compute T and S in the tbl
      !===============================
      CALL isf_tbl_avg( misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav, ze3, ts(:,:,:,jp_tem,Kmm), & ! <<== in
         &                                                                                        zttbl  ) ! ==>> out
      !
      CALL isf_tbl_avg( misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav, ze3, ts(:,:,:,jp_sal,Kmm), & ! <<== in
         &                                                                                        zstbl  ) ! ==>> out
      !
      !==========================================
      ! 2.: compute velocity in the tbl if needed
      !==========================================
      IF ( TRIM(cn_gammablk) == 'vel_stab' .OR. TRIM(cn_gammablk) == 'vel' ) THEN
         ! compute velocity in tbl
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zvel(ji,jj,jk) = 0.5_wp * ( uu(ji-1,jj,jk,Kmm) + uu(ji,jj,jk,Kmm) )
         END_3D
         CALL isf_tbl_avg( misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav, ze3, zvel, & ! <<== in
            &                                                                        zutbl  ) ! ==>> out
         !
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zvel(ji,jj,jk) = 0.5_wp * ( vv(ji,jj-1,jk,Kmm) + vv(ji,jj,jk,Kmm) )
         END_3D
         CALL isf_tbl_avg( misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav, ze3, zvel, & ! <<== in
            &                                                                        zvtbl  ) ! ==>> out
         !
         ! compute ustar (AD15 eq. 27)
         DO_2D( 0, 0, 0, 0 )         
            zustar(ji,jj) = SQRT( rCd0_top(ji,jj) * ( ( zutbl(ji,jj)*zutbl(ji,jj) + zvtbl(ji,jj)*zvtbl(ji,jj) ) + r_ke0_top ) &
               &                ) * mskisf_cav(ji,jj)
         END_2D
         !
         !
      ELSE
         zustar(:,:) = 0._wp ! not used
      ENDIF
      !
      !==============================
      ! 3.: compute ice shelf melting
      !==============================
      ! --- initialisation --- !
      IF ( TRIM(cn_gammablk) == 'vel_stab' ) THEN
         !
         DO_2D( 0, 0, 0, 0 )
            zqoce(ji,jj) = -pfwf(ji,jj) * rLfusisf        * mskisf_cav(ji,jj)
            zqh_b(ji,jj) =  ptsc(ji,jj,jp_tem) * rho0_rcp * mskisf_cav(ji,jj) ! last time step total heat fluxes (to speed up convergence)

            ikt = mikt(ji,jj)
            ! compute Rc number (as done in zdfric.F90)
!!gm better to do it like in the new zdfric.F90   i.e. avm weighted Ri computation
            zcoef = 0.5_wp / e3w(ji,jj,ikt+1,Kmm)
            !                                            ! shear of horizontal velocity
            zdku = zcoef * (  ( uu(ji-1,jj  ,ikt  ,Kmm) + uu(ji,jj,ikt  ,Kmm) )  &   ! add () for NP repro
               &            - ( uu(ji-1,jj  ,ikt+1,Kmm) + uu(ji,jj,ikt+1,Kmm) )  )
            zdkv = zcoef * (  ( vv(ji  ,jj-1,ikt  ,Kmm) + vv(ji,jj,ikt  ,Kmm) )  &   ! add () for NP repro
               &            - ( vv(ji  ,jj-1,ikt+1,Kmm) + vv(ji,jj,ikt+1,Kmm) )  )
            !                                            ! richardson number (minimum value set to zero)
            zRc(ji,jj) = MAX( rn2(ji,jj,ikt+1), 1.e-20_wp ) / MAX( zdku*zdku + zdkv*zdkv, 1.e-20_wp )
         END_2D
         !
      ENDIF
      !
      ! Calculate freezing temperature (overwritten only in case cn_isfcav_mlt="3eq")
      CALL eos_fzp( zstbl(:,:), ztfrz(:,:), risfdep(:,:), kbnd=0 )
      !
      ! --- iteration loop --- !
      iconv = 0
      l_converged(:,:) = .FALSE.
      ! Convergence calculated until all sub-domain grid points have converged
      ! Calculations keep going for all grid points until sub-domain convergence (vectorisation optimisation)
      ! but values are not taken into account (results independant of MPI partitioning)
      !
      DO WHILE ( ( .NOT. ALL (l_converged(:,:)) ) .AND. iconv < iconv_max )
         iconv = iconv + 1
         !
         ! compute gammat everywhere (2d) - useless if melt specified
         IF ( TRIM(cn_isfcav_mlt) .NE. 'spe' ) THEN
            CALL isfcav_gammats( Kmm, l_converged, zustar, zttbl, zstbl, zqoce, pfwf, zRc,  & ! <<== in
               &                                                          zgammat, zgammas  ) ! ==>> inout
         END IF
         !
         ! compute tfrz, latent heat and melt (2d)
         CALL isfcav_mlt( kt, l_converged, zgammat, zgammas, zttbl, zstbl, & ! <<== in
            &                                    zqhc, zqoce, pfwf, ztfrz  ) ! ==>> inout
         !
         ! define if we need to iterate
         SELECT CASE ( cn_gammablk )
         CASE ( 'spe', 'vel' )
            l_converged(:,:) = .TRUE.                ! no convergence needed
         CASE ( 'vel_stab' )
            DO_2D( 0, 0, 0, 0 )                      ! compute error between 2 iterations
               IF ( .NOT. l_converged(ji,jj) ) THEN
                  zerr = ABS( zqhc(ji,jj) + zqoce(ji,jj) - zqh_b(ji,jj) )
                  IF( zerr <= 0.01_wp ) THEN         ! convergence is reached
                     l_converged(ji,jj) = .TRUE.
                  ELSE                               ! converge is not yet reached
                     l_converged(ji,jj) = .FALSE.
                     zqh_b(ji,jj) = zqhc(ji,jj) + zqoce(ji,jj)
                  ENDIF
               ENDIF
            END_2D
            !
         END SELECT
         !
      END DO
      IF( iconv == iconv_max )   CALL ctl_stop( 'STOP', 'isf_cav: vel_stab gamma formulation had too many iterations ...' )   ! too many iterations
      !
      !
      DO_2D( 0, 0, 0, 0 )
         ! compute heat content flux ( > 0 from isf to oce)
         zqlat(ji,jj) = - pfwf(ji,jj) * rLfusisf    ! 2d latent heat flux (W/m2)
         !
         ! set temperature content
         ptsc(ji,jj,jp_tem) = ( zqhc(ji,jj) + zqoce(ji,jj) ) * r1_rho0_rcp
      END_2D
      !
      ! output fluxes
      CALL isf_diags_flx( Kmm, misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav, 'cav', pfwf, zqoce, zqlat, zqhc )
      !
      ! --- output --- !
      CALL iom_put( 'ttbl_cav' , zttbl  (:,:) * mskisf_cav(A2D(0)) )
      CALL iom_put( 'stbl'     , zstbl  (:,:) * mskisf_cav(A2D(0)) )
      CALL iom_put( 'isfgammat', zgammat(:,:) * mskisf_cav(A2D(0)) )
      CALL iom_put( 'isfgammas', zgammas(:,:) * mskisf_cav(A2D(0)) )
      IF ( TRIM(cn_gammablk) == 'vel_stab' .OR. TRIM(cn_gammablk) == 'vel' ) THEN
         CALL iom_put( 'isfustar', zustar              )
         CALL iom_put( 'utbl'    , zutbl  * mskisf_cav(A2D(0)) )
         CALL iom_put( 'vtbl'    , zvtbl  * mskisf_cav(A2D(0)) )
      ENDIF
      !
      CALL iom_put('isftfrz_cav'    ,                ztfrz(:,:)   * mskisf_cav(A2D(0)) )   ! freezing point at the interface
      CALL iom_put('isfthermald_cav', ( zttbl(:,:) - ztfrz(:,:) ) * mskisf_cav(A2D(0)) )   ! thermal driving at the interface
      !
      IF( cn_isfcav_mlt == '3eq' ) THEN                                                 ! conductive heat flux through the ice
         ALLOCATE( ztmp(A2D(0)) )
         DO_2D( 0, 0, 0, 0 )
            ztmp(ji,jj) = rhoisf * rcpisf * rkappa / MAX( risfdep(ji,jj), 1.e-20 ) * ( ztfrz(ji,jj) - rtsurf ) * mskisf_cav(ji,jj)
         END_2D
         CALL iom_put('qconisf', ztmp )
         DEALLOCATE( ztmp )
      ENDIF
      !
      ! --- debug --- !
      IF ( ln_isfdebug ) THEN
         IF(lwp) WRITE(numout,*) ''
         CALL debug( 'isfcav: gammaT ', zgammat(:,:) )
         CALL debug( 'isfcav: gammaS ', zgammas(:,:) )
         CALL debug( 'isfcav: ptsc T ', ptsc(:,:,1)  )
         CALL debug( 'isfcav: ptsc S ', ptsc(:,:,2)  )
         CALL debug( 'isfcav: qhc    ', zqhc (:,:)   )
         CALL debug( 'isfcav: qoce   ', zqoce(:,:)   )
         CALL debug( 'isfcav: fwf    ', pfwf (:,:)   )
         IF(lwp) WRITE(numout,*) ''
      END IF
      !
   END SUBROUTINE isf_cav

   SUBROUTINE isf_cav_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_cav_init ***
      !!
      !! ** Purpose : initialisation of variable needed to compute melt under an ice shelf
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ierr
      INTEGER ::   ji, jj     ! dummy loop indices
      !!---------------------------------------------------------------------
      !
      !==============
      ! 0: allocation
      !==============
      CALL isf_alloc_cav()
      !
      !==================
      ! 1: initialisation
      !==================
      DO_2D( 1, 1, 1, 1 )
         ! top and bottom level of the 'top boundary layer'
         misfkt_cav(ji,jj)    = mikt(ji,jj)
         misfkb_cav(ji,jj)    = 1
         !
         ! thickness of 'tbl' and fraction of bottom cell affected by 'tbl'
         rhisf_tbl_cav(ji,jj) = 0.0_wp
         rfrac_tbl_cav(ji,jj) = 0.0_wp
         !
         ! cavity mask
         mskisf_cav(ji,jj) = ( 1._wp - tmask(ji,jj,1) ) * ssmask(ji,jj)
         !
      END_2D
      !
#if ! defined key_RK3
      !================
      ! 2: activate restart
      !================
      !
      !================
      ! 3: read restart
      !================
      !
      ! MLF: read cav variable from restart
      IF ( ln_rstart ) CALL isfrst_read('cav', risf_cav_tsc, fwfisf_cav, risf_cav_tsc_b, fwfisf_cav_b)
      !
#endif
      !==========================================
      ! 3: specific allocation and initialisation (depending of scheme choice)
      !==========================================
      !
      SELECT CASE ( TRIM(cn_isfcav_mlt) )
      CASE( 'spe' )

         ALLOCATE( sf_isfcav_fwf(1), STAT=ierr )
         ALLOCATE( sf_isfcav_fwf(1)%fnow(jpi,jpj,1), sf_isfcav_fwf(1)%fdta(jpi,jpj,1,2) )
         CALL fld_fill( sf_isfcav_fwf, (/ sn_isfcav_fwf /), cn_isfdir, 'isf_cav_init', 'read fresh water flux isf data', 'namisf' )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '  ==>> The ice shelf melt inside the cavity is read from forcing files'

      CASE( '2eq' )
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '  ==>> The original ISOMIP melt formulation is used to compute melt under the ice shelves'

      CASE( '3eq' )
         ! coeficient for linearisation of potential tfreez
         ! Crude approximation for pressure (but commonly used)
         IF ( ln_teos10 ) THEN   ! linearisation from Jourdain et al. (2017)
            risf_lamb1 =-0.0564_wp
            risf_lamb2 = 0.0773_wp
            risf_lamb3 =-7.8633e-8 * grav * rho0
         ELSE                  ! linearisation from table 4 (Asay-Davis et al., 2015)
            risf_lamb1 =-0.0573_wp
            risf_lamb2 = 0.0832_wp
            risf_lamb3 =-7.5300e-8 * grav * rho0
         ENDIF

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '  ==>> The 3 equations melt formulation is used to compute melt under the ice shelves'

      CASE DEFAULT
         CALL ctl_stop(' cn_isfcav_mlt method unknown (spe, 2eq, 3eq), check namelist')
      END SELECT
      !
   END SUBROUTINE isf_cav_init

END MODULE isfcav
