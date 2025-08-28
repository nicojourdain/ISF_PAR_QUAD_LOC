MODULE isfstp
   !!======================================================================
   !!                       ***  MODULE  isfstp  ***
   !! Ice Shelves :  compute iceshelf load, melt and heat flux
   !!======================================================================
   !! History :  3.2  !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X  !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4  !  2013-03  (P. Mathiot) Merging + parametrization
   !!            4.1  !  2019-09  (P. Mathiot) Split param/explicit ice shelf and re-organisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isfstp       : compute iceshelf melt and heat flux
   !!----------------------------------------------------------------------
   USE isf_oce                                      ! isf variables
   USE isfpar                                       ! ice shelf parametrisation
   USE isfcav                                       ! ice shelf cavity
   USE isfrst , ONLY: isfrst_write                  ! ice shelf restart read/write subroutine
   USE isfload, ONLY: isf_load                      ! ice shelf load
   USE isftbl , ONLY: isf_tbl_lvl                   ! ice shelf boundary layer
   USE isfcpl , ONLY: isfcpl_rst_write, isfcpl_init ! isf variables

   USE oce    , ONLY: ssh
   USE par_oce        ! ocean space and time domain
   USE dom_oce        ! ocean space and time domain
   USE zdfdrg,  ONLY: r_Cdmin_top, r_ke0_top            ! vertical physics: top/bottom drag coef.
   !
   USE lib_mpp, ONLY: ctl_stop, ctl_nam
   USE in_out_manager ! I/O manager
   USE lbclnk
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   isf_stp, isf_init, isf_nam  ! routine called in sbcmod and divhor

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
 
   SUBROUTINE isf_stp( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_stp  ***
      !!
      !! ** Purpose : compute total heat flux and total fwf due to ice shelf melt
      !!
      !! ** Method  : For each case (parametrisation or explicity cavity) :
      !!              - define the before fields
      !!              - compute top boundary layer properties
      !!                (in case of parametrisation, this is the 
      !!                 depth range model array used to compute mean far fields properties)
      !!              - compute fluxes
      !!              - write restart variables
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time step
      INTEGER, INTENT(in) ::   Kmm   ! ocean time level index
      !
      INTEGER ::   ji, jj, jk, ikt                     ! loop index
      REAL(wp), DIMENSION(A2D(1))     ::   zhtmp  ! temporary array for thickness
      REAL(wp), DIMENSION(A2D(1),jpk) ::   ze3t   ! 3D workspace for key_qco
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('isf')
      !
      ! temporary arrays for key_qco
      DO_2D( 1, 1, 1, 1 )
         zhtmp(ji,jj) = ht(ji,jj,Kmm)
      END_2D
      DO_3D( 1, 1, 1, 1, 1, jpk )
         ze3t(ji,jj,jk) = e3t(ji,jj,jk,Kmm)
      END_3D
      !
      !=======================================================================
      ! 1.: compute melt and associated heat fluxes in the ice shelf cavities
      !=======================================================================
      !
      IF ( ln_isfcav_mlt ) THEN
         !
         ! --- before time step --- ! 
#if ! defined key_RK3
         IF ( kt /= nit000 ) THEN         ! MLF : need risf_cav_tsc_b update 
            DO_2D( 2, 2, 2, 2 )
               fwfisf_cav_b(ji,jj) = fwfisf_cav(ji,jj)
            END_2D
            DO_3D( 0, 0, 0, 0, 1, jpts )
               risf_cav_tsc_b(ji,jj,jk) = risf_cav_tsc(ji,jj,jk)
            END_3D
         END IF
#endif
         !
         ! --- deepest level (misfkb), thickness (rhisf) & fraction of deepest cell affected by tbl (rfrac) --- !
         DO_2D( 1, 1, 1, 1 )
            ! limit the tbl to water depth and to the top level thickness
            ikt = misfkt_cav(ji,jj)  ! tbl top indices
            rhisf_tbl_cav(ji,jj) = MAX( MIN( rn_htbl * mskisf_cav(ji,jj), zhtmp(ji,jj) ), ze3t(ji,jj,ikt) )
         END_2D

         CALL isf_tbl_lvl( ze3t, misfkt_cav, rhisf_tbl_cav, &  ! <<== in
            &                    misfkb_cav, rfrac_tbl_cav )   ! ==>> out
         !
         ! --- ice shelf melt (fwfisf) and temperature trend (risf) --- !
         CALL isf_cav( kt, Kmm, risf_cav_tsc, fwfisf_cav(A2D(0)) )    ! <<==>> inout
         !
      END IF
      ! 
      !=================================================================================
      ! 2.: compute melt and associated heat fluxes for not resolved ice shelf cavities
      !=================================================================================
      !
      IF ( ln_isfpar_mlt ) THEN
         !
         ! --- before time step --- ! 
#if ! defined key_RK3
         IF ( kt /= nit000 ) THEN          ! MLF : need risf_par_tsc_b update
            DO_2D( 2, 2, 2, 2 )
               fwfisf_par_b(ji,jj)   = fwfisf_par(ji,jj)
            END_2D
            DO_3D( 0, 0, 0, 0, 1, jpts )
               risf_par_tsc_b(ji,jj,jk) = risf_par_tsc(ji,jj,jk)
            END_3D
         END IF
#endif
         !
         ! --- deepest level (misfkb), thickness (rhisf) & fraction of deepest cell affected by tbl (rfrac) --- !
         ! by simplicity, we assume the top level where param applied do not change with time (done in init part)
         !      limit the tbl to water depth and to the top level thickness
         DO_2D( 1, 1, 1, 1 )
            ikt = misfkt_par(ji,jj)  ! tbl top indices
            rhisf_tbl_par(ji,jj) = MAX( MIN( rhisf0_tbl_par(ji,jj), zhtmp(ji,jj) ), ze3t(ji,jj,ikt) )
         END_2D

         CALL isf_tbl_lvl( ze3t, misfkt_par, rhisf_tbl_par, &  ! <<== in
            &                    misfkb_par, rfrac_tbl_par )   ! ==>> out
         !
         ! --- ice shelf melt (fwfisf) and temperature trend (risf) --- !
         CALL isf_par( kt, Kmm, risf_par_tsc, fwfisf_par(A2D(0)) )    ! <<==>> inout
         !
      END IF
      !
      !
      CALL lbc_lnk( 'isf_stp', fwfisf_par  , 'T', 1.0_wp, fwfisf_cav  , 'T', 1.0_wp )
      CALL lbc_lnk( 'isf_stp', risf_par_tsc(:,:,jp_tem), 'T', 1.0_wp, risf_par_tsc(:,:,jp_sal), 'T', 1.0_wp)

      !clem: these lbc are needed since we calculate fwf only in the interior
      IF( ln_isfcpl ) THEN
         CALL lbc_lnk( 'isf_stp', fwfisf_par  , 'T', 1.0_wp, fwfisf_cav  , 'T', 1.0_wp, &
!!clem#if ! defined key_RK3
!!            &                     fwfisf_par_b, 'T', 1.0_wp, fwfisf_cav_b, 'T', 1.0_wp, &
!!#endif
            &                     risfcpl_ssh, 'T', 1.0_wp, risfcpl_cons_ssh, 'T', 1.0_wp ) ! needed in dynspg_ts, stp2d
         CALL lbc_lnk( 'isf_stp', risfcpl_vol, 'T', 1.0_wp, risfcpl_cons_vol, 'T', 1.0_wp ) ! needed in dynspg_ts, stp2d, sshwzv, dynatf
      ELSE
         CALL lbc_lnk( 'isf_stp', fwfisf_par  , 'T', 1.0_wp, fwfisf_cav  , 'T', 1.0_wp  &
!!clem#if ! defined key_RK3
!!            &                   , fwfisf_par_b, 'T', 1.0_wp, fwfisf_cav_b, 'T', 1.0_wp  &
!!#endif
            &        )     
      ENDIF
      !
      !==================
      ! 3.: write restart
      !==================
#if ! defined key_RK3
      ! MLF: write restart variables (qoceisf, qhcisf, fwfisf for now and before)
      IF( ln_isfcav_mlt .AND. lrst_oce )   CALL isfrst_write( kt, 'cav', risf_cav_tsc , fwfisf_cav )
      ! MLF: write restart variables (qoceisf, qhcisf, fwfisf for now and before)
      IF( ln_isfpar_mlt .AND. lrst_oce )   CALL isfrst_write( kt, 'par', risf_par_tsc , fwfisf_par )
#endif
      IF( ln_isfcpl     .AND. lrst_oce )   CALL isfcpl_rst_write( kt, Kmm )
      !
      IF( ln_timing )   CALL timing_stop('isf')
      !
   END SUBROUTINE isf_stp

   
   SUBROUTINE isf_init( Kbb, Kmm, Kaa )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isfstp_init  ***
      !!
      !! ** Purpose :   Initialisation of the ice shelf public variables
      !!
      !! ** Method  :   Read the namisf namelist, check option compatibility and set derived parameters
      !!
      !! ** Action  : - read namisf parameters
      !!              - allocate memory
      !!              - output print
      !!              - ckeck option compatibility
      !!              - call cav/param/isfcpl init routine
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm, Kaa   ! ocean time level indices
      INTEGER ::   ji, jj, jk                  ! loop index
      !!----------------------------------------------------------------------
      !
      ! constrain: l_isfoasis need to be known
      !
      CALL isf_nam()                                              ! Read namelist
      !
      IF( .NOT. ln_isf ) RETURN
      !
      CALL isf_alloc()                                            ! Allocate public array
      !
      ! initalisation of fwf and tsc array to 0
      DO_2D( 0, 0, 0, 0 )
         fwfisf_oasis(ji,jj) = 0._wp
      END_2D
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         risfload    (ji,jj) = 0._wp
         fwfisf_par  (ji,jj) = 0._wp
         fwfisf_cav  (ji,jj) = 0._wp
      END_2D
      DO_3D( 0, 0, 0, 0, 1, jpts )
         risf_cav_tsc(ji,jj,jk) = 0._wp
         risf_par_tsc(ji,jj,jk) = 0._wp
      END_3D
#if ! defined key_RK3
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         fwfisf_par_b(ji,jj) = 0._wp
         fwfisf_cav_b(ji,jj) = 0._wp
      END_2D
      DO_3D( 0, 0, 0, 0, 1, jpts )
         risf_cav_tsc_b(ji,jj,jk) = 0._wp
         risf_par_tsc_b(ji,jj,jk) = 0._wp
      END_3D
#endif
      !
      CALL isf_ctl()                                              ! check option compatibility
      !
      IF( ln_isfcav )   CALL isf_load( Kmm, risfload )            ! compute ice shelf load
      !
      ! terminate routine now if no ice shelf melt formulation specify
      IF( ln_isf ) THEN
         !
         IF( ln_isfcav_mlt )   CALL isf_cav_init()                ! initialisation melt in the cavity
         !
         IF( ln_isfpar_mlt )   CALL isf_par_init()                ! initialisation parametrised melt
         !
         IF( ln_isfcpl     )   CALL isfcpl_init( Kbb, Kmm, Kaa )  ! initialisation ice sheet coupling
         !
      END IF
         
  END SUBROUTINE isf_init

  
  SUBROUTINE isf_ctl()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_ctl  ***
      !!
      !! ** Purpose :   output print and option compatibility check
      !!
      !!----------------------------------------------------------------------
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'isf_init : ice shelf initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namisf :'
         !
         WRITE(numout,*) '   ice shelf cavity (open or parametrised)  ln_isf = ', ln_isf
         WRITE(numout,*)
         !
         IF ( ln_isf ) THEN
#if defined key_qco && ! defined key_isf 
            IF (ln_isfcav ) CALL ctl_stop( 'STOP', 'isf_ctl: ice shelf requires both ln_isf=T AND key_isf activated' ) 
#endif 
            WRITE(numout,*) '      Add debug print in isf module           ln_isfdebug     = ', ln_isfdebug
            WRITE(numout,*)
            WRITE(numout,*) '      melt inside the cavity                  ln_isfcav_mlt   = ', ln_isfcav_mlt
            IF ( ln_isfcav_mlt) THEN
               WRITE(numout,*) '         melt formulation                         cn_isfcav_mlt= ', TRIM(cn_isfcav_mlt)
               WRITE(numout,*) '         thickness of the top boundary layer      rn_htbl      = ', rn_htbl
               WRITE(numout,*) '         gamma formulation                        cn_gammablk  = ', TRIM(cn_gammablk) 
               IF ( TRIM(cn_gammablk) .NE. 'spe' ) THEN 
                  WRITE(numout,*) '         gammat coefficient                       rn_gammat0   = ', rn_gammat0  
                  WRITE(numout,*) '         gammas coefficient                       rn_gammas0   = ', rn_gammas0  
                  WRITE(numout,*) '         top background ke used (from namdrg_top) rn_ke0       = ', r_ke0_top
                  WRITE(numout,*) '         top drag coef.    used (from namdrg_top) rn_Cd0       = ', r_Cdmin_top
               END IF
            END IF
            WRITE(numout,*) ''
            !
            WRITE(numout,*) '      ice shelf melt parametrisation          ln_isfpar_mlt   = ', ln_isfpar_mlt
            IF ( ln_isfpar_mlt ) THEN
               WRITE(numout,*) '         isf parametrisation formulation         cn_isfpar_mlt   = ', TRIM(cn_isfpar_mlt)
            END IF
            WRITE(numout,*) ''
            !
            WRITE(numout,*) '      Coupling to an ice sheet model          ln_isfcpl       = ', ln_isfcpl
            IF ( ln_isfcpl ) THEN
               WRITE(numout,*) '         conservation activated ln_isfcpl_cons     = ', ln_isfcpl_cons
               WRITE(numout,*) '         number of call of the extrapolation loop  = ', nn_drown
            ENDIF
            WRITE(numout,*) ''
            !
         ELSE
            !
            IF ( ln_isfcav ) THEN
               WRITE(numout,*) ''
               WRITE(numout,*) '   W A R N I N G: ice shelf cavities are open BUT no melt will be computed or read from file !'
               WRITE(numout,*) ''
            END IF
            !
         END IF

         IF (ln_isfcav) THEN
            WRITE(numout,*) '      Ice shelf load method                   cn_isfload        = ', TRIM(cn_isfload)
            WRITE(numout,*) '         Temperature used to compute the ice shelf load            = ', rn_isfload_T
            WRITE(numout,*) '         Salinity    used to compute the ice shelf load            = ', rn_isfload_S
         END IF
         WRITE(numout,*) ''
         FLUSH(numout)

      END IF
      !

      !---------------------------------------------------------------------------------------------------------------------
      ! sanity check  ! issue ln_isfcav not yet known as well as l_isfoasis  => move this call in isf_stp ?
      ! melt in the cavity without cavity
      IF ( ln_isfcav_mlt .AND. (.NOT. ln_isfcav) ) &
          &   CALL ctl_stop('ice shelf melt in the cavity activated (ln_isfcav_mlt) but no cavity detected in domcfg (ln_isfcav), STOP' )
      !
      ! ice sheet coupling without cavity
      IF ( ln_isfcpl .AND. (.NOT. ln_isfcav) ) &
         &   CALL ctl_stop('coupling with an ice sheet model detected (ln_isfcpl) but no cavity detected in domcfg (ln_isfcav), STOP' )
      !
      IF ( ln_isfcpl .AND. ln_isfcpl_cons .AND. lk_linssh ) &
         &   CALL ctl_stop( 'The coupling between NEMO and an ice sheet model with the conservation option does not work with the linssh option' )
      !
      IF ( l_isfoasis .AND. .NOT. ln_isf ) CALL ctl_stop( ' OASIS send ice shelf fluxes to NEMO but NEMO does not have the isf module activated' )
      !
      IF ( l_isfoasis .AND. ln_isf ) THEN
         !
         ! NEMO coupled to ATMO model with isf cavity need oasis method for melt computation 
         IF ( ln_isfcav_mlt .AND. TRIM(cn_isfcav_mlt) /= 'oasis' ) CALL ctl_stop( 'cn_isfcav_mlt = oasis is the only option availble if fwf send by oasis' )
         IF ( ln_isfpar_mlt .AND. TRIM(cn_isfpar_mlt) /= 'oasis' ) CALL ctl_stop( 'cn_isfpar_mlt = oasis is the only option availble if fwf send by oasis' )
         !
         ! oasis melt computation with cavity open and cavity parametrised (not coded)
         IF ( ln_isfcav_mlt .AND. ln_isfpar_mlt ) THEN
            IF ( TRIM(cn_isfpar_mlt) == 'oasis' .AND. TRIM(cn_isfcav_mlt) == 'oasis' ) CALL ctl_stop( 'cn_isfpar_mlt = oasis and cn_isfcav_mlt = oasis not coded' )
         END IF
         !
      ENDIF
   END SUBROUTINE isf_ctl

   
   SUBROUTINE isf_nam
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_nam  ***
      !!
      !! ** Purpose :   Read ice shelf namelist cfg and ref
      !!
      !!----------------------------------------------------------------------
      INTEGER               :: ios                  ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      NAMELIST/namisf/ ln_isf        ,                                                           & 
         &             cn_gammablk   , rn_gammat0    , rn_gammas0    , rn_htbl, sn_isfcav_fwf,   &
         &             ln_isfcav_mlt , cn_isfcav_mlt , sn_isfcav_fwf , nn_isfpar_basin,          &
         &             ln_isfpar_mlt , cn_isfpar_mlt , sn_isfpar_fwf , sn_isfpar_basin,          &
         &             sn_isfpar_zmin, sn_isfpar_zmax, sn_isfpar_Leff, sn_isfpar_area,           &
         &             ln_isfcpl     , nn_drown      , ln_isfcpl_cons, ln_isfdebug,              &
         &             cn_isfload    , rn_isfload_T  , rn_isfload_S  , cn_isfdir  ,              &
         &             rn_isfpar_bg03_gt0, rn_isfpar_Kcoeff
      !!----------------------------------------------------------------------
      !
      READ_NML_REF(numnam,namisf)
      READ_NML_CFG(numnam,namisf)
      IF(lwm) WRITE ( numond, namisf )

   END SUBROUTINE isf_nam
   !!
   !!======================================================================
END MODULE isfstp
