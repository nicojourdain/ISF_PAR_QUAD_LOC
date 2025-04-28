MODULE isfparmlt
   !!======================================================================
   !!                       ***  MODULE  isfparmlt  ***
   !! Ice shelf parametrisation module :  update surface ocean boundary condition under ice
   !!                   shelf using an ice shelf melt parametrisation
   !!======================================================================
   !! History :  4.0  !                        original code
   !!            4.2  ! 2025-04  (N. Jourdain) Interactive quadratic parameterisation
   !!----------------------------------------------------------------------

   USE isf_oce                  ! ice shelf
   USE isftbl , ONLY: isf_tbl   ! ice shelf depth average
   USE isfutils,ONLY: debug     ! debug subroutine

   USE dom_oce                        ! ocean space and time domain
   USE oce    , ONLY: ts              ! ocean dynamics and tracers
   USE phycst , ONLY: rcp, rho0, grav ! physical constants
   USE eosbn2 , ONLY: eos_fzp         ! equation of state

   USE in_out_manager              ! I/O manager
   USE iom        , ONLY: iom_put  ! I/O library
   USE fldread    , ONLY: fld_read, FLD, FLD_N !
   USE lib_fortran, ONLY: glob_sum !
   USE lib_mpp    , ONLY: ctl_stop !

   IMPLICIT NONE

   PRIVATE

   PUBLIC  isfpar_mlt 


   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcisf.F90 10536 2019-01-16 19:21:09Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

! -------------------------------------------------------------------------------------------------------
! -------------------------------- PUBLIC SUBROUTINE ----------------------------------------------------
! -------------------------------------------------------------------------------------------------------

   SUBROUTINE isfpar_mlt( kt, Kmm, pqhc, pqoce, pqfwf )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isfpar_mlt  ***
      !!
      !! ** Purpose : Compute Salt and Heat fluxes related to ice_shelf 
      !!              melting and freezing 
      !!
      !! ** Method  :  2 parameterizations are available according
      !!                        1 : Specified melt flux
      !!                        2 : Beckmann & Goose parameterization
      !!
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pqfwf, pqoce, pqhc  ! fresh water, ice-ocean heat and heat content fluxes
      !!-------------------------- IN  -------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      !!---------------------------------------------------------------------
      !
      ! Choose among the available ice shelf parametrisation
      SELECT CASE ( cn_isfpar_mlt )
      CASE ( 'spe' )          ! specified runoff in depth (Mathiot et al., GMD, 2017)
         CALL isfpar_mlt_spe(kt, Kmm, pqhc, pqoce, pqfwf)
      CASE ( 'bg03' )         ! Linear parametrisation (Beckmann and Goosse, Oce. Mod., 2003)
         CALL isfpar_mlt_bg03(kt, Kmm, pqhc, pqoce, pqfwf)
      CASE ( 'quad_loc' ) ! Quadratic local melt parameterisation (Burgard et al., TC, 2022)
         CALL isfpar_mlt_quad_loc(kt, Kmm, pqhc, pqoce, pqfwf)
      CASE ( 'oasis' )        ! Climate model (ice-shelf melt prescribed together with iceberg melt to conserve the ice sheet mass)
         CALL isfpar_mlt_oasis( kt, Kmm, pqhc, pqoce, pqfwf)
      CASE DEFAULT
         CALL ctl_stop('STOP', 'unknown isf melt formulation : cn_isfpar (should not see this)')
      END SELECT
      !
      IF (ln_isfdebug) THEN
         IF(lwp) WRITE(numout,*) ''
         CALL debug( 'isfpar_mlt qhc  :', pqhc (:,:) )
         CALL debug( 'isfpar_mlt qoce :', pqoce(:,:) )
         CALL debug( 'isfpar_mlt qfwf :', pqfwf(:,:) )
         IF(lwp) WRITE(numout,*) ''
      END IF
      !
   END SUBROUTINE isfpar_mlt

! -------------------------------------------------------------------------------------------------------
! -------------------------------- PRIVATE SUBROUTINE ---------------------------------------------------
! -------------------------------------------------------------------------------------------------------

   SUBROUTINE isfpar_mlt_spe(kt, Kmm, pqhc, pqoce, pqfwf)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isfpar_mlt_spe  ***
      !!
      !! ** Purpose : prescribed ice shelf melting in case ice shelf cavities are closed.
      !!              data read into a forcing files.
      !!
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pqhc, pqfwf, pqoce  ! fresh water and ice-ocean heat fluxes
      !!-------------------------- IN  -------------------------------------
      INTEGER,  INTENT(in) :: kt
      INTEGER,  INTENT(in) :: Kmm    !  ocean time level index
      !!--------------------------------------------------------------------
      INTEGER :: jk
      REAL(wp), DIMENSION(jpi,jpj,jpk)  :: ztfrz3d
      REAL(wp), DIMENSION(jpi,jpj)      :: ztfrz
      !!--------------------------------------------------------------------
      !
      ! 0. ------------Read specified fwf from isf to oce
      CALL fld_read ( kt, 1, sf_isfpar_fwf   )
      !
      ! compute ptfrz
      ! 1. ------------Mean freezing point
      DO jk = 1,jpk
         CALL eos_fzp(ts(:,:,jk,jp_sal,Kmm), ztfrz3d(:,:,jk), gdept(:,:,jk,Kmm))
      END DO
      CALL isf_tbl(Kmm, ztfrz3d, ztfrz, 'T', misfkt_par, rhisf_tbl_par, misfkb_par, rfrac_tbl_par )
      !
      pqfwf(:,:) =   sf_isfpar_fwf(1)%fnow(:,:,1)      ! fresh water flux from the isf (fwfisf <0 mean melting)       ( > 0 from isf to oce)
      pqoce(:,:) = - pqfwf(:,:) * rLfusisf             ! ocean/ice shelf flux assume to be equal to latent heat flux  ( > 0 from isf to oce)
      pqhc (:,:) =   pqfwf(:,:) * ztfrz(:,:) * rcp     ! heat content flux                                            ( > 0 from isf to oce)
      !
      CALL iom_put('isftfrz_par', ztfrz(:,:) * mskisf_par(:,:) )
      !
   END SUBROUTINE isfpar_mlt_spe

   SUBROUTINE isfpar_mlt_bg03(kt, Kmm, pqhc, pqoce, pqfwf)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isfpar_mlt_bg03  ***
      !!
      !! ** Purpose : compute an estimate of ice shelf melting and 
      !!              latent, ocean-ice and heat content heat fluxes
      !!              in case cavities are closed based on the far fields T and S properties. 
      !!
      !! ** Method  : The ice shelf melt is computed as proportional to the differences between the 
      !!              mean temperature and mean freezing point in front of the ice shelf averaged 
      !!              over the ice shelf min ice shelf draft and max ice shelf draft and the freezing point
      !!
      !! ** Reference : Beckmann and Goosse (2003), "A parameterization of ice shelf-ocean
      !!                interaction for climate models", Ocean Modelling 5(2003) 157-170.
      !!----------------------------------------------------------------------
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pqhc, pqfwf, pqoce  ! fresh water and ice-ocean heat fluxes
      !!-------------------------- IN  -------------------------------------
      INTEGER,  INTENT(in) :: kt
      INTEGER,  INTENT(in) :: Kmm    !  ocean time level index
      !!--------------------------------------------------------------------
      INTEGER :: jk
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztfrz3d        ! freezing point
      REAL(wp), DIMENSION(jpi,jpj)     :: ztfrz          ! freezing point
      REAL(wp), DIMENSION(jpi,jpj)     :: ztavg          ! temperature avg
      !!----------------------------------------------------------------------
      !
      ! 0. ------------Mean freezing point
      DO jk = 1,jpk
         CALL eos_fzp(ts(:,:,jk,jp_sal,Kmm), ztfrz3d(:,:,jk), gdept(:,:,jk,Kmm))
      END DO
      CALL isf_tbl(Kmm, ztfrz3d, ztfrz, 'T', misfkt_par, rhisf_tbl_par, misfkb_par, rfrac_tbl_par )
      !
      ! 1. ------------Mean temperature
      CALL isf_tbl(Kmm, ts(:,:,:,jp_tem,Kmm), ztavg, 'T', misfkt_par, rhisf_tbl_par, misfkb_par, rfrac_tbl_par )
      !
      ! 2. ------------Net heat flux and fresh water flux due to the ice shelf
      pqfwf(:,:) =   rho0 * rcp * rn_isfpar_bg03_gt0 * risfLeff(:,:) * e1t(:,:) * (ztavg(:,:) - ztfrz(:,:) ) * r1_e1e2t(:,:) / rLfusisf  ! ( > 0 from isf to oce)
      pqoce(:,:) = - pqfwf(:,:) * rLfusisf             ! ocean/ice shelf flux assume to be equal to latent heat flux  ( > 0 from isf to oce)
      pqhc (:,:) =   pqfwf(:,:) * ztfrz(:,:) * rcp     ! heat content flux                                            ( > 0 from isf to oce)
      !
      ! 3. ------------BG03 output
      ! output ttbl
      CALL iom_put('ttbl_par', ztavg(:,:) * mskisf_par(:,:) )
      !
      ! output thermal driving
      CALL iom_put('isfthermald_par',( ztavg(:,:) - ztfrz(:,:) ) * mskisf_par(:,:))
      !
      ! output freezing point used to define the thermal driving and heat content fluxes
      CALL iom_put('isftfrz_par', ztfrz(:,:) * mskisf_par(:,:) )
      !
   END SUBROUTINE isfpar_mlt_bg03

   SUBROUTINE isfpar_mlt_quad_loc(kt, Kmm, pqhc, pqoce, pqfwf)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isfpar_mlt_quad_loc  ***
      !!
      !! ** Purpose : compute an estimate of ice shelf melting and 
      !!              latent, ocean-ice and heat content heat fluxes
      !!              for non-resolved parts of cavities based on nearby T and S properties. 
      !!
      !! ** Method  : The input files consist of: zmin and zmax, the minimum and maximal depths
      !!              defining the interfacial water colums that are used to calculate the thermal
      !!              forcing and to inject the parameterised melt; a map of basin/sector defined 
      !!              by specific numbers; the non-resolved (i.e. parameterised) ice-shelf area
      !!              per basin and per vertical level.
      !!             
      !!              A mean profile of thermal forcing (averaged over interfacial water columns) 
      !!              is calculated in each basin and used to calculate a melt profile, which is 
      !!              then uniformly distributed over the active water columns.
      !!
      !! ** Reference : Burgard, C., Jourdain, N. C., Reese, R., Jenkins, A., and Mathiot, P. (2022). 
      !!                An assessment of basal melt parameterisations for Antarctic ice shelves, 
      !!                The Cryosphere, 16, 4931–4975, doi:10.5194/tc-16-4931-2022 
      !!                
      !!----------------------------------------------------------------------
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pqhc, pqfwf, pqoce  ! fresh water and ice-ocean heat fluxes
      !!-------------------------- IN  -------------------------------------
      INTEGER,  INTENT(in) :: kt
      INTEGER,  INTENT(in) :: Kmm    !  ocean time level index
      !!--------------------------------------------------------------------
      REAL(wp) ::   zfillvalue = 1.e20_wp    ! missing values in outputs
      REAL(wp) ::   zf         = 1.4e-4_wp   ! mean Coriolis parameter [s^-1]
      REAL(wp) ::   zbeta      = 7.8e-4_wp   ! salt contraction coefficient [1.e3]
      REAL(wp) ::   zsin_theta = 2.9e-3_wp   ! assuming a representative "Antarctic slope" [1]
      !!--------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)     :: ztfrz         ! mean freezing temperature in interfacial water columns [degC]
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztf3d         ! 3d thermal forcing [degC]
      !REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztftfs3d      ! TF*|TF|*Sloc*e1t*e2t (where TF = thermal forcing)
      REAL(wp), DIMENSION(jpk)         :: ztftfs        ! TF*|TF|*Sloc*e1t*e2t (where TF = thermal forcing)
      REAL(wp), DIMENSION(nn_isfpar_basin,jpk) :: zmelt ! Parameterised melt per basin and per level [kg s^-1]
      REAL(wp), DIMENSION(nn_isfpar_basin,jpk) :: zzztf2s ! TF*|TF|*Sloc for interfacial cells  [degC^2 1.e-3]
      INTEGER  :: ji, jj, jk, kbasin         ! dummy loop indices
      REAL(wp) :: zcoef                      ! Constant coefficient used in the param [kg m^-2 s^-1 degC^-2 1.e3]
      !!----------------------------------------------------------------------
      !
      ! Coefficient for the melting parameterization:
      !   see eq. 14 of Burgard et al. (2022) 
      !   NB: their melt is in meters of ice per second while we here use kg/m2/s
      zcoef = rho0 * ( rcp / rLfusisf )**2 * zbeta * grav * zsin_theta * 0.5_wp / zf
      !
      pqfwf(:,:) = 0._wp
      !
      ! ztf3d is the freezing temperature at all levels:
      DO jk = 1,jpk
         CALL eos_fzp(ts(:,:,jk,jp_sal,Kmm), ztf3d(:,:,jk), gdept(:,:,jk,Kmm))
      END DO
      ! ztfrz is the mean freezing temperature in the interfacial water columns, i.e., between
      ! specified zmin and zmax (only used for the heat content flux):
      CALL isf_tbl(Kmm, ztf3d, ztfrz, 'T', misfkt_par, rhisf_tbl_par, misfkb_par, rfrac_tbl_par )
      !
      DO kbasin=1,nn_isfpar_basin
        !
        IF ( ln_exchg(kbasin) ) THEN
          !
          ! Calculate ztftfs3d as TF*|TF|*Sloc*e1t*e2t (where TF = thermal forcing) [degC^2 1.e-3 m^2]:
          DO jk = 1,jpk
            ztftfs(jk) = SUM( ( ts(:,:,jk,jp_tem,Kmm) - ztf3d(:,:,jk) ) * abs( ts(:,:,jk,jp_tem,Kmm) - ztf3d(:,:,jk) ) &
            &                  * ts(:,:,jk,jp_sal,Kmm) * e1e2t(:,:) * mskisf_exchg(:,:,kbasin) * tmask(:,:,jk) ) 
          ENDDO
          CALL mpp_sum( 'isf_par_mlt_quad_loc', ztftfs(:) )
          zzztf2s(kbasin,:) = ztftfs(:) / area_exchg(kbasin,:)
          !DO jk = 1,jpk
          !  ztftfs3d(:,:,jk) = ( ts(:,:,jk,jp_tem,Kmm) - ztf3d(:,:,jk) ) * abs( ts(:,:,jk,jp_tem,Kmm) - ztf3d(:,:,jk) ) &
          !  &                  * ts(:,:,jk,jp_sal,Kmm) * e1e2t(:,:) * mskisf_exchg(:,:,kbasin) * tmask(:,:,jk)
          !END DO
          !! Average profile of TF*|TF|*Sloc in the interfacial ocean grid cells [degC^2 1.e-3]
          !! (no problem if area_exchg=0 at some levels, these values will be replaced in the next loop)
          !zzztf2s(kbasin,:) = glob_sum( 'isfparmlt', ztftfs3d(:,:,:) ) / area_exchg(kbasin,:)
          !
          DO jk = 1,jpk
            !
            ! Melt per basin per vertical level [kg s^-1] :
            ! NB1: rn_isfpar_Kcoeff is defined in eq. 16 of Burgard et al. (2022)
            !      and was calibrated at 1.16e-4 in that paper (here specified by user).
            ! NB2: vertical extrapolation of levels with area_exchg=0 is done here using jk_exchg.
            zmelt(kbasin,jk) = rn_isfpar_Kcoeff * zcoef * zzztf2s(kbasin,jk_exchg(kbasin,jk)) * risf_par_area(kbasin,jk)
            !
            ! 2D net fresh water flux due to ice shelf melting ( > 0 from isf to oce) [kg m^-2 s^-1]
            ! (will be redistributed between zmin and zmax by subroutines isf_hdiv_mlt and tra_isf_mlt)
            ! NB: here, jk_exchg is used to put the meltwater in ocean cells of the exchange zone.
            pqfwf(:,:) = pqfwf(:,:) + zmelt(kbasin,jk) * mskisf_exchg(:,:,kbasin) &
            &                         * tmask(:,:,jk_exchg(kbasin,jk)) / area_effec(kbasin,jk_exchg(kbasin,jk))
            !
          END DO
        ELSE
          !
          zzztf2s(kbasin,:) = zfillvalue 
          zmelt(kbasin,:) = zfillvalue
        ENDIF
        !
      ENDDO ! kbasin
      !
      ! output per basin and per vertical level (to possibly redistribute per depth to finer-scale ice shelf draft):
      CALL iom_put('melt_bas_isf_par', zmelt )    ! parameterised melt [kg/s] -> to be redistributed to ice-sheet model
      CALL iom_put('tf2s_bas_isf_par', zzztf2s )  ! mean ( TF*|TF|*Sloc ) [degC^2 1.e-3]    -> debug/check diagnostic
      !
      ! Associated heat flux and heat content
      pqoce(:,:) = - pqfwf(:,:) * rLfusisf          ! ocean/ice shelf latent heat flux ( > 0 from isf to oce)
      pqhc (:,:) =   pqfwf(:,:) * ztfrz(:,:) * rcp  ! heat content flux                ( > 0 from isf to oce)
      !
   END SUBROUTINE isfpar_mlt_quad_loc

   SUBROUTINE isfpar_mlt_oasis(kt, Kmm, pqhc , pqoce, pqfwf )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE isfpar_mlt_oasis  ***
      !!
      !! ** Purpose    : scale the fwf read from input file by the total amount received by the sbccpl interface
      !!
      !! ** Purpose    : - read ice shelf melt from forcing file and scale it by the input file total amount => pattern
      !!                 - compute total amount of fwf given by sbccpl (fwfisf_oasis)
      !!                 - scale fwf and compute heat fluxes
      !!
      !!---------------------------------------------------------------------
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: pqhc, pqoce, pqfwf  ! heat content, latent heat and fwf fluxes
      !!-------------------------- IN  -------------------------------------
      INTEGER                     , INTENT(in   ) :: kt                  ! current time step
      INTEGER                     , INTENT(in   ) :: Kmm                 !  ocean time level index
      !!--------------------------------------------------------------------
      INTEGER                           :: jk                            ! loop index
      REAL(wp)                          :: zfwf_fld, zfwf_oasis          ! total fwf in the forcing fields (pattern) and from the cpl interface (amount)
      REAL(wp), DIMENSION(jpi,jpj)      :: ztfrz                         ! tbl freezing temperature
      REAL(wp), DIMENSION(jpi,jpj)      :: zfwf                          ! 2d fwf map after scaling
      REAL(wp), DIMENSION(jpi,jpj,jpk)  :: ztfrz3d
      !!--------------------------------------------------------------------
      !
      ! 0. ------------Read specified runoff
      CALL fld_read ( kt, 1, sf_isfpar_fwf   )
      !
      ! 1. ------------Mean freezing point (needed for heat content flux)
      DO jk = 1,jpk
         CALL eos_fzp(ts(:,:,jk,jp_sal,Kmm), ztfrz3d(:,:,jk), gdept(:,:,jk,Kmm))
      END DO
      CALL isf_tbl(Kmm, ztfrz3d, ztfrz, 'T', misfkt_par, rhisf_tbl_par, misfkb_par, rfrac_tbl_par )
      !
      ! 2. ------------Scale isf melt pattern with total amount from oasis
      ! ice shelf 2d map of fwf from isf to oce
      zfwf(:,:) = sf_isfpar_fwf(1)%fnow(:,:,1)
      !
      ! compute glob sum from input file
      ! (PM) should we consider delay sum as in fwb ? (it will offset by 1 time step if I understood well)
      zfwf_fld = glob_sum('isfcav_mlt', e1e2t(:,:) * zfwf(:,:))
      !
      ! compute glob sum from atm->oce ice shelf fwf
      ! (PM) should we consider delay sum as in fwb ?
      zfwf_oasis = glob_sum('isfcav_mlt', e1e2t(:,:) * fwfisf_oasis(:,:))
      !
      ! scale fwf
      zfwf(:,:) = zfwf(:,:) * zfwf_oasis / zfwf_fld
      ! 
      ! 3. -----------Define fwf and qoce
      ! ocean heat flux is assume to be equal to the latent heat
      pqfwf(:,:) =   zfwf(:,:)                         ! fwf                ( > 0 from isf to oce)
      pqoce(:,:) = - pqfwf(:,:) * rLfusisf             ! ocean heat flux    ( > 0 from isf to oce) (assumed to be the latent heat flux)
      pqhc (:,:) =   pqfwf(:,:) * ztfrz(:,:) * rcp     ! heat content flux  ( > 0 from isf to oce)
      !
      CALL iom_put('isftfrz_par', ztfrz )
      !
   END SUBROUTINE isfpar_mlt_oasis

END MODULE isfparmlt
