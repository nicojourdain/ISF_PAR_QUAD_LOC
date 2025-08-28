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
   USE isftbl                   ! ice shelf depth average

   USE par_oce                  ! ocean space and time domain
   USE dom_oce                  ! ocean space and time domain
   USE oce    , ONLY: ts        ! ocean dynamics and tracers
   USE phycst , ONLY: rcp, rho0 ! physical constants
   USE eosbn2 , ONLY: eos_fzp   ! equation of state

   USE xios                                    ! Issue with the output at format (nbasin,jpk) with iom_put (tilling issue)
   USE in_out_manager                          ! I/O manager
   USE iom        , ONLY: iom_put, iom_use     ! I/O library
   USE fldread    , ONLY: fld_read, FLD, FLD_N !
   USE lib_fortran, ONLY: local_2Dsum, glob_2Dsum
   USE lib_mpp    , ONLY: mpp_sum, mpp_max, ctl_stop    !
   USE lbclnk

   IMPLICIT NONE

   PRIVATE

   PUBLIC  isfpar_mlt 

   TYPE(FLD), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     :: sf_isfpar_fwf

   !! * Substitutions
#  include "domzgr_substitute.h90"   
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE isfpar_mlt( kt, Kmm, ptfrz, ptavg, pqhc, pqoce, pfwf )
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
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ptfrz, ptavg       ! tbl freezing and averaged temperatures
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pfwf, pqoce, pqhc  ! fresh water, ice-ocean heat and heat content fluxes
      !!---------------------------------------------------------------------
      !
      ! Choose among the available ice shelf parametrisation
      SELECT CASE ( cn_isfpar_mlt )
      CASE ( 'spe' )      ! specified runoff in depth (Mathiot et al., 2017)
         !
         CALL isfpar_mlt_spe(   kt, Kmm, ptfrz, pqhc, pqoce, pfwf )
         !
      CASE ( 'bg03' )     ! Beckmann and Goosse parametrisation 
         !
         CALL isfpar_mlt_bg03(  kt, Kmm, ptfrz, ptavg, pqhc, pqoce, pfwf )
         !
      CASE ( 'quad_loc' ) ! Quadratic local melt parameterisation (Burgard et al., TC, 2022)
         !
         CALL isfpar_mlt_quad_loc(kt, Kmm, pqhc, pqoce, pfwf)
         !
      CASE ( 'oasis' )    ! Climate model
         !
         CALL isfpar_mlt_oasis( kt, Kmm, ptfrz, pqhc, pqoce, pfwf )
         !
      CASE DEFAULT
         CALL ctl_stop('STOP', 'unknown isf melt formulation : cn_isfpar (should not see this)')
      END SELECT
      !
   END SUBROUTINE isfpar_mlt


   SUBROUTINE isfpar_mlt_spe( kt, Kmm, ptfrz, pqhc, pqoce, pfwf )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isfpar_mlt_spe  ***
      !!
      !! ** Purpose : prescribed ice shelf melting in case ice shelf cavities are closed.
      !!              data read into a forcing files.
      !!
      !!--------------------------------------------------------------------
      INTEGER,  INTENT(in) ::   kt
      INTEGER,  INTENT(in) ::   Kmm    !  ocean time level index
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ptfrz              ! tbl freezing temp
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pfwf, pqoce  ! fresh water and ice-ocean heat fluxes
      !!
      INTEGER ::   ji, jj     ! dummy loop indices
      !!--------------------------------------------------------------------
      !
      ! Read specified fwf from isf to oce
      CALL fld_read ( kt, 1, sf_isfpar_fwf )
      !
      DO_2D( 0, 0, 0, 0 )
         pfwf(ji,jj) =   sf_isfpar_fwf(1)%fnow(ji,jj,1)    * mskisf_par(ji,jj)  ! fresh water flux from the isf (fwfisf <0 mean melting)       ( > 0 from isf to oce)
         pqoce(ji,jj) = - pfwf(ji,jj) * rLfusisf           * mskisf_par(ji,jj)  ! ocean/ice shelf flux assume to be equal to latent heat flux  ( > 0 from isf to oce)
         pqhc (ji,jj) =   pfwf(ji,jj) * ptfrz(ji,jj) * rcp * mskisf_par(ji,jj)  ! heat content flux                                            ( > 0 from isf to oce)
      END_2D
      !
      !
   END SUBROUTINE isfpar_mlt_spe

   SUBROUTINE isfpar_mlt_bg03( kt, Kmm, ptfrz, ptavg, pqhc, pqoce, pfwf )
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
      INTEGER,  INTENT(in) ::   kt
      INTEGER,  INTENT(in) ::   Kmm    !  ocean time level index
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ptfrz, ptavg       ! tbl freezing and averaged temp
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pfwf, pqoce  ! fresh water and ice-ocean heat fluxes
      !!
      INTEGER ::   ji, jj     ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      ! Net heat flux and fresh water flux due to the ice shelf
      DO_2D( 0, 0, 0, 0 )
         pfwf (ji,jj) =  rho0 * rcp * rn_isfpar_bg03_gt0 * risfLeff(ji,jj) * e1t(ji,jj) * ( ptavg(ji,jj) - ptfrz(ji,jj) ) &
            &                 * r1_e1e2t(ji,jj) / rLfusisf * mskisf_par(ji,jj)   ! ( > 0 from isf to oce)
         pqoce(ji,jj) = - pfwf(ji,jj) * rLfusisf           * mskisf_par(ji,jj)   ! ocean/ice shelf flux assume to be equal to latent heat flux  ( > 0 from isf to oce)
         pqhc (ji,jj) =   pfwf(ji,jj) * ptfrz(ji,jj) * rcp * mskisf_par(ji,jj)   ! heat content flux                                            ( > 0 from isf to oce)
      END_2D
      !
   END SUBROUTINE isfpar_mlt_bg03

   SUBROUTINE isfpar_mlt_quad_loc(kt, Kmm, pqhc, pqoce, pfwf)
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
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) :: pqhc, pfwf, pqoce  ! fresh water and ice-ocean heat fluxes
      !!-------------------------- IN  -------------------------------------
      INTEGER,  INTENT(in) :: kt
      INTEGER,  INTENT(in) :: Kmm    !  ocean time level index
      !!--------------------------------------------------------------------
      COMPLEX(dp), DIMENSION(nbasins_glo*jpk) :: ctmp  ! 1d version of ztf2s for the mppsum
      COMPLEX(dp), DIMENSION(nbasins_glo,jpk) :: ctf2s ! TF*|TF|*Sloc for interfacial cells  [degC^2 1.e-3]
      REAL(wp)   , DIMENSION(nbasins_glo*jpk) :: ztmp  ! 1d version for the mppsum
      REAL(wp), DIMENSION(A2D(0))         :: ztfrz    ! mean freezing temperature in interfacial water columns [degC]
      REAL(wp), DIMENSION(A2D(0),jpk)     :: ztfrz3d  ! 3d freezing temperature [degC]
      REAL(wp), DIMENSION(A2D(0),jpk)     :: ztmp3d   ! tmp array for e3 or depth for eos_fzp and isf_tbl_avg
      REAL(wp), DIMENSION(A2D(0))         :: ztf      ! 2d thermal forcing [degC]
      REAL(wp), DIMENSION(A2D(0))         :: ztf2sa   ! 2d thermal forcing ^2 * s * area [degC^2 1e-3 m^2]
      REAL(wp), DIMENSION(nbasins_glo,jpk) :: zmelt    ! Parameterised melt per basin and per level [kg s^-1]
      REAL(wp), DIMENSION(nbasins_glo,jpk) :: ztf2s    ! TF*|TF|*Sloc for interfacial cells  [degC^2 1.e-3]
      INTEGER  :: ji, jj, jk, jke, jb_loc, jb_glo      ! dummy loop indices
      REAL(wp) :: zcoef, zdum                          ! Constant coefficient used in the param [kg m^-2 s^-1 degC^-2 1.e3]
      !!----------------------------------------------------------------------
      !
      ! 1: freezing point definition
      ! ----------------------------
      ! 
      ! ztfrz3d is the freezing temperature at all levels:
      DO_3D( 0 ,0, 0, 0, 1, jpk )
         ztmp3d(ji,jj,jk) = gdept(ji,jj,jk,Kmm)
      END_3D
      CALL eos_fzp(ts(A2D(0),:,jp_sal,Kmm), ztfrz3d(:,:,:), ztmp3d(:,:,:), 0)

      ! ztfrz is the mean freezing temperature in the interfacial water columns, i.e., between
      ! specified zmin and zmax (only used for the heat content flux):
      !
      DO_3D( 0 ,0, 0, 0, 1, jpk )
         ztmp3d(ji,jj,jk) = e3t(ji,jj,jk,Kmm)
      END_3D
      CALL isf_tbl_avg(misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par, ztmp3d, ztfrz3d, ztfrz )

      DO_3D( 0 ,0, 0, 0, 1, jpk )
         ztmp3d(ji,jj,jk) = gdept(ji,jj,jk,Kmm)
      END_3D
      !
      ! 2: Main computation
      ! ------------------
      !
      ! 2.1: Thermal forcing computation
      ! --------------------------------
      ! Calculate ztftfs3d as TF*|TF|*Sloc*e1t*e2t (where TF = thermal forcing) [degC^2 1.e-3 m^2]:
      IF ( ANY(ln_exchg(:)) ) THEN
         ! Main loops: outer on jk for better memory access on ctf2s(:, jk)
         DO jk = 1, jpk
            ! Compute the product once into ztf
            DO_2D( 0 ,0, 0, 0)
               ztf(ji,jj) = ( ts(ji,jj,jk,jp_tem,Kmm) - ztfrz3d(ji,jj,jk) )
               ztf2sa(ji,jj) = ztf(ji,jj) * ABS(ztf(ji,jj)) * ts(ji,jj,jk,jp_sal,Kmm) * e1e2t(ji,jj) * tmask(ji,jj,jk)
            END_2D
      
            ! Perform local sum if necessary
            DO jb_glo = 1, nbasins_glo
               IF (ln_exchg(jb_glo)) THEN
                  jb_loc = idx_basin_glo_to_loc(jb_glo)
                  ctf2s(jb_glo, jk) = local_2Dsum(ztf2sa(:,:) * mskisf_exchg(:,:,jb_loc) * r1_area_exchg(jb_loc, jk))
               ELSE
                  ctf2s(jb_glo, jk) = CMPLX(0.e0, 0.e0, dp)
               END IF
            END DO
         END DO
      ELSE
         ctf2s(:,:) = CMPLX(0.e0, 0.e0, dp)
      END IF
      ! 
      ! Average profile of TF*|TF|*Sloc in the interfacial ocean grid cells [degC^2 1.e-3]
      ! because mpp_sum is only 0d or 1d, we need to reshape
      ctmp(:) = RESHAPE(ctf2s, [nbasins_glo * jpk] )
      CALL mpp_sum( 'isfparmlt', ctmp(:) )
      ztf2s(:,:) = REAL(RESHAPE(ctmp(:), [nbasins_glo, jpk]), wp)
      !
      !
      ! 2.2: Melt computation
      ! --------------------------------
      !
      pfwf(:,:) = 0._wp
      zmelt(:,:) = -HUGE(1._wp)
      IF ( ANY(ln_exchg(:)) ) THEN
         DO jk = 1, jpk
            DO jb_glo = 1, nbasins_glo
               IF (ln_exchg(jb_glo)) THEN
                  jb_loc = idx_basin_glo_to_loc(jb_glo)
                  jke = jk_exchg(jb_loc, jk)
                  ! Melt per basin per vertical level [kg s^-1] :
                  ! NB2: vertical extrapolation of levels with area_exchg=0 is done here using jk_exchg.
                  zmelt(jb_glo,jk) = rtf2s_to_melt_loc(jb_loc,jk) * ztf2s(jb_glo,jke)
                  pfwf(:,:) = pfwf(:,:) + zmelt(jb_glo,jk) * mskisf_exchg(:,:,jb_loc) * tmask(A2D(0),jke) * r1_area_exchg(jb_loc,jke)
               END IF
            END DO
         END DO
      END IF
      !
      ! 3: Diagnostics
      ! output per basin and per vertical level (to possibly redistribute per depth to finer-scale ice shelf draft):
      ! use mpp_max, as the same basin can be spread on multiple subdomain and in this case zmelt is already the total and the same on each subdomain 
      IF ( iom_use('melt_bas_isf_par') ) THEN
         ztmp(:) = RESHAPE(zmelt, [nbasins_glo * jpk] )
         CALL mpp_max( 'isfparmlt', ztmp(:) )
         zmelt(:,:) = RESHAPE(ztmp(:), [nbasins_glo, jpk])
         CALL xios_send_field( 'melt_bas_isf_par', zmelt(:,:) ) ! parameterised melt [kg/s] -> to be redistributed to ice-sheet model
      END IF
      
      IF ( iom_use('tf2s_bas_isf_par') ) THEN
         CALL xios_send_field( 'tf2s_bas_isf_par', ztf2s(:,:) ) ! mean ( TF*|TF|*Sloc ) [degC^2 1.e-3]    -> debug/check diagnostic
      END IF
      !
      ! 4: Output of the subroutine
      ! Associated heat flux and heat content
      pqoce(:,:) = - pfwf(:,:) * rLfusisf          ! ocean/ice shelf latent heat flux ( > 0 from isf to oce)
      pqhc (:,:) =   pfwf(:,:) * ztfrz(:,:) * rcp  ! heat content flux                ( > 0 from isf to oce)
      !
   END SUBROUTINE isfpar_mlt_quad_loc

   SUBROUTINE isfpar_mlt_oasis(kt, Kmm, ptfrz, pqhc , pqoce, pfwf )
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
      INTEGER                    , INTENT(in   ) ::   kt                 ! current time step
      INTEGER                    , INTENT(in   ) ::   Kmm                !  ocean time level index
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ptfrz              ! tbl freezing temp
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pqoce, pfwf  ! heat content, latent heat and fwf fluxes
      !!
      INTEGER  ::   ji, jj           ! dummy loop indices
      REAL(wp), DIMENSION(A2D(0),2) ::   ztmp
      REAL(wp), DIMENSION(2)        ::   zbg
      !!--------------------------------------------------------------------
      !
      ! Read specified fwf from isf to oce
      CALL fld_read ( kt, 1, sf_isfpar_fwf )
      !
      DO_2D( 0, 0, 0, 0 )
         ! ice shelf 2d map
         pfwf(ji,jj) = sf_isfpar_fwf(1)%fnow(ji,jj,1)
         ztmp(ji,jj,1) = e1e2t(ji,jj) * pfwf(ji,jj)
         ztmp(ji,jj,2) = e1e2t(ji,jj) * fwfisf_oasis(ji,jj)
      END_2D
      !
      ! compute glob sum from input file and from atm->oce ice shelf fwf
      zbg = glob_2Dsum( 'isfpar_mlt', ztmp, cdelay = 'isfmlt_tag' )
      !
      ! Define fwf and qoce
      ! ocean heat flux is assume to be equal to the latent heat
      DO_2D( 0, 0, 0, 0 )
         pfwf (ji,jj) =   pfwf(ji,jj) * zbg(2) / zbg(1)    * mskisf_par(ji,jj)   ! scale fwf ( > 0 from isf to oce)
         pqoce(ji,jj) = - pfwf(ji,jj) * rLfusisf           * mskisf_par(ji,jj)   ! ocean heat flux    ( > 0 from isf to oce) (assumed to be the latent heat flux)
         pqhc (ji,jj) =   pfwf(ji,jj) * ptfrz(ji,jj) * rcp * mskisf_par(ji,jj)   ! heat content flux  ( > 0 from isf to oce)
      END_2D
      !
   END SUBROUTINE isfpar_mlt_oasis

END MODULE isfparmlt
