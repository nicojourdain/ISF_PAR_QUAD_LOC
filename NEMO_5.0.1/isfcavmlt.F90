MODULE isfcavmlt
   !!======================================================================
   !!                       ***  MODULE  isfcavmlt  ***
   !! ice shelf module :  update surface ocean boundary condition under ice
   !!                   shelves
   !!======================================================================
   !! History :  4.0  !  2019-09  (P. Mathiot) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isfcav_mlt    : compute or read ice shelf fwf/heat fluxes from isf 
   !!                   to oce
   !!----------------------------------------------------------------------
   USE isf_oce                  ! ice shelf

   USE par_oce                            ! ocean space and time domain
   USE dom_oce                            ! ocean space and time domain
   USE phycst , ONLY: rcp, rho0, rho0_rcp ! physical constants

   USE in_out_manager              ! I/O manager
   USE fldread    , ONLY: fld_read, FLD, FLD_N !
   USE lib_fortran, ONLY: glob_2Dsum
   USE lib_mpp    , ONLY: ctl_stop

   IMPLICIT NONE
   PRIVATE

   PUBLIC   isfcav_mlt
   !
   TYPE(FLD)  , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: sf_isfcav_fwf

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE isfcav_mlt( kt, l_converged, pgt, pgs , pttbl, pstbl, &
      &                                    pqhc, pqoce, pfwf, ptfrz )
      !!----------------------------------------------------------------------
      !!
      !!                          ***  ROUTINE isfcav_mlt  ***
      !!
      !! ** Purpose    : compute or read ice shelf fwf/heat fluxes in the ice shelf cavity
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      LOGICAL , DIMENSION(A2D(0)), INTENT(in   ) ::   l_converged         ! true when melting converges (per grid point)
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   pgt, pgs            ! gamma t and gamma s
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   pttbl, pstbl        ! top boundary layer tracer
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pqoce, pfwf   ! heat and fwf fluxes
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   ptfrz               ! tbl freezing temperature
      !!---------------------------------------------------------------------
      !
      ! compute latent heat and melt (2d)
      SELECT CASE ( cn_isfcav_mlt )
      CASE ( 'spe' )   ! ice shelf melt specified (read input file, and heat fluxes derived from
         !
         CALL isfcav_mlt_spe  ( kt, l_converged,                          pqhc, pqoce, pfwf, ptfrz )
         !
      CASE ( '2eq' )   !  ISOMIP  formulation (2 equations) for volume flux (Hunter et al., 2006)
         !
         CALL isfcav_mlt_2eq  (     l_converged, pgt,       pttbl,        pqhc, pqoce, pfwf, ptfrz )
         !
      CASE ( '3eq' )   ! ISOMIP+ formulation (3 equations) for volume flux (Asay-Davis et al., 2015)
         !
         CALL isfcav_mlt_3eq  (     l_converged, pgt, pgs , pttbl, pstbl, pqhc, pqoce, pfwf, ptfrz )
         !
      CASE ( 'oasis' ) ! fwf pass trough oasis
         !
         CALL isfcav_mlt_oasis( kt, l_converged,                          pqhc, pqoce, pfwf, ptfrz )
         !
      CASE DEFAULT
         CALL ctl_stop('STOP', 'unknown isf melt formulation : cn_isfcav (should not see this)')
      END SELECT
      !
   END SUBROUTINE isfcav_mlt


   SUBROUTINE isfcav_mlt_spe( kt, l_converged,                          &  ! <<== in
      &                                        pqhc, pqoce, pfwf, ptfrz )  ! ==>> inout
      !!----------------------------------------------------------------------
      !!
      !!                          ***  ROUTINE isfcav_mlt_spe  ***
      !!
      !! ** Purpose    : - read ice shelf melt from forcing file
      !!                 - compute ocea-ice heat flux (assuming it is equal to latent heat)
      !!                 - compute heat content flux
      !!---------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kt                 ! current time step
      LOGICAL , DIMENSION(A2D(0)), INTENT(in   ) ::   l_converged        ! true when melting converges (per grid point)
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pqoce, pfwf  ! heat content, latent heat and fwf fluxes
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ptfrz              ! tbl freezing temp
      !!
      INTEGER ::   ji, jj     ! dummy loop indices
      !!--------------------------------------------------------------------
      !
      ! read input file of fwf (from isf to oce; ie melt)
      CALL fld_read ( kt, 1, sf_isfcav_fwf )
      !
      ! define fwf and qoce (ocean heat flux is assume to be equal to the latent heat)
      DO_2D( 0, 0, 0, 0 )
         IF ( .NOT. l_converged(ji,jj) ) THEN
            pfwf (ji,jj) =   sf_isfcav_fwf(1)%fnow(ji,jj,1)   * mskisf_cav(ji,jj)   ! fwf                ( > 0 from isf to oce)
            pqoce(ji,jj) = - pfwf(ji,jj) * rLfusisf           * mskisf_cav(ji,jj)   ! ocean heat flux    ( > 0 from isf to oce)
            pqhc (ji,jj) =   pfwf(ji,jj) * ptfrz(ji,jj) * rcp * mskisf_cav(ji,jj)   ! heat content flux  ( > 0 from isf to oce)
         ENDIF
      END_2D
      !
   END SUBROUTINE isfcav_mlt_spe

   SUBROUTINE isfcav_mlt_2eq( l_converged, pgt , pttbl,             &  ! <<== in
      &                                    pqhc, pqoce, pfwf, ptfrz )  ! ==>> inout
      !!----------------------------------------------------------------------
      !!
      !!                          ***  ROUTINE isfcav_mlt_2eq  ***
      !!
      !! ** Purpose    : Compute ice shelf fwf/heqt fluxes using ISOMIP formulation (Hunter et al., 2006)
      !!
      !! ** Method     : The ice shelf melt latent heat is defined as being equal to the ocean/ice heat flux.
      !!                 From this we can derived the fwf, ocean/ice heat flux and the heat content flux as being :
      !!                    fwf  = Gammat * rho0 * Cp * ( Tw - Tfrz ) / Lf 
      !!                   qhoce = qlat
      !!                   qhc   = fwf * Cp * Tfrz
      !!
      !! ** Reference  : Hunter,  J.  R.:  Specification  for  test  models  of  ice  shelf  cavities,  
      !!                 Tech.  Rep.  June,  Antarctic  Climate  &  Ecosystems  Cooperative  Research  Centre,  available  at:  
      !!                 http://staff.acecrc.org.au/~bkgalton/ISOMIP/test_cavities.pdf (last access: 21 July 2016), 2006.
      !!
      !!--------------------------------------------------------------------
      LOGICAL , DIMENSION(A2D(0)), INTENT(in   ) ::   l_converged        ! true when melting converges (per grid point)
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   pgt                ! temperature exchange coeficient
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   pttbl              ! temperature and salinity in top boundary layer
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pqoce, pfwf  ! hean content, ocean-ice heat and fwf fluxes
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ptfrz              ! tbl freezing temp
      !!
      INTEGER  ::   ji, jj     ! dummy loop indices
      !!--------------------------------------------------------------------
      !
      DO_2D( 0, 0, 0, 0 )
         !
         ! compute ocean-ice heat flux and then derive fwf assuming that ocean heat flux equal latent heat
         IF ( .NOT. l_converged(ji,jj) ) THEN
            pfwf (ji,jj) =   pgt (ji,jj) * rho0_rcp * ( pttbl(ji,jj) - ptfrz(ji,jj) ) / rLfusisf * mskisf_cav(ji,jj)  ! fresh water flux  ( > 0 from isf to oce)
            pqoce(ji,jj) = - pfwf(ji,jj) * rLfusisf                                              * mskisf_cav(ji,jj)  ! ocea-ice flux     ( > 0 from isf to oce)
            pqhc (ji,jj) =   pfwf(ji,jj) * ptfrz(ji,jj) * rcp                                    * mskisf_cav(ji,jj)  ! heat content flux ( > 0 from isf to oce)
            !
         ENDIF
      END_2D
      !
   END SUBROUTINE isfcav_mlt_2eq

   SUBROUTINE isfcav_mlt_3eq( l_converged, pgt, pgs , pttbl, pstbl, &  ! <<== in
      &                                    pqhc, pqoce, pfwf, ptfrz )  ! ==>> inout
      !!----------------------------------------------------------------------
      !!
      !!                          ***  ROUTINE isfcav_mlt_3eq  ***
      !!
      !! ** Purpose    : Compute ice shelf fwf/heqt fluxes using the 3 equation formulation 
      !!
      !! ** Method     : The melt rate is determined considering the heat balance, the salt balance
      !!                 at the phase change interface and a linearisation of the equation of state.
      !!
      !! ** Reference  : - Holland, D. M. and Jenkins, A.,
      !!                   Modeling Thermodynamic Ice-Ocean Interactions at the Base of an Ice Shelf,
      !!                   J. Phys. Oceanogr., 29, 1999.
      !!                 - Asay-Davis, X. S., Cornford, S. L., Durand, G., Galton-Fenzi, B. K., Gladstone, 
      !!                   R. M., Gudmundsson, G. H., Hattermann, T., Holland, D. M., Holland, D., Holland, 
      !!                   P. R., Martin, D. F., Mathiot, P., Pattyn, F., and Seroussi, H.:
      !!                   Experimental design for three interrelated marine ice sheet and ocean model intercomparison projects: 
      !!                   MISMIP v. 3 (MISMIP +), ISOMIP v. 2 (ISOMIP +) and MISOMIP v. 1 (MISOMIP1), 
      !!                   Geosci. Model Dev., 9, 2471-2497, https://doi.org/10.5194/gmd-9-2471-2016, 2016. 
      !!
      !!--------------------------------------------------------------------
      LOGICAL , DIMENSION(A2D(0)), INTENT(in   ) ::   l_converged        ! true when melting converges (per grid point)
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   pgt  , pgs         ! heat/salt exchange coeficient
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   pttbl, pstbl       ! mean temperature and salinity in top boundary layer
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pqoce, pfwf  ! latent heat and fwf fluxes
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   ptfrz              ! tbl freezing temp
      !!
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::   zeps1, zeps2, zeps3, zeps4, zeps6, zeps7        ! dummy local scalar for quadratic equation resolution
      REAL(wp) ::   zaqe, zbqe, zcqe, zaqer, zdis, zsfrz, zcfac     ! dummy local scalar for quadratic equation resolution
      !!--------------------------------------------------------------------
      !
      ! compute upward heat flux zhtflx and upward water flux zwflx
      ! Resolution of a 3d equation from equation 24, 25 and 26 (note conduction through the ice has been added to Eq 24)
      DO_2D( 0, 0, 0, 0 )
         !
         IF ( .NOT. l_converged(ji,jj) ) THEN
            ! compute coeficient to solve the 2nd order equation
            zeps1 = rho0_rcp * pgt(ji,jj)
            zeps2 = rLfusisf * rho0 * pgs(ji,jj)
            zeps3 = rhoisf * rcpisf * rkappa / MAX(risfdep(ji,jj),1.e-20)
            zeps4 = risf_lamb2 + risf_lamb3 * risfdep(ji,jj)
            zeps6 = zeps4 - pttbl(ji,jj)
            zeps7 = zeps4 - rtsurf
            !
            ! solve the 2nd order equation to find zsfrz
            zaqe  = risf_lamb1 * (zeps1 + zeps3)
            zaqer = 0.5_wp / MIN(zaqe,-1.e-20)
            zbqe  = zeps1 * zeps6 + zeps3 * zeps7 - zeps2
            zcqe  = zeps2 * pstbl(ji,jj)
            zdis  = zbqe * zbqe - 4.0_wp * zaqe * zcqe
            !
            ! Presumably zdis can never be negative because gammas is very small compared to gammat
            zsfrz=(-zbqe - SQRT(zdis)) * zaqer
            IF ( zsfrz < 0.0_wp ) zsfrz=(-zbqe + SQRT(zdis)) * zaqer  ! check this if this if is needed
            !
            ! compute t freeze (eq. 25)
            ptfrz(ji,jj) = zeps4 + risf_lamb1 * zsfrz
            !
            ! compute the upward water and heat flux (eq. 24 and eq. 26)
            pfwf (ji,jj) = - rho0     * pgs(ji,jj) * ( zsfrz - pstbl(ji,jj) ) / MAX(zsfrz,1.e-20) * mskisf_cav(ji,jj) ! fresh water flux    ( > 0 from isf to oce)
            pqoce(ji,jj) = - rho0_rcp * pgt(ji,jj) * ( pttbl(ji,jj) - ptfrz(ji,jj) )              * mskisf_cav(ji,jj) ! ocean-ice heat flux ( > 0 from isf to oce)
            pqhc (ji,jj) =   rcp      * pfwf(ji,jj) * ptfrz(ji,jj)                                * mskisf_cav(ji,jj) ! heat content   flux ( > 0 from isf to oce)
            !
         ENDIF
         !
      END_2D
      !
   END SUBROUTINE isfcav_mlt_3eq

   SUBROUTINE isfcav_mlt_oasis( kt, l_converged,                          &  ! <<== in
      &                                          pqhc, pqoce, pfwf, ptfrz )  ! ==>> inout
      !!----------------------------------------------------------------------
      !!                          ***  ROUTINE isfcav_mlt_oasis  ***
      !!
      !! ** Purpose    : scale the fwf read from input file by the total amount received by the sbccpl interface
      !!
      !! ** Purpose    : - read ice shelf melt from forcing file => pattern
      !!                 - total amount of fwf is given by sbccpl (fwfisf_cpl)
      !!                 - scale fwf and compute heat fluxes
      !!
      !!---------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kt                 ! current time step
      LOGICAL , DIMENSION(A2D(0)), INTENT(in   ) ::   l_converged        ! true when melting converges (per grid point)
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pqhc, pqoce, pfwf  ! heat content, latent heat and fwf fluxes
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ptfrz              ! tbl freezing temp
      !!
      INTEGER  ::   ji, jj
      REAL(wp), DIMENSION(A2D(0),2) ::   ztmp
      REAL(wp), DIMENSION(2)        ::   zbg
      !!--------------------------------------------------------------------
      !
      ! read input file of fwf from isf to oce
      CALL fld_read ( kt, 1, sf_isfcav_fwf )
      !
      DO_2D( 0, 0, 0, 0 )
         ! ice shelf 2d map
         IF ( .NOT. l_converged(ji,jj) )   pfwf(ji,jj) = sf_isfcav_fwf(1)%fnow(ji,jj,1)
         ztmp(ji,jj,1) = e1e2t(ji,jj) * pfwf(ji,jj)
         ztmp(ji,jj,2) = e1e2t(ji,jj) * fwfisf_oasis(ji,jj)
      END_2D
      !
      ! compute glob sum from input file and from atm->oce ice shelf fwf
      ! (PM) should consider delay sum as in fwb (1 time step offset if I well understood)
      zbg(:) = glob_2Dsum( 'isfcav_mlt', ztmp, cdelay = 'fwfisf' )
      !
      ! compute the upward water and heat flux (ocean heat flux is assume to be equal to the latent heat)
      DO_2D( 0, 0, 0, 0 )
         IF ( .NOT. l_converged(ji,jj) ) THEN
            pfwf (ji,jj) =   pfwf(ji,jj) * zbg(2) / zbg(1)    * mskisf_cav(ji,jj)   ! scale fwf
            pqoce(ji,jj) = - pfwf(ji,jj) * rLfusisf           * mskisf_cav(ji,jj)   ! ocean heat flux    ( > 0 from isf to oce)
            pqhc (ji,jj) =   pfwf(ji,jj) * ptfrz(ji,jj) * rcp * mskisf_cav(ji,jj)   ! heat content flux  ( > 0 from isf to oce)
         ENDIF
      END_2D
      !
   END SUBROUTINE isfcav_mlt_oasis

END MODULE isfcavmlt
