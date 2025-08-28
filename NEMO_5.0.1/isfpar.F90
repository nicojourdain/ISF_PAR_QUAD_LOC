MODULE isfpar
   !!======================================================================
   !!                       ***  MODULE  isfpar  ***
   !! ice shelf module :  update ocean boundary condition under ice
   !!                   shelf
   !!======================================================================
   !! History :  3.2  !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X  !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4  !  2013-03  (P. Mathiot) Merging + parametrization
   !!            4.1  !  2019-09  (P. Mathiot) Restructuration
   !!            4.2  !  2021-05  (C. Ethe   ) Test and fix oasis case
   !!            4.2  !  2025-04  (N. Jourdain) Interactive quadratic param.
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isfpar       : compute ice shelf melt using a prametrisation of ice shelf cavities
   !!----------------------------------------------------------------------
   USE isf_oce        ! ice shelf
   !
   USE isfrst   , ONLY: isfrst_read               ! ice shelf restart read/write subroutine
   USE isftbl                                     ! ice shelf top boundary layer properties subroutine
   USE isfparmlt, ONLY: sf_isfpar_fwf, isfpar_mlt ! ice shelf melt formulation subroutine
   USE isfdiags , ONLY: isf_diags_flx             ! ice shelf diags subroutine
   USE isfutils , ONLY: debug, read_2dcstdta      ! ice shelf debug subroutine
   !
   USE oce      , ONLY: ts             ! ocean dynamics and tracers
   USE dom_oce                         ! ocean space and time domain
   USE par_oce                         ! ocean space and time domain
   USE phycst                          ! physical constants
   USE eosbn2   , ONLY: eos_fzp        ! equation of state
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE fldread        ! read input field at current time step
   USE lib_fortran, ONLY: local_2Dsum, glob_2Dsum
   USE lib_mpp    , ONLY: mpp_sum, ctl_stop !
   USE lbclnk

   IMPLICIT NONE
   PRIVATE

   PUBLIC   isf_par, isf_par_init

   TYPE(FLD_N), PUBLIC :: sn_isfpar_fwf      !: information about the isf melting file to be read
   TYPE(FLD_N), PUBLIC :: sn_isfpar_zmax     !: information about the grounding line depth file to be read
   TYPE(FLD_N), PUBLIC :: sn_isfpar_zmin     !: information about the calving   line depth file to be read
   TYPE(FLD_N), PUBLIC :: sn_isfpar_Leff     !: information about the effective length     file to be read
   TYPE(FLD_N), PUBLIC :: sn_isfpar_basin    !: information about the ice-shelf basins     file to be read
   TYPE(FLD_N), PUBLIC :: sn_isfpar_area     !: information on non-resolved ice-shelf area file to be read

   !! * Substitutions   
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
 
   SUBROUTINE isf_par( kt, Kmm, ptsc, pfwf )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE isf_par ***      
      !!
      !! ** Purpose : compute the heat and fresh water due to ice shelf melting/freezing using a parametrisation 
      !!
      !! ** Comment : in isf_par and all its call tree, 
      !!              'tbl' means parametrisation layer (ie how the far field temperature/salinity is computed) 
      !!              instead of in a proper top boundary layer as at the ice shelf ocean interface
      !!              as the action to compute the properties of the tbl or the parametrisation layer are the same,
      !!              (ie average T/S over a specific depth (can be across multiple levels))
      !!              the name tbl was kept.
      !!
      !! ** Convention : all fluxes are from isf to oce
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt                                     ! ocean time step
      INTEGER, INTENT(in) ::   Kmm                                    ! ocean time level index
      REAL(wp), DIMENSION(A2D(0))     , INTENT(inout) ::   pfwf
      REAL(wp), DIMENSION(A2D(0),jpts), INTENT(inout) ::   ptsc
      !!
      INTEGER ::   ji, jj, jk
      REAL(wp), DIMENSION(A2D(0)) ::   ztfrz, ztavg                   ! tbl freezing and averaged temperatures
      REAL(wp), DIMENSION(A2D(0)) ::   zqoce, zqhc, zqlat
      REAL(wp), DIMENSION(A2D(0),jpk) ::   ztfrz3d, ztmp
      !!---------------------------------------------------------------------
      !
      ! Mean freezing point
      CALL eos_fzp( ts, Kmm, ztfrz3d, 0 )
!!st      DO_3D( 0 ,0, 0, 0, 1, jpk )
!!st         ztmp(ji,jj,jk) = gdept(ji,jj,jk,Kmm)
!!st      END_3D
!!st      CALL eos_fzp( ts(A2D(0),:,jp_sal,Kmm), ztfrz3d(:,:,:), ztmp )
      !
      DO_3D( 0 ,0, 0, 0, 1, jpk )
         ztmp(ji,jj,jk) = e3t  (ji,jj,jk,Kmm)
      END_3D
      CALL isf_tbl_avg( misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par, ztmp, ztfrz3d, & ! <<== in
         &                                                                            ztfrz  ) ! ==>> out

      ! Mean temperature (only for bg03)
      CALL isf_tbl_avg( misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par, ztmp, ts(:,:,:,jp_tem,Kmm), & ! <<== in
         &                                                                                         ztavg  ) ! ==>> out

      ! compute heat content, latent heat and melt fluxes (2d)
      CALL isfpar_mlt( kt, Kmm, ztfrz, ztavg, zqhc, zqoce, pfwf )
      !
      DO_2D( 0, 0, 0, 0 )
         ! compute latent heat flux (from isf to oce)
         zqlat(ji,jj) = - pfwf(ji,jj) * rLfusisf    ! 2d latent heat flux (W/m2)
         !
         ! set temperature content
         ptsc(ji,jj,jp_tem) = ( zqhc (ji,jj) + zqoce(ji,jj) ) * r1_rho0_rcp
      END_2D
      !
      ! output fluxes
      CALL isf_diags_flx( Kmm, misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par, 'par', pfwf, zqoce, zqlat, zqhc )
      !
      ! outputs
      CALL iom_put('isftfrz_par', ztfrz(:,:) * mskisf_par(:,:) ) ! freezing temperature
      IF( cn_isfpar_mlt == 'bg03' ) THEN
         CALL iom_put('ttbl_par',         ztavg(:,:)                * mskisf_par(:,:) )      ! ttbl
         CALL iom_put('isfthermald_par',( ztavg(:,:) - ztfrz(:,:) ) * mskisf_par(:,:) )      ! thermal driving
      ENDIF
      !
      ! debugs
      IF ( ln_isfdebug ) THEN
         IF(lwp) WRITE(numout,*)
         CALL debug('isf_par: ptsc T', ptsc (:,:,1))
         CALL debug('isf_par: ptsc S', ptsc (:,:,2))
         CALL debug( 'isfpar: qhc   ', zqhc (:,:)  )
         CALL debug( 'isfpar: qoce  ', zqoce(:,:)  )
         CALL debug( 'isfpar: fwf   ', pfwf (:,:)  )
         IF(lwp) WRITE(numout,*) ''
      END IF
      !
   END SUBROUTINE isf_par

   SUBROUTINE isf_par_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_par_init  ***
      !!
      !! ** Purpose : initialisation of the variable needed for the parametrisation of ice shelf melt
      !!
      !!------------------------------------------------------------------------------
      INTEGER                          :: ierr, inum, jb_loc, jb_glo, ji, jj, jk, indx
      REAL(wp)                         :: zepsln = 1.e-20_wp  
      INTEGER                          :: kbest, kb_used, klen, ialloc
      REAL(wp), DIMENSION(jpk)         :: ztmp
      COMPLEX(dp), DIMENSION(jpk)      :: ctmp
      INTEGER , ALLOCATABLE, DIMENSION(:)   :: klvl
      REAL(wp), DIMENSION(A2D(1))     :: ztblmax, ztblmin, zid
      CHARACTER(1024)                  :: cinfo
      !!--------------------------------------------------------------------
      ! quad_loc
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zisf_par_area_glo, zarea_exchg  ! area distribution per isf basin (from input file)
      REAL(wp) ::   zf         = 1.4e-4   ! mean Coriolis parameter [s^-1]
      REAL(wp) ::   zbeta      = 7.8e-4   ! salt contraction coefficient [1.e3]
      REAL(wp) ::   zsin_theta = 2.9e-3   ! assuming a representative "Antarctic slope" [1]
      REAL(wp) ::   zeps       = 1.0e-20
      REAL(wp) ::   zcoef
      !!------------------------------------------------------------------------------
      !
      !==============
      ! 0: allocation
      !==============
      CALL isf_alloc_par()
      !
      !==================
      ! 1: initialisation
      !==================
      DO_2D( 1, 1, 1, 1 )
         misfkt_par   (ji,jj) = 1
         misfkb_par   (ji,jj) = 1
         rhisf_tbl_par(ji,jj) = 1e-20
         rfrac_tbl_par(ji,jj) = 0.0_wp
      END_2D
      !
      ! define isf tbl tickness, top and bottom indice
      CALL read_2dcstdta( TRIM(sn_isfpar_zmax%clname), TRIM(sn_isfpar_zmax%clvar), ztblmax )
      CALL read_2dcstdta( TRIM(sn_isfpar_zmin%clname), TRIM(sn_isfpar_zmin%clvar), ztblmin )
      !
      DO_2D( 1, 1, 1, 1 )
         ! mask ice shelf parametrisation location
         ztblmax(ji,jj) = ztblmax(ji,jj) * ssmask(ji,jj)
         ztblmin(ji,jj) = ztblmin(ji,jj) * ssmask(ji,jj)
         !
         ! if param used under an ice shelf overwrite ztblmin by the ice shelf draft
         IF( risfdep(ji,jj) > 0._wp .AND. ztblmin(ji,jj) > 0._wp )   ztblmin(ji,jj) = risfdep(ji,jj)
         !
         ! enforce zmin to be above the bottom (gdepw_0 of the bottom wet cell)
         IF (ztblmin(ji,jj) > gdepw_0(ji,jj,mbkt(ji,jj))) ztblmin(ji,jj) = gdepw_0(ji,jj,mbkt(ji,jj)) 
         !
         ! ensure ztblmax <= bathy
         ztblmax(ji,jj) = MIN( ztblmax(ji,jj), bathy(ji,jj) )
      END_2D
      !
      ! initial tbl thickness
      DO_2D( 1, 1, 1, 1 )
         rhisf0_tbl_par(ji,jj) = ztblmax(ji,jj) - ztblmin(ji,jj)
      END_2D

      !
      ! define iceshelf parametrisation mask
      mskisf_par = 0
      WHERE ( rhisf0_tbl_par(A2D(0)) > 0._wp )
         mskisf_par(:,:) = 1
      END WHERE

      ! compute ktop and update ztblmin to gdepw_0 at misfkt_par
      CALL isf_tbl_ktop(ztblmin, misfkt_par) !   out: misfkt_par
      !                                      ! inout: ztblmin
      !
#if ! defined key_RK3
      !================
      ! 2: read restart
      !================
      ! MLF: read par variable from restart
      IF ( ln_rstart ) CALL isfrst_read( 'par', risf_par_tsc, fwfisf_par, risf_par_tsc_b, fwfisf_par_b )
#endif
      !
      !==========================================
      ! 3: specific allocation and initialisation (depending of scheme choice)
      !==========================================
      SELECT CASE ( TRIM(cn_isfpar_mlt) )
         !
      CASE ( 'spe' )
         !
         ALLOCATE( sf_isfpar_fwf(1), STAT=ierr )
         ALLOCATE( sf_isfpar_fwf(1)%fnow(A2D(0),1), sf_isfpar_fwf(1)%fdta(A2D(0),1,2) )
         CALL fld_fill( sf_isfpar_fwf, (/ sn_isfpar_fwf /), cn_isfdir, 'isf_par_init', 'read fresh water flux isf data', 'namisf' )
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>   ice shelf melt rate read from forcing field (cn_isfmlt_par = spe)'
         !
      CASE ( 'bg03' )
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>   bg03 parametrisation (cn_isfmlt_par = bg03)'
         !
         ! read effective length
         CALL read_2dcstdta(TRIM(sn_isfpar_Leff%clname), TRIM(sn_isfpar_Leff%clvar), risfLeff)
         risfLeff = risfLeff*1000.0_wp           !: convertion in m
         !
      CASE ( 'quad_loc' )
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>   Ice shelf melt rate calculated through'
         IF(lwp) WRITE(numout,*) '              the quadratic-local parametrisation (cn_isfmlt_par = quad_loc)'
         !
         IF ( TRIM(cn_isfpar_mlt) == 'quad_loc' ) THEN
            IF(lwp) WRITE(numout,*) '            basin file sn_isfpar_basin%name = ', TRIM(sn_isfpar_basin%clname)
            IF(lwp) WRITE(numout,*) '            basin map variable name sn_isfpar_basin% = ', TRIM(sn_isfpar_basin%clvar)
            IF(lwp) WRITE(numout,*) '            area distribution variable name sn_isfpar_area% = ', TRIM(sn_isfpar_area%clvar)
            IF(lwp) WRITE(numout,*) '            number of basin nn_isfpar_basin = ', nn_isfpar_basin
            IF(lwp) WRITE(numout,*) '            Tuning coeficient rn_isfpar_Kcoeff = ', rn_isfpar_Kcoeff
         END IF
         !
         nbasins_glo = nn_isfpar_basin
         !
         CALL isf_alloc_par_quad('glo')
         !
         ! read basin ID (2d map)
         CALL iom_open (TRIM(sn_isfpar_basin%clname), inum)
         !
         ! get basin map
         CALL iom_get  ( inum, jpdom_global , TRIM(sn_isfpar_basin%clvar), zid)
         id_basin_isfpar(:,:) = NINT(zid(A2D(0)))
         !
         ! get basin array
         CALL iom_get  ( inum, jpdom_unknown, 'basin', rbasisf_num(1:nbasins_glo))
         CALL iom_close( inum )
         !
         ! read non-resolved (i.e. parameterised) ice-shelf area [m2] per basin and per vertical level
         ALLOCATE(zisf_par_area_glo(nbasins_glo,jpk), STAT=ialloc)
         IF( ialloc > 0 ) THEN
            CALL ctl_stop( 'isfpar: unable to allocate zisf_par_area_glo' )   ;   RETURN
         ENDIF

         CALL iom_open( TRIM(sn_isfpar_area%clname), inum )
         CALL iom_get( inum, jpdom_unknown, TRIM(sn_isfpar_area%clvar), zisf_par_area_glo )
         CALL iom_close(inum) 

         !---------------------------------------------------------
         ! 2. Identify which basins are present in this subdomain
         !    !!!! Basin id need to be between 1 and nbasins_glo => to be improve by reading the basin_id variable
         idx_basin_glo_to_loc(:) = 0
         kb_used = 0
         DO jb_glo = 1, nbasins_glo
            IF ( ANY( id_basin_isfpar(:,:) == jb_glo) ) THEN
               kb_used = kb_used + 1
               idx_basin_glo_to_loc(jb_glo) = kb_used
            ENDIF
         END DO
         nbasins_loc = kb_used

         !---------------------------------------------------------
         ! 3. Allocate array that have a local size
         CALL isf_alloc_par_quad('loc')
         !
         ALLOCATE(zarea_exchg(nbasins_loc,jpk), STAT=ialloc)
         IF( ialloc > 0 ) THEN
            CALL ctl_stop( 'isfpar: unable to allocate zarea_exchg' )   ;   RETURN
         ENDIF
         zarea_exchg(:,:) = 0._wp

         !---------------------------------------------------------
         ! 4. compute flag if there is is to parametrized 
         ln_exchg(:) = .FALSE.
         DO jb_glo = 1, nbasins_glo
            jb_loc = idx_basin_glo_to_loc(jb_glo)
            IF ( jb_loc  > 0 ) THEN
               idx_basin_loc_to_glo( jb_loc ) = jb_glo
               ln_exchg(jb_glo) = .TRUE.
            ENDIF
         END DO
         !
         !---------------------------------------------------------
         ! 5. Compute mask, area, jk [nbasin,jpk] for the param
         !
         ! Find interfacial water columns in individual basins, save corresponding mask and area:
         !    mskisf_exchg(ji,jj,kbasin) is a horizontal 2d mask defining the exchange zone for individual basins.
         DO jb_loc = 1, nbasins_loc
            jb_glo = idx_basin_loc_to_glo(jb_loc)
            !
            ! c=MERGE(a,b,l) is a F90 feature eq to IF ( l ) c=a ; ELSE c=b 
            DO_2D( 0, 0, 0, 0 )
               mskisf_exchg(ji,jj,jb_loc) = MERGE(REAL(mskisf_par(ji,jj),wp), 0.0_wp, id_basin_isfpar(ji,jj) == jb_glo)
            END_2D
            !
         END DO
         !
         DO jb_glo = 1, nbasins_glo
            !
            IF ( ln_exchg(jb_glo) ) THEN
               !
               jb_loc = idx_basin_glo_to_loc(jb_glo)
               !
               ! Area of the exchange zone per basin and per vertical level [m^2]
               DO jk = 1,jpk
                  ctmp(jk) = local_2Dsum( e1e2t(A2D(0)) * tmask(A2D(0),jk) * mskisf_exchg(:,:,jb_loc) )
               END DO
               CALL mpp_sum( 'isf_par_init', ctmp(:) )
               zarea_exchg(jb_loc,:) = REAL(ctmp(:),wp)

               IF (SUM(zarea_exchg(jb_loc,:)) == 0._wp) THEN
                  WRITE(cinfo, '(A,I3,A)') 'Basin ', INT(rbasisf_num(jb_glo)), ' is not represented. Check zmin wrt bathymetry'
                  CALL ctl_stop('STOP',cinfo)
               END IF
               !
               ! Compute all the valid wet level (jk_exchg) from area
               DO jk = 1, jpk
                  jk_exchg(jb_loc,jk) = MERGE(jk, 0, zarea_exchg(jb_loc,jk) >= zepsln )
               ENDDO
               !
               ! Fill jk for depth below and above valid data at the front to compute melt in case draft of GL above or below the
               ! valid wet level 
               ! method: nearest neighbor forward/backward fill
               ztmp(:) = REAL(jk_exchg(jb_loc,:), wp)

               ! Collect indices of valid (non-zero) values
               klen = COUNT(jk_exchg(jb_loc,:) /= 0)
               ALLOCATE(klvl(klen), STAT=ialloc)
               IF( ialloc > 0 ) THEN
                  CALL ctl_stop( 'isfpar: unable to allocate klvl' )   ;   RETURN
               ENDIF
               klvl = PACK([(jj, jj = 1, jpk)], MASK = jk_exchg(jb_loc,:) /= 0)

               DO jk = 1, jpk
                  IF (jk_exchg(jb_loc,jk) == 0) THEN

                     ! Find closest valid index to current jk
                     kbest = MINLOC(ABS(klvl - jk), 1)
 
                     ! Fill with nearest valid vertical level
                     ztmp(jk) = jk_exchg(jb_loc,klvl(kbest))
                  END IF
               END DO
               DEALLOCATE(klvl)
               !
               jk_exchg(jb_loc,:) = INT(ztmp(:))
               !
               ! compute 1 / area_exchag as this is the what is used after:
               PRINT *, 'zarea_exchg : ',MAXVAL(zarea_exchg), MINVAL(zarea_exchg)
               r1_area_exchg(:,:) = 1.0_wp / MAX(zarea_exchg(:,:),zepsln)
               !
               ! compute scale factor to convert tf2s into melt
               ! Coefficient (zcoef) for the melting parameterization:
               !   see eq. 14 of Burgard et al. (2022)
               !   NB: their melt is in meters of ice per second while we here use kg/m2/s
               !   NB1: rn_isfpar_Kcoeff is defined in eq. 16 of Brugard et al. (2022)
               !        and was calibrated at 1.16e-4 in that paper (here specified by user).
               zcoef = rn_isfpar_Kcoeff * rho0 * ( rcp / rLfusisf )**2 * zbeta * grav * zsin_theta * 0.5_wp / zf
               rtf2s_to_melt_loc(jb_loc,:) = zcoef * zisf_par_area_glo(jb_glo,:)
               !
            ELSE
               ! need to do the mpp_sum here to match the one done above by subdomain that contain the basin jb_glo
               ctmp(:) = CMPLX( 0.e0, 0.e0, dp )
               CALL mpp_sum( 'isf_par_init', ctmp(:) )  
            END IF
            !
         END DO
         DEALLOCATE(zisf_par_area_glo)
         !
      CASE ( 'oasis' )
         !
         ALLOCATE( sf_isfpar_fwf(1), STAT=ierr )
         ALLOCATE( sf_isfpar_fwf(1)%fnow(A2D(0),1), sf_isfpar_fwf(1)%fdta(A2D(0),1,2) )
         CALL fld_fill( sf_isfpar_fwf, (/ sn_isfpar_fwf /), cn_isfdir, 'isf_par_init', 'read fresh water flux isf data', 'namisf' )
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>    isf melt provided by OASIS (cn_isfmlt_par = oasis)'
         !
      CASE DEFAULT
         CALL ctl_stop( 'sbc_isf_init: wrong value of nn_isf' )
      END SELECT
      !
   END SUBROUTINE isf_par_init

END MODULE isfpar
