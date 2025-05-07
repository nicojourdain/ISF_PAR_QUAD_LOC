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
   USE isfrst   , ONLY: isfrst_write, isfrst_read ! ice shelf restart read/write subroutine
   USE isftbl   , ONLY: isf_tbl_ktop, isf_tbl_lvl ! ice shelf top boundary layer properties subroutine
   USE isfparmlt, ONLY: sf_isfpar_fwf, isfpar_mlt ! ice shelf melt formulation subroutine
   USE isfdiags , ONLY: isf_diags_flx             ! ice shelf diags subroutine
   USE isfutils , ONLY: debug, read_2dcstdta      ! ice shelf debug subroutine
   !
   USE dom_oce  , ONLY: bathy          ! ocean space and time domain
   USE par_oce  , ONLY: jpi,jpj        ! ocean space and time domain
   USE phycst   , ONLY: r1_rho0_rcp    ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE fldread        ! read input field at current time step
   USE lib_fortran, ONLY: local_sum
   USE lib_mpp    , ONLY: mpp_sum, ctl_stop !

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
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcisf.F90 10536 2019-01-16 19:21:09Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
 
   SUBROUTINE isf_par( kt, Kmm, ptsc, pqfwf )
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
      !!-------------------------- OUT --------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)     , INTENT(inout) :: pqfwf
      REAL(wp), DIMENSION(jpi,jpj,jpts), INTENT(inout) :: ptsc
      !!-------------------------- IN  --------------------------------------
      INTEGER, INTENT(in) ::   kt                                     ! ocean time step
      INTEGER, INTENT(in) ::   Kmm                                    ! ocean time level index
      !!---------------------------------------------------------------------
      INTEGER ::   ji, jj
      REAL(wp), DIMENSION(jpi,jpj) :: zqoce, zqhc, zqlat, zqh
      !!---------------------------------------------------------------------
      !
      ! compute heat content, latent heat and melt fluxes (2d)
      CALL isfpar_mlt( kt, Kmm, zqhc, zqoce, pqfwf  )
      !
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         ! compute heat and water flux (from isf to oce)
         pqfwf(ji,jj) = pqfwf(ji,jj) * mskisf_par(ji,jj)
         zqoce(ji,jj) = zqoce(ji,jj) * mskisf_par(ji,jj)
         zqhc (ji,jj) = zqhc(ji,jj)  * mskisf_par(ji,jj)
         !
         ! compute latent heat flux (from isf to oce)
         zqlat(ji,jj) = - pqfwf(ji,jj) * rLfusisf    ! 2d latent heat flux (W/m2)
         !
         ! total heat flux (from isf to oce)
         zqh(ji,jj) = ( zqhc (ji,jj) + zqoce(ji,jj) )
         !
         ! set temperature content
         ptsc(ji,jj,jp_tem) = zqh(ji,jj) * r1_rho0_rcp
      END_2D
      !
      ! output fluxes
      CALL isf_diags_flx( Kmm, misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par, 'par', pqfwf, zqoce, zqlat, zqhc)
      !
      ! write restart variables (qoceisf, qhcisf, fwfisf for now and before)
      IF (lrst_oce) CALL isfrst_write(kt, 'par', ptsc, pqfwf)
      !
      IF ( ln_isfdebug ) THEN
         IF(lwp) WRITE(numout,*)
         CALL debug('isf_par: ptsc T',ptsc(:,:,1))
         CALL debug('isf_par: ptsc S',ptsc(:,:,2))
         CALL debug('isf_par: pqfwf fwf',pqfwf(:,:))
         IF(lwp) WRITE(numout,*)
      END IF
      !
   END SUBROUTINE isf_par

   SUBROUTINE isf_par_init
      !!------------------------------------------------------------------------------
      !!                  ***  ROUTINE isf_par_init  ***
      !!
      !! ** Purpose : initialisation of the variable needed for the parametrisation of ice shelf melt
      !!
      !!------------------------------------------------------------------------------
      INTEGER                          :: ierr, inum, jb_loc, jb_glo, ji, jj, jk, indx
      REAL(wp)                         :: epsln = 1.e-20_wp  
      INTEGER                          :: kbest, kb_used, klen, ialloc
      REAL(wp), DIMENSION(jpk)         :: ztmp
      COMPLEX(dp), DIMENSION(jpk)      :: ctmp
      INTEGER , ALLOCATABLE, DIMENSION(:)   :: klvl
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zisf_par_area_glo  ! area distribution per isf basin (from input file)
      REAL(wp), DIMENSION(jpi,jpj)     :: ztblmax, ztblmin, zzid
      !
      !!------------------------------------------------------------------------------
      !
      ! allocation
      CALL isf_alloc_par()
      !
      ! initialisation
      misfkt_par(:,:)     = 1         ; misfkb_par(:,:)       = 1         
      rhisf_tbl_par(:,:)  = 1e-20     ; rfrac_tbl_par(:,:)    = 0.0_wp
      !
      ! define isf tbl tickness, top and bottom indice
      CALL read_2dcstdta(TRIM(sn_isfpar_zmax%clname), TRIM(sn_isfpar_zmax%clvar), ztblmax)
      CALL read_2dcstdta(TRIM(sn_isfpar_zmin%clname), TRIM(sn_isfpar_zmin%clvar), ztblmin)
      !
      ! mask ice shelf parametrisation location
      ztblmax(:,:) = ztblmax(:,:) * ssmask(:,:)
      ztblmin(:,:) = ztblmin(:,:) * ssmask(:,:)
      !
      ! if param used under an ice shelf overwrite ztblmin by the ice shelf draft
      WHERE ( risfdep > 0._wp .AND. ztblmin > 0._wp )
         ztblmin(:,:) = risfdep(:,:)
      END WHERE
      !
      ! ensure ztblmax <= bathy
      WHERE ( ztblmax(:,:) > bathy(:,:) )
         ztblmax(:,:) = bathy(:,:)
      END WHERE
      !
      ! compute index of top ocean level and update ztblmin to gdepw_0(misfkt_par) 
      CALL isf_tbl_ktop(ztblmin, misfkt_par) !   out: misfkt_par
      !                                      ! inout: ztblmin
      !
      ! initial tbl thickness
      rhisf0_tbl_par(:,:) = ztblmax(:,:) - ztblmin(:,:)
      !
      ! define iceshelf parametrisation mask
      mskisf_par = 0
      WHERE ( rhisf0_tbl_par(:,:) > 0._wp )
         mskisf_par(:,:) = 1._wp
      END WHERE
      !
      ! read par variable from restart
      IF ( ln_rstart ) CALL isfrst_read('par', risf_par_tsc, fwfisf_par, risf_par_tsc_b, fwfisf_par_b)
      !
      SELECT CASE ( TRIM(cn_isfpar_mlt) )
         !
      CASE ( 'spe' )
         !
         ALLOCATE( sf_isfpar_fwf(1), STAT=ierr )
         ALLOCATE( sf_isfpar_fwf(1)%fnow(jpi,jpj,1), sf_isfpar_fwf(1)%fdta(jpi,jpj,1,2) )
         CALL fld_fill( sf_isfpar_fwf, (/ sn_isfpar_fwf /), cn_isfdir, 'isf_par_init', 'read fresh water flux isf data', 'namisf' )
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>   ice shelf melt rate read from forcing field (cn_isfmlt_par = spe)'
         !
      CASE ( 'bg03' )
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>   ice shelf melt rate calculated through the' 
         IF(lwp) WRITE(numout,*) '              Beckman and Goose 2003 parametrisation (cn_isfmlt_par = bg03)'
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
         ! get basin map
         CALL iom_get  ( inum, jpdom_global , TRIM(sn_isfpar_basin%clvar), zzid)
         id_basin_isfpar(:,:) = NINT(zzid(:,:))
         ! get basin array
         CALL iom_get  ( inum, jpdom_unknown, 'basin', rbasisf_num(1:nbasins_glo))
         CALL iom_close( inum )
         rbasisf_num(:) = 0._wp
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
            IF ( ANY( id_basin_isfpar(:,:) * tmask_i(:,:) == jb_glo) ) THEN
               kb_used = kb_used + 1
               idx_basin_glo_to_loc(jb_glo) = kb_used
            ENDIF
         END DO
         nbasins_loc = kb_used

         !---------------------------------------------------------
         ! 3. Allocate array that have a local size
         CALL isf_alloc_par_quad('loc')

         ln_exchg(:) = .FALSE.
         DO jb_glo = 1, nbasins_glo
            jb_loc = idx_basin_glo_to_loc(jb_glo)
            IF ( jb_loc  > 0 ) THEN
               idx_basin_loc_to_glo( jb_loc ) = jb_glo
               ln_exchg(jb_glo) = .TRUE.
            ENDIF
         END DO
         !
         ! Find interfacial water columns in individual basins, save corresponding mask and area:
         !    mskisf_exchg(ji,jj,kbasin) is a horizontal 2d mask defining the exchange zone for individual basins.
         !mskisf_exchg(:,:,:) = 0._wp
         DO jb_loc = 1, nbasins_loc
            jb_glo = idx_basin_loc_to_glo(jb_loc)
            !
            ! c=MERGE(a,b,l) is a F90 feature eq to IF ( l ) c=a ; ELSE c=b 
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               mskisf_exchg(ji,jj,jb_loc) = MERGE(1.0_wp, 0.0_wp, id_basin_isfpar(ji,jj) == jb_glo)
            END_2D
            !
         END DO
         !
         ! likely not needed as all the cell filled later.
         !jk_exchg(:,:) = 0
         DO jb_glo = 1, nbasins_glo
            !
            IF ( ln_exchg(jb_glo) ) THEN
               !
               jb_loc = idx_basin_glo_to_loc(jb_glo)
               !
               ! Area of the exchange zone per basin and per vertical level [m^2]
               DO jk = 1,jpk
                  ctmp(jk) = local_sum( e1e2t(:,:) * tmask(:,:,jk) * mskisf_exchg(:,:,jb_loc) )
               END DO
               CALL mpp_sum( 'isf_par_init', ctmp(:) )
               area_exchg(jb_loc,:) = REAL(ctmp(:),wp)
               !
               ! Compute all the valid wet level (jk_exchg) from area
               DO jk = 1, jpk
                  jk_exchg(jb_loc,jk) = MERGE(jk, 0, area_exchg(jb_loc,jk) >= epsln )
               ENDDO
               !
               ! Extract ice shelf draft area for only the local basin
               risf_par_area_loc(jb_loc,:) = zisf_par_area_glo(jb_glo,:)
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
         ALLOCATE( sf_isfpar_fwf(1)%fnow(jpi,jpj,1), sf_isfpar_fwf(1)%fdta(jpi,jpj,1,2) )
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
