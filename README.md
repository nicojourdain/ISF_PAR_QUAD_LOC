# Interactive ice-shelf melt parameterisation for NEMO

Implementation of a qudratic local formulation to interactively calculate melt rates beneath ice shelves that are not represented on the NEMO grid. It can be used to parameterise entire ice shelf cavities or parts of ice shelf cavities (e.g. the ocean in the vicinity of the grounding zone).

Calculations are made per glacial drainage basin to avoid computation of nearest neighbours when the ice shelf grounding line evolves.

The required inputs are zmin,zmax (as previously for prescribed melt and Beckman & Goose 2003), and **isfpar\_basin** the map of basin numbers, and **isfpar\_area** the non-resolved ice shelf area (m^2) per basin and per vertical level.

There are new inputs in the ocean namelist:
```bash
      ln_isfpar_mlt = .true.   ! ice shelf melting parametrised
         cn_isfpar_mlt = 'quad_loc'  ! ice shelf melting parametrisation (spe/bg03/oasis)
         !                           ! spe      = fwfisf is read from a forcing field ( melt > 0; freezing < 0 )
         !                           ! quad_loc = melt computed using a quadratic local parameterisation ( melt > 0; freezing < 0 )
         !                           ! bg03     = melt computed using Beckmann and Goosse parametrisation
         !                           ! oasis    = fwfisf is given by oasis and pattern by file sn_isfpar_fwf
         !
         !* 'quad_loc' case
         nn_isfpar_basin = 155       ! max basin ID in sn_isfpar_basin
         rn_isfpar_Kcoeff = 7.20e-3  ! Tuning coef for the parameterization
         !
         !*** File definition ***
         !
         !* all cases
         !______________!___________________!___________________!___________!_____________!_________!___________!__________!__________!_______________!
         !              !     file name     ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
         !              !                   !  (if <0  months)  !   name    !  (logical)  !  (T/F)  ! 'monthly' ! filename ! pairing  ! filename      !
         sn_isfpar_zmax = 'isfpar_zmin_zmax',        0.         ,  'zmax'   ,  .false.    , .true.  , 'yearly'  ,    ''    ,   ''     ,    ''
         sn_isfpar_zmin = 'isfpar_zmin_zmax',        0.         ,  'zmin'   ,  .false.    , .true.  , 'yearly'  ,    ''    ,   ''     ,    ''
         !
         !* 'quad_loc' case
         sn_isfpar_basin = 'isfpar_basin'   ,        0.         ,  'basin'  ,  .false.    , .true.  , 'yearly'  ,    ''    ,   ''     ,    ''
         sn_isfpar_area  = 'isfpar_area'    ,        0.         ,  'area'   ,  .false.    , .true.  , 'yearly'  ,    ''    ,   ''     ,    ''
         !
```

The melt rates are saved in usual files as a (x,y,time) field. In addition, a new grid is defined in grid\_def\_nemo.xml:
```bash
  <grid id="isfbas" >
     <axis axis_ref="nisfbas" />
```
so that the melt is also saved in the (nbasin,z) dimensions. This can be used by a coupled ice sheet model to directly redistribute melt rates per basin and per vertical level.


