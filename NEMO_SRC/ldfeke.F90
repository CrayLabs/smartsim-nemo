MODULE ldfeke
   !!======================================================================
   !!                       ***  MODULE  ldfeke  ***
   !! Ocean physics:  Eddy induced velocity coefficient according to the
   !!                 GEOMETRIC parameterization scheme
   !!=====================================================================
   !! History :  4.0  !  2017-11  (J. Mak, G. Madec) original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_eke       : time step depth-integrated EKE and update aeiw (and
   !!                   from that aeiu and aeiv) according to the GEOMETRIC
   !!                   parameterization scheme
   !!   ldf_eke_init  : initialization, namelist read, and parameters control
   !!   ldf_eke_rst   : read/write eke restart in ocean restart file
   !!----------------------------------------------------------------------
   USE oce            ! ocean: dynamics and active tracers variables
   USE phycst         ! physical constants
   USE dom_oce        ! domain: ocean
   USE ssdb           ! SmartSim database client
   USE ldfslp         ! lateral physics: slope of iso-neutral surfaces
   USE ldftra         ! lateral physics: eddy coefficients
   USE zdfmxl         ! vertical physics: mixed layer

   USE geo2ocean,  ONLY: rot_rep         ! rotation from east-north to i-j

   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_eke        ! routine called in stprk3.F90
   PUBLIC   ldf_eke_eiv    ! routine called in stprk3.F90
   PUBLIC   ldf_eke_init   ! routine called in nemogcm.F90


   !                                 !!** Namelist  namldf_eke  **
   REAL(wp) ::   rn_ekedis                  !  dissipation time scale of EKE                  [days]
   REAL(wp) ::   rn_geom                    !  geometric parameterization master coefficient     [-]
   REAL(wp) ::   rn_eke_lap                 !  diffusion of EKE                               [m2/s]
   REAL(wp) ::   rn_eke_init                !  initial value of total EKE                    [m2/s2]
   REAL(wp) ::   rn_eke_min                 !  background value of total EKE                 [m2/s2]
   REAL(wp) ::   rn_ross_min                !  tapering based of minimum Rossby radius           [m]
   REAL(wp) ::   rn_aeiv_min, rn_aeiv_max   !  min and max bounds of aeiv coefficient         [m2/s]
   REAL(wp) ::   rn_SFmin, rn_SFmax         !  min and max bounds of Structure Function          [-]
   REAL(wp) ::   zf0, zbeta                 !  f0 and beta for computing Rossby speed
   !
   INTEGER, PUBLIC  ::   nn_eke_opt         !  option for what term to include in eddy energy budget
   INTEGER          ::   nn_eke_dis         !  option for taking constant or spatially varying linear dissipation
   !
   LOGICAL  ::   ln_adv_wav                 !  option for having advection by Rossby wave or not
   LOGICAL  ::   ln_beta_plane              !  option for computing long Rossby phase speed
   !
   CHARACTER(len=32) :: cn_model_key = "ml-eke" !  key where ML model for EKE is stored


   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   eke_geom              ! vertical sum of total Eddy Kinetic Energy [m3/s2]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   r1_ekedis             ! linear dissipation rate (= 1/rn_ekedis)   [  /s ]

   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   eke_keS       ! source term of eke due to KE dissipation

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_eke_eiv( kt, Kbb )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ldf_eke_eiv  ***
      !!
      !! ** Purpose :   compute the eddy induced velocity coefficient (aeiu/v)
      !!              from the eddy kinetic energy (eke_geom)
      !!
      !! ** Notes   :   If nn_aei_ijk_t = 32 then aeiuw/v is updated
      !!
      !! ** Method  :   GEOMETRIC calculates the Gent-McWilliams / eddy induced
      !!              velocity coefficient according to
      !!
      !!                    aeiw = alpha * (\int E dz) / (\int S M^2 / N dz),
      !!
      !!              where (\int E dz) is the depth-integrated eddy energy
      !!              (at the previous time level), informed by a parameterized
      !!              depth-integrated eddy energy, where
      !!
      !! ** Action  : * calculate aeiw
      !!
      !! References : Marshall, Maddison, Berloff JPO 2012   ; Mak et al JPO 2018
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      INTEGER, INTENT(in) ::   Kbb      ! time-level indicators
      !
      INTEGER  ::   ji, jj, jk          ! dummy loop arguments
      INTEGER  ::   ik                  ! local integer
      REAL(wp) ::   z1_f20, zfw, ze3w   ! local scalar
      REAL(wp) ::   zn2                 !   -      -
      REAL(wp) ::   zck, zwslpi, zwslpj !   -      -
      REAL(wp) ::   zn2_min = 1.e-8_wp  ! minimum value of N2 used to compute structure function
      REAL(wp), DIMENSION(A2D(1),jpk) ::   zwrk3d      ! 3D workspace (structure function and then 3D EIV coeff)
      REAL(wp), DIMENSION(A2D(1))     ::   zn_slp      ! 2D workspace, PE-KE conversion
      REAL(wp), DIMENSION(A2D(1))     ::   zn, zross   !  -     -      tapering near coasts
      REAL(wp), DIMENSION(A2D(1))     ::   zaeiw       ! 2D EIV coefficient
      !!--------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('ldf_eke_eiv')
      !
      IF( kt == nit000 .AND. lwp) THEN       !* Control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_eke_eiv : GEOMETRIC parameterization: eddy induced velocity coefficient'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !                         work out the 3D structure function (SF) here and
      !                         store in 3d workspace
      !
      ! current: Ferreria et al, SF = N^2 / N^2_ref, W-points
      !                          capped between rn_SFmax and rn_SFmin
      zwrk3d(:,:,:) = 0._wp
      DO_3D( 1, 1, 1, 1, 2, jpkm1 )
         IF( jk <= nmln(ji,jj) ) THEN ! above and at mixed layer
            zwrk3d(ji,jj,jk) = rn_SFmax
         ELSE
            ik   = MIN( nmln(ji,jj), mbkt(ji,jj) ) + 1     ! one level below the mixed layer (MIN in case ML depth is the ocean depth)
            zwrk3d(ji,jj,jk) = MAX( 0._wp , rn2b(ji,jj,jk) ) / MAX( zn2_min , rn2b(ji,jj,ik) )      ! Structure Function : N^2 / N^2_ref
            zwrk3d(ji,jj,jk) = MAX(  rn_SFmin , MIN( zwrk3d(ji,jj,jk) , rn_SFmax )  )
         ENDIF
      END_3D
      !                                                    ! vertical integration
      zn    (:,:) = 0._wp
      zn_slp(:,:) = 0._wp
      DO_3D( 1, 1, 1, 1, 2, jpkm1 )
         zn2     = MAX( 0._wp , rn2b(ji,jj,jk) )
         !
         ze3w      = e3w(ji,jj,jk,Kbb) * tmask(ji,jj,jk)
         !
         zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * ze3w                         ! for working out taper at small rossby radius regions later
         !
         zck     =   ( umask(ji,jj,jk) + umask(ji-1,jj,jk) )   &            ! taken from ldfslp, undo the slope reduction
            &      * ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) * 0.25_wp      !   near topographic features
         zwslpi  = wslpi(ji,jj,jk) / MAX( zck, 0.1_wp)                      ! (just to avoid dividing by zeros)
         zwslpj  = wslpj(ji,jj,jk) / MAX( zck, 0.1_wp)
         !
         zn_slp(ji,jj) = zn_slp(ji,jj) + ze3w * zwrk3d(ji,jj,jk)   &        ! note this >=0 and structure function weighted
            &                          * SQRT( zn2 * ( zwslpi * zwslpi + zwslpj * zwslpj )  )
      END_3D
      !
      IF( l_eke_eiv ) THEN
         !                    !==  resulting EIV coefficient  ==!
         !
         !                          ! surface value
         z1_f20 = 1._wp / (  2._wp * omega * sin( rad * 20._wp )  )
         DO_2D( 1, 1, 1, 1 )
            ! Rossby radius at w-point, Ro = .5 * sum_jpk(N) / f
            zfw = MAX( ABS( ff_t(ji,jj) ) , 1.e-10_wp )
            zross(ji,jj) = 0.5_wp * zn(ji,jj) / zfw
            !
            zaeiw(ji,jj) = rn_geom * eke_geom(ji,jj) / MAX( 1.e-10_wp , zn_slp(ji,jj) ) * tmask(ji,jj,1)   ! zn_slp has SF multiplied to it
            zaeiw(ji,jj) = MIN(  rn_aeiv_max, zaeiw(ji,jj)  )                                      ! bound aeiv from above
            zaeiw(ji,jj) = zaeiw(ji,jj)      &
               ! tanh taper to deal with some some large values near coast
               &           * 0.5_wp * (   1._wp - TANH(  ( -ht_0(ji,jj) + 800._wp     ) / 300._wp  )   )   &
               ! tanh taper of aeiv on internal Rossby radius
               &           * 0.5_wp * (   1._wp + TANH(  ( zross(ji,jj) - rn_ross_min ) * 0.5_wp   )   )
            zaeiw(ji,jj) = zaeiw(ji,jj) * MIN(  1._wp, ABS( ff_t(ji,jj) * z1_f20 )  )              ! tropical decrease
            zaeiw(ji,jj) = MAX(  rn_aeiv_min, zaeiw(ji,jj)  )                                      ! bound aeiv from below
         END_2D
         !
         !                          ! inner value
         ! bottom value is already set to zero, use the un-masked zaeiw(ji,jj) with the structure function to
         ! set interior values. Re-purpose 3D workspace, replacing structure function (SF) with 3D eiv coefficient
         zwrk3d(:,:,1) = zaeiw(:,:)
         DO_3D( 1, 1, 1, 1, 2, jpkm1 )
            zwrk3d(ji,jj,jk) = zaeiw(ji,jj) * zwrk3d(ji,jj,jk) * wmask(ji,jj,jk)   ! SF has already been capped
         END_3D
         !                          ! aei at u- and v-points
         aeiu(:,:,jpk) = 0._wp
         aeiv(:,:,jpk) = 0._wp
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            aeiu(ji,jj,jk) =    (  ( zwrk3d(ji,jj,jk  ) + zwrk3d(ji+1,jj,jk  ) )    &   ! add () for NP repro
               &                 + ( zwrk3d(ji,jj,jk+1) + zwrk3d(ji+1,jj,jk+1) )  ) &   ! add () for NP repro
               &           / MAX(  wmask(ji,jj,jk  ) + wmask(ji+1,jj,jk  )        &
               &                 + wmask(ji,jj,jk+1) + wmask(ji+1,jj,jk+1) , 1._wp ) * umask(ji,jj,jk)
            aeiv(ji,jj,jk) =    (  ( zwrk3d(ji,jj,jk  ) + zwrk3d(ji,jj+1,jk  ) )    &   ! add () for NP repro
               &                 + ( zwrk3d(ji,jj,jk+1) + zwrk3d(ji,jj+1,jk+1) )  ) &   ! add () for NP repro
               &           / MAX(  wmask(ji,jj,jk  ) + wmask(ji,jj+1,jk  )        &
               &                 + wmask(ji,jj,jk+1) + wmask(ji,jj+1,jk+1) , 1._wp ) * vmask(ji,jj,jk)
         END_3D
         !
         CALL lbc_lnk( 'ldfeke', aeiu , 'U', 1._wp , aeiv , 'V', 1._wp )
      END IF
      !
      !                    !==  output EKE related variables  ==!
      !
      !                    !--  diagnostics  --!
      CALL iom_put( "aeiv_geom" , zwrk3d(A2D(0),jpk) )      ! eddy induced coefficient from GEOMETRIC param
      CALL iom_put( "rossby_rad", zross(A2D(0))      )      ! internal Rossby deformation radius
!!
!! jm: modified ldf_tra in ldftra.F90 to include the ln_eke_equ flag
!!     into consideration, otherwise iom_put is called twice aeiu and aeiv
!!
      CALL iom_put( "aeiu_2d", aeiu(:,:,1) )      ! surface u-EIV coeff.
      CALL iom_put( "aeiv_2d", aeiv(:,:,1) )      ! surface v-EIV coeff.
      CALL iom_put( "aeiu_3d", aeiu(:,:,:) )      ! 3D      u-EIV coeff.
      CALL iom_put( "aeiv_3d", aeiv(:,:,:) )      ! 3D      v-EIV coeff.
      !
      IF( ln_timing )  CALL timing_stop('ldf_eke_eiv')
      !
   END SUBROUTINE ldf_eke_eiv


   SUBROUTINE ldf_eke( kt, Kbb )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ldf_eke  ***
      !!
      !! ** Purpose :   compute eddy kinetic energy wrt GEOMETRIC parametrization
      !!
      !! ** Notes   :   If nn_aei_ijk_t = 32 then eke and aeiw are BOTH updated
      !!                If ln_eke_equ = .true. in namtra_ldfeiv but nn_aei_ijk_t
      !!              is something else, then ONLY eddy equation is updated
      !!              (but the eddy energy is passive and doesn't do anything)
      !!
      !! ** Method  :   GEOMETRIC calculates the Gent-McWilliams / eddy induced
      !!              velocity coefficient according to
      !!
      !!                    aeiw = alpha * (\int E dz) / (\int S M^2 / N dz),
      !!
      !!              where (\int E dz) is the depth-integrated eddy energy
      !!              (at the previous time level), informed by a parameterized
      !!              depth-integrated eddy energy, where
      !!
      !!    nn_eke_opt    =  0 => default: just PE->EKE growth and linear dissipation
      !!                  !
      !!                  =  1 => default + advection by depth-averaged flow
      !!                  !
      !!                  =  2 => default + advection + contribution to EKE from
      !!                          momentum dissipation
      !!                  !
      !!                  = 88 => ONLY advection
      !!                  !
      !!                  = 99 => ONLY Laplacian diffusion
      !!
      !!              S is a structure function, and M and N are horizontal and
      !!              vertical buoyancy frequencies
      !!
      !!                  linear dissipation may be specified by
      !!
      !!    nn_eke_dus    =  0 => constant
      !!                  !
      !!                  =-20 => read in a geom_diss_2D.nc field (give it in days)
      !!
      !! ** Action  : * time step depth-integrated eddy energy budget
      !!
      !! References : Marshall, Maddison, Berloff JPO 2012   ; Mak et al JPO 2018
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      INTEGER, INTENT(in) ::   Kbb      ! time-level indicators
      !
      INTEGER  ::   ji, jj, jk          ! dummy loop arguments
      INTEGER  ::   rc                  ! return code
      REAL(wp) ::   z1_f20, ze3w        ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj      !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj      !   -      -
      REAL(wp) ::   zn_slp2, zn2        !   -      -
      REAL(wp) ::   zmsku, zaeiu_w      !   -      -
      REAL(wp) ::   zmskv, zaeiv_w      !   -      -
      REAL(wp) ::   zen, zed            !   -      -
      REAL(wp) ::   zeke_rhs            !   -      -
      REAL(wp) ::   zc_rosu, zc_rosv    !   -      -
      REAL(wp) ::   ze3wd, zck          !   -      -
      REAL(wp) ::   zfw, ztmpx, ztmpy   !   -      -
      REAL(wp) ::   zwslpi, zwslpj      !   -      -
      REAL(wp), DIMENSION(A2D(0))     ::   zeke_peS               ! 2D workspace, PE-KE conversion
      REAL(wp), DIMENSION(A2D(0))     ::   zadv_ubt               !  -     -      barotropic advection
      REAL(wp), DIMENSION(A2D(0))     ::   zlap                   !  -     -      diffusion
      REAL(wp), DIMENSION(A2D(0))     ::   zdis                   !  -     -      linear dissipation
      REAL(wp), DIMENSION(A2D(0))     ::   zadv_wav               !  -     -      wave advection
      REAL(wp), DIMENSION(A2D(1))     ::   zn                     !  -     -      tapering near coasts
      REAL(wp), DIMENSION(A2D(1))     ::   zwx, zwy, zeke_ht      !  -     -      barotropic advection
      REAL(wp), DIMENSION(A2D(1))     ::   zc_rosi, zc_rosj       !  -     -      wave advection
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zc1, zc_ros      !  1st baroclinic mode and long Rossby speed [m /s ]
      REAL(wp), DIMENSION(A2D(1)) ::   zwrk2d, zn_slp      ! 2D workspace (to calculate vertically averaged isopycnal slope)
      REAL(wp), DIMENSION(SIZE(zwrk2d),4) :: z_feat_vec
      REAL(wp), DIMENSION(SIZE(zwrk2d)) :: z_pred_vec
      REAL(wp), DIMENSION(A2D(0)) :: z_rv, z_ross, zn_mke
      CHARACTER(len=25), DIMENSION(1) :: zc_feat_key, zc_pred_key
      !!--------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('ldf_eke')
      !
      IF( kt == nit000 .AND. lwp) THEN       !* Control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_eke : GEOMETRIC parameterization (total EKE time evolution equation)'
         WRITE(numout,*) '~~~~~~~'
      ENDIF
      !                    !==  EIV mean conversion to EKE (ah N^2 slp2) & N slp  ==!
      !
      !                         !*  parametrized PE_EKE conversion due to eddy induced velocity
      zeke_peS(:,:) = 0._wp
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         zmsku = wmask(ji,jj,jk) / MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)          &
            &                           + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk) , 1._wp  )
         zmskv = wmask(ji,jj,jk) / MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)          &
            &                           + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk) , 1._wp  )
         !
         zaeiu_w = (   ( aeiu(ji  ,jj,jk-1) + aeiu(ji-1,jj,jk) )          &   !  add () for NP repro
            &        + ( aeiu(ji-1,jj,jk-1) + aeiu(ji  ,jj,jk) )  ) * zmsku   !  add () for NP repro
         zaeiv_w = (   ( aeiv(ji,jj  ,jk-1) + aeiv(ji,jj-1,jk) )          &   !  add () for NP repro
            &        + ( aeiv(ji,jj-1,jk-1) + aeiv(ji,jj  ,jk) )  ) * zmskv   !  add () for NP repro
         !
         zn_slp2 = (   ( zaeiu_w * wslpi(ji,jj,jk) * wslpi(ji,jj,jk) )   &    ! (slope ** 2) * aeiv + add () for NP repro
            &        + ( zaeiv_w * wslpj(ji,jj,jk) * wslpj(ji,jj,jk) )   )    ! JM 28 Jun: undo slope reduction here too?
         zn2     = MAX( 0._wp , rn2b(ji,jj,jk) )
         !
         ze3w      = e3w(ji,jj,jk,Kbb) * tmask(ji,jj,jk)
         !
         zeke_peS(ji,jj) = zeke_peS(ji,jj) + ze3w * zn2 * zn_slp2           ! note this is >=0
      END_3D
      !
      ! temporary depth-averaged eke variable for use in advection and diffusion
      DO_2D( 1, 1, 1, 1 )
         zeke_ht(ji,jj) = eke_geom(ji,jj) / MAX( ht(ji,jj,Kbb), 1._wp ) * tmask(ji,jj,1)
      END_2D
      !
      !                         !*  upstream advection with initial mass fluxes & intermediate update
      !                               !* upstream tracer flux in the i and j direction
      DO_2D( 1, 0, 1, 0 )
         ! upstream scheme by depth-averaged velocity (which should be un_adv * r1_hu etc.)
         ! however, zwx with depth-averaged eke should then multiplied by hu to restore dimensions, so no r1_h[uv] factor here
         !
         zfp_ui = un_adv(ji,jj) + ABS( un_adv(ji,jj) )
         zfm_ui = un_adv(ji,jj) - ABS( un_adv(ji,jj) )
         zfp_vj = vn_adv(ji,jj) + ABS( vn_adv(ji,jj) )
         zfm_vj = vn_adv(ji,jj) - ABS( vn_adv(ji,jj) )
         ! advection but with dimensions restored (would have a h[uv] factor normally, but cancelled from above)
         !
         zwx(ji,jj) = 0.5_wp * e2u(ji,jj) * ( zfp_ui * zeke_ht(ji  ,jj  )   &
                                            + zfm_ui * zeke_ht(ji+1,jj  ) )
         zwy(ji,jj) = 0.5_wp * e1v(ji,jj) * ( zfp_vj * zeke_ht(ji  ,jj  )   &
                                            + zfm_vj * zeke_ht(ji  ,jj+1) )
      END_2D
      !                               !* divergence of ubt advective fluxes
      DO_2D( 0, 0, 0, 0 )
         zadv_ubt(ji,jj) = - (  ( zwx(ji,jj) - zwx(ji-1,jj  ) )   &
            &                 + ( zwy(ji,jj) - zwy(ji  ,jj-1) )  ) * r1_e1e2t(ji,jj)
      END_2D
      !
      IF( lk_linssh ) THEN            !* linear ssh : add advective fluxes through z=0
         DO_2D( 0, 0, 0, 0 )
            zadv_ubt(ji,jj) = - ww(ji,jj,1) * eke_geom(ji,jj) / ( ht_0(ji,jj) + 1._wp-tmask(ji,jj,1) )   ! jm (03 Mar 18): was eke_n
         END_2D
      ENDIF
      !
      !                         !* same as above but for advection by Rossby waves
      zadv_wav(:,:) = 0._wp
      IF( ln_adv_wav ) THEN
         ALLOCATE( zc1(A2D(1)) , zc_ros(A2D(1)) )
         !
         zn(:,:) = 0._wp                        ! for working out taper at small rossby radius regions
         DO_3D( 1, 1, 1, 1, 2, jpkm1 )
            zn2       = MAX( 0._wp , rn2b(ji,jj,jk) )
            ze3w      = e3w(ji,jj,jk,Kbb) * tmask(ji,jj,jk)
            zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * ze3w
         END_3D
         !
         zc1   (:,:) = 0._wp
         zc_ros(:,:) = 0._wp
         DO_2D( 1, 1, 1, 1 )
            ! compute only for deep enough places
            IF( ht_0(ji,jj) > 300._wp ) THEN   ! jm: should use something like ht_b really...
               ! compute vertical mode phase speed on T point
               zc1(ji,jj) = MIN( 10._wp, zn(ji,jj) / rpi )
               ! compute long Rossby phase speed on T point (minus sign later)
               IF( ln_beta_plane ) THEN
                  zc_ros(ji,jj) = zc1(ji,jj) * zc1(ji,jj) * ABS(zbeta) / (zf0 * zf0)
               ELSE
                  zc_ros(ji,jj) = zc1(ji,jj) * zc1(ji,jj) * COS( rad * gphit(ji,jj) )   &
                                / (  ra * ff_t(ji,jj) * SIN( rad * gphit(ji,jj) )        &
                                + rsmall  )
               ENDIF
               ! cap the Rossby phase speeds by the fastest equatorial Rossby wave speed
               ! the minus sign for westward propagation goes here
               zc_ros(ji,jj) = -MIN( zc1(ji,jj) / 3._wp, zc_ros(ji,jj) )
            ENDIF
         END_2D
         !
         zwx(:,:) = 0._wp ! wipe the advective contributions from above
         !
         ! do the rotation of the advective velocity (using the zeroed arrays)
         CALL rot_rep(zc_ros, zwx, 'T', 'en->i', zc_rosi)
         CALL rot_rep(zc_ros, zwx, 'T', 'en->j', zc_rosj)  ! re-use the zero array
         !
         DO_2D( 1, 0, 1, 0 )
            ! average onto grid
            zc_rosu    = 0.5_wp * ( zc_rosi(ji,jj) + zc_rosi(ji+1,jj  ) ) * umask(ji,jj,1)
            zc_rosv    = 0.5_wp * ( zc_rosj(ji,jj) + zc_rosj(ji  ,jj+1) ) * vmask(ji,jj,1)
            ! upstream advection
            zfp_ui     = zc_rosu + ABS( zc_rosu )
            zfm_ui     = zc_rosu - ABS( zc_rosu )
            zfp_vj     = zc_rosv + ABS( zc_rosv )
            zfm_vj     = zc_rosv - ABS( zc_rosv )
            !
            ! create trend with appropriate averaging (will be divided e1e2t later)
            zwx(ji,jj) = 0.5_wp * e2u(ji,jj) * hu(ji,jj,Kbb) * ( zfp_ui * zeke_ht(ji  ,jj  )   &
                                                               + zfm_ui * zeke_ht(ji+1,jj  ) )
            zwy(ji,jj) = 0.5_wp * e1v(ji,jj) * hv(ji,jj,Kbb) * ( zfp_vj * zeke_ht(ji  ,jj  )   &
                                                               + zfm_vj * zeke_ht(ji  ,jj+1) )
         END_2D
         !
         !                            !* divergence of wav advective fluxes
         z1_f20 = 1._wp / (  2._wp * omega * SIN( rad * 20._wp )  )

         DO_2D( 0, 0, 0, 0 )
            zadv_wav(ji,jj) = - (  ( zwx(ji,jj) - zwx(ji-1,jj  ) )   &
               &                 + ( zwy(ji,jj) - zwy(ji  ,jj-1) )  ) * r1_e1e2t(ji,jj)
            zadv_wav(ji,jj) = zadv_wav(ji,jj) * MIN(  1._wp, ABS( ff_t(ji,jj) * z1_f20 )  )   ! tropical decrease
         END_2D
      ENDIF
      !
      !                         !* divergence of diffusive fluxes
      !
      ! use depth-averaged eke here; zeke_ht already exists from above
      !
      IF( rn_eke_lap >= 0._wp ) THEN
         DO_2D( 1, 0, 1, 0 )
            zwx(ji,jj) = rn_eke_lap * e2_e1u(ji,jj) * hu(ji,jj,Kbb) * ( zeke_ht(ji+1,jj  ) -  zeke_ht(ji,jj) ) * umask(ji,jj,1)   ! rn_eke_lap is constant (for now) and NOT masked
            zwy(ji,jj) = rn_eke_lap * e1_e2v(ji,jj) * hv(ji,jj,Kbb) * ( zeke_ht(ji  ,jj+1) -  zeke_ht(ji,jj) ) * vmask(ji,jj,1)   !      before it is pahu and pahv which IS masked
         END_2D
         DO_2D( 0, 0, 0, 0 )
            zlap(ji,jj) = (  ( zwx(ji,jj) - zwx(ji-1,jj  ) )   &
               &           + ( zwy(ji,jj) - zwy(ji  ,jj-1) )  ) * r1_e1e2t(ji,jj)
         END_2D
      ELSE
         zlap(:,:) = 0._wp
      ENDIF
                                !* form the trend for linear dissipation
      DO_2D( 0, 0, 0, 0 )
         zdis(ji,jj) = - r1_ekedis(ji,jj) * (eke_geom(ji,jj) - rn_eke_min) * tmask(ji,jj,1)
      END_2D
      !
      !                    !==  time stepping of EKE Eq.  ==!
      !
      ! note: the rn_eke_min term is a forcing in eddy equation, thus damping in mean equation
      !       added to prevent overshoots and oscillations of energy from exponential growth/decay
      !       maintains a background (depth-integrated) eddy energy level
      !
      !
      SELECT CASE( nn_eke_opt )   ! Specification of what to include in EKE budget
      !
      CASE(   0  )  !  default: just PE->EKE growth and linear dissipation
         DO_2D( 0, 0, 0, 0 )
            zeke_rhs =                                     zeke_peS(ji,jj) + zdis(ji,jj) + zlap(ji,jj)
            eke_geom(ji,jj) = eke_geom(ji,jj) + rDt * zeke_rhs * ssmask(ji,jj)
         END_2D
      CASE(   1  )  !  as default but with full advection
         DO_2D( 0, 0, 0, 0 )
            zeke_rhs = zadv_ubt(ji,jj) + zadv_wav(ji,jj) + zeke_peS(ji,jj) + zdis(ji,jj) + zlap(ji,jj)
            eke_geom(ji,jj) = eke_geom(ji,jj) + rDt * zeke_rhs * ssmask(ji,jj)
         END_2D
      CASE(   2  )  !  full thing with additional KE->EKE growth
         DO_2D( 0, 0, 0, 0 )
            zeke_rhs = zadv_ubt(ji,jj) + zadv_wav(ji,jj) + zeke_peS(ji,jj) + zdis(ji,jj) + zlap(ji,jj) + eke_keS(ji,jj)
            eke_geom(ji,jj) = eke_geom(ji,jj) + rDt * zeke_rhs * ssmask(ji,jj)
         END_2D
      CASE(  77  )  !  ONLY advection by rossby waves
         DO_2D( 0, 0, 0, 0 )
            zeke_rhs = zadv_wav(ji,jj)
            eke_geom(ji,jj) = eke_geom(ji,jj) + rDt * zeke_rhs * ssmask(ji,jj)
         END_2D
      CASE(  88  )  !  ONLY advection by mean flow
         DO_2D( 0, 0, 0, 0 )
            zeke_rhs = zadv_ubt(ji,jj)
            eke_geom(ji,jj) = eke_geom(ji,jj) + rDt * zeke_rhs * ssmask(ji,jj)
         END_2D
      CASE(  99  )  !  ONLY diffusion
         DO_2D( 0, 0, 0, 0 )
            zeke_rhs = zlap(ji,jj)
            eke_geom(ji,jj) = eke_geom(ji,jj) + rDt * zeke_rhs * ssmask(ji,jj)
         END_2D
      CASE(  100  ) ! Calculate via machine-learning model
         !! Calculate features used to predict EKE

         z_rv(:,:) = 0.
         ! Relative vorticity
         DO_2D( 0, 0, 0, 0 )
            z_rv(ji,jj) = (   e2v(ji+1,jj  ) * vv(ji+1,jj  ,1,Kbb) - e2v(ji,jj) * vv(ji,jj,1,Kbb)    &
            &               - e1u(ji  ,jj+1) * uu(ji  ,jj+1,1,Kbb) + e1u(ji,jj) * uu(ji,jj,1,Kbb)  ) * r1_e1e2f(ji,jj)
         END_2D

         ! Vertically averaged isopycnal slope
         zwrk2d(:,:) = 0._wp
         DO_3D( 1, 1, 1, 1, 2, jpkm1 )
            IF( jk > nmln(ji,jj) ) THEN ! above and at mixed layer
               zwrk2d(ji,jj) = zwrk2d(ji,jj) + e3w(ji,jj,jk,Kbb) * tmask(ji,jj,jk)
            ENDIF
         END_3D
         !                                                    ! vertical average
         zn_slp(:,:) = 0._wp
         DO_3D( 1, 1, 1, 1, 2, jpkm1 )
            IF( jk > nmln(ji,jj) ) THEN ! above and at mixed layer
               ze3wd      = e3w(ji,jj,jk,Kbb) * tmask(ji,jj,jk)
               zn2       = MAX( 0._wp , rn2b(ji,jj,jk) )
               zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * ze3w                         ! Needed for Rossby radius calculation
               !
               zck     =   ( umask(ji,jj,jk) + umask(ji-1,jj,jk) )   &            ! taken from ldfslp, undo the slope reduction
                  &      * ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) * 0.25_wp      !   near topographic features
               zwslpi  = wslpi(ji,jj,jk) / MAX( zck, 0.1_wp)                      ! (just to avoid dividing by zeros)
               zwslpj  = wslpj(ji,jj,jk) / MAX( zck, 0.1_wp)
               !
               zn_slp(ji,jj) = zn_slp(ji,jj) + ze3w / zwrk2d(ji,jj)   &           ! Update average weighted by fraction of column height
                  &                          * SQRT( zwslpi * zwslpi + zwslpj * zwslpj )
            ENDIF
         END_3D
         ! Rossby Radius over grid spacing
         DO_2D( 0, 0, 0, 0 )
            ! Rossby radius at w-point, Ro = .5 * sum_jpk(N) / f
            zfw = MAX( ABS( ff_t(ji,jj) ) , 1.e-10_wp )
            z_ross(ji,jj) = 0.5_wp * zn(ji,jj) / zfw
            z_ross(ji,jj) = z_ross(ji,jj) / MAX(e1u(ji,jj), e2v(ji,jj))
         END_2D
         ! Resolved surface kinetic energy
         DO_2D( 0, 0, 0, 0 )
            ztmpx = uu(ji-1,jj  ,1,Kbb) + uu(ji,jj,jk,Kbb)
            ztmpy = vv(ji  ,jj-1,1,Kbb) + vv(ji,jj,jk,Kbb)
            zn_mke(ji,jj) = 0.25_wp * ( ztmpx*ztmpx + ztmpy*ztmpy )
         END_2D

         z_feat_vec(:,1) = PACK(zn_mke, .true.)
         z_feat_vec(:,2) = PACK(zn_slp, .true.)
         z_feat_vec(:,3) = PACK(z_rv, .true.)
         z_feat_vec(:,4) = PACK(z_ross, .true.)

         zc_feat_key(1) = TRIM(client_prefix)//"features"
         zc_pred_key(1) = TRIM(client_prefix)//"prediction"

         rc = sr_client%put_tensor( zc_feat_key(1), z_feat_vec, SHAPE(z_feat_vec) )
         rc = sr_client%run_model( cn_model_key, zc_feat_key, zc_pred_key )
         rc = sr_client%unpack_tensor( zc_pred_key(1), z_pred_vec, SHAPE(z_pred_vec))
         eke_geom(:,:) = RESHAPE(z_pred_vec, SHAPE(eke_geom))

      END SELECT

      !
      CALL lbc_lnk( 'ldfeke', eke_geom(:,:), 'T', 1._wp )   ! Lateral boundary conditions on eke_geom  (unchanged sign)
      !
      IF( lrst_oce .AND. kt == nitrst )   CALL ldf_eke_rst( kt, 'WRITE' )   ! write Kaa since time level
      !                                                                                ! indicators will be shuffled before the
      !                                                                                ! rest of the restart is written
      !                    !==  output EKE related variables  ==!
      !
      CALL iom_put( "eke"            , eke_geom )           ! parameterized total EKE (EPE+ EKE)
      CALL iom_put( "trd_eke_adv_ubt", zadv_ubt )           ! ubt    advective trend of EKE(LHS)
      CALL iom_put( "trd_eke_dis"    , zdis )               ! dissipative trend of EKE     (RHS)
      CALL iom_put( "trd_eke_lap"    , zlap )               ! diffusive trend of EKE       (RHS)
      CALL iom_put( "trd_eke_peS"    , zeke_peS )           ! PE to EKE source trend       (RHS)
      CALL iom_put( "trd_eke_keS"    ,  eke_keS )           ! KE to EKE source trend       (RHS)
      IF( ln_adv_wav ) THEN
         CALL iom_put( "trd_eke_adv_wav", zadv_wav )        ! rossby advective trend of EKE(LHS)
         CALL iom_put( "c1_vert"        , zc1(A2D(0)) )     ! 1st baroclinic mode phase speed
         CALL iom_put( "c_ros"          , zc_ros(A2D(0)) )  ! long Rossby phase speed (should be innreer)
         DEALLOCATE( zc1 , zc_ros )
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('ldf_eke')
      !
   END SUBROUTINE ldf_eke


   SUBROUTINE ldf_eke_init( Kbb )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eke_init  ***
      !!
      !! ** Purpose :   Initialization of the depth-integrated eke, vertical
      !!              structure function and max/min caps of aeiw when using the
      !!              GEOMETRIC parameterization
      !!
      !! ** Method  :   Read the namldf_eke namelist and check the parameters
      !!              called at the first timestep (nit000)
      !!
      !! ** input   :   Namlist namldf_eke
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   Kbb   ! time level indicator
      INTEGER ::   ios, inum
      INTEGER ::   ji, jj
      INTEGER ::   ierr
      !!
      NAMELIST/namldf_eke/  rn_ekedis  , nn_eke_dis ,   & ! GEOM master params (lambda and option, alpha, e0)
         &                  rn_geom    , rn_eke_init,   &
         &                  rn_eke_lap ,                & ! coeff of Laplacian diffusion of EKE
         &                  rn_eke_min ,                & ! background value of (depth-integrated) EKE
         &                  rn_ross_min,                & ! taper aeiv based on Rossby internal radius
         &                  rn_aeiv_max, rn_aeiv_min,   & ! caps and options on result aeiv
         &                  rn_SFmin,    rn_SFmax,      & ! caps on vertical structure function
         &                  nn_eke_opt ,                & ! option for EKE budget terms
         &                  ln_adv_wav ,                & ! option for advection by Rossby wave
         &                  ln_beta_plane               & ! beta plane option for computing Rossby wave speed
         &                  cn_model_key                  ! key where ML model for EKE is stored
      !!----------------------------------------------------------------------
      !
      !
      READ_NML_REF(numnam,namldf_eke)
      READ_NML_(numnam_cfg,cfg,namldf_eke,.TRUE.)
      IF(lwm) WRITE ( numond, namldf_eke )
      !
      IF(lwp) THEN                    !* Control print
         WRITE(numout,*)
         WRITE(numout,*)    'ldf_eke_init : GEOMETRIC parameterization (total EKE time evolution equation)'
         WRITE(numout,*)    '~~~~~~~~~~~~'
         WRITE(numout,*)    '   Namelist namldf_eke : set the GEOMETRIC parameters'
         WRITE(numout,*)    '      aeiw updated according to GEOMETRIC             l_eke_eiv   = ', l_eke_eiv
         WRITE(numout,*)    '      dissipation time scale of EKE in days           rn_ekedis   = ', rn_ekedis, ' days'
         WRITE(numout,*)    '         option for dissipation field                 nn_eke_dis  = ', nn_eke_dis
         SELECT CASE( nn_eke_dis )
         !
         CASE(   0  )
            WRITE(numout,*) '         spatially constant                                         '
         CASE(   2  )
            WRITE(numout,*) '         read in a 2d varying file geom_diss_2D.nc                  '
         END SELECT
         !
         WRITE(numout,*)    '      geometric parameterization master coefficient   rn_geom     = ', rn_geom
         WRITE(numout,*)    '      coeff of Lapican diffusion of EKE               rn_eke_lap  = ', rn_eke_lap
         WRITE(numout,*)    '      initial total EKE value                         rn_eke_init = ', rn_eke_init
         WRITE(numout,*)    '      background values of total EKE                  rn_eke_min  = ', rn_eke_min
         WRITE(numout,*)    '      taper aeiv subject to min Rossby radius         rn_ross_min = ', rn_ross_min
         WRITE(numout,*)    '      maximum bound of aeiv coeff                     rn_aeiv_max = ', rn_aeiv_max
         WRITE(numout,*)    '      minimum bound of aeiv coeff                     rn_aeiv_min = ', rn_aeiv_min
         WRITE(numout,*)    '      minimum bound of Structure Function             rn_SFmin    = ', rn_SFmin
         WRITE(numout,*)    '      maximum bound of Structure Function             rn_SFmax    = ', rn_SFmax
         WRITE(numout,*)    '         [set rnSFmin = rnSFmax = 1 to use aeiv = aeiv(x,y,t)]      '
         WRITE(numout,*)    '      option for terms in eddy energy budget          nn_eke_opt  = ', nn_eke_opt
         WRITE(numout,*)    '         default: just PE->EKE growth and linear dissipation        '
         SELECT CASE( nn_eke_opt )
         !
         CASE(   0  )
            WRITE(numout,*) '         default                                                     '
         CASE(   1  )
            WRITE(numout,*) '         default + advection                                         '
         CASE(   2  )
            WRITE(numout,*) '         default + advection + KE->EKE                               '
         CASE(  77  )
            WRITE(numout,*) '         ONLY advection by rossby waves (no growth or dissipation)   '
         CASE(  88  )
            WRITE(numout,*) '         ONLY advection by z-avg mean flow (no growth or dissipation)'
         CASE(  99  )
            WRITE(numout,*) '         ONLY Laplacian diffusion          (no growth or dissipation)'
         CASE(  100  )
            WRITE(numout,*) '         use ML model to calculate EKE'
         CASE DEFAULT
            CALL ctl_stop('ldf_eke: wrong choice nn_eke_opt, set at least to 0 (default)')
         END SELECT
         !
         WRITE(numout,*)    '      advection of energy by Rossby waves             ln_adv_wav  =  ', ln_adv_wav
         !
      ENDIF
      !                                ! allocate eke arrays
      ALLOCATE( eke_geom(A2D(nn_hls)) , eke_keS(A2D(0)) , r1_ekedis(A2D(0)) , STAT=ierr )
      !
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ldf_eke_init: failed to allocate arrays')
      !
      eke_keS(:,:) = 0._wp                    ! To avoid unset values in unused array if nn_eke_opt /= 2
      !                                       ! (array is unused but may still be mistakeningly requested in XML output)
      !
      SELECT CASE( nn_eke_dis )               ! Specification of linear dissipation
      !
      CASE(   0  )      !==  constant  ==!
         IF(lwp) WRITE(numout,*) '      linear EKE dissipation coef. = constant = ', rn_ekedis, ' days'
         r1_ekedis(:,:) = 1._wp / (rn_ekedis * rday)
         !
      CASE( -20  )      !== fixed horizontal shape read in file  ==!
         IF(lwp) WRITE(numout,*) '      linear EKE dissipation coef. = F(i,j) read in geom_diss_2D.nc file'
         CALL iom_open ( 'geom_diss_2D.nc', inum )
         CALL iom_get  ( inum, jpdom_auto, 'r1_ekedis', r1_ekedis, cd_type = 'T', psgn = 1._wp ) ! read in as time-scale in days...
         CALL iom_close( inum )
         DO_2D( 0, 0, 0, 0 )
            r1_ekedis(ji,jj) = 1._wp / (r1_ekedis(ji,jj) * rday)   ! ...convert rate in per second
         END_2D
         !
      CASE DEFAULT
         CALL ctl_stop('ldf_eke_init: wrong choice for nn_eke_dis, the option for linear dissipation in GEOMETRIC')
      END SELECT
      !
      CALL ldf_eke_rst( nit000, 'READ', Kbb )   !* read or initialize all required files
      !
      ! if using beta_plane, compute beta and f0 using values from the 1st domain
      IF( ln_beta_plane ) THEN
        IF( narea == 1 ) THEN
           zf0 = ff_f(1,1)
           zbeta = (ff_f(1,2) - zf0) / e2t(1,1)
        ELSE
           zf0   = HUGE(1._wp)
           zbeta = HUGE(1._wp)
        ENDIF
        CALL mpp_min( 'ldfeke', zf0 )
        CALL mpp_min( 'ldfeke', zbeta )
        IF(lwp) WRITE(numout,*)    '         beta plane option is    ln_beta_plane =', ln_beta_plane
        IF(lwp) WRITE(numout,*)    '         f0   = ', zf0,   '    s-1'
        IF(lwp) WRITE(numout,*)    '         beta = ', zbeta, 'm-1 s-1'
      ENDIF

      IF ( nn_eke_opt == 100 ) THEN
         IF ( .NOT. ln_ssdb ) CALL ctl_stop('ldf_eke_init: EKE to be predicted by an ML model, but ln_ssdb is not set')
      ENDIF

      !
      !
   END SUBROUTINE ldf_eke_init


   SUBROUTINE ldf_eke_rst( kt, cdrw, Kbb )
     !!---------------------------------------------------------------------
     !!                   ***  ROUTINE eke_rst  ***
     !!
     !! ** Purpose :   Read or write EKE file (eke) in restart file
     !!
     !! ** Method  :   use of IOM library
     !!              if the restart does not contain EKE, eke is set
     !!              according to rn_eke_init, and aeiu = aeiv = 10 m s^-2
     !!----------------------------------------------------------------------
     INTEGER                   , INTENT(in) ::   kt       ! ocean time-step
     INTEGER         , OPTIONAL, INTENT(in) ::   Kbb      ! time level indicator (needed for 'ht')
     CHARACTER(len=*)          , INTENT(in) ::   cdrw     ! "READ"/"WRITE" flag
     !
     INTEGER ::   jit, jk      ! dummy loop indices
     INTEGER ::   id1          ! local integer
     !!----------------------------------------------------------------------
     !
     IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise
        !                                   ! ---------------
        IF( ln_rstart ) THEN                   !* Read the restart file
           id1 = iom_varid( numror, 'eke_b'   , ldstop = .FALSE. )
           IF(lwp) write(numout,*) ' ldfeke init restart'
        ELSE                                   !* Start from rest with a non zero value (required)
           id1 = -1
        ENDIF
        !
        IF( id1 > 0 ) THEN        ! all required arrays exist
           IF(lwp) write(numout,*) ' ldfeke init from restart file'
           CALL iom_get( numror, jpdom_auto, 'eke_b', eke_geom(:,:), cd_type = 'T', psgn = 1._wp )
        ELSE                      ! one at least array is missing
           IF(lwp) write(numout,*) ' ldfeke init restart from namelist'
           eke_geom(:,:)  = rn_eke_init * ht(:,:,Kbb) * ssmask(:,:)
        ENDIF
        !
     ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
        !                                   ! -------------------
        IF(lwp) WRITE(numout,*) '---- eke-rst ----'
        CALL iom_rstput( kt, nitrst, numrow, 'eke_b', eke_geom(:,:)  )
        !
     ENDIF
     !
   END SUBROUTINE ldf_eke_rst

   !!======================================================================
END MODULE ldfeke
