      MODULE vorticity_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes relative (s-1) and  potential (m-1 s-1)       !
!  vorticity for an adiabatic Boussinesq fluid where the potential     !
!  density is conserved:                                               !
!                                                                      !
!    pvor = 1/rho0 dot_product(avor, grad(pden))                       !
!                                                                      !
!  where "avor" is the absolute (relative plus planetary) vorticity    !
!  and "pden" is the potential density (a conserved quantity).         !
!                                                                      !
!    avor = rvor + f                                                   !
!                                                                      !
!  In curvilinear coordinates, the vertical component of relative      !
!  vorticity and potential vorticity are:                              !
!  are:                                                                !
!                                                                      !
!    rvor = mn * [d(v/n)/d(xi) - d(u/m)/d(eta)]                        !
!                                                                      !
!    pvor = mn/rho0 * [f/mn +                                          !
!                                                                      !
!                      d(v/n)/d(xi) - d(u/m)/d(eta)] * d(pden)/d(z) +  !
!                                                                      !
!           1/rho0 * [1/n d(pden)/d(eta) d(u)/d(z) -                   !
!                                                                      !
!                     1/m d(pden)/d(xi)  d(v)/d(z)]                    !
!                                                                      !
!  In addition, the vertically integrated (shallow water) relative     !
!  and potential vorticity are computed.                               !
!                                                                      !
!  The relative and potential vorticity is discretized at horizontal   !
!  PSI-points and vertical RHO-points.                                 !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: vorticity
      PUBLIC  :: vorticity_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE vorticity (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_average
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 5)
      CALL vorticity_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng), nrhs(ng),                          &
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % fomn,                             &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % om_u,                             &
     &                     GRID(ng) % on_v,                             &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % zice,                             &
     &                     GRID(ng) % z_r,                              &
     &                     OCEAN(ng) % pden,                            &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
     &                     OCEAN(ng) % zeta,                            &
     &                     AVERAGE(ng) % avgpvor3d,                     &
     &                     AVERAGE(ng) % avgrvor3d,                     &
     &                     AVERAGE(ng) % avgpvor2d,                     &
     &                     AVERAGE(ng) % avgrvor2d)
      CALL wclock_off (ng, iNLM, 5)
      RETURN
      END SUBROUTINE vorticity
!
!***********************************************************************
      SUBROUTINE vorticity_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kout, nout,                            &
     &                           pmask, umask, vmask,                   &
     &                           fomn, h, om_u, on_v, pm, pn,           &
     &                           zice,                                  &
     &                           z_r, pden, u, v,                       &
     &                           ubar, vbar, zeta,                      &
     &                           pvor, rvor,                            &
     &                           pvor_bar, rvor_bar)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_p2d_tile
      USE exchange_3d_mod, ONLY : exchange_p3d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kout, nout
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: fomn(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: pden(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(out) :: pvor_bar(LBi:,LBj:)
      real(r8), intent(out) :: rvor_bar(LBi:,LBj:)
      real(r8), intent(out) :: pvor(LBi:,LBj:,:)
      real(r8), intent(out) :: rvor(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: k, k1, k2
      real(r8) :: cff
      real(r8) :: dVdx_p, dUde_p
      real(r8) :: dRde_pr, dRdx_pr, dRdz_pr, dUdz_pr, dVdz_pr
      real(r8) :: orho0
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dRde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dRdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dUde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dVdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dRdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dUdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dVdz
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Compute 3D relative and potential vorticity.
!-----------------------------------------------------------------------
!
!  Compute horizontal and vertical gradients.  Notice the recursive
!  blocking sequence for vertical placement of the gradients is:
!
!      dRdz,dUdz,dVdz(:,:,k1) k-1/2   W-points
!      dRdz,dUdz,dVdz(:,:,k2) k+1/2   W-points
!
      orho0=1.0_r8/rho0
      k2=1
      K_LOOP : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.gt.0) THEN
          DO j=Jstr-1,JendR
            DO i=Istr,IendR
              cff=0.5_r8*(pm(i,j)+pm(i-1,j))
              cff=cff*umask(i,j)
              dRdx(i,j)=cff*(pden(i  ,j,k)-                             &
     &                       pden(i-1,j,k))
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=Istr-1,IendR
              cff=0.5_r8*(pn(i,j)+pn(i,j-1))
              cff=cff*vmask(i,j)
              dRde(i,j)=cff*(pden(i,j  ,k)-                             &
     &                       pden(i,j-1,k))
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=Istr,IendR
              dUde(i,j)=om_u(i,j  )*u(i,j  ,k,nout)-                    &
     &                  om_u(i,j-1)*u(i,j-1,k,nout)
              dUde(i,j)=dUde(i,j)*pmask(i,j)
              dVdx(i,j)=on_v(i  ,j)*v(i  ,j,k,nout)-                    &
     &                  on_v(i-1,j)*v(i-1,j,k,nout)
              dVdx(i,j)=dVdx(i,j)*pmask(i,j)
            END DO
          END DO
        END IF
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,JendR
            DO i=Istr-1,IendR
              dRdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstr-1,JendR
            DO i=Istr,IendR
              dUdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=Istr-1,IendR
              dVdz(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,JendR
            DO i=Istr-1,IendR
              cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              dRdz(i,j,k2)=cff*(pden(i,j,k+1)-                          &
     &                          pden(i,j,k  ))
            END DO
          END DO
          DO j=Jstr-1,JendR
            DO i=Istr,IendR
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)-z_r(i-1,j,k)+          &
     &                            z_r(i  ,j,k+1)-z_r(i  ,j,k)))
              dUdz(i,j,k2)=cff*(u(i,j,k+1,nout)-                        &
     &                          u(i,j,k  ,nout))
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=Istr-1,IendR
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)-z_r(i,j-1,k)+          &
     &                            z_r(i,j  ,k+1)-z_r(i,j  ,k)))
              dVdz(i,j,k2)=cff*(v(i,j,k+1,nout)-                        &
     &                          v(i,j,k  ,nout))
            END DO
          END DO
        END IF
!
!  Compute relative vorticity (second-1) and potential vorticity
!  (meter-1 second-1) at horizontal PSI-points and vertical RHO-points.
!
        IF (k.gt.0) THEN
          DO j=Jstr,JendR
            DO i=Istr,IendR
              cff=pm(i,j)*pn(i,j)
              dRde_pr=dRde(i-1,j  )+dRde(i,j)
              dRdx_pr=dRdx(i  ,j-1)+dRdx(i,j)
              dRdz_pr=0.125_r8*(dRdz(i-1,j-1,k1)+dRdz(i-1,j-1,k2)+      &
     &                          dRdz(i  ,j-1,k1)+dRdz(i  ,j-1,k2)+      &
     &                          dRdz(i-1,j  ,k1)+dRdz(i-1,j  ,k2)+      &
     &                          dRdz(i  ,j  ,k1)+dRdz(i  ,j  ,k2))
              dUdz_pr=dUdz(i  ,j-1,k1)+dUdz(i  ,j-1,k2)+                &
     &                dUdz(i  ,j  ,k1)+dUdz(i  ,j  ,k2)
              dVdz_pr=dVdz(i-1,j  ,k1)+dVdz(i-1,j  ,k2)+                &
     &                dVdz(i  ,j  ,k1)+dVdz(i  ,j  ,k2)
              rvor(i,j,k)=cff*(dVdx(i,j)-dUde(i,j))
              pvor(i,j,k)=orho0*                                        &
     &                    (cff*dRdz_pr*(fomn(i,j)+                      &
     &                                  dVdx(i,j)-dUde(i,j))+           &
     &                     0.125_r8*(dUdz_pr*dRde_pr-dVdz_pr*dRdx_pr))
            END DO
          END DO
        END IF
      END DO K_LOOP
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          pvor)
        CALL exchange_p3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rvor)
      END IF
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pvor,                                         &
     &                    rvor)
!
!-----------------------------------------------------------------------
!  Compute 2D relative and potential vorticity.
!-----------------------------------------------------------------------
!
!  Compute vertically-integrated relative vorticity (second-1) and
!  potential vorticity (meter-1 second-1) at PSI-points.
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          cff=pm(i,j)*pn(i,j)
          dVdx_p=on_v(i  ,j)*vbar(i  ,j,kout)-                          &
     &           on_v(i-1,j)*vbar(i-1,j,kout)
          dVdx_p=dVdx_p*pmask(i,j)
          dUde_p=om_u(i,j  )*ubar(i,j  ,kout)-                          &
     &           om_u(i,j-1)*ubar(i,j-1,kout)
          dUde_p=dUde_p*pmask(i,j)
          rvor_bar(i,j)=cff*(dVdx_p-dUde_p)
          pvor_bar(i,j)=cff*((fomn(i,j)+dVdx_p-dUde_p)/                 &
     &                       (h(i,j)-ABS(zice(i,j))+zeta(i,j,kout)))
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pvor_bar)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rvor_bar)
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pvor_bar,                                     &
     &                    rvor_bar)
      RETURN
      END SUBROUTINE vorticity_tile
      END MODULE vorticity_mod
