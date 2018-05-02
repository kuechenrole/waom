      MODULE lmd_skpp_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine determines the depth of surface  oceanic boundary      !
!  layer, hsbl, as the shallowest depth where the bulk Richardson      !
!  number is equal to the critical value, Ric.                         !
!                                                                      !
!  Then,  it computes the vertical mixing coefficients  within the     !
!  boundary layer. They depend on surface forcing and the magnitude    !
!  and gradient of interior mixing below  the boundary layer.  The     !
!  ocean interior is allowed to force the boundary layer through a     !
!  dependence of the nondimensional vertical shape function G(sigma)   !
!  and its vertical derivative at  sigma=1  on the interior  mixing    !
!  coefficients, and it vertical derivative at d=hsbl. The boundary    !
!  layer mixing coefficients are computed by matching these values.    !
!                                                                      !
! Reference:                                                           !
!                                                                      !
!  Large, W.G., J.C. McWilliams, and S.C. Doney, 1994: A Review        !
!    and model with a nonlocal boundary layer parameterization,        !
!    Reviews of Geophysics, 32,363-403.                                !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: lmd_skpp
      CONTAINS
!
!***********************************************************************
      SUBROUTINE lmd_skpp (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
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
      CALL lmd_skpp_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng),                                     &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v,                                &
     &                    OCEAN(ng) % pden,                             &
     &                    FORCES(ng) % srflx,                           &
     &                    FORCES(ng) % stflx,                           &
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
     &                    FORCES(ng) % sustr,                           &
     &                    FORCES(ng) % svstr,                           &
     &                    MIXING(ng) % alpha,                           &
     &                    MIXING(ng) % beta,                            &
     &                    MIXING(ng) % bvf,                             &
     &                    MIXING(ng) % ghats,                           &
     &                    MIXING(ng) % Akt,                             &
     &                    MIXING(ng) % Akv,                             &
     &                    MIXING(ng) % hsbl,                            &
     &                    MIXING(ng) % ksbl)
      RETURN
      END SUBROUTINE lmd_skpp
!
!***********************************************************************
      SUBROUTINE lmd_skpp_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
     &                          rmask,                                  &
     &                          f, Hz, z_r, z_w,                        &
     &                          u, v, pden,                             &
     &                          srflx, stflx,                           &
     &                          bustr, bvstr, sustr, svstr,             &
     &                          alpha,                                  &
     &                          beta,                                   &
     &                          bvf,                                    &
     &                          ghats,                                  &
     &                          Akt, Akv, hsbl, ksbl)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE shapiro_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: pden(LBi:,LBj:,:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: alpha(LBi:,LBj:)
      real(r8), intent(in) :: beta(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: hsbl(LBi:,LBj:)
      integer,  intent(out) :: ksbl(LBi:,LBj:)
      real(r8), intent(out) :: ghats(LBi:,LBj:,0:,:)
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
      real(r8), parameter :: small = 1.0E-20_r8
      real(r8) :: Gm, Gt, Gs, K_bl, Ribot, Ritop, Rk
      real(r8) :: Uk, Ustarb, Ustar3, Vk, Vtc
      real(r8) :: a1, a2, a3, cff, cff1, cff2,cff_up, cff_dn
      real(r8) :: depth, dK_bl, hekman, hmonob, sigma, zbl
      real(r8) :: zetahat, zetapar
      real(r8), dimension (IminS:ImaxS) :: Rref
      real(r8), dimension (IminS:ImaxS) :: Uref
      real(r8), dimension (IminS:ImaxS) :: Vref
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Bflux
      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: dV
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Bo
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Bosol
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Bfsfc
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gm1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gt1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gs1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Ustar
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGm1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGt1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGs1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: f1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: sl_dpth
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: swdk
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: wm
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: ws
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: zgrid
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
!  Initialize relevant parameters.
!-----------------------------------------------------------------------
!
      Vtc=lmd_Cv*SQRT(-lmd_betaT)/(SQRT(lmd_cs*lmd_epsilon)*            &
     &                             lmd_Ric*vonKar*vonKar)
!
!-----------------------------------------------------------------------
!  Get approximation of surface layer depth using "lmd_eps" and
!  boundary layer depth from previous time step.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          sl_dpth(i,j)=lmd_epsilon*(z_w(i,j,N(ng))-hsbl(i,j))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute turbulent friction velocity (m/s) "Ustar" from wind stress
!  at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Ustar(i,j)=SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+   &
     &                         (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
          Ustar(i,j)=Ustar(i,j)*rmask(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute surface turbulent buoyancy forcing "Bo" (m2/s3). Remove
!  incoming solar shortwave radiation because this contribution is
!  included in "Bosol".  Compute surface radiative buoyancy forcing
!  "Bosol" (m2/s3).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Bo(i,j)=g*(alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))-          &
     &               beta (i,j)*stflx(i,j,isalt))
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux (m2/s3) at W-points.  Notice that the
!  radiative bouyancy flux is distributed vertically using decay
!  function, swdk. Begin computation of nonlocal transport, "ghats".
!-----------------------------------------------------------------------
!
      DO k=0,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            zgrid(i,j)=z_w(i,j,N(ng))-z_w(i,j,k)
          END DO
        END DO
        CALL lmd_swfrac_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        -1.0_r8, zgrid, swdk)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Bflux(i,j,k)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
            Bflux(i,j,k)=Bflux(i,j,k)*rmask(i,j)
            cff=1.0_r8-(0.5_r8+SIGN(0.5_r8,Bflux(i,j,k)))
            ghats(i,j,k,itemp)=-cff*(stflx(i,j,itemp)-srflx(i,j)+       &
     &                               srflx(i,j)*(1.0_r8-swdk(i,j)))
            ghats(i,j,k,isalt)=cff*stflx(i,j,isalt)
          END DO
        END DO
      END DO
!
!=======================================================================
!  Compute bulk Richardson number "Rib" and then find depth of the
!  oceanic surface boundary layer "hsbl", such that Rib(hsbl)=Ric.
!=======================================================================
!
      DO j=Jstr,Jend
!
! Construct parabolic splines for vertical derivatives of potential
! density and velocity components at W-points.  FC is a scratch array.
!
        DO i=Istr,Iend
          FC(i,0)=0.0_r8
          dR(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*Hz(i,j,k+1)
            dR(i,k)=cff*(6.0_r8*(pden(i,j,k+1)-pden(i,j,k))-            &
     &                   Hz(i,j,k)*dR(i,k-1))
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+     &
     &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
     &                   Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+     &
     &                           v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))-    &
     &                   Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO
        DO i=Istr,Iend
          dR(i,N(ng))=0.0_r8
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            dR(i,k)=dR(i,k)-FC(i,k)*dR(i,k+1)
            dU(i,k)=dU(i,k)-FC(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-FC(i,k)*dV(i,k+1)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute bulk Richardson number "Rib" and then find depth of oceanic
!  surface boundary layer "hsbl".
!
!                  [Br - B(d)] * d
!     Rib(d) = ----------------------- ;     Rib(hsbl)=Ric     (1)
!              |Vr - V(d)|^2 + Vt(d)^2
!
!  where "Br" and "Vr" are the surface reference buoyancy and velocity
!  while "B(d)" and "V(d)" are the bouyancy and velocity at depth "d".
!
!  In the code below, the criterion "Rib(hsbl)=Ric" is reformulated
!  as follows:
!
!     Rib(d)       Ritop(d)
!     ------ = --------------- = 1                             (2)
!      Ric      Ric * Ribot(d)
!
!  where "Ritop" and "Ribot" are numerator and denominator in Eq. (1).
!  In its turn, Eq. (2) is rewritten in the following form:
!
!     FC(d) = Ritop(d) - Ric * Ribot(d) = 0                    (3)
!
!  That is, the planetary boundary layer extends to the depth where
!  the critical function "FC(d)" changes its sign.
!-----------------------------------------------------------------------
!
!  Compute potential density and velocity components surface reference
!  values.
!
        cff1=1.0_r8/3.0_r8
        cff2=1.0_r8/6.0_r8
        DO i=Istr,Iend
          Rref(i)=pden(i,j,N(ng))+                                      &
     &            Hz(i,j,N(ng))*(cff1*dR(i,N(ng))+cff2*dR(i,N(ng)-1))
          Uref(i)=0.5_r8*(u(i,j,N(ng),nstp)+u(i+1,j,N(ng),nstp))+       &
     &            Hz(i,j,N(ng))*(cff1*dU(i,N(ng))+cff2*dU(i,N(ng)-1))
          Vref(i)=0.5_r8*(v(i,j,N(ng),nstp)+v(i,j+1,N(ng),nstp))+       &
     &            Hz(i,j,N(ng))*(cff1*dV(i,N(ng))+cff2*dV(i,N(ng)-1))
        END DO
!
!  Compute turbulent velocity scales for momentum (wm) and tracers (ws).
!  Then, compute critical function (FC) for bulk Richardson number.
!
        DO i=Istr,Iend
          FC(i,N(ng))=0.0_r8
          DO k=N(ng),1,-1
            depth=z_w(i,j,N(ng))-z_w(i,j,k-1)
            IF (Bflux(i,j,k-1).lt.0.0_r8) THEN
              sigma=MIN(sl_dpth(i,j),depth)
            ELSE
              sigma=depth
            END IF
            Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
            zetahat=vonKar*sigma*Bflux(i,j,k-1)
            zetapar=zetahat/(Ustar3+small)
            IF (zetahat.ge.0.0_r8) THEN                         ! stable
              wm(i,j)=vonKar*Ustar(i,j)/(1.0_r8+5.0_r8*zetapar)
              ws(i,j)=wm(i,j)
            ELSE                                              ! unstable
              IF (zetapar.gt.lmd_zetam) THEN
                wm(i,j)=vonKar*Ustar(i,j)*                              &
     &                  (1.0_r8-16.0_r8*zetapar)**0.25_r8
              ELSE
                wm(i,j)=vonKar*(lmd_am*Ustar3-lmd_cm*zetahat)**r3
              END IF
              IF (zetapar.gt.lmd_zetas) THEN
                ws(i,j)=vonKar*Ustar(i,j)*                              &
     &                  (1.0_r8-16.0_r8*zetapar)**0.5_r8
              ELSE
                ws(i,j)=vonKar*(lmd_as*Ustar3-lmd_cs*zetahat)**r3
              END IF
            END IF
!
            Rk=pden(i,j,k)-                                             &
     &         Hz(i,j,k)*(cff1*dR(i,k-1)+cff2*dR(i,k))
            Uk=0.5_r8*(u(i,j,k,nstp)+u(i+1,j,k,nstp))-                  &
     &         Hz(i,j,k)*(cff1*dU(i,k-1)+cff2*dU(i,k))
            Vk=0.5_r8*(v(i,j,k,nstp)+v(i,j+1,k,nstp))-                  &
     &         Hz(i,j,k)*(cff1*dV(i,k-1)+cff2*dV(i,k))
!
            Ritop=-gorho0*(Rref(i)-Rk)*depth
            Ribot=(Uref(i)-Uk)**2+(Vref(i)-Vk)**2+                      &
     &            Vtc*depth*ws(i,j)*SQRT(ABS(bvf(i,j,k-1)))
            FC(i,k-1)=Ritop-lmd_Ric*Ribot
          END DO
        END DO
!
! Linearly interpolate to find "hsbl" where Rib/Ric=1.
!
        DO i=Istr,Iend
          ksbl(i,j)=1
          hsbl(i,j)=z_w(i,j,1)
        END DO
        DO k=N(ng),2,-1
          DO i=Istr,Iend
            IF ((ksbl(i,j).eq.1).and.(FC(i,k-1).gt.0.0_r8)) THEN
              hsbl(i,j)=(z_w(i,j,k)*FC(i,k-1)-z_w(i,j,k-1)*FC(i,k))/    &
     &                  (FC(i,k-1)-FC(i,k))
              ksbl(i,j)=k
            END IF
          END DO
        END DO
      END DO
!
!  Compute total buoyancy flux at surface boundary layer depth,
!  "Bfsfc".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hsbl(i,j)
          zgrid(i,j)=zgrid(i,j)*rmask(i,j)
        END DO
      END DO
      CALL lmd_swfrac_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      -1.0_r8, zgrid, swdk)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Bfsfc(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
          Bfsfc(i,j)=Bfsfc(i,j)*rmask(i,j)
        END DO
      END DO
!
!  Under neutral and stable conditions, the depth of the surface
!  boundary layer is required to be less than Ekman and Monin-Obukov
!  depths.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF ((Ustar(i,j).gt.0.0_r8).and.(Bfsfc(i,j).gt.0.0_r8)) THEN
            hekman=lmd_cekman*Ustar(i,j)/MAX(ABS(f(i,j)),eps)
            hmonob=lmd_cmonob*Ustar(i,j)*Ustar(i,j)*Ustar(i,j)/         &
     &             MAX(vonKar*Bfsfc(i,j),eps)
            hsbl(i,j)=(z_w(i,j,N(ng))-                                  &
     &                 MIN(hekman,hmonob,z_w(i,j,N(ng))-hsbl(i,j)))
          END IF
          hsbl(i,j)=MIN(hsbl(i,j),z_w(i,j,N(ng)))
          hsbl(i,j)=MAX(hsbl(i,j),z_w(i,j,0))
          hsbl(i,j)=hsbl(i,j)*rmask(i,j)
        END DO
      END DO
!
!  Apply gradient or periodic boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  hsbl)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    hsbl)
!
!  Apply Shapiro filter
!
      CALL shapiro2d_tile (ng, tile, iNLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     rmask,                                       &
     &                     hsbl)
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          hsbl(i,j)=MIN(hsbl(i,j),z_w(i,j,N(ng)))
          hsbl(i,j)=MAX(hsbl(i,j),z_w(i,j,0))
          hsbl(i,j)=hsbl(i,j)*rmask(i,j)
        END DO
      END DO
!
!  Apply gradient or periodic boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  hsbl)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    hsbl)
!
!  Find new boundary layer index "ksbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ksbl(i,j)=1
          DO k=N(ng),2,-1
            IF ((ksbl(i,j).eq.1).and.(z_w(i,j,k-1).lt.hsbl(i,j))) THEN
              ksbl(i,j)=k
            END IF
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux at final surface boundary layer depth.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hsbl(i,j)
          zgrid(i,j)=zgrid(i,j)*rmask(i,j)
        END DO
      END DO
      CALL lmd_swfrac_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      -1.0_r8, zgrid, swdk)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Bfsfc(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
          Bfsfc(i,j)=Bfsfc(i,j)*rmask(i,j)
        END DO
      END DO
!
!=======================================================================
!  Compute vertical mixing coefficients within the planetary boundary
!  layer.
!=======================================================================
!
!  Compute tubulent velocity scales (wm,ws) at "hsbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          sl_dpth(i,j)=lmd_epsilon*(z_w(i,j,N(ng))-hsbl(i,j))
          IF (Bfsfc(i,j).gt.0.0_r8) THEN
            cff=1.0_r8
          ELSE
            cff=lmd_epsilon
          END IF
          sigma=cff*(z_w(i,j,N(ng))-hsbl(i,j))
          Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
          zetahat=vonKar*sigma*Bfsfc(i,j)
          zetapar=zetahat/(Ustar3+small)
          IF (zetahat.ge.0.0_r8) THEN                           ! stable
            wm(i,j)=vonKar*Ustar(i,j)/(1.0_r8+5.0_r8*zetapar)
            ws(i,j)=wm(i,j)
          ELSE                                                ! unstable
            IF (zetapar.gt.lmd_zetam) THEN
              wm(i,j)=vonKar*Ustar(i,j)*                                &
     &                (1.0_r8-16.0_r8*zetapar)**0.25_r8
            ELSE
              wm(i,j)=vonKar*(lmd_am*Ustar3-lmd_cm*zetahat)**r3
            END IF
            IF (zetapar.gt.lmd_zetas) THEN
              ws(i,j)=vonKar*Ustar(i,j)*                                &
     &                (1.0_r8-16.0_r8*zetapar)**0.5_r8
            ELSE
              ws(i,j)=vonKar*(lmd_as*Ustar3-lmd_cs*zetahat)**r3
            END IF
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute nondimensional shape function Gx(sigma) in terms of the
!  interior diffusivities at sigma=1 (Gm1, Gs1, Gt1) and its vertical
!  derivative evaluated "hsbl" via interpolation.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          f1(i,j)=5.0_r8*MAX(0.0_r8,Bfsfc(i,j))*vonKar/                 &
     &            (Ustar(i,j)*Ustar(i,j)*Ustar(i,j)*Ustar(i,j)+eps)
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zbl=z_w(i,j,N(ng))-hsbl(i,j)
          IF (hsbl(i,j).gt.z_w(i,j,1)) THEN
            k=ksbl(i,j)
            cff=1.0_r8/(z_w(i,j,k)-z_w(i,j,k-1))
            cff_dn=cff*(hsbl(i,j)-z_w(i,j,k-1))
            cff_up=cff*(z_w(i,j,k)-hsbl(i,j))
!
!  Compute nondimensional shape function for viscosity "Gm1" and its
!  vertical derivative "dGm1dS" evaluated at "hsbl".
!
            K_bl=cff_dn*Akv(i,j,k)+cff_up*Akv(i,j,k-1)
            dK_bl=cff*(Akv(i,j,k)-Akv(i,j,k-1))
            Gm1(i,j)=K_bl/(zbl*wm(i,j)+eps)
            Gm1(i,j)=Gm1(i,j)*rmask(i,j)
            dGm1dS(i,j)=MIN(0.0_r8,-dK_bl/(wm(i,j)+eps)-K_bl*f1(i,j))
!
!  Compute nondimensional shape function for diffusion of temperature
!  "Gt1" and its vertical derivative "dGt1dS" evaluated at "hsbl".
!
            K_bl=cff_dn*Akt(i,j,k,itemp)+cff_up*Akt(i,j,k-1,itemp)
            dK_bl=cff*(Akt(i,j,k,itemp)-Akt(i,j,k-1,itemp))
            Gt1(i,j)=K_bl/(zbl*ws(i,j)+eps)
            Gt1(i,j)=Gt1(i,j)*rmask(i,j)
            dGt1dS(i,j)=MIN(0.0_r8,-dK_bl/(ws(i,j)+eps)-K_bl*f1(i,j))
!
!  Compute nondimensional shape function for diffusion of salinity
!  "Gs1" and its vertical derivative "dGs1dS" evaluated at "hsbl".
!
            K_bl=cff_dn*Akt(i,j,k,isalt)+cff_up*Akt(i,j,k-1,isalt)
            dK_bl=cff*(Akt(i,j,k,isalt)-Akt(i,j,k-1,isalt))
            Gs1(i,j)=K_bl/(zbl*ws(i,j)+eps)
            Gs1(i,j)=Gs1(i,j)*rmask(i,j)
            dGs1dS(i,j)=MIN(0.0_r8,-dK_bl/(ws(i,j)+eps)-K_bl*f1(i,j))
          ELSE
!
!  If the surface boundary layer extends to the bottom, assume that
!  the neutral boundary layer similarity theory holds at the bottom.
!
            ksbl(i,j)=0
!
!  Compute nondimensional bottom shape function for viscosity.
!
            Ustarb=SQRT(SQRT((0.5_r8*(bustr(i,j)+bustr(i+1,j)))**2+     &
     &                       (0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))**2))
            Ustarb=Ustarb*rmask(i,j)
            dK_bl=vonKar*Ustarb
            K_bl=dK_bl*(hsbl(i,j)-z_w(i,j,0))
            Gm1(i,j)=K_bl/(zbl*wm(i,j)+eps)
            Gm1(i,j)=Gm1(i,j)*rmask(i,j)
            dGm1dS(i,j)=MIN(0.0_r8,-dK_bl/(wm(i,j)+eps)-K_bl*f1(i,j))
!
!  Compute nondimensional bottom shape function for diffusion of
!  temperature.
!
            Gt1(i,j)=K_bl/(zbl*ws(i,j)+eps)
            Gt1(i,j)=Gt1(i,j)*rmask(i,j)
            dGt1dS(i,j)=MIN(0.0_r8,-dK_bl/(ws(i,j)+eps)-K_bl*f1(i,j))
!
!  Compute nondimensional bottom shape function for diffusion of
!  salinity.
!
            Gs1(i,j)=Gt1(i,j)
            dGs1dS(i,j)=dGt1dS(i,j)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute surface boundary layer mixing coefficients.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            zbl=z_w(i,j,N(ng))-hsbl(i,j)
            IF (k.gt.ksbl(i,j)) THEN
!
!  Compute turbulent velocity scales at vertical W-points.
!
              depth=z_w(i,j,N(ng))-z_w(i,j,k)
              IF (Bflux(i,j,k).lt.0.0_r8) THEN
                sigma=MIN(sl_dpth(i,j),depth)
              ELSE
                sigma=depth
              END IF
              Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
              zetahat=vonKar*sigma*Bflux(i,j,k)
              zetapar=zetahat/(Ustar3+small)
              IF (zetahat.ge.0.0_r8) THEN                       ! stable
                wm(i,j)=vonKar*Ustar(i,j)/(1.0_r8+5.0_r8*zetapar)
                ws(i,j)=wm(i,j)
              ELSE                                            ! unstable
                IF (zetapar.gt.lmd_zetam) THEN
                  wm(i,j)=vonKar*Ustar(i,j)*                            &
     &                    (1.0_r8-16.0_r8*zetapar)**0.25_r8
                ELSE
                  wm(i,j)=vonKar*(lmd_am*Ustar3-lmd_cm*zetahat)**r3
                END IF
                IF (zetapar.gt.lmd_zetas) THEN
                  ws(i,j)=vonKar*Ustar(i,j)*                            &
     &                    (1.0_r8-16.0_r8*zetapar)**0.5_r8
                ELSE
                  ws(i,j)=vonKar*(lmd_as*Ustar3-lmd_cs*zetahat)**r3
                END IF
              END IF
!
!  Set polynomial coefficients for shape function.
!
              sigma=depth/(zbl+eps)
              sigma=sigma*rmask(i,j)
              a1=sigma-2.0_r8
              a2=3.0_r8-2.0_r8*sigma
              a3=sigma-1.0_r8
!
!  Compute nondimesional shape functions.
!
              Gm=a1+a2*Gm1(i,j)+a3*dGm1dS(i,j)
              Gt=a1+a2*Gt1(i,j)+a3*dGt1dS(i,j)
              Gs=a1+a2*Gs1(i,j)+a3*dGs1dS(i,j)
!
!  Compute boundary layer mixing coefficients, combine them
!  with interior mixing coefficients.
!
              Akv(i,j,k)=depth*wm(i,j)*(1.0_r8+sigma*Gm)
              Akt(i,j,k,itemp)=depth*ws(i,j)*(1.0_r8+sigma*Gt)
              AKt(i,j,k,isalt)=depth*ws(i,j)*(1.0_r8+sigma*Gs)
!
!  Compute boundary layer nonlocal transport (m/s2).
!
              cff=lmd_Cg*(1.0_r8-(0.5_r8+SIGN(0.5_r8,Bflux(i,j,k))))/   &
     &            (zbl*ws(i,j)+eps)
              ghats(i,j,k,itemp)=cff*ghats(i,j,k,itemp)
              ghats(i,j,k,isalt)=cff*ghats(i,j,k,isalt)
!
!  Set vertical mixing coefficients to fit neutral log layer
!  similarity theory.
!
            ELSE
              ghats(i,j,k,itemp)=0.0_r8
              ghats(i,j,k,isalt)=0.0_r8
            END IF
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE lmd_skpp_tile
      END MODULE lmd_skpp_mod
