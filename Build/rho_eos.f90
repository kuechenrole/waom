      MODULE rho_eos_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rho_eos
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rho_eos (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
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
      CALL wclock_on (ng, iNLM, 14)
      CALL rho_eos_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   OCEAN(ng) % t,                                 &
     &                   COUPLING(ng) % rhoA,                           &
     &                   COUPLING(ng) % rhoS,                           &
     &                   MIXING(ng) % bvf,                              &
     &                   MIXING(ng) % alpha,                            &
     &                   MIXING(ng) % beta,                             &
     &                   MIXING(ng) % alfaobeta,                        &
     &                   OCEAN(ng) % pden,                              &
     &                   OCEAN(ng) % rho)
      CALL wclock_off (ng, iNLM, 14)
      RETURN
      END SUBROUTINE rho_eos
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         rmask,                                   &
     &                         Hz,                                      &
     &                         z_r, z_w, t,                             &
     &                         rhoA, rhoS,                              &
     &                         bvf,                                     &
     &                         alpha, beta,                             &
     &                         alfaobeta,                               &
     &                         pden,                                    &
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_eoscoef
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k
      real(r8) :: SedDen, Tp, Tpr10, Ts, Tt, sqrtTs
      real(r8) :: bulk_dn, bulk_up, den_dn, den_up
      real(r8) :: cff, cff1, cff2
      real(r8), dimension(0:9) :: C
      real(r8), dimension(0:9) :: dCdT(0:9)
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: wrk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den1
      real (r8) :: pp, rz0, rz1, rz2, rz3, sau, ss, ctu, zu, deltas
      real (r8) :: sa, ct, p, r1
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
!=======================================================================
!  Nonlinear equation of state.  Notice that this equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!=======================================================================
!
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
!
            Tt=MAX(-5.0_r8,t(i,j,k,nrhs,itemp))
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
!
            dCdT(0)=Q01+Tt*(2.0_r8*Q02+Tt*(3.0_r8*Q03+Tt*(4.0_r8*Q04+   &
     &                      Tt*5.0_r8*Q05)))
            dCdT(1)=U01+Tt*(2.0_r8*U02+Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            dCdT(2)=V01+Tt*2.0_r8*V02
!
            den1(i,k)=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)
!
!  Compute d(den1)/d(S) and d(den1)/d(T) derivatives used in the
!  computation of thermal expansion and saline contraction
!  coefficients.
!
            Dden1DS(i,k)=C(1)+1.5_r8*C(2)*sqrtTs+2.0_r8*W00*Ts
            Dden1DT(i,k)=dCdT(0)+Ts*(dCdT(1)+sqrtTs*dCdT(2))
!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
!
            dCdT(3)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            dCdT(4)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            dCdT(5)=D01+Tt*2.0_r8*D02
            dCdT(6)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            dCdT(7)=F01+Tt*2.0_r8*F02
            dCdT(8)=G02+Tt*2.0_r8*G03
            dCdT(9)=H01+Tt*2.0_r8*H02
!
            bulk0(i,k)=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1(i,k)=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2(i,k)=C(8)+Ts*C(9)
            bulk (i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))
!
!  Compute d(bulk)/d(S) and d(bulk)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
            DbulkDS(i,k)=C(4)+sqrtTs*1.5_r8*C(5)-                       &
     &                   Tp*(C(7)+sqrtTs*1.5_r8*G00-Tp*C(9))
            DbulkDT(i,k)=dCdT(3)+Ts*(dCdT(4)+sqrtTs*dCdT(5))-           &
     &                   Tp*(dCdT(6)+Ts*dCdT(7)-                        &
     &                       Tp*(dCdT(8)+Ts*dCdT(9)))
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0_r8/(bulk(i,k)+Tpr10)
            den(i,k)=den1(i,k)*bulk(i,k)*cff
            den(i,k)=den(i,k)-1000.0_r8
            den(i,k)=den(i,k)*rmask(i,j)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  (rhoS) used in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=den(i,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points:
!
!                  bvf = - g/rho d(rho)/d(z).
!
!  The density anomaly difference is computed by lowering/rising the
!  water parcel above/below adiabatically at W-point depth "z_w".
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrT,IendT
            bulk_up=bulk0(i,k+1)-                                       &
     &              z_w(i,j,k)*(bulk1(i,k+1)-                           &
     &                          bulk2(i,k+1)*z_w(i,j,k))
            bulk_dn=bulk0(i,k  )-                                       &
     &              z_w(i,j,k)*(bulk1(i,k  )-                           &
     &                          bulk2(i,k  )*z_w(i,j,k))
            cff1=1.0_r8/(bulk_up+0.1_r8*z_w(i,j,k))
            cff2=1.0_r8/(bulk_dn+0.1_r8*z_w(i,j,k))
            den_up=cff1*(den1(i,k+1)*bulk_up)
            den_dn=cff2*(den1(i,k  )*bulk_dn)
            bvf(i,j,k)=-g*(den_up-den_dn)/                              &
     &                 (0.5_r8*(den_up+den_dn)*                         &
     &                  (z_r(i,j,k+1)-z_r(i,j,k)))
          END DO
        END DO
        DO i=IstrT,IendT
          bvf(i,j,0)=0.0_r8
          bvf(i,j,N(ng))=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
            Tpr10=0.1_r8*z_r(i,j,k)
!
!  Compute thermal expansion and saline contraction coefficients.
!
            cff=bulk(i,k)+Tpr10
            cff1=Tpr10*den1(i,k)
            cff2=bulk(i,k)*cff
            wrk(i,k)=(den(i,k)+1000.0_r8)*cff*cff
            Tcof(i,k)=-(DbulkDT(i,k)*cff1+                              &
     &                  Dden1DT(i,k)*cff2)
            Scof(i,k)= (DbulkDS(i,k)*cff1+                              &
     &                  Dden1DS(i,k)*cff2)
            alfaobeta(i,j,k)=Tcof(i,k)/Scof(i,k)
          END DO
          IF (k.eq.N(ng)) THEN
            DO i=IstrT,IendT
              cff=1.0_r8/wrk(i,N(ng))
              alpha(i,j)=cff*Tcof(i,N(ng))
              beta (i,j)=cff*Scof(i,N(ng))
            END DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Load "in situ" density anomaly (kg/m3 - 1000) and potential
!  density anomaly (kg/m3 - 1000) referenced to the surface into global
!  arrays. Notice that this is done in a separate (i,k) DO-loops to
!  facilitate the adjoint.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
            rho(i,j,k)=den(i,k)
            pden(i,j,k)=(den1(i,k)-1000.0_r8)
            pden(i,j,k)=pden(i,j,k)*rmask(i,j)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rho)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          pden)
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          alfaobeta)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          beta)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoS)
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          bvf)
      END IF
!
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rho, pden)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alfaobeta)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alpha, beta)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhoA, rhoS)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bvf)
      RETURN
      END SUBROUTINE rho_eos_tile
      END MODULE rho_eos_mod
