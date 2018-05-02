      MODULE set_tides_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group        Robert Hetland   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine adds tidal elevation (m) and tidal currents (m/s) to   !
!  sea surface height and 2D momentum climatologies, respectively.     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_tides
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_tides (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_tides
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
      CALL wclock_on (ng, iNLM, 11)
      CALL set_tides_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     NTC(ng),                                     &
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
     &                     TIDES(ng) % SSH_Tamp,                        &
     &                     TIDES(ng) % SSH_Tphase,                      &
     &                     TIDES(ng) % UV_Tangle,                       &
     &                     TIDES(ng) % UV_Tphase,                       &
     &                     TIDES(ng) % UV_Tmajor,                       &
     &                     TIDES(ng) % UV_Tminor,                       &
     &                     TIDES(ng) % Tperiod)
      CALL wclock_off (ng, iNLM, 11)
      RETURN
      END SUBROUTINE set_tides
!
!***********************************************************************
      SUBROUTINE set_tides_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NTC,                                   &
     &                           angler,                                &
     &                           rmask, umask, vmask,                   &
     &                           SSH_Tamp, SSH_Tphase,                  &
     &                           UV_Tangle, UV_Tphase,                  &
     &                           UV_Tmajor, UV_Tminor,                  &
     &                           Tperiod)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_ncparam
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_boundary
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variables declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: NTC
!
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: Tperiod(MTC)
      real(r8), intent(in) :: SSH_Tamp(LBi:,LBj:,:)
      real(r8), intent(in) :: SSH_Tphase(LBi:,LBj:,:)
      real(r8), intent(in) :: UV_Tangle(LBi:,LBj:,:)
      real(r8), intent(in) :: UV_Tmajor(LBi:,LBj:,:)
      real(r8), intent(in) :: UV_Tminor(LBi:,LBj:,:)
      real(r8), intent(in) :: UV_Tphase(LBi:,LBj:,:)
!
!  Local variables declarations.
!
      logical :: update
      integer :: ILB, IUB, JLB, JUB
      integer :: i, itide, j
      real(r8) :: Cangle, Cphase, Sangle, Sphase
      real(r8) :: angle, cff, phase, omega, ramp
      real(r8) :: bry_cor, bry_pgr, bry_str, bry_val
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Etide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Utide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vtide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
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
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      ILB=BOUNDS(ng)%LBi(-1)
      IUB=BOUNDS(ng)%UBi(-1)
      JLB=BOUNDS(ng)%LBj(-1)
      JUB=BOUNDS(ng)%UBj(-1)
!
!=======================================================================
!  Process tidal forcing used at the lateral boundaries.
!=======================================================================
!
!  In refinement nesting applications, tidal forcing is only necessary
!  in the coarser large scale grid.
!
      NEEDED : IF (.not.(RefinedGrid(ng).and.                           &
     &                    RefineScale(ng).gt.0)) THEN
!
!  Set time-ramping parameter.
!
        ramp=TANH((tdays(ng)-dstart)/1.0_r8)
!
!-----------------------------------------------------------------------
!  Add tidal elevation (m) to sea surface height climatology.
!-----------------------------------------------------------------------
!
      Etide(:,:)=0.0_r8
      cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
      DO itide=1,NTC
        IF (Tperiod(itide).gt.0.0_r8) THEN
          omega=cff/Tperiod(itide)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Etide(i,j)=Etide(i,j)+                                    &
     &                   ramp*SSH_Tamp(i,j,itide)*                      &
     &                   COS(omega-SSH_Tphase(i,j,itide))
              Etide(i,j)=Etide(i,j)*rmask(i,j)
            END DO
          END DO
        END IF
      END DO
!
!  If appropriate, load tidal forcing into boundary arrays.  The "zeta"
!  boundary arrays are important for the Flather or reduced physics
!  boundary conditions for 2D momentum. To avoid having two boundary
!  points for these arrays, the values of "zeta_west" and "zeta_east"
!  are averaged at u-points.  Similarly, the values of "zeta_south"
!  and "zeta_north" is averaged at v-points. Noticed that these
!  arrays are also used for the clamped conditions for free-surface.
!  This averaging is less important for that type ob boundary
!  conditions.
!
        IF (LBC(iwest,isFsur,ng)%acquire.or.                            &
     &      LBC(iwest,isUbar,ng)%acquire.or.                            &
     &      LBC(iwest,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrR,JendR
              BOUNDARY(ng)%zeta_west(j)=BOUNDARY(ng)%zeta_west(j)+      &
     &                                  0.5_r8*(Etide(Istr-1,j)+        &
     &                                          Etide(Istr  ,j))
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_west)
        END IF
!
        IF (LBC(ieast,isFsur,ng)%acquire.or.                            &
     &      LBC(ieast,isUbar,ng)%acquire.or.                            &
     &      LBC(ieast,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrR,JendR
              BOUNDARY(ng)%zeta_east(j)=BOUNDARY(ng)%zeta_east(j)+      &
     &                                  0.5_r8*(Etide(Iend  ,j)+        &
     &                                          Etide(Iend+1,j))
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_east)
        END IF
!
        IF (LBC(isouth,isFsur,ng)%acquire.or.                           &
     &      LBC(isouth,isUbar,ng)%acquire.or.                           &
     &      LBC(isouth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrR,IendR
              BOUNDARY(ng)%zeta_south(i)=BOUNDARY(ng)%zeta_south(i)+    &
     &                                   0.5_r8*(Etide(i,Jstr-1)+       &
     &                                           Etide(i,Jstr  ))
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_south)
        END IF
!
        IF (LBC(inorth,isFsur,ng)%acquire.or.                           &
     &      LBC(inorth,isUbar,ng)%acquire.or.                           &
     &      LBC(inorth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrR,IendR
              BOUNDARY(ng)%zeta_north(i)=BOUNDARY(ng)%zeta_north(i)+    &
     &                                   0.5_r8*(Etide(i,Jend  )+       &
     &                                           Etide(i,Jend+1))
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_north)
        END IF
!
!-----------------------------------------------------------------------
!  Add tidal currents (m/s) to 2D momentum climatologies.
!-----------------------------------------------------------------------
!
        Utide(:,:)=0.0_r8
        Vtide(:,:)=0.0_r8
        cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
        DO itide=1,NTC
          IF (Tperiod(itide).gt.0.0_r8) THEN
            omega=cff/Tperiod(itide)
            DO j=MIN(JstrR,Jstr-1),JendR
              DO i=MIN(IstrR,Istr-1),IendR
                angle=UV_Tangle(i,j,itide)-angler(i,j)
                Cangle=COS(angle)
                Sangle=SIN(angle)
                phase=omega-UV_Tphase(i,j,itide)
                Cphase=COS(phase)
                Sphase=SIN(phase)
                Uwrk(i,j)=UV_Tmajor(i,j,itide)*Cangle*Cphase-           &
     &                    UV_Tminor(i,j,itide)*Sangle*Sphase
                Vwrk(i,j)=UV_Tmajor(i,j,itide)*Sangle*Cphase+           &
     &                    UV_Tminor(i,j,itide)*Cangle*Sphase
              END DO
            END DO
            DO j=JstrR,JendR
              DO i=Istr,IendR
                Utide(i,j)=Utide(i,j)+                                  &
     &                     ramp*0.5_r8*(Uwrk(i-1,j)+Uwrk(i,j))
                Utide(i,j)=Utide(i,j)*umask(i,j)
              END DO
            END DO
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                Vtide(i,j)=(Vtide(i,j)+                                 &
     &                      ramp*0.5_r8*(Vwrk(i,j-1)+Vwrk(i,j)))
                Vtide(i,j)=Vtide(i,j)*vmask(i,j)
              END DO
            END DO
          END IF
        END DO
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
        IF (Lm2CLM(ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              CLIMA(ng)%ubarclm(i,j)=CLIMA(ng)%ubarclm(i,j)+Utide(i,j)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              CLIMA(ng)%vbarclm(i,j)=CLIMA(ng)%vbarclm(i,j)+Vtide(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              CLIMA(ng)%ubarclm)
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              CLIMA(ng)%vbarclm)
          END IF
          CALL mp_exchange2d (ng, tile, iNLM, 2,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        CLIMA(ng)%ubarclm,                        &
     &                        CLIMA(ng)%vbarclm)
        END IF
!
!  If appropriate, load tidal forcing into boundary arrays.
!
        IF (LBC(iwest,isUbar,ng)%acquire.and.                           &
     &      LBC(iwest,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrR,JendR
              BOUNDARY(ng)%ubar_west(j)=BOUNDARY(ng)%ubar_west(j)+      &
     &                                  Utide(Istr,j)
            END DO
            DO j=Jstr,JendR
              BOUNDARY(ng)%vbar_west(j)=BOUNDARY(ng)%vbar_west(j)+      &
     &                                  Vtide(Istr-1,j)
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_west)
          CALL mp_boundary (ng, iNLM, Jstr,  JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_west)
        END IF
!
        IF (LBC(ieast,isUbar,ng)%acquire.and.                           &
     &      LBC(ieast,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrR,JendR
              BOUNDARY(ng)%ubar_east(j)=BOUNDARY(ng)%ubar_east(j)+      &
     &                                  Utide(Iend+1,j)
            END DO
            DO j=Jstr,JendR
              BOUNDARY(ng)%vbar_east(j)=BOUNDARY(ng)%vbar_east(j)+      &
     &                                  Vtide(Iend+1,j)
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_east)
          CALL mp_boundary (ng, iNLM, Jstr,  JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_east)
        END IF
!
        IF (LBC(isouth,isUbar,ng)%acquire.and.                          &
     &      LBC(isouth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,IendR
              BOUNDARY(ng)%ubar_south(i)=BOUNDARY(ng)%ubar_south(i)+    &
     &                                   Utide(i,Jstr-1)
            END DO
            DO i=IstrR,IendR
              BOUNDARY(ng)%vbar_south(i)=BOUNDARY(ng)%vbar_south(i)+    &
     &                                   Vtide(i,Jstr)
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, Istr,  IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_south)
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_south)
        END IF
!
        IF (LBC(inorth,isUbar,ng)%acquire.and.                          &
     &      LBC(inorth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,IendR
              BOUNDARY(ng)%ubar_north(i)=BOUNDARY(ng)%ubar_north(i)+    &
     &                                   Utide(i,Jend+1)
            END DO
            DO i=IstrR,IendR
              BOUNDARY(ng)%vbar_north(i)=BOUNDARY(ng)%vbar_north(i)+    &
     &                                   Vtide(i,Jend+1)
            END DO
            update=.TRUE.
          END IF
          CALL mp_boundary (ng, iNLM, Istr,  IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_north)
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_north)
        END IF
      END IF NEEDED
      RETURN
      END SUBROUTINE set_tides_tile
      END MODULE set_tides_mod
