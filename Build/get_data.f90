      SUBROUTINE get_data (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and other data from      !
!  NetCDF files.  If there is more than one time-record,  data is      !
!  loaded into global  two-time  record arrays. The interpolation      !
!  is carried elsewhere.                                               !
!                                                                      !
!  Currently, this routine is only executed in serial mode by the      !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)
      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, ic, my_tile
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks time dependent data.
!-----------------------------------------------------------------------
!
!  Point Source/Sink vertically integrated mass transport.
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        CALL get_ngfld (ng, iNLM, idRtra, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % QbarG(:,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iNLM, idRtrc(i), SSF(ng)%ncid,            &
     &                    1, SSF(ng), update(1),                        &
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Surface wind stress components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUsms, ncFRCid(idUsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % umask,                                 &
     &                FORCES(ng) % sustrG)
      IF (exit_flag.ne.NoError) RETURN
      CALL get_2dfld (ng, iNLM, idVsms, ncFRCid(idVsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % vmask,                                 &
     &                FORCES(ng) % svstrG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTsur(itemp),                          &
     &                ncFRCid(idTsur(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % stflxG(:,:,:,itemp))
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: sea surface temperature
!  (SST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSTc, ncFRCid(idSSTc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % sstG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: heat flux sensitivity to
!  SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, iddQdT, ncFRCid(iddQdT,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % dqdtG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idsfwf, ncFRCid(idsfwf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % stflxG(:,:,:,isalt))
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSSc, ncFRCid(idSSSc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % sssG)
      IF (exit_flag.ne.NoError) RETURN
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.  In
!  grid refinement, only the coarser grid (RefineScale(ng)=0) open
!  boundary conditions data is processed and needed.
!=======================================================================
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(iwest,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        DO i=1,NT(ng)
          IF (LBC(iwest,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(iwest,i), BRY(ng)%ncid,    &
     &                      1, BRY(ng), update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
        DO i=1,NT(ng)
          IF (LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(ieast,i), BRY(ng)%ncid,    &
     &                      1, BRY(ng), update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
        DO i=1,NT(ng)
          IF (LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(isouth,i), BRY(ng)%ncid,   &
     &                      1, BRY(ng), update(1),                      &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
        DO i=1,NT(ng)
          IF (LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(inorth,i), BRY(ng)%ncid,   &
     &                      1, BRY(ng), update(1),                      &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_north(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
      END IF
!
!=======================================================================
!  Read in data from Climatology NetCDF file.
!=======================================================================
!
!  Free-surface.
!
      IF (LsshCLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idSSHc, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  CLIMA(ng) % sshG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  2D momentum.
!
      IF (Lm2CLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idUbcl, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % umask,                               &
     &                  CLIMA(ng) % ubarclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_2dfld (ng, iNLM, idVbcl, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % vmask,                               &
     &                  CLIMA(ng) % vbarclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  3D momentum.
!
      IF (Lm3CLM(ng)) THEN
        CALL get_3dfld (ng, iNLM, idUclm, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  GRID(ng) % umask,                               &
     &                  CLIMA(ng) % uclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_3dfld (ng, iNLM, idVclm, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  GRID(ng) % vmask,                               &
     &                  CLIMA(ng) % vclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Tracers.
!
      ic=0
      DO i=1,NT(ng)
        IF (LtracerCLM(i,ng)) THEN
          ic=ic+1
          CALL get_3dfld (ng, iNLM, idTclm(i), CLM(ng)%ncid,            &
     &                    1, CLM(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
     &                    GRID(ng) % rmask,                             &
     &                    CLIMA(ng) % tclmG(:,:,:,:,ic))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
      RETURN
      END SUBROUTINE get_data
