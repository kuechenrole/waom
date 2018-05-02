      MODULE mod_clima
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sea surface height fields.                                          !
!                                                                      !
!   ssh         Climatology for sea surface height (m).                !
!   sshG        Latest two-time snapshots of input "ssh" grided        !
!                 data used for interpolation.                         !
!   zeta_ads    Sensitivity functional for sea surface height.         !
!   zeta_adsF   Latest two-time snapshots of input "zeta_ads" grided   !
!                 data used fot interpolation.                         !
!                                                                      !
!  2D momentum fields.                                                 !
!                                                                      !
!   ubarclm     Vertically integrated U-momentum climatology (m/s).    !
!   ubarclmG    Latest two-time snapshots of input "ubarclm" grided    !
!                 data used for interpolation.                         !
!   ubar_ads    Sensitivity functional for vertically integrated       !
!                 U-momentum.                                          !
!   ubar_adsG   Latest two-time snapshots of input "ubar_ads" grided   !
!                 data used for interpolation.                         !
!   vbarclm     Vertically integrated V-momentum climatology (m/s).    !
!   vbarclmG    Latest two-time snapshots of input "vbarclm" grided    !
!                 data used for interpolation.                         !
!   vbar_ads    Sensitivity functional for vertically integrated       !
!                 V-momentum.                                          !
!   vbar_adsG   Latest two-time snapshots of input "vbar_ads" grided   !
!                 data used for interpolation.                         !
!                                                                      !
!  Tracer fields.                                                      !
!                                                                      !
!   tclm        Climatology for tracer type variables (usually,        !
!                 temperature: degC; salinity: PSU).                   !
!   tclmG       Latest two-time snapshots of input "tclm" grided       !
!                 data used for interpolation.                         !
!   t_ads       Sensitivity functional for tracer type variables.      !
!   t_adsG      Latest two-time snapshots of input "t_ads" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  3D momentum climatology.                                            !
!                                                                      !
!   uclm        3D U-momentum climatology (m/s).                       !
!   uclmG       Latest two-time snapshots of input "uclm" grided       !
!                 data used for interpolation.                         !
!   u_ads       Sensitivity functional for 3D U-momentum.              !
!   u_adsG      Latest two-time snapshots of input "u_ads" grided      !
!                 data used for interpolation.                         !
!   vclm        3D V-momentum climatology (m/s).                       !
!   vclmG       Latest two-time snapshots of input "vclm" grided       !
!                 data used for interpolation.                         !
!   v_ads       Sensitivity functional for 3D V-momentum.              !
!   v_adsG      Latest two-time snapshots of input "v_ads" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Nudging variables.                                                  !
!                                                                      !
!   M2nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 2D momentum data.                                    !
!   M3nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 3D momentum data.                                    !
!   Tnudgcof    Time-scale (1/sec) coefficients for nudging towards    !
!                 tracer data.                                         !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_CLIMA
!
!  Climatology/Nudging arrays.
!
          real(r8), pointer :: ssh(:,:)
          real(r8), pointer :: sshG(:,:,:)
          real(r8), pointer :: ubarclm(:,:)
          real(r8), pointer :: vbarclm(:,:)
          real(r8), pointer :: ubarclmG(:,:,:)
          real(r8), pointer :: vbarclmG(:,:,:)
          real(r8), pointer :: uclm(:,:,:)
          real(r8), pointer :: vclm(:,:,:)
          real(r8), pointer :: uclmG(:,:,:,:)
          real(r8), pointer :: vclmG(:,:,:,:)
          real(r8), pointer :: tclm(:,:,:,:)
          real(r8), pointer :: tclmG(:,:,:,:,:)
!
!  Nudging coefficient arrays.
!
          real(r8), pointer :: M2nudgcof(:,:)
          real(r8), pointer :: M3nudgcof(:,:,:)
          real(r8), pointer :: Tnudgcof(:,:,:,:)
        END TYPE T_CLIMA
        TYPE (T_CLIMA), allocatable :: CLIMA(:)
      CONTAINS
      SUBROUTINE allocate_clima (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( CLIMA(Ngrids) )
!
!  Climatology/Nudging arrays.
!
      IF (LsshCLM(ng)) THEN
        allocate ( CLIMA(ng) % ssh(LBi:UBi,LBj:UBj) )
        allocate ( CLIMA(ng) % sshG(LBi:UBi,LBj:UBj,2) )
      END IF
!
      IF (Lm2CLM(ng)) THEN
        allocate ( CLIMA(ng) % ubarclm(LBi:UBi,LBj:UBj) )
        allocate ( CLIMA(ng) % vbarclm(LBi:UBi,LBj:UBj) )
        allocate ( CLIMA(ng) % ubarclmG(LBi:UBi,LBj:UBj,2) )
        allocate ( CLIMA(ng) % vbarclmG(LBi:UBi,LBj:UBj,2) )
      END IF
!
      IF (Lm3CLM(ng)) THEN
        allocate ( CLIMA(ng) % uclm(LBi:UBi,LBj:UBj,N(ng)) )
        allocate ( CLIMA(ng) % vclm(LBi:UBi,LBj:UBj,N(ng)) )
        allocate ( CLIMA(ng) % uclmG(LBi:UBi,LBj:UBj,N(ng),2) )
        allocate ( CLIMA(ng) % vclmG(LBi:UBi,LBj:UBj,N(ng),2) )
      END IF
!
      IF (ANY(LtracerCLM(:,ng)).or.ANY(LnudgeTCLM(:,ng))) THEN
        allocate ( CLIMA(ng) % tclm(LBi:UBi,LBj:UBj,N(ng),NTCLM(ng)) )
        allocate ( CLIMA(ng) % tclmG(LBi:UBi,LBj:UBj,N(ng),2,           &

     &                               NTCLM(ng)) )
      END IF
!
!  Nudging coefficient arrays.
!
      IF (LnudgeM2CLM(ng)) THEN
        allocate ( CLIMA(ng) % M2nudgcof(LBi:UBi,LBj:UBj) )
      END IF
      IF (LnudgeM3CLM(ng)) THEN
        allocate ( CLIMA(ng) % M3nudgcof(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        allocate ( CLIMA(ng) % Tnudgcof(LBi:UBi,LBj:UBj,N(ng),NTCLM(ng)) )
      END IF
      RETURN
      END SUBROUTINE allocate_clima
      SUBROUTINE initialize_clima (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      integer :: itrc, k
      real(r8), parameter :: IniVal = 0.0_r8
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
!  Set array initialization range.
!
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Climatology/Nudging arrays.
!
      IF (LsshCLM(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            CLIMA(ng) % ssh(i,j) = IniVal
            CLIMA(ng) % sshG(i,j,1) = IniVal
            CLIMA(ng) % sshG(i,j,2) = IniVal
          END DO
        END DO
      END IF
!
      IF (Lm2CLM(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            CLIMA(ng) % ubarclm(i,j) = IniVal
            CLIMA(ng) % vbarclm(i,j) = IniVal
            CLIMA(ng) % ubarclmG(i,j,1) = IniVal
            CLIMA(ng) % ubarclmG(i,j,2) = IniVal
            CLIMA(ng) % vbarclmG(i,j,1) = IniVal
            CLIMA(ng) % vbarclmG(i,j,2) = IniVal
          END DO
        END DO
      END IF
!
      IF (Lm3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              CLIMA(ng) % uclm(i,j,k) = IniVal
              CLIMA(ng) % vclm(i,j,k) = IniVal
              CLIMA(ng) % uclmG(i,j,k,1) = IniVal
              CLIMA(ng) % uclmG(i,j,k,2) = IniVal
              CLIMA(ng) % vclmG(i,j,k,1) = IniVal
              CLIMA(ng) % vclmG(i,j,k,2) = IniVal
            END DO
          END DO
        END DO
      END IF
!
      IF (ANY(LtracerCLM(:,ng)).or.ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                CLIMA(ng) % tclm(i,j,k,itrc) = IniVal
                CLIMA(ng) % tclmG(i,j,k,1,itrc) = IniVal
                CLIMA(ng) % tclmG(i,j,k,2,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
!
!  Nudging coefficient arrays.
!
      IF (LnudgeM2CLM(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            CLIMA(ng) % M2nudgcof(i,j) = IniVal
          END DO
        END DO
      END IF
!
      IF (LnudgeM3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              CLIMA(ng) % M3nudgcof(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                CLIMA(ng) % Tnudgcof(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_clima
      END MODULE mod_clima
