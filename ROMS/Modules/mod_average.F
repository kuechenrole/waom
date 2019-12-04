      MODULE mod_average
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  The strategy here is  to define all possible pointers in the        !
!  time-averaged structure and allocate only those requested by        !
!  the user. This will facilitate a better management of memory.       !
!                                                                      !
!  Time-averaged state variables for output purposes.                  !
!                                                                      !
!  avgu2d     2D velocity  component (m/s) in the XI-direction.        !
!  avgv2d     2D velocity  component (m/s) in the ETA-direction.       !
!  avgu2dE    2D Eastward  component (m/s) at RHO-points.              !
!  avgv2dN    2D Northward component (m/s) at RHO-points.              !
!  avgzeta    Free surface (m).                                        !
!  avgssflxSOF    Ice shelf melt rate (m/s).
!  avgshflxSOF    Ice shelf melt rate (m/s).
!  avgismr    Ice shelf melt rate (m/s).
!  avgustar   friction velocity (m/s).
!  avgTstar   thermal driving (degC).
!  avgu3d     3D velocity  component (m/s) in the XI-direction.        !
!  avgv3d     3D velocity  component (m/s) in the ETA-direction.       !
!  avgu3dE    3D Eastward  component (m/s) at RHO-points.              !
!  avgv3dN    3D Northward component (m/s) at RHO-points.              !
!  avgw3d     S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!  avgwvel    3D "true" vertical velocity (m/s).                       !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  avgAKt     Vertical diffusion of temperature (m2/s).                !
!  avgAKv     Vertical viscosity (m2/s).                               !
!  avgAKs     Vertical diffusion of Salinity (m2/s).                   !
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
!                                                                      !
!  Time-averaged surface and bottom fluxes.                            !
!                                                                      !
!  avgsus     Surface u-momentum stress (N/m2).                        !
!  avgsvs     Surface v-momentum stress (N/m2).                        !
!  avgbus     Bottom u-momentum stress (N/m2).                         !
!  avgbvs     Bottom v-momentum stress (N/m2).                         !
!  avgstf     Surface net heat flux (W/m2).                            !
!  avgswf     Surface net freshwater flux (kg/m2/s).                   !
!  avgsrf     Shortwave radiation flux (W/m2).                         !
!                                                                      !
!  Time-averaged quadratic fields.                                     !
!                                                                      !
!  avgZZ      Quadratic term <zeta*zeta> for free-surface.             !
!  avgU2      Quadratic term <ubar*ubar> for 2D momentum at U-points.  !
!  avgV2      Quadratic term <vbar*vbar> for 2D momentum at V-points.  !
!  avgUU      Quadratic term <u*u> for 3D momentum at U-points.        !
!  avgVV      Quadratic term <v*v> for 3D momentum at V-points.        !
!  avgUV      Quadratic term <u*v> for 3D momentum at RHO-points.      !
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                         !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                         !
!  avgTT      Quadratic term <t*t> for tracers.                        !
!  avgUT      Quadratic term <u*t> for potential temperature and       !
!               salinity at U-points.                                  !
!  avgVT      Quadratic term <v*t> for potential temperature and       !
!               salinity at V-points.                                  !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).             !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).             !
!                                                                      !
!  Time-averages vorticity fields.                                     !
!                                                                      !
!  avgpvor2d  2D, vertically integrated, potential vorticity.          !
!  avgrvor2d  2D, vertically integrated, relative vorticity.           !
!  avgpvor3d  3D potential vorticity.                                  !
!  rvorvor2d  3D relative vorticity.                                   !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_AVERAGE
!
!  Time-averaged state variables.
!
          real(r8), pointer :: avgzeta(:,:)
          real(r8), pointer :: avgssflxSOF(:,:)
          real(r8), pointer :: avgshflxSOF(:,:)
          real(r8), pointer :: avgismr(:,:)
          real(r8), pointer :: avgustar(:,:)
          real(r8), pointer :: avgTstar(:,:)
          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)
          real(r8), pointer :: avgu2dE(:,:)
          real(r8), pointer :: avgv2dN(:,:)
          real(r8), pointer :: avgu3d(:,:,:)
          real(r8), pointer :: avgv3d(:,:,:)
          real(r8), pointer :: avgu3dE(:,:,:)
          real(r8), pointer :: avgv3dN(:,:,:)
          real(r8), pointer :: avgw3d(:,:,:)
          real(r8), pointer :: avgwvel(:,:,:)
          real(r8), pointer :: avgrho(:,:,:)
          real(r8), pointer :: avgt(:,:,:,:)
          real(r8), pointer :: avgAKv(:,:,:)
          real(r8), pointer :: avgAKt(:,:,:)
          real(r8), pointer :: avgAKs(:,:,:)
          real(r8), pointer :: avghsbl(:,:)
          real(r8), pointer :: avghbbl(:,:)
!
!  Time-averaged surface and bottom fluxes.
!
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)
          real(r8), pointer :: avgbus(:,:)
          real(r8), pointer :: avgbvs(:,:)
          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
          real(r8), pointer :: avgsrf(:,:)
!
!  Time-averaged quadratic fields.
!
          real(r8), pointer :: avgZZ(:,:)
          real(r8), pointer :: avgU2(:,:)
          real(r8), pointer :: avgV2(:,:)
          real(r8), pointer :: avgUU(:,:,:)
          real(r8), pointer :: avgUV(:,:,:)
          real(r8), pointer :: avgVV(:,:,:)
          real(r8), pointer :: avgHuon(:,:,:)
          real(r8), pointer :: avgHvom(:,:,:)
          real(r8), pointer :: avgTT(:,:,:,:)
          real(r8), pointer :: avgUT(:,:,:,:)
          real(r8), pointer :: avgVT(:,:,:,:)
          real(r8), pointer :: avgHuonT(:,:,:,:)
          real(r8), pointer :: avgHvomT(:,:,:,:)
!
!  Time-averaged vorticity fields.
!
          real(r8), pointer :: avgpvor2d(:,:)
          real(r8), pointer :: avgrvor2d(:,:)
          real(r8), pointer :: avgpvor3d(:,:,:)
          real(r8), pointer :: avgrvor3d(:,:,:)
        END TYPE T_AVERAGE
        TYPE (T_AVERAGE), allocatable :: AVERAGE(:)
      CONTAINS
      SUBROUTINE allocate_average (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
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
      IF (ng.eq.1 ) allocate ( AVERAGE(Ngrids) )
!
!  Time-averaged state variables.
!
      IF (Aout(idFsur,ng)) THEN
        allocate ( AVERAGE(ng) % avgzeta(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idssflxSOF,ng)) THEN
         allocate ( AVERAGE(ng) % avgssflxSOF(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idshflxSOF,ng)) THEN
         allocate ( AVERAGE(ng) % avgshflxSOF(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idismr,ng)) THEN
         allocate ( AVERAGE(ng) % avgismr(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idustar,ng)) THEN
         allocate ( AVERAGE(ng) % avgustar(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idTstar,ng)) THEN
         allocate ( AVERAGE(ng) % avgTstar(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUbar,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbar,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idu2dE,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2dE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idv2dN,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2dN(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idVvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idu3dE,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3dE(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idv3dN,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3dN(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idOvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgw3d(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idWvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgwvel(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idDano,ng)) THEN
        allocate ( AVERAGE(ng) % avgrho(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (Aout(idVvis,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKv(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idTdif,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKt(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idSdif,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKs(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idHsbl,ng)) THEN
        allocate ( AVERAGE(ng) % avghsbl(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idHbbl,ng)) THEN
        allocate ( AVERAGE(ng) % avghbbl(LBi:UBi,LBj:UBj) )
      END IF
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        allocate ( AVERAGE(ng) % avgsus(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVsms,ng)) THEN
        allocate ( AVERAGE(ng) % avgsvs(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUbms,ng)) THEN
        allocate ( AVERAGE(ng) % avgbus(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbms,ng)) THEN
        allocate ( AVERAGE(ng) % avgbvs(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idTsur(itemp),ng)) THEN
        allocate ( AVERAGE(ng) % avgstf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idTsur(isalt),ng)) THEN
        allocate ( AVERAGE(ng) % avgswf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idSrad,ng)) THEN
        allocate ( AVERAGE(ng) % avgsrf(LBi:UBi,LBj:UBj) )
      END IF
!
!  Time-averaged quadratic fields.
!
      IF (Aout(idZZav,ng)) THEN
        allocate ( AVERAGE(ng) % avgU2(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idU2av,ng)) THEN
        allocate ( AVERAGE(ng) % avgV2(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idV2av,ng)) THEN
        allocate ( AVERAGE(ng) % avgZZ(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUUav,ng)) THEN
        allocate ( AVERAGE(ng) % avgUU(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idVVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgVV(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idUVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgUV(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idHUav,ng)) THEN
        allocate ( AVERAGE(ng) % avgHuon(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idHVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgHvom(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (ANY(Aout(idTTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgTT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (ANY(Aout(idUTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgUT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (ANY(Aout(idVTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgVT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (ANY(Aout(iHUTav(:),ng))) THEN
        allocate ( AVERAGE(ng)% avgHuonT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (ANY(Aout(iHVTav(:),ng))) THEN
        allocate ( AVERAGE(ng)% avgHvomT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
!
!  Time-averaged vorticity fields.
!
      IF (Aout(id2dPV,ng)) THEN
        allocate ( AVERAGE(ng) % avgpvor2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(id2dRV,ng)) THEN
        allocate ( AVERAGE(ng) % avgrvor2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(id3dPV,ng)) THEN
        allocate ( AVERAGE(ng) % avgpvor3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(id3dRV,ng)) THEN
        allocate ( AVERAGE(ng) % avgrvor3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      RETURN
      END SUBROUTINE allocate_average
      SUBROUTINE initialize_average (ng, tile)
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
      USE mod_ncparam
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
!  Time-averaged state variables.
!
      IF (Aout(idFsur,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgzeta(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idssflxSOF,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgssflxSOF(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idshflxSOF,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgshflxSOF(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idismr,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgismr(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idustar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgustar(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idTstar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgTstar(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idu2dE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2dE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idv2dN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2dN(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUvel,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVvel,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idu3dE,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3dE(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idv3dN,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3dN(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idOvel,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgw3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWvel,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgwvel(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idDano,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgrho(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgt(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVvis,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKv(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTdif,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKt(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idSdif,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKs(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHsbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghsbl(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idHbbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghbbl(i,j) = IniVal
          END DO
        END DO
      END IF
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsvs(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbvs(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idTsur(itemp),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgstf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idTsur(isalt),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgswf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idSrad,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsrf(i,j) = IniVal
          END DO
        END DO
      END IF
!
!  Time-averaged quadratic fields.
!
      IF (Aout(idZZav,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgU2(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idU2av,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgV2(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idV2av,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgZZ(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUUav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgUU(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgVV(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idUVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgUV(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHUav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgHuon(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgHvom(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idTTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgTT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idUTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgUT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idVTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgVT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(iHUTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgHuonT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(iHVTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgHvomT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
!
!  Time-averaged vorticity fields.
!
      IF (Aout(id2dPV,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgpvor2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(id2dRV,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgrvor2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(id3dPV,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgpvor3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(id3dRV,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgrvor3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_average
      END MODULE mod_average
