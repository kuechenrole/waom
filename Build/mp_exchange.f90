      MODULE mp_exchange_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Distributed-memory tile exchange:                                   !
!                                                                      !
!  This routine updates the I,J tile overlap halo of NV variables.     !
!  It exchanges the specified number of "ghost-points".  In  order     !
!  to minimize the number send and receive calls, the ghost-points     !
!  are included in the buffers.  Therefore, the order of the pack,     !
!  send, receive, and unpack is crucial.                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number.                                  !
!     model       Calling model identifier.                            !
!     Nvar        Number of variables for aggregated exchanges.        !
!     Istr        Starting tile index in the I-direction.              !
!     Iend        Ending   tile index in the I-direction.              !
!     Jstr        Starting tile index in the J-direction.              !
!     Jend        Ending   tile index in the J-direction.              !
!     LBi         I-dimension Lower bound.                             !
!     UBi         I-dimension Upper bound.                             !
!     LBj         J-dimension Lower bound.                             !
!     UBj         J-dimension Upper bound.                             !
!     LBk         K-dimension Lower bound.                             !
!     UBk         K-dimension Upper bound.                             !
!     LBt         T-dimension Lower bound.                             !
!     UBt         T-dimension Upper bound.                             !
!     Nghost      Number of ghost-points in the halo region.           !
!     EW_periodic Switch indicating EW periodicity exchanges.          !
!     NS_periodic Switch indicating NS periodicity exchanges.          !
!     A           2D tiled array to process.                           !
!     B           2D tiled array (optional) to process.                !
!     C           2D tiled array (optional) to process.                !
!     D           2D tiled array (optional) to process.                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A           Updated tiled array.                                 !
!     B           Updated tiled array (optional).                      !
!     C           Updated tiled array (optional).                      !
!     D           Updated tiled array (optional).                      !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!  mp_exchange2d         2D variables tile exchanges                   !
!  mp_exchange2d_bry     2D boundary variables tile exchanges          !
!  mp_exchange3d         3D variables tile exchanges                   !
!  mp_exchange3d_bry     3D boundary variables tile exchanges          !
!  mp_exchange4d         4D variables tile exchanges                   !
!                                                                      !
!  ad_mp_exchange2d      2D variables tile adjoint exchanges           !
!  ad_mp_exchange2d_bry  2D boundary variables tile adjoint exchanges  !
!  ad_mp_exchange3d      3D variables tile adjoint exchanges           !
!  ad_mp_exchange3d_bry  3D boundary variables tile adjoint exchanges  !
!  ad_mp_exchange4d      4D variables tile adjoint exchanges           !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,  &
     &                           GrecvW, GsendW, Wtile, Wexchange,      &
     &                           GrecvE, GsendE, Etile, Eexchange,      &
     &                           GrecvS, GsendS, Stile, Sexchange,      &
     &                           GrecvN, GsendN, Ntile, Nexchange)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic
      integer, intent(in) :: ng, Nghost
      logical, intent(out) :: Wexchange, Eexchange
      logical, intent(out) :: Sexchange, Nexchange
      integer, intent(out) :: GrecvW, GsendW, Wtile
      integer, intent(out) :: GrecvE, GsendE, Etile
      integer, intent(out) :: GrecvS, GsendS, Stile
      integer, intent(out) :: GrecvN, GsendN, Ntile
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: MyRankI, MyRankJ, Null_Value, rank
      integer, dimension(-1:NtileI(ng),-1:NtileJ(ng)) :: table
!
!-----------------------------------------------------------------------
!  Set tile partition table for looking up adjacent processes.
!-----------------------------------------------------------------------
!
!  Notice that a null value is used in places that data transmition is
!  not required.
!
      Null_Value=MPI_PROC_NULL
      DO j=-1,NtileJ(ng)
        DO i=-1,NtileI(ng)
          table(i,j)=Null_Value
        END DO
      END DO
      rank=0
      DO j=0,NtileJ(ng)-1
        DO i=0,NtileI(ng)-1
          table(i,j)=rank
          IF (MyRank.eq.rank) THEN
            MyRankI=i
            MyRankJ=j
          END IF
          rank=rank+1
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Determine the rank of Western and Eastern tiles.  Then, determine
!  the number of ghost-points to send and receive in the West- and
!  East-directions.
!-----------------------------------------------------------------------
!
!  This logic only works for two and three ghost points. The number of
!  ghost-points changes when periodic boundary condition are activated.
!  The periodicity is as follows:
!
!  If two ghost-points:
!
!                      Lm-2  Lm-1  Lm   Lm+1  Lm+2
!                      -2    -1     0    1     2
!
!  If three ghost-points:
!
!                      Lm-2  Lm-1  Lm   Lm+1  Lm+2   Lm+3
!                      -2    -1     0    1     2      3
!
      IF (EW_periodic) THEN
        IF ((table(MyRankI-1,MyRankJ).eq.Null_Value).and.               &
     &      (NtileI(ng).gt.1)) THEN
          Wtile=table(NtileI(ng)-1,MyRankJ)
          Etile=table(MyRankI+1,MyRankJ)
          GsendW=Nghost
          GsendE=Nghost
          IF (NghostPoints.eq.3) THEN
            GrecvW=Nghost
          ELSE
            GrecvW=Nghost+1
          END IF
          GrecvE=Nghost
        ELSE IF ((table(MyRankI+1,MyRankJ).eq.Null_Value).and.          &
     &           (NtileI(ng).gt.1)) THEN
          Wtile=table(MyRankI-1,MyRankJ)
          Etile=table(0,MyRankJ)
          GsendW=Nghost
          IF (NghostPoints.eq.3) THEN
            GsendE=Nghost
          ELSE
            GsendE=Nghost+1
          END IF
          GrecvW=Nghost
          GrecvE=Nghost
        ELSE
          Wtile=table(MyRankI-1,MyRankJ)
          Etile=table(MyRankI+1,MyRankJ)
          GsendW=Nghost
          GsendE=Nghost
          GrecvW=Nghost
          GrecvE=Nghost
        END IF
      ELSE
        Wtile=table(MyRankI-1,MyRankJ)
        Etile=table(MyRankI+1,MyRankJ)
        GsendW=Nghost
        GsendE=Nghost
        GrecvW=Nghost
        GrecvE=Nghost
      END IF
!
!  Determine exchange switches.
!
      IF (Wtile.eq.Null_Value) THEN
        Wexchange=.FALSE.
      ELSE
        Wexchange=.TRUE.
      END IF
      IF (Etile.eq.Null_Value) THEN
        Eexchange=.FALSE.
      ELSE
        Eexchange=.TRUE.
      END IF
!
!-----------------------------------------------------------------------
!  Determine the rank of Southern and Northern tiles.  Then, determine
!  the number of ghost-points to send and receive in the South- and
!  North-directions.
!-----------------------------------------------------------------------
!
!  This logic only works for two and three ghost-points. The number of
!  ghost-points changes when periodic boundary condition are activated.
!  The periodicity is as follows:
!
!  If two ghost-points:
!
!                      Mm-2  Mm-1  Mm   Mm+1  Mm+2
!                      -2    -1     0    1     2
!
!  If three ghost-points:
!
!                      Mm-2  Mm-1  Mm   Mm+1  Mm+2  Mm+3
!                      -2    -1     0    1     2     3
!
      IF (NS_periodic) THEN
        IF ((table(MyRankI,MyRankJ-1).eq.Null_Value).and.               &
     &      (NtileJ(ng).gt.1)) THEN
          Stile=table(MyRankI,NtileJ(ng)-1)
          Ntile=table(MyRankI,MyRankJ+1)
          GsendS=Nghost
          GsendN=Nghost
          IF (NghostPoints.eq.3) THEN
            GrecvS=Nghost
          ELSE
            GrecvS=Nghost+1
          END IF
          GrecvN=Nghost
        ELSE IF ((table(MyRankI,MyRankJ+1).eq.Null_Value).and.          &
     &           (NtileJ(ng).gt.1)) then
          Stile=table(MyRankI,MyRankJ-1)
          Ntile=table(MyRankI,0)
          GsendS=Nghost
          IF (NghostPoints.eq.3) THEN
            GsendN=Nghost
          ELSE
            GsendN=Nghost+1
          END IF
          GrecvS=Nghost
          GrecvN=Nghost
        ELSE
          Stile=table(MyRankI,MyRankJ-1)
          Ntile=table(MyRankI,MyRankJ+1)
          GsendS=Nghost
          GsendN=Nghost
          GrecvS=Nghost
          GrecvN=Nghost
        END IF
      ELSE
        Stile=table(MyRankI,MyRankJ-1)
        Ntile=table(MyRankI,MyRankJ+1)
        GsendS=Nghost
        GsendN=Nghost
        GrecvS=Nghost
        GrecvN=Nghost
      END IF
!
!  Determine exchange switches.
!
      IF (Stile.eq.Null_Value) THEN
        Sexchange=.FALSE.
      ELSE
        Sexchange=.TRUE.
      END IF
      IF (Ntile.eq.Null_Value) THEN
        Nexchange=.FALSE.
      ELSE
        Nexchange=.TRUE.
      END IF
      RETURN
      END SUBROUTINE tile_neighbors
!
!***********************************************************************
      SUBROUTINE mp_exchange2d (ng, tile, model, Nvar,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Nghost, EW_periodic, NS_periodic,       &
     &                          A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic
      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
      real(r8), intent(inout), optional :: B(LBi:,LBj:)
      real(r8), intent(inout), optional :: C(LBi:,LBj:)
      real(r8), intent(inout), optional :: D(LBi:,LBj:)
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange
      integer :: i, icS, icN, ioff, Imin, Imax, Ilen
      integer :: j, jcW, jcE, joff, Jmin, Jmax, Jlen
      integer :: m, mc, Ierror, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN
      integer, dimension(MPI_STATUS_SIZE,4) :: status
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeJ(ng)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloSizeI(ng)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloSizeI(ng)) :: recvS, recvN
      character (len=MPI_MAX_ERROR_STRING) :: string
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
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 40)
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*Jlen
      NSsize=Nvar*(Nghost+pp)*Ilen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE2D - communication buffer too small, ',  &
     &          a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*Jlen
          i=Istr+m-1
          DO j=Jmin,Jmax
            sizeW=sizeW+1
            jcW=1+(j-Jmin)+mc
            sendW(jcW)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
              sendW(jcW)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
              sendW(jcW)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
              sendW(jcW)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*Jlen
          i=Iend-GsendE+m
          DO j=Jmin,Jmax
            sizeE=sizeE+1
            jcE=1+(j-Jmin)+mc
            sendE(jcE)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
              sendE(jcE)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
              sendE(jcE)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
              sendE(jcE)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE2D - error during ',a,                 &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*Jlen
          i=Istr-m
          DO j=Jmin,Jmax
            jcW=1+(j-Jmin)+mc
            A(i,j)=recvW(jcW)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
              B(i,j)=recvW(jcW)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
              C(i,j)=recvW(jcW)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
              D(i,j)=recvW(jcW)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvE
          mc=(m-1)*Jlen
          i=Iend+m
          DO j=Jmin,Jmax
            jcE=1+(j-Jmin)+mc
            A(i,j)=recvE(jcE)
          ENDDO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
              B(i,j)=recvE(jcE)
            ENDDO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
              C(i,j)=recvE(jcE)
            ENDDO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
              D(i,j)=recvE(jcE)
            ENDDO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*Ilen
          j=Jstr+m-1
          DO i=Imin,Imax
            sizeS=sizeS+1
            icS=1+(i-Imin)+mc
            sendS(icS)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
              sendS(icS)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
              sendS(icS)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
              sendS(icS)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*Ilen
          j=Jend-GsendN+m
          DO i=Imin,Imax
            sizeN=sizeN+1
            icN=1+(i-Imin)+mc
            sendN(icN)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
              sendN(icN)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
              sendN(icN)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
              sendN(icN)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*Ilen
          j=Jstr-m
          DO i=Imin,Imax
            icS=1+(i-Imin)+mc
            A(i,j)=recvS(icS)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
              B(i,j)=recvS(icS)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
              C(i,j)=recvS(icS)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
              D(i,j)=recvS(icS)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvN
          mc=(m-1)*Ilen
          j=Jend+m
          DO i=Imin,Imax
            icN=1+(i-Imin)+mc
            A(i,j)=recvN(icN)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
              B(i,j)=recvN(icN)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
              C(i,j)=recvN(icN)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
              D(i,j)=recvN(icN)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 40)
      RETURN
      END SUBROUTINE mp_exchange2d
!
!***********************************************************************
      SUBROUTINE mp_exchange2d_bry (ng, tile, model, Nvar, boundary,    &
     &                              LBij, UBij,                         &
     &                              Nghost, EW_periodic, NS_periodic,   &
     &                              A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic
      integer, intent(in) :: ng, tile, model, Nvar, boundary
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: Nghost
!
      real(r8), intent(inout) :: A(LBij:)
      real(r8), intent(inout), optional :: B(LBij:)
      real(r8), intent(inout), optional :: C(LBij:)
      real(r8), intent(inout), optional :: D(LBij:)
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange
      integer :: i, icS, icN
      integer :: j, jcW, jcE
      integer :: m, Ierror, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN
      integer, dimension(MPI_STATUS_SIZE,4) :: status
!
      real(r8), dimension(Nvar*HaloBry(ng)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloBry(ng)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloBry(ng)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloBry(ng)) :: recvS, recvN
      character (len=MPI_MAX_ERROR_STRING) :: string
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
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 43)
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Adjust exchange swiches according to boundary edge to process.
!
      Wexchange=Wexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Eexchange=Eexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Sexchange=Sexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
      Nexchange=Nexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)
      NSsize=Nvar*(Nghost+pp)
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE2D_BRY - communication buffer too ',     &
     &          'small, ',a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        jcW=0
        sizeW=0
        DO m=1,GsendW
          i=Istr+m-1
          sizeW=sizeW+1
          jcW=jcW+1
          sendW(jcW)=A(i)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            sizeW=sizeW+1
            jcW=jcW+1
            sendW(jcW)=B(i)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            sizeW=sizeW+1
            jcW=jcW+1
            sendW(jcW)=C(i)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            sizeW=sizeW+1
            jcW=jcW+1
            sendW(jcW)=D(i)
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        jcE=0
        sizeE=0
        DO m=1,GsendE
          i=Iend-GsendE+m
          sizeE=sizeE+1
          jcE=jcE+1
          sendE(jcE)=A(i)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            sizeE=sizeE+1
            jcE=jcE+1
            sendE(jcE)=B(i)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            sizeE=sizeE+1
            jcE=jcE+1
            sendE(jcE)=C(i)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            sizeE=sizeE+1
            jcE=jcE+1
            sendE(jcE)=D(i)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE2D_BRY - error during ',a,             &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
        jcW=0
        DO m=GrecvW,1,-1
          i=Istr-m
          jcW=jcW+1
          A(i)=recvW(jcW)
        END DO
        IF (PRESENT(B)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            jcW=jcW+1
            B(i)=recvW(jcW)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            jcW=jcW+1
            C(i)=recvW(jcW)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            jcW=jcW+1
            D(i)=recvW(jcW)
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        jcE=0
        DO m=1,GrecvE
          i=Iend+m
          jcE=jcE+1
          A(i)=recvE(jcE)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GrecvE
            i=Iend+m
            jcE=jcE+1
            B(i)=recvE(jcE)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GrecvE
            i=Iend+m
            jcE=jcE+1
            C(i)=recvE(jcE)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GrecvE
            i=Iend+m
            jcE=jcE+1
            D(i)=recvE(jcE)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        icS=0
        sizeS=0
        DO m=1,GsendS
          j=Jstr+m-1
          sizeS=sizeS+1
          icS=icS+1
          sendS(icS)=A(j)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            sizeS=sizeS+1
            icS=icS+1
            sendS(icS)=B(j)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            sizeS=sizeS+1
            icS=icS+1
            sendS(icS)=C(j)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            sizeS=sizeS+1
            icS=icS+1
            sendS(icS)=D(j)
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        icN=0
        sizeN=0
        DO m=1,GsendN
          j=Jend-GsendN+m
          sizeN=sizeN+1
          icN=icN+1
          sendN(icN)=A(j)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            sizeN=sizeN+1
            icN=icN+1
            sendN(icN)=B(j)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            sizeN=sizeN+1
            icN=icN+1
            sendN(icN)=C(j)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            sizeN=sizeN+1
            icN=icN+1
            sendN(icN)=D(j)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        icS=0
        DO m=GrecvS,1,-1
          j=Jstr-m
          icS=icS+1
          A(j)=recvS(icS)
        END DO
        IF (PRESENT(B)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            icS=icS+1
            B(j)=recvS(icS)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            icS=icS+1
            C(j)=recvS(icS)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            icS=icS+1
            D(j)=recvS(icS)
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        icN=0
        DO m=1,GrecvN
          j=Jend+m
          icN=icN+1
          A(j)=recvN(icN)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GrecvN
            j=Jend+m
            icN=icN+1
            B(j)=recvN(icN)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GrecvN
            j=Jend+m
            icN=icN+1
            C(j)=recvN(icN)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GrecvN
            j=Jend+m
            icN=icN+1
            D(j)=recvN(icN)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 43)
      RETURN
      END SUBROUTINE mp_exchange2d_bry
!
!***********************************************************************
      SUBROUTINE mp_exchange3d (ng, tile, model, Nvar,                  &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Nghost, EW_periodic, NS_periodic,       &
     &                          A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic
      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
      real(r8), intent(inout), optional :: B(LBi:,LBj:,LBk:)
      real(r8), intent(inout), optional :: C(LBi:,LBj:,LBk:)
      real(r8), intent(inout), optional :: D(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange
      integer :: i, ikS, ikN, ioff, Imin, Imax, Ilen, IKlen
      integer :: j, jkW, jkE, joff, Jmin, Jmax, Jlen, JKlen
      integer :: k, kc,  m, mc, Ierror, Klen, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN
      integer, dimension(MPI_STATUS_SIZE,4) :: status
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                          (UBk-LBk+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                          (UBk-LBk+1)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                          (UBk-LBk+1)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                          (UBk-LBk+1)) :: recvS, recvN
      character (len=MPI_MAX_ERROR_STRING) :: string
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
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 41)
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IKlen=Ilen*Klen
      JKlen=Jlen*Klen
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*JKlen
      NSsize=Nvar*(Nghost+pp)*IKlen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE3D - communication buffer too small, ',  &
     &          a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*JKlen
          i=Istr+m-1
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jkW=1+(j-Jmin)+kc
              sendW(jkW)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
                sendW(jkW)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
                sendW(jkW)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
                sendW(jkW)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*JKlen
          i=Iend-GsendE+m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jkE=1+(j-Jmin)+kc
              sendE(jkE)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
                sendE(jkE)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
                sendE(jkE)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
                sendE(jkE)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*JKlen
          i=Istr-m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              jkW=1+(j-Jmin)+kc
              A(i,j,k)=recvW(jkW)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
                B(i,j,k)=recvW(jkW)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
                C(i,j,k)=recvW(jkW)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
                D(i,j,k)=recvW(jkW)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE3D - error during ',a,                 &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvE
          mc=(m-1)*JKlen
          i=Iend+m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              jkE=1+(j-Jmin)+kc
              A(i,j,k)=recvE(jkE)
            END DO
          ENDDO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
                B(i,j,k)=recvE(jkE)
              END DO
            ENDDO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
                C(i,j,k)=recvE(jkE)
              END DO
            ENDDO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
                D(i,j,k)=recvE(jkE)
              END DO
            ENDDO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*IKlen
          j=Jstr+m-1
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              sizeS=sizeS+1
              ikS=1+(i-Imin)+kc
              sendS(ikS)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
                sendS(ikS)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
                sendS(ikS)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
                sendS(ikS)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*IKlen
          j=Jend-GsendN+m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              sizeN=sizeN+1
              ikN=1+(i-Imin)+kc
              sendN(ikN)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
                sendN(ikN)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
                sendN(ikN)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
                sendN(ikN)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*IKlen
          j=Jstr-m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              ikS=1+(i-Imin)+kc
              A(i,j,k)=recvS(ikS)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
                B(i,j,k)=recvS(ikS)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
                C(i,j,k)=recvS(ikS)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
                D(i,j,k)=recvS(ikS)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvN
          mc=(m-1)*IKlen
          j=Jend+m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              ikN=1+(i-Imin)+kc
              A(i,j,k)=recvN(ikN)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
                B(i,j,k)=recvN(ikN)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
                C(i,j,k)=recvN(ikN)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
                D(i,j,k)=recvN(ikN)
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 41)
      RETURN
      END SUBROUTINE mp_exchange3d
!
!***********************************************************************
      SUBROUTINE mp_exchange3d_bry (ng, tile, model, Nvar, boundary,    &
     &                              LBij, UBij, LBk, UBk,               &
     &                              Nghost, EW_periodic, NS_periodic,   &
     &                              A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic
      integer, intent(in) :: ng, tile, model, Nvar, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk
      integer, intent(in) :: Nghost
!
      real(r8), intent(inout) :: A(LBij:,LBk:)
      real(r8), intent(inout), optional :: B(LBij:,LBk:)
      real(r8), intent(inout), optional :: C(LBij:,LBk:)
      real(r8), intent(inout), optional :: D(LBij:,LBk:)
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange
      integer :: i, ikS, ikN, ioff
      integer :: j, jkW, jkE, joff
      integer :: k, m, mc, Ierror, Klen, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN
      integer, dimension(MPI_STATUS_SIZE,4) :: status
!
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: recvS, recvN
      character (len=MPI_MAX_ERROR_STRING) :: string
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
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 43)
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Adjust exchange swiches according to boundary edge to process.
!
      Wexchange=Wexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Eexchange=Eexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Sexchange=Sexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
      Nexchange=Nexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Klen=UBk-LBk+1
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*Klen
      NSsize=Nvar*(Nghost+pp)*Klen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE3D_BRY - communication buffer too ',     &
     &          'small, ', a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*Klen
          i=Istr+m-1
          DO k=LBk,UBk
            sizeW=sizeW+1
            jkW=1+(k-LBk)+mc
            sendW(jkW)=A(i,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
              sendW(jkW)=B(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
              sendW(jkW)=C(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
              sendW(jkW)=D(i,k)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*Klen
          i=Iend-GsendE+m
          DO k=LBk,UBk
            sizeE=sizeE+1
            jkE=1+(k-LBk)+mc
            sendE(jkE)=A(i,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
              sendE(jkE)=B(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
              sendE(jkE)=C(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
              sendE(jkE)=D(i,k)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*Klen
          i=Istr-m
          DO k=LBk,UBk
            jkW=1+(k-LBk)+mc
            A(i,k)=recvW(jkW)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
              B(i,k)=recvW(jkW)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
              C(i,k)=recvW(jkW)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
              D(i,k)=recvW(jkW)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE3D_BRY - error during ',a,             &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvE
          mc=(m-1)*Klen
          i=Iend+m
          DO k=LBk,UBk
            jkE=1+(k-LBk)+mc
            A(i,k)=recvE(jkE)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
              B(i,k)=recvE(jkE)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
              C(i,k)=recvE(jkE)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
              D(i,k)=recvE(jkE)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*Klen
          j=Jstr+m-1
          DO k=LBk,UBk
            sizeS=sizeS+1
            ikS=1+(k-LBk)+mc
            sendS(ikS)=A(j,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
              sendS(ikS)=B(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
              sendS(ikS)=C(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
              sendS(ikS)=D(j,k)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*Klen
          j=Jend-GsendN+m
          DO k=LBk,UBk
            sizeN=sizeN+1
            ikN=1+(k-LBk)+mc
            sendN(ikN)=A(j,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
              sendN(ikN)=B(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
              sendN(ikN)=C(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
              sendN(ikN)=D(j,k)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*Klen
          j=Jstr-m
          DO k=LBk,UBk
            ikS=1+(k-LBk)+mc
            A(j,k)=recvS(ikS)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
              B(j,k)=recvS(ikS)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
              C(j,k)=recvS(ikS)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
              D(j,k)=recvS(ikS)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvN
          mc=(m-1)*Klen
          j=Jend+m
          DO k=LBk,UBk
            ikN=1+(k-LBk)+mc
            A(j,k)=recvN(ikN)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
              B(j,k)=recvN(ikN)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
              C(j,k)=recvN(ikN)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
              D(j,k)=recvN(ikN)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 43)
      RETURN
      END SUBROUTINE mp_exchange3d_bry
!
!***********************************************************************
      SUBROUTINE mp_exchange4d (ng, tile, model, Nvar,                  &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt, &
     &                          Nghost, EW_periodic, NS_periodic,       &
     &                          A, B, C)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic
      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt
      integer, intent(in) :: Nghost
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBt:)
      real(r8), intent(inout), optional :: B(LBi:,LBj:,LBk:,LBt:)
      real(r8), intent(inout), optional :: C(LBi:,LBj:,LBk:,LBt:)
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange
      integer :: i, ikS, ikN, ioff, Imin, Imax, Ilen, IKlen, IKTlen
      integer :: j, jkW, jkE, joff, Jmin, Jmax, Jlen, JKlen, JKTlen
      integer :: k, kc, m, mc, Ierror, Klen, Lstr, Tlen, pp
      integer :: l, lc
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN
      integer, dimension(MPI_STATUS_SIZE,4) :: status
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: recvS, recvN
      character (len=MPI_MAX_ERROR_STRING) :: string
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
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 42)
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      Tlen=UBt-LBt+1
      IKlen=Ilen*Klen
      JKlen=Jlen*Klen
      IKTlen=IKlen*Tlen
      JKTlen=JKlen*Tlen
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*JKTlen
      NSsize=Nvar*(Nghost+pp)*IKTlen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE4D - communication buffer too small, ',  &
     &          a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*JKTlen
          i=Istr+m-1
          DO l=LBt,UBt
           lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=1+(j-Jmin)+kc
                sendW(jkW)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKTlen
            i=Istr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeW=sizeW+1
                  jkW=joff+1+(j-Jmin)+kc
                  sendW(jkW)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKTlen
            i=Istr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeW=sizeW+1
                  jkW=joff+1+(j-Jmin)+kc
                  sendW(jkW)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*JKTlen
          i=Iend-GsendE+m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=1+(j-Jmin)+kc
                sendE(jkE)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKTlen
            i=Iend-GsendE+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeE=sizeE+1
                  jkE=joff+1+(j-Jmin)+kc
                  sendE(jkE)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKTlen
            i=Iend-GsendE+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeE=sizeE+1
                  jkE=joff+1+(j-Jmin)+kc
                  sendE(jkE)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE4D - error during ',a,                 &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*JKTlen
          i=Istr-m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                jkW=1+(j-Jmin)+kc
                A(i,j,k,l)=recvW(jkW)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKTlen
            i=Istr-m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkW=joff+1+(j-Jmin)+kc
                  B(i,j,k,l)=recvW(jkW)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKTlen
            i=Istr-m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkW=joff+1+(j-Jmin)+kc
                  C(i,j,k,l)=recvW(jkW)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvE
          mc=(m-1)*JKTlen
          i=Iend+m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                jkE=1+(j-Jmin)+kc
                A(i,j,k,l)=recvE(jkE)
              END DO
            END DO
          ENDDO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKTlen
            i=Iend+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkE=joff+1+(j-Jmin)+kc
                  B(i,j,k,l)=recvE(jkE)
                END DO
              END DO
            ENDDO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKTlen
            i=Iend+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkE=joff+1+(j-Jmin)+kc
                  C(i,j,k,l)=recvE(jkE)
                END DO
              END DO
            ENDDO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*IKTlen
          j=Jstr+m-1
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=1+(i-Imin)+kc
                sendS(ikS)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKTlen
            j=Jstr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeS=sizeS+1
                  ikS=ioff+1+(i-Imin)+kc
                  sendS(ikS)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKTlen
            j=Jstr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeS=sizeS+1
                  ikS=ioff+1+(i-Imin)+kc
                  sendS(ikS)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*IKTlen
          j=Jend-GsendN+m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=1+(i-Imin)+kc
                sendN(ikN)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKTlen
            j=Jend-GsendN+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeN=sizeN+1
                  ikN=ioff+1+(i-Imin)+kc
                  sendN(ikN)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKTlen
            j=Jend-GsendN+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeN=sizeN+1
                  ikN=ioff+1+(i-Imin)+kc
                  sendN(ikN)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*IKTlen
          j=Jstr-m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                ikS=1+(i-Imin)+kc
                A(i,j,k,l)=recvS(ikS)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKTlen
            j=Jstr-m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikS=ioff+1+(i-Imin)+kc
                  B(i,j,k,l)=recvS(ikS)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKTlen
            j=Jstr-m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikS=ioff+1+(i-Imin)+kc
                  C(i,j,k,l)=recvS(ikS)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
        DO m=1,GrecvN
          mc=(m-1)*IKTlen
          j=Jend+m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                ikN=1+(i-Imin)+kc
                A(i,j,k,l)=recvN(ikN)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKTlen
            j=Jend+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikN=ioff+1+(i-Imin)+kc
                  B(i,j,k,l)=recvN(ikN)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKTlen
            j=Jend+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikN=ioff+1+(i-Imin)+kc
                  C(i,j,k,l)=recvN(ikN)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 42)
      RETURN
      END SUBROUTINE mp_exchange4d
      END MODULE mp_exchange_mod
