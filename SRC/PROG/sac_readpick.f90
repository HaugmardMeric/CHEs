! programme principal sac_readpick
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! ------- readpick : lit les fichiers sac avec les pick.sh       --------   .
! ------- et ecrit un fichier (phases list) lisible par le programme CHE    .
! -------                                                        --------   .
! ------- septembre 2014                                         --------   .
! -------                                                        --------   .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr           --------   .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
!  This program is distributed for research purposes and in the hope        !
!  that it will be useful however without any warranties.                   !
!  Use it on your own risk. Please, report bugs/improvements/... .          !
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
!                                                                           !
!
program readpick
  ! --------------------------------------------------------------------    .
  use modparam
  use typetemps
  use time
  use sac_i_o
  ! -------                                                      --------   .
  implicit none
  ! ---------------------------------------------------------------------   . pour SAC
  real(KIND=4), dimension(:), allocatable :: sacfile
  real(KIND=4) :: delta1, b1, e1
  integer(KIND=4) :: NN, npts1
  character(LEN=112) :: file1
  ! ---------------------------------------------------------------------   . autres
  integer(KIND=wi) :: JD
  type(dataone) :: Adata
  ! --------------------------------------------------------------------    . lecture du fichier sac
  NN = IARGC()
  if ((NN < 1).or.(NN > 1)) then
    write(*,'(a)') "usage:  ./readpick.exe sacfile"
    write(*,'(a)') "        sacfile - input sac file"
    stop
  endif
  call GETARG(1, file1)
  ! --------------------------------------------------------------------    .
  call rbsac(file1,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,        &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
    nzsec,nzmsec,nvhdr,norid,nevid,npts1,internal4,nwfid,nxsize,nysize,unused8,  &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
    kinst,sacfile)
  ! --------------------------------------------------------------------    .
  if (nvhdr /= 6) then
    write(*,*) "ERROR - File: '", TRIM(adjustl(file1)), "' appears to be of non-native &
    &byte-order or is not a SAC file."
    stop
  endif
  ! --------------------------------------------------------------------    .
  Adata%sta%staname=(TRIM(adjustl(kstnm))//'000')                           ! nom de la station
  ! --------------------------------------------------------------------    . initialise
  ! si phase est près du début :  10 sec
  if (T1.le.(b1+10.0_4)) T1=-12345.0_4
  if (T2.le.(b1+10.0_4)) T2=-12345.0_4
  if (T3.le.(b1+10.0_4)) T3=-12345.0_4
  if (T4.le.(b1+10.0_4)) T4=-12345.0_4
  if (T5.le.(b1+10.0_4)) T5=-12345.0_4
  if (T6.le.(b1+10.0_4)) T6=-12345.0_4
  if (T7.le.(b1+10.0_4)) T7=-12345.0_4
  if (T8.le.(b1+10.0_4)) T8=-12345.0_4
  T0=-12345.0_4
  a=-12345.0_4
  o=-12345.0_4
  T9=-12345.0_4
  ! --------------------------------------------------------------------    .
  ! ONDES DIRECTES                                                          .
  ! --------------------------------------------------------------------    .
  ! T1 : temps arrivée de l'onde Pg
  ! T2 : incertitudes absolue sur l'onde Pg                                 ! -> calcul du delta après
  ! --------------------------------------------------------------------    .
  ! T3 : temps arrivée de l'onde Sg
  ! T4 : incertitudes absolue sur l'onde Sg
  ! ---------------------------------------------------------------------   .
  if ((real(T1,wr).ne.-12345.0_wr).and.(real(T2,wr).ne.-12345.0_wr)) then
    Adata%typeonde='G'
    Adata%coefP=0 ! modification a posteriori si necessaire
    Adata%coefS=0
    ! -------                                                  --------    . date
    Adata%tpsR%date%year=int(nzyear,wi)-2000
    ! -------                                                  --------    . jour julien
    Adata%tpsR%date%month=1
    Adata%tpsR%date%day=1
    call JDATE (JD,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
    call GDATE (JD+int(nzjday,wi)-1, Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
    ! -------                                                  --------    .
    Adata%tpsR%date%hour=nzhour
    Adata%tpsR%date%min=nzmin
    Adata%tpsR%secP=real(nzsec,wr)+ real(T1,wr) + real(nzmsec,wr)/1000.0_wr
    ! -------                                                  --------    .
    Adata%sigP=max(abs(real(T1-T2,wr)),real(delta1,wr)/2.0_wr)
    ! ------------------------------------------------------------------    .
    if ((real(T3,wr).ne.-12345.0_wr).and.(real(T4,wr).ne.-12345.0_wr)) then
      Adata%andS='S'
      Adata%tpsR%secS=real(nzsec,wr)+ real(T3,wr) + real(nzmsec,wr)/1000.0_wr
    else
      Adata%andS='X'
      Adata%tpsR%secS=0.0_wr
    endif
    ! -------                                                  --------    .
    Adata%sigS=max(abs(real(T3-T4,wr)),real(delta1,wr)/2.0_wr)
    ! -------                                                  --------    .
    call basetime(Adata%tpsR)                                              ! respect de la base 60, 24, ...
    ! -------                                                  --------    .
    if (Adata%andS=='S') then
      write(*,1000)Adata%sta%staname,'P',Adata%typeonde,Adata%coefP,Adata%tpsR%date%year, &
        Adata%tpsR%date%month,Adata%tpsR%date%day,Adata%tpsR%date%hour,Adata%tpsR%date%min, &
        Adata%tpsR%secP,Adata%tpsR%secS,Adata%andS,Adata%coefS,Adata%sigP,Adata%sigS
    else
      write(*,1001)Adata%sta%staname,'P',Adata%typeonde,Adata%coefP,Adata%tpsR%date%year, &
        Adata%tpsR%date%month,Adata%tpsR%date%day,Adata%tpsR%date%hour,Adata%tpsR%date%min, &
        Adata%tpsR%secP,Adata%sigP
    endif
  endif
  ! --------------------------------------------------------------------    .
  ! ONDES REFRACTEES                                                        .
  ! --------------------------------------------------------------------    .
  ! T5 : temps arrivée de l'onde Pn
  ! T6 : incertitudes absolue sur l'onde Pn
  ! --------------------------------------------------------------------    .
  ! T7 : temps arrivée de l'onde Sn
  ! T8 : incertitudes absolue sur l'onde Sn
  ! ---------------------------------------------------------------------   .
  if ((real(T5,wr).ne.-12345.0_wr).and.(real(T6,wr).ne.-12345.0_wr)) then
    Adata%typeonde='N'
    Adata%coefP=0 ! modification a posteriori si necessaire
    Adata%coefS=0
    ! -------                                                  --------    . date
    Adata%tpsR%date%year=int(nzyear,wi)-2000
    ! -------                                                  --------    . jour julien
    Adata%tpsR%date%month=1
    Adata%tpsR%date%day=1
    call JDATE (JD,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
    call GDATE (JD+int(nzjday,wi)-1, Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
    ! -------                                                  --------    .
    Adata%tpsR%date%hour=nzhour
    Adata%tpsR%date%min=nzmin
    Adata%tpsR%secP=real(nzsec,wr)+ real(T5,wr) + real(nzmsec,wr)/1000.0_wr
    ! -------                                                  --------    .
    Adata%sigP=max(abs(real(T5-T6,wr)),real(delta1,wr)/2.0_wr)
    ! ------------------------------------------------------------------    .
    if ((real(T7,wr).ne.-12345.0_wr).and.(real(T8,wr).ne.-12345.0_wr)) then
      Adata%andS='S'
      Adata%tpsR%secS=real(nzsec,wr)+ real(T7,wr) + real(nzmsec,wr)/1000.0_wr
    else
      Adata%andS='X'
      Adata%tpsR%secS=0.0_wr
    endif
    ! -------                                                  --------    .
    Adata%sigS=max(abs(real(T7-T8,wr)),real(delta1,wr)/2.0_wr)
    ! -------                                                  --------    .
    call basetime(Adata%tpsR)                                              ! respect de la base 60, 24, ...
    ! -------                                                  --------    .
    if (Adata%andS=='S') then
      write(*,1000)Adata%sta%staname,'P',Adata%typeonde,Adata%coefP,Adata%tpsR%date%year, &
        Adata%tpsR%date%month,Adata%tpsR%date%day,Adata%tpsR%date%hour,Adata%tpsR%date%min, &
        Adata%tpsR%secP,Adata%tpsR%secS,Adata%andS,Adata%coefS,Adata%sigP,Adata%sigS
    else
      write(*,1001)Adata%sta%staname,'P',Adata%typeonde,Adata%coefP,Adata%tpsR%date%year, &
        Adata%tpsR%date%month,Adata%tpsR%date%day,Adata%tpsR%date%hour,Adata%tpsR%date%min, &
        Adata%tpsR%secP,Adata%sigP
    endif
  endif
  ! --------------------------------------------------------------------    . format de lecture des données
1000 format(a4,2a1,1x,i1,1x,5i2.2,f6.3,5x,f7.3,1x,a1,i1.1,3x,f6.3,4x,f6.3)
1001 format(a4,2a1,1x,i1,1x,5i2.2,f6.3,18x,f6.3)
  ! --------------------------------------------------------------------    .

end program readpick

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .

