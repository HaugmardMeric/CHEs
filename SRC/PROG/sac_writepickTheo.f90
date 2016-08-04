! programme principal sac_writepick
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! ------- writepick : lit les fichiers (phases list) lisible     --------   .
! ------- par le programme CHE et ecrit un fichier               --------   .
! ------- sac avec pick.sh                                       --------   .
! -------                                                        --------   .
! ------- fevrier 2015                                           --------   .
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
program writepickTheo
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
  character(LEN=4) :: astation
  character(LEN=112) :: file1
  character(LEN=112) :: fileppk
  ! ---------------------------------------------------------------------   . autres
  real(KIND=wr) :: alpha,beta
  integer(KIND=wi) :: JD, ok
  type(dataone) :: Adata
  ! --------------------------------------------------------------------    . lecture
  NN = IARGC()
  if ((NN < 2).or.(NN > 2)) then
    write(*,'(a)') "usage:  ./readpick.exe sacfile ppkdatafile"
    write(*,'(a)') "        sacfile - input sac file (.bin)"
    write(*,'(a)') "        ppkdatafile - input pick file (.txt)"
    stop
  endif
  call GETARG(1, file1)
  call GETARG(2, fileppk)
  ! -------                                                      --------   . lecture du fichier sac
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
  if (nvhdr /= 6) then
    write(*,*) "ERROR 1 - File: '", TRIM(adjustl(file1)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif

  ! --------------------------------------------------------------------    . temps du début de la trace, -> tpsTh <-
  Adata%tpsTh%date%year=int(nzyear,wi)
  ! -------                                                      --------   . jour julien
  Adata%tpsTh%date%month=1
  Adata%tpsTh%date%day=1
  call JDATE (JD,Adata%tpsTh%date%year,Adata%tpsTh%date%month,Adata%tpsTh%date%day)
  call GDATE (JD+int(nzjday,wi)-1, Adata%tpsTh%date%year,Adata%tpsTh%date%month,Adata%tpsTh%date%day)
  call JDATE (Adata%tpsTh%date%Jday,Adata%tpsTh%date%year,Adata%tpsTh%date%month,Adata%tpsTh%date%day)
  ! -------                                                      --------   .
  Adata%tpsTh%date%hour=nzhour
  Adata%tpsTh%date%min=nzmin
  Adata%tpsTh%secP=real(nzsec,wr) + real(nzmsec,wr)/1000.0_wr
  Adata%tpsTh%secS=Adata%tpsTh%secP
  call basetime(Adata%tpsTh)
  ! --------------------------------------------------------------------    . initialise
  T0=-12345.0_4
  a=-12345.0_4
  o=-12345.0_4
  ! --------------------------------------------------------------------    . lecture du fichier données picks, -> tpsR <-
  ok=0
  open(111, FILE = TRIM(adjustl(fileppk)),status='old',iostat = ok)
  if (ok .ne. 0) then
    write(*,*)'problème : le fichier '//TRIM(adjustl(fileppk))//' n''existe pas '
    stop
  endif
  ! -------                                                      --------   .
  read(111,1000,iostat = ok)Adata%sta%staname,Adata%typeonde,           &
    Adata%coefP,Adata%tpsR%date%year,Adata%tpsR%date%month,             &
    Adata%tpsR%date%day,Adata%tpsR%date%hour,Adata%tpsR%date%min,       &
    Adata%tpsR%secP,Adata%tpsR%secS,Adata%andS,Adata%coefS,             &
    Adata%sigP,Adata%sigS
  if (ok .ne. 0) then
    write(*,*)'fichier '//TRIM(adjustl(fileppk))//' vide ...'
    stop
  endif
  close(111)
  ! -------                                                      --------   .
  if((Adata%sigP.lt.0.0_wr).and.(IsNaN(Adata%sigP))) then
    write(*,*)'problème dans lectdata : les incertitudes sur les données P n''existent pas ',Adata%sigP
    stop
  endif
  if ((Adata%andS.ne."S").and.(Adata%sigS.lt.0.0_wr).and.(IsNaN(Adata%sigS))) then
    write(*,*)'problème dans lectdata : les incertitudes sur les données S n''existent pas ',Adata%sigS
    stop
  endif
  ! -------                                                      --------   .
  Adata%tpsR%date%year=Adata%tpsR%date%year+2000
  ! -------                                                      --------   . respect du decoupage des années en mois et jours avec prise en compte des années bisextiles
  call basetime(Adata%tpsR)
  call JDATE(Adata%tpsR%date%Jday,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
  call GDATE (Adata%tpsR%date%Jday,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
  call basetime(Adata%tpsR)                                                 ! respect du decoupage temps dans la base composite 60/12/365 ...
  ! -------                                                      --------   . verif station
  astation=TRIM(adjustl(kstnm))//'000'
  if (Adata%sta%staname.ne.astation) then
    write(*,*)'problème : les stations sont différentes : ',Adata%sta%staname,(TRIM(adjustl(kstnm))//'000')
    stop
  endif
  ! --------------------------------------------------------------------    .
  call difftime(alpha,beta,Adata%tpsR,Adata%tpsTh)
  ! --------------------------------------------------------------------    .
  ! ONDES DIRECTES THEORIQUES                                               .
  ! --------------------------------------------------------------------    .
  if (Adata%typeonde.eq.'G') then
    if(Adata%andS.eq.'S') then
      a=real(alpha,4)
      o=real(beta,4)
      if (.not.(o.gt.(0.0_4+b1))) then
        o=-12345.0_4
      endif
    else
      a=real(alpha,4)
    endif
    if (.not.(a.gt.(0.0_4+b1))) then
      a=-12345.0_4
      o=-12345.0_4
    endif
  endif
  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .
  call wbsac(file1,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,        &
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
    write(*,*) "ERROR 2 - File: '", TRIM(adjustl(file1)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif
  ! --------------------------------------------------------------------    . format de lecture des données
1000 format(a4,1x,a1,1x,i1,1x,5i2.2,f6.3,5x,f7.3,1x,a1,i1.1,3x,f6.3,4x,f6.3)
  ! --------------------------------------------------------------------    .

end program writepickTheo

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .

