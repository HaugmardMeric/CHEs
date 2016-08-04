! programme principal sac_writepick
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! ------- writepick : lit le catalogue                           --------   .
! ------- inscrit les pointés théoriques dans les fichiers sac   --------   .
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
program writepickCata
  ! --------------------------------------------------------------------    .
  use modparam
  use typetemps
  use time
  use sac_i_o
  use datalecture
  use pb_direct
  ! -------                                                      --------   .
  implicit none
  ! ---------------------------------------------------------------------   . pour SAC
  real(KIND=4), dimension(:), allocatable :: sacfile
  real(KIND=4) :: delta1, b1, e1
  integer(KIND=4) :: NN, npts1
  character(LEN=112) :: file1
  ! ---------------------------------------------------------------------   .
  type(parametre) :: param                                                  ! paramètres théoriques du séisme
  type(seismes) :: theseisme(2)                                             ! paramètres du catalogue
  integer(KIND=wi) :: find                                                  ! ce séisme est présent dans le catalogue, 0-1-2 fois
  integer(KIND=wi) :: i,j,count,nbsta,nbtps(nbseismes)
  type(stations), allocatable, dimension(:) :: datasta
  type(dataone) :: datatemps                                                ! données
  logical :: critique                                                       ! .true. si distance hypo + 5 km < distance hypo critique pour la réfraction
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  ! ---------------------------------------------------------------------   . autres
  real(KIND=wr) :: alpha,beta
  integer(KIND=wi) :: JD, ok
  type(dataone) :: Adata
  ! --------------------------------------------------------------------    . lecture
  NN = IARGC()
  if ((NN < 1).or.(NN > 1)) then
    write(*,'(a)') "usage:  ./readpick.exe sacfile ppkdatafile"
    write(*,'(a)') "        sacfile - input sac file (.bin)"
    stop
  endif
  call GETARG(1, file1)
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
  ! ------- si moho non tabulaire                                --------    .
  acentroid%lonC=moho_lon ; acentroid%latC=moho_lat
  acentroid%NS=moho_NS ; acentroid%EO=moho_EO
  call alph2vect(acentroid); call vect2alph(acentroid)
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
  Adata%sta%staname=TRIM(adjustl(kstnm))//'000'
  ! --------------------------------------------------------------------    .
  ! lecture du catalogue
  ! --------------------------------------------------------------------    . d'après la géol :
  param%VC=6.0_wr
  param%VM=8.0_wr
  param%Zmoho=30.0_wr
  param%VpVs=1.71_wr
  ! --------------------------------------------------------------------    . au pif :
  param%Zhypo=5.0_wr
  param%lon=-2.0_wr
  param%lat=47.5_wr
  ! --------------------------------------------------------------------    .
  param%Tzero%date=Adata%tpsTh%date
  param%Tzero%sec=Adata%tpsTh%secP+60.0_wr                                  ! la trace commence 60 s avant le séisme du catalogue
  call basetime(param%Tzero)
  ! -------                                                      --------   .
  call catalogue(param,theseisme,find)
  if (find.lt.1) then                                                       ! deux séismes max !
    write(*,*)'problème dans writepick : séismes au catalogue'
    stop
  elseif (find.gt.2) then
    write(*,*)'attention dans writepick : plusieurs séismes au catalogue',find
    write(*,*)theseisme(1)
    write(*,*)theseisme(2)
  endif
  ! -------                                                      --------   .
  param%Zhypo=theseisme(1)%pfd
  param%lon=theseisme(1)%lon
  param%lat=theseisme(1)%lat
  param%Tzero%date = theseisme(1)%tps_init%date
  param%Tzero%sec = theseisme(1)%tps_init%sec
  call basetime(param%Tzero)
  ! --------------------------------------------------------------------    .
  ! lire la station
  ! --------------------------------------------------------------------    .
  call lectnbdata(nbsta,nbtps)
  allocate(datasta(nbsta))
  ! -------                                                      --------   .
  open(503, FILE = 'DATA/sta.d',status='old',iostat = ok)
  if (ok .ne. 0) then
    open(503, FILE = 'sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème : le fichier DATA/sta.d n''existe pas '
      stop
    endif
  endif
  ! -------                                                      --------   .
  do i=1,nbsta
    read(503,*,iostat = ok) datasta(i)
  enddo
  close(503)
  ! -------                                                      --------   .
  count=0
  do j=1,nbsta
    if(Adata%sta%staname.eq.datasta(j)%staname) then
      datatemps%sta=datasta(j)                                              ! attribution d'une station
      count=count+1
    endif
  enddo
  ! -------                                                      --------   .
  if(count.gt.1) then                                                       ! vérification de l'absence de doublons
    write(*,*)'problème : station ',datatemps%sta%staname,' en double in file : DATA/sta.d'
    stop
  elseif(count.eq.0) then
    write(*,*)'problème : station ',datatemps%sta%staname,' non répertoriée'
    stop
  endif
  ! -------                                                      --------   .
  deallocate(datasta)
  ! --------------------------------------------------------------------    .
  ! distance épi, problème direct
  ! --------------------------------------------------------------------    .
  datatemps%typeonde='G'
  datatemps%andS='S'
  call tempsTheoDirectone(param,datatemps,critique,acentroid)
  ! --------------------------------------------------------------------    .
  ! tps arrivée, tpsR = Tzero (catalogue) + parcours
  ! --------------------------------------------------------------------    .
  Adata%tpsR%date=param%Tzero%date
  Adata%tpsR%secP=param%Tzero%sec+datatemps%tpsparcP
  Adata%tpsR%secS=param%Tzero%sec+datatemps%tpsparcS
  ! -------                                                      --------   . respect du decoupage des années en mois et jours avec prise en compte des années bisextiles
  call basetime(Adata%tpsR)
  call JDATE(Adata%tpsR%date%Jday,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
  call GDATE (Adata%tpsR%date%Jday,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
  call basetime(Adata%tpsR)                                                 ! respect du decoupage temps dans la base composite 60/12/365 ...
  ! --------------------------------------------------------------------    .
  ! ondes P et S directes
  ! --------------------------------------------------------------------    .
  call difftime(alpha,beta,Adata%tpsR,Adata%tpsTh)
  ! --------------------------------------------------------------------    .
  ! ONDES DIRECTES THEORIQUES                                               .
  ! --------------------------------------------------------------------    .
  a=real(alpha,4)
  o=real(beta,4)
  if (.not.(o.gt.(0.0_4+b1))) then
    o=-12345.0_4
  endif
  ! -------                                                      --------   .
  if (.not.(a.gt.(0.0_4+b1))) then
    a=-12345.0_4
    o=-12345.0_4
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

end program writepickCata

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .

