! programme principal sac_coda
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! ------- calcul de la magnitude Md d'un séisme                  --------   .
! ------- d'après la durée de la coda (formule de Lee et al., 1972) ----    .
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
program coda
  ! --------------------------------------------------------------------    .
  use modparam
  use typetemps
  use time
  use statistiques
  use sac_i_o
  ! -------                                                      --------   .
  implicit none
  ! -------                                                      --------   .
  type(dataone) :: Adata,pick,ref
  type(date_secPS) :: tps0
  ! ---------------------------------------------------------------------   .
  real(KIND=wr), parameter :: dmax=300.0_wr                                 ! distance maximale de prise en compte
  ! ---------------------------------------------------------------------   . pour SAC
  real(KIND=4), dimension(:), allocatable :: sacfile
  real(KIND=4) :: delta1, b1, e1
  integer(KIND=4) :: NN, ios, npts1, j
  character(LEN=112) :: file1,ofile
  ! ---------------------------------------------------------------------   . autres
  type(date_sec) :: tpsref
  real(KIND=wr) :: ttp,tts,deltatps
  real(KIND=wr) :: AmpliTot,AmpliMin,AmpliMax,AmpliNORM
  real(KIND=wr) :: ML,moyP,ecP,moyS
  integer(KIND=wi) :: i, n, jD, numero,nP, nS
  character (LEN=5) :: nbseisme
  logical :: existe1
  ! --------------------------------------------------------------------    . lecture du fichier sac
  NN = IARGC()
  if ((NN > 2).or.(NN < 2)) then
    write(*,'(a)') "usage:  sac_coda sacfile"
    write(*,'(a)') "        sacfile - input sac file"
    write(*,'(a)') "        ofile - output xy file"
    stop
  endif
  call GETARG(1, file1)
  open(UNIT=1,FILE=file1,STATUS='OLD',IOSTAT=ios)
  if (ios > 0) then
    write(*,*) "ERROR - Input file: '", TRIM(adjustl(file1)), "' does not exist ..."
    close(1)
    stop
  endif
  close(1)
  call GETARG(2, ofile)
  ! --------------------------------------------------------------------    . lecture des données : numéro du séisme, Pg et Sg, d_épi
  read(*,*)numero,pick%tpsR,tpsref,pick%depi,pick%dhypo,AmpliNORM
  tps0%date=tpsref%date
  tps0%secP=tpsref%sec
  tps0%secS=0.0_wr
  call basetime(tps0)
  call basetime(pick%tpsR)
  ! --------------------------------------------------------------------    .
  if (pick%depi.lt.dmax) then                                               ! distance maximale
    ! ------------------------------------------------------------------    .
    open(UNIT=2,FILE=ofile,STATUS='REPLACE',IOSTAT=ios)
    if (ios > 0) then
      write(*,*) "ERROR - Input file: '", TRIM(adjustl(ofile)), "' does not exist ..."
      stop
    endif
    CALL rbsac(file1,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,        &
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
      write(*,*) "ERROR - File: '", TRIM(adjustl(file1)), "' appears to be of non-native &
      &byte-order or is not a SAC file."
      stop
    endif
    pick%sta%staname=(TRIM(adjustl(kstnm))//'000')                          ! nom de la station
    ! ------------------------------------------------------------------    .
    ! ------------------------------------------------------------------    . normalisation amplitude
    AmpliMin=1.e9_wr
    AmpliMax=-1.e9_wr
    do j=0,(npts1-1)
      sacfile(j+1)=sqrt(sacfile(j+1))
      if (AmpliMin.gt.real(sacfile(j+1),wr)) AmpliMin=real(sacfile(j+1),wr)
      if (AmpliMax.lt.real(sacfile(j+1),wr)) AmpliMax=real(sacfile(j+1),wr)
    enddo
    AmpliTot=max(abs(AmpliMin),abs(AmpliMax))
    do j=0,(npts1-1)
      sacfile(j+1)=sacfile(j+1)*real(AmpliNORM/AmpliTot,4)
    enddo
    ! -------                                                    --------   . temps initial
    Adata%tpsR%date%year=int(nzyear,4)
    ! -------                                                   --------    . jour julien
    Adata%tpsR%date%month=1
    Adata%tpsR%date%day=1
    call jDATE (jD,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
    call GDATE (jD+int(nzjday,4)-1,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
    ! -------                                                   --------    .
    Adata%tpsR%date%hour=int(nzhour,4)
    Adata%tpsR%date%min=int(nzmin,4)
    Adata%tpsR%secP=real(nzsec,wr)+real(b1,wr)
    Adata%tpsR%secS=0.0_wr
    call jDATE (Adata%tpsR%date%jday,Adata%tpsR%date%year,Adata%tpsR%date%month,Adata%tpsR%date%day)
    call basetime(Adata%tpsR)
    ! ------------------------------------------------------------------    . moyenne du bruit, avant la P
    moyP=0.0_wr
    nP=0
    ref=Adata
    ref%tpsR%secP = ref%tpsR%secP - real(delta1,wr)
    do j=0,(npts1-1)
      ref%tpsR%secP = ref%tpsR%secP + real(delta1,wr)
      call basetime(ref%tpsR)
      call difftime(ttp,tts,pick%tpsR,ref%tpsR)
      if ((j.gt.int(35.0_wr/real(delta1,wr),wi)).and.(ttp.ge.1.0_wr).and.(ttp.le.600.0_wr)) then
        ! 35 sec apres début signal (taper) ; au moins 1 sec . avant la P (Pg ou Pn) ; moins de 600 secondes avant la P
        np=nP+1
        moyP=moyP+real(abs(sacfile(j+1)),wr)
        endif
    enddo
    if (nP.gt.(100)) then
      moyP=moyP/real(nP,wr)
      ! ----------------------------------------------------------------    . ecartype sur la moyenne du bruit, avant la P
      nP=0
      ecP=0.0_wr
      ref=Adata
      ref%tpsR%secP = ref%tpsR%secP - real(delta1,wr)
      do j=0,(npts1-1)
        ref%tpsR%secP = ref%tpsR%secP + real(delta1,wr)
        call basetime(ref%tpsR)
        call difftime(ttp,tts,pick%tpsR,ref%tpsR)
        if ((j.gt.int(35.0_wr/real(delta1,wr),wi)).and.(ttp.ge.1.0_wr).and.(ttp.le.600.0_wr)) then
          ! 35 sec apres début signal (taper) ; au moins 1 sec . avant la P (Pg ou Pn) ; moins de 600 secondes avant la P
          np=nP+1
          ecP=ecP+(moyP-real(abs(sacfile(j+1)),wr))**2.0_wr
        endif
      enddo
      ecP=sqrt(ecP/real(nP,wr))
      ! ----------------------------------------------------------------    . moyenne glissante sur la coda
      nS=0
      ref=Adata
      ref%tpsR%secP = ref%tpsR%secP - real(delta1,wr)
      boucles : do j=0,(npts1-50)
        ref%tpsR%secP = ref%tpsR%secP + real(delta1,wr)
        ref%tpsR%secS = ref%tpsR%secP
        call basetime(ref%tpsR)
        call difftime(ttp,tts,pick%tpsR,ref%tpsR)
        if ((ttp.le.0.0_wr).and.(j.gt.int(35.0_wr/real(delta1,wr),wi)))then ! après la première P et 35 sec apres le début du signal (taper)
          ns=ns+1
          moyS=0.0_wr
          n=0
          ! moyenne sur les 10 sec d'avant
          do i=-int(10.0_wr/real(delta1,wr),wi),0
            n=n+1
            moyS=moyS+real(abs(sacfile(j+1+i)),wr)
          enddo
          moyS=moyS/real(n,wr)
          call difftime(deltatps,Ml,ref%tpsR,tps0)
          write(2,*) pick%dhypo+moyS,deltatps
          if ((tts.le.0.0_wr).and.(moyS.lt.(2.0_wr*(moyP+0.0_wr*ecP)))) then ! après la première S
            ! magnitude durée (Md), la fin de la coda = la moitié du bruit avant le séisme, cf :
            ! [Kayal, j.R. (2008) : Microearthquake Seismology and Seismotectonics of South Asia , springer (§3.15)]
            exit boucles
          endif
        endif
      enddo boucles
      write(2,*)
      close(2)
      ! ----------------------------------------------------------------    .
      if ((ns.gt.int(1.0_wr/real(delta1,wr),wi))) then                      ! au moins 1.0 sec de signal
        write(nbseisme(1:5),'(i5)')numero
          inquire (file="OUTPUT/files/mag-"//trim(adjustl(nbseisme))//".d",exist=existe1)
          if (existe1) then
            open(UNIT=101,FILE="OUTPUT/files/mag-"//trim(adjustl(nbseisme))//".d", &
            status="old",position="append",action="write")
          else
            open(UNIT=101,FILE="OUTPUT/files/mag-"//trim(adjustl(nbseisme))//".d",STATUS='new')
          endif
          ! formule de Lee et al. (1972)
          Ml = -0.87_wr + 2.0_wr*log10(real(ns,wr)*real(delta1,wr)) + 0.0035_wr*pick%depi
          write(101,10011)kstnm,Ml,real(ns,wr)*real(delta1,wr),pick%depi
          close(101)
      endif
      ! ----------------------------------------------------------------    .
    endif
    ! ------------------------------------------------------------------    .
  endif
  ! --------------------------------------------------------------------    .
10011 format (a4,3(1x,f19.6))
  ! --------------------------------------------------------------------    .

end program coda

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .



