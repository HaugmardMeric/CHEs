! programme principal sac_stalta_kurtosis
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! ------- stalta_kurtosis :                                      --------   .
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
! CONTAINS : subroutine Moment
!
program stalta_kurtosis
  ! --------------------------------------------------------------------    .
  use modparam
  use typetemps
  use time
  use sac_i_o
  ! -------                                                      --------   .
  implicit none
  ! ---------------------------------------------------------------------   . pour SAC
  real(k), dimension(:), allocatable :: sacfile1,sacfile2,sacfile3
  real(k), dimension(:), allocatable :: sacfilestalta, sacfilekurt, sacfilekurtdervivcentre
  real(k), dimension(:), allocatable :: F2,F3,F4
  real(k) :: delta1, b1, e1
  real(k) :: delta2, b2, e2
  real(k) :: delta3, b3, e3
  integer(k) :: NN, npts1, npts2, npts3
  character(LEN=112) :: file1,file2,file3
  character(LEN=112) :: outfile
  integer(k) :: anbpts
  ! ---------------------------------------------------------------------   . autres
  integer(KIND=wi) :: JD, inc, i,j
  integer(KIND=wi) :: sta, lta
  real(KIND = wr) :: nsec, ave, sdev, skew, curt1, curt2, curt3, nextmax
  real(KIND = wr), dimension(:), allocatable :: Asta, Alta, FC
  real(KIND=wr) :: deltaP1,deltaS1,deltaP2,deltaS2
  type(dataone) :: Adata1,Adata2,Adata3
  logical :: existe1
  ! --------------------------------------------------------------------    . lecture
  NN = IARGC()
  if ((NN < 1).or.(NN > 3).or.(NN == 2)) then
    write(*,'(a)') "usage:  ./sac_stalta_kurtosis sacfile (1 or 3 files : Z, N , E)"
    write(*,'(a)') "        sacfile - input sac file (.bin)",NN
    stop
  endif
  call GETARG(1, file1)
  inquire (file=file1,exist=existe1) ! option différente selon compilo !
  if (existe1) then
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
    kinst,sacfile1)
  if (nvhdr /= 6) then
    write(*,*) "ERROR 1 - File 1 : '", TRIM(adjustl(file1)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif
  allocate(FC(npts1))
  ! --------------------------------------------------------------------    . temps du début de la trace, -> tpsTh <-
  Adata1%tpsTh%date%year=int(nzyear,wi)
  ! -------                                                      --------   . jour julien
  Adata1%tpsTh%date%month=1
  Adata1%tpsTh%date%day=1
  call JDATE (JD,Adata1%tpsTh%date%year,Adata1%tpsTh%date%month,Adata1%tpsTh%date%day)
  call GDATE (JD+int(nzjday,wi)-1, Adata1%tpsTh%date%year,Adata1%tpsTh%date%month,Adata1%tpsTh%date%day)
  call JDATE (Adata1%tpsTh%date%Jday,Adata1%tpsTh%date%year,Adata1%tpsTh%date%month,Adata1%tpsTh%date%day)
  ! -------                                                      --------   .
  Adata1%tpsTh%date%hour=nzhour
  Adata1%tpsTh%date%min=nzmin
  Adata1%tpsTh%secP=real(nzsec,wr) + real(nzmsec,wr)/1000.0_wr
  Adata1%tpsTh%secS=Adata1%tpsTh%secP
  call basetime(Adata1%tpsTh)
  ! --------------------------------------------------------------------    .

  ! --------------------------------------------------------------------    .
  if (NN==3) then
    call GETARG(2,file2)
    call GETARG(3,file3)
    ! -------                                                    --------   . lecture du fichier sac
    call rbsac(file2,delta2,depmin,depmax,scale,odelta,b2,e2,o,a,internal1,        &
      t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
      resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
      user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
      internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
      unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
      nzsec,nzmsec,nvhdr,norid,nevid,npts2,internal4,nwfid,nxsize,nysize,unused8,  &
      iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
      imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
      unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
      kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
      kinst,sacfile2)
    if (nvhdr /= 6) then
      write(*,*) "ERROR 1 - File 2 : '", TRIM(adjustl(file2)), "' appears to be of non-native byte-order or is not a SAC file."
      stop
    endif
    ! ------------------------------------------------------------------    . temps du début de la trace, -> tpsTh <-
    Adata2%tpsTh%date%year=int(nzyear,wi)
    ! -------                                                    --------   . jour julien
    Adata2%tpsTh%date%month=1
    Adata2%tpsTh%date%day=1
    call JDATE (JD,Adata2%tpsTh%date%year,Adata2%tpsTh%date%month,Adata2%tpsTh%date%day)
    call GDATE (JD+int(nzjday,wi)-1, Adata2%tpsTh%date%year,Adata2%tpsTh%date%month,Adata2%tpsTh%date%day)
    call JDATE (Adata2%tpsTh%date%Jday,Adata2%tpsTh%date%year,Adata2%tpsTh%date%month,Adata2%tpsTh%date%day)
    ! -------                                                    --------   .
    Adata2%tpsTh%date%hour=nzhour
    Adata2%tpsTh%date%min=nzmin
    Adata2%tpsTh%secP=real(nzsec,wr) + real(nzmsec,wr)/1000.0_wr
    Adata2%tpsTh%secS=Adata2%tpsTh%secP
    call basetime(Adata2%tpsTh)
    ! ------------------------------------------------------------------    .

    ! -------                                                    --------   . lecture du fichier sac
    call rbsac(file3,delta3,depmin,depmax,scale,odelta,b3,e3,o,a,internal1,        &
      t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
      resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
      user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
      internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
      unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
      nzsec,nzmsec,nvhdr,norid,nevid,npts3,internal4,nwfid,nxsize,nysize,unused8,  &
      iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
      imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
      unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
      kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
      kinst,sacfile3)
    if (nvhdr /= 6) then
      write(*,*) "ERROR 1 - File 3 : '", TRIM(adjustl(file3)), "' appears to be of non-native byte-order or is not a SAC file."
      stop
    endif
    ! ------------------------------------------------------------------    . temps du début de la trace, -> tpsTh <-
    Adata3%tpsTh%date%year=int(nzyear,wi)
    ! -------                                                    --------   . jour julien
    Adata3%tpsTh%date%month=1
    Adata3%tpsTh%date%day=1
    call JDATE (JD,Adata3%tpsTh%date%year,Adata3%tpsTh%date%month,Adata3%tpsTh%date%day)
    call GDATE (JD+int(nzjday,wi)-1, Adata3%tpsTh%date%year,Adata3%tpsTh%date%month,Adata3%tpsTh%date%day)
    call JDATE (Adata3%tpsTh%date%Jday,Adata3%tpsTh%date%year,Adata3%tpsTh%date%month,Adata3%tpsTh%date%day)
    ! -------                                                    --------   .
    Adata3%tpsTh%date%hour=nzhour
    Adata3%tpsTh%date%min=nzmin
    Adata3%tpsTh%secP=real(nzsec,wr) + real(nzmsec,wr)/1000.0_wr
    Adata3%tpsTh%secS=Adata3%tpsTh%secP
    call basetime(Adata3%tpsTh)
    ! ------------------------------------------------------------------    . verification 1
    call difftime(deltaP1,deltaS1,Adata1%tpsTh,Adata2%tpsTh)
    call difftime(deltaP2,deltaS2,Adata1%tpsTh,Adata3%tpsTh)
    if ((deltaP1+deltaS2).gt.(real(delta1,wr)*4.0_wr)) then
      write(*,*)'problème dans stalta_kurtosis, fichiers non homogènes en temps '
      write(*,*)Adata1%tpsTh
      write(*,*)Adata2%tpsTh,deltaP1
      write(*,*)Adata3%tpsTh,deltaP2
      write(*,*)delta1
      stop
    endif
    ! ------------------------------------------------------------------    . verification 2
    if ((npts1.ne.npts2).or.(npts1.ne.npts3)) then
      write(*,*)'problème dans stalta_kurtosis, fichiers non homogènes en durée '
      write(*,*)npts1,npts2,npts3
      stop
    endif
    ! ------------------------------------------------------------------    . verification 3
    if ((delta1.ne.delta2).or.(delta1.ne.delta3)) then
      write(*,*)'problème dans stalta_kurtosis, fichiers non homogènes en taux échantillonnage'
      write(*,*)delta1,delta2,delta3
      stop
    endif
  endif
  ! --------------------------------------------------------------------    .




  ! --------------------------------------------------------------------    .
  ! KURTOSIS
  ! cf : Baillard, W. C. Crawford, V. Ballu, C. Hibert, & A. Mangeney (2014) :
  !      An Automatic Kurtosis-Based P- and S-Phase Picker Designed for Local Seismic Networks.
  !      BSSA, Vol. 104, No. 1, 16p.
  ! --------------------------------------------------------------------    .
  allocate(sacfilekurt(npts1))
  allocate(sacfilekurtdervivcentre(npts1))
  allocate(F2(npts1),F3(npts1),F4(npts1))
  nsec=2.5_wr                                                               ! 1.0 secondes
  inc=int(nsec/real(delta1,wr)+1.0_wr,wi)
  sacfilekurt=real(0.0,k)
  anbpts=int(inc,k)
  ! -------                                                      --------   .
  ! The central moment of order d at sample k can be written as
  ! -------                                                      --------   .
  do i=anbpts,npts1
    call Moment(sacfile1(i-anbpts+1:i),anbpts,ave,sdev,skew,curt1)
    if (NN==3) then
      call Moment(sacfile2(i-anbpts+1:i),anbpts,ave,sdev,skew,curt2)
      call Moment(sacfile3(i-anbpts+1:i),anbpts,ave,sdev,skew,curt3)
      sacfilekurt(i)=real((curt1+curt2+curt3)/3.0_wr,k)
    else
      sacfilekurt(i)=real(curt1,k)
    endif
  enddo
  ! -------                                                      --------   .
  ! The first transformation essentially cleans the initial CF
  ! of all strictly negative gradients (Fig. 2c), because only
  ! positive gradients characterize the transition from noise to
  ! a coherent signal.
  ! -------                                                      --------   .
  do i=1,npts1-1
   if ((sacfilekurt(i+1)-sacfilekurt(i)).gt.0.0_k) then
    F2(i+1)=F2(i)+(sacfilekurt(i+1)-sacfilekurt(i))
   else
    F2(i+1)=F2(i)
   endif
  enddo
  F2(npts1)=F2(npts1-1)
  ! -------                                                      --------   .
  ! The second transformation removes a linear trend from
  ! F2, so that the first and last values equal zero.
  ! In this way, the onsets become local minima.
  ! -------                                                      --------   .
  do i=1+1,npts1
    F3(i)=F2(i)-((F2(npts1)-F2(1))/(real(npts1-1,k))*real(i-1,k)+F2(1))
  enddo
  ! -------                                                      --------   .
  ! The final transformation makes the amplitude of the
  ! minima amplitude scale with the total change in the kurtosis
  ! that follows, so that the greatest minima correspond to the
  ! greatest onset strengths
  ! -------                                                      --------   .
  do i=1,npts1-2
    ! prochain maxima
    do j=i+1,npts1
      if (F3(j+1).lt.F3(j)) then
        nextmax=F3(j)
        exit
      endif
    enddo
    if ((F3(i)-nextmax).lt.0.0_k) then
      F4(i)=F3(i)-real(nextmax,k)
    else
      F4(i)=0.0_k
    endif
  enddo
  F4(npts1)=0.0_k
  F4(npts1-1)=0.0_k
  F4(npts1-2)=0.0_k
  ! -------                                                       --------   .
  !Dérivée centré d'ordre 2
  !do i=2,npts1-1
  !sacfilekurtdervivcentre(i)=(sacfilekurt(i+1)-sacfilekurt(i-1))/(2.0_k*delta1)
  !enddo
  !sacfilekurtdervivcentre(1)=sacfilekurtdervivcentre(2)
  !sacfilekurtdervivcentre(npts1)=sacfilekurtdervivcentre(npts1)
  ! --------------------------------------------------------------------    .
  deallocate(F2,F3)
  ! --------------------------------------------------------------------    .



  ! --------------------------------------------------------------------    .
  ! STA / LTA                                                               .
  ! --------------------------------------------------------------------    .
  if (NN==3) then
    ! -------                                                    --------   . fonction caractéristique
    FC=(real(sacfile1,wr)**2.0_wr+real(sacfile2,wr)**2.0_wr+real(sacfile3,wr)**2.0_wr)/3.0_wr
  else
    ! -------                                                    --------   . fonction caractéristique : energie
    FC=real(sacfile1,wr)**2.0_wr
   endif
  ! --------------------------------------------------------------------    .
  allocate(sacfilestalta(npts1))
  allocate(Asta(npts1))
  allocate(Alta(npts1))
  Asta=0.1e-9_wr
  Alta=0.1e-9_wr
  nsec=3.0_wr                                                               ! 2.5 secondes pour sta
  sta=int(nsec/real(delta1,wr)+1.0_wr,wi)
  nsec=30.0_wr                                                              ! 25 secondes pour lta
  lta=int(nsec/real(delta1,wr)+1.0_wr,wi)
  sacfilestalta=real(0.0,k)
  ! --------------------------------------------------------------------    . initialisation STA
  Asta(lta+sta+1)=0.00001_wr
  do i=1,lta
   Alta(lta+sta+1)=Alta(lta+sta+1)+FC(i)
  enddo
  Alta(lta+sta+1)=Alta(lta+sta+1)/real(lta,wr)
  ! --------------------------------------------------------------------    . initialisation LTA
  Alta(lta+sta+1)=0.00001_wr
  do i=lta+1,lta+sta
    Asta(lta+sta+1)=Asta(lta+sta+1)+FC(i)
  enddo
  Asta(lta+sta+1)=Asta(lta+sta+1)/real(sta,wr)
  ! --------------------------------------------------------------------    . STA - LTA
  do i=lta+sta+sta,npts1-1
    ! démarre à + 2 sta sinon instable
    Alta(i)= FC(i-sta-1)/real(lta,wr) + (1.0_wr-1.0_wr/real(lta,wr))*Alta(i-1)
    Asta(i)= FC(i-1)/real(sta,wr) + (1.0_wr-1.0_wr/real(sta,wr))*Asta(i-1)
    sacfilestalta(i)=real(Asta(i)/Alta(i),k)
  enddo

  deallocate(Asta,Alta,FC)
  ! --------------------------------------------------------------------    .

  ! --------------------------------------------------------------------    .
  ! ecriture des résultats                                                   .
  ! --------------------------------------------------------------------    .
  outfile=TRIM(adjustl(file1))//'.kurt'
  ! --------------------------------------------------------------------    .
  call wbsac(outfile,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,      &
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
    kinst,F4)

  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .
  outfile=TRIM(adjustl(file1))//'.stalta'
  ! --------------------------------------------------------------------    .
  call wbsac(outfile,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,      &
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
    kinst,sacfilestalta)
  ! --------------------------------------------------------------------    .
  if (nvhdr /= 6) then
    write(*,*) "ERROR 2 - File: '", TRIM(adjustl(file1)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif
  ! --------------------------------------------------------------------    .
  deallocate(sacfilestalta,sacfilekurt,sacfilekurtdervivcentre,F4)
  else
    write(*,*) "NO FILE ... ",file1
  endif
  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .

CONTAINS

  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .

  subroutine Moment(Adata1,n,ave,sdev,skew,curt)
    ! ------------------------------------------------------------------    .
    ! calcul de la moyenne (ave), de l'écart-type (sdev)
    ! du coefficient de dissymétrie (skew) et
    ! du coefficient d’aplatissement de Pearson (curt pour kurtosis)
    ! ------------------------------------------------------------------    .
    implicit none
    ! ------------------------------------------------------------------    .
    integer(k), intent (in) :: n
    real(KIND=k), intent (in) :: Adata1(n)
    real(KIND=wr), intent (out) ::   ave,sdev,skew,curt
    ! -------                                                   --------    .
    integer(KIND=wi) :: j
    real(KIND=wr) :: var,p,s
    ! ------------------------------------------------------------------    .
    if(n.le.1) then
      write(*,*)'problème dans Moment : n < 2 ! ',n
      stop
    end if
    s=0.0_wr
    do j=1,n
      s=s+real(Adata1(j),wr)
    end do
    ave=s/real(n,wr)                                                        ! moyenne
    ! ------------------------------------------------------------------    .
    var=0.0_wr
    skew=0.0_wr
    curt=0.0_wr
    do j=1,n
      s=real(Adata1(j),wr)-ave
      p=s*s
      var=var+p
      p=p*s
      skew=skew+p
      p=p*s
      curt=curt+p
    end do
    ! ------------------------------------------------------------------    .
    !var=var/real(n-1,wr)                                                   ! variance echantillon
    var=var/real(n,wr)                                                      ! variance population
    sdev=sqrt(var)                                                          ! ecart-type
    if(var.ne.0.0_wr) then
      skew=skew/(real(n,wr)*sdev**3.0_wr)                                   ! coefficient de dissymétrie
      curt=curt/(real(n,wr)*var**2.0_wr)-3.0_wr                             ! coefficient d’aplatissement de Pearson
    else
      skew=-1.0e9_wr
      curt=-1.0e9_wr
    end if
    ! ------------------------------------------------------------------    .
    ! write(*,*) n,ave,sdev,skew,curt
    ! ------------------------------------------------------------------    .
  end subroutine Moment

  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .

end program stalta_kurtosis

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .

