! programme principal sac_ZNE_2_LQT
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! -------                                                        --------   .
! ------- mars 2015                                              --------   .
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
program ZNE_2_LQT
  ! --------------------------------------------------------------------    .
  use modparam
  use distance_epi
  use sac_i_o
  ! -------                                                      --------   .
  implicit none
  real (kind=wr) :: dlat1,dlat2,dlon1,dlon2,d,alti,pfd
  real (kind=wr) :: a_baz                                                   ! backazimuth measured clockwise from north
  real (kind=wr) :: alpha                                                   ! angle of incidence, measured from vertical
  real (kind=wr) :: M3D(3,3), val(7)
  integer(KIND=wi) :: ok
  character(LEN=4) :: Asta
  logical :: found
  ! ---------------------------------------------------------------------   . pour SAC
  real(k), dimension(:), allocatable :: sacZ, sacN, sacE, sacL, sacQ, sacT
   CHARACTER(LEN=8) :: kstnm1,kstnm2,kstnm3
  real(k) :: delta1, b1, e1
  real(k) :: delta2, b2, e2
  real(k) :: delta3, b3, e3
  integer(k) :: NN, npts1, npts2, npts3
  character(LEN=112) :: file1,file2,file3
  character(LEN=112) :: outfile
  ! --------------------------------------------------------------------    .
  ! Plesinger, A., M. Hellweg and D. Seidl (1986) : Interactive high-resolution polarization analysis of broadband seismograms. J. Geophysics, 59, p. 129-139. 
  ! (http://service.iris.edu/irisws/rotation/docs/1/help/)
  ! --------------------------------------------------------------------    .
  ! Z , E , and N represents the 3 seismograms with original orientations
  ! --------------------------------------------------------------------    .
  ! L, Q, and T represent the three seismograms that are output as below.
  ! L – Aligned in direction of P wave propagation
  ! Q – Aligned in the direction of the SV phase movement
  ! T – Aligned in the direction of the SH phase movement
  ! --------------------------------------------------------------------    .
  ! ------- Lecture des données                                 --------    .
  ! --------------------------------------------------------------------    .
  NN = IARGC()
  if ((NN > 3).or.(NN < 3)) then
    write(*,'(a)') "usage:  ./ZNE_2_LQT.exe sacfileZ sacfileN sacfileE latEvent LonEvent pdfEvent"
    write(*,'(a)') "        sacfile - input sac file (.bin)"
    stop
  endif
  call GETARG(1,file1)
  call GETARG(2,file2)
  call GETARG(3,file3)
  ! --------------------------------------------------------------------    . lecture de l'hypocentre
  write(*,*)"hypocentre : lat(deg),lon(deg),pfd(m)"
  read(*,*)dlat2,dlon2,pfd

  !  dlat2=47.657_wr fixe ?
  !  dlon2=-2.8_wr
  !  pfd=4.0_wr

  ! --------------------------------------------------------------------    . lecture du fichier sac Z
  call rbsac(file1,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,        &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
    nzsec,nzmsec,nvhdr,norid,nevid,npts1,internal4,nwfid,nxsize,nysize,unused8,  &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm1,khole,ko,ka,kt0,kt1,&
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
    kinst,sacZ)
  ! -------                                                      --------   .
  if (nvhdr /= 6) then
    write(*,*) "ERROR 1 - File Z : '", TRIM(adjustl(file1)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif
  ! -------                                                      --------   . verification de la composante
  if((kcmpnm(3:3).ne."Z").and.(kcmpnm(3:3).ne."1")) then
    write(*,*)'problème mauvaise composante -> Z,1'
    write(*,*)kcmpnm
    stop
  endif
  ! --------------------------------------------------------------------    . lecture du fichier sac N
  call rbsac(file2,delta2,depmin,depmax,scale,odelta,b2,e2,o,a,internal1,        &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
    nzsec,nzmsec,nvhdr,norid,nevid,npts2,internal4,nwfid,nxsize,nysize,unused8,  &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm2,khole,ko,ka,kt0,kt1,&
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
    kinst,sacN)
  ! -------                                                      --------   .
  if (nvhdr /= 6) then
    write(*,*) "ERROR 1 - File N : '", TRIM(adjustl(file2)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif
  ! -------                                                      --------   . verification de la composante
  if((kcmpnm(3:3).ne."N").and.(kcmpnm(3:3).ne."2")) then
    write(*,*)'problème mauvaise composante -> N,2'
    write(*,*)kcmpnm
    stop
  endif
  ! --------------------------------------------------------------------    .  lecture du fichier sac E
  call rbsac(file3,delta3,depmin,depmax,scale,odelta,b3,e3,o,a,internal1,        &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
    nzsec,nzmsec,nvhdr,norid,nevid,npts3,internal4,nwfid,nxsize,nysize,unused8,  &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm3,khole,ko,ka,kt0,kt1,&
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
    kinst,sacE)
  ! -------                                                      --------   .
  if (nvhdr /= 6) then
    write(*,*) "ERROR 1 - File E : '", TRIM(adjustl(file3)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif
  ! -------                                                      --------   . verification de la composante
  if((kcmpnm(3:3).ne."E").and.(kcmpnm(3:3).ne."3")) then
    write(*,*)'problème mauvaise composante -> E,3'
    write(*,*)kcmpnm
    stop
  endif
  ! --------------------------------------------------------------------    . verification 1
  if ((npts1.ne.npts2).or.(npts1.ne.npts3)) then
    write(*,*)'problème dans stalta_ZNE_2_LQT, fichiers non homogènes en durée '
    write(*,*)npts1,npts2,npts3
    stop
  endif
  ! --------------------------------------------------------------------    . verification 2
  if ((delta1.ne.delta2).or.(delta1.ne.delta3)) then
    write(*,*)'problème dans stalta_ZNE_2_LQT, fichiers non homogènes en taux échantillonnage'
    write(*,*)delta1,delta2,delta3
    stop
  endif
  ! --------------------------------------------------------------------    . verification 3
  if ((kstnm1.ne.kstnm2).or.(kstnm1.ne.kstnm3)) then
    write(*,*)'problème dans stalta_ZNE_2_LQT, differentes stations'
    write(*,*)kstnm1,kstnm2,kstnm3
    stop
  endif
  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .

  ! --------------------------------------------------------------------    .
  ! ------- station                                             --------    .
  ! --------------------------------------------------------------------    .
  found=.false.
  open(503, FILE = 'sta.d',status='old',iostat = ok)
  if (ok .ne. 0) then
    write(*,*)'problème : le fichier sta.d n''existe pas '
    stop
  endif
  ! --------------------------------------------------------------------    .
  do while(ok.eq.0)
    read(503,*,iostat = ok) Asta,val

    if ((TRIM(adjustl(Asta))==TRIM(adjustl(kstnm1(1:4)))).or.(TRIM(adjustl(Asta))==TRIM(adjustl(kstnm1(1:4)))//"0")) then ! station en trois letters + "0"
      if(found) then
        write(*,*)'problème : la station '//TRIM(adjustl(kstnm1))//' est présente deux fois dans sta.d'
        stop
      else
        dlat1=val(1)
        dlon1=val(2)
        alti=val(3)
        found=.true.
      endif
    endif
  enddo
  close(503)
  if (.not.found) then
    write(*,*)'problème : la station '//TRIM(adjustl(kstnm1))//' n''existe pas '
    stop
  endif

  ! --------------------------------------------------------------------    .
  ! ------- rotation                                            --------    .
  ! --------------------------------------------------------------------    .
  call dellipsgc(dlat2,dlon2,dlat1,dlon1,d,a_baz)
  ! write(*,*)d,a_baz,alpha/pi*180.0_wr
  ! --------------------------------------------------------------------    .
  ! incidence is the angle from vertical at which an incoming ray arrives. 
  ! A ray arriving from directly below the station would have an incidence of 0 deg.
  alpha=atan(d/(pfd+alti/1000.0_wr))
  ! --------------------------------------------------------------------    .

  ! --------------------------------------------------------------------    .
  M3D(1,1)=cos(alpha)
  M3D(1,2)=-sin(alpha)*sin(a_baz)
  M3D(1,3)=-sin(alpha)*cos(a_baz)
  M3D(2,1)=sin(alpha)
  M3D(2,2)=cos(alpha)*sin(a_baz)
  M3D(2,3)=cos(alpha)*cos(a_baz)
  M3D(3,1)=0.0_wr
  M3D(3,2)=-cos(a_baz)
  M3D(3,3)=sin(a_baz)
  ! --------------------------------------------------------------------    .
  allocate(sacL(npts1),sacQ(npts1),sacT(npts1))

  sacL = M3D(1,1)*sacZ + M3D(1,2)*sacE + M3D(1,3)*sacN
  sacQ = M3D(2,1)*sacZ + M3D(2,2)*sacE + M3D(2,3)*sacN
  sacT = M3D(3,1)*sacZ + M3D(3,2)*sacE + M3D(3,3)*sacN

  ! --------------------------------------------------------------------    .
  ! ------- Écriture des données                                --------    .
  ! --------------------------------------------------------------------    .
  outfile=TRIM(adjustl(file1))//'.L'
  kcmpnm='HHL'
  call wbsac(outfile,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,      &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
    nzsec,nzmsec,nvhdr,norid,nevid,npts1,internal4,nwfid,nxsize,nysize,unused8,  &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm1,khole,ko,ka,kt0,kt1,&
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
    kinst,sacL)
  ! --------------------------------------------------------------------    .
  outfile=TRIM(adjustl(file2))//'.Q'
  kcmpnm='HHQ'
  call wbsac(outfile,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,      &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
    nzsec,nzmsec,nvhdr,norid,nevid,npts1,internal4,nwfid,nxsize,nysize,unused8,  &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm2,khole,ko,ka,kt0,kt1,&
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
    kinst,sacQ)
  ! --------------------------------------------------------------------    .
  outfile=TRIM(adjustl(file3))//'.T'
  kcmpnm='HHT'
  call wbsac(outfile,delta1,depmin,depmax,scale,odelta,b1,e1,o,a,internal1,      &
    t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
    resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
    user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
    internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
    unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
    nzsec,nzmsec,nvhdr,norid,nevid,npts1,internal4,nwfid,nxsize,nysize,unused8,  &
    iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
    imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
    unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm3,khole,ko,ka,kt0,kt1,&
    kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
    kinst,sacT)
  ! --------------------------------------------------------------------    .
  deallocate(sacZ,sacN,sacE,sacL,sacQ,sacT)
  ! --------------------------------------------------------------------    .
end program ZNE_2_LQT

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .

