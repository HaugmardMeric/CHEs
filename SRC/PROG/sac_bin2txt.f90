! programme principal sac_bin2txt
! ***********************************************************************    .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! -------                                                        --------    .
! ------- sac2txt                                                --------    .
! -------                                                        --------    .
! ------- septembre 2014                                         --------    .
! -------                                                        --------    .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr           --------    .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
!  This program is distributed for research purposes and in the hope         !
!  that it will be useful however without any warranties.                    !
!  Use it on your own risk. Please, report bugs/improvements/... .           !
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
!                                                                            !
!    transforme un fichier sac en .txt -> hodochrone (dis hypo, km)          !
!                                                                            !
program bin2txt
  ! --------------------------------------------------------------------    .
  use sac_i_o
  use modparam
  use typetemps
  use time
  ! -------                                                      --------   .
  implicit none
  ! -------                                                      --------   . pour SAC
  real(KIND=4), dimension(:), allocatable :: sacfile
  real(KIND=4) :: delta1, b1, e1
  integer(KIND=4) :: NN, ios, npts1
  character(LEN=112) :: file1,ofile
  ! ---------------------------------------------------------------------   . autres
  type(date_sec) :: d,ref,tzero
  real(KIND=wr) :: AmpliTot,AmpliMin,AmpliMax,AmpliNORM, tt
  real(KIND=wr) :: Disthypo
  integer(KIND=wi) :: jourj, j
  ! ---------------------------------------------------------------------   . lecture du fichier sac
  NN = IARGC()
  if ((NN < 2).or.(NN > 2)) THEN
    write(*,'(a)') "usage:  sac2xy sacfile ofile"
    write(*,'(a)') "        sacfile - input sac file"
    write(*,'(a)') "        ofile - output xy file"
    stop
  endif
  call GETARG(1, file1)
  open(UNIT=1,FILE=file1,STATUS='old',IOSTAT=ios)
  if (ios > 0) THEN
    write(*,*) "ERROR - Input file: '", TRIM(adjustl(file1)), "' does not exist ..."
    CLOSE(1)
    stop
  endif
  close(1)
  call GETARG(2,ofile)
  ! ---------------------------------------------------------------------   .
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
      write(*,*) "ERROR - File: '", TRIM(adjustl(file1)), "' appears to be of non-native &
      &byte-order or is not a SAC file."
      stop
    endif
  ! ---------------------------------------------------------------------   . lecture
  read(*,*)tzero,Disthypo,AmpliNORM
  ! ---------------------------------------------------------------------   . normalisation amplitude
  AmpliMin=1.e9_wr
  AmpliMax=-1.e9_wr
  do j=0,(int(npts1,wi)-1)
    if (AmpliMin.gt.real(sacfile(j+1),wr)) AmpliMin=real(sacfile(j+1),wr)
    if (AmpliMax.lt.real(sacfile(j+1),wr)) AmpliMax=real(sacfile(j+1),wr)
  enddo
  AmpliTot=max(abs(AmpliMin),abs(AmpliMax))
  ! ---------------------------------------------------------------------   . compute jour julien 
  d%date%year=int(nzyear,wi)
  d%date%hour=int(nzhour,wi)
  d%date%min=int(nzmin,wi)
  d%sec=real(nzsec,wr)+ real(b1,wr) + real(nzmsec,wr)/1000.0_wr
  ref%date%year=d%date%year
  ref%date%month=int(1,wi)
  jourj=int(nzjday,wi)
  ref%date%day=jourj
  ref%date%hour=int(0,wi)
  ref%date%min=int(0,wi)
  ref%sec=1.0_wr
  call basetime(ref)
  call jDATE(ref%date%jday,ref%date%year,ref%date%month,ref%date%day)
  call GDATE(ref%date%jday,ref%date%year,ref%date%month,ref%date%day)
  call basetime(ref)
  d%date%year=ref%date%year
  d%date%jday=ref%date%jday
  d%date%month=ref%date%month
  d%date%day=ref%date%day
  call basetime(d)
  ! ---------------------------------------------------------------------   . écriture
  open(UNIT=2,FILE=ofile,STATUS='REPLACE')
  ref=d
  ref%sec = ref%sec - real(delta1,wr)
  do j=0,(int(npts1,wi)-1)
    ref%sec = ref%sec + real(delta1,wr)
    call basetime(ref)
    call difftime(tt,ref,tzero)
    if (tt.gt.-30.0_wr) then
      write(2,'(2f18.6)') tt,Disthypo+real(sacfile(j+1),wr)*AmpliNORM/AmpliTot
    endif
  enddo
  close(2)
  ! ---------------------------------------------------------------------   .

end program bin2txt


! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .



