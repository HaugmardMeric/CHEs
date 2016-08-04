! programme principal sac_spectre
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! ------- spectre : lit les fichiers (phases list) lisible       --------   .
! ------- par le programme CHE et trace les spectres             --------   .
! ------- autour des arrivées de Pg et Sg                        --------   .
! -------                                                        --------   .
! ------- septembre/octobre 2015                                 --------   .
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
program spectre
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
  ! ---------------------------------------------------------------------   . fft
  integer :: i,j,iter,t,nn1,nn2
  double precision, allocatable :: spdata(:),temps(:),freq(:)
  double precision :: thedelta, amin, amax
  ! --------------------------------------------------------------------    .
  double precision :: plusoumoins=0.75
  logical :: plot=.false.
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
  ! --------------------------------------------------------------------    .
  thedelta=real(delta1)
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
  ! ONDES DIRECTES                                                          .
  ! --------------------------------------------------------------------    .
  ! T1 : temps arrivée de l'onde Pg
  ! T2 : incertitudes absolue sur l'onde Pg                                 ! -> calcul du delta après
  ! --------------------------------------------------------------------    .
  ! T3 : temps arrivée de l'onde Sg
  ! T4 : incertitudes absolue sur l'onde Sg
  ! ---------------------------------------------------------------------   .
  if (Adata%typeonde.eq.'G') then
    if(Adata%andS.eq.'S') then
      T1=real(alpha,4)
      T2=T1+real(Adata%sigP,4)
      T3=real(beta,4)
      T4=T3+real(Adata%sigS,4)
      if (.not.((T3.gt.(0.0_4+b1)).and.(T4.lt.(real(npts1,4)*delta1+b1)))) then
        T3=-12345.0_4
        T4=-12345.0_4
      endif
      if (.not.((T4.gt.(0.0_4+b1)).and.(T5.lt.(real(npts1,4)*delta1+b1)))) then
        T3=-12345.0_4
        T4=-12345.0_4
      endif
    else
      T1=real(alpha,4)
      T2=T1+real(Adata%sigP,4)
    endif
    if (.not.((T1.gt.(0.0_4+b1)).and.(T2.lt.(real(npts1,4)*delta1+b1)))) then
      T1=-12345.0_4
      T2=-12345.0_4
      T3=-12345.0_4
      T4=-12345.0_4
    endif
    if (.not.((T2.gt.(0.0_4+b1)).and.(T3.lt.(real(npts1,4)*delta1+b1)))) then
      T1=-12345.0_4
      T2=-12345.0_4
      T3=-12345.0_4
      T4=-12345.0_4
    endif
    ! ------------------------------------------------------------------    .
    ! ONDES REFRACTEES                                                      .
    ! ------------------------------------------------------------------    .
    ! T5 : temps arrivée de l'onde Pn
    ! T6 : incertitudes absolue sur l'onde Pn
    ! ------------------------------------------------------------------    .
    ! T7 : temps arrivée de l'onde Sn
    ! T8 : incertitudes absolue sur l'onde Sn
    ! -------------------------------------------------------------------   .
  elseif (Adata%typeonde.eq.'N') then                                       !
    if(Adata%andS.eq.'S') then
      T5=real(alpha,4)
      T6=T5+real(Adata%sigP,4)
      T7=real(beta,4)
      T8=T7+real(Adata%sigS,4)
      if (.not.((T7.gt.(0.0_4+b1)).and.(T8.lt.(real(npts1,4)*delta1+b1)))) then
        T7=-12345.0_4
        T8=-12345.0_4
      endif
      if (.not.((T8.gt.(0.0_4+b1)).and.(T9.lt.(real(npts1,4)*delta1+b1)))) then
        T8=-12345.0_4
        T7=-12345.0_4
      endif
    else
      T5=real(alpha,4)
      T6=T5+real(Adata%sigP,4)
    endif
    if (.not.((T5.gt.(0.0_4+b1)).and.(T6.lt.(real(npts1,4)*delta1+b1)))) then
      T5=-12345.0_4
      T6=-12345.0_4
      T8=-12345.0_4
      T7=-12345.0_4
    endif
    if (.not.((T6.gt.(0.0_4+b1)).and.(T7.lt.(real(npts1,4)*delta1+b1)))) then
      T6=-12345.0_4
      T5=-12345.0_4
      T8=-12345.0_4
      T7=-12345.0_4
    endif
  ! --------------------------------------------------------------------    .
  else
    write(*,*)'problème : onde ni directe ni réfractée ... ? ', Adata%typeonde
    stop
  endif
  ! --------------------------------------------------------------------    .
  if (nvhdr /= 6) then
    write(*,*) "ERROR 2 - File: '", TRIM(adjustl(file1)), "' appears to be of non-native byte-order or is not a SAC file."
    stop
  endif
  ! --------------------------------------------------------------------    .


  t=int(plusoumoins/thedelta)

  ! --------------------------------------------------------------------    .
  ! ONDE Pg :
  ! --------------------------------------------------------------------    .
  do i=1,npts1
    if (real(real(i)*thedelta).lt.(real(T1))) then
      iter=i
    endif
  enddo
  ! --------------------------------------------------------------------    .
  nn1=min(iter+t,npts1)-max(iter-t,1)
  nn2=2**max0(int(dlog(dble(nn1))/dlog(2.0d0))+1,3)                         ! prochaine puissance de 2 (fft)
  ! --------------------------------------------------------------------    .
  allocate(spdata(2*nn2),temps(2*nn2),freq(2*nn2))
  ! --------------------------------------------------------------------    .
  freq(1)=0.0d0
  freq(2)=0.0d0
  do i=3,nn2-1,2
    freq(i)=real((i+1)/2-1)/real(nn2*thedelta)
    freq(i+1)=freq(i)
  enddo
  freq(nn2+1)=1.0/real(2*thedelta)
  freq(nn2+2)=freq(nn2+1)
  do i=nn2+2,2*nn2-2,2
    freq(i)=real((-(i+1)+nn2)/2-1)/real(nn2*thedelta)
    freq(i+1)=freq(i)
  enddo
  freq(2*nn2-1)=-1.0/real(nn2*thedelta)
  freq(2*nn2)=freq(2*nn2-1)
  ! --------------------------------------------------------------------    .
  spdata(:)=0.0d0
  temps(:)=0.0d0
  ! --------------------------------------------------------------------    .
  j=max(iter-t,1)
  do i=1,2*nn2-2,2
    if(j<min(iter+t,npts1)) then
      spdata(i)=sacfile(j)
    else
      spdata(i)=0.0d0
    endif
    j=j+1
  enddo
  ! --------------------------------------------------------------------    .
  do i=3,2*nn2-2,2
    temps(i)=temps(i-1)+thedelta
    temps(i+1)=temps(i-1)+thedelta
  enddo
  ! --------------------------------------------------------------------    .
  if (plot) then
  open(100, FILE = 'data.xy1',status='replace')
  do i=1,2*nn1,2
    write(100,*)temps(i),spdata(i)
  enddo
  close(100)
  endif
  ! --------------------------------------------------------------------    .
  amin=999.999d0
  amax=-999.999d0
  do i=3,nn2-2,2
    if (log10(sqrt(spdata(i)**2.0d0+spdata(i+1)**2.0d0))<amin) amin=log10(sqrt(spdata(i)**2.0d0+spdata(i+1)**2.0d0))
    if (log10(sqrt(spdata(i)**2.0d0+spdata(i+1)**2.0d0))>amax) amax=log10(sqrt(spdata(i)**2.0d0+spdata(i+1)**2.0d0))
  enddo
  ! --------------------------------------------------------------------    .
  call four1(spdata,nn2,1)
  open(102, FILE = 'OUTPUT/GMT/spectre_'//TRIM(adjustl(kstnm))//'_'//TRIM(adjustl(kcmpnm))//'.ph',status='replace')
  open(101, FILE = 'OUTPUT/GMT/spectre_'//TRIM(adjustl(kstnm))//'_'//TRIM(adjustl(kcmpnm))//'.am',status='replace')
  do i=3,nn2-2,2
    write(101,*)freq(i),log10(sqrt(spdata(i)**2.0d0+spdata(i+1)**2.0d0))/(amax-amin)*1000.d0-amin+0.1d0
    write(102,*)freq(i),atan2(spdata(i+1),spdata(i))
  enddo
  close(101)
  close(102)
  ! --------------------------------------------------------------------    .
  if (plot) then
  call four1(spdata,nn2,-1)
  open(103, FILE = 'data.xy2',status='replace')
  do i=1,2*nn1,2
    spdata(i)=spdata(i)/real(nn2)
    write(103,*)temps(i),spdata(i)
  enddo
  close(103)
  endif
  ! --------------------------------------------------------------------    .
  deallocate(spdata,temps,freq)

  ! --------------------------------------------------------------------    . format de lecture des données
  1000 format(a4,1x,a1,1x,i1,1x,5i2.2,f6.3,5x,f7.3,1x,a1,i1.1,3x,f6.3,4x,f6.3)
  ! --------------------------------------------------------------------    .


  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .
  CONTAINS
  ! --------------------------------------------------------------------    .
  ! --------------------------------------------------------------------    .


  subroutine four1(adata,nn,isign)
    ! -------------------------------------------------------------------   .
    implicit none
    ! -------------------------------------------------------------------   .
    ! FFT routine :
    ! Press W. H., Teukolsky S.A., Vetterling W.T.
    ! Numerical Recipes in Fortran 77 -- The Art of Scientific Computing --
    ! Second Edition - Volume 1 of Fortran Numerical Recipes1
    ! -------------------------------------------------------------------   .
    integer, intent (in) :: isign,nn
    double precision, intent (inout) :: adata(2*nn)
    ! -------                                                   --------    .
    integer i,istep,j,m,mmax,n
    double precision tempi,tempr,theta,wi,wpi,wpr,wr,wtemp
    ! -------------------------------------------------------------------   .
    n=2*nn
    j=1
    do i=1,n,2                                                              ! This is the bit-reversal section of the routine.
      if(j.gt.i)then
        tempr=adata(j)                                                      ! Exchange the two complex numbers.
        tempi=adata(j+1)
        adata(j)=adata(i)
        adata(j+1)=adata(i+1)
        adata(i)=tempr
        adata(i+1)=tempi
      endif
      m=n/2
      10001 if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        goto 10001
      endif
      j=j+m
    enddo
    ! -------------------------------------------------------------------   .
    mmax=2                                                                  ! Here begins the Danielson-Lanczos section of the routine.
    ! -------                                                   --------    .
    10002 if (n.gt.mmax) then                                               ! Outer loop executed log2 nn times.
      istep=2*mmax
      theta=6.28318530717959d0/real(isign*mmax)                             ! Initialize for the trigonometric recurrence
      wpr=-2.d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      do m=1,mmax,2                                                         ! Here are the two nested inner loops.
        do i=m,n,istep
          j=i+mmax                                                          ! This is the Danielson-Lanczos formula:
          tempr=sngl(wr)*adata(j)-sngl(wi)*adata(j+1)
          tempi=sngl(wr)*adata(j+1)+sngl(wi)*adata(j)
          adata(j)=adata(i)-tempr
          adata(j+1)=adata(i+1)-tempi
          adata(i)=adata(i)+tempr
          adata(i+1)=adata(i+1)+tempi
        enddo
        wtemp=wr                                                            ! Trigonometric recurrence.
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
      enddo
      mmax=istep
      goto 10002                                                            ! Not yet done.
    endif                                                                   ! All done.
    ! -------                                                   --------    .
  end subroutine four1
  ! --------------------------------------------------------------------    .


end program spectre

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .

