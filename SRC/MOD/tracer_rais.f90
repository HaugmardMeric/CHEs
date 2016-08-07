! octobre 2015
! tracer de rais 1D, couches tabulaires et homogènes
! module non utilise pour la recherche du modèle de Terre dans CHE
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE ray_tracing

    use modparam

    implicit none

    private

    public  :: tracerays

    ! -----------------------------------------------------------------    .
    real (kind=wr), parameter :: eps = 0.000001_wr                         ! le chouïa
    ! -----------------------------------------------------------------    .

! ---------------------------------------------------------------------    .
CONTAINS
! ---------------------------------------------------------------------    .

  subroutine tracerays(distancepi,distancehypo,dcritiqueH,pfdseisme,altidudesta,lon, &
                         lat,pdfmoho,tdirectP,trefP,tdirectS,trefS,modele)
    ! -------                                                  --------    . mh
    ! trace les rais sismiques des ondes Pg,Sg,Pn,Sn pour un couple station-séisme,
    ! avec calcul du temps de parcours et un modèle de Terre donnée (modele={'A';'S';'F';'C'})
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) :: distancepi,pfdseisme                     ! pour un couple station-séisme
    real (kind=wr), intent(in) :: altidudesta,lon,lat                      ! pour un couple station-séisme
    real(KIND=wr), intent (out) :: pdfmoho
    real (kind=wr), intent(out) :: tdirectP,trefP,tdirectS,trefS           ! temps parcours ondes directes et réfractées
    real (kind=wr), intent(out) :: distancehypo,dcritiqueH
    character (LEN=1), intent(in), optional :: modele
    ! -----------------------------------------------------------------    .
    integer (kind=wi) :: nbl                                               ! nombre couche du modèle
    integer (kind=wi) :: nc                                                ! couche dans laquelle est le séisme
    ! -------                                                  --------    .
    integer (kind=wi) :: i
    ! -------                                                  --------    .
    real (kind=wr), allocatable :: vP(:),vS(:)                             ! vitesse P et S dans chaque couche
    real (kind=wr), allocatable :: h(:)                                    ! profondeur cumulés du toit la ieme couche
    real (kind=wr), allocatable :: hpuiss(:)                               ! puissance la ieme couche
    real (kind=wr), allocatable :: X(:),X2(:,:)                            ! distance, projetée horizontalement, parcourue dans chaque couche
    logical, parameter :: plotfigures=.false.
    ! -----------------------------------------------------------------    . choix d'un modèle de terre :
    if(present(modele)) then
      ! -------                                                --------    .
      if (modele=='A') then
        ! call ModTerrArroucau(nbl,vP,vS,h,lon,lat,forcelettre='a')
        call ModTerrArroucau(nbl,vP,vS,h,lon,lat,pdfmoho)
      elseif (modele=='S') then
        call ModTerrSiHex_Haslach(nbl,vP,vS,h,pdfmoho)
      elseif (modele=='F') then
        call ModTerr_fun(nbl,vP,vS,h,pdfmoho)
      elseif (modele=='C') then
        call ModTerr_fun(nbl,vP,vS,h,pdfmoho)
      else
        write(*,*)'problème dans tracerays2, le modèle de terre : ',modele,' n''existe pas '
        stop
      endif
      ! -------                                                --------    . par défaut
    else
      call ModTerrSiHex_Haslach(nbl,vP,vS,h,pdfmoho)
    endif
    ! -----------------------------------------------------------------    . VÉRIF modèle de Terre
    do i=2,nbl
      if(h(i).le.h(i-1)) then
        write(*,*)'problème dans tracerays2 : ! ordre profondeur ', h
        stop
      endif
    enddo
    ! -------                                                  --------    .
    if (altidudesta.gt.h(2)) then
      write(*,*)'problème dans tracerays2 : ! altidude de la station non adapté au modèle de terre'
      stop
    endif
    ! -------                                                  --------    . dans quelle couche est le séisme ?
    do i=1,nbl
      if(pfdseisme.gt.h(i)) nc=i
    enddo
    if (nc.ge.nbl) then
      ! write(*,*)'problème dans tracerays2 : le séisme est trop profond'
      tdirectP=-1.e9_wr
      trefP=1.e9_wr
      tdirectS=1.e9_wr
      trefS=1.e9_wr
      distancehypo=0.0_wr
      dcritiqueH=0.0_wr
    else
      ! -------                                                --------    .
      allocate(hpuiss(nbl),X(nbl))
      ! ---------------------------------------------------------------    . puissance de chaque couche
      hpuiss=0.0_wr
      do i=1,nc
        if(i==nc) then
          hpuiss(i)=pfdseisme-h(i)
        else
          hpuiss(i)=h(i+1)-h(i)
        endif
      enddo
      hpuiss(1)=hpuiss(1)-altidudesta
      h(1)=altidudesta
      ! -------                                                --------    .
      distancehypo=0.0_wr
      ! ---------------------------------------------------------------    . pour onde Pg
      call traceondeGKIM(distancepi,pfdseisme,altidudesta,nbl,nc,vP,hpuiss,X,distancehypo,tdirectP)
      if ((tdirectP.gt.0.0_wr).and.(plotfigures)) then
        call printAray(int(distancepi*100.0_wr),nbl,nc,h,X,X2,pfdseisme,distancepi,altidudesta,lettre='Pg')
      endif
      ! -------                                                --------    . pour onde Pn
      call traceondeN(distancepi,pfdseisme,altidudesta,nbl,nc,vP,h,X,X2,distancehypo,dcritiqueH,trefP)
      if ((trefP.gt.0.0_wr).and.(plotfigures)) then
        call printAray(int(distancepi*100.0_wr),nbl,nc,h,X,X2,pfdseisme,distancepi,altidudesta,lettre='Pn')
      endif
      ! -------                                                --------    . pour onde Sg
      call traceondeGKIM(distancepi,pfdseisme,altidudesta,nbl,nc,vS,hpuiss,X,distancehypo,tdirectS)
      if ((tdirectS.gt.0.0_wr).and.(plotfigures)) then
        call printAray(int(distancepi*100.0_wr),nbl,nc,h,X,X2,pfdseisme,distancepi,altidudesta,lettre='Sg')
      endif
      ! -------                                                --------    .  pour onde Sn
      call traceondeN(distancepi,pfdseisme,altidudesta,nbl,nc,vS,h,X,X2,distancehypo,dcritiqueH,trefS)
      if ((trefS.gt.0.0_wr).and.(plotfigures)) then
        call printAray(int(distancepi*100.0_wr),nbl,nc,h,X,X2,pfdseisme,distancepi,altidudesta,lettre='Sn')
      endif
      ! -------                                                --------    .
      deallocate(hpuiss,X,X2)
    endif
    ! -----------------------------------------------------------------    .
    deallocate(vP,vS,h)
    ! -----------------------------------------------------------------    .
  end subroutine tracerays

! ---------------------------------------------------------------------    .

  subroutine traceondeDichotomie(distancepi,nbl,nc,v,hpuiss,X,theta)
    ! -------                                                  --------    . mh
    ! trace onde directe, recherche de theta par dichotomie 
    ! (on retire ici un quart seulement de l'espace, pas la moitié)
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) :: distancepi                               ! pour un couple station-séisme
    real (kind=wr), intent(in) :: hpuiss(nbl)
    integer (kind=wi), intent(in) :: nbl,nc
    real (kind=wr), intent(in) :: v(nbl)
    ! -------                                                  --------    .
    real (kind=wr), intent(out) :: X(nbl)
    real (kind=wr), intent(out) :: theta
    ! -----------------------------------------------------------------    .
    integer (kind=wi) :: i
    ! -------                                                  --------    .
    real (kind=wr) :: thetamin,thetamax                                    ! takeoff angle (clockwise from the upward vertical direction)
    real (kind=wr) :: sommeX,sommeXmin,sommeXmax                           ! somme des Xi
    ! -----------------------------------------------------------------    .
    thetamin=0.0_wr*pi/180.0_wr+eps/10.0_wr
    thetamax=90.0_wr*pi/180.0_wr-eps/10.0_wr
    ! -----------------------------------------------------------------    .
    i=0
    do while ((abs(sommeX-distancepi).gt.eps).and.(i.lt.20000).and.(thetamin.ne.thetamax))
      i=i+1
      theta=(thetamin+thetamax)/2._wr
      call CalcX(nbl,nc,hpuiss,(thetamin+theta)/2._wr,X,v,sommeXmin)
      call CalcX(nbl,nc,hpuiss,(theta+thetamax)/2._wr,X,v,sommeXmax)
      call CalcX(nbl,nc,hpuiss,theta,X,v,sommeX)
      if(abs(sommeXmin-distancepi).le.abs(sommeXmax-distancepi)) then
        thetamax=(thetamax+theta)/2._wr
      else
        thetamin=(theta+thetamin)/2._wr
      endif
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine traceondeDichotomie

! ---------------------------------------------------------------------    .

  subroutine traceondeGKIM(distancepi,pfdseisme,altidudesta,nbl,nc,v,hpuiss,X,distancehypo,temps)
    ! -------                                                  --------    . mh
    ! ray tracing : d'après Kim and Baag (2002) :
    ! Rapid and Accurate Two-Point Ray Tracing Based on
    ! a Quadratic Equation of Takeoff Angle in Layered Media with Constant (BSSA)
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) :: distancepi,pfdseisme,altidudesta         ! pour un couple station-séisme
    real (kind=wr), intent(in) :: hpuiss(nbl)
    integer (kind=wi), intent(in) :: nbl,nc
    real (kind=wr), intent(inout) :: v(nbl)
    ! -------                                                  --------    .
    real (kind=wr), intent(out) :: X(nbl),distancehypo
    real (kind=wr), intent(out) :: temps                                   ! temps de parcours des ondes
    ! -----------------------------------------------------------------    .
    integer (kind=wi) :: i,j
    ! -------                                                  --------    .
    real (kind=wr) :: vx(nbl),vt(nbl),vxtot                                ! vérif ...
    real (kind=wr) :: thetainit,theta                                      ! takeoff angle (clockwise from the upward vertical direction)
    real (kind=wr) :: deltatheta1,deltatheta2                              ! incrément de theta, angle (2 solutions équation 2nd degré)
    real (kind=wr) :: sommeX                                               ! somme des Xi
    ! -----------------------------------------------------------------    .
    sommeX=99999.999_wr
    ! -----------------------------------------------------------------    .
    call initialvalue1(nbl,nc,v,hpuiss,thetainit,distancepi,pfdseisme,altidudesta)
    ! -----------------------------------------------------------------    .
    ! test pour plusieurs valeurs initales si besoin, avec incrément tous les dégrés, jusqu'à convergence
    j=-1
    do while ((abs(sommeX-distancepi).gt.eps).and.(j.lt.180))
      ! ---------------------------------------------------------------    .
      j=j+1
      theta=thetainit+real(j,wr)/0.5_wr*pi/180.0_wr
      ! ---------------------------------------------------------------    .
      i=1
      call CalcX(nbl,nc,hpuiss,theta,X,v,sommeX)
      do while ((abs(sommeX-distancepi).gt.eps).and.(i.lt.20))
        i=i+1
        call CalcDeltaTheta(nbl,v,hpuiss,theta,X,deltatheta1,deltatheta2,distancepi,nc)
        theta=theta+deltatheta2
        do while((theta*180.0_wr/pi).ge.89.9_wr)
          theta=theta-1.0_wr
        enddo
        call CalcX(nbl,nc,hpuiss,theta,X,v,sommeX)
        !write(*,*)i,'DELTA : -- ',sommeX-distancepi,deltatheta1*180.0_wr/pi,theta*180.0_wr/pi
      enddo
      ! ---------------------------------------------------------------    .
      if (isnan(sommeX))sommeX=99999.999_wr
      ! ---------------------------------------------------------------    .
      if (j.gt.179) then
        call traceondeDichotomie(distancepi,nbl,nc,v,hpuiss,X,theta)       ! trace onde directe, recherche de theta par dichotomie
        call CalcX(nbl,nc,hpuiss,theta,X,v,sommeX)
        j=200
      endif
    enddo
    ! -----------------------------------------------------------------    .VERIF
    temps=0.0_wr
    vxtot=0.0_wr
    distancehypo=0.0_wr
    vt(nc)=theta
    do i=nc,1,-1
      if ((i-1).gt.0) vt(i-1)=asin(sin(vt(i))*v(i-1)/v(i))
      vx(i)=tan(vt(i))*hpuiss(i)
      vxtot=vxtot+vx(i)
      temps=temps+(hpuiss(i)/cos(vt(i)))/v(i)
      distancehypo=distancehypo+hpuiss(i)/cos(vt(i))
    enddo
    ! -----------------------------------------------------------------    .pb cos(x)
    if (isnan(temps))then
      v=0.0_wr
      X=0.0_wr
      distancehypo=0.0_wr
      temps=0.0_wr
    endif
    ! -----------------------------------------------------------------    .
    ! write(*,'(a,f15.1,2f15.9,f15.1,f20.10)')'K',distancepi,vxtot-distancepi,sommeX-distancepi,temps,theta*180.0_wr/pi
    ! -----------------------------------------------------------------    .
  end subroutine traceondeGKIM

! ---------------------------------------------------------------------    .

  subroutine traceondeN(distancepi,pfdseisme,altidudesta,nbl,nc,v,h,X,Xplot,distancehypo,dcritiqueH,temps)
    ! -------                                                  --------    . mh
    ! ray tracing : onde réfracté au moho (toujours la dernière couche du modèle de terre !)
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) :: distancepi,pfdseisme,altidudesta         ! pour un couple station-séisme
    real (kind=wr), intent(in) :: h(nbl)
    integer (kind=wi), intent(in) :: nbl,nc
    real (kind=wr), intent(inout) :: v(nbl)
    ! -------                                                  --------    .
    real (kind=wr), intent(inout), allocatable :: Xplot(:,:)
    real (kind=wr), intent(out) :: X(nbl),distancehypo,dcritiqueH
    real (kind=wr), intent(out) :: temps                                   ! temps de parcours des ondes
    ! -----------------------------------------------------------------    .
    integer (kind=wi) :: i,j
    integer (kind=wi) :: nm                                                ! couche de la réfraction
    ! -------                                                  --------    .
    real (kind=wr) :: theta(nbl),hpuiss(nbl),hequi(nbl),xcumule
    ! -----------------------------------------------------------------    .
    nm=nbl
    ! -----------------------------------------------------------------    . modèle équivalent pour la première réfracté
    ! ne prend pas encore en compte : altidudesta
    ! -----------------------------------------------------------------    .
    hequi=0.0_wr
    hpuiss=0.0_wr
    do i=1,nbl-1
      hpuiss(i)=h(i+1)-h(i)
    enddo
    ! -------                                                  --------    .
    distancehypo=0.0_wr
    ! -------                                                  --------    .
    do i=1,nm-1
      if (i.lt.nc) then
        hequi(i)=hpuiss(i)
      elseif(i==nc) then
        hequi(i)=hpuiss(i)+(h(nc+1)-pfdseisme)
      else ! (i.gt.nc)
        hequi(i)=2.0_wr*hpuiss(i)
      endif
    enddo
    ! -----------------------------------------------------------------    .
    theta(nm)=pi/2.0_wr
    do i=nm-1,1,-1
      theta(i)=asin(sin(theta(i+1))*v(i)/v(i+1))
    enddo
    ! -----------------------------------------------------------------    .
    xcumule=0.0_wr
    dcritiqueH=0.0_wr
    X(nm)=0.0_wr
    do i=1,nm-1
      dcritiqueH=dcritiqueH+hequi(i)/cos(theta(i))
      X(i)=tan(theta(i))*hequi(i)
      X(nm)= X(nm)+X(i)
    enddo
    X(nm)= distancepi-X(nm)
    ! -----------------------------------------------------------------    .
    if (allocated(Xplot)) deallocate(Xplot)
    allocate(Xplot(2*int(nbl,4)-nc,2))
    Xplot=0.0_wr
    j=1
    Xplot(j,1)=tan(theta(nc))*(h(nc+1)-pfdseisme)
    Xplot(j,2)=h(nc+1)
    do i=nc+1,nm-1
      j=j+1
      Xplot(j,1)=Xplot(j-1,1)+tan(theta(i))*hpuiss(i)
      Xplot(j,2)=h(i+1)
    enddo
    j=j+1
    Xplot(j,1)=Xplot(j-1,1)+X(nm)
    Xplot(j,2)=Xplot(j-1,2)
    do i=nm-1,1,-1
      j=j+1
      Xplot(j,1)=Xplot(j-1,1)+tan(theta(i))*hpuiss(i)
      Xplot(j,2)=h(i)
    enddo
    ! -----------------------------------------------------------------    .
    distancehypo=dcritiqueH+X(nm)
    temps=0.0_wr
    do i=1,nm-1
      temps=temps+(hequi(i)/cos(theta(i)))/v(i)
    enddo
    if (distancehypo.ge.dcritiqueH) then
        temps=temps+X(nm)/v(nm)
    else
      temps=0.0_wr
      Xplot=0.0_wr
      X=0.0_wr
    endif
    ! -----------------------------------------------------------------    .pb
    if (isnan(temps))then
      temps=0.0_wr
      Xplot=0.0_wr
      X=0.0_wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine traceondeN

! ---------------------------------------------------------------------    .

  subroutine initialvalue1(nbl,nc,vP,h,theta,distancepi,pfdseisme,altidudesta)
    ! -------                                                  --------    . mh
    ! inital value 1, primary estimation : kinematic properties of rays
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) ::vP(nbl),h(nbl),distancepi,pfdseisme,altidudesta
    integer (KIND=wi), intent(in) :: nbl,nc
    real (kind=wr), intent(out) :: theta
    ! -------                                                  --------    .
    integer (KIND=wi) :: i
    real (kind=wr) :: thetaI,thetaH,thetaM
    real (kind=wr) :: vmoy,vmax,S,si(nc)
    ! -----------------------------------------------------------------    .
    thetaH=atan(distancepi/abs(pfdseisme-altidudesta))
    ! -------                                                  --------    . vitesse maximale
    vmax=-99.99_wr
    do i=1,nc
      if(vmax.le.vP(i)) vmax=vP(i)
    enddo
    !if(vmax.ne.vP(nc)) write(*,*)'problème dans initialvalue : vitesse on max',vmax,vP(nc)
    ! -------                                                  --------    .
    !vmax=vP(nbl)
    ! -------                                                  --------    .
    vmoy=0.0_wr
    S=0.0_wr
    si=0.0_wr
    ! -------                                                  --------    .
    do i=1,nc
      thetaI=asin(vP(i)/vmax*sin(thetaH))
      si(i)=h(i)/cos(thetaI)
      S=S+si(i)
      vmoy=vmoy+vP(i)*si(i)
    enddo
    vmoy=vmoy/S
    ! -----------------------------------------------------------------    .
    if((vmax/vmoy*sin(thetaH)).le.1.0_wr) then
      thetaM=asin(vmax/vmoy*sin(thetaH))
    else
      thetaM=thetaH
    endif
    ! -------                                                  --------    .
    theta=(thetaH+thetaM)/2.0_wr
    ! -----------------------------------------------------------------    .
    if (isnan(theta))theta=thetaH*1.1_wr
    ! -------                                                  --------    .
    do while((theta*180.0_wr/pi).ge.90.0_wr)
      theta=theta-1.5_wr*pi/180.0_wr
    enddo
    !write(*,*)'initial value 1 : ',thetainit*180.0_wr/pi
    ! -----------------------------------------------------------------    .
  end subroutine initialvalue1

! ---------------------------------------------------------------------    .

  subroutine CalcDeltaTheta(nbl,vP,h,theta,X,deltatheta1,deltatheta2,distancepi,nc)
    ! -------                                                  --------    . mh
    ! CalcDeltaTheta
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) ::vP(nbl),h(nbl),distancepi,theta,X(nbl)
    integer (KIND=wi), intent(in) :: nbl,nc
    real (kind=wr), intent(out) :: deltatheta1,deltatheta2
    ! -------                                                  --------    .
    integer (KIND=wi) :: i
    real (kind=wr) :: A,B,C
    real (kind=wr) :: lambda(nbl)
    real (kind=wr) :: racine
    real (kind=wr) :: sommeXi
    ! -----------------------------------------------------------------    .
    A=0.0_wr
    B=0.0_wr
    sommeXi=0.0_wr
    ! -------                                                  --------    .
    do i=1,nc
      ! -------                                                --------    .
      lambda(i)=vP(i)/vP(nc)
      racine=sqrt(1.0_wr-(lambda(i)**2.0_wr)*(sin(theta)**2.0_wr))
      ! -------                                                --------    .
      A=A+0.5_wr*h(i)*(3.0_wr*(lambda(i)**3.0_wr)*(cos(theta)**2.0_wr)*sin(theta)/(racine**5.0_wr) &
        - lambda(i)*sin(theta)/(racine**3.0_wr))
      ! -------                                                --------    .
      B=B+h(i)*lambda(i)*cos(theta)/(racine**3.0_wr)
      ! -------                                                --------    .
      sommeXi=sommeXi+X(i)
      C=sommeXi-distancepi
      ! -------                                                --------    .
    enddo
    ! -------                                                  --------    .
    if((B**2.0_wr-4.0_wr*A*C).gt.0.0_wr) then
      deltatheta1=(-B-sqrt(B**2.0_wr-4.0_wr*A*C))/(2.0_wr*A)
      deltatheta2=(-B+sqrt(B**2.0_wr-4.0_wr*A*C))/(2.0_wr*A)
    else
      deltatheta1=0.0_wr
      deltatheta2=0.0_wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine CalcDeltaTheta

! ---------------------------------------------------------------------    .

  subroutine CalcX(nbl,nc,h,theta,X,v,sommeX)
    ! -------                                                  --------    . mh
    ! CalcX
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) ::h(nbl),v(nbl),theta
    integer (KIND=wi), intent(in) :: nbl,nc
    real (kind=wr), intent(out) :: X(nbl),sommeX
    ! -------                                                  --------    .
    integer (KIND=wi) :: i
    real (kind=wr) :: lambda(nbl),racine                                   !
    ! -----------------------------------------------------------------    .
    sommeX=0.0_wr
    X=0.0_wr
    do i=1,nc
      lambda(i)=v(i)/v(nc)
      racine=sqrt(1.0_wr-(lambda(i)**2.0_wr)*(sin(theta)**2.0_wr))
      X(i)=h(i)*lambda(i)*sin(theta)/racine
      sommeX=sommeX+X(i)
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine CalcX

! ---------------------------------------------------------------------    .

  subroutine printAray(k,nbl,nc,h,X,Xplot,pfdseisme,distancepi,altidudesta,lettre)
    ! -------                                                  --------    . mh
    ! printAray
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(in) :: h(nbl),X(nbl),pfdseisme,distancepi,altidudesta
    integer (KIND=wi), intent(in) :: k,nbl,nc
    character (LEN=2), intent(in), optional :: lettre
    real (kind=wr), intent(in), allocatable :: Xplot(:,:)
    ! -------                                                  --------    .
    real (kind=wr) :: X2(nbl)
    integer (KIND=wi) :: i,ok
    character (LEN=5) :: numberchaine
    ! -----------------------------------------------------------------    .
    write(numberchaine(1:5),'(i5.2)')k
    ! -----------------------------------------------------------------    .
    if(present(lettre)) then
      open(105,FILE='ray-'//lettre//'-'//trim(adjustl(numberchaine))//'.d',status='replace',iostat = ok)
    else
      open(105,FILE='ray-'//trim(adjustl(numberchaine))//'.d',status='replace',iostat = ok)
    endif
    if (ok .ne. 0) then
      write(*,*)'problème dans printAray : le fichier n''existe pas '
      stop
    endif
    ! -------                                                  --------    .
    write(105,'(2f15.7)')0.0_wr,pfdseisme
    ! -------                                                  --------    .
    if(present(lettre)) then
    if(lettre(2:)=='n') then
      ! -------                                                --------    .
      do i=1,2*int(nbl-nc,4)
        write(105,*)Xplot(i,1),Xplot(i,2)
      enddo
      ! -------                                                --------    .
      elseif(lettre(2:)=='g') then
        ! -------                                              --------    .
        X2(nbl)=X(nbl)
        do i=nbl-1,1,-1
          X2(i)=X2(i+1)+X(i)
        enddo
        do i=nc,1,-1
          write(105,*)x2(i),h(i)
        enddo
        ! -------                                              --------    .
      else
        write(*,*)'problème dans printAray : onde ni "g" ni "n" -> ',lettre
        stop
       ! -------                                              --------    .
      endif
    else
      ! -------                                                --------    .
      X2(nbl)=X(nbl)
      do i=nbl-1,1,-1
        X2(i)=X2(i+1)+X(i)
      enddo
      do i=nc,1,-1
        write(105,*)x2(i),h(i)
      enddo
      ! -------                                                --------    .
    endif
    ! -------                                                  --------    .
    close(105)
    ! -----------------------------------------------------------------    .
    !if(k==0) then
      open(99, FILE = 'mod.d',status='replace',iostat = ok)
      do i=1,nbl
        write(99,*)">"
        write(99,*)-distancepi-10.0_wr,h(i)
        write(99,*)distancepi+10.0_wr,h(i)
      enddo
      close(99)
      ! ---------------------------------------------------------------    .
      open(101, FILE = 'plot.sh',status='replace',iostat = ok)
      write(101,'(a)')"gmtset LABEL_FONT_SIZE 15"
      write(101,'(a)')'gmtset HEADER_FONT_SIZE 15'
      write(101,'(a)')"gmtset ANNOT_FONT_PRIMARY Times-Roman"
      write(101,'(a)')"gmtset ANNOT_FONT_SECONDARY Times-Roman"
      write(101,'(a)')"gmtset PAPER_MEDIA A3"
      write(101,'(a)')"gmtset TIME_LANGUAGE FR"
      write(101,'(a)')"gmtset CHAR_ENCODING ISOLatin1+"
      ! ---------------------------------------------------------------    .
      write(101,'(a)')"geoproj=-JX12i/-1.2i"
      write(101,'(a)')"geozone=-R-55/550/-7/77"
      write(101,'(2a)')"psxy $geozone $geoproj mod.d -Wthinnest,green -Ba100f50/a10f2NsWeg100  -m -K > file.ps"
      write(101,*)"echo '",0.0_wr,pfdseisme,"' | psxy $geozone $geoproj -Sa0.25 -Gorange -Wthinnest,red -O -K >> file.ps"
      write(101,*)"echo '",distancepi,altidudesta,"' | psxy $geozone $geoproj -St0.25 -Gblue -Wthinnest,blue -O -K -N >> file.ps"
      write(101,'(2a)')"psxy $geozone $geoproj ray-Pg-*.d -Wthinnest,blue,-- -O -K >> file.ps"
      write(101,'(2a)')"psxy $geozone $geoproj ray-Sg-*.d -Wthinnest,red,-- -O -K >> file.ps"
      write(101,'(2a)')"psxy $geozone $geoproj ray-Pn-*.d -Wthinnest,blue,-. -O -K >> file.ps"
      write(101,'(2a)')"psxy $geozone $geoproj ray-Sn-*.d -Wthinnest,red,-. -O >> file.ps"
      write(101,*)"ps2raster file.ps -Tf -A -P "
      ! ---------------------------------------------------------------    .
      close(101)
    !endif
    ! -----------------------------------------------------------------    .
  end subroutine printAray

! ---------------------------------------------------------------------    .

  subroutine ModTerrArroucau(nbl,vP,vS,h,lon,lat,pdfmoho,forcelettre)
    ! -------                                                  --------    . mh
    ! Modèles de Terre Arroucau, phD 2006
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(inout), allocatable :: h(:),vP(:),vS(:)
    integer (KIND=wi), intent(out) :: nbl
    real(KIND=wr), intent (out) :: pdfmoho
    ! -------                                                  --------    .
    real (kind=wr), intent(in) :: lon,lat
    character (LEN=1), intent(in), optional :: forcelettre
    ! -------                                                  --------    .
    character (LEN=1) :: lettre
    ! -----------------------------------------------------------------    .
    if (allocated(vP)) deallocate(vP)
    if (allocated(vS)) deallocate(vS)
    if (allocated(h)) deallocate(h)
    ! -----------------------------------------------------------------    .
    nbl=9
    allocate(vP(nbl),vS(nbl),h(nbl))
    ! -----------------------------------------------------------------    .
    if(present(forcelettre)) then
      lettre=forcelettre
    else
      if ((lon.ge.-6.0_wr).and.(lon.le.-1.5_wr).and.(lat.ge.45.0_wr).and.(lat.le.47.0_wr)) then
        lettre='a'
      elseif ((lon.ge.-6.0_wr).and.(lon.le.-3.0_wr).and.(lat.ge.47.0_wr).and.(lat.le.48.0_wr)) then
        lettre='b'
      elseif ((lon.ge.-6.0_wr).and.(lon.le.-3.0_wr).and.(lat.ge.48.0_wr).and.(lat.le.50.0_wr)) then
        lettre='c'
      elseif ((lon.ge.-3.0_wr).and.(lon.le.0.0_wr).and.(lat.ge.49.0_wr).and.(lat.le.50.0_wr)) then
        lettre='d'
      elseif ((lon.ge.-1.5_wr).and.(lon.le.0.0_wr).and.(lat.ge.48.0_wr).and.(lat.le.49.0_wr)) then
        lettre='d'
      elseif ((lon.ge.-3.0_wr).and.(lon.le.-1.5_wr).and.(lat.ge.48.0_wr).and.(lat.le.49.0_wr)) then
        lettre='e'
      elseif ((lon.ge.-3.0_wr).and.(lon.le.-1.5_wr).and.(lat.ge.47.0_wr).and.(lat.le.48.0_wr)) then
        lettre='f'
      elseif ((lon.ge.-1.0_wr).and.(lon.le.0.0_wr).and.(lat.ge.47.0_wr).and.(lat.le.48.0_wr)) then
        lettre='g'
      elseif ((lon.ge.-1.0_wr).and.(lon.le.0.0_wr).and.(lat.ge.46.0_wr).and.(lat.le.47.0_wr)) then
        lettre='h'
      elseif ((lon.ge.-1.0_wr).and.(lon.le.1.0_wr).and.(lat.ge.45.0_wr).and.(lat.le.46.0_wr)) then
        lettre='i'
      elseif ((lon.ge.0.0_wr).and.(lon.le.1.0_wr).and.(lat.ge.46.0_wr).and.(lat.le.47.0_wr)) then
        lettre='j'
      elseif ((lon.ge.0.0_wr).and.(lon.le.1.0_wr).and.(lat.ge.47.0_wr).and.(lat.le.50.0_wr)) then
        lettre='k'
      else
        !write(*,*)'problème dans ModTerrArroucau, le séisme est hors zone (a-k)'
        lettre='z'
      endif
    endif
    ! -----------------------------------------------------------------    .
    !write(*,*) 'LETTRE   ',lettre
    ! -----------------------------------------------------------------    .
    h(1)=0.0_wr
    h(2)=4.0_wr
    h(3)=8.0_wr
    h(4)=12.0_wr
    h(5)=16.0_wr
    h(6)=20.0_wr
    h(7)=24.0_wr
    h(8)=28.0_wr
    h(9)=32.0_wr
    ! -------                                                  --------    .
    if (lettre=='a') then
      vP(1)=6.0_wr
      vP(2)=6.1_wr
      vP(3)=6.1_wr
      vP(4)=6.1_wr
      vP(5)=6.1_wr
      vP(6)=6.8_wr
      vP(7)=7.0_wr
      vP(8)=7.8_wr
      vP(9)=8.0_wr
    elseif (lettre=='b') then
      vP(1)=6.0_wr
      vP(2)=6.0_wr
      vP(3)=6.0_wr
      vP(4)=6.1_wr
      vP(5)=6.1_wr
      vP(6)=6.2_wr
      vP(7)=7.0_wr
      vP(8)=7.4_wr
      vP(9)=8.1_wr
    elseif (lettre=='c') then
      vP(1)=6.0_wr
      vP(2)=6.0_wr
      vP(3)=6.1_wr
      vP(4)=6.1_wr
      vP(5)=6.2_wr
      vP(6)=6.2_wr
      vP(7)=6.9_wr
      vP(8)=7.3_wr
      vP(9)=8.1_wr
    elseif (lettre=='d') then
      vP(1)=6.1_wr
      vP(2)=6.1_wr
      vP(3)=6.1_wr
      vP(4)=6.1_wr
      vP(5)=6.3_wr
      vP(6)=6.5_wr
      vP(7)=7.1_wr
      vP(8)=7.8_wr
      vP(9)=8.2_wr
    elseif (lettre=='e') then
      vP(1)=6.0_wr
      vP(2)=6.0_wr
      vP(3)=6.1_wr
      vP(4)=6.2_wr
      vP(5)=6.3_wr
      vP(6)=6.4_wr
      vP(7)=6.8_wr
      vP(8)=7.3_wr
      vP(9)=8.1_wr
    elseif (lettre=='f') then
      vP(1)=6.0_wr
      vP(2)=6.0_wr
      vP(3)=6.1_wr
      vP(4)=6.1_wr
      vP(5)=6.2_wr
      vP(6)=6.4_wr
      vP(7)=6.8_wr
      vP(8)=7.4_wr
      vP(9)=8.1_wr
    elseif (lettre=='g') then
      vP(1)=6.0_wr
      vP(2)=6.1_wr
      vP(3)=6.2_wr
      vP(4)=6.2_wr
      vP(5)=6.2_wr
      vP(6)=6.4_wr
      vP(7)=6.9_wr
      vP(8)=7.5_wr
      vP(9)=7.9_wr
    elseif (lettre=='h') then
      vP(1)=6.0_wr
      vP(2)=6.1_wr
      vP(3)=6.1_wr
      vP(4)=6.1_wr
      vP(5)=6.2_wr
      vP(6)=6.5_wr
      vP(7)=6.7_wr
      vP(8)=7.3_wr
      vP(9)=8.0_wr
    elseif (lettre=='i') then
      vP(1)=5.2_wr
      vP(2)=6.0_wr
      vP(3)=6.1_wr
      vP(4)=6.1_wr
      vP(5)=6.2_wr
      vP(6)=6.3_wr
      vP(7)=6.7_wr
      vP(8)=7.3_wr
      vP(9)=8.0_wr
    elseif (lettre=='j') then
      vP(1)=6.0_wr
      vP(2)=6.1_wr
      vP(3)=6.1_wr
      vP(4)=6.2_wr
      vP(5)=6.3_wr
      vP(6)=6.4_wr
      vP(7)=6.9_wr
      vP(8)=7.5_wr
      vP(9)=8.0_wr
    elseif (lettre=='k') then
      vP(1)=6.0_wr
      vP(2)=6.1_wr
      vP(3)=6.2_wr
      vP(4)=6.2_wr
      vP(5)=6.2_wr
      vP(6)=6.5_wr
      vP(7)=6.8_wr
      vP(8)=7.4_wr
      vP(9)=8.0_wr
    elseif (lettre=='z') then
      !write(*,*)'nouveau modèle : ModTerrSiHex_Haslach'
      lettre='z'
    else
      write(*,*)'problème dans ModTerrArroucau, le modèle : ',lettre,' n''existe pas '
      stop
    endif
    ! -------                                                  --------    .
    vS=vP/1.68_wr
    ! -------                                                  --------    .
    if (lettre=='z') then
      Call ModTerrSiHex_Haslach(nbl,vP,vS,h,pdfmoho)
    endif
    ! -------                                                  --------    .
    pdfmoho=h(nbl)
    ! -----------------------------------------------------------------    .
  end subroutine ModTerrArroucau

! ---------------------------------------------------------------------    .

  subroutine ModTerrSiHex_Haslach(nbl,vP,vS,h,pdfmoho)
    ! -------                                                  --------    . mh
    ! Modèles de Terre SiHex Haslach
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(inout), allocatable :: h(:),vP(:),vS(:)
    integer (KIND=wi), intent(out) :: nbl
    real(KIND=wr), intent (out) :: pdfmoho
    ! -----------------------------------------------------------------    .
    if (allocated(vP)) deallocate(vP)
    if (allocated(vS)) deallocate(vS)
    if (allocated(h)) deallocate(h)
    ! -----------------------------------------------------------------    .
    nbl=3
    allocate(vP(nbl),vS(nbl),h(nbl))
    ! -----------------------------------------------------------------    .
    vP(1)=5.9_wr
    vP(2)=6.5_wr
    vP(3)=8.2_wr
    ! -------                                                  --------    .
    vS(1)=3.4_wr
    vS(2)=3.7_wr
    vS(3)=4.4_wr
    ! -------                                                  --------    .
    h(1)=0.0_wr
    h(2)=20.0_wr
    h(3)=30.0_wr
    ! -------                                                  --------    .
    pdfmoho=h(nbl)
    ! -----------------------------------------------------------------    .
  end subroutine ModTerrSiHex_Haslach

! ---------------------------------------------------------------------    .

  subroutine ModTerrCEA(nbl,vP,vS,h,pdfmoho)
    ! -------                                                  --------    . mh
    ! Modèles de Terre urilisé au CEA, globale à toute la france.
    ! Le moho est à 25.9 km.
    ! pour les inversions, et pour chaque station : le modèle 1D est moyenné à partir du modèle 3D CRUST2.0.
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(inout), allocatable :: h(:),vP(:),vS(:)
    integer (KIND=wi), intent(out) :: nbl
    real(KIND=wr), intent (out) :: pdfmoho
    ! -----------------------------------------------------------------    .
    if (allocated(vP)) deallocate(vP)
    if (allocated(vS)) deallocate(vS)
    if (allocated(h)) deallocate(h)
    ! -----------------------------------------------------------------    .
    nbl=3
    allocate(vP(nbl),vS(nbl),h(nbl))
    ! -----------------------------------------------------------------    .
    vP(1)=3.00_wr
    vP(2)=6.03_wr
    vP(3)=8.16_wr
    ! -------                                                  --------    .
    vS(1)=1.73_wr
    vS(2)=3.56_wr
    vS(3)=4.65_wr
    ! -------                                                  --------    .
    h(1)=0.0_wr
    h(2)=0.9_wr
    h(3)=25.0_wr
    ! -------                                                  --------    .
    pdfmoho=h(nbl)
    ! -----------------------------------------------------------------    .
  end subroutine ModTerrCEA

! ---------------------------------------------------------------------    .

  subroutine ModTerr_fun(nbl,vP,vS,h,pdfmoho)
    ! -------                                                  --------    . mh
    ! Modèles de Terre fun
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (kind=wr), intent(inout), allocatable :: h(:),vP(:),vS(:)
    integer (KIND=wi), intent(out) :: nbl
    real(KIND=wr), intent (out) :: pdfmoho
    ! -------                                                  --------    .
    integer (KIND=wi) :: i
    ! -----------------------------------------------------------------    .
    if (allocated(vP)) deallocate(vP)
    if (allocated(vS)) deallocate(vS)
    if (allocated(h)) deallocate(h)
    ! -----------------------------------------------------------------    .
    nbl=12
    allocate(vP(nbl),vS(nbl),h(nbl))
    ! -----------------------------------------------------------------    .
    vP(1)=5.2_wr
    vS(1)=vP(1)/1.65_wr
    h(1)=0.0
    do i=2,nbl
        vP(i)=vP(i-1)+1._wr
        h(i)=5.0_wr+h(i-1)
        vS(i)=vP(i)/(1.60+h(i)/100._wr)
    enddo
    ! -------                                                  --------    .
    vP(7)=vP(1)
    vP(10)=vP(7)
    !vS(1)=vP(1)/1.60
    ! -------                                                  --------    .
    pdfmoho=h(nbl)
    ! -----------------------------------------------------------------    .
end subroutine ModTerr_fun

END MODULE ray_tracing



! *********************************************************************    .
! *********************************************************************    .


