! Librairie de subroutines permettant des regressions linéaires
! octobre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE statistiques

    use modparam

    implicit none

    private

    public  :: correlationaffpond, correlationpond
    public  :: inv_normal_cumulative_distrib_func
    public  :: Rpcalc
    public  :: autovariance

    ! -------                                                  --------    .
    interface Rpcalc
      module procedure Rpcalc_bis, Rpcalc_ter   ! différents arguments
    end interface Rpcalc
    ! -------                                                  --------    .

CONTAINS

! ---------------------------------------------------------------------    .

  subroutine correlationaffpond(a,R2,dph,XY)
    ! ------- calcul correlation ponderé affine                --------    .mh
    ! calcul le coeficient directeur a pour XY(:,1) = X ; XY(:,2) = Y et XY(:,3) = pondération [0;1]
    ! R2 correspond au chi2
    ! ordonnée à l'origine nulle
    ! -------                                                  --------    .
    implicit none
    integer (KIND=wi), intent (in) :: dph
    real (KIND=wr), intent (in)  ::  XY(dph,3)
    real (KIND=wr), intent (out) ::  R2,a
    ! -------                                                  --------    .
    real (KIND=wr) :: sxi,syi,b
    integer (KIND=wi) :: i
    ! -----------------------------------------------------------------    .
    if (dph.gt.1) then
      sxi=0.0_wr
      syi=0.0_wr
      do i=1,dph
        sxi=sxi+XY(i,2)*XY(i,1)*XY(i,3)
        syi=syi+XY(i,1)*XY(i,3)*XY(i,1)
      enddo
      a=(sxi)/(syi)
      b=0.0_wr
      call chi2calc(XY,a,b,dph,R2)
      if(R2.gt.0.0_wr) then
        if(IsNaN(a)) then
          write(*,*)'problème dans correlationaffpond : IsNaN(a)',a
          stop
        endif
        if(IsNaN(b)) then
          write(*,*)'problème dans correlationaffpond : IsNaN(a)',b
          stop
        endif
      else
        R2=-1._wr
        a=-1._wr
      endif
    else
      R2=-1._wr
      a=-1._wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine correlationaffpond

! ---------------------------------------------------------------------    .

  subroutine correlationpond(a,b,R2,dph,XY)
    ! ------- calcul correlation ponderé                       --------    .mh
    ! calcul le coeficient directeur a et ordonnée à l'origine b
    ! pour XY(:,1) = X ; XY(:,2) = Y et XY(:,3) = pondération [0;1]
    ! R2 correspond au chi2
    ! -------                                                  --------    .
    implicit none
    integer (KIND=wi), intent (in) :: dph
    real (KIND=wr), intent (in)  ::  XY(dph,3)
    real (KIND=wr), intent (out) ::  R2,a,b
    real (KIND=wr) :: sxi,syi,sxiyi,sxi2,spond
    integer (KIND=wi) :: i
    ! -----------------------------------------------------------------    .
    if (dph.gt.1) then
      a=0.0_wr
      b=0.0_wr
      sxi=0.0_wr
      syi=0.0_wr
      sxiyi=0.0_wr
      sxi2=0.0_wr
      spond=0.0_wr
      do i=1,dph
        spond=spond + XY(i,3)
        sxi=sxi+XY(i,1)*XY(i,3)
        syi=syi+XY(i,2)*XY(i,3)
        sxiyi=sxiyi+XY(i,1)*XY(i,2)*XY(i,3)
        sxi2=sxi2+XY(i,1)*XY(i,1)*XY(i,3)
      enddo
      a=(spond*sxiyi-sxi*syi)/(spond*sxi2-sxi*sxi)
      b=(syi*sxi2-sxi*sxiyi)/(spond*sxi2-sxi*sxi)
      call chi2calc(XY,a,b,dph,R2)
      if(R2.gt.0.0_wr) then
        if(IsNaN(a)) then
          write(*,*)'problème dans correlationpond : IsNaN(a)',a
          stop
        endif
        if(IsNaN(b)) then
          write(*,*)'problème dans correlationpond : IsNaN(b)',b
          stop
        endif
      else
        R2=-1._wr
        a=-1._wr
      endif
    else
      R2=-1._wr
      a=-1._wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine correlationpond

! ---------------------------------------------------------------------    .

  subroutine chi2calc(XY,a,b,dph,chi2)
    ! -------                                                  --------    .mh
    ! calcul du chi2
    ! -------                                                  --------    .
    implicit none
    integer (KIND=wi), intent (in) :: dph
    real (KIND=wr), intent (in)  :: XY(dph,3),a,b
    real (KIND=wr), intent (out) :: chi2
    ! -------                                                  --------    .
    real (KIND=wr) ::  yth(dph), test
    integer (KIND=wi) i, ik
    ! -----------------------------------------------------------------    .
    chi2=0.0_wr
    do ik=1,dph
      yth(ik)=a*XY(ik,1)+b
      chi2=chi2+(XY(ik,2)-yth(ik))**(2.0_wr)
    enddo
    chi2=chi2/real(dph,wr)
    if(IsNaN(chi2)) then
      write(*,*)'problème dans chi2calc : IsNaN(chi2)',a,b,dph,chi2
      test=0.0_wr
      do i=1,dph
        ! write(*,*)i,XY(i,1),XY(i,2),XY(i,3)
        test=test+XY(i,3)
      enddo
      if (test.ge.0.00000001_wr) then
        stop
      else
        chi2=-999.99_wr
      endif
    endif
    ! -----------------------------------------------------------------    .
  end subroutine chi2calc

! ---------------------------------------------------------------------    .

  subroutine Rpcalc_bis(XY,dph,nb,Rp)
    ! -------                                                  --------    .mh
    ! calcul du Rp (abs), Coefficient de corrélation linéaire de Bravais-Pearson
    ! -------                                                  --------    .
    implicit none
    integer (KIND=wi), intent (in) :: dph,nb
    real (KIND=wr), intent (in)  :: XY(nb,3)
    real (KIND=wr), intent (out) :: Rp
    ! -------                                                  --------    .
    real (KIND=wr) ::  A, B, C, x, y
    integer (KIND=wi) i
    ! -----------------------------------------------------------------    .
    if (dph.gt.2) then
      ! ---------------------------------------------------------------    .
      if (nb.lt.dph) then
        write(*,*)'problème dans Rpcalc_bis : i < dph ',nb,' < ',dph
        stop
      endif
      ! -------                                                --------    .
      A = 0.0_wr
      B = 0.0_wr
      C = 0.0_wr
      x = 0.0_wr
      y = 0.0_wr
      ! -------                                                --------    .
      do i=1,dph
        if (IsNaN(X)) then
          write(*,*)'problème dans Rpcalc_bis : IsNaN(x)',x
          stop
        endif
        if (IsNaN(y)) then
          write(*,*)'problème dans Rpcalc_bis : IsNaN(y)',y
          stop
        endif
        x = x + XY(i,1)
        y = y + XY(i,2)
      enddo
      x = x/real(dph,wr)
      y = y/real(dph,wr)
      ! -------                                                --------    .
      do i=1,dph
        A = A + (XY(i,1) - x) * (XY(i,2) - y)
        B = B + (XY(i,1) - x)**2.0_wr
        C = C + (XY(i,2) - y)**2.0_wr
      enddo
      ! -------                                                --------    .
      Rp = abs(A/sqrt(B*C))
      ! -------                                                --------    .
      if (IsNaN(Rp)) then
        write(*,*)A, B, C, x, y, dph,nb
        write(*,*)'problème dans Rpcalc_bis : IsNaN(Rp) : ',Rp,dph,nb
        write(*,*)XY(1,:)
        write(*,*)XY(2,:)
        write(*,*)
        Rp=0.0_wr ! si que pg
        !stop
      endif
      ! ---------------------------------------------------------------    .
    else
      ! ---------------------------------------------------------------    .
      Rp=0.0_wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine Rpcalc_bis

! ---------------------------------------------------------------------    .

  subroutine Rpcalc_ter(vX,vY,nb,Rp)
    ! -------                                                  --------    .mh
    ! calcul du Rp, Coefficient de corrélation linéaire de Bravais-Pearson
    ! -------                                                  --------    .
    implicit none
    integer (KIND=wi), intent (in) :: nb
    real (KIND=wr), intent (in)  :: vX(nb), vY(nb)
    real (KIND=wr), intent (out) :: Rp
    ! -------                                                  --------    .
    real (KIND=wr) ::  A, B, C, x, y
    integer (KIND=wi) i
    ! -----------------------------------------------------------------    .
    if (nb.gt.2) then
      ! ---------------------------------------------------------------    .
      A = 0.0_wr
      B = 0.0_wr
      C = 0.0_wr
      x = 0.0_wr
      y = 0.0_wr
      ! -------                                                --------    .
      do i=1,nb
        if (IsNaN(X)) then
          write(*,*)'problème dans Rpcalc_ter : IsNaN(x)',x
          stop
        endif
        if (IsNaN(y)) then
          write(*,*)'problème dans Rpcalc_ter : IsNaN(y)',y
          stop
        endif
        x = x + vX(i)
        y = y + vY(i)
      enddo
      x = x/real(nb,wr)
      y = y/real(nb,wr)
      ! -------                                                --------    .
      do i=1,nb
        A = A + (vX(i) - x) * (vY(i) - y)
        B = B + (vX(i) - x)**2.0_wr
        C = C + (vY(i) - y)**2.0_wr
      enddo
      ! -------                                                --------    .
      Rp = A/sqrt(B*C)
      ! -------                                                --------    .
      if (IsNaN(Rp)) then
        write(*,*)A, B, C, x, y, nb
        write(*,*)'problème dans Rpcalc_ter : IsNaN(Rp)',Rp,nb
        write(*,*)
        stop
      endif
      ! ---------------------------------------------------------------    .
    else
      ! ---------------------------------------------------------------    .
      Rp=0.0_wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine Rpcalc_ter

! ---------------------------------------------------------------------    .

  subroutine inv_normal_cumulative_distrib_func(p,dinvnorm)
    ! -------                                                  --------    .mh
    ! free, fast, and accurate way of computing the inverse normal cumulative distribution function.
    ! algorithm with a relative error less than 1.15 x 10-9 in the entire region
    ! on sait : p(dinvnorm)=0.5_wr*erfc(-dinvnorm/sqrt(2.0_wr))
    ! on cherche : dinvnorm(p)=?
    ! Peter John Acklam (pjacklam@online.no)
    ! http://home.online.no/~pjacklam/notes/invnorm
    ! -------                                                  --------    .
    implicit none
    real (KIND=wr), intent (in)  :: p
    real (KIND=wr), intent (out) :: dinvnorm
    ! -------                                                  --------    .
    real (KIND=wr) :: p_low,p_high
    real (KIND=wr) :: a1,a2,a3,a4,a5,a6
    real (KIND=wr) :: b1,b2,b3,b4,b5
    real (KIND=wr) :: c1,c2,c3,c4,c5,c6
    real (KIND=wr) :: d1,d2,d3,d4
    real (KIND=wr) :: z,q,r
    ! -----------------------------------------------------------------    .
    a1=-39.6968302866538_wr
    a2=220.946098424521_wr
    a3=-275.928510446969_wr
    a4=138.357751867269_wr
    a5=-30.6647980661472_wr
    a6=2.50662827745924_wr
    b1=-54.4760987982241_wr
    b2=161.585836858041_wr
    b3=-155.698979859887_wr
    b4=66.8013118877197_wr
    b5=-13.2806815528857_wr
    c1=-0.00778489400243029_wr
    c2=-0.322396458041136_wr
    c3=-2.40075827716184_wr
    c4=-2.54973253934373_wr
    c5=4.37466414146497_wr
    c6=2.93816398269878_wr
    d1=0.00778469570904146_wr
    d2=0.32246712907004_wr
    d3=2.445134137143_wr
    d4=3.75440866190742_wr
    p_low=0.02425_wr
    p_high=1.0_wr-p_low
    ! -------                                                  --------    .
    if ((p.le.0.0_wr).or.(p.ge.1.0_wr)) then
      write(*,*)'problème dans inv_normal_cumulative_distrib_func : p= ',p
    endif
    ! -------                                                  --------    .
    if(p.lt.p_low) then
      q=sqrt(-2.0_wr*log(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_wr)
      dinvnorm=z
    else
      if(p.le.p_high) then
        q=p-0.5_wr
        r=q*q
        z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_wr)
        dinvnorm=z
      elseif(p.lt.1.0_wr) then
        q=sqrt(-2.0_wr*log(1.0_wr-p))
        z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_wr)
        dinvnorm=z
      endif
    endif
    ! -------                                                  --------    .
    if (IsNaN(dinvnorm)) then
      write(*,*)'problème dans inv_normal_cumulative_distrib_func : IsNaN(dinvnorm)',dinvnorm
      stop
    endif
    ! -----------------------------------------------------------------    .
  end subroutine inv_normal_cumulative_distrib_func

    ! -----------------------------------------------------------------    .

  subroutine autovariance(vec,itermax,k,nom)
    ! -------                                                  --------    .
    ! Calcul de la fonction d’autocovariance, Ck (mesure la covariance entre une variable
    ! et cette même variable à des dates différentes, pour un délai k)
    ! La fonction d’autocovariance (Ck) est normalisée en fonction d’autocorrélation (rk)
    !
    ! La représentation des fonctions d’autocorrélation des paramètres permet de donner
    ! des indications sur le nombre d’itérations requis pour que les valeurs des paramètres
    ! échantillonnées par l’algorithme soient décorrélées, ainsi que sur la valeur de l’écarttype
    ! de la gaussienne à employer (e.g. Drilleau, 2013).

    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(in) :: itermax,k
    real(KIND=wr), intent(in) :: vec(itermax)
    character(len=50), intent(in) :: nom
    ! -------                                                  --------    .
    integer(kind=wi) :: i, t, ok
    real(KIND=wr) :: moy,C0,Ck,rk
    integer(kind=wi) :: num=0
    ! -------                                                  --------    .
    num=num+1
    ! -------                                                  --------    .
    moy=0.0_wr
    do i=1,itermax
      moy=moy+vec(i)
    enddo
    moy=moy/real(itermax,wr)
    ! -------                                                  --------    .
    C0=0.0_wr
    do t=1,itermax-1
      C0=C0+(vec(t)-moy)**2.0_wr
    enddo
    C0=C0/real(itermax-1,wr)
    ! -------                                                  --------    .
    ok=0
    open(unit=5000+num,file=nom,STATUS="replace",iostat=ok)
    if (ok .ne. 0) then
      write(*,*)"problème dans autovariance : le fichier "//nom//" n''existe pas "
      stop
    endif
    ! -------                                                  --------    .
    do i=0,k,1
      Ck=0.0_wr
      do t=1,itermax-i-k
        Ck=Ck+(vec(t)-moy)*(vec(t+i)-moy)
      enddo
      Ck=Ck/real(itermax-i-k,wr)
      rk=Ck/C0
      write(5000+num,*)i,rk
    enddo
    close(5000+num)
    ! -------                                                  --------    .
  end subroutine autovariance

END MODULE statistiques



! *********************************************************************    .
! *********************************************************************    .


