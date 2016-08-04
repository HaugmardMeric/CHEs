! module permetant le calcul des points d'intersection entre deux cercles
! quelconques sur une sphère.
! fevrier 2014
! *********************************************************************    .mh
! ------- Philippe Cance philippe.cance@univ-nantes.fr         --------    .pc
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .mh
! ------- Ianis Gaudot ianis.gaudot@univ-nantes.fr             --------    .ig
! *********************************************************************    .
! ---------------------------------------------------------------------    .


MODULE dist_cercle

    use modparam

    implicit none

    private

    public  :: dist2c

    real(KIND=wr), parameter :: d2r=pi/180._wr
    real(KIND=wr), parameter :: r2d=180._wr/pi


CONTAINS

! ---------------------------------------------------------------------    .


  subroutine dist2c(t1,t2,p1,p2,a1,a2,t_1,p_1,t_2,p_2)
    ! -------                                                  --------    .mh
    ! calcul des points d'intersection entre deux cercles quelconques sur une sphère
    ! t (theta) : latitudes [-90;90]
    ! p (phi) : longitudes [-180;180]
    ! t et p représente le centre du cecle, à la surface (Rt=6371 km)
    ! a : demi-angle d'ouveture du cône [0;180], depuis le centre
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real(KIND=wr), intent(in) :: t1,t2,p1,p2,a1,a2
    real(KIND=wr), intent(out) :: t_1,p_1,t_2,p_2
    ! -------                                                  --------    .
    real(KIND=wr) :: t11,t21,p11,p21
    ! -----------------------------------------------------------------    . conversion lat -> colat
    t11=90.0_wr-t1
    t21=90.0_wr-t2
    ! -----------------------------------------------------------------    . conversion lon -> colon
    p11=p1
    p21=p2
    if (p1.lt.0.0_wr) then
      p11=360.0_wr+p1
    endif
    if (p2.lt.0.0_wr) then
      p21=360.0_wr+p2
    endif
    ! -----------------------------------------------------------------    . calcul des points d'intersection
    call dist_2sc(t11,t21,p11,p21,a1,a2,t_1,p_1,t_2,p_2)
    ! -----------------------------------------------------------------    . if Okay
    if ((p_1.lt.1000.0_wr).and.(p_2.lt.1000.0_wr).and.(t_1.lt.1000.0_wr).and.(t_2.lt.1000.0_wr)) then
      ! -------                                                --------    . conversion colat -> lat
      t_1=90.0_wr-t_1
      t_2=90.0_wr-t_2
      ! -------                                                --------    . conversion colon -> lon
      if (p_1.gt.180.0_wr) then
        p_1=p_1-360.0_wr
      endif
      if (p_2.gt.180.0_wr) then
        p_2=p_2-360.0_wr
      endif
      ! -------                                                --------    . else
      if ((p_1.lt.-170.0_wr).or.(p_2.lt.-170.0_wr).or.(p_1.gt.170.0_wr).or.(p_2.gt.170.0_wr)) then
        write(*,*)'problème dans dist2c : longitude actuelle trop près des coutures !',p_1,p_2
        stop
      endif
      if ((t_1.lt.-80.0_wr).or.(t_2.lt.-80.0_wr).or.(t_1.gt.80.0_wr).or.(t_2.gt.80.0_wr)) then
        write(*,*)'problème dans dist2c : latitude actuelle trop près des coutures !',t1,t2
        stop
      endif
      ! -------                                                --------    .
    endif
    ! -----------------------------------------------------------------    .
  end subroutine dist2c

! ---------------------------------------------------------------------    .

    subroutine dist_2sc(theta1,theta2,phi1,phi2,alpha1,alpha2,thetaA,phiA,thetaB,phiB)
    ! -------                                                  --------    .pc
    ! calcul des points d'intersection entre deux cercles quelconques sur une sphère
    ! Philippe Cance
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real(KIND=wr), intent(in) :: theta1, theta2, phi1, phi2, alpha1, alpha2
    real(KIND=wr), intent(out) :: thetaA, phiA, thetaB, phiB
    ! -------                                                  --------    .
    real(KIND=wr) :: heisenberg
    real(KIND=wr) :: v1(3), v2(3), v3(3), M0(3), M1(3)
    real(KIND=wr) :: a(2),b(2),c(2),alpha(2),tmp1(2)
    real(KIND=wr) :: aA,bB,cC,alphaA,alphaB,alphaC
    real(KIND=wr) :: detla, t
    ! -----------------------------------------------------------------    .
    heisenberg=real(3,wr)*spacing(pi)
    ! -----------------------------------------------------------------    .
    if ((abs(theta1-theta2).lt.heisenberg).and.(abs(sin(phi1-phi2)).lt.heisenberg)) then
      write(*,*)'problème dans dist_2sc : cercles paralleles !'
      thetaA=10000.0_wr
      phiA=10000.0_wr
      thetaB=thetaA
      phiB=phiA
      stop
    endif
    ! -----------------------------------------------------------------    .
    v1 = (/ cos(d2r*phi1)*sin(d2r*theta1), sin(d2r*phi1)*sin(d2r*theta1), cos(d2r*theta1) /) ! v1=vec(OO1)
    v2 = (/ cos(d2r*phi2)*sin(d2r*theta2), sin(d2r*phi2)*sin(d2r*theta2), cos(d2r*theta2) /)! v2=vec(OO2)
    call vec_prod(v1, v2, v3) ! v3=vec(OO1)^vec(OO2)
    a=(/v1(1),v2(1)/)
    b=(/v1(2),v2(2)/)
    c=(/v1(3),v2(3)/)
    alpha=(/cos(d2r*alpha1),cos(d2r*alpha2)/)
    aA=det2(b,c)
    bB=det2(c,a)
    cC=det2(a,b)
    alphaA=det2(alpha,a)
    alphaB=det2(alpha,b)
    alphaC=det2(alpha,c)
    ! -----------------------------------------------------------------    .
    v1 = (/ alphaA, alphaB, alphaC /)
    v2 = (/ aA, bB, cC /)
    call vec_prod(v1,v2,M0)
    M0 = M0 / dot_product(v2,v2) ! M0 = P1 inter P2 inter (O,O1,O2)
    detla=dot_product(M0,M0)
    ! -----------------------------------------------------------------    .
    if (detla.gt.(1.0_wr)) then
      ! write(*,*)'problème dans dist_2sc : cercles non secants !',detla
      thetaA=10000.0_wr
      phiA=10000.0_wr
      thetaB=thetaA
      phiB=phiA
      ! ---------------------------------------------------------------    . un unique point d'intersection
    else if (abs(detla-1.0_wr).lt.heisenberg) then
      M1 = M0
      thetaA = acos(M1(3))
      tmp1 = M1(1:2) / sin(thetaA)
      thetaA = r2d * thetaA
      phiA= r2d * acos(tmp1(1))
      if (tmp1(2).lt.0.0_wr) phiA= phiA+ 180.0_wr
      thetaB=thetaA
      phiB=phiA
      ! ---------------------------------------------------------------    . deux points d'intersection
    else
      ! -------                                                --------    . solution 1
      t = sqrt((1.0_wr - detla) / dot_product(v3,v3))
      M1 = M0 + t * v3
      thetaA = acos(M1(3))
      tmp1 = M1(1:2) / sin(thetaA)
      thetaA = r2d * thetaA
      phiA= r2d * acos(tmp1(1))
      if (abs(tmp1(2)) .gt. heisenberg) then
        if (tmp1(2).lt.0.0_wr) phiA= 360.0_wr - phiA
      endif
      ! -------                                                --------    . solution 2
      M1 = M0 - t * v3
      thetaB = acos(M1(3))
      tmp1 = M1(1:2) / sin(thetaB)
      thetaB = r2d * thetaB
      phiB= r2d * acos(tmp1(1))
      if (abs(tmp1(2)) .gt. heisenberg) then
        if (tmp1(2).lt.0.0_wr) phiB= 360.0_wr - phiB
      endif
    endif
    ! -----------------------------------------------------------------    .
end subroutine dist_2sc

! ---------------------------------------------------------------------    .

  function det2(a1,a2) result(res)
    ! -------                                                  --------    .pc
    ! determinant
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real(KIND=wr),dimension(2),intent(in) :: a1,a2
    real(KIND=wr) :: res
    ! -----------------------------------------------------------------    .
    res = a1(1)*a2(2)-a1(2)*a2(1)
    ! -----------------------------------------------------------------    .
  end function det2

! ---------------------------------------------------------------------    .

  subroutine vec_prod(u1,u2,res)
    ! -------                                                  --------    .pc
    ! produit vectoriel
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real(KIND=wr),dimension(3),intent(in):: u1,u2
    real(KIND=wr),dimension(3),intent(out)::res
    real(KIND=wr),dimension(2) :: tmp1, tmp2
    ! -----------------------------------------------------------------    .
    tmp1 = u1(2:3)
    tmp2 = u2(2:3)
    res(1) = det2(tmp1,tmp2)
    tmp1 = (/ u1(3), u1(1) /)
    tmp2 = (/ u2(3), u2(1) /)
    res(2) = det2(tmp1,tmp2)
    tmp1 = u1(1:2)
    tmp2 = u2(1:2)
    res(3) = det2(tmp1,tmp2)
    ! -----------------------------------------------------------------    .
  end subroutine vec_prod

! ---------------------------------------------------------------------    .

END MODULE dist_cercle


! *********************************************************************    .
! *********************************************************************    .




