! novembre 2013
! calcul la distance épicentrale et le back-azimuth (en degrés)
! *********************************************************************    .mh
! ------- Eric Beucler eric.beucler@univ-nantes.fr (12.07.2005) -------    .eb
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .mh
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE distance_epi

    use modparam

    implicit none

    private

    public  :: dellipsgc


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine dellipsgc(dlat1,dlon1,dlat2,dlon2,d,bazfin)
    ! -------                                                  --------    .eb
    ! calcul plus complexe de la distance épicentrale et du back-azimuth
    ! -------                                                  --------    .mh
    implicit none
    real (kind=wr), intent(in) :: dlat1,dlat2,dlon1,dlon2                  ! coordonnées des points [en degrés]
    real (kind=wr), intent(out), optional :: bazfin                        ! back-azimuth [en degrés :0;360]
    real (kind=wr), intent(out) :: d                                       ! distance épicentrale [en kilomètres]
    ! -------                                                  --------    .
    real (kind=wr) :: val, a, f, az, baz
    real (kind=wr) :: alph0,alpha,b,bb,d2r,delta,dphi,eps,i1,i2,m
    real (kind=wr) :: nx,ny,nz,p,rlat1,rlat2,rlatn,rlon1,rlon2
    ! -----------------------------------------------------------------    .
    val=abs(dlat1-dlat2)+abs(dlon1-dlon2)
    if(val.gt.(real(2,wr)*spacing(dlat1))) then                            ! pb si distance épicentre => zéro
    ! -----------------------------------------------------------------    .
      az = 0.0_wr
      baz = 0.0_wr
      d = 0.0_wr
      eps = real(2,wr)*spacing(dlat1)
      d2r = pi/180.0_wr
      if ((abs(dlat1-dlat2) < eps).and.&
           ((abs(dlon1-dlon2) < eps).or.(abs(dlon1-dlon2) == 360.0_wr))) then
        write(*,*)'problème dans dellipsgc : mauvaise coordonnées'
      endif
      a = rT * 1000.0_wr
      f = 0.0_wr                                                           ! sphère sans aplatissement au pôles
      b = a*(1._wr-f)
      ! ---------------------------------------------------------------    .
      ! Conversion to radians and verification of non nullity:
      ! Latitude  1:
      rlat1 = dlat1*d2r
      if (abs(rlat1) < eps) then
        if (rlat1 < 0.0_wr) rlat1 = -1.0_wr*eps
        if (rlat1 >= 0.0_wr) rlat1 = eps
      end if
      ! Longitude 1:
      if (dlon1 > (180.0_wr+eps)) then
        rlon1 = (dlon1-360.0_wr)*d2r
      else if (dlon1 < (-180.0_wr-eps)) then
        rlon1 = (dlon1+360.0_wr)*d2r
      else
        rlon1 = dlon1*d2r
      end if
      if (abs(rlon1) < eps) then
        if (rlon1 < 0.0_wr) rlon1 = -1.0_wr*eps
        if (rlon1 >= 0.0_wr) rlon1 = eps
      end if
      if ((dlat1 == 90.0_wr).or.(dlat1 == -90.0_wr)) rlon1 = eps
      ! Latitude  2:
      rlat2 = dlat2*d2r
      if (abs(rlat2) < eps) then
        if (rlat2 < 0.0_wr) rlat2 = -1.0_wr*eps
        if (rlat2 >= 0.0_wr) rlat2 = eps
      end if
      ! Longitude 1:
      if (dlon2 > (180.0_wr+eps)) then
        rlon2 = (dlon2-360.0_wr)*d2r
      else if (dlon2 < (-180.0_wr-eps)) then
        rlon2 = (dlon2+360.0_wr)*d2r
      else
        rlon2 = dlon2*d2r
      end if
      if (abs(rlon2) < eps) then
        if (rlon2 < 0.0_wr) rlon2 = -1.0_wr*eps
        if (rlon2 >= 0.0_wr) rlon2 = eps
      end if
      if ((dlat2 == 90.0_wr).or.(dlat2 == -90.0_wr)) rlon2 = eps
      dphi = rlon2-rlon1
      ! ---------------------------------------------------------------    .
      ! Computes path angle (identical to the spherical case):
      if (((cos(rlat1)*cos(rlat2)*cos(rlon1-rlon2))+ &
            (sin(rlat1)*sin(rlat2))) >= 1.0_wr) then
        delta = eps
      elseif (((cos(rlat1)*cos(rlat2)*cos(rlon1-rlon2))+ &
            (sin(rlat1)*sin(rlat2))) <= -1.0_wr) then
        delta = pi-eps
      else
        delta = acos((cos(rlat1)*cos(rlat2)*cos(rlon1-rlon2))+ &
            (sin(rlat1)*sin(rlat2)))
      end if
      ! ---------------------------------------------------------------    .
      ! Computes azimuth (angle between north and path 1->2):
      if (((sin(dphi)*cos(rlat2))/sin(delta)) > 1.0_wr) then
        az = pi/2.0_wr
      elseif (((sin(dphi)*cos(rlat2))/sin(delta)) < -1.0_wr) then
        az = -1.0_wr*pi/2.0_wr
      else
        az = asin((sin(dphi)*cos(rlat2))/sin(delta))
        if (((sin(rlat2)-(sin(rlat1)*cos(delta)))/&
          (cos(rlat1)*sin(delta))) < 0.0_wr) az = pi-az
      end if
      if ((abs(dphi) < eps).and.(rlat1 > rlat2)) az = 180.0_wr*d2r
      if (abs(az) < eps) then
        az = 0.0_wr
      elseif (az >= (2.0_wr*pi)) then
        az = (2.0_wr*pi)-az
      elseif (az < 0.0_wr) then
        az = (2.0_wr*pi)+az
      end if
      ! ---------------------------------------------------------------    .
      ! Computes back-azimuth (angle between north and path 2->1):
      dphi = -1.0_wr*dphi
      if (((sin(dphi)*cos(rlat1))/sin(delta)) > 1.0_wr) then
        baz = pi/2.0_wr
      elseif (((sin(dphi)*cos(rlat1))/sin(delta)) < -1.0_wr) then
        baz = -1.0_wr*pi/2.0_wr
      else
        baz = asin(sin(dphi)*cos(rlat1)/sin(delta))
        if (((sin(rlat1)-(sin(rlat2)*cos(delta)))/&
          (cos(rlat2)*sin(delta))) < 0.0_wr) baz = pi-baz
      end if
      if ((abs(dphi) < eps).and.(rlat2 > rlat1)) baz = 180.0_wr*d2r
      if (abs(baz) < eps) then
        baz = 0.0_wr
      elseif (baz >= (2.0_wr*pi)) then
        baz = (2.0_wr*pi)-baz
      elseif (baz < 0.0_wr) then
        baz = (2.0_wr*pi)+baz
      end if
      ! ---------------------------------------------------------------    .
      ! Computes the angle between great circle and the equatorial plane:
      nx = (cos(rlat1)*sin(rlon1)*sin(rlat2))- &
       (sin(rlat1)*cos(rlat2)*sin(rlon2))
      ny = (sin(rlat1)*cos(rlat2)*cos(rlon2))- &
       (cos(rlat1)*cos(rlon1)*sin(rlat2))
      nz = (cos(rlat1)*cos(rlon1)*cos(rlat2)*sin(rlon2))- &
       (cos(rlat1)*sin(rlon1)*cos(rlat2)*cos(rlon2))
      if ((sqrt(nx**2+ny**2) < eps).and.(abs(nz) > eps)) then
        rlatn = pi/2.0_wr
      elseif ((sqrt(nx**2+ny**2) < eps).and.(abs(nz) < eps)) then
        rlatn = 0.0_wr
      else
        rlatn = atan(nz/sqrt(nx**2+ny**2))
      end if
      rlatn = (pi/2.0_wr)+rlatn
      if (rlatn > (pi/2.0_wr)) rlatn = pi-rlatn
      ! ---------------------------------------------------------------    .
      ! Computes the new length of the half axe:
      bb = (a*b)/sqrt(((b*cos(rlatn))**2)+((a*sin(rlatn))**2))
      if (bb < b) bb = b
      if (bb > a) bb = a
      ! ---------------------------------------------------------------    .
      ! Computes the initial angle in the great circle plane:
      if (abs(rlatn) < eps) then
        alph0 = 0.0_wr
      elseif (abs(rlatn-pi) < eps) then
        alph0 = pi
      else
        if (rlat1 <= rlat2) then
          if (abs(sin(rlat1)/sin(rlatn)) >= 1.0_wr) then
            if ((rlat1*rlatn) < 0.0_wr) then
              alph0 = -1.0_wr*(pi/2.0_wr)
            else
              alph0 = pi/2.0_wr
            end if
          else
            alph0 = asin(sin(rlat1)/sin(rlatn))
          end if
        else
          if (abs(sin(rlat2)/sin(rlatn)) >= 1.0_wr) then
            if ((rlat2*rlatn) < 0.0_wr) then
              alph0 = -1.0_wr*(pi/2.0_wr)
            else
              alph0 = pi/2.0_wr
            end if
          else
           alph0 = asin(sin(rlat2)/sin(rlatn))
          end if
        end if
      end if
      ! ---------------------------------------------------------------    .
      ! Computes the curvilign integral along the elliptic path:
      m = 1.0_wr-(a**2/bb**2)
      if (abs(m) < eps) m = 0.0_wr
      p = 2.0_wr*pi*sqrt(0.5_wr*(a**2+bb**2)) ! Perimeter of the ellipse.
      ! ---------------------------------------------------------------    .
      ! Second integral between 0 and alph0:
      i2 = 0.0_wr
      alpha = alph0
      if (abs(1.0_wr-sin(alpha)**2) < eps) then
        i2 = p/4.0_wr
        if (alpha < 0.0_wr) i2 = -i2
      else
        if (alpha > (pi/2.0_wr)) then
          i2 = p/2.0_wr
          alpha = pi-alpha
          i2 = i2-(a*sqrt(1.0_wr-sin(alpha)**2)*asin(sin(alpha))/&
               sqrt((1.0_wr-(m*sin(alpha)**2))*(1.0_wr-sin(alpha)**2)))
        else
          i2 = a*sqrt(1.0_wr-sin(alpha)**2)*asin(sin(alpha))/&
                sqrt((1.0_wr-(m*sin(alpha)**2))*(1.0_wr-sin(alpha)**2))
        end if
      end if
      ! ---------------------------------------------------------------    .
      ! First integral between 0 and (alph0+delta):
      i1 = 0.0_wr
      alpha = alph0+delta
      if (abs(1.0_wr-sin(alpha)**2) < eps) then
        i1 = p/4.0_wr
        if (alpha < 0.0_wr) i1 = -i1
      else
        if (alpha > (pi/2.0_wr)) then
          i1 = p/2.0_wr
          alpha = pi-alpha
          i1 = i1-(a*sqrt(1.0_wr-sin(alpha)**2)*asin(sin(alpha))/&
                sqrt((1.0_wr-(m*sin(alpha)**2))*(1.0_wr-sin(alpha)**2)))
        else
          i1 = a*sqrt(1.0_wr-sin(alpha)**2)*asin(sin(alpha))/&
                sqrt((1.0_wr-(m*sin(alpha)**2))*(1.0_wr-sin(alpha)**2))
        end if
      end if
      d = i1-i2
      ! ---------------------------------------------------------------    .
      ! Back to kilometers and degrees:
      d = d/1000.0_wr
      az = az/d2r
      baz = baz/d2r
    else
      ! write(*,*)" problème dans dellipsgc : l'épicentre est trop porche de la station ! ",val
      d = 0.000000001_wr
      az = 0.000000001_wr
    endif
    ! -----------------------------------------------------------------    . BAZ, output in present
    if(present(bazfin)) bazfin=baz
    ! -----------------------------------------------------------------    .
  end subroutine dellipsgc

END MODULE distance_epi



! *********************************************************************    .
! *********************************************************************    .


