! Librairie de subroutines permettant le tirage aléatoire des paramètres
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE tirage

    use modparam
    use typetemps, only : parametresinv, date_sec
    use mt19937

    implicit none

    private

    public  :: reflexion
    public  :: tirage_T, tirage_H, tirageUN


CONTAINS

    ! -----------------------------------------------------------------    .

  subroutine reflexion(valin,min,max,valout)
    ! -------                                                  --------    .mh
    ! verifie si la valeur est dans le prior
    ! sinon : effectue une réflexion
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr),intent(in) :: valin,min,max
    real(KIND=wr),intent(out) :: valout
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    ! -----------------------------------------------------------------    .
    ! -------                                                  --------    .
    if (max.lt.min) write(*,*)'problème dans reflexion : ',max,'<',min
    i=0
    valout=valin
    do while ((valout.ge.max).or.(valout.le.min))                          ! tant que hors des bornes
      i=i+1
      ! -------                                                --------    . réflexion supérieure
      if (valout.ge.max) then
        valout=valout-2.0_wr*abs(valout-max)
      end if
      ! -------                                                --------    . réflexion inférieure
      if (valout.le.min) then
        valout=valout+2.0_wr*abs(min-valout)
      end if
      ! -------                                                --------    .
      if(i.gt.250) then
        write(*,*)'problème dans reflexion : problème tirage hors bornes', valin,min,max,valout
        stop
      endif
      ! -------                                                --------    . fin tant que
    enddo
    if (IsNaN(valout)) write(*,*)'problème dans reflexion : IsNaN(valout)', valin,min,max,valout
    ! -----------------------------------------------------------------    .
  end subroutine reflexion

! ---------------------------------------------------------------------    .

  subroutine tirage_norm(vNew,vOld,vmini,vmaxi,vecartype)
    ! -------                                                  --------    .mh
    ! tirage d'un paramètre selon une loi normale bornée
    ! modif : non ! 
    ! c'est plutot un changement de variable ! pour tiage normal
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr),intent(in) :: vOld,vmini,vmaxi,vecartype
    real(KIND=wr),intent(out) :: vNew
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    ! -----------------------------------------------------------------    .
    vNew = normal(vOld,vecartype)                                          ! tirage aléatoire
    ! -------                                                  --------    .
    i=0
    do while ((vNew.ge.vmaxi).or.(vNew.le.vmini))                          ! tant que hors des bornes
      i=i+1
      ! -------                                                --------    . réflexion supérieure
      if (vNew.ge.vmaxi) then
        vNew=vNew-2.0_wr*abs(vNew-vmaxi)
      end if
      ! -------                                                --------    . réflexion inférieure
      if (vNew.le.vmini) then
        vNew=vNew+2.0_wr*abs(vmini-vNew)
      end if
      ! -------                                                --------    .
      if(i.eq.200) then
        write(*,*)'avertissement dans tirage_norm : tirage(param) -> iter > 100',vNew,vOld,vmini,vmaxi,vecartype
        vNew = normal(vOld,vecartype)
      endif
      ! -------                                                --------    .
      if(i.gt.250) then
        write(*,*)'problème dans tirage_norm : problème tirage_norm',vNew,vOld,vmini,vmaxi,vecartype
        write(*,*)'souvent : si pas de contraintes sur moho ; il remonte et Zhypo est bien au delà !!!!! '
        vNew=vmini+0.5_wr*abs(vmini-vmaxi)
        write(*,*)'RECENTRE !',vNew
        exit
      endif
      ! -------                                                --------    . fin tant que
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine tirage_norm

! ---------------------------------------------------------------------    .

  subroutine tirage_norm_lon(vNew,vOld,vmini,vmaxi,vecartype,lat,X,Y,R,dist_et)
    ! -------                                                  --------    .mh
    ! tirage de la longitude selon une loi normale bornée -> le prior de l'épicentre est ici un cercle !
    ! -------                                                  --------    .
    use dist_cercle
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr),intent(inout) :: vOld,vmini,vmaxi,vecartype
    real(KIND=wr),intent(in) :: X,Y,R,lat,dist_et
    real(KIND=wr),intent(out) :: vNew
    integer(KIND=wi), save :: x_0=0
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    real(KIND=wr) :: nR,Lo1,Lo2,La1,La2,alpha,eclat
    ! -----------------------------------------------------------------    . calcul de l'écart-type
    eclat  = dist_et * 360.0_wr / ( 2.0_wr * pi * rT)
    vecartype = eclat / sin((90.0_wr-lat)/180.0_wr*pi)
    ! -----------------------------------------------------------------    . calcul des bornes min et max
    ! petit cercle : centre(X,Y) et rayon (R, en km)
    nR=180.0_wr*R/pi/rT
    ! petit cercle : de latitude (lat) : 2pts {x,lat}, pour tous x
    alpha=90._wr-lat
    ! -------                                                  --------    .
    call dist2c(Y,90.0_wr,X,0.0_wr,nR,alpha,La1,Lo1,La2,Lo2)
    ! -----------------------------------------------------------------    .
    if ((abs(Lo1)+abs(La1)+abs(Lo2)+abs(La2)).lt.1000._wr) then            ! une ou deux intersections
      vmini=min(Lo1,Lo2,vOld)
      vmaxi=max(Lo1,Lo2,vOld)
      ! ---------------------------------------------------------------    .
      if(isnan(vmini))stop
      if(isnan(vmaxi))stop
      ! ---------------------------------------------------------------    .
      vNew = normal(vOld,vecartype)                                        ! tirage aléatoire
      ! -------                                                --------    .
      i=0
      do while ((vNew.ge.vmaxi).or.(vNew.le.vmini))                        ! tant que hors des bornes
        i=i+1
        ! -------                                              --------    . réflexion supérieure
        if (vNew.ge.vmaxi) then
          vNew=vNew-2.0_wr*abs(vNew-vmaxi)
        end if
        ! -------                                              --------    . réflexion inférieure
        if (vNew.le.vmini) then
          vNew=vNew+2.0_wr*abs(vmini-vNew)
        end if
        ! -------                                              --------    .
        if(i.eq.200) then
          write(*,*)'avertissement dans tirage_norm_lon : tirage(param) -> iter > 100', &
          vNew,vOld,vmini,vmaxi,vecartype
          vNew = normal(vOld,vecartype)
        endif
        ! -------                                              --------    .
        if(i.gt.250) then
          if (abs(vmini-vmaxi).lt.0.0007_wr) then                          ! quelques mètres de différence  (~50m)
            write(*,*)'avertissement dans tirage_norm_lon : quelques mètres de différence (~50m)',vmini,vmaxi
            vNew = vOld
          else
            write(*,*)'problème dans tirage_norm_lon : problème tirage_norm', &
              vNew,vOld,vmini,vmaxi,vecartype
            stop
          endif
        endif
        ! -------                                              --------    . fin tant que
      enddo
    ! -----------------------------------------------------------------    .
    else                                                                   ! pas d'intersections
      vNew=vOld
      x_0=x_0+1
      if(x_0.gt.10*nbseismes) then
        write(*,*)'problème dans tirage_norm_lon : pas d''intersection'
        write(*,*)vNew,vOld,vmini,vmaxi
        write(*,*)Y,90.0_wr,X,0.0_wr,nR,alpha,La1,Lo1,La2,Lo2
      elseif(x_0.gt.250*nbseismes) then
        write(*,*)'problème dans tirage_norm_lon : pas d''intersection > 250 fois'
        stop
      endif
    endif
    ! -----------------------------------------------------------------    .
  end subroutine tirage_norm_lon

! ---------------------------------------------------------------------    .

  subroutine tirage_norm_lat(vNew,vOld,vmini,vmaxi,vecartype,lon,X,Y,R)
    ! -------                                                  --------    .mh
    ! tirage de la latitude selon une loi normale bornée -> le prior de l'épicentre est ici un cercle !
    ! -------                                                  --------    .
    use dist_cercle
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr),intent(inout) :: vOld,vmini,vmaxi,vecartype
    real(KIND=wr),intent(in) :: X,Y,R,lon
    real(KIND=wr),intent(out) :: vNew
    integer(KIND=wi), save :: x_0=0
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    real(KIND=wr) :: nR,Lo1,Lo2,La1,La2,phi
    ! -----------------------------------------------------------------    . calcul des bornes min et max
    ! petit cercle : centre(X,Y) et rayon (R, en km)
    nR=180.0_wr*R/pi/rT ! deg
    ! grand cercle : de longitude (lon) : 2pts {lon,x}, pour tous x
    ! -------                                                  --------    .
    phi=lon+90.0_wr
    call dist2c(0.0_wr,Y,phi,X,90.0_wr,nr,La1,Lo1,La2,Lo2)
    ! -----------------------------------------------------------------    .
    if ((abs(Lo1)+abs(La1)+abs(Lo2)+abs(La2)).lt.1000._wr) then            ! une ou deux intersections
      vmini=min(La1,La2,vOld)
      vmaxi=max(La1,La2,vOld)
      ! ---------------------------------------------------------------    .
      if(isnan(vmini))stop
      if(isnan(vmaxi))stop
      ! ---------------------------------------------------------------    .
      vNew = normal(vOld,vecartype)                                        ! tirage aléatoire
      ! -------                                                --------    .
      i=0
      do while ((vNew.ge.vmaxi).or.(vNew.le.vmini))                        ! tant que hors des bornes
        i=i+1
        ! -------                                              --------    . réflexion supérieure
        if (vNew.ge.vmaxi) then
          vNew=vNew-2.0_wr*abs(vNew-vmaxi)
        end if
        ! -------                                              --------    . réflexion inférieure
        if (vNew.le.vmini) then
          vNew=vNew+2.0_wr*abs(vmini-vNew)
        end if
        ! -------                                              --------    .
        if(i.eq.200) then
          write(*,*)'avertissement dans tirage_norm_lat : tirage(param) -> iter > 250',vNew,vOld,vmini,vmaxi,vecartype
          vNew = normal(vOld,vecartype)
        endif
        ! -------                                              --------    .
        if(i.gt.250) then
          if (abs(vmini-vmaxi).lt.0.0005_wr) then                          ! quelques mètres de différence (~50m)
            write(*,*)'avertissement dans tirage_norm_lat : quelques mètres de différence (~50m)',vmini,vmaxi
            vNew = vOld
          else
            write(*,*)'problème dans tirage_norm_lat : problème tirage_norm',vNew,vOld,vmini,vmaxi,vecartype
            stop
          endif
        endif
        ! -------                                              --------    . fin tant que
      enddo
    ! -----------------------------------------------------------------    .
    else                                                                   ! pas d'intersections
      vNew=vOld
      x_0=x_0+1
      if(x_0.gt.10*nbseismes) then
        write(*,*)'problème dans tirage_norm_lat : pas d''intersection'
        write(*,*)vNew,vOld,vmini,vmaxi
        write(*,*)0.0_wr,Y,phi,X,90.0_wr,nr,La1,Lo1,La2,Lo2
      elseif(x_0.gt.250*nbseismes) then
        write(*,*)'problème dans tirage_norm_lat : pas d''intersection > 100 fois'
        stop
      endif
    endif
    ! -----------------------------------------------------------------    .
  end subroutine tirage_norm_lat

! ---------------------------------------------------------------------    .

  subroutine tirage_inv_norm(vNew,vOld,vmini,vmaxi,vecartype)
    ! -------                                                  --------    .mh
    ! tirage de l'inverse d'un paramètre selon une densité de probabilité de type loi normale
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr),intent(in) :: vOld,vmini,vmaxi,vecartype
    real(KIND=wr),intent(out) :: vNew
    ! -------                                                  --------    .
    real(KIND=wr) :: ivOld,ivmini,ivmaxi,ivecartype
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    ! -----------------------------------------------------------------    .
    ivOld=1.0_wr/vOld
    ivmini=1.0_wr/vmaxi
    ivmaxi=1.0_wr/vmini
    ivecartype=vecartype/vOld/vOld
    vNew = 1.0_wr/normal(ivOld,ivecartype)
    ! -------                                                  --------    .
    i=0
    do while ((vNew.ge.vmaxi).or.(vNew.le.vmini))                          ! tant que hors des bornes
      i=i+1
      ! -------                                                --------    . réflexion supérieure
      if (vNew.ge.vmaxi) then
        vNew=vNew-2.0_wr*abs(vNew-vmaxi)
      end if
      ! -------                                                --------    . réflexion inférieure
      if (vNew.le.vmini) then
        vNew=vNew+2.0_wr*abs(vmini-vNew)
      end if
      ! -------                                                --------    .
      if(i.eq.200) then
        write(*,*)'avertissement dans tirage_inv_norm  : tirage_log(param) -> iter > 100',vNew
        vNew = 1.0_wr/normal(ivOld,ivecartype)
      endif
      ! -------                                                --------    .
      if(i.gt.250) then
        write(*,*)'problème dans tirage_inv_norm : tirage_log_norm',vNew
        stop
      endif
      ! -------                                                --------    . fin tant que
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine tirage_inv_norm

! ---------------------------------------------------------------------    .

  subroutine tirage_log_norm(vNew,vOld,vmini,vmaxi,vecartype)
    ! -------                                                  --------    .mh
    ! tirage d'un paramètre selon une densité de probabilité identique pour un paramètre x et son inverse 1/x
    ! attention : ce n'est ni un tirage type loi normale, ni un tirage type loi log-normale
    ! log -> ln -> logarithme naturel (à base e = 2,71828... )
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr),intent(in) :: vOld,vmini,vmaxi,vecartype
    real(KIND=wr),intent(out) :: vNew
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    real(KIND=wr) :: sigma
    ! -----------------------------------------------------------------    .
    sigma=log((vOld+vecartype)/vOld)
    vNew = exp(normal(log(vOld),sigma))                                    ! pour une distribution normale
    ! -------                                                  --------    .
    i=0
    do while ((vNew.ge.vmaxi).or.(vNew.le.vmini))                          ! tant que hors des bornes
      i=i+1
      ! -------                                                --------    . réflexion supérieure
      if (vNew.ge.vmaxi) then
        vNew=vNew-2.0_wr*abs(vNew-vmaxi)
      end if
      ! -------                                                --------    . réflexion inférieure
      if (vNew.le.vmini) then
        vNew=vNew+2.0_wr*abs(vmini-vNew)
      end if
      ! -------                                                --------    .
      if(i.eq.200) then
        write(*,*)'avertissement dans tirage_log_norm  : tirage_log(param) -> iter > 100',vNew
        vNew = exp(normal(log(vOld),sigma))
      endif
      ! -------                                                --------    .
      if(i.gt.250) then
        write(*,*)'problème dans tirage_log_norm : tirage_log_norm',vNew
        stop
      endif
      ! -------                                                --------    . fin tant que
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine tirage_log_norm

! ---------------------------------------------------------------------    .

  subroutine tirage_norm_time(vNew,vOld,vmini,vmaxi,vecartype)
    ! -------                                                  --------    .mh
    ! tirage d'un paramètre spécial (date) selon une loi normale
    ! -------                                                  --------    .
    use time
    implicit none
    type(date_sec),intent(in) :: vOld,vmini,vmaxi,vecartype
    type(date_sec),intent(out) :: vNew
    ! -------                                                  --------    .
    type(date_sec) :: vnull
    integer(KIND=wi) :: i
    real(KIND=wr) :: delta,deltamin,deltamax
    ! -----------------------------------------------------------------    .
    vNew=vOld                                                              ! on conserve la même date
    ! -------                                                  --------    .
    call tempszero(vnull%date)
    vnull%sec=0.0_wr
    call difftime(delta,vnull,vecartype)                                   ! écartype en sec
    ! -------                                                  --------    .
    vNew%sec = normal(vOld%sec,abs(delta))                                 ! tirage aléatoire
    ! -------                                                  --------    .
    call basetime(vNew)                                                    ! reste en base 60/12/365 ...
    ! -------                                                  --------    .
    i=0
    call difftime(deltamin,vNew,vmini)
    call difftime(deltamax,vNew,vmaxi)
    do while ((deltamax.ge.0.0_wr).or.(deltamin.le.0.0_wr))                ! tant que hors des bornes
      i=i+1
      ! -------                                                --------    . réflexion supérieure
      if (deltamax.ge.0.0_wr) then
        vNew%sec=vNew%sec-2.0_wr*abs(deltamax)
        call basetime(vNew)                                                ! reste en base 60/12/365 ...
        call difftime(deltamax,vNew,vmaxi)
      end if
          ! -------                                            --------    . réflexion inférieure
      if (deltamin.le.0.0_wr) then
        vNew%sec=vNew%sec+2.0_wr*abs(deltamin)
        call basetime(vNew)                                                ! reste en base 60/12/365 ...
        call difftime(deltamin,vNew,vmini)
      end if
      ! -------                                                --------    .
      if(i.eq.10) then
        write(*,*)'avertissement dans tirage_norm_time : tirage(date) -> iter > 10',vNew%sec
        vNew=vOld
        vNew%sec = normal(vOld%sec,vecartype%sec)
        call basetime(vNew)                                                ! reste en base 60/12/365 ...
      endif
      ! -------                                                --------    .
      if(i.gt.50) then
        write(*,*)'problème dans tirage_norm_time : tirage_norm_time'
        stop
      endif
      ! -------                                                --------    . fin tant que
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine tirage_norm_time

! ---------------------------------------------------------------------    .

  subroutine tirage_H(p,i,all)
    ! -------                                                  --------    .mh
    ! tirage des tous les paramètres Hypocentraux simultanément ou un seul pour un séisme
    ! les paramètres sont libres ou fixes (FLAGhypofixe)
    ! -------                                                  --------    .
    ! all = true : tire tous les parametres
    ! all = false : tire au hasard, un des parametres
    ! -------                                                  --------    . (FLAGhypofixe : MOD/modparam.f90)
    ! FLAGhypofixe = true : tirage aléatoire selon une loi normale dont la moyenne est fixe
    ! FLAGhypofixe = false : tirage aléatoire selon une loi normale 
    !   dont la moyenne correspond à l'iteration précédente (vrai McMC)
    ! -------                                                  --------    .
    use time
    ! -------                                                  --------    .
    implicit none
    type(parametresinv), intent(inout) :: p
    integer(KIND=wi), intent(in) :: i
    logical, intent(in) :: all
    ! -------                                                  --------    .
    type(date_sec) :: tpsref
    ! -------                                                  --------    .
    real(KIND=wr) :: val, val2
    integer(KIND=wi) :: pourcentage
    logical :: P1,P2,P3,P4
    ! -----------------------------------------------------------------    . tous les parametres ou un seul ?
    if (all) then
      P1=.true.
      P2=.true.
      P3=.true.
      P4=.true.
    else
      pourcentage=int(genrand_real1()*100._wr)                             ! aléatoire de 0 à 99
      P1=.false.
      P2=.false.
      P3=.false.
      P4=.false.
      select case (pourcentage)
        case (0:24)
          P1=.true.
        case (25:49)
          P2=.true.
        case (50:74)
          P3=.true.
        case (75:99)
          P4=.true.
      end select
    endif
    ! -----------------------------------------------------------------    .
    if (FLAGhypofixe) then
      ! ---------------------------------------------------------------    . paramètres hypocentraux fixes
      val=p%mini%Lat(i)+(p%maxi%Lat(i)-p%mini%Lat(i))/2.0_wr
      if (P1) call tirage_norm(p%valNew%Lat(i),val,p%mini%Lat(i),p%maxi%Lat(i),p%ecartype%Lat(i))
      ! -------                                                --------    .
      val=p%mini%lon(i)+(p%maxi%lon(i)-p%mini%lon(i))/2.0_wr
      if (P2) call tirage_norm(p%valNew%lon(i),val,p%mini%lon(i),p%maxi%lon(i),p%ecartype%lon(i))
      ! -------                                                --------    .
      tpsref=p%mini%Tzero(i)
      tpsref%sec=tpsref%sec+3.0_wr*p%ecartype%Tzero(i)%sec
      if (P3) call tirage_norm_time(p%valNew%Tzero(i),tpsref,p%mini%Tzero(i),p%maxi%Tzero(i),p%ecartype%Tzero(i))
      ! -------                                                --------    .
      val2=p%maxi%Zhypo(i)
      if(p%maxi%Zhypo(i).gt.(p%valNew%Zmoho-0.1_wr)) then                  ! hypocentre dans la croûte
        val2=p%valNew%Zmoho-0.1_wr
      endif
      val=p%mini%Zhypo(i)+(p%maxi%Zhypo(i)-p%mini%Zhypo(i))/2.0_wr
      if (P4) call tirage_norm(p%valNew%Zhypo(i),val,p%mini%Zhypo(i),val2,p%ecartype%Zhypo(i))
      ! ---------------------------------------------------------------    .
    else
      ! ---------------------------------------------------------------    . paramètres hypocentraux libres
      if (P1) call tirage_norm_time(p%valNew%Tzero(i),p%valOld%Tzero(i),p%mini%Tzero(i),p%maxi%Tzero(i),p%ecartype%Tzero(i))
      ! -------                                                --------    .
      if (P2) call tirage_norm_lat(p%valNew%Lat(i),p%valOld%Lat(i),p%mini%Lat(i),p%maxi%Lat(i), &
        p%ecartype%Lat(i),p%valNew%Lon(i),p%centreX(i),p%centreY(i),p%Rayon(i))
      ! -------                                                --------    .
      if (P3) call tirage_norm_lon(p%valNew%Lon(i),p%valOld%Lon(i),p%mini%Lon(i),p%maxi%Lon(i), &
        p%ecartype%Lon(i),p%valNew%Lat(i),p%centreX(i),p%centreY(i),p%Rayon(i),p%ec_horizontal)
      ! -------                                                --------    .
      val2=p%maxi%Zhypo(i)
      if(p%maxi%Zhypo(i).gt.(p%valNew%Zmoho-0.1_wr)) then                  ! hypocentre dans la croûte
        val2=p%valNew%Zmoho-0.1_wr
      endif
      if (P4) call tirage_norm(p%valNew%Zhypo(i),p%valOld%Zhypo(i),p%mini%Zhypo(i),val2,p%ecartype%Zhypo(i))
      ! ---------------------------------------------------------------    .
    endif
    ! -----------------------------------------------------------------    .
  end subroutine tirage_H

! ---------------------------------------------------------------------    .

  subroutine tirage_T(p,all,vpvs)
    ! -------                                                  --------    .mh
    ! tirage des tous les paramètres de Terre simultanément ou un seul
    ! les paramètres sont libres ou fixes (FLAGterrefixe)
    ! -------                                                  --------    .
    ! all = true : tire tous les parametres
    ! all = false : tire au hasard, un des parametres
    ! vpvs = true : tire ou peut tirer VpVs
    ! vpvs = false : ne peut pas tirer VpVs
    ! -------                                                  --------    . (FLAGterrefixe : MOD/modparam.f90)
    ! FLAGterrefixe = true : tirage aléatoire selon une loi normale dont la moyenne est fixe
    ! FLAGterrefixe = false : tirage aléatoire selon une loi normale
    !   dont la moyenne correspond à l'iteration précédente (vrai McMC)
    ! -------                                                  --------    .
    implicit none
    type(parametresinv), intent(inout) :: p
    logical, intent(in) :: all, vpvs
    ! -------                                                  --------    .
    real(KIND=wr) :: val,val2
    integer(KIND=wi) :: i,pourcentage
    logical :: P1,P2,P3,P4
    ! -----------------------------------------------------------------    . tous les parametres ou un seul ?
    if (all) then
      P1=.true.
      P2=.true.
      P3=.true.
      if (vpvs) then
        P4=.true.
      else
        P4=.false.
      endif
    else
      pourcentage=int(genrand_real1()*100._wr)                             ! aléatoire de 0 à 99
      P1=.false.
      P2=.false.
      P3=.false.
      P4=.false.
      select case (pourcentage)
        !case (0:24)
        case (0:44)
          P1=.true.
        !case (25:49)
        case (45:89)
          P2=.true.
        !case (50:74)
        case (90:94)
          P3=.true.
        !case (75:99)
        case (95:99)
          if (vpvs) then
            P4=.true.
          else
            P4=.false.
            pourcentage=int(genrand_real1()*100._wr)                       ! aléatoire de 0 à 99
            select case (pourcentage)
              case (0:32)
                P1=.true.
              case (33:65)
                P2=.true.
              case (66:99)
                P3=.true.
            end select
          endif
      end select
    endif
    ! -----------------------------------------------------------------    .
    if (FLAGterrefixe) then
      ! ---------------------------------------------------------------    . paramètres de terres fixes
      if (P1) then
        val=p%mini%Zmoho
        do i=1,nbseismes
          if(val.lt.(p%valNew%Zhypo(i)+0.1_wr)) then                       ! hypocentre dans la croûte
            val=p%valNew%Zhypo(i)+0.1_wr
          endif
        enddo
        val2=val+(p%maxi%Zmoho-val)/2.0_wr
        call tirage_norm(p%valNew%Zmoho,val2,val,p%maxi%Zmoho,p%ecartype%Zmoho)
      endif
      ! -------                                                --------    .
      val2=p%maxi%VC
      if(p%maxi%VC.gt.(p%valNew%VM-0.1_wr)) then                           ! respecte VM > VC
        val2=p%valNew%VM-0.1_wr
      endif
      val=p%mini%VC+(p%maxi%VC-p%mini%VC)/2.0_wr
      if (P2) call tirage_log_norm(p%valNew%VC,val,p%mini%VC,val2,p%ecartype%VC)
      ! -------                                                --------    .
      if (P3) then
        val2=p%mini%VM
        if(p%mini%VM.lt.(p%valNew%VC+0.1_wr)) then                         ! respecte VM > VC
          val2=p%valNew%VC+0.1_wr
        endif
        val=p%mini%VM+(p%maxi%VM-p%mini%VM)/2.0_wr
        call tirage_log_norm(p%valNew%VM,val,val2,p%maxi%VM,p%ecartype%VM)
      endif
      ! -------                                                --------    .
      val=p%mini%VpVs+(p%maxi%VpVs-p%mini%VpVs)/2.0_wr
      if (P4) call tirage_norm(p%valNew%VpVs,val,p%mini%VpVs,p%maxi%VpVs,p%ecartype%VpVs)
      ! ---------------------------------------------------------------    .
    else
      ! ---------------------------------------------------------------    . paramètres de terres libres
      if (P1) then
        val=p%mini%Zmoho
        do i=1,nbseismes
          if(val.lt.(p%valNew%Zhypo(i)+0.1_wr)) then                       ! hypocentre dans la croûte
            val=p%valNew%Zhypo(i)+0.1_wr
          endif
        enddo
        call tirage_norm(p%valNew%Zmoho,p%valOld%Zmoho,val,p%maxi%Zmoho,p%ecartype%Zmoho)
      endif
      ! -------                                                --------    .
      if (P2) call tirage_log_norm(p%valNew%VC,p%valOld%VC,p%mini%VC,p%maxi%VC,p%ecartype%VC)
      ! -------                                                --------    .
      if (P3) then
        val=p%mini%VM
        if(p%mini%VM.lt.(p%valNew%VC+0.1_wr)) then                         ! respecte VM > VC
          val=p%valNew%VC+0.1_wr
        endif
        call tirage_log_norm(p%valNew%VM,p%valOld%VM,val,p%maxi%VM,p%ecartype%VM)
      endif
      ! -------                                                --------    .
      if (P4) call tirage_norm(p%valNew%VpVs,p%valOld%VpVs,p%mini%VpVs,p%maxi%VpVs,p%ecartype%VpVs)
    endif
    ! -----------------------------------------------------------------    .
  end subroutine tirage_T

! ---------------------------------------------------------------------    .

  subroutine tirageUN(p,ap1,i)
    ! -------                                                  --------    .mh
    ! tirage d'un seul paramètres : ap(1-8) pour le ième séisme
    ! avec ap = 1->VC,2->VM,3->Zmoho,4->VpVs,5->Lat,6->Lon,7->Zhypo,8->Tzero
    ! -------                                                  --------    .
    implicit none
    type(parametresinv), intent(inout) :: p
    integer(KIND=wi), intent(in) :: ap1                                    ! quel paramètre ?
    integer(KIND=wi), intent(in) :: i                                      ! quel séisme ?
    ! -------                                                  --------    .
    real(KIND=wr) :: val,val2
    integer(KIND=wi) :: j,ap
    ! -----------------------------------------------------------------    .
    if (ap1.gt.8) then
      ap=4+mod(ap1-5,4)+1
    else
      ap=ap1
    endif
    ! -----------------------------------------------------------------    . VC -> 1
    if (ap==1) then
      val=p%maxi%VC
      if(p%maxi%VC.gt.(p%valNew%VM-0.1_wr)) then                           ! respecte VM > VC
        val=p%valNew%VM-0.1_wr
      endif
      call tirage_log_norm(p%valNew%VC,p%valOld%VC,p%mini%VC,val,p%ecartype%VC)
    endif
    ! -----------------------------------------------------------------    . VM -> 2
    if (ap==2) then
      val=p%mini%VM
      if(p%mini%VM.lt.(p%valNew%VC+0.1_wr)) then                           ! respecte VM > VC
        val=p%valNew%VC+0.1_wr
      endif
      call tirage_log_norm(p%valNew%VM,p%valOld%VM,val,p%maxi%VM,p%ecartype%VM)
    endif
    ! -----------------------------------------------------------------    . Zmoho -> 3
    if (ap==3) then
      val=p%mini%Zmoho
      do j=1,nbseismes
        if(val.lt.(p%valNew%Zhypo(j)+0.1_wr)) then                         ! hypocentre dans la croûte
          val=p%valNew%Zhypo(j)+0.1_wr
        endif
      enddo
      call tirage_norm(p%valNew%Zmoho,p%valOld%Zmoho,val,p%maxi%Zmoho,p%ecartype%Zmoho)
    endif
    ! -----------------------------------------------------------------    . VpVs -> 4
    if (ap==4) call tirage_norm(p%valNew%VpVs,p%valOld%VpVs,p%mini%VpVs,p%maxi%VpVs,p%ecartype%VpVs)
    ! -----------------------------------------------------------------    . Lat(i) -> 5
    if (ap==5) call tirage_norm_lat(p%valNew%Lat(i),p%valOld%Lat(i),p%mini%Lat(i),p%maxi%Lat(i), &
      p%ecartype%Lat(i),p%valNew%Lon(i),p%centreX(i),p%centreY(i),p%Rayon(i))
    ! -----------------------------------------------------------------    . Lon(i) -> 6
    if (ap==6) call tirage_norm_lon(p%valNew%Lon(i),p%valOld%Lon(i),p%mini%Lon(i),p%maxi%Lon(i), &
      p%ecartype%Lon(i),p%valNew%Lat(i),p%centreX(i),p%centreY(i),p%Rayon(i),p%ec_horizontal)
    ! -----------------------------------------------------------------    . Zhypo(i) -> 7
    if (ap==7) then
      val2=p%maxi%Zhypo(i)
      if(p%maxi%Zhypo(i).gt.(p%valNew%Zmoho-0.1_wr)) then                  ! hypocentre dans la croûte
        val2=p%valNew%Zmoho-0.1_wr
      endif
      call tirage_norm(p%valNew%Zhypo(i),p%valOld%Zhypo(i),p%mini%Zhypo(i), &
        val2,p%ecartype%Zhypo(i))
    endif
    ! -----------------------------------------------------------------    . Tzero(i) -> 8
    if (ap==8) call tirage_norm_time(p%valNew%Tzero(i),p%valOld%Tzero(i),p%mini%Tzero(i), &
      p%maxi%Tzero(i),p%ecartype%Tzero(i))
    ! -----------------------------------------------------------------    .
  end subroutine tirageUN



END MODULE tirage



! *********************************************************************    .
! *********************************************************************    .


