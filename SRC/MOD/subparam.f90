! Librairie de subroutines concernant les parametres d'inversion
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE sub_param

    use modparam

    implicit none

    private

    public  :: lect_prior
    public  :: moycoldruns
    public  :: init_div
    public  :: calc_accept
    public  :: lectparam
    public  :: nb_mod_selec
    public  :: lect_mod_select
    public  :: moy_mod_select
    public  :: moy_ec, mediane
    public  :: dist_apriori
    public  :: paramfixe
    public  :: inR, outR


CONTAINS

    ! -----------------------------------------------------------------    .

  subroutine lect_prior(param_p,param_init,CorH)
    ! -------                                                  --------    .mh
    ! lecture du prior dans PARAM/priorIn.d
    ! -------                                                  --------    .
    use typetemps, only : parametresinv, parametres
    use time
    ! -------                                                  --------    .
    implicit none
    type(parametresinv), intent(out) :: param_p                            ! paramètres de l'inversion
    type(parametres), intent(in) :: param_init                             ! un jeu de paramètres
    character (LEN=1), intent (in) :: CorH                                 ! cold or Hot runs
    ! -------                                                  --------    .
    integer(KIND=wi) :: i, ok
    real(KIND=wr) :: dist, dist_et                                         ! disance de recherche
    real(KIND=wr), parameter :: dcherche = 350.0_wr                        ! à priori : épicentre est à moins de 350 km de la premiere station !!!
    real(KIND=wr) :: val
    ! -----------------------------------------------------------------    .
    param_p%Rayon=dcherche
    ! -----------------------------------------------------------------    .
    ok=0
    param_p%valNew=param_init
    param_p%valOld=param_init
    param_p%mini=param_init
    param_p%maxi=param_init
    param_p%ecartype=param_init
    ! -------                                                  --------    . min max pour Lon et Lat
    do i=1, nbseismes
      if ((param_p%valNew%Lat(i).gt.89._wr).or. &
        (param_p%valNew%Lat(i).lt.-89._wr).or. &
        (param_p%valNew%Lon(i).gt.179._wr).or. &
        (param_p%valNew%Lon(i).lt.-179._wr)) then
        write(*,*)'problème dans lect_prior 1 : longitude ou latitude actuelle trop près des coutures !'
        write(*,*)'pas prévu',i,param_p%valNew%Lat(i),param_p%valNew%Lon(i)
        stop
      endif
    enddo
    ! -------                                                  --------    .
    ok=0
    if (CorH.eq."C") then
      open(980, FILE = 'PARAM/priorIn_COLD.d',status='old',iostat = ok)
    elseif (CorH.eq."H") then
      open(980, FILE = 'PARAM/priorIn_HOT.d',status='old',iostat = ok)
    else
      write(*,*)'problème dans lect_prior : ni Coldruns ni Hotruns'
    endif
    if (ok .ne. 0) then
      write(*,*)'problème dans lect_prior : le fichier PARAM/priorIn.d n''existe pas '
       stop
    endif
    ! -------                                                  --------    . PRIOR pour le modèle de terre
    read(980,*)param_p%mini%Vc,param_p%maxi%Vc,param_p%ecartype%Vc         ! vitesse croute (en km/s)
    read(980,*)param_p%mini%Vm,param_p%maxi%Vm,param_p%ecartype%Vm         ! vitesse manteau (en km/s)
    read(980,*)param_p%mini%Zmoho,param_p%maxi%Zmoho,param_p%ecartype%Zmoho ! moho (km), positif vers le bas, 0 = niveau de la mer
    read(980,*)param_p%mini%VpVs,param_p%maxi%VpVs,param_p%ecartype%VpVs   ! Vp/Vs (sans unités), croute et manteau
    ! profondeur de hypoventre (km), positif vers le bas, négatif si relief, 0 = niveau de la mer
    read(980,*)param_p%mini%Zhypo(1),param_p%maxi%Zhypo(1),param_p%ecartype%Zhypo(1)
    do i=2, nbseismes
        param_p%mini%Zhypo(i)=param_p%mini%Zhypo(1)
        param_p%maxi%Zhypo(i)= param_p%maxi%Zhypo(1)
        param_p%ecartype%Zhypo(i)=param_p%ecartype%Zhypo(1)
    enddo
    ! -------                                                  --------    . séisme dans la croute
    if ((param_p%maxi%Zmoho+0.1_wr).le.param_p%maxi%Zhypo(1)) then
      param_p%maxi%Zmoho=param_p%maxi%Zhypo(1)+0.1_wr
    endif
    ! -------                                                  --------    .
    read(980,*)dist_et                                                     ! distance (km) correspondant à l'écart type pour Lon et Lat
    param_p%ec_horizontal=dist_et
    ! -------                                                  --------    .
    call tempszero(param_p%ecartype%Tzero(1)%date)
    read(980,*)param_p%ecartype%Tzero(1)%sec                               ! écart type (secondes) pour le temps initial
    call basetime(param_p%ecartype%Tzero(1))                               ! reste en base 60/12/365 ...
    close(980)
    ! -------                                                  --------    . on cherche la solution dans un rayon de 'dcherche' km (pas au-delà !!!!)
    do i=1, nbseismes
      dist = dcherche * 360.0_wr / ( 2.0_wr * pi * rT)
      param_p%mini%Lat(i) = param_init%Lat(i) - dist
      param_p%maxi%Lat(i) = param_init%Lat(i) + dist
      dist = dcherche * 360.0_wr / ( 2.0_wr * pi * rT) / &
        sin((90.0_wr-param_init%Lat(i))/180.0_wr*pi)
      param_p%mini%Lon(i) = param_init%Lon(i) - dist
      param_p%maxi%Lon(i) = param_init%Lon(i) + dist
      ! -------                                                --------    . écart-types lon et lat homogènes (lu en km : dist_et)
      param_p%ecartype%Lat(i) = dist_et * 360.0_wr / ( 2.0_wr * pi * rT)
      param_p%ecartype%Lon(i) = param_p%ecartype%Lat(i) /  &
        sin((90.0_wr-param_init%Lat(i))/180.0_wr*pi)
      ! -------                                                --------    .
      if ((param_p%maxi%Lat(i).gt.89.9_wr).or. &
            (param_p%mini%Lat(i).lt.-89.9_wr).or. &
            (param_p%maxi%Lon(i).gt.179.9_wr).or. &
            (param_p%mini%Lon(i).lt.-179.9_wr)) then
        write(*,*)'problème dans lect_prior 2 : longitude ou latitude min/max trop près des coutures !'
        write(*,*)'pas prévu',i,param_p%valNew%Lat(i),param_p%valNew%Lon(i)
        stop
      endif
    enddo
    ! -------                                                  --------    . min max pour le temps initial
    do i=1, nbseismes
      param_p%maxi%Tzero(i)%sec = param_init%Tzero(i)%sec + 1.0_wr * 60.0_wr ! plus ou moins 1 minutes
      call basetime(param_p%maxi%Tzero(i))                                 ! reste en base 60/12/365 ...
      param_p%mini%Tzero(i)%sec = param_init%Tzero(i)%sec - 1.0_wr * 60.0_wr ! plus ou moins 1 minutes
      call basetime(param_p%mini%Tzero(i))                                 ! reste en base 60/12/365 ...
      ! -------                                                --------    . au cas où ...
      call basetime(param_p%valNew%Tzero(i))                               ! reste en base 60/12/365 ...
      call basetime(param_p%valOld%Tzero(i))                               ! reste en base 60/12/365 ...
      param_p%ecartype%Tzero(i)=param_p%ecartype%Tzero(1)
    enddo
    ! -------                                                  --------    . verification au cas ou ...
    ! coldruns appartient bien au nouveau prior ? :
    if (CorH.eq."H") then
      do i=1, nbseismes
        if (param_p%valNew%lon(i).le.param_p%mini%lon(i)) param_p%valNew%lon(i)=param_p%mini%lon(i)
        if (param_p%valNew%lon(i).ge.param_p%maxi%lon(i)) param_p%valNew%lon(i)=param_p%maxi%lon(i)
        if (param_p%valNew%lat(i).le.param_p%mini%lat(i)) param_p%valNew%lat(i)=param_p%mini%lat(i)
        if (param_p%valNew%lat(i).ge.param_p%maxi%lat(i)) param_p%valNew%lat(i)=param_p%maxi%lat(i)
        if (param_p%valNew%Zhypo(i).le.param_p%mini%Zhypo(i)) param_p%valNew%Zhypo(i)=param_p%mini%Zhypo(i)
        if (param_p%valNew%Zhypo(i).ge.param_p%maxi%Zhypo(i)) param_p%valNew%Zhypo(i)=param_p%maxi%Zhypo(i)
        call difftime(val,param_p%valNew%Tzero(i),param_p%maxi%Tzero(i))
        if (val.gt.0.0_wr) param_p%valNew%Tzero(i)=param_p%maxi%Tzero(i)
        call difftime(val,param_p%valNew%Tzero(i),param_p%mini%Tzero(i))
        if (val.lt.0.0_wr) param_p%valNew%Tzero(i)=param_p%mini%Tzero(i)
      enddo
      if (param_p%valNew%VC.le.param_p%mini%VC) param_p%valNew%VC=param_p%mini%VC
      if (param_p%valNew%VC.ge.param_p%maxi%VC) param_p%valNew%VC=param_p%maxi%VC
      if (param_p%valNew%VM.le.param_p%mini%VM) param_p%valNew%VM=param_p%mini%VM
      if (param_p%valNew%VM.ge.param_p%maxi%VM) param_p%valNew%VM=param_p%maxi%VM
      if (param_p%valNew%Zmoho.le.param_p%mini%Zmoho) param_p%valNew%Zmoho=param_p%mini%Zmoho
      if (param_p%valNew%Zmoho.ge.param_p%maxi%Zmoho) param_p%valNew%Zmoho=param_p%maxi%Zmoho
      if (param_p%valNew%VpVs.le.param_p%mini%VpVs) param_p%valNew%VpVs=param_p%mini%VpVs
      if (param_p%valNew%VpVs.ge.param_p%maxi%VpVs) param_p%valNew%VpVs=param_p%maxi%VpVs
      param_p%valOld=param_p%valNew
      ! -------                                                --------    .
      ok=0
      if (param_init%VC.le.param_p%mini%VC) ok=-1
      if (param_init%VC.ge.param_p%maxi%VC) ok=-1
      if (param_init%VM.le.param_p%mini%VM) ok=-1
      if (param_init%VM.ge.param_p%maxi%VM) ok=-1
      if (param_init%Zmoho.le.param_p%mini%Zmoho) ok=-1
      if (param_init%Zmoho.ge.param_p%maxi%Zmoho) ok=-1
      if (param_init%VpVs.le.param_p%mini%VpVs) ok=-1
      if (param_init%VpVs.ge.param_p%maxi%VpVs) ok=-1
      do i=1,nbseismes
        if (param_init%Lon(i).le.param_p%mini%Lon(i)) ok=-1
        if (param_init%Lon(i).ge.param_p%maxi%Lon(i)) ok=-1
        if (param_init%Lat(i).le.param_p%mini%Lat(i)) ok=-1
        if (param_init%Lat(i).ge.param_p%maxi%Lat(i)) ok=-1
        if (param_init%Zhypo(i).le.param_p%mini%Zhypo(i)) ok=-1
        if (param_init%Zhypo(i).ge.param_p%maxi%Zhypo(i)) ok=-1
        call difftime(val,param_init%Tzero(i),param_p%maxi%Tzero(i))
        if (val.gt.0.0_wr) ok=-1
        call difftime(val,param_init%Tzero(i),param_p%mini%Tzero(i))
        if (val.lt.0.0_wr) ok=-1
      enddo
      ! -------                                                --------    .
      if ((ok==-1).and.(.not.FLAGterrefixe)) then
        write(*,*)'problème dans lect_prior : rectifier les priors pour hotruns'
        write(*,*) '(trop resserré par rapport aux coldruns)'
        write(*,*)
        write(*,*)param_init
        write(*,*)param_p%mini
        write(*,*)param_p%maxi
        stop
      endif
      ! -------                                                --------    .
    endif
    ! -----------------------------------------------------------------    .
  end subroutine lect_prior

    ! -----------------------------------------------------------------    .

  subroutine paramfixe(param_p)
    ! -------                                                  --------    .mh
    ! lecture des paramtres fixe si nécessaires
    ! -------                                                  --------    .
    use typetemps, only : parametresinv, date_sec
    use time
    ! -------                                                  --------    .
    implicit none
    type(parametresinv), intent(inout) :: param_p                          ! paramètres de l'inversion
    ! -------                                                  --------    .
    type(date_sec) :: tpsref,tps
    real(KIND=wr) :: moy, ec
    integer(KIND=wi) :: i
    integer(KIND=wi) :: ok
    ! -----------------------------------------------------------------    . lecture des paramètres de terre si fixes
    ok=0
    if(FLAGterre) then
      open(unit=50,file="PARAM/paramTerre.d",STATUS="old",iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans paramfixe : le fichier PARAM/paramTerre.d n''existe pas '
        stop
      endif
      read(50,*)moy,ec
      param_p%valNew%VC=moy
      param_p%valOld%VC=moy
      param_p%mini%VC=moy-3.0_wr*ec
      param_p%maxi%VC=moy+3.0_wr*ec
      param_p%ecartype%VC=ec
      read(50,*)moy,ec
      param_p%valNew%VM=moy
      param_p%valOld%VM=moy
      param_p%mini%VM=moy-3.0_wr*ec
      param_p%maxi%VM=moy+3.0_wr*ec
      param_p%ecartype%VM=ec
      read(50,*)moy,ec
      param_p%valNew%Zmoho=moy
      param_p%valOld%Zmoho=moy
      param_p%mini%Zmoho=moy-3.0_wr*ec
      param_p%maxi%Zmoho=moy+3.0_wr*ec
      param_p%ecartype%Zmoho=ec
      read(50,*)moy,ec
      param_p%valNew%VpVs=moy
      param_p%valOld%VpVs=moy
      param_p%mini%VpVs=moy-3.0_wr*ec
      param_p%maxi%VpVs=moy+3.0_wr*ec
      param_p%ecartype%VpVs=ec
      close(50)
    endif
    ! -----------------------------------------------------------------    . lecture des paramètres hypocentraux si fixes
    ok=0
    if(FLAGhypo) then
      open(unit=51,file="PARAM/paramHypo.d",STATUS="old",iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans paramfixe : le fichier PARAM/paramHypo.d n''existe pas '
        stop
      endif
      do i=1,nbseismes
        read(51,*)moy,ec
        param_p%valNew%lon(i)=moy
        param_p%valOld%lon(i)=moy
        param_p%mini%lon(i)=moy-3.0_wr*ec
        param_p%maxi%lon(i)=moy+3.0_wr*ec
        param_p%ecartype%lon(i)=ec
        read(51,*)moy,ec
        param_p%valNew%lat(i)=moy
        param_p%valOld%lat(i)=moy
        param_p%mini%lat(i)=moy-3.0_wr*ec
        param_p%maxi%lat(i)=moy+3.0_wr*ec
        param_p%ecartype%lat(i)=ec
        read(51,*)moy,ec
        param_p%valNew%Zhypo(i)=moy
        param_p%valOld%Zhypo(i)=moy
        param_p%mini%Zhypo(i)=moy-3.0_wr*ec
        param_p%maxi%Zhypo(i)=moy+3.0_wr*ec
        param_p%ecartype%Zhypo(i)=ec
        read(51,*)tpsref
        read(51,*)moy,ec
        tpsref%sec=moy
        call basetime(tpsref)
        param_p%valNew%Tzero(i)=tpsref
        param_p%valOld%Tzero(i)=tpsref
        tps=tpsref
        tps%sec=tps%sec-3.0_wr*ec
        param_p%mini%Tzero(i)=tps
        tps=tpsref
        tps%sec=tps%sec+3.0_wr*ec
        param_p%maxi%Tzero(i)=tps
        call tempszero(tps%date)
        tps%sec=ec
        param_p%ecartype%Tzero(i)=tps
      enddo
      close(51)
    endif
    ! -----------------------------------------------------------------    .
  end subroutine paramfixe

    ! -----------------------------------------------------------------    .

  subroutine moycoldruns(nbChaineMV,param_best,misfit,nbChaineMVhot,dc)
    ! -------                                                  --------    .mh
    ! fait la moyenne (et écart-type) des meilleurs modèles pour les coldruns
    ! -------                                                  --------    .
    use typetemps, only : parametres, fcout, coldmoy, coldmoyval
    use time
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: nbChaineMV                            ! nombre de coldruns
    integer(KIND=wi), intent (in) :: nbChaineMVhot                         ! nombre de hotrun
    type(parametres), intent (in) :: param_best(nbChaineMV)                ! triés par misfit (triparam)
    type(fcout), intent (in) :: misfit(nbChaineMV)
    type(coldmoy), intent (out) :: dc                                      ! moyennes et écarts-types des modèles du coldrun
    ! -------                                                  --------    .
    type(coldmoyval) :: cval                                               ! modèles du coldrun
    integer(KIND=wi) :: i, j
    real(KIND=wr) :: vecT(nbChaineMV)
    real(KIND=wr) :: vecS(nbChaineMVhot)
    ! -----------------------------------------------------------------    . allocation dynamique
    allocate(cval%Tmis(nbChaineMV),cval%TVC(nbChaineMV),cval%TVM(nbChaineMV), &
    cval%TZmoho(nbChaineMV),cval%TVpVs(nbChaineMV),cval%TLat(nbChaineMV,nbseismes), &
    cval%TLon(nbChaineMV,nbseismes),cval%TZhypo(nbChaineMV,nbseismes), &
    cval%TTzero(nbChaineMV,nbseismes))
    allocate(cval%Smis(nbChaineMVhot),cval%SVC(nbChaineMVhot), &
    cval%SVM(nbChaineMVhot),cval%SZmoho(nbChaineMVhot), &
    cval%SVpVs(nbChaineMVhot),cval%SLat(nbChaineMVhot,nbseismes), &
    cval%SLon(nbChaineMVhot,nbseismes),cval%SZhypo(nbChaineMVhot,nbseismes), &
    cval%STzero(nbChaineMVhot,nbseismes))
    ! -------                                                  --------    .
    do i=1, nbseismes
      dc%tempsrefcold(i) = param_best(1)%Tzero(i)
      dc%tempsrefcold(i)%sec = 0.0_wr
      call tempszero(dc%moytot%par%Tzero(i)%date)
      call tempszero(dc%ectot%par%Tzero(i)%date)
      call tempszero(dc%moyselect%par%Tzero(i)%date)
      call tempszero(dc%ecselect%par%Tzero(i)%date)
    enddo
    ! -------                                                  --------    . toutes les chaînes
    do i=1, nbChaineMV
      cval%Tmis(i) = misfit(i)%best
      cval%TVC(i) = param_best(i)%VC
      cval%TVM(i) = param_best(i)%VM
      cval%TZmoho(i) = param_best(i)%Zmoho
      cval%TVpVs(i) = param_best(i)%VpVs
      do j=1, nbseismes
        cval%TLat(i,j) = param_best(i)%Lat(j)
        cval%TLon(i,j) = param_best(i)%Lon(j)
        cval%TZhypo(i,j) = param_best(i)%Zhypo(j)
        call difftime(cval%TTzero(i,j),param_best(i)%Tzero(j),dc%tempsrefcold(j))
      enddo
    enddo
    ! -------                                                  --------    . les chaînes sélectionnées
    do i =1,nbChaineMVhot
      cval%Smis(i) = misfit(i)%best
      cval%SVC(i) = param_best(i)%VC
      cval%SVM(i) = param_best(i)%VM
      cval%SZmoho(i) = param_best(i)%Zmoho
      cval%SVpVs(i) = param_best(i)%VpVs
      do j=1, nbseismes
        cval%SLat(i,j) = param_best(i)%Lat(j)
        cval%SLon(i,j) = param_best(i)%Lon(j)
        cval%SZhypo(i,j) = param_best(i)%Zhypo(j)
        call difftime(cval%STzero(i,j),param_best(i)%Tzero(j),dc%tempsrefcold(j))
      enddo
    enddo
    ! -------                                                  --------    . calcul moyenne pour tous les meilleurs modeles des coldruns
    call moy_ec(cval%Tmis,nbChaineMV,nbChaineMV,dc%moytot%mis,dc%ectot%mis)
    call moy_ec(cval%TVC,nbChaineMV,nbChaineMV,dc%moytot%par%VC,dc%ectot%par%VC)
    call moy_ec(cval%TVM,nbChaineMV,nbChaineMV,dc%moytot%par%VM,dc%ectot%par%VM)
    call moy_ec(cval%TZmoho,nbChaineMV,nbChaineMV,dc%moytot%par%Zmoho,dc%ectot%par%Zmoho)
    call moy_ec(cval%TVpVs,nbChaineMV,nbChaineMV,dc%moytot%par%VpVs,dc%ectot%par%VpVs)
    do j=1, nbseismes
      vecT=cval%TLon(:,j)
      call moy_ec(vecT,nbChaineMV,nbChaineMV,dc%moytot%par%Lon(j),dc%ectot%par%Lon(j))
      vecT=cval%TLat(:,j)
      call moy_ec(vecT,nbChaineMV,nbChaineMV,dc%moytot%par%Lat(j),dc%ectot%par%Lat(j))
      vecT=cval%TZhypo(:,j)
      call moy_ec(vecT,nbChaineMV,nbChaineMV,dc%moytot%par%Zhypo(j),dc%ectot%par%Zhypo(j))
      vecT=cval%TTzero(:,j)
      call moy_ec(vecT,nbChaineMV,nbChaineMV,dc%moytot%par%Tzero(j)%sec,dc%ectot%par%Tzero(j)%sec)
    enddo
    ! -------                                                  --------    . calcul moyenne pour tous les meilleurs modeles des coldruns sélectionnés
    call moy_ec(cval%Smis,nbChaineMVhot,nbChaineMVhot,dc%moyselect%mis,dc%ecselect%mis)
    call moy_ec(cval%SVC,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%VC,dc%ecselect%par%VC)
    call moy_ec(cval%SVM,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%VM,dc%ecselect%par%VM)
    call moy_ec(cval%SZmoho,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%Zmoho,dc%ecselect%par%Zmoho)
    call moy_ec(cval%SVpVs,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%VpVs,dc%ecselect%par%VpVs)
    do j=1, nbseismes
      vecS(:)=cval%SLon(:,j)
      call moy_ec(vecS,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%Lon(j),dc%ecselect%par%Lon(j))
      vecS(:)=cval%SLat(:,j)
      call moy_ec(vecS,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%Lat(j),dc%ecselect%par%Lat(j))
      vecS(:)=cval%SZhypo(:,j)
      call moy_ec(vecS,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%Zhypo(j),dc%ecselect%par%Zhypo(j))
      vecS(:)=cval%STzero(:,j)
      call moy_ec(vecS,nbChaineMVhot,nbChaineMVhot,dc%moyselect%par%Tzero(j)%sec,dc%ecselect%par%Tzero(j)%sec)
    enddo
    ! -------                                                  --------    .
    deallocate(cval%Tmis,cval%TVC,cval%TVM,cval%TZmoho,cval%TVpVs,cval%TLat,cval%TLon, &
    cval%TZhypo,cval%TTzero,cval%Smis,cval%SVC,cval%SVM,cval%SZmoho,cval%SVpVs,cval%SLat, &
    cval%SLon,cval%SZhypo,cval%STzero)
    ! -----------------------------------------------------------------    . vérif
    do i=1,nbseismes
      if (dc%moytot%par%lon(i).gt.179._wr) dc%moytot%par%lon(i)=0.0_wr
      if (dc%moytot%par%lon(i).lt.-179._wr) dc%moytot%par%lon(i)=0.0_wr
      if (dc%ectot%par%lon(i).gt.179._wr) dc%ectot%par%lon(i)=0.0_wr
      if (dc%ectot%par%lon(i).lt.-179._wr) dc%ectot%par%lon(i)=0.0_wr
      if (dc%moytot%par%lat(i).gt.89._wr) dc%moytot%par%lat(i)=0.0_wr
      if (dc%moytot%par%lat(i).lt.-89._wr) dc%moytot%par%lat(i)=0.0_wr
      if (dc%ectot%par%lat(i).gt.89._wr) dc%ectot%par%lat(i)=0.0_wr
      if (dc%ectot%par%lat(i).lt.-89._wr) dc%ectot%par%lat(i)=0.0_wr
      if (dc%moyselect%par%lon(i).gt.179._wr) dc%moyselect%par%lon(i)=0.0_wr
      if (dc%moyselect%par%lon(i).lt.-179._wr) dc%moyselect%par%lon(i)=0.0_wr
      if (dc%ecselect%par%lon(i).gt.179._wr) dc%moyselect%par%lon(i)=0.0_wr
      if (dc%ecselect%par%lon(i).lt.-179._wr) dc%moyselect%par%lon(i)=0.0_wr
      if (dc%moyselect%par%lat(i).gt.89._wr) dc%ecselect%par%lat(i)=0.0_wr
      if (dc%moyselect%par%lat(i).lt.-89._wr) dc%ecselect%par%lat(i)=0.0_wr
      if (dc%ecselect%par%lat(i).gt.89._wr) dc%ecselect%par%lat(i)=0.0_wr
      if (dc%ecselect%par%lat(i).lt.-89._wr) dc%ecselect%par%lat(i)=0.0_wr
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine moycoldruns

    ! -----------------------------------------------------------------    .

  subroutine init_div(misfit,acceptance)
    ! -------                                                  --------    .mh
    ! initialise quelques variables
    ! -------                                                  --------    .
    use typetemps, only : fcout, accept
    ! -------                                                  --------    .
    implicit none
    type(fcout),intent(out) :: misfit                                      !    fonction coût
    type(accept),intent(out) :: acceptance                                 !       acceptance
    ! -----------------------------------------------------------------    .
    misfit%old=100000.0_wr
    misfit%best=100000.0_wr
    acceptance%N=int(0,wl)                                                 ! modèle non accepté
    acceptance%O=int(0,wl)                                                 ! modèle accepté car meilleur que précendant
    acceptance%NO=int(0,wl)                                                ! modèle accepté car repêché par le Metropolis
    acceptance%val=0.0_wr
    ! -----------------------------------------------------------------    .
  end subroutine init_div

    ! -----------------------------------------------------------------

  subroutine calc_accept(acceptance)
    ! -------                                                  --------    .mh
    ! calcul de l'acceptance
    ! -------                                                  --------    .
    use typetemps, only : accept
    ! -------                                                  --------    .
    implicit none
    type(accept),intent(inout) :: acceptance                               !       acceptance (%)
    ! -----------------------------------------------------------------    .
    ! acceptance = ( accepté du premier coup + repêchés ) / total
    acceptance%val = real(acceptance%O+acceptance%NO,wr)/ &
                 real(acceptance%NO+acceptance%N+acceptance%O,wr)*100.0_wr
    ! -----------------------------------------------------------------    .
  end subroutine calc_accept

    ! -----------------------------------------------------------------    .

  subroutine lectparam(nbChaineMV1,nbChaineMV2,maxiter2,maxiter1,chut)
    ! -------                                                  --------    .mh
    ! lecture du nombre de chaîne de Markox (nbChaineMV)
    ! lecture du nombre d'itération par chaîne de Markox (maxiter)
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi),intent(out) :: nbChaineMV1,nbChaineMV2                ! cold puis hot runs
    integer(KIND=wi),intent(out) :: maxiter1,maxiter2                      ! cold puis hot runs
    integer(KIND=wi) :: ok
    logical, intent(in), optional :: chut
    ! -------                                                  --------    .
    logical :: printtest
    ! -----------------------------------------------------------------    .
    ok = 0
    open(unit=981,file="PARAM/iteration.d",STATUS="old",iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans lectparam : le fichier PARAM/iteration.d n''existe pas '
      stop
    endif
    read(981,*)nbChaineMV1,maxiter1
    read(981,*)nbChaineMV2,maxiter2
    close(981)
    ! -------                                                  --------    . verif
    if (nbChaineMV2 .gt. nbChaineMV1) then
        write(*,*)'problème dans lectparam : nbChaineMV hot > nbChaineMV cold '
        stop
    endif
    ! -------                                                  --------    .
    printtest=.true.
    if (present(chut)) then
      if (chut) printtest=.false.
    endif
    ! -------                                                  --------    .
    if(printtest) then
      write(*,*)'nombre de séismes :                    ',nbseismes
      write(*,*)'nombre de chaînes de Markov (cold) :   ',nbChaineMV1
      write(*,*)'nombre de chaînes de Markov (hot) :    ',nbChaineMV2
      write(*,*)'nombre d''itérations par chaîne (cold) :',maxiter1
      write(*,*)'nombre d''itérations par chaîne (hot) : ',maxiter2
      write(*,*)'nombre modèles testés          :       ',maxiter1*nbChaineMV1+maxiter2*nbChaineMV2
    endif
    ! -----------------------------------------------------------------    .
  end subroutine lectparam

    ! -----------------------------------------------------------------    .

  subroutine nb_mod_selec(nbparam,nbChaineMV)
    ! -------                                                  --------    .mh
    ! relecture des modeles sélectionnés (distrib. a posteriori), définition du nombre de modèles sélectionnés après hotruns
    ! -------                                                  --------    .
    use typetemps, only : paramisfit
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(out) :: nbparam
    integer(KIND=wi), intent(in) :: nbChaineMV
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,k,ok,noctet
    type(paramisfit) :: pm
    character (LEN=5) :: numberchaine
    ! -----------------------------------------------------------------    .
    nbparam = 0
    ! -------                                                  --------    .
    do i=1,nbChaineMV
      ok = 0
      write(numberchaine(1:5),'(i5)')i
      inquire ( iolength = noctet ) pm
      open(unit=1000+i,file="OUTPUT/files/"//trim(adjustl(numberchaine))//".bin",STATUS="old",access='direct',RECL=noctet,iostat=ok)
      if (ok .ne. 0) then
        write(*,*)"problème dans nb_mod_selec : le fichier OUTPUT/files/"//trim(adjustl(numberchaine))//".bin n''existe pas "
        stop
      endif
      k=0
      do while(ok .eq. 0)                                                  ! boucle pour compter le nombre de lignes du fichier
        k=k+1
        read(1000+i,REC=k,iostat = ok)pm%mis,pm%par
        if (ok .eq. 0) nbparam = nbparam + 1
      enddo
      close(1000+i)
    enddo
    write(*,*)'nombre de modèles retenus      :       ',nbparam
    ! -----------------------------------------------------------------    .
  end subroutine nb_mod_selec

    ! -----------------------------------------------------------------    .

  subroutine lect_mod_select(p,dp,nbChaineMV,mis,p_best)
    ! -------                                                  --------    .mh
    ! relecture des modeles sélectionnés (distrib. a posteriori) après hotruns
    ! -------                                                  --------    .
    use statistiques
    use typetemps
    use time
    ! -------                                                  --------    .
    implicit none
    type(parametresinv), intent(in) ::  p                                  ! paramètres
    integer(KIND=wi), intent(in) :: nbChaineMV                             ! nombre de chaînes
    type(parametres), intent(in):: p_best(nbChaineMV)                      ! meilleur jeu de parametre pour chaque chaîne (en liens avec "mis")
    type(fcout), intent(in) :: mis(nbChaineMV)                             ! meilleur fonction coût pour chaque chaîne
    type(densityplot), intent(inout) :: dp                                 ! vecteur des modèles sélectionnés
    ! -------                                                  --------    .
    type(paramisfit) :: pm
    integer(KIND=wi) :: ok,i,j,k,l,m,n
    integer(KIND=wi) :: noct
    character (LEN=5) :: numberchaine
    character (LEN=5) :: numberseisme
    type(date_sec) :: temps_mod
    real(KIND=wr) :: minmis,size,gap
    real(KIND=wr) :: minmax_5,minmax_1
    real(KIND=wr) :: minmax_2(nbseismes),minmax_3(nbseismes),minmax_4(nbseismes)
    real(KIND=wr) :: minmax_6(nbseismes),minmax_7(nbseismes),minmax_8(nbseismes)
    character(len=50) :: namefile
    ! -----------------------------------------------------------------    .
    call nb_mod_selec(dp%nbparam,nbChaineMV)                               ! cherche le nombre de modèles sélectionnés
    ! -------                                                  --------    .
    allocate(dp%mis%vec(dp%nbparam))                                       ! allocation dynamique du 'type' dp
    allocate(dp%VC%vec(dp%nbparam))
    allocate(dp%VM%vec(dp%nbparam))
    allocate(dp%Zmoho%vec(dp%nbparam))
    allocate(dp%VpVs%vec(dp%nbparam))
    do i=1,nbseismes
      allocate(dp%Lat(i)%vec(dp%nbparam))
      allocate(dp%Lon(i)%vec(dp%nbparam))
      allocate(dp%Zhypo(i)%vec(dp%nbparam))
      allocate(dp%Tzero(i)%vec(dp%nbparam))
    enddo
    ! -------                                                  --------    .
    ! sélectionne le modèle dont le misfit est le plus bas, afin de définir une réference de temps en Année, Mois, Jour, heure et min
    ! pour la suite du programme, le parametre Temps_zero correspondra à ce temps_ref (min) + dp%Tzero%vec(n) (secondes relatives)
    ! -------                                                  --------    .
    do j=1,nbseismes
      minmis=1.e9_wr
      do i=1,nbChaineMV
        if(mis(i)%best.lt.minmis) then
            minmis=mis(i)%best
            dp%temps_ref(j)=p_best(i)%Tzero(j)
          endif
      enddo
      dp%temps_ref(j)%sec=0.0_wr
      call basetime(dp%temps_ref(j))
    enddo
    ! -------                                                  --------    . lecture des modèles sélectionnés
    n=0
    l=1
    do i=1,nbChaineMV
      ok = 0
      write(numberchaine(1:5),'(i5)')i
      inquire ( iolength = noct ) pm
      open(unit=1250+i,file="OUTPUT/files/"//trim(adjustl(numberchaine))//".bin",STATUS="old",access='direct',RECL=noct,iostat=ok)
      ! -------                                                --------    .
      if (ok .ne. 0) then
        write(*,*)"problème dans lect_mod_select : le fichier OUTPUT/files/"//trim(adjustl(numberchaine))//".bin n''existe pas "
        stop
      endif
      ! -------                                                --------    .
      k=0
      do while(ok.eq.0)
        k=k+1
        read(1250+i,REC=k,iostat = ok)pm%mis,pm%par
        if (ok .eq. 0) then
          n=n+1
          dp%mis%vec(n)=pm%mis
          dp%VC%vec(n)=pm%par%VC
          dp%VM%vec(n)=pm%par%VM
          dp%Zmoho%vec(n)=pm%par%Zmoho
          dp%VpVs%vec(n)=pm%par%VpVs
          do j=1, nbseismes
            dp%Lat(j)%vec(n)=pm%par%Lat(j)
            dp%Lon(j)%vec(n)=pm%par%Lon(j)
            dp%Zhypo(j)%vec(n)=pm%par%Zhypo(j)
            temps_mod=pm%par%Tzero(j)
            call difftime(dp%Tzero(j)%vec(n),temps_mod,dp%temps_ref(j))    ! dp%Tzero : temps relatif à temps_ref
          enddo
        endif
      enddo
      ! ---------------------------------------------------------------    .
      ! calcul des fonction d'autovariance pour chaque parametre et chaque chaîne
      ! ---------------------------------------------------------------    .
      m=min(autocorr,k/3)
      ! -------                                                --------    .
      namefile="OUTPUT/GMT/autovar_VC"//trim(adjustl(numberchaine))//".txt"
      call autovariance(dp%VC%vec(l:n),k,m,namefile)
      namefile="OUTPUT/GMT/autovar_VM"//trim(adjustl(numberchaine))//".txt"
      call autovariance(dp%VM%vec(l:n),k,m,namefile)
      namefile="OUTPUT/GMT/autovar_Zmoho"//trim(adjustl(numberchaine))//".txt"
      call autovariance(dp%Zmoho%vec(l:n),k,m,namefile)
      namefile="OUTPUT/GMT/autovar_VpVs"//trim(adjustl(numberchaine))//".txt"
      call autovariance(dp%VpVs%vec(l:n),k,m,namefile)
      do j=1, nbseismes
        write(numberseisme(1:5),'(i5)')j
        namefile="OUTPUT/GMT/autovar_lon_"//trim(adjustl(numberseisme))//"_"//trim(adjustl(numberchaine))//".txt"
        call autovariance(dp%Lon(j)%vec(l:n),k,m,namefile)
        namefile="OUTPUT/GMT/autovar_lat_"//trim(adjustl(numberseisme))//"_"//trim(adjustl(numberchaine))//".txt"
        call autovariance(dp%Lat(j)%vec(l:n),k,m,namefile)
        namefile="OUTPUT/GMT/autovar_Zhypo_"//trim(adjustl(numberseisme))//"_"//trim(adjustl(numberchaine))//".txt"
        call autovariance(dp%Zhypo(j)%vec(l:n),k,m,namefile)
        namefile="OUTPUT/GMT/autovar_Tzero_"//trim(adjustl(numberseisme))//"_"//trim(adjustl(numberchaine))//".txt"
        call autovariance(dp%Tzero(j)%vec(l:n),k,m,namefile)
      enddo
      ! ---------------------------------------------------------------    .
      l=n+1
      ! -------                                                --------    .
      close(1250+i)
    enddo
    ! -------                                                  --------    .
    if (n.ne.dp%nbparam) then
      write(*,*)"problème dans lect_mod_select : mauvaise lecture des modeles sélectionnés"
    endif
    ! ------- nom des légendes GMT                             --------    . format GMT
    dp%mis%char     ='       fonction co\373t       '
    dp%VC%char      =' V@-C @-\050km\056s@+-1@+\051 '
    dp%VM%char      =' V@-M @-\050km\056s@+-1@+\051 '
    dp%Zmoho%char   ='       moho \050km\051        '
    dp%VpVs%char    ='       V@-P @- / V@-S @-      '
    do i=1,nbseismes
      dp%Lat(i)%char     ='     latitude \050\260\051    '
      dp%Lon(i)%char     ='    longitude \050\260\051    '
      dp%Zhypo(i)%char   ='    Hypocentre \050km\051     '
      dp%Tzero(i)%char   ='  temps initial \050s\051  '
    enddo
    ! ------- nom des fichiers d'entrée GMT                    --------    . format GMT
    dp%mis%name="mis"
    dp%VC%name="_vc"
    dp%VM%name="_vm"
    dp%Zmoho%name="_zm"
    dp%VpVs%name="vps"
    do i=1,nbseismes
      dp%Lat(i)%name="lat"
      dp%Lon(i)%name="lon"
      dp%Zhypo(i)%name="_zh"
      dp%Tzero(i)%name="_to"
    enddo
    ! -------                                                  --------    .
    ! cherche les bornes max et min du prior pour les paramètres de structures + pour la profondeur du séisme (Zhypo)
    ! -------                                                  --------    .
    dp%VC%themin=p%mini%VC
    dp%VM%themin=p%mini%VM
    dp%Zmoho%themin=p%mini%Zmoho
    dp%VpVs%themin=p%mini%VpVs
    do i=1,nbseismes
      dp%Zhypo(i)%themin=p%mini%Zhypo(i)
    enddo
    dp%VC%themax=p%maxi%VC
    dp%VM%themax=p%maxi%VM
    dp%Zmoho%themax=p%maxi%Zmoho
    dp%VpVs%themax=p%maxi%VpVs
    do i=1,nbseismes
      dp%Zhypo(i)%themax=p%maxi%Zhypo(i)
    enddo
    ! -------                                                  --------    .
    ! cherche les bornes max et min des modéles sélectionnés pour les paramètres de hypocenraux (lon,lat,temps initial) et la fonction coût
    ! -------                                                  --------    .
    minmax_1=100000.0_wr
    minmax_2=100000.0_wr
    minmax_3=100000.0_wr
    minmax_4=100000.0_wr
    minmax_5=-100000.0_wr
    minmax_6=-100000.0_wr
    minmax_7=-100000.0_wr
    minmax_8=-100000.0_wr
    do i=1,dp%nbparam
      if(dp%mis%vec(i).lt.minmax_1) then
        minmax_1=dp%mis%vec(i)
        dp%mis%themin=minmax_1
      endif
      if(dp%mis%vec(i).gt.minmax_5) then
        minmax_5=dp%mis%vec(i)
        dp%mis%themax=minmax_5
      endif
      do j=1,nbseismes
        if(dp%Lon(j)%vec(i).lt.minmax_2(j)) then
          minmax_2(j)=dp%Lon(j)%vec(i)
          dp%Lon(j)%themin=minmax_2(j)
        endif
        if(dp%Lat(j)%vec(i).lt.minmax_3(j)) then
          minmax_3(j)=dp%Lat(j)%vec(i)
          dp%Lat(j)%themin=minmax_3(j)
        endif
        if(dp%Tzero(j)%vec(i).lt.minmax_4(j)) then
          minmax_4(j)=dp%Tzero(j)%vec(i)
          dp%Tzero(j)%themin=minmax_4(j)
        endif
        if(dp%Lon(j)%vec(i).gt.minmax_6(j)) then
          minmax_6(j)=dp%Lon(j)%vec(i)
          dp%Lon(j)%themax=minmax_6(j)
        endif
        if(dp%Lat(j)%vec(i).gt.minmax_7(j)) then
          minmax_7(j)=dp%Lat(j)%vec(i)
          dp%Lat(j)%themax=minmax_7(j)
        endif
        if(dp%Tzero(j)%vec(i).gt.minmax_8(j)) then
          minmax_8(j)=dp%Tzero(j)%vec(i)
          dp%Tzero(j)%themax=minmax_8(j)
        endif
      enddo
    enddo
    ! -------                                                  --------    .
    gap=0.05_wr                                                            ! ajoute un gap de 5% autour du min et du max pour affichage
    do j=1,nbseismes
      size=dp%Lon(j)%themax-dp%Lon(j)%themin
      dp%Lon(j)%themax=dp%Lon(j)%themax+gap*size
      dp%Lon(j)%themin=dp%Lon(j)%themin-gap*size
      size=dp%Lat(j)%themax-dp%Lat(j)%themin
      dp%Lat(j)%themax=dp%Lat(j)%themax+gap*size
      dp%Lat(j)%themin=dp%Lat(j)%themin-gap*size
    enddo
    gap=0.01_wr                                                            ! ajoute un gap de 1% autour du min et du max pour affichage
    do j=1,nbseismes
      size=dp%Tzero(j)%themax-dp%Tzero(j)%themin
      dp%Tzero(j)%themax=dp%Tzero(j)%themax+gap*size
      dp%Tzero(j)%themin=dp%Tzero(j)%themin-gap*size
    enddo
    size=dp%mis%themax-dp%mis%themin
    dp%mis%themax=dp%mis%themax+gap*size
    dp%mis%themin=dp%mis%themin-gap*size
    ! -----------------------------------------------------------------    .
  end subroutine lect_mod_select

    ! -----------------------------------------------------------------    .

  subroutine moy_mod_select(dp,nbChaineMV,param_best,misfit)
    ! -------                                                  --------    .mh
    ! calcul moy(s), ecart-type(s) et mode des modèles sélectionnés, sur tous les paramètres
    ! -------                                                  --------    .
    use typetemps, only : densityplot, parametres, fcout
    use time, only : difftime
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    type(densityplot), intent(inout) :: dp
    type(parametres), intent(in) :: param_best(nbChaineMV)
    type(fcout), intent(in) :: misfit(nbChaineMV)
    integer(KIND=wi), intent(in) :: nbChaineMV
    ! -------                                                  --------    .
    real(KIND=wr) :: vec(nbChaineMV)
    integer :: i,j
    ! -----------------------------------------------------------------    .
    if(dp%nbparam.gt.11000) then                                            ! si au moins 11000 modèles (on travail sur des grands nombres, tout de même)
    ! calcul des moyennes :
      call moy_mod_select_one (dp%mis,dp%mis,dp%nbparam,dp%deltaxy,dp%mis%vec10000)
      call moy_mod_select_one (dp%VC,dp%mis,dp%nbparam,dp%deltaxy,dp%VC%vec10000)
      call moy_mod_select_one (dp%VM,dp%mis,dp%nbparam,dp%deltaxy,dp%VM%vec10000)
      call moy_mod_select_one (dp%Zmoho,dp%mis,dp%nbparam,dp%deltaxy,dp%Zmoho%vec10000)
      call moy_mod_select_one (dp%VpVs,dp%mis,dp%nbparam,dp%deltaxy,dp%VpVs%vec10000)
      do i=1,nbseismes
        call moy_mod_select_one (dp%Lat(i),dp%mis,dp%nbparam,dp%deltaxy,dp%Lat(i)%vec10000)
        call moy_mod_select_one (dp%Lon(i),dp%mis,dp%nbparam,dp%deltaxy,dp%Lon(i)%vec10000)
        call moy_mod_select_one (dp%Zhypo(i),dp%mis,dp%nbparam,dp%deltaxy,dp%Zhypo(i)%vec10000)
        call moy_mod_select_one (dp%Tzero(i),dp%mis,dp%nbparam,dp%deltaxy,dp%Tzero(i)%vec10000)
      enddo
    else
      write(*,*)'problème dans moy_mod_select : pas assez de modèles, d''un point de vue statistique'
      stop
    endif

    ! -----------------------------------------------------------------    . 
    ! ecriture des moyenne et écartyes dans paramHypo.d et paramTerre.d
    ! -> fichiers ustilisables comme INPUT dans une autre execution ...
    ! -----------------------------------------------------------------    .
    open(unit=50,file="OUTPUT/input/paramTerre_new.d",STATUS="replace")
      write(50,*)dp%VC%moy_1000,dp%VC%ec_1000
      write(50,*)dp%VM%moy_1000,dp%VM%ec_1000
      write(50,*)dp%Zmoho%moy_1000,dp%Zmoho%ec_1000
      write(50,*)dp%VpVs%moy_1000,dp%VpVs%ec_1000
    close(50)
    open(unit=51,file="OUTPUT/input/paramHypo_new.d",STATUS="replace")
      do i=1,nbseismes
        write(51,*)dp%lon(i)%moy_1000,dp%lon(i)%ec_1000
        write(51,*)dp%lat(i)%moy_1000,dp%lat(i)%ec_1000
        write(51,*)dp%Zhypo(i)%moy_1000,dp%Zhypo(i)%ec_1000
        write(51,*)dp%temps_ref(i)
        write(51,*)dp%Tzero(i)%moy_1000,dp%Tzero(i)%ec_1000
      enddo
    close(51)
    ! -----------------------------------------------------------------    .
    ! -------                                                  --------    .
    ! calcul des moyennes sur l'ensemble du meilleur modèle de chaque chaîne
    do i=1,nbChaineMV
      vec(i)=misfit(i)%best
    enddo
    call moy_ec(vec,nbChaineMV,nbChaineMV,dp%mis%moy_bestchaine,dp%mis%ec_bestchaine)
    ! -------                                                  --------    .
    do i=1,nbChaineMV
      vec(i)=param_best(i)%VC
    enddo
    call moy_ec(vec,nbChaineMV,nbChaineMV,dp%VC%moy_bestchaine,dp%VC%ec_bestchaine)
    ! -------                                                  --------    .
    do i=1,nbChaineMV
      vec(i)=param_best(i)%VM
    enddo
    call moy_ec(vec,nbChaineMV,nbChaineMV,dp%VM%moy_bestchaine,dp%VM%ec_bestchaine)
    ! -------                                                  --------    .
    do i=1,nbChaineMV
      vec(i)=param_best(i)%Zmoho
    enddo
    call moy_ec(vec,nbChaineMV,nbChaineMV,dp%Zmoho%moy_bestchaine,dp%Zmoho%ec_bestchaine)
    ! -------                                                  --------    .
    do i=1,nbChaineMV
      vec(i)=param_best(i)%VpVs
    enddo
    call moy_ec(vec,nbChaineMV,nbChaineMV,dp%VpVs%moy_bestchaine,dp%VpVs%ec_bestchaine)
    ! -------                                                  --------    .
    do j=1,nbseismes
      do i=1,nbChaineMV
        vec(i)=param_best(i)%Lat(j)
      enddo
      call moy_ec(vec,nbChaineMV,nbChaineMV,dp%Lat(j)%moy_bestchaine,dp%Lat(j)%ec_bestchaine)
      ! -------                                                --------    .
      do i=1,nbChaineMV
        vec(i)=param_best(i)%Lon(j)
      enddo
      call moy_ec(vec,nbChaineMV,nbChaineMV,dp%Lon(j)%moy_bestchaine,dp%Lon(j)%ec_bestchaine)
      ! -------                                                --------    .
      do i=1,nbChaineMV
        vec(i)=param_best(i)%Zhypo(j)
      enddo
      call moy_ec(vec,nbChaineMV,nbChaineMV,dp%Zhypo(j)%moy_bestchaine,dp%Zhypo(j)%ec_bestchaine)
      ! -------                                                --------    .
      do i=1,nbChaineMV
        call difftime(vec(i),param_best(i)%Tzero(j),dp%temps_ref(j))
      enddo
      call moy_ec(vec,nbChaineMV,nbChaineMV,dp%Tzero(j)%moy_bestchaine,dp%Tzero(j)%ec_bestchaine)
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine moy_mod_select

    ! -----------------------------------------------------------------    .

  subroutine moy_mod_select_one (a_param,mis_param,nbparam,deltaxy,VEC_10000_tri)
    ! -------                                                  --------    .mh
    ! calcul des moyennes, ecart-types et modes des modèles sélectionnés, pour un unique paramètre
    ! là, j'aurais pu être plus élégant au niveau syntaxique,
    ! en utilisant les fonctions de la norme f90 
    ! au lieu de boucler, comme un cochon, sur tout.
    ! mais pour ~500 000 modèles, ca reste rapide ...
    ! -------                                                  --------    .
    use typetemps , only : densityplot_one
    use cpt_temps
    use tri , only : tri_bulle
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    type(densityplot_one), intent(inout) :: a_param
    type(densityplot_one), intent(in) :: mis_param
    integer(KIND=wi), intent(in) :: nbparam
    integer(KIND=wi), intent(in) :: deltaxy
    real(KIND=wr), intent(out)  :: VEC_10000_tri(10000,2)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,n,size
    integer(KIND=wi) :: vmode(deltaxy)
    real(KIND=wr) :: VEC_10000(10000,2), max, med
    integer(KIND=wi), save :: count = 0
    character (LEN=30) :: chaine
    logical :: ok
    ! -----------------------------------------------------------------    .
    a_param%delta=(a_param%themax-a_param%themin)/real(deltaxy,wr)         ! calcul le pas (discretisation pour le mode et diagramme de densité)
    ! -----------------------------------------------------------------    .
    ! vecteur avec les 10 000 meilleurs modèles (les modèles sont tous différents ; c.-à.-d. sans doublons)
    ! -----------------------------------------------------------------    .
    do i=1,10000                                                           ! initialisation
      VEC_10000(i,1)=-1.0_wr
      VEC_10000(i,2)=10000.0_wr+real(i*100,wr)
    enddo
    ! -------                                                  --------    .
    bestmod : do i=1,nbparam
      write(chaine(1:30),*)"paramètre : ",a_param%name                     ! chaîne pour la barre d'avancement
      count=count+1
      call progress(count,nbparam*(5+4*nbseismes),chaine)                  ! barre d'avancement : il existe (4 + 4*nbseismes) paramètres + (1) la fonction coût
      ! -------                                                --------    . cherche les 10 000 meilleurs modèles (non uniques)
      max=-1.0_wr
      piremod : do j=1,10000                                               ! cherche le pire modèles (possédant le plus grand misfit)
        ok=.true.
        if ((VEC_10000(j,2)==mis_param%vec(i)).and.(VEC_10000(j,1)==a_param%vec(i))) then
          ok=.false.                                                       ! modèle existe déjà, pas de doublons
          exit piremod
        else
          if(VEC_10000(j,2).gt.max) then
            max=VEC_10000(j,2)
            n=j
          endif
        endif
      enddo piremod
      if((ok).and.(VEC_10000(n,2).gt.mis_param%vec(i))) then               ! remplace le pire modèle
        VEC_10000(n,1)=a_param%vec(i)
        VEC_10000(n,2)=mis_param%vec(i)
      endif
    enddo bestmod
    ! -----------------------------------------------------------------    .
    n=10000
    call tri_bulle(VEC_10000,n,VEC_10000_tri)                              ! tri croissant des meilleurs modèles (sans doublons)
    if (VEC_10000_tri(10000,1)==-1.0_wr) then
      write(*,*)'problème dans moy_mod_select_one : pas assez de modèles différents (< 10 000)'
      stop
    endif
    ! -----------------------------------------------------------------    .
    a_param%best=VEC_10000_tri(1,1)                                        ! modèle possédant la fonction coût la plus basse
    ! -----------------------------------------------------------------    .
    do i=1,deltaxy                                                         ! initialise, calcul du mode
      vmode(i)=0
    enddo
    ! -------                                                  --------    . calcul de la médiane
    med=2.0_wr
    call mediane(med,a_param%vec,nbparam,a_param%mediane)
    ! -------                                                  --------    . calcul du mode
    do i=1,nbparam
     if ((a_param%vec(i).lt.a_param%themin).or.(a_param%vec(i).gt.a_param%themax)) then
        write(*,*)'problème dans moy_mod_select_one 1 : calcul du mode impossible',i,a_param%name, &
        a_param%vec(i),a_param%themin,a_param%themax,a_param%delta,j
        ! stop
      else
        j = int((a_param%vec(i)-a_param%themin)/a_param%delta)+1
        if ((j.gt.0).and.(j.le.deltaxy)) then
          vmode(j) = vmode(j) + 1
        else
          write(*,*)'problème dans moy_mod_select_one 2 : calcul du mode impossible',i,a_param%name, &
          a_param%vec(i),a_param%themin,a_param%themax,a_param%delta,j
          stop
        endif

      endif
    enddo
    ! -------                                                  --------    .
    j=0
    max=-1.0_wr
    do i=1,deltaxy                                                         ! mode : maximum
      if(real(vmode(i),wr).gt.max) then
        max = real(vmode(i),wr)
        j=i
      endif
    enddo
    a_param%mode = a_param%themin + (real(j,wr)-0.5_wr) * a_param%delta
    ! -------                                                  --------    . calcul des moyennes et des écart-types
    call moy_ec(a_param%vec,nbparam,nbparam,a_param%moy_tot,a_param%ec_tot)! pour les nbparam modèles (tous les sélectionnés)
    n=10000
    size=n
    call moy_ec(VEC_10000_tri(:,1),size,n,a_param%moy_10000,a_param%ec_10000) ! pour les 10000 meilleurs modèles
    n=1000
    call moy_ec(VEC_10000_tri(:,1),size,n,a_param%moy_1000,a_param%ec_1000) ! pour les 1000 meilleurs modèles
    n=100
    call moy_ec(VEC_10000_tri(:,1),size,n,a_param%moy_100,a_param%ec_100) ! pour les 100 meilleurs modèles
    ! -----------------------------------------------------------------    .
  end subroutine moy_mod_select_one

    ! -----------------------------------------------------------------    .

  subroutine mediane(val,tab,bufLength,med)
    ! -------                                                --------    .mh
    ! calcul de la médiane (si val=2.0)
    ! calcul du premier quatrile (si val=4.0)
    ! calcul du second quatrile (si val=4.0/3.0)
    ! -------                                                --------    .
    ! Modified algorithm according to http://www.geocities.com/zabrodskyvlada/3alg.html
    ! Contributed by Heinz Klar
    ! -------                                                --------    .
    implicit none
    integer (KIND=wi), intent(in) :: bufLength
    real (KIND=wr), intent(in) :: val
    real (KIND=wr), intent(in) :: tab(bufLength)
    real (KIND=wr), intent(out) :: med
    ! -------                                                --------    .
    real (KIND=wr) :: buf(bufLength)
    integer (KIND=wi):: i,j,n,l,m
    real (KIND=wr) :: dum
    ! ---------------------------------------------------------------    .
    if(bufLength.lt.750000) then
      buf=tab
      m=bufLength-1
      n=int(real(bufLength,8)/val)
      med=buf(n)
      l=1
      ! -------                                              --------    .
      do while (l.lt.m)
        i=l
        j=m
        do while ((j.ge.n).and.(i.le.n))
          do while (buf(i).lt.med)
            i=i+1
          enddo
          do while (med.lt.buf(j))
            j=j-1
          enddo
          dum=buf(j)
          buf(j)=buf(i)
          buf(i)=dum
          i=i+1
          j=j-1
        enddo
        if (j.lt.n) l=i
        if (n.lt.i) m=j
        med=buf(n)
      enddo
    else
      write(*,*)'problème dans mediane : trop de modèles ',bufLength
      med=0.0_wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine mediane

    ! -----------------------------------------------------------------    .

  subroutine inR(D,R,nbtps,nbstaR,mis)
    ! -------                                                  --------    .mh
    ! lecture des résidus dans D, à chaque itération  et stock dans R
    ! dans la suite R est écrit dans un fichier (cf outR) | (si FLAGresSTA=.true.)
    ! -------                                                  --------    .
    use typetemps, only : dataall, residus
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent(in) :: nbtps(nbseismes), nbstaR
    real(KIND=wr), intent(in) :: mis
    type(dataall), intent(in) :: D(nbseismes)
    type(residus), intent(inout) :: R(nbstaR)                              ! résidus
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k
    ! -----------------------------------------------------------------    .
    do i=1,nbseismes
      do j=1,nbtps(i)
        do k=1,nbstaR
          if (D(i)%datatps(j)%sta%staname==R(k)%staname) then
            if (D(i)%datatps(j)%typeonde=='G') then
              R(k)%nbPg = R(k)%nbPg + 1
              R(k)%resPg(R(k)%nbPg,1) = D(i)%datatps(j)%dTP
              R(k)%resPg(R(k)%nbPg,2) = mis
              if (D(i)%datatps(j)%andS=='S') then
                R(k)%nbSg = R(k)%nbSg + 1
                R(k)%resSg(R(k)%nbSg,1) = D(i)%datatps(j)%dTS
                R(k)%resSg(R(k)%nbSg,2) = mis
              endif
            elseif (D(i)%datatps(j)%typeonde=='N') then
              R(k)%nbPn = R(k)%nbPn + 1
              R(k)%resPn(R(k)%nbPn,1) = D(i)%datatps(j)%dTP
              R(k)%resPn(R(k)%nbPn,2) = mis
              if (D(i)%datatps(j)%andS=='S') then
                R(k)%nbSn = R(k)%nbSn + 1
                R(k)%resSn(R(k)%nbSn,1) = D(i)%datatps(j)%dTS
                R(k)%resSn(R(k)%nbSn,2) = mis
              endif
            else
              write(*,*)'problème dans inR : onde ni directe ni réfractée ... ? '
              stop
            endif
          endif
        enddo
      enddo
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine inR

   ! -----------------------------------------------------------------    .

  subroutine outR(R,nbstaR)
    ! -------                                                  --------    .mh
    ! écriture des résidus dans des fichiers | (si FLAGresSTA=.true.)
    ! -------                                                  --------    .
    use typetemps, only : residus
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent(in) :: nbstaR
    type(residus), intent(inout) :: R(nbstaR)                              ! résidus (si FLAGresSTA=.true.)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,noctet,ok
    real(KIND=wr) :: min,max,val
    ! -----------------------------------------------------------------    .
    ok = 0
    i = 0
    do while (ok == 0)
      i=i+1
      if (R(i)%nbPg.gt.1) then                                             ! il existe des ondes directes
        inquire (iolength = noctet)R(i)%resPg(1,1)
        ok=1
      elseif (R(i)%nbPn.gt.1) then
        inquire (iolength = noctet)R(i)%resPn(1,1)                         ! ou sinon, au moins des réfractées
        ok=1
      endif
      if (i.gt.nbstaR) then
        write(*,*)'problème dans outR', i
        stop
      endif
    enddo
    ! -------                                                  --------    .
    do i=1,nbstaR
      ! ---------------------------------------------------------------    . Pg
      if (R(i)%nbPg.gt.0) then
        open(unit=1000,file="OUTPUT/files/STA/"//R(i)%staname//"-PG.bin",access='direct',RECL=(noctet))
        min=1.e9_wr
        max=-1.e9_wr
        do j=1,R(i)%nbPg
          if(min.gt.R(i)%resPg(j,2)) min=R(i)%resPg(j,2)
          if(max.lt.R(i)%resPg(j,2)) max=R(i)%resPg(j,2)
        enddo
        ! -------                                              --------    . 
        ! on garde dans les 10 meilleurs %
        ! c'est pas très propre, mais ca a le mérite d'être rapide !
        val=min+(max-min)/10.0_wr
        ! -------                                              --------    .
        do j=1,R(i)%nbPg
          if(R(i)%resPg(j,1).le.val) write(1000,rec=R(i)%nbPgT+j)R(i)%resPg(j,1)
        enddo
        close(1000)
        R(i)%nbPgT=R(i)%nbPgT+R(i)%nbPg
        R(i)%nbPg=0
      endif
      ! ---------------------------------------------------------------    . Pn
      if (R(i)%nbPn.gt.0) then
        open(unit=1000,file="OUTPUT/files/STA/"//R(i)%staname//"-PN.bin",access='direct',RECL=(noctet))
        min=1.e9_wr
        max=-1.e9_wr
        do j=1,R(i)%nbPn
          if(min.gt.R(i)%resPn(j,2)) min=R(i)%resPn(j,2)
          if(max.lt.R(i)%resPn(j,2)) max=R(i)%resPn(j,2)
        enddo
        ! -------                                              --------    . 
        ! on garde dans les 10 meilleurs %
        ! c'est pas très propre, mais ca a le mérite d'être rapide !
        val=min+(max-min)/10.0_wr
        ! -------                                              --------    .
        do j=1,R(i)%nbPn
          if(R(i)%resPn(j,1).le.val) write(1000,rec=R(i)%nbPnT+j)R(i)%resPn(j,1)
        enddo
        close(1000)
        R(i)%nbPnT=R(i)%nbPnT+R(i)%nbPn
        R(i)%nbPn=0
      endif
      ! ---------------------------------------------------------------    . Sg
      if (R(i)%nbSg.gt.0) then
        open(unit=1000,file="OUTPUT/files/STA/"//R(i)%staname//"-SG.bin",access='direct',RECL=(noctet))
        min=1.e9_wr
        max=-1.e9_wr
        do j=1,R(i)%nbSg
          if(min.gt.R(i)%resSg(j,2)) min=R(i)%resSg(j,2)
          if(max.lt.R(i)%resSg(j,2)) max=R(i)%resSg(j,2)
        enddo
        ! -------                                              --------    . 
        ! on garde dans les 10 meilleurs %
        ! c'est pas très propre, mais ca a le mérite d'être rapide !
        val=min+(max-min)/10.0_wr
        ! -------                                              --------    .
        do j=1,R(i)%nbSg
          if(R(i)%resSg(j,1).le.val) write(1000,rec=R(i)%nbSgT+j)R(i)%resSg(j,1)
        enddo
        close(1000)
        R(i)%nbSgT=R(i)%nbSgT+R(i)%nbSg
        R(i)%nbSg=0
      endif
      ! ---------------------------------------------------------------    . Sn
      if (R(i)%nbSn.gt.0) then
        open(unit=1000,file="OUTPUT/files/STA/"//R(i)%staname//"-SN.bin",access='direct',RECL=(noctet))
        min=1.e9_wr
        max=-1.e9_wr
        do j=1,R(i)%nbSn
          if(min.gt.R(i)%resSn(j,2)) min=R(i)%resSn(j,2)
          if(max.lt.R(i)%resSn(j,2)) max=R(i)%resSn(j,2)
        enddo
        ! -------                                              --------    . 
        ! on garde dans les 10 meilleurs %
        ! c'est pas très propre, mais ca a le mérite d'être rapide !
        val=min+(max-min)/10.0_wr
        ! -------                                              --------    .
        do j=1,R(i)%nbSn
          if(R(i)%resSn(j,1).le.val) write(1000,rec=R(i)%nbSnT+j)R(i)%resSn(j,1)
        enddo
        close(1000)
        R(i)%nbSnT=R(i)%nbSnT+R(i)%nbSn
        R(i)%nbSn=0
      endif
      ! ---------------------------------------------------------------    .
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine outR

    ! -----------------------------------------------------------------    .

  subroutine moy_ec(param,size,n,moy,ec)
    ! -------                                                  --------    .mh
    ! calcul moy et ecart-type des n meilleurs modèles
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent(in) :: size, n
    real(KIND=wr), intent(in) :: param(size)
    real(KIND=wr), intent(out) :: moy,ec
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    ! -----------------------------------------------------------------    .
    if (n.le.size) then
      moy = 0.0_wr
      ec = 0.0_wr
      do i = 1,n
        moy =  moy + param(i)
      enddo
      moy =  moy / real(n,wr)
      do i = 1,n
        ec = ec + (param(i)-moy)**2.0_wr
      enddo
      ec = sqrt(ec/real(n,wr))
    else
      write(*,*)'problème dans moy_ec : calcul de moyenne impossible'
      stop
    endif
    ! -----------------------------------------------------------------    .
  end subroutine moy_ec

    ! -----------------------------------------------------------------    .

  subroutine dist_apriori(j,rang,nbtps,D,nbm_ap,temps_ref,pEpis,nb,acentroid)
    ! -------                                                  --------    .mh
    !  calcul une densité de probabilité a priori pour chaque paramètre
    !  c'est à dire un point de départ de la Chaine de Markov
    !  de plus, la méthode Geiger est appliqué sur chaque séismes avec des modèles aléatoirs de terre
    ! -----------------------------------------------------------------    .
    use typetemps, only : dataall, date_sec, parametre, parametres, priorEPI, amoho_centroid, mvPall_2_P1
    use time, only : basetime,difftime
    use tri, only : triparam
    use mt19937
    use cpt_temps
    use invGEIGER
    use rechercheepi
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: j
    integer(KIND=wi), intent (in) :: rang
    integer(KIND=wi), intent(in) :: nbtps(nbseismes)
    type(dataall), intent(inout) :: D(nbseismes)
    integer(KIND=wi), intent(in) :: nbm_ap
    type(date_sec), intent(in) :: temps_ref(nbseismes)
    type(priorEPI), intent (inout) :: pEpis(nbseismes)
    integer(KIND=wi), intent (in) :: nb                                    ! au carré, nb de cases possible
    type (amoho_centroid), intent (in) :: acentroid
    ! -------                                                  --------    .
    type(parametres), allocatable :: papriori(:)
    type(parametre), allocatable :: one_param(:)
    real(KIND=wr), allocatable :: valmis(:)
    type(parametre) :: pgeiger,pgeigertest,pgeigermoy,pgeigerec
    integer(KIND=wi) :: i,n,l,k
    real(KIND=wr) :: val,moy,ec,moyMIS,ecMIS
    character(len=5) :: x
    character(LEN=30) :: one_string,two_string
    real(KIND=wr) :: onemisfit
    logical :: conv                                                        ! "true" si geiger convergent
    ! -----------------------------------------------------------------    .
    character (LEN=1) :: mod
    ! -----------------------------------------------------------------    .
    allocate(papriori(nbm_ap))
    ! -----------------------------------------------------------------    .   
    ! on ne sauve les modèles de terre que pour un rang :
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(one_string(1:30),'(a30)')"                                                      "
    if (rang==0) then
      open(unit=992,file="OUTPUT/GMT/aprio_vc-1.bin",STATUS="replace",access='direct',RECL=8)
      open(unit=993,file="OUTPUT/GMT/aprio_vm-1.bin",STATUS="replace",access='direct',RECL=8)
      open(unit=994,file="OUTPUT/GMT/aprio_zm-1.bin",STATUS="replace",access='direct',RECL=8)
      open(unit=995,file="OUTPUT/GMT/apriovps-1.bin",STATUS="replace",access='direct',RECL=8)
    endif
    ! -------                                                  --------    .
    do i=1,nbm_ap
      if (rang==0) write(one_string(1:30),'(a20,1x,i4,a1,i4)')" densité apriori:   ",j,"/",nbseismes
      if (rang==0) call progress(i,nbm_ap,one_string)
      ! -------                                                --------    .
      call initparam(nbtps,D,papriori(i),pEpis,nb)
      ! -------                                                --------    . apriori
      if (rang==0) then
        write(992,REC=i)real(papriori(i)%VC,8)
        write(993,REC=i)real(papriori(i)%VM,8)
        write(994,REC=i)real(papriori(i)%Zmoho,8)
        write(995,REC=i)real(papriori(i)%VpVs,8)
      endif
    enddo
    ! -------                                                  --------    . fin fichiers pour le modele de terre
    if (rang==0) then
      close(992)
      close(993)
      close(994)
      close(995)
    endif
    ! -----------------------------------------------------------------    . apriori
    write(x(1:5),'(i5)')j
    open(unit=996,file="OUTPUT/GMT/aprio_zh-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    open(unit=997,file="OUTPUT/GMT/apriolon-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    open(unit=998,file="OUTPUT/GMT/apriolat-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    open(unit=999,file="OUTPUT/GMT/aprio_to-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    do i=1,nbm_ap
      write(996,REC=i)real(papriori(i)%Zhypo(j),8)
      write(997,REC=i)real(papriori(i)%Lon(j),8)
      write(998,REC=i)real(papriori(i)%Lat(j),8)
      call difftime(val,papriori(i)%Tzero(j),temps_ref(j))
      write(999,REC=i)real(val,8)
    enddo
    close(996)
    close(997)
    close(998)
    close(999)
    write(two_string(1:30),'(a30)')"                                                      "
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(x(1:5),'(i5)')j
    ! -------                                                  --------    .
    open(unit=1996,file="OUTPUT/GMT/geiger_zh-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    open(unit=1997,file="OUTPUT/GMT/geigerlon-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    open(unit=1998,file="OUTPUT/GMT/geigerlat-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    open(unit=1999,file="OUTPUT/GMT/geiger_to-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=8)
    do i=1,nbm_ap/50
      ! -------                                                --------    .
      if (rang==0) then
        write(two_string(1:30),'(a20,1x,i4,a1,i4)')" méthode Geiger : 1 ",j,"/",nbseismes
        call progress(i,nbm_ap/50,two_string)
      endif
      ! -------                                              --------    . geiger
      mod='I'
      call mvPall_2_P1(pgeiger,papriori(i),j)
      call dogeigerone(j,nbtps(j),D(j)%datatps,pgeiger,acentroid,mod,onemisfit,chut=.true.,con=conv)
      ! -------                                              --------    .
      write(1996,REC=i)real(pgeiger%Zhypo,8)
      write(1997,REC=i)real(pgeiger%Lon,8)
      write(1998,REC=i)real(pgeiger%Lat,8)

      call difftime(val,pgeiger%Tzero,temps_ref(j))
      write(1999,REC=i)real(val,8)
      ! -------                                              --------    .
    enddo
    close(1996)
    close(1997)
    close(1998)
    close(1999)
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    n=5                                                                    ! nb de modèles * toutes pfd ...
    ! -----------------------------------------------------------------    .
    ! avec les modèles d'Arroucau                                          !
    ! -----------------------------------------------------------------    .
    open(unit=2999,file="OUTPUT/GMT/Arroucau_geiger-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=40)
    l=0
    do i=1,n
      ! -------                                                --------    .
      if (rang==0) then
        write(two_string(1:30),'(a20,1x,i4,a1,i4)')" méthode Geiger : 2 ",i,"/",n
        call progress(i,n,two_string)
      endif
      ! -------                                                --------    . geiger
      call mvPall_2_P1(pgeiger,papriori(i),j)
      mod='A'
      pgeiger%Zhypo=genrand_real3()
      do while(pgeiger%Zhypo.lt.(32.0_wr-0.1_wr))
        pgeigertest=pgeiger
        call dogeigerone(j,nbtps(j),D(j)%datatps,pgeigertest,acentroid,mod,onemisfit,chut=.true.,con=conv)
        if (conv) then
          l=l+1
          call difftime(val,pgeigertest%Tzero,temps_ref(j))
          write(2999,REC=l)real(onemisfit,8),real(pgeigertest%Zhypo,8),real(pgeigertest%Lon,8),real(pgeigertest%Lat,8),real(val,8)
        endif
        pgeiger%Zhypo=pgeiger%Zhypo+5.0_wr*genrand_real3()
      enddo

      ! -------                                                --------    .
    enddo
    ! -------                                                  --------    . relecture
    if (l.ge.3) then
      ! -------                                                --------    .
      allocate(valmis(l),one_param(l))
      ! -------                                                --------    .
      do i=1,l
        read(2999,REC=i)onemisfit,one_param(i)%Zhypo,one_param(i)%Lon,one_param(i)%Lat,val
        valmis(i)=onemisfit
        one_param(i)%Tzero=temps_ref(j)
        one_param(i)%Tzero%sec = one_param(i)%Tzero%sec + val
        call basetime(one_param(i)%Tzero)
        one_param(i)%VC=0.0_wr
        one_param(i)%VM=0.0_wr
        one_param(i)%Zmoho=0.0_wr
        one_param(i)%VpVs=0.0_wr
      enddo
      ! -------                                                --------    . tri croissant
      call triparam(l,valmis,one_param)
      ! -------                                                --------    . calcul moy et ecart-type des i meilleurs modèles
      i=max(2,3*l/4)
      call moy_ec(valmis,l,i,moyMIS,ecMIS)
      call moy_ec(one_param(:)%Zhypo,l,i,moy,ec)
      pgeigermoy%Zhypo=moy
      pgeigerec%Zhypo=ec
      call moy_ec(one_param(:)%Lat,l,i,moy,ec)
      pgeigermoy%Lat=moy
      pgeigerec%Lat=ec
      call moy_ec(one_param(:)%Lon,l,i,moy,ec)
      pgeigermoy%Lon=moy
      pgeigerec%Lon=ec
      valmis=0.0_wr
      open(unit=299,file="OUTPUT/GMT/ArrALL-"//trim(adjustl(x))//".txt",STATUS="replace")
      open(unit=298,file="OUTPUT/GMT/ArrALLt-"//trim(adjustl(x))//".txt",STATUS="replace")
      do k=1,i
        call difftime(valmis(k),one_param(k)%Tzero,temps_ref(j))
        write(299,*)real(one_param(k)%Lon,8),real(one_param(k)%Lat,8)
        write(298,*)real(one_param(k)%Lon,8),real(one_param(k)%Lat,8),"10 0 5 LM Arroucau"
      enddo
      close(298)
      close(299)
      call moy_ec(valmis,l,i,moy,ec)
      deallocate(valmis,one_param)
      ! -------                                                --------    .
      open(unit=2998,file="OUTPUT/GMT/Arroucau_f-"//trim(adjustl(x))//".d",STATUS="replace")
      write(2998,*)i,moyMIS,ecMIS
      write(2998,*)pgeigermoy%Zhypo,pgeigermoy%Lon,pgeigermoy%Lat,moy
      write(2998,*)pgeigerec%Zhypo,pgeigerec%Lon,pgeigerec%Lat,ec
      close(2998)
      open(unit=2997,file="OUTPUT/GMT/Arroucau_all-"//trim(adjustl(x))//".tex",STATUS="replace")
      write(2997,*)"\begin{table}[!ht] \scriptsize \centering  \renewcommand{\arraystretch}{1.5}"
      write(2997,*)"\begin{tabular}{|>{\centering}m{2.5cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
      write(2997,*)"|>{\centering}m{2.2cm}|>{\centering}m{2.2cm}|m{2.1cm}<{\centering}|}"
      write(2997,*)"\hline "
      write(2997,*)"{\bf \large mod\`ele} & {fonction co\^ut} & {Z$_{hypo}$ (km)} & {longitude (\degree)} & ", &
        " {latitude (\degree)} & {T$_{z\acute{e}ro}$ (s)} \\"
      write(2997,*)"\hline "
      write(2997,*)"\hline "
      write(2997,5007)"mod\`eles de terre de Arroucau & \np{",real(moyMIS,wr), &
        "} $\pm$ \np{",2.0_wr*real(ecMIS,wr),"} & \np{",pgeigermoy%Zhypo,"} $\pm$ \np{",2.0_wr*pgeigerec%Zhypo, &
        "} & \np{",pgeigermoy%Lon,"} $\pm$ \np{",2.0_wr*pgeigerec%Lon,"} & \np{",pgeigermoy%Lat, &
        "} $\pm$ \np{",2.0_wr*pgeigerec%Lat,"} & \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} \\ "
      write(2997,*)"\hline "
      write(2997,*)"\end{tabular}"
      write(2997,*)"\end{table}"
      close(2997)
    else
      open(unit=2998,file="OUTPUT/GMT/Arroucau_f-"//trim(adjustl(x))//".d",STATUS="replace")
      write(2998,*)"0 0.0 0.0"
      write(2998,*)"0.0 0.0 0.0 0.0"
      write(2998,*)"0.0 0.0 0.0 0.0"
      close(2998)
      open(unit=2997,file="OUTPUT/GMT/Arroucau_all-"//trim(adjustl(x))//".tex",STATUS="replace")
      write(2997,*)" "
      close(2997)
      open(unit=299,file="OUTPUT/GMT/ArrALL-"//trim(adjustl(x))//".txt",STATUS="replace")
      write(299,*)"0 0.0 0.0"
      close(299)
      open(unit=298,file="OUTPUT/GMT/ArrALLt-"//trim(adjustl(x))//".txt",STATUS="replace")
      write(298,*)"0 0.0 0.0 10 0 5 LM Arroucau"
      close(298)
    endif
    close(2999)

    ! -----------------------------------------------------------------    .
    ! avec le modèle Si-HEx                                                !
    ! -----------------------------------------------------------------    .

    open(unit=3999,file="OUTPUT/GMT/SiHex_geiger-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=40)
    l=0
    do i=1,n
      ! -------                                                --------    .
      if (rang==0) then
        write(two_string(1:30),'(a20,1x,i4,a1,i4)')" méthode Geiger : 3 ",i,"/",n
        call progress(i,n,two_string)
      endif
      ! -------                                                --------    . geiger
      call mvPall_2_P1(pgeiger,papriori(i),j)
      mod='S'
      pgeiger%Zhypo=genrand_real3()
      do while(pgeiger%Zhypo.le.(30.0_wr-0.1_wr))
        pgeigertest=pgeiger
        call dogeigerone(j,nbtps(j),D(j)%datatps,pgeigertest,acentroid,mod,onemisfit,chut=.true.,con=conv)
        if (conv) then
          l=l+1
          call difftime(val,pgeigertest%Tzero,temps_ref(j))
        write(3999,REC=l)real(onemisfit,8),real(pgeigertest%Zhypo,8),real(pgeigertest%Lon,8),real(pgeigertest%Lat,8),real(val,8)
        endif
        pgeiger%Zhypo=pgeiger%Zhypo+5.0_wr*genrand_real3()
      enddo
      ! -------                                                --------    .
    enddo
    ! -------                                                  --------    . relecture
    if (l.ge.3) then
      ! -------                                                --------    .
      allocate(valmis(l),one_param(l))
      ! -------                                                --------    .
      do i=1,l
        read(3999,REC=i)onemisfit,one_param(i)%Zhypo,one_param(i)%Lon,one_param(i)%Lat,val
        valmis(i)=onemisfit
        one_param(i)%Tzero=temps_ref(j)
        one_param(i)%Tzero%sec = one_param(i)%Tzero%sec + val
        call basetime(one_param(i)%Tzero)
        one_param(i)%VC=0.0_wr
        one_param(i)%VM=0.0_wr
        one_param(i)%Zmoho=0.0_wr
        one_param(i)%VpVs=0.0_wr
      enddo
      ! -------                                                --------    . tri croissant
      call triparam(l,valmis,one_param)
      ! -------                                                --------    . calcul moy et ecart-type des i meilleurs modèles
      i=max(2,3*l/4)
      call moy_ec(valmis,l,i,moyMIS,ecMIS)
      call moy_ec(one_param(:)%Zhypo,l,i,moy,ec)
      pgeigermoy%Zhypo=moy
      pgeigerec%Zhypo=ec
      call moy_ec(one_param(:)%Lat,l,i,moy,ec)
      pgeigermoy%Lat=moy
      pgeigerec%Lat=ec
      call moy_ec(one_param(:)%Lon,l,i,moy,ec)
      pgeigermoy%Lon=moy
      pgeigerec%Lon=ec
      valmis=0.0_wr
      open(unit=399,file="OUTPUT/GMT/SiHexALL-"//trim(adjustl(x))//".txt",STATUS="replace")
      open(unit=398,file="OUTPUT/GMT/SiHexALLt-"//trim(adjustl(x))//".txt",STATUS="replace")
      do k=1,i
        call difftime(valmis(k),one_param(k)%Tzero,temps_ref(j))
        write(399,*)real(one_param(k)%Lon,8),real(one_param(k)%Lat,8)
        write(398,*)real(one_param(k)%Lon,8),real(one_param(k)%Lat,8),"10 0 5 LM SiHex"
      enddo
      close(399)
      close(398)
      call moy_ec(valmis,l,i,moy,ec)
      deallocate(valmis,one_param)
      ! -------                                                --------    . 
      open(unit=3998,file="OUTPUT/GMT/SiHex_f-"//trim(adjustl(x))//".d",STATUS="replace")
      write(3998,*)i,moyMIS,ecMIS
      write(3998,*)pgeigermoy%Zhypo,pgeigermoy%Lon,pgeigermoy%Lat,moy
      write(3998,*)pgeigerec%Zhypo,pgeigerec%Lon,pgeigerec%Lat,ec
      close(3998)
      open(unit=3996,file="OUTPUT/GMT/SiHex_all-"//trim(adjustl(x))//".tex",STATUS="replace")
      write(3996,*)"\begin{table}[!ht] \scriptsize \centering  \renewcommand{\arraystretch}{1.5}"
      write(3996,*)"\begin{tabular}{|>{\centering}m{2.5cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
      write(3996,*)"|>{\centering}m{2.2cm}|>{\centering}m{2.2cm}|m{2.1cm}<{\centering}|}"
      write(3996,*)"\hline "
      write(3996,*)"\hline "
      write(3996,5007)"mod\`eles de terre \textsc{Si--Hex} & \np{",real(moyMIS,wr), &
        "} $\pm$ \np{",2.0_wr*real(ecMIS,wr),"} & \np{",pgeigermoy%Zhypo,"} $\pm$ \np{",2.0_wr*pgeigerec%Zhypo, &
        "} & \np{",pgeigermoy%Lon,"} $\pm$ \np{",2.0_wr*pgeigerec%Lon,"} & \np{",pgeigermoy%Lat, &
        "} $\pm$ \np{",2.0_wr*pgeigerec%Lat,"} & \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} \\ "
      write(3996,*)"\hline "
      write(3996,*)"\end{tabular}"
      write(3996,*)"\end{table}"
      close(3996)
    else
      open(unit=3998,file="OUTPUT/GMT/SiHex_f-"//trim(adjustl(x))//".d",STATUS="replace")
      write(3998,*)"0 0.0 0.0"
      write(3998,*)"0.0 0.0 0.0 0.0"
      write(3998,*)"0.0 0.0 0.0 0.0"
      close(3998)
      open(unit=3996,file="OUTPUT/GMT/SiHex_all-"//trim(adjustl(x))//".tex",STATUS="replace")
      write(3996,*)" "
      close(3996)
      open(unit=399,file="OUTPUT/GMT/SiHexALL-"//trim(adjustl(x))//".txt",STATUS="replace")
      write(399,*)"0 0.0 0.0"
      close(399)
      open(unit=398,file="OUTPUT/GMT/SiHexALLt-"//trim(adjustl(x))//".txt",STATUS="replace")
      write(398,*)"0 0.0 0.0 10 0 5 LM SiHex"
      close(398)
    endif
    close(3999)

    ! -----------------------------------------------------------------    .
    ! avec le modèle du CEA                                                !
    ! -----------------------------------------------------------------    .

    open(unit=4999,file="OUTPUT/GMT/CEA_geiger-"//trim(adjustl(x))//".bin",STATUS="replace",access='direct',RECL=40)
    l=0
    do i=1,n
      ! -------                                                --------    .
      if (rang==0) then
        write(two_string(1:30),'(a20,1x,i4,a1,i4)')" méthode Geiger : 4 ",i,"/",n
        call progress(i,n,two_string)
      endif
      ! -------                                                --------    . geiger
      call mvPall_2_P1(pgeiger,papriori(i),j)
      mod='C'
      pgeiger%Zhypo=genrand_real3()
      do while(pgeiger%Zhypo.le.(25.0_wr-0.1_wr))
        pgeigertest=pgeiger
        call dogeigerone(j,nbtps(j),D(j)%datatps,pgeigertest,acentroid,mod,onemisfit,chut=.true.,con=conv)
        if (conv) then
          l=l+1
          call difftime(val,pgeigertest%Tzero,temps_ref(j))
          write(4999,REC=l)real(onemisfit,8),real(pgeigertest%Zhypo,8),real(pgeigertest%Lon,8),real(pgeigertest%Lat,8),real(val,8)
        endif
        pgeiger%Zhypo=pgeiger%Zhypo+5.0_wr*genrand_real3()
      enddo
      ! -------                                                --------    .
    enddo
    ! -------                                                  --------    . relecture
    if (l.ge.3) then
      ! -------                                                --------    .
      allocate(valmis(l),one_param(l))
      ! -------                                                --------    .
      do i=1,l
        read(4999,REC=i)onemisfit,one_param(i)%Zhypo,one_param(i)%Lon,one_param(i)%Lat,val
        valmis(i)=onemisfit
        one_param(i)%Tzero=temps_ref(j)
        one_param(i)%Tzero%sec = one_param(i)%Tzero%sec + val
        call basetime(one_param(i)%Tzero)
        one_param(i)%VC=0.0_wr
        one_param(i)%VM=0.0_wr
        one_param(i)%Zmoho=0.0_wr
        one_param(i)%VpVs=0.0_wr
      enddo
      ! -------                                                --------    . tri croissant
      call triparam(l,valmis,one_param)
      ! -------                                                --------    . calcul moy et ecart-type des i meilleurs modèles
      i=max(2,3*l/4)
      call moy_ec(valmis,l,i,moyMIS,ecMIS)
      call moy_ec(one_param(:)%Zhypo,l,i,moy,ec)
      pgeigermoy%Zhypo=moy
      pgeigerec%Zhypo=ec
      call moy_ec(one_param(:)%Lat,l,i,moy,ec)
      pgeigermoy%Lat=moy
      pgeigerec%Lat=ec
      call moy_ec(one_param(:)%Lon,l,i,moy,ec)
      pgeigermoy%Lon=moy
      pgeigerec%Lon=ec
      valmis=0.0_wr
      open(unit=499,file="OUTPUT/GMT/CEAALL-"//trim(adjustl(x))//".txt",STATUS="replace")
      open(unit=498,file="OUTPUT/GMT/CEAALLt-"//trim(adjustl(x))//".txt",STATUS="replace")
      do k=1,i
        call difftime(valmis(k),one_param(k)%Tzero,temps_ref(j))
        write(499,*)real(one_param(k)%Lon,8),real(one_param(k)%Lat,8)
        write(498,*)real(one_param(k)%Lon,8),real(one_param(k)%Lat,8),"10 0 5 LM CEA"
      enddo
      close(499)
      close(498)
      call moy_ec(valmis,l,i,moy,ec)
      deallocate(valmis,one_param)
      ! -------                                                --------    .
      open(unit=4998,file="OUTPUT/GMT/CEA_f-"//trim(adjustl(x))//".d",STATUS="replace")
      write(4998,*)i,moyMIS,ecMIS
      write(4998,*)pgeigermoy%Zhypo,pgeigermoy%Lon,pgeigermoy%Lat,moy
      write(4998,*)pgeigerec%Zhypo,pgeigerec%Lon,pgeigerec%Lat,ec
      close(4998)
      open(unit=4996,file="OUTPUT/GMT/CEA_all-"//trim(adjustl(x))//".tex",STATUS="replace")
      write(4996,*)"\begin{table}[!ht] \scriptsize \centering  \renewcommand{\arraystretch}{1.5}"
      write(4996,*)"\begin{tabular}{|>{\centering}m{2.5cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
      write(4996,*)"|>{\centering}m{2.2cm}|>{\centering}m{2.2cm}|m{2.1cm}<{\centering}|}"
      write(4996,*)"\hline "
      write(4996,*)"\hline "
      write(4996,5007)"mod\`eles de terre \textsc{CEA} & \np{",real(moyMIS,wr), &
      "} $\pm$ \np{",2.0_wr*real(ecMIS,wr),"} & \np{",pgeigermoy%Zhypo,"} $\pm$ \np{",2.0_wr*pgeigerec%Zhypo, &
      "} & \np{",pgeigermoy%Lon,"} $\pm$ \np{",2.0_wr*pgeigerec%Lon,"} & \np{",pgeigermoy%Lat, &
      "} $\pm$ \np{",2.0_wr*pgeigerec%Lat,"} & \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} \\ "
      write(4996,*)"\hline "
      write(4996,*)"\end{tabular}"
      write(4996,*)"\end{table}"
      close(4996)
    else
      open(unit=4998,file="OUTPUT/GMT/CEA_f-"//trim(adjustl(x))//".d",STATUS="replace")
      write(4998,*)"0 0.0 0.0"
      write(4998,*)"0.0 0.0 0.0 0.0"
      write(4998,*)"0.0 0.0 0.0 0.0"
      close(4998)
      open(unit=4996,file="OUTPUT/GMT/CEA_all-"//trim(adjustl(x))//".tex",STATUS="replace")
      write(4996,*)" "
      close(4996)
      open(unit=499,file="OUTPUT/GMT/CEAALL-"//trim(adjustl(x))//".txt",STATUS="replace")
      write(499,*)"0 0.0 0.0"
      close(499)
      open(unit=498,file="OUTPUT/GMT/CEAALLt-"//trim(adjustl(x))//".txt",STATUS="replace")
      write(498,*)"0 0.0 0.0 10 0 5 LM CEA"
      close(498)
    endif
    close(4999)
    ! -----------------------------------------------------------------    .

    ! -----------------------------------------------------------------    .
5007 format (a,f15.2,a,f15.2,a,f9.2,a,f9.1,a,f9.4,a,f9.3,a,f9.4,a,f9.3,a,f9.2,a,f9.2,a)
    ! -----------------------------------------------------------------    .
    deallocate(papriori)
    ! -----------------------------------------------------------------    .
  end subroutine dist_apriori

END MODULE sub_param



! *********************************************************************    .
! *********************************************************************    .


