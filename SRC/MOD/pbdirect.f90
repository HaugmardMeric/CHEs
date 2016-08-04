! Librairie de subroutines pour le problème direct
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE pb_direct

    use modparam

    implicit none

    private

    public  :: tempsTheoDirectone, tempsTheoDirect, tempsTheoDirectone_AUTRE
    public  :: Wadatiplot, chatelainplot
    public  :: pPn_sSn, reflechie2, reflechie, refracte, refracte_mohovar, directe ! à completer ...


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine tempsTheoDirect(nbtps,p,D,critique,acentroid)
    ! -----------------------------------------------------------------    .
    ! Calcul les temps théoriques des arrivées des ondes pour tous les séismes
    ! -----------------------------------------------------------------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: nbtps(nbseismes)
    type(parametres), intent (inout) :: p                                  ! paramètres d'inv.
    type(dataall), intent (inout) :: D(nbseismes)                          ! données
    logical, intent (out) :: critique                                      ! .true. si distance hypo + 5 km < distance hypo critique pour la réfraction
    type (amoho_centroid), intent (in) :: acentroid
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j
    type(parametre) :: param
    ! -----------------------------------------------------------------    .
    do i=1, nbseismes
      do j=1,nbtps(i)
        call mvPall_2_P1(param,p,i)
        call tempsTheoDirectone(param,D(i)%datatps(j),critique,acentroid)
        call mvP1_2_Pall(p,param,i)
      enddo
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine tempsTheoDirect

  ! -------------------------------------------------------------------    .

  subroutine tempsTheoDirectone(param,datatemps,critique,acentroid)
    ! -------                                                  --------    .mh
    ! Calcul les temps théoriques des arrivées des ondes
    ! pour un couple staion-hypocentre sur un seul séisme
    ! -------                                                  --------    .
    use typetemps
    use distance_epi
    use time
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent (in) :: param                                  ! paramètres d'inv.
    type(dataone), intent (inout) :: datatemps                             ! données
    logical, intent (out) :: critique                                      ! .true. si distance hypo + 5 km < distance hypo critique pour la réfraction
    type (amoho_centroid), intent (in) :: acentroid
    ! -----------------------------------------------------------------    . calcul de distance épicentrale
    call dellipsgc(datatemps%sta%lat,datatemps%sta%lon,param%Lat,param%Lon,datatemps%depi,datatemps%baz)
    ! -------                                                  --------    .
    ! calcul de distance hypocentrale avec prise en compte de l'altitude de la station
    ! -------                                                  --------    . niveau zéro : celui de la mer
    datatemps%dhypo = sqrt(datatemps%depi**2.0_wr + (param%Zhypo+datatemps%sta%alti/1000.0_wr)**2.0_wr )
    ! -------                                                  --------    .
    ! calcul de temps de parcours des ondes théoriques et les temps théoriques d'arrivées des ondes
    ! -------                                                  --------    .
    critique=.false.                                                       ! cas normal
    if (datatemps%typeonde.eq.'G') then                                    ! ondes directes
      ! -------                                                --------    .
      call directe(param,datatemps%dhypo,datatemps%tpsparcS,datatemps%tpsparcP)
      ! -------                                                --------    .
    elseif (datatemps%typeonde.eq.'N') then                                ! ondes réfractées
      ! if(datatemps%depi.gt.2500.0_wr) datatemps%coefS=max(3,datatemps%coefS)
      ! car Vp/Vs constant pour le profil, mais peut diverger pour les réfracté à longue distance
      ! -------                                                --------    .
      if(FLAG_non_tabulaire) then                                          ! moho incliné
        call refracte_mohovar(acentroid,param,datatemps%sta,datatemps%dhypo,datatemps%tpsparcS,datatemps%tpsparcP, &
        alti = datatemps%sta%alti, dcritique = datatemps%dcritiqueH)
      else
        call refracte(param,datatemps%dhypo,datatemps%tpsparcS,datatemps%tpsparcP, &
        alti = datatemps%sta%alti, dcritique = datatemps%dcritiqueH)       ! moho non incliné
      endif
      ! -------                                                --------    .
      if((datatemps%dhypo+5.0_wr).lt.datatemps%dcritiqueH) then
        critique=.true. ! distance hypo + 5 km < distance hypo critique pour la réfraction
      endif
    else
      write(*,*)'problème dans tempsTheoDirectone : onde ni directe ni réfractée ... ? ',datatemps%typeonde
      stop
    endif
    ! -------                                                  --------    .
    ! calcul du temps d'arrivée théorique absolu
    ! -------                                                  --------    .
    datatemps%tpsTh%date=param%Tzero%date                                  ! même date
    datatemps%tpsTh%secP = param%Tzero%sec + datatemps%tpsparcP            ! change les secondes
    datatemps%tpsTh%secS = param%Tzero%sec + datatemps%tpsparcS            ! change les secondes
    call basetime(datatemps%tpsTh)                                         ! reste en base 60/12/365 ... pour P et S
    ! -------                                                  --------    .
    ! calcul de la difference de temps d'arrivée entre données et modèle
    ! -------                                                  --------    .
    call difftime(datatemps%dTP,datatemps%dTS,datatemps%tpsR,datatemps%tpsTh)
    if(datatemps%andS.eq.'X') datatemps%dTS=0.0_wr                         ! si il n'existe pas d'ondes S
    if(IsNaN(datatemps%dTP)) then
      write(*,*)'problème dans tempsTheoDirectone 1 : datatemps%dTP = NaN'
      write(*,*)param
      write(*,*)datatemps
      stop
    endif
    if(IsNaN(datatemps%dTS)) then
      write(*,*)'problème dans tempsTheoDirectone 2 : datatemps%dTS = NaN'
      write(*,*)param
      write(*,*)datatemps
      stop
    endif
    ! -----------------------------------------------------------------    .
  end subroutine tempsTheoDirectone

    ! -----------------------------------------------------------------    .

  subroutine tempsTheoDirectone_AUTRE(param,datatemps,pdfmoho,critique,m)
    ! -------                                                  --------    .mh
    ! Calcul les temps théoriques des arrivées des ondes
    ! pour un couple staion-hypocentre sur un seul séisme
    ! selon un modèle de terre tabulaire à n couches
    ! -------                                                  --------    .
    use typetemps
    use distance_epi
    use time
    use ray_tracing
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent (in) :: param                                  ! paramètres d'inv.
    type(dataone), intent (inout) :: datatemps                             ! données
    logical, intent (out) :: critique                                      ! .true. si distance hypo + 5 km < distance hypo critique pour la réfraction
    character (LEN=1), intent(in), optional :: m
    real(KIND=wr), intent (out) :: pdfmoho
    ! -------                                                  --------    .
    real (kind=wr) :: altidudesta,tdirectP,trefP,tdirectS,trefS
    ! -----------------------------------------------------------------    .
    altidudesta=-datatemps%sta%alti/1000._wr
    ! -------                                                  --------    . calcul de distance épicentrale
    call dellipsgc(datatemps%sta%lat,datatemps%sta%lon,param%Lat,param%Lon,datatemps%depi,datatemps%baz)
    ! -------                                                  --------    .
    ! calcul de temps de parcours des ondes théoriques et les temps théoriques d'arrivées des ondes
    ! -------                                                  --------    .
    if(present(m)) then
      call tracerays(datatemps%depi,datatemps%dhypo,datatemps%dcritiqueH,param%Zhypo,altidudesta, &
        param%Lon,param%Lat,pdfmoho,tdirectP,trefP,tdirectS,trefS,m)
    else
      call tracerays(datatemps%depi,datatemps%dhypo,datatemps%dcritiqueH,param%Zhypo,altidudesta, &
        param%Lon,param%Lat,pdfmoho,tdirectP,trefP,tdirectS,trefS)
    endif
    ! -------                                                  --------    .
    critique=.false.                                                       ! cas normal
    ! -------                                                  --------    .
    if (datatemps%typeonde.eq.'G') then                                    ! ondes directes
      datatemps%tpsparcP=tdirectP
      datatemps%tpsparcS=tdirectS
    elseif (datatemps%typeonde.eq.'N') then                                ! ondes réfractées
      ! if(datatemps%depi.gt.2500.0_wr) datatemps%coefS=max(3,datatemps%coefS)
      ! car Vp/Vs constant pour le profil, mais peut diverger pour les réfracté à longue distance
      datatemps%tpsparcP=trefP
      datatemps%tpsparcS=trefS
      if((datatemps%dhypo+5.0_wr).lt.datatemps%dcritiqueH) then
        critique=.true. ! distance hypo + 5 km < distance hypo critique pour la réfraction
      endif
    else
      write(*,*)'problème dans tempsTheoDirectone_AUTRE : onde ni directe ni réfractée ... ? '
      stop
    endif
    ! -------                                                  --------    .
    ! calcul du temps d'arrivée théorique absolu
    ! -------                                                  --------    .
    datatemps%tpsTh%date=param%Tzero%date                                  ! même date
    datatemps%tpsTh%secP = param%Tzero%sec + datatemps%tpsparcP            ! change les secondes
    datatemps%tpsTh%secS = param%Tzero%sec + datatemps%tpsparcS            ! change les secondes
    call basetime(datatemps%tpsTh)                                         ! reste en base 60/12/365 ... pour P et S
    ! -------                                                  --------    .
    ! calcul de la difference de temps d'arrivée entre données et modèle
    ! -------                                                  --------    .
    call difftime(datatemps%dTP,datatemps%dTS,datatemps%tpsR,datatemps%tpsTh)
    ! -----------------------------------------------------------------    .
    if(datatemps%andS.eq.'X') datatemps%dTS=0.0_wr                         ! si il n'existe pas d'ondes S
    if(IsNaN(datatemps%dTP)) then
      write(*,*)'problème dans tempsTheoDirectone_AUTRE 1 : datatemps%dTP = NaN'
      write(*,*)param
      write(*,*)datatemps
      stop
    endif
    if(IsNaN(datatemps%dTS)) then
      write(*,*)'problème dans tempsTheoDirectone_AUTRE 2 : datatemps%dTS = NaN'
      write(*,*)param
      write(*,*)datatemps
      stop
    endif
    ! -----------------------------------------------------------------    .
  end subroutine tempsTheoDirectone_AUTRE

    ! -----------------------------------------------------------------    .

  subroutine Wadatiplot(nbtps,D,param_best,vpvs,a,R2,XY,nb,sig,OK,atype)
    ! -------                                                  --------    .mh
    ! Calcul la regression sur le Wadati plot (Ts-Tp en fonction de Tp-To)
    ! Calcul une bonne estimation de VpVS si les autres paramètres sont déja pas mauvais
    ! wadati (1933)
    ! -------                                                  --------    .
    use typetemps
    use time
    use statistiques, only : correlationaffpond
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(in) :: nbtps(nbseismes)
    type(dataall), intent(in) :: D(nbseismes)
    type(parametres), intent(in) :: param_best
    ! -------                                                  --------    .
    integer (KIND=wi), parameter :: taille=2500
    real (KIND=wr), intent(out), optional :: vpvs
    real (KIND=wr), intent(out), optional :: a
    real (KIND=wr), intent(out), optional :: R2
    real (KIND=wr), intent(out), optional :: XY(taille,3)
    real (KIND=wr), intent(out), optional :: sig(taille,2)
    integer (KIND=wi), intent(out), optional :: nb
    logical, intent(out), optional :: OK
    character(len=1), intent(in), optional :: atype
    ! -------                                                  --------    .
    integer (KIND=wi) :: i,j,n
    type(date_sec) :: one_tps_1,one_tps_2
    real (KIND=wr) :: XY_coef(taille,3),a_coef,R2_coef
    real (KIND=wr), dimension(:,:), allocatable :: XYbis
    real (KIND=wr) :: sig1(taille,2)
    logical :: test,tropval
    ! -----------------------------------------------------------------    . initialisation
    do i=1,taille
      XY_coef(i,1)=-1.0_wr
      XY_coef(i,2)=-1.0_wr
      XY_coef(i,3)=-1.0_wr
      sig1(i,1)=-1.0_wr
      sig1(i,2)=-1.0_wr
    enddo
    ! -----------------------------------------------------------------    .
    n=0
    tropval=.true.
    do j=1,nbseismes
      do i=1,nbtps(j)
        if(tropval) then                                                   ! n < taille=2500
           if(D(j)%datatps(i)%andS.eq.'S') then
            ! -------                                          --------    . onde G, N ou les 2
            if(present(atype)) then
              if (D(j)%datatps(i)%typeonde.eq.atype) then
                test = .true.
              else
                test=.false.
              endif
            else
              test = .true.
            endif
            if (test) then
              n=n+1
              ! -------                                        --------    . onde P
              one_tps_1%date = D(j)%datatps(i)%tpsR%date
              one_tps_1%sec = D(j)%datatps(i)%tpsR%secP
              call difftime(XY_coef(n,1),one_tps_1,param_best%Tzero(j))
              ! -------                                        --------    . onde S
              one_tps_2%date = D(j)%datatps(i)%tpsR%date
              one_tps_2%sec = D(j)%datatps(i)%tpsR%secS
              call difftime(XY_coef(n,2),one_tps_2,one_tps_1)
              ! -------                                        --------    . incertitudes sur les données -> pour les figures
              sig1(n,1)=D(j)%datatps(i)%sigP
              sig1(n,2)=D(j)%datatps(i)%sigS                               ! GMT_wadati prend en compte la propagation de erreurs
              ! -------                                        --------    .
              XY_coef(n,3)= D(j)%datatps(i)%ws * D(j)%datatps(i)%wp
              ! -------                                        --------    .
              if (n.eq.(taille)) then
                !do o=1,taille
                  !write(*,*)o,XY_coef(o,1),XY_coef(o,2),sig1(o,1), sig1(o,2)
                !enddo
                write(*,*)'problème dans Wadatiplot : trop de données : ',j,i,n
                tropval=.false.
              endif
              ! -------                                        --------    .
            endif
          endif
        endif
      enddo
    enddo
    ! -------                                                  --------    .
    allocate(XYbis(n,3))
    do i=1,n
      XYbis(i,1) = XY_coef(i,1)
      XYbis(i,2) = XY_coef(i,2)
      XYbis(i,3) = XY_coef(i,3)
    enddo
    ! -------                                                  --------    .
    call correlationaffpond(a_coef,R2_coef,n,XYbis)
    ! -------                                                  --------    .
    if (present(a)) a = a_coef
    if (present(R2)) R2 = R2_coef
    if (present(vpvs)) vpvs = 1.0_wr+a_coef
    if (present(XY)) XY=XY_coef
    if (present(nb)) nb = n
    if (present(sig)) sig = sig1
    ! -------                                                  --------    .
    deallocate(XYbis)
    if (n.gt.2) then
      if (present(OK)) OK = .false.
    else
      if (present(OK)) OK = .true.
    endif
    ! -------                                                  --------    .
    if (present(vpvs)) then
      if (IsNaN(vpvs)) then
        write(*,*)'problème dans Wadatiplot : IsNaN(vpvs)',a_coef
      endif
    endif
    ! -----------------------------------------------------------------    .
  end subroutine Wadatiplot

    ! -----------------------------------------------------------------    .

  subroutine chatelainplot(nbtps,D,vpvs,a,R2,XY,nb,sig,OK,atype,nom_sta)
    ! -------                                                  --------    .mh
    ! Calcul la regression sur le châtelain plot
    ! avec : Ts1-Ts2 en fonction de Tp1-Tp2
    ! -------                                                  --------    .
    ! le diagramme de Châtelain (ou Wadati modifié) est indépendant des parametres
    ! et peux ainsi identifier une erreur de pointé sur les ondes
    ! -------                                                  --------    .
    ! Châtelain (1978) : "Etude fine de la sismicité en zone de collision continentale au moyen
    ! d'un réseau de stations portables : la région Hindu-Kush Pamir" these de doctorat
    ! -------                                                  --------    .
    use typetemps
    use time
    use statistiques, only : correlationaffpond
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(in) :: nbtps(nbseismes)
    type(dataall), intent(in) :: D(nbseismes)
    ! -------                                                  --------    .
    integer (KIND=wi), parameter :: taille=5000
    ! -------                                                  --------    .
    real (KIND=wr), intent(out), optional :: vpvs
    real (KIND=wr), intent(out), optional :: a
    real (KIND=wr), intent(out), optional :: R2
    real (KIND=wr), intent(out), optional :: XY(taille,3)
    real (KIND=wr), intent(out), optional :: sig(taille,2)
    integer (KIND=wi), intent(out), optional :: nb
    logical, intent(out), optional :: OK
    character(len=1), intent(in), optional :: atype
    character(len=4), intent(out), optional :: nom_sta(taille,2)
    ! -------                                                  --------    .
    integer (KIND=wi) :: i,j,k,l,m,n
    type(date_sec) :: one_tps_1,one_tps_2
    real (KIND=wr) :: val1,val2,XY_coef(taille,3),a_coef,R2_coef
    real (KIND=wr), dimension(:,:), allocatable :: XYbis
    real (KIND=wr) :: sig1(taille,2)
    logical :: test, deja, tropval
    ! -----------------------------------------------------------------    .
    if(present(nom_sta)) nom_sta(:,:)='123_'
    ! -------                                                  --------    .
    do i=1,taille
      XY_coef(i,1)=-1.0_wr
      XY_coef(i,2)=-1.0_wr
      XY_coef(i,3)=-1.0_wr
      sig1(i,1)=-1.0_wr
      sig1(i,2)=-1.0_wr
    enddo
    n=0
    tropval=.true.
    do j=1,nbseismes
      do i=1,nbtps(j)
        do k=1,nbtps(j)
          if(tropval) then                                                 ! n < taille=1000

            ! -------                                          --------    . onde G, N ou les 2
            if(present(atype)) then
              if ((D(j)%datatps(i)%typeonde.eq.atype).and.(D(j)%datatps(k)%typeonde.eq.atype)) then
                test = .true.
              else
                test=.false.
              endif
            else
              if (D(j)%datatps(i)%typeonde.eq.D(j)%datatps(k)%typeonde) then
                test = .true.
              else
                test = .false.
              endif
            endif
            if ((test).and.(D(j)%datatps(i)%andS=='S').and.(D(j)%datatps(k)%andS=='S')) then
              ! -------                                        --------    . onde P et S
              n=n+1
              one_tps_1%date = D(j)%datatps(i)%tpsR%date
              one_tps_1%sec = D(j)%datatps(i)%tpsR%secP
              one_tps_2%date = D(j)%datatps(k)%tpsR%date
              one_tps_2%sec = D(j)%datatps(k)%tpsR%secP
              call difftime(val1,one_tps_1,one_tps_2)
              one_tps_1%date = D(j)%datatps(i)%tpsR%date
              one_tps_1%sec = D(j)%datatps(i)%tpsR%secS
              one_tps_2%date = D(j)%datatps(k)%tpsR%date
              one_tps_2%sec = D(j)%datatps(k)%tpsR%secS
              call difftime(val2,one_tps_1,one_tps_2)
              ! -------                                        --------    . déja présent ?
              deja=.true.
              do l=1,n
                if((abs(val1)==XY_coef(l,1)).and.(abs(val2)==XY_coef(l,2))) then
                  deja=.false.
                endif
              enddo
              if(deja)then
                XY_coef(n,1)=abs(val1)
                XY_coef(n,2)=abs(val2)
                ! -------                                    --------    . incertitudes sur les données -> pour les figures
                if((abs(val1)==0.0_wr).and.(abs(val2)==0.0_wr)) then
                  sig1(n,1)=0.0_wr
                  sig1(n,2)=0.0_wr
                else
                  sig1(n,1)=sqrt(D(j)%datatps(i)%sigP**2.0_wr+D(j)%datatps(k)%sigP**2.0_wr)
                  sig1(n,2)=sqrt(D(j)%datatps(i)%sigS**2.0_wr+D(j)%datatps(k)%sigS**2.0_wr)
                endif
                ! -------                                      --------    .
                XY_coef(n,3)=(D(j)%datatps(i)%ws+D(j)%datatps(i)%wp+D(j)%datatps(k)%ws+D(j)%datatps(k)%wp)/4.0_wr
                if(present(nom_sta)) then
                  nom_sta(n,1)=D(j)%datatps(i)%sta%staname
                  nom_sta(n,2)=D(j)%datatps(k)%sta%staname
                endif
              else
                n=n-1
              endif
              ! -------                                        --------    .
              if (n.eq.(taille)) then
                write(*,*)'problème dans chatelainplot : trop de données : ',j,i,m,k,n
                !do l=1,taille
                  !write(*,*)l,XY_coef(l,1),XY_coef(l,2),sig1(l,1), sig1(l,2)
                !enddo
                tropval=.false.
              endif
              ! -------                                        --------    .
            endif
          endif
        enddo
      enddo
    enddo
    ! -------                                                  --------    .
    allocate(XYbis(n,3))
    do i=1,n
      XYbis(i,1) = XY_coef(i,1)
      XYbis(i,2) = XY_coef(i,2)
      XYbis(i,3) = XY_coef(i,3)
    enddo
    ! -------                                                  --------    .
    call correlationaffpond(a_coef,R2_coef,n,XYbis)
    ! -------                                                  --------    .
    if (present(a)) a = a_coef
    if (present(R2)) R2 = R2_coef
    if (present(vpvs)) vpvs = a_coef
    if (present(XY)) XY=XY_coef
    if (present(nb)) nb = n
    if (present(sig)) sig = sig1
    ! -------                                                  --------    .
    deallocate(XYbis)
    if (n.gt.2) then
      if (present(OK)) OK = .false.
    else
      if (present(OK)) OK = .true.
    endif
    ! -------                                                  --------    .
    if (present(vpvs)) then
      if (IsNaN(vpvs)) then
        write(*,*)'problème dans chatelainplot : IsNaN(vpvs)',a_coef
      endif
    endif
    ! -----------------------------------------------------------------    .
  end subroutine chatelainplot

    ! -----------------------------------------------------------------    .

  subroutine directe(param,dishypo,Tps,Tpp)
    ! -------                                                  --------    .mh
    ! calcul du temps d'arrivée des ondes P et S directes (Pg et Sg)
    ! pour 2 distances hypocentrales
    ! modèle tabulaire
    ! VpVs constant le long du profile
    ! pas de séismes sous le moho
    ! -------                                                  --------    .
    use typetemps, only : parametre
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent(in) :: param                                   ! parametres du modèle
    real(kind=wr), intent(in) :: dishypo                                   ! distance hypocentrale
    real(kind=wr), intent(out) :: Tps,Tpp                                  ! temps des ondes Pg et Sg
    ! -----------------------------------------------------------------    .
    if(param%Zhypo.lt.(param%Zmoho-0.1_wr)) then                           ! séisme au dessus du moho !
      Tps= dishypo/param%VC*param%VpVs
      Tpp= dishypo/param%VC
    else
      Tps= 1.e9_wr
      Tpp= 1.e9_wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine directe

    ! -----------------------------------------------------------------    .

  subroutine refracte(param,dishypo,Tps,Tpp,alti,dcritique)
    ! -------                                                  --------    .mh
    ! calcul du temps d'arrivée des ondes P et S réfractées (Pn et Sn)
    ! pour une distance hypocentrale
    ! modèle tabulaire
    ! VpVs constant le long du profile
    ! -------                                                  --------    .
    use typetemps, only : parametre
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent(in) :: param                                   ! parametres du modèle
    real(kind=wr), intent(in) :: dishypo                                   ! distance hypocentrale
    real(kind=wr), intent(out) :: Tps,Tpp                                  ! temps des ondes Pn et Sn
    real(kind=wr), intent(in), optional :: alti                            ! altitude de la station (m)
    real(kind=wr), intent(out), optional :: dcritique                      ! distance hypocentrale critique (à partir de laquelle les premières réfractions ont lieu)
    ! -------                                                  --------    .
    real(kind=wr) :: anglei,depi,dist_c,dist_m,altista,dc
    ! -----------------------------------------------------------------    .
    if(param%Zhypo.lt.(param%Zmoho-0.1_wr)) then                           ! séisme au dessus du moho !
      anglei = asin(param%VC/param%VM)
      if(IsNaN(anglei)) then
        write(*,*)'problème dans refracte : paramètre de rai = NaN car VC > VM '
        stop
      endif
      ! -------                                                --------    .
      if (present(alti)) then
        altista = alti/1000.0_wr
      else
        altista = 0.0_wr
      endif
      depi = sqrt(dishypo**2.0_wr - (param%Zhypo+altista)**2.0_wr)
      ! -------                                                --------    .
      dc = (param%Zmoho-param%Zhypo)*tan(anglei) + &
         (param%Zmoho+altista)*tan(anglei)
      dist_m = depi - dc
      dist_c = (param%Zmoho-param%Zhypo)/cos(anglei) + &
             (param%Zmoho+altista)/cos(anglei)
      Tpp = dist_c/param%VC + dist_m/param%Vm
      Tps = param%VpVs*Tpp
      dc = sqrt(dc**2.0_wr + (param%Zhypo+altista)**2.0_wr )
      ! -------                                                --------    .
      if (present(dcritique)) dcritique = dc
    else
      Tps= 1.e9_wr
      Tpp= 1.e9_wr
      if (present(dcritique)) dcritique = 1.e9_wr
    endif
    ! -----------------------------------------------------------------    .
  end subroutine refracte

    ! -----------------------------------------------------------------    .

  subroutine refracte_mohovar(mohocentroid,param,sta,dishypo,Tps,Tpp,alti,dcritique,pfd)
    ! -------                                                  --------    .mh
    ! calcul du temps d'arrivée des ondes P et S réfractées (Pn et Sn)
    ! pour une distance hypocentrale
    ! - > modèle non tabulaire avec mohocentroid
    ! VpVs constant le long du profile
    ! -------                                                  --------    .
    use typetemps
    use distance_epi
    ! -------                                                  --------    .
    implicit none
    type (amoho_centroid), intent (in) :: mohocentroid
    type(parametre), intent(in) :: param                                   ! parametres du modèle
    type(stations), intent (in) :: sta
    real(kind=wr), intent(in) :: dishypo                                   ! distance hypocentrale
    real(kind=wr), intent(out) :: Tps,Tpp                                  ! temps des ondes Pn et Sn
    real(kind=wr), intent(in), optional :: alti                            ! altitude de la station (m)
    real(kind=wr), intent(out), optional :: dcritique                      ! distance hypocentrale critique (à partir de laquelle les premières réfractions ont lieu)
    real(kind=wr), intent(out), optional :: pfd                            ! profondeur du moho sous le séisme
    ! -------                                                  --------    .
    real(kind=wr) :: anglei,depi,dist_c,dist_m,altista,dc
    real(kind=wr) :: pfdmohoSTA,pfdmohoEvent,angleap,pfd1,pfd2
    real(kind=wr) :: dSTAlon,dSTAlat,dEVlon,dEVlat
    ! -----------------------------------------------------------------    .
    anglei = asin(param%VC/param%VM)
    if(IsNaN(anglei)) then
      write(*,*)'problème dans refracte_mohovar : paramètre de rai = NaN car VC > VM '
      stop
    endif
    ! -------                                                  --------    .
    if (present(alti)) then
      altista = alti/1000.0_wr
    else
      altista = 0.0_wr
    endif
    depi = sqrt(dishypo**2.0_wr - (param%Zhypo+altista)**2.0_wr)
    ! -------                                                  --------    .
    ! vecteur normal du moho (lon, lat, z) en km
    call dellipsgc(mohocentroid%latC,mohocentroid%lonC,sta%lat,mohocentroid%lonC,dSTAlat)
    if (sta%Lat.gt.mohocentroid%latC) dSTAlat=-dSTAlat
    call dellipsgc(mohocentroid%latC,mohocentroid%lonC,mohocentroid%latC,sta%lon,dSTAlon)
    if (mohocentroid%lonC.gt.sta%Lon) dSTAlon=-dSTAlon
    ! prod scalaire : sta(dSTAlon,dSTAlat,pfdmohoSTA) . mohocentroid(alpha,beta,gamma) = 0
    pfdmohoSTA=param%Zmoho-dSTAlon*mohocentroid%alph/mohocentroid%gamma-dSTAlat*mohocentroid%beta/mohocentroid%gamma
    call dellipsgc(mohocentroid%latC,mohocentroid%lonC,param%Lat,mohocentroid%lonC,dEVlat)
    if (param%Lat.gt.mohocentroid%latC) dEVlat=-dEVlat
    call dellipsgc(mohocentroid%latC,mohocentroid%lonC,mohocentroid%latC,param%Lon,dEVlon)
    if (mohocentroid%lonC.gt.param%Lon) dEVlon=-dEVlon
    ! prod scalaire : event(dEVlon,dEVlat,0) . mohocentroid(alpha,beta,gamma) = 0
    pfdmohoEvent=param%Zmoho-dEVlon*mohocentroid%alph/mohocentroid%gamma-dEVlat*mohocentroid%beta/mohocentroid%gamma
    ! -------                                                  --------    .
    angleap=atan(abs(pfdmohoEvent-pfdmohoSTA)/depi)                        ! pendage apparent du moho
    ! -------                                                  --------    .
    pfd1=(pfdmohoEvent-param%Zhypo)*cos(angleap)                           ! distance entre moho et séisme, perpendiculaire au moho
    pfd2=(pfdmohoSTA+altista)*cos(angleap)                                 ! distance entre moho et station, perpendiculaire au moho
    dc = pfd1*tan(anglei) + pfd2*tan(anglei)                               ! distance parcourue dans la croûte projetée sur le moho
    dist_m = depi/cos(angleap) - dc                                        ! distance parcourue dans le manteau
    dist_c = pfd1/cos(anglei) + pfd2/cos(anglei)                           ! distance parcourue dans la croûte
    Tpp = dist_c/param%VC + dist_m/param%Vm
    Tps = param%VpVs*Tpp
    dc = sqrt(dc**2.0_wr + (param%Zhypo+altista)**2.0_wr) * cos(anglei)    ! distance dcritique pour les premieres réfractées
    ! -------                                                  --------    .
    if (present(dcritique)) dcritique = dc
    if (present(pfd)) pfd = pfdmohoEvent
    ! -----------------------------------------------------------------    .
  end subroutine refracte_mohovar

    ! -----------------------------------------------------------------    .

  subroutine reflechie(param,dishypo,Tps,Tpp,alti)
    ! -------                                                  --------    .mh
    ! calcul du tmps d'arrivée des ondes P et S réflechies (PmP et SmS)
    ! pour une distance hypocentrale
    ! modèle tabulaire
    ! VpVs constant le long du profile
    ! -------                                                  --------    .
    use typetemps, only : parametre
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent(in) :: param                                   ! parametres du modèle
    real(kind=wr), intent(in) :: dishypo                                   ! distance hypocentrale
    real(kind=wr), intent(out) :: Tps,Tpp                                  ! temps des ondes PmP et SmS
    real(kind=wr), intent(in), optional :: alti                            ! altitude de la station (m)
    ! -------                                                  --------    .
    real(kind=wr) :: depi,dc1,dc2, altista
    ! -----------------------------------------------------------------    .
    if (present(alti)) then
      altista = alti/1000.0_wr
    else
      altista = 0.0_wr
    endif
    ! -------                                                  --------    .
    if((param%Zhypo+altista).lt.dishypo) then
      depi = sqrt(dishypo**2.0_wr-(param%Zhypo+altista)**2.0_wr)
      ! -------                                                --------    .
      dc1 = depi * 1.0_wr / &
            ((param%Zmoho-param%Zhypo) / (param%Zmoho+altista) + 1.0_wr)
      dc2 = depi - dc1
      dc1 = sqrt(dc1**2.0_wr + (param%Zmoho-param%Zhypo)**2.0_wr)
      dc2 = sqrt(dc2**2.0_wr + (param%Zmoho+altista)**2.0_wr)
      Tpp = (dc1 + dc2) / param%VC
      Tps= Tpp*param%VpVs
    else
      write(*,*)'problème dans reflechie : réflexion impossible '
    endif
    ! -----------------------------------------------------------------    .
  end subroutine reflechie

    ! -----------------------------------------------------------------    .

  subroutine reflechie2(param,dishypo,Tps,Tpp,alti)
    ! -------                                                  --------    .mh
    ! calcul du tmps d'arrivée des ondes P et S réflechies (PmP2 et SmS2)
    ! pour une distance hypocentrale
    ! modèle tabulaire
    ! VpVs constant le long du profile
    ! -------                                                  --------    .
    use typetemps, only : parametre
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent(in) :: param                                   ! parametres du modèle
    real(kind=wr), intent(in) :: dishypo                                   ! distance hypocentrale
    real(kind=wr), intent(out) :: Tps,Tpp                                  ! temps des ondes PmP2 et SmS2
    real(kind=wr), intent(in), optional :: alti                            ! altitude de la station (m)
      ! -------                                                --------    .
    real(kind=wr) :: depi,dc1,dc2,altista
    ! -----------------------------------------------------------------    .
    if (present(alti)) then
      altista = alti/1000.0_wr
    else
      altista = 0.0_wr
    endif
    ! -------                                                  --------    .
    if((param%Zhypo+altista).lt.dishypo) then
      depi = sqrt(dishypo**2.0_wr-(param%Zhypo+altista)**2.0_wr)
      ! -------                                                --------    .
      dc1 = depi * 1.0_wr / ((param%Zmoho-param%Zhypo) / &
            (3.0_wr*(param%Zmoho+altista))+1.0_wr)
      dc2 = (depi - dc1)/3.0_wr
      dc1 = sqrt(dc1**2.0_wr+(param%Zmoho-param%Zhypo)**2.0_wr)
      dc2 = sqrt(dc2**2.0_wr+(param%Zmoho+altista)**2.0_wr)
      Tpp = (dc1 + 3.0_wr*dc2)/param%VC
      Tps= Tpp*param%VpVs
    else
      write(*,*)'problème dans reflechie2 : réflexion impossible '
    endif
    ! -----------------------------------------------------------------    .
  end subroutine reflechie2

    ! -----------------------------------------------------------------    .

  subroutine pPn_sSn(param,dishypo,Tps,Tpp,alti,dcritique)
    ! -------                                                  --------    .mh
    ! calcul du temps d'arrivée des ondes P et S réfractées (pPn et sSn)
    ! pour une distance hypocentrale
    ! modèle tabulaire
    ! VpVs constant le long du profile
    ! -------                                                  --------    .
    use typetemps, only : parametre
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent(in) :: param                                   ! parametres du modèle
    real(kind=wr), intent(inout) :: dishypo                                ! distance hypocentrale
    real(kind=wr), intent(out) :: Tps,Tpp                                  ! temps des ondes pPn et sSn
    real(kind=wr), intent(in), optional :: alti                            ! altitude de la station (m)
    real(kind=wr), intent(out), optional :: dcritique                      ! distance hypocentrale critique (à partir de laquelle les premières réfractions ont lieu)
    ! -------                                                  --------    .
    real(kind=wr) :: anglei,depi,dist_c,dist_m,altista,dc
    ! -----------------------------------------------------------------    .
    if (present(alti)) then
      altista = alti/1000.0_wr
    else
      altista = 0.0_wr
    endif
    ! -----------------------------------------------------------------    .
    anglei = asin(param%VC/param%VM)
    if(IsNaN(anglei)) then
      write(*,*)'problème dans pPn_pSn : paramètre de rai = NaN car VC > VM '
      stop
    endif
      ! -------                                                --------    .
    if (dishypo==0.0_wr) then
      depi = tan(anglei)*param%Zhypo + (2.0_wr*param%Zmoho+1.5_wr*altista)*tan(anglei)
      dishypo = sqrt(depi**2.0_wr + param%Zhypo**2.0_wr)
    else
      depi = sqrt(dishypo**2.0_wr-(param%Zhypo)**2.0_wr)
    endif
      ! -------                                                --------    .
    dc = tan(anglei)*param%Zhypo + (2.0_wr*param%Zmoho+1.5_wr*altista)*tan(anglei)
    dist_m = depi  - dc
    dist_c = (param%Zhypo)/cos(anglei) + (2.0_wr*param%Zmoho+1.5_wr*altista)/cos(anglei)
    Tpp = dist_c/param%VC + dist_m/param%Vm
    Tps = param%VpVs*Tpp
    dc = sqrt(dc**2.0_wr + (param%Zhypo+altista)**2.0_wr )
    if (present(dcritique)) dcritique = dc
    ! -----------------------------------------------------------------    .
  end subroutine pPn_sSn

END MODULE pb_direct



! *********************************************************************    .
! *********************************************************************    .


