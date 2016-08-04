! déclaration de l'ensemble des structures dérivées utilisées dans le programme CHE
! septembre 2013
! **********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         ---------    .
! **********************************************************************    .
! ----------------------------------------------------------------------    .


MODULE typetemps

    use modparam

    implicit none

    private

    ! type :
    public  :: date_min, date_sec, date_secPS
    public  :: stations
    public  :: dataone, dataall
    public  :: pond
    public  :: parametre, parametres
    public  :: parametresinv, paramisfit
    public  :: densityplot_one, densityplot
    public  :: fcout
    public  :: accept
    public  :: seismes
    public  :: coldmoyval,coldmoy
    public  :: ellip
    public  :: residus
    public  :: priorEPI
    public  :: amoho_centroid

    ! subroutines :
    public  :: mvP1_2_Pall, mvPall_2_P1
    public  :: vect2alph, alph2vect

! ----------------------------------------------------------------------    .

  ! -------                                                     --------    ! date
    TYPE date_min
        integer(KIND=wi) :: Jday,year,month,day,hour,min                    ! jours Julien, années, mois, jours du mois, heures et minutes en base 60, 24, 365, 12 ...
    END TYPE date_min
  ! -------                                                     --------    .
  ! -------                                                     --------    ! date en secondes
    TYPE date_sec
        type(date_min) :: date                                              ! en minute
        real(KIND=wr) :: sec                                                ! en secondes
    END TYPE date_sec
  ! -------                                                     --------    .
  ! -------                                                     --------    ! date en secondes (avec arrivées des ondes S et P)
    TYPE date_secPS
        type(date_min) :: date                                              ! en minute
        real(KIND=wr) :: secP,secS                                          ! en secondes
    END TYPE date_secPS
  ! -------                                                     --------    .
  ! -------                                                     --------    ! caractéristiques du réseau sismologique utilisé
    TYPE stations
        character(LEN=4) :: staname                                         ! nom station
        real(KIND=wr) :: lat,lon,alti                                       ! coordonnées station (degrée et m.a.s.l.)
        real(KIND=wr) :: res_Pg,res_Pn,res_Sg,res_Sn                        ! résidus à la station, onde P et S
    END TYPE stations
  ! -------                                                     --------    .
  ! -------                                                     --------    ! les DONNÉES, phase list (temps d'arrivées des ondes) pour 1 station, 1 séisme
    TYPE dataone
        type(stations) :: sta                                               ! données sur la station
        type(date_secPS) :: tpsR, tpsTh                                     ! temps d'arrivées des ondes R (réels) et Th (théorique)
        character(LEN=1) :: typeonde                                        ! "G" pour onde directe, "N" pour réfractée au moho
        integer(KIND=wi) :: coefP,coefS                                     ! coef qualité de lecture (0(best)-4(worst))
        character(LEN=1) :: andS                                            ! "S" si ondes S lues
        real(KIND=wr) :: dTP,dTS                                            ! différence de temps d'arrivées des ondes R (réels) et Th (théorique), ondes P et S
        real(KIND=wr) :: ws, wp                                             ! pondération des ondes P et S [0;1] (qualité et distance)
        real(KIND=wr) :: tpsparcP,tpsparcS                                  ! temps de parcours théorique des ondes P et S
        real(KIND=wr) :: depi,dhypo                                         ! distance épi- et hypocentrale en km
        real(KIND=wr) :: dcritiqueH                                         ! distance hypocentrale à partir de laquelle la réfraction commence, en km
        real(KIND=wr) :: baz                                                ! back-azimuth en degrée
        real(KIND=wr) :: sigP,sigS                                          ! écart-type sur les données en secondes
    END TYPE dataone
  ! -------                                                     --------    .
  ! -------                                                     --------    ! les DONNÉES (temps d'arrivées des ondes) pour tous les séismes
    TYPE dataall
        type(dataone), dimension(:), allocatable :: datatps                 ! les DONNÉES (temps d'arrivées des ondes) pour 1 séisme
    END TYPE dataall
  ! -------                                                     --------    .
  ! -------                                                     --------    ! coeficients de pondération pour chaque données
    TYPE pond
        real(KIND=wr) :: S_Pg, S_Sg, S_Pn, S_Sn                             ! somme des coef.
        integer(KIND=wi) :: nPg, nPn ,nSg, nSn                              ! nombre de données pour chaque
    END TYPE pond
  ! -------
  ! -------                                                     --------    ! ensemble des PARAMETRES de l'inversion pour plusieurs séismes
    TYPE parametres
        real(KIND=wr) :: VC,VM,Zmoho,VpVs                                   ! paramètres de terre
        real(KIND=wr) :: Lat(nbseismes),Lon(nbseismes),Zhypo(nbseismes)     ! paramètres des hypocentres
        type(date_sec) :: Tzero(nbseismes)                                  ! temps initial
    END TYPE parametres
  ! -------                                                     --------    .
  ! -------                                                     --------    ! ensemble des parametres de l'inversion pour 1 seul séisme
    TYPE parametre
        real(KIND=wr) :: VC,VM,Zmoho,VpVs,Lat,Lon,Zhypo                     !
        type(date_sec) :: Tzero                                             ! temps initial
    END TYPE parametre
  ! -------                                                     --------    .
  ! -------                                                     --------    . valeurs des parametres au cours de l'inversion
    TYPE parametresinv
        type(parametres) :: valNew,valOld,mini,maxi,ecartype                ! à l'itération actuelle, précedente, valeur minimale, maximale et écart type de la gaussienne d'échantillonage
        real(KIND=wr) :: centreX(nbseismes),centreY(nbseismes),Rayon(nbseismes),ec_horizontal
    END TYPE parametresinv
  ! -------                                                     --------    .
  ! -------                                                     --------    . un modèle (i.e. misfit + jeu de paramètres)
    TYPE paramisfit
        type(parametres) :: par                                             ! un jeu de paramètre
        real(KIND=wr) :: mis                                                ! valeur de la fonction coût
    END TYPE paramisfit
  ! -------                                                     --------    .
  ! -------                                                     --------    ! relecture modèles sélectionnés et calculs a posteriori
    TYPE densityplot_one
        real(KIND=wr), dimension(:), allocatable :: vec                     ! ensemble de valeurs prise lors de l'inverssion pour ce parametre
        character (LEN=30) :: char                                          ! chaîne de caractères pour les légendes GMT
        character (LEN=3) :: name                                           ! nom en trois lettre pour les nom des fichiers
        real(KIND=wr)themax,themin                                          ! min et max ou prior
        real(KIND=wr) :: delta                                              ! valeur de la disrcétisation pour le diagramme de densité, les histogrammes et le calcul du mode
        real(KIND=wr) :: best,mode,mediane                                  ! meilleur modèle, mode et médiane
        real(KIND=wr) :: moy_tot,moy_100,moy_1000,moy_10000                 ! moyenne totale puis des 100, 1000, et 10000 meilleurs modèles (plus petite fonction coût)
        real(KIND=wr) :: ec_tot,ec_100,ec_1000,ec_10000                     ! 1 écart type
        real(KIND=wr) :: moy_bestchaine,ec_bestchaine                       ! moyenne et écart type de l'ensemble du meilleur modèle de chaque chaîne
        real(KIND=wr) :: vec10000(10000,2)                                  ! stocke les 10000 meilleurs modèles (valeur du paramètre, valeur de la fonction coût)
    END TYPE densityplot_one
  ! -------                                                     --------    .
  ! -------                                                     --------    ! relecture modèles sélectionnés et calculs a posteriori
    TYPE densityplot
        type(densityplot_one) :: mis,VC,VM,Zmoho,VpVs                       ! ensemble des paramètres
        type(densityplot_one) :: Lat(nbseismes),Lon(nbseismes)
        type(densityplot_one) :: Zhypo(nbseismes),Tzero(nbseismes)
        type(date_sec) :: temps_ref(nbseismes)                              ! temps zéro en minutes du séisme
        integer(KIND=wi) :: nbparam                                         ! nombre de modéles sélectionnées au cours de l'inversion McMC
        integer(KIND=wi) :: deltaxy                                         ! nombre de disrcétisation pour le diagramme de densité, les histogrammes et le calcul du mode
    END TYPE densityplot
  ! -------                                                     --------    .
  ! -------                                                     --------    ! fonction coût au cours de l'inversion
    TYPE fcout
        real(KIND=wr) :: old,new,best                                       ! à l'itération précedente, actuelle, meilleur pour la chaîne
    END TYPE fcout
  ! -------                                                     --------    .
  ! -------                                                     --------    ! acceptance coût au cours de l'inversion
    TYPE accept
        integer(KIND=wl) :: N,O,NO                                          ! nombre de modèles non acceptés, acceptés et repêchés
        real(KIND=wr) :: val                                                ! valeurs en pourcentage
    END TYPE accept
  ! -------                                                     --------    .
  ! -------                                                     --------    ! événements sismiques du catalogue
    TYPE seismes
        type(date_sec) :: tps_init                                          ! date
        real(KIND=wr) :: mag                                                ! mL
        real(KIND=wr) :: lon, lat, pfd                                      ! coordonnées (degrés), profondeur hypocentre (km)
        real(KIND=wr) :: d_t, d_epi, d_p                                    ! différences en temps (s), distance épicentrale (km) et profondeur (km) entre le catalogue et le meilleur modèle rencontré lors de l'inversion
        character(LEN=20) :: name                                           ! nom du bulletin
    END TYPE seismes
  ! -------                                                     --------    !
  ! -------                                                     --------    ! pour des statistiques sur les coldruns (T : toutes les chaines ; S : les chaines sélectionnées)
    TYPE coldmoyval
        real(KIND=wr), dimension(:), allocatable :: Tmis,TVC,TVM,TZmoho,TVpVs
        real(KIND=wr), dimension(:,:), allocatable :: TLat,TLon,TZhypo,TTzero
        real(KIND=wr), dimension(:), allocatable :: Smis,SVC,SVM,SZmoho,SVpVs
        real(KIND=wr), dimension(:,:), allocatable :: SLat,SLon,SZhypo,STzero
    END TYPE coldmoyval
  ! -------                                                     --------    !
  ! -------                                                     --------    !
    TYPE coldmoy
        type(date_sec) :: tempsrefcold(nbseismes)
        type(paramisfit) :: moytot,ectot,moyselect,ecselect                 ! statistiques sur les coldruns
    END TYPE coldmoy
  ! -------                                                     --------    !
    TYPE ellip
        real(KIND=wr) :: ang, axeA, axeB                                    ! ellipse (autour de l'épicentre) +- 1 sigma des 1000 meilleurs modeles
    END TYPE ellip
  ! -------                                                     --------    !
    TYPE residus
        character(LEN=4) :: staname
        real(KIND=wr), dimension(:,:), allocatable :: resPg,resSg,resPn,resSn ! résidus à la station, onde P et S
        integer(KIND=wi) :: nbPg,nbSg,nbPn,nbSn,nbPgT,nbSgT,nbPnT,nbSnT
    END TYPE residus
  ! -------                                                     --------    !
    TYPE apriorEPI
        real(KIND=wr) lat,lon,distcarre                                     ! rechercher epicentre initial
    END TYPE apriorEPI
  ! -------                                                     --------    !
    TYPE priorEPI
        integer(KIND=wi) nb                                                 !  nb de cellules (mailles)
        type(apriorEPI), dimension(:), allocatable :: pEpi
    END TYPE priorEPI
    ! -------                                                    --------    ! amoho_centroid
    TYPE amoho_centroid
        real(KIND=wr) :: lonC,latC                                           ! centroïde
        real(KIND=wr) :: alph,beta,gamma                                     ! vecteur normal du moho (lon, lat, z)
        real(KIND=wr) :: NS,EO                                               ! angle apparant du moho (lon,lat) au centroïde (0 degrés : horizontal ; ~ 90 degrés : vertical)
    END TYPE amoho_centroid
  ! ---------------------------------------------------------------------   .


CONTAINS

  ! ---------------------------------------------------------------------   .

  subroutine mvPall_2_P1(p1,p2,j)
    ! -------                                                   --------    .mh
    ! parametres pour tous les séismes -> pour le jieme séisme
    ! -------                                                   --------    .
    implicit none
    ! -------                                                   --------    .
    integer(KIND=wi), intent (in) :: j
    type(parametres), intent (in) :: p2
    type(parametre), intent (inout) :: p1
    ! ------------------------------------------------------------------    .
    p1%VC=p2%VC
    p1%VM=p2%VM
    p1%Zmoho=p2%Zmoho
    p1%VpVs=p2%VpVs
    p1%Lat=p2%Lat(j)
    p1%Lon=p2%Lon(j)
    p1%Zhypo=p2%Zhypo(j)
    p1%Tzero=p2%Tzero(j)
    ! ------------------------------------------------------------------    .
  end subroutine mvPall_2_P1

  ! ---------------------------------------------------------------------   .

  subroutine mvP1_2_Pall(p2,p1,j)
    ! -------                                                   --------    .mh
    ! parametre pour le jieme séisme -> pour tous les séismes
    ! -------                                                   --------    .
    implicit none
    ! -------                                                   --------    .
    integer(KIND=wi), intent (in) :: j
    type(parametres), intent (inout) :: p2
    type(parametre), intent (in) :: p1
    ! ------------------------------------------------------------------    .
    p2%VC=p1%VC
    p2%VM=p1%VM
    p2%Zmoho=p1%Zmoho
    p2%VpVs=p1%VpVs
    p2%Lat(j)=p1%Lat
    p2%Lon(j)=p1%Lon
    p2%Zhypo(j)=p1%Zhypo
    p2%Tzero(j)=p1%Tzero
    ! ------------------------------------------------------------------    .
  end subroutine mvP1_2_Pall

  ! ---------------------------------------------------------------------   .

  subroutine alph2vect(acentroid)
    ! -------                                                   --------    .mh
    ! angle apparant du moho -> vecteur normal, moho non tabulaire
    ! NS et OE : angle entre 0 (horizontal) et ~90 (vertical), penche positivement vers l'Est et le Sud
    ! -------                                                   --------    .
    implicit none
    type (amoho_centroid), intent (inout) :: acentroid
    ! ------------------------------------------------------------------    .
    acentroid%alph=tan(acentroid%EO*pi/180._wr)
    acentroid%beta=tan(acentroid%NS*pi/180._wr)
    acentroid%gamma=-1.0_wr
    ! ------------------------------------------------------------------    .
  end subroutine alph2vect

  ! --------------------------------------------------------------------    .

  subroutine vect2alph(acentroid)
    ! -------                                                   --------    .mh
    ! vecteur normal -> angle apparant du moho, moho non tabulaire
    ! NS et OE : angle entre 0 (horizontal) et ~90 (vertical), penche positivement vers l'Est et le Sud
    ! -------                                                    --------    .
    implicit none
    type (amoho_centroid), intent (inout) :: acentroid
    ! ------------------------------------------------------------------    .
    acentroid%NS=atan(acentroid%beta/acentroid%gamma)/pi*180._wr
    acentroid%EO=atan(acentroid%alph/acentroid%gamma)/pi*180._wr
    ! ------------------------------------------------------------------    .
  end subroutine vect2alph

END MODULE typetemps



! **********************************************************************    .
! **********************************************************************    .


