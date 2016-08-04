! Librairies pour le programme CHEs
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE modparam

    implicit none

    private

    public :: wi, wl, wr
    public :: nbseismes
    public :: pi, rT, Tminsec
    public :: FLAGterre, FLAGterrefixe
    public :: FLAGhypo, FLAGhypofixe
    public :: FLAGresSTA
    public :: FLAGcercles
    public :: plotgraph,plotposteriori
    public :: libre
    public :: tracessac
    public :: autocorr
    public :: FLAG_non_tabulaire,moho_lon,moho_lat,moho_NS,moho_EO
    public :: printnbseismes

    ! -------   parametres ForTran precision numérique :       --------    .
    integer, parameter :: intg = selected_int_kind(r=9)                  ! |int| < 10^9
    integer, parameter :: long = selected_int_kind(r=18)                 ! |int| < 10^18
    integer, parameter :: dobl = selected_real_kind(p=15,r=307)          ! |int| < 10^307, 15 chiffres significatifs au moins
    integer, parameter :: wi = intg
    integer, parameter :: wl = long
    integer, parameter :: wr = dobl
    ! -----------------------------------------------------------------    .

    ! *****************************************************************    .
    ! -----------------------------------------------------------------    .
    ! nombre de séismes
    integer(KIND=wi), parameter :: nbseismes=6
    ! -----------------------------------------------------------------    .
    ! *****************************************************************    .
    ! -----------------------------------------------------------------    .
    ! -------                                                  --------    . paramètres de terre
    ! tirage aléatoire selon une loi normale dont la moyenne est fixe :
    logical, parameter :: FLAGterrefixe=.false.
    ! prior resserré : (lu dans PARAM/paramTerre.d)
    logical, parameter :: FLAGterre=.false.
    ! -------                                                  --------    . paramètres hypocentraux
    ! tirage aléatoire selon une loi normale dont la moyenne est fixe :
    logical, parameter :: FLAGhypofixe=.false.
    ! prior resserré : (lu dans PARAM/paramHypo.d)
    logical, parameter :: FLAGhypo=.false.
    ! -------                                                  --------    . résidus
    ! calcul résidus aux stations (pour chaque modèle sauvé en hotruns,
    ! les différence temps théoriques-réels sont stockés par types d'onde et
    ! par stations afin de voir un décalage, ou résidus à la station)
    logical, parameter :: FLAGresSTA=.true.
    ! -------                                                  --------    . nombre aléatoire
    ! générateur aléatoire fixé (libre=vrai)
    logical, parameter :: libre=.true.
    ! -------                                                  --------    .
    ! plot les traces sac sur hodochrones et calcul de la magnitude Ml
    logical, parameter :: tracessac=.true.
    ! -------                                                  --------    .
    ! prise en compte des cercles de pondération xmin et xmax
    ! attention : parfois l'écipentre se déplace et aucune station se rouve dans le cercle 
    logical, parameter :: FLAGcercles=.true.
    ! -------                                                  --------    .
    ! plot les graphs (chaine + autocorr + ... ) pour chaque param
    logical, parameter :: plotgraph=.true.
    ! -------                                                  --------    .
    ! plot les graphs (posteriori )
    logical, parameter :: plotposteriori=.true.
    ! -------                                                  --------    .
    ! écart minimal de temps entre deux stations pour la recherche initiale
    real(KIND=wr), parameter :: Tminsec=0.50_wr
    ! -------                                                  --------    .
    ! vrai si moho penche ; faux sinon
    logical, parameter :: FLAG_non_tabulaire=.false.
    ! centre ou est estimée la profondeur du moho
    real(KIND=wr), parameter :: moho_lon=-2.499999999_wr
    real(KIND=wr), parameter :: moho_lat=47.499999999_wr
    ! angle entre 0 (horizontal) et ~90 (vertical), le moho penche positivement vers l'Est et le Sud
    real(KIND=wr), parameter :: moho_NS=0.0_wr
    real(KIND=wr), parameter :: moho_EO=0.0_wr

    ! -------                                                  --------    .
    integer(KIND=wi), parameter :: autocorr=10000                          ! nombre du modèles maximal pour le calcul de l'autocorrelation

    ! -----------------------------------------------------------------    .
    ! *****************************************************************    .
    ! -----------------------------------------------------------------    .

    real(KIND=wr), parameter :: pi = 3.141592653589793238_wr
    real(KIND=wr), parameter :: rT = 6371.0_wr                             ! rayon terrestre moyen

    ! -----------------------------------------------------------------    .
    ! *****************************************************************    .
    ! -----------------------------------------------------------------    .

CONTAINS

  ! ---------------------------------------------------------------------   .

  subroutine printnbseismes
    ! -------                                                   --------    .mh
    ! écriture du nombre de séismes
    ! -------                                                   --------    .
    implicit none
    integer(KIND=wi) :: ok
    ! -------                                                   --------    .
    open(1999, FILE = 'OUTPUT/GMT/nbseisme.d',status='replace',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans printnbseismes : le fichier OUTPUT/GMT/nbseisme.d n''existe pas '
      stop
    endif
    ! -------                                                   --------    .
    write(1999,*)nbseismes
    ! -------                                                   --------    .
    close(1999)
    ! ------------------------------------------------------------------    .
  end subroutine printnbseismes

END MODULE modparam



! *********************************************************************    .
! *********************************************************************    .


