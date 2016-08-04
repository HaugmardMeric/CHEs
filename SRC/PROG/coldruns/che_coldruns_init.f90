! programme principal I
! ***********************************************************************    .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! -------                                                        --------    .
! ------- CHE2013_coldruns version 2.2                           --------    .
! ------- octobre 2013 - décembre 2014                           --------    .
! -------                                                        --------    .
! ------- Prog. basé uniquement sur des méthodes non linéaires   --------    .
! -------                                                        --------    .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr           --------    .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
!  This program is distributed for research purposes and in the hope         !
!  that it will be useful however without any warranties.                    !
!  Use it on your own risk. Please, report bugs/improvements/... .           !
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! avec la participation de Ianis Gaudot, Éric Beucler et Philippe Cance
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
!                                                                            !
!         " Une réponse approximative à la bonne question, qui est           !
!            souvent mal posée, est bien meilleure que la réponse            !
!           exacte à une mauvaise question que l’on peut toujours            !
!                        formuler de façon précise … "                       !
!                                                                            !
!       Tukey, J.W. (1962) : The Future of Data Analysis.                    !
!       The Annals of Mathematical Statistics, Vol. 33, No. 1, p. 1-67.      !
!                                                                            !
!                                                                            !
! -----------------------------------------------------------------------    .
! version 1.1 : inversion d'un séisme (jan 2014)                             !
! version 1.2 : inversion de plusieurs séismes (juillet 2014)                !
! version 1.3 : parallelisation OpenMP (septembre 2014)                      !
! version 1.4 : parallelisation MPI (octobre 2014)                           !
! version 1.5 : initialisation du prior auto (novembre 2014)                 !
! version 1.6 : parallelisation full MPI (décembre 2014)                     !
! version 1.7 : ajout d'une notice (décembre 2014)                           !
! version 1.8 : compilation avec ifort et gfortran (janvier 2015)            !
! version 1.9 : moho non tabulaire (fevrier 2015)                            !
! version 2.0 : test 50 séismes (septembre 2015)                             !
! version 2.1 : gestion des carrières, d'après Pascal Guterman (octobre 2015)!
! version 2.2 : ajout de modèle de terre différents pour le problème directe (novembre 2015)!
! version 2.3 : calcules a posteriori (janvier 2016)                         !
! -----------------------------------------------------------------------    .
! Initialisation des Coldruns                                                !
!                                                                            !
program che2013_coldruns_init
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use mt19937
  use datalecture
  use affiche
  use sub_param
  use rechercheepi
  use time
  use figure_posteriori
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none
  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(parametres), dimension(:), allocatable :: param_init                  ! paramètres d'inv.
  type(fcout) :: misfit                                                      ! fonction coût
  type(accept) :: acceptance                                                 ! acceptance
  type(parametresinv) :: p                                                   ! paramètres d'inv.
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  type(date_secPS) :: midi20
  ! -------                                                      --------    .
  real(KIND=wr) :: xmin(nbseismes), xmax(nbseismes)                          ! cercles pond.
  real(KIND=wr) :: deltaP,deltaS
  ! -------                                                      --------    .
  integer(KIND=wi) :: i,j,k
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: nbChaineMVhot,nbChaineMVcold                           ! nombre chaînes
  integer(KIND=wi) :: maxiterhot,maxitercold                                 ! nombre d'itérations
  integer(KIND=wi) :: nbsta, nbtps(nbseismes)                                ! nombre station et nombre de données de temps par séismes
  integer(KIND=wi) :: nbmod
  ! -------                                                      --------    .
  character (LEN=5) :: numchaine
  character (LEN=20) :: nomfichier

  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  ! ------- initialisation                                       --------    .
  call print_mess_1
  call initseed(libre)                                                       ! aléatoire calé sur temps CPU
  call printnbseismes
  ! ------- si moho non tabulaire                                --------    .
  acentroid%lonC=moho_lon ; acentroid%latC=moho_lat
  acentroid%NS=moho_NS ; acentroid%EO=moho_EO
  call alph2vect(acentroid); call vect2alph(acentroid)
  ! ------- lecture des données                                  --------    . phases and stations list
  call lectnbdata(nbsta,nbtps)                                               ! nombre de données
  do i=1,nbseismes
    allocate(D(i)%datatps(nbtps(i)))                                         ! alloue par séisme
  enddo
  call lectdata(nbsta,nbtps,D)                                               ! temps d'arrivés par séisme
  ! ------- cas synthétiques (si besoin)                         --------    .
  call mksynth(nbtps,D,acentroid)
  ! ------- réduction du prior pour les paramètres épicentraux   --------    .
  print*,'initialisation du prior pour les paramètres épicentraux'
  call zoneRecherche(nbtps,D,pEpis,mb)
  ! ---------------------------------------------------------------------    .
  ! -------  nb de chaines de Markov et d'iterations par chaine  --------    .
  call lectparam(nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold)
  allocate(param_init(nbChaineMVcold))
  ! ---------------------------------------------------------------------    .

  do i=1,nbChaineMVcold
    ! -------------------------------------------------------------------    .
    ! ------- initialisation du modèle de terre puis des paramètres hypocentraux
    call initparam(nbtps,D,param_init(i),pEpis,mb)                           ! par la méthode des hémisphères -> méthode non linéaire
    ! ------- lecture du prior                                   --------    .
    call lect_prior(p,param_init(i),"C")                                     ! C -> coldruns
    ! ------- initialisation des bornes épicentales              --------    .
    p%centreY(:)=param_init(i)%lat(:)
    p%centreX(:)=param_init(i)%lon(:)
    p%Rayon=500.0_wr                                                         ! recherche de l'épicentre dans ce rayon, pas au delas
    ! -------------------------------------------------------------------    .
  enddo
  ! ------- recherche les distances de pondération               --------    .
  call cerclespond(nbtps,D,param_init(1),xmin,xmax,acentroid)
  ! ------- initialisation des autres paramètres                 --------    .
  call init_div(misfit,acceptance)
  ! ---------------------------------------------------------------------    . si certains parametres fixes
  call paramfixe(p)
  ! ---------------------------------------------------------------------    .


  ! ------------------------------------------------------------------------ .
  ! Carriere
  ! ------------------------------------------------------------------------ .
  !do i=1,nbseismes
  !  midi20=D(i)%datatps(1)%tpsR
  !  midi20%date%hour=11                                                      ! midi TU, pour Vannes en hiver
  !  midi20%date%min=20
  !  call basetime(midi20)
  !  call difftime(deltaP,deltaS,midi20,D(i)%datatps(1)%tpsR)
  !  if ((abs(deltaP).lt.1.00_wr).or.(i==45)) then              ! une heure
  !    write(*,*)'seisme ',i,' : tire de carriere ?'
  !    write(*,*)'--- > Zhypo max = 6 km'
  !    p%maxi%Zhypo(i) = 6.d0
  !    p%ecartype%Zhypo(i) =  p%ecartype%Zhypo(i)/4.0_wr
  !  endif
  !  ! ---------------------------------------------------------------------- .
  !enddo
  ! ------------------------------------------------------------------------- .


  ! ------------------------------------------------------------------------- .
  ! étude des gradients sur la fonction coût
  ! ------------------------------------------------------------------------- .
  if (plotposteriori) then
    nbmod=1
    nomfichier='POST_COLD_i'
    !call PosterioriExploration(p,nbmod,pEpis,nbtps,D,acentroid,xmin,xmax,mb,nomfichier)
  endif
  ! ------------------------------------------------------------------------- .


  ! production de fichier
  ! parce que dans la norme MPI-Fortran 2008,
  ! faire un call MPI_BCAST avec des types dérivés, c'est pas prévu !
  ! --------------------  un fichier par cold runs  ---------------------    .
  do k=1,nbChaineMVcold
    write(numchaine(1:5),'(i5)')k
    open(unit=100,file="OUTPUT/files/Cold/In_"//trim(adjustl(numchaine))//".bin", &
        status="replace",form="unformatted",access="sequential")
      write(100)nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold
      write(100)nbsta,nbtps
      write(100)p
      do i=1,nbseismes
        do j=1,nbtps(i)
          write(100)D(i)%datatps(j)
        enddo
      enddo
      write(100)misfit
      write(100)acceptance
      write(100)param_init(k)
      write(100)xmin,xmax
      do i=1,nbseismes
        write(100)pEpis(i)%nb
        do j=1,pEpis(i)%nb
          write(100)pEpis(i)%pEpi(j)
        enddo
      enddo
      write(100)mb
      write(100)acentroid
    close(100)
  enddo
  ! ---------------------  fin du programme -----------------------------    .
  deallocate(param_init)
  do i=1,nbseismes
    deallocate(D(i)%datatps,pEpis(i)%pEpi)
  enddo
  ! ---------------------------------------------------------------------    .

end program che2013_coldruns_init



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .


