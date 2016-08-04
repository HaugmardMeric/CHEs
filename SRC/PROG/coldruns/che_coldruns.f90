! programme principal I bis
! ***********************************************************************    .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! -------                                                        --------    .
! ------- CHE2013_coldruns version 1.5                           --------    .
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
!      Les voyages initiatiques de Ernesto CHE Gevara de 1951 à 1954 :       !
!       recherche à grande longeur d'onde d'un esprit révolutionnaire.       !
!                                                                            !
! -----------------------------------------------------------------------    .
! Période de rodage                                                          !
!                                                                            !
program che2013_coldruns
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use mt19937
  use cpt_temps
  use algo_metropolis
  use tirage
  use pb_direct
  use affiche
  use sub_param
  use sub_misfit
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none

    include 'mpif.h'

  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(parametres) :: param_init, param_best                                 ! paramètres d'inv.
  type(fcout) :: misfit                                                      ! fonction coût
  type(accept) :: acceptance                                                 ! acceptance
  type(parametresinv) :: p                                                   ! paramètres d'inv.
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  ! -------                                                      --------    .
  real(KIND=wr) :: xmin(nbseismes), xmax(nbseismes)                          ! cercles pond.
  ! -------                                                      --------    .
  integer(KIND=wi) :: i,j,k,l,ok
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: nbChaineMVhot,nbChaineMVcold                           ! nombre chaînes
  integer(KIND=wi) :: maxiterhot,maxitercold                                 ! nombre d'itérations
  integer(KIND=wi) :: nbsta, nbtps(nbseismes)                                ! nombre station et nombre de données de temps par séismes
  ! -------                                                      --------    .
  logical :: critique                                                        ! si moho trop bas, onde refacté observé mais non prédite (vrai)
  logical :: savemod                                                         ! modèle sauvé (vrai)
  logical :: accepte                                                         ! modèle accepté (vrai) ou rejeté (faux)
  logical :: div                                                             ! chaîne divergente
  ! -------                                                      --------    .
  character (LEN=30) :: chaine
  character (LEN=5) :: numbchaine
  ! -------                                                      --------    . MPI :
  integer :: nb_procs,rang,code,err
  integer, dimension(:), allocatable :: allseed
  integer :: seed
  ! -------                                                      --------    .
  logical, parameter :: plotmisfit=.true.
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    . MPI_BEGIN
  call MPI_INIT(code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  ! ---------------------------------------------------------------------    . MPI_BEGIN
  ! ------- initialisation de la graine                          --------    .
  if (rang==0) then
    call print_mess_2
    call initseed(libre)                                                     ! aléatoire calé sur temps CPU
    ! ------- clacul de nb_procs graines                         --------    .
    allocate(allseed(nb_procs))
    do i=1,nb_procs
      allseed(i)=abs(genrand_int31())+10000
      do while (allseed(i) > 10000000)
        seed = int(10000123.0_wr*genrand_real1())
        allseed(i) = allseed(i) - seed
      end do
    enddo
  endif
  ! ------- partage des graines                                  --------    .
  call MPI_SCATTER(allseed,1,MPI_INTEGER,seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
  call init_genrand (int(seed*(rang+1),wi))                                  ! initialisation des graines pour chaque processus
  ! write(*,*)'générateur de nombre aléatoire :       ',seed*(int(rang+1))
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  ! ------- lecture paramètres et données                        --------    .
  write(numbchaine(1:5),'(i5)')rang+1
  ok=0
  open(unit=100,file="OUTPUT/files/Cold/In_"//trim(adjustl(numbchaine))//".bin", &
      status="old",form="unformatted",access="sequential",iostat = ok)
  read(100)nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold
  read(100)nbsta,nbtps
  read(100)p
  do i=1,nbseismes
    allocate(D(i)%datatps(nbtps(i)))
    do j=1,nbtps(i)
      read(100)D(i)%datatps(j)
    enddo
  enddo
  read(100)misfit
  read(100)acceptance
  read(100)param_init
  read(100)xmin,xmax
  do i=1,nbseismes
    read(100)pEpis(i)%nb
    allocate(pEpis(i)%pEpi(pEpis(i)%nb))
    do j=1,pEpis(i)%nb
      read(100)pEpis(i)%pEpi(j)
    enddo
  enddo
  read(100)mb
  read(100)acentroid
  close(100)
  ! -------                                                      --------    .
  if (ok .ne. 0) then
    write(*,*)"problème dans che2013_coldruns : le fichier ", &
      "OUTPUT/files/Cold/In_"//trim(adjustl(numbchaine))//".bin n''existe pas "
    call MPI_ABORT(MPI_COMM_WORLD,err,code)
  endif

  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .
  ! ------- début MCMc                                           --------    . COLDRUNS
  ! ---------------------------------------------------------------------    .
  if(plotmisfit) open(unit=10,file="OUTPUT/files/Cold/mis_"//trim(adjustl(numbchaine))//".txt",status="replace")
  ! ---------------------------------------------------------------------    .
  i=rang+1
  ! ---------------------------------------------------------------------    .
  if (rang==0) call print_messchaine(i,nbChaineMVcold)

    ! -------------------------------------------------------------------    .
    ! ------- début d'une chaîne                                 --------    .
    ! -------------------------------------------------------------------    .
    savemod=.true.
    div=.false.
    unechaine : do j=1,maxitercold
      ! -----------------------------------------------------------------    . progression en %
      if (i==1) then
        if ((j.gt.1000).or.(j.ne.1).or.(j.ne.maxitercold)) then
          if(mod(j,100)==0) write(chaine(1:30),'(a5,i12,a13)')" 1 - ",j," modèles   "
        else
          write(chaine(1:30),'(a5,i12,a13)')" 1 - ",j," modèles     "
        endif
        call progress(j,maxitercold,chaine)
      endif
      ! -----------------------------------------------------------------    .
      k=0
      critique=.true.
      do while(critique)                                                     ! tant que distance hypocentral << distance critique pour la réfraction
        k=k+1

          l=int(genrand_real1()*real(nbseismes,wr)+1.0_wr)                   ! aléatoire de 1 à nbseismes, un seul séisme
          call tirage_H(p,l,all=.true.)                                      ! tirage des tous les paramètres hypocentraux pour le 'l'ième séisme
          call tirage_T(p,all=.true.,vpvs=.true.)                            ! tirage des tous les paramètres de terre simultanément

        ! ------- problème direct pour le jeu de paramètre tiré et chaque donnée
        call tempsTheoDirect(nbtps,p%valNew,D,critique,acentroid)
        if(k==5) then
          critique=.false.                                                   ! pour sortir apres 5 essais
          ! ------- une donnée réfractée observée ne peux pas etre inférieur à la distance hypocentrale minimale pour la refraction
          write(*,*)'problème dans che_coldruns : distance épi + 5 km < distance hypocentrale critique pour la réfraction '
        endif
      enddo
      ! ------- calcul de la fonction coût                       --------    .
      call compute_misfit(nbtps,D,misfit%new,xmin,xmax,'C',div)
      ! -----------------------------------------------------------------    .
      if (div) then                                                          ! fin de chaine si divergent
        if (rang==0) write(chaine(1:30),'(a5,i12,a13)')" 1 - ",j," modèles   "
        if (rang==0) call progress(maxitercold,maxitercold,chaine)
        exit unechaine
      endif
      ! ------- Metropolis (acceptation et rejet des modèles)    --------    .
      call metropolis(p,param_best,misfit,acceptance,savemod,accepte)        !  tout le script est là ...
      ! -----------------------------------------------------------------    .
      if(plotmisfit) write(10,*)j,misfit%new
      ! -----------------------------------------------------------------    .
    enddo unechaine
    ! -------------------------------------------------------------------    .
    ! ------- fin d'une chaîne                                   --------    .
    ! -------------------------------------------------------------------    .
    call calc_accept(acceptance)

  ! call MPI_Barrier(MPI_COMM_WORLD,err)
  ! call print_mess_finchainemin(misfit%best,acceptance%val)

  if(plotmisfit) close(10)
  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .
  ! ------- fin MCMc                                             --------    . COLDRUNS
  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .

  ! ----------------  un fichier OUTPUT par cold runs  ------------------    .
  open(unit=100,file="OUTPUT/files/Cold/Out_"//trim(adjustl(numbchaine))//".bin", &
      status="replace",form="unformatted",access="sequential")
    write(100)nbsta,nbtps
    write(100)p
    do i=1,nbseismes
      do j=1,nbtps(i)
        write(100)D(i)%datatps(j)
      enddo
    enddo
    write(100)misfit
    write(100)acceptance
    write(100)param_init
    write(100)param_best
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

  ! ---------------------  fin du programme -----------------------------    .
  do i=1,nbseismes
    deallocate(D(i)%datatps,pEpis(i)%pEpi)
  enddo
  ! ---------------------------------------------------------------------    . MPI_END
  call MPI_FINALIZE(code)
  ! ---------------------------------------------------------------------    . MPI_END

end program che2013_coldruns



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .


