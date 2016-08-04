! programme principal II bis
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
!                                                                            !
!             L'expérience cubaine, puis au Congo ou en Bolivie :            !
!                     différentes réalisations concrètes                     !
!                 mettant à jour le visage de la Révolution.                 !
!                                                                            !
! -----------------------------------------------------------------------    .
! Période des réalisations de densités a posteriori                          !
!                                                                            !
program che2013_hotruns
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use time
  use mt19937
  use cpt_temps
  use algo_metropolis
  use tirage
  use datalecture
  use pb_direct
  use affiche
  use sub_param
  use sub_misfit
  use distance_epi
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none

    include 'mpif.h'

  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(parametres) :: param_init,param_best                                  ! paramètres d'inv.
  type(fcout) :: misfit                                                      ! fonction coût
  type(accept) :: acceptance                                                 ! acceptance
  type(parametresinv) :: p                                                   ! paramètres d'inv.
  type(coldmoy) :: dc                                                        ! modèles du coldrun
  type(residus), dimension(:), allocatable :: R                              ! résidus (si FLAGresSTA=.true.)
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  ! -------                                                      --------    .
  real(KIND=wr) :: xmin(nbseismes), xmax(nbseismes)                          ! cercles pond.
  ! real(KIND=wr) :: VPVSch                                                    ! ratio VpVs défini par le diagramme de Châpelain
  ! -------                                                      --------    .
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: i,j,k,l,ok,pourcentage
  integer(KIND=wi) :: noctet, nbauto
  integer(KIND=wi) :: nbChaineMVhot,nbChaineMVcold                           ! nombre chaînes
  integer(KIND=wi) :: maxiterhot,maxitercold                                 ! nombre d'itérations
  integer(KIND=wi) :: nbsta,nbstaR,nbtps(nbseismes)                          ! nombre station possible, de phases et reel nombre de données de temps
  integer(KIND=wi) :: nmod                                                   ! nombre modèles sélectionnés par chaine
  ! -------                                                      --------    .
  logical :: critique                                                        ! si moho trop bas, onde refacté observé mais non prédite (vrai)
  logical :: savemod                                                         ! modèle sauvé (vrai)
  logical :: initauto                                                        ! modèle sauvé (vrai)
  logical :: accepte                                                         ! modèle accepté (vrai) ou rejeté (faux)
  ! -------                                                      --------    .
  character (LEN=5) :: numberchaine
  character (LEN=30) :: chaine
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
    call print_mess_2bis
    call initseed(libre)                                                     ! aléatoire calé sur temps CPU
    ! ------- clacul de nb_procs graines                         --------    .

    allocate(allseed(nb_procs))
    do i=1,nb_procs
      allseed(i)=abs(genrand_int31())+10000
      do while (allseed(i) > 10000000)
        seed = int(10000123.0_wr*genrand_real1())
        allseed(i) = allseed(i) - seed
      enddo
    enddo
  endif
  ! ------- partage des graines                                  --------    .
  call MPI_SCATTER(allseed,1,MPI_INTEGER,seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
  call init_genrand (int(seed*(rang+1),wi))                                  ! initialisation des graines pour chaque processus
  !write(*,*)'générateur de nombre aléatoire :       ',seed*(int(rang+1))
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  ! ------- lecture des données                                  --------    .
  call lectnbdata(nbsta,nbtps)
  ! ---------------------------------------------------------------------    .
  ! -------  nb de chaines de Markov et d'iterations par chaine  --------    .
  call lectparam(nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold,chut=.true.)
  ! ---------------------------------------------------------------------    . lecture des coldruns
  ok = 0
  write(numberchaine(1:5),'(i5)')rang+1
  open(unit=100,file="OUTPUT/files/Hot/In_"//trim(adjustl(numberchaine))//".bin", &
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
  read(100)dc
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
  if (ok.ne.0) then
    write(*,*)"problème dans che_hotruns : le fichier OUTPUT/files/Hot/In_"//trim(adjustl(numberchaine))//".bin n''existe pas "
    call MPI_ABORT(MPI_COMM_WORLD,err,code)
  endif
  p%valOLd = param_init
  p%valNew = param_init
  ! ---------------------------------------------------------------------    . calul résidus aux stations
  if (FLAGresSTA) call initR(D,R,maxiterhot,nbtps,nbsta,nbstaR)
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .
  ! ------- début MCMc                                           --------    . HOTRUNS
  ! ---------------------------------------------------------------------    .
  if((rang==0).and.(plotmisfit)) open(unit=10,file="OUTPUT/files/Hot/misH_1.txt",status="replace")
  if((rang==0).and.(plotmisfit)) open(unit=11,file="OUTPUT/files/Hot/misH_2.txt",status="replace")
  ! ---------------------------------------------------------------------    .
  i=rang+1
  ! ---------------------------------------------------------------------    .
    if (rang==0) call print_messchaine(i,nbChaineMVhot)
    ! ------- initialisation                                     --------    .
    nmod = 0                                                                 ! nombre de modèle sélectionnés pdt une chaîne
    ! -------------------------------------------------------------------    .
    ! ------- création d'un fichier par chaîne                   --------    .
    write(numberchaine(1:5),'(i5)')i
    inquire ( iolength = noctet ) misfit%new,p%valNew
    open(unit=205+i,file="OUTPUT/files/"//trim(adjustl(numberchaine))//".bin",status="replace",access='direct',RECL=(noctet))
    initauto=.false.                                                         ! les premiers modèles ne sont pas sauvés
    savemod=.false.
    ! -------------------------------------------------------------------    . calul du ratio VpVs
    !call Chatelainplot(nbtps,D,vpvs=VPVSch)
    !write(*,'(a,i3.3)')' VP/VS, Chatelain : 1,',int(VPVSch*1000.0_wr,wi)-1000
    !call Wadatiplot(nbtps,D,p%valNew,vpvs=VPVSch)
    !write(*,'(a,i3.3)')' VP/VS, Wadati : 1,',int(VPVSch*1000.0_wr,wi)-1000
    ! -------------------------------------------------------------------    .
    ! ------- début d'une chaîne                                 --------    .
    ! -------------------------------------------------------------------    .
    unechaine  : do k=1,maxiterhot
      ! -----------------------------------------------------------------    . progression en %
      if (i==1) then
        if ((k.gt.1000).or.(k.ne.1).or.(k.ne.maxiterhot)) then
          if(mod(k,100)==0) write(chaine(1:30),'(a5,i12,a13)')" 2 - ",k," modèles   "
        else
          write(chaine(1:30),'(a5,i12,a13)')" 2 - ",k," modèles     "
          call progress(k,maxiterhot,chaine)
        endif
        call progress(k,maxiterhot,chaine)
      endif
      ! -----------------------------------------------------------------    .
      j=0
      critique=.true.
      do while(critique)                                                     ! tant que distance hypocentral << distance critique pour la réfraction
        j=j+1

        ! ------- tirage au sort des paramètres dans le prior      --------  .
        if(k.ne.1) then
          pourcentage=int(genrand_real1()*100._wr)                           ! aléatoire de 0 à 99
          l=int(genrand_real1()*real(nbseismes,wr)+1.0_wr)                   ! aléatoire de 1 à nbseismes
          select case (pourcentage)
            case (0:24)
              call tirage_H(p,l,all=.true.)                                  ! 25 % : tirage gaussien de de tous les paramètres hypocentraux pour le 'l'ième séisme
            case (25:49)
              call tirage_H(p,l,all=.false.)                                 ! 25 % : tirage gaussien de 1 des 4 paramètres hypocentraux pour le 'l'ième séisme
            case (50:74)
              call tirage_T(p,all=.true.,vpvs=.true.)                        ! 25 % : tirage gaussien de de tous les paramètres de terre (dont VpVs)
            case (75:99)
              call tirage_T(p,all=.false.,vpvs=.true.)                       ! 25 % : tirage gaussien de 1 des 4 paramètres de terre (dont VpVs)
          end select
        endif

        ! -----------------------------------------------------------------  .
        ! ------- problème direct pour le jeu de paramètre tiré et chaque donnée
        call tempsTheoDirect(nbtps,p%valNew,D,critique,acentroid)
        if(j==5) then
          critique=.false.                                                   ! pour sortir apres 5 essais
          ! ------- une donnée réfractée observée ne peux pas etre inférieur à la distance hypocentrale minimale pour la refraction
          write(*,*)'problème dans che_hotruns : distance épi + 5 km < distance épi critique pour la réfraction '
        endif
      enddo
      ! ------- calcul de la fonction coût                       --------    .
      call compute_misfit(nbtps,D,misfit%new,xmin,xmax,'H')
      ! ------- modele sauvé ?                                   --------    .
      if (.not.initauto) then
        if (k.gt.(maxiterhot/10)) initauto=.true.                            ! on ne garde pas les 10 premier %
        nbauto=(maxiterhot/10)+max(10,int(normal(real(8*nbseismes,wr),real(nbseismes,wr)))) ! procchain modèle sauvé
      else
        if(k==nbauto) then
          savemod=.true.                                                     ! le modele est sauvé
          nmod=nmod+1
          nbauto=k+max(10,int(normal(real(8*nbseismes,wr),real(nbseismes,wr)))) ! procchain modèle sauvé
        else
          savemod=.false.                                                    ! le modele n'est pas sauvé
        endif
      endif
      ! ------- Metropolis (acceptation et rejet des modèles)    --------    .
      call metropolis(p,param_best,misfit,acceptance,savemod,accepte)
      ! -------                                                  --------    .
      if((rang==0).and.(plotmisfit)) then
        if(initauto) then
          if (savemod) write(11,*)k+maxitercold,misfit%new
        else
         write(10,*)k+maxitercold,misfit%new
        endif
      endif
      ! ------- écriture dans les fichiers                       --------    .
      if(savemod) then
        write(205+i,rec=nmod)misfit%new,p%valNew                             ! écriture du modèle sélectioné dans un fichier files/k.bin
        if (FLAGresSTA) call inR(D,R,nbtps,nbstaR,misfit%new)                ! sauve les résidus dans R -> evaluation des retard aux stations
      endif
      ! -----------------------------------------------------------------    .
    enddo unechaine
    call calc_accept(acceptance)

    ! -------------------------------------------------------------------    .
    ! ------- fin d'une chaîne                                   --------    .
    ! -------------------------------------------------------------------    .
    if (FLAGresSTA) call outR(R,nbstaR)                                      ! sauve les résidus dans des fichiers .bin
    ! -------------------------------------------------------------------    .
    close(205+i)
    ! -------------------------------------------------------------------    .
    !call MPI_Barrier(MPI_COMM_WORLD,err)
    !call print_mess_finchainemin(misfit%best,acceptance%val)

  if((rang==0).and.(plotmisfit)) close(10)
  if((rang==0).and.(plotmisfit)) close(11)
  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .
  ! ------- fin MCMc                                             --------    . HOTRUNS
  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .

  ! ----------------  un fichier OUTPUT par cold runs  ------------------    .
  open(unit=100,file="OUTPUT/files/Hot/Out_"//trim(adjustl(numberchaine))//".bin", &
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
  write(100)dc
  write(100)nmod
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
    deallocate (D(i)%datatps,pEpis(i)%pEpi)
  enddo
  if(allocated(R)) deallocate(R)
  ! ---------------------------------------------------------------------    . MPI_END
  call MPI_FINALIZE(code)
  ! ---------------------------------------------------------------------    . MPI_END

end program che2013_hotruns



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .

