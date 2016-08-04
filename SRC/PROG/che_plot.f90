! programme principal III
! ***********************************************************************    .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! -------                                                        --------    .
! ------- CHE2013_plots version 1.5                              --------    .
! ------- octobre 2013 - décembre 2014                           --------    .
! -------                                                        --------    .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr           --------    .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
!  This program is distributed for research purposes and in the hope         !
!  that it will be useful however without any warranties.                    !
!  Use it on your own risk. Please, report bugs/improvements/... .           !
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
!                                                                            !
!         Mort de Ernesto CHE Gevara en 1967 (La Higuera, Bolivie),          !
!        il devient alors un symbole de la Révolution bien identifié,        !
!     mise en place d'une analyse a posteriori de la Révolution cubaine.     !
!                                                                            !
! -----------------------------------------------------------------------    .
!                                                                            !
!                    Alea jacta est ; les dés sont jetés.                    !
!                                                                            !
! -----------------------------------------------------------------------    .
! période d'exploitation des résultats                                       !
!                                                                            !
program che2013_plot
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use latexscript
  use datalecture
  use affiche
  use sub_param
  use figure_GMT
  use mt19937
  use figure_posteriori
  use time
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none
  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(parametres), dimension(:), allocatable :: param_best                  ! paramètres d'inv.
  type(fcout), dimension(:), allocatable :: misfit                           ! fonction coût
  type(accept), dimension(:), allocatable :: acceptance                      ! acceptance
  type(parametresinv) :: p                                                   ! paramètres d'inv.
  type(densityplot) :: dp
  type(coldmoy) :: dc                                                        ! modèles du coldrun
  type(ellip) :: E(nbseismes)                                                ! éllipses
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  type(parametres), dimension(:), allocatable :: modelesIN
  ! -------                                                      --------    .
  real(KIND=wr) :: xmin(nbseismes), xmax(nbseismes)                          ! cercles pond.
  ! -------                                                      --------    .
  integer(KIND=wi) :: i,j,k,ok
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: nbChaineMVhot,nbChaineMVcold                           ! nombre chaînes
  integer(KIND=wi) :: maxiterhot,maxitercold                                 ! nombre d'itérations
  integer(KIND=wi) :: nbsta, nbtps(nbseismes)                                ! nombre station et nombre de données de temps
  integer(KIND=wi), dimension(:), allocatable :: nmod                        ! nombre modèles sélectionnés par chaine
  integer(KIND=wi) :: nbmod
  ! -------                                                      --------    .
  character (LEN=20) :: nomfichier
  character (LEN=5) :: nbdeseismes
  character(LEN=4), dimension(:), allocatable :: nomsta
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  call initseed(libre)                                                       ! aléatoire calé sur temps CPU

  ! ---------------------------------------------------------------------    .
  ! -------  relecture de quelques varaiables                    --------    .
  ! ---------------------------------------------------------------------    .
  ! ------- lecture des données                                  --------    .
  call lectnbdata(nbsta,nbtps)
  do i=1,nbseismes
    allocate(D(i)%datatps(nbtps(i)))
  enddo
  ! -------  nb de chaines de Markov et d'iterations par chaine  --------    .
  call lectparam(nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold,chut=.true.)
  ! -------                                                      --------    .
  allocate(param_best(nbChaineMVhot))
  allocate(misfit(nbChaineMVhot),acceptance(nbChaineMVhot))
  allocate(nmod(nbChaineMVhot))
  ! ---------------------------------------------------------------------    .
  ! lecture des hotruns
  ok=0
  open(unit=400,file="OUTPUT/files/passHot2Plot.bin",status="old",form="unformatted",access="sequential",iostat = ok)
  do i=1,nbseismes
    read(400)D(i)%datatps
  enddo
  read(400)nmod
  read(400)p
  read(400)misfit
  read(400)param_best
  read(400)xmin
  read(400)xmax
  read(400)acceptance
  read(400)dc
  do i=1,nbseismes
    read(400)pEpis(i)%nb
    allocate(pEpis(i)%pEpi(pEpis(i)%nb))
    do j=1,pEpis(i)%nb
      read(400)pEpis(i)%pEpi(j)
    enddo
  enddo
  read(400)mb
  read(400)acentroid
  close(400)
  if (ok.ne.0) write(*,*)'problème dans che_plot : le fichier OUTPUT/files/passHot2Plot.txt n''existe pas '

  ! ---------------------------------------------------------------------    .
  ! -------  traitement a posteriori des modèles sélectionnés    --------    .
  ! ---------------------------------------------------------------------    .
  call lect_mod_select(p,dp,nbChaineMVhot,misfit,param_best)                 ! relecture des modèles sélectionnés dans les fichier de sortie McMC
  ! ------- tri des meilleurs modèles                            --------    .
  call print_mess_4
  dp%deltaxy=150                                                             ! maillage X et Y du diagramme de densité
  call moy_mod_select(dp,nbChaineMVhot,param_best,misfit)                    ! calcule les moy, et et mode pour chaque paramètre

  ! ---------------------------------------------------------------------    .
  ! -------  géneration des densité à priori                     --------    .
  ! ---------------------------------------------------------------------    . ecriture pour che_apriori.90
  do i=1,nbseismes
    write(nbdeseismes(1:5),'(i5)')i
    open(unit=500,file="OUTPUT/files/Plot/apriori_"//trim(adjustl(nbdeseismes))//".bin", &
      status="replace",form="unformatted",access="sequential")
    write(500)i,nbtps,dp%temps_ref,mb
    do k=1,nbseismes
      write(500)D(k)%datatps
    enddo
    do k=1,nbseismes
      write(500)pEpis(k)%nb
      do j=1,pEpis(k)%nb
        write(500)pEpis(k)%pEpi(j)
      enddo
    enddo
    write(500)acentroid
    close(500)
  enddo

  ! ---------------------------------------------------------------------    .
  ! ------ moyenne 100 meilleurs modèles avec toutes les stations -------    .
  ! ---------------------------------------------------------------------    .
  call mksynthallsta(acentroid,dp)


  ! ------------------------------------------------------------------------- .
  ! étude des gradients sur la fonction coût
  ! ------------------------------------------------------------------------- .
  !if (plotposteriori) then
  !  nbmod=5000
  !  allocate(modelesIN(nbmod))
  !  do i=1,nbmod
  !    modelesIN(i)%VC=dp%VC%vec10000(i,1)
  !    modelesIN(i)%VM=dp%VM%vec10000(i,1)
  !    modelesIN(i)%Zmoho=dp%Zmoho%vec10000(i,1)
  !    modelesIN(i)%VpVs=dp%VpVs%vec10000(i,1)
  !    modelesIN(i)%Lat(:)=dp%Lat(:)%vec10000(i,1)
  !    modelesIN(i)%Lon(:)=dp%Lon(:)%vec10000(i,1)
  !    modelesIN(i)%Zhypo(:)=dp%Zhypo(:)%vec10000(i,1)
  !    modelesIN(i)%Tzero(:)=dp%temps_ref(:)
  !    modelesIN(i)%Tzero(:)%sec=modelesIN(i)%Tzero(:)%sec+dp%Tzero(:)%vec10000(i,1)
  !    do j=1,nbseismes
  !      call basetime(modelesIN(i)%Tzero(j))
  !    enddo
  !  enddo
  !  nomfichier='POST_HOTS_i'
    !call PosterioriExploration(p,nbmod,pEpis,nbtps,D,acentroid,xmin,xmax,mb,nomfichier,modelesIN)
  !  deallocate(modelesIN)
  !endif
  ! ------------------------------------------------------------------------- .




  ! ---------------------------------------------------------------------    .
  ! -------  géneration des scripts                              --------    .
  ! ---------------------------------------------------------------------    .
  ! ------- diagramme de densité, scripts GMT                    --------    .
  call print_mess_5
  call GMTfull(dp,nmod,nbChaineMVhot,xmin,xmax,nbtps,nbsta,D,E,nomsta,acentroid)
  ! ------- scriptes LaTeX                                       --------    .
  call latexfull(dc,dp,xmin,xmax,nbChaineMVhot,acceptance,param_best,misfit,E,nomsta,acentroid,nbtps,D)
  ! ---------------------------------------------------------------------    .

  ! ---------------------  fin du programme -----------------------------    .
  deallocate(dp%mis%vec,dp%VC%vec,dp%VM%vec,dp%Zmoho%vec,dp%VpVs%vec)
  do i=1,nbseismes
    deallocate (D(i)%datatps,pEpis(i)%pEpi)
    deallocate(dp%Lat(i)%vec,dp%Lon(i)%vec,dp%Zhypo(i)%vec,dp%Tzero(i)%vec)
  enddo
  deallocate(param_best,misfit,acceptance,nmod)
  if(allocated(nomsta)) deallocate(nomsta)
  ! ---------------------------------------------------------------------    .
  call print_line
  ! ---------------------------------------------------------------------    .

end program che2013_plot



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .



