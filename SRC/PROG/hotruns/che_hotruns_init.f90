! programme principal II
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
! -----------------------------------------------------------------------    .
! Initialisation des Hotruns                                                 !
!                                                                            !
program che2013_hotruns
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use mt19937
  use tirage
  use time
  use affiche
  use sub_param
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none
  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(parametres), dimension(:), allocatable :: param_init                  ! paramètres d'inv.
  type(fcout) :: misfit                                                      ! fonction coût
  type(accept) :: acceptance                                                 ! acceptance
  type(parametresinv) :: p                                                   ! paramètres d'inv.
  type(coldmoy) :: dc                                                        ! modèles du coldrun
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  ! type(date_secPS) :: midi20
  ! -------                                                      --------    .
  real(KIND=wr) :: xmin(nbseismes), xmax(nbseismes)                          ! cercles pond.
  real(KIND=wr) :: val1, val2
  ! real(KIND=wr) :: deltaP,deltaS
  ! -------                                                      --------    .
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: i,j,k,ok
  integer(KIND=wi) :: nbChaineMVhot,nbChaineMVcold                           ! nombre chaînes
  integer(KIND=wi) :: maxiterhot,maxitercold                                 ! nombre d'itérations
  integer(KIND=wi) :: nbsta,nbtps(nbseismes)                                 ! nombre station possible et reel nombre de données de temps
  ! -------                                                      --------    .
  character (LEN=5) :: numberchaine
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  ! -------  nb de chaines de Markov et d'iterations par chaine  --------    .
  call lectparam(nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold,chut=.true.)
  allocate(param_init(nbChaineMVcold))
  ! ---------------------------------------------------------------------    . lecture des coldruns
  ok = 0
  open(unit=200,file="OUTPUT/files/passCold2Hot.bin",status="old",form="unformatted",access="sequential",iostat = ok)
    read(200)nbsta,nbtps
    do i=1,nbseismes
      allocate(D(i)%datatps(nbtps(i)))
      do j=1,nbtps(i)
        read(200)D(i)%datatps(j)
      enddo
    enddo
    read(200)param_init
    read(200)xmin,xmax
    read(200)dc
    do i=1,nbseismes
      read(200)pEpis(i)%nb
      allocate(pEpis(i)%pEpi(pEpis(i)%nb))
      do j=1,pEpis(i)%nb
        read(200)pEpis(i)%pEpi(j)
      enddo
    enddo
    read(200)mb
    read(200)acentroid
  close(200)
  if (ok.ne.0) then
    write(*,*)'problème dans che_hotruns : le fichier OUTPUT/files/passCold2Hot.bin n''existe pas '
  endif
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .
  ! ------- début MCMc                                           --------    . HOTRUNS
  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .

    ! ------- lecture du prior                                   --------    .
    do i=1,nbChaineMVhot
      call lect_prior(p,param_init(i),"H")                                   ! H -> hotruns
    enddo

    ! ------- initialisation des autres paramètres               --------    .
    call init_div(misfit,acceptance)
    ! -------------------------------------------------------------------    . parametres hypocentraux
    ! ------- initialisation des bornes temporelles              --------    . 
    !
    ! centre du prior correspond à la moyenne des coldruns sélectionnées,
    ! plus ou moin 3 écart-types (plus une demi seconde)
    !
    do j=1,nbseismes
      p%maxi%Tzero(j) = dc%tempsrefcold(j)                                   ! borne sup.
      p%maxi%Tzero(j)%sec = dc%moyselect%par%Tzero(j)%sec + max(3.00_wr*dc%ecselect%par%Tzero(j)%sec+0.75_wr,0.75_wr)  ! 3 ecartypes superieurs minimum : 0,75 seconde
      p%maxi%Tzero(j)%sec = real(int(p%maxi%Tzero(j)%sec)+1,wr)
      call basetime(p%maxi%Tzero(j))
      p%mini%Tzero(j) = dc%tempsrefcold(j)                                   ! borne inf.
      p%mini%Tzero(j)%sec =  dc%moyselect%par%Tzero(j)%sec - max(3.00_wr*dc%ecselect%par%Tzero(j)%sec+0.75_wr,0.75_wr) ! 3 ecartypes inférieurs minimum : 0,75 seconde
      p%mini%Tzero(j)%sec = real(int(p%mini%Tzero(j)%sec),wr)
      call basetime(p%mini%Tzero(j))
    enddo
    ! ------- initialisation des bornes épicentales              --------    . centre du prior correspond à la moyenne des coldruns sélectionnées, plus ou moin 3 écart-types
    p%centreY(:) = dc%moyselect%par%lat(:)
    p%centreX(:) = dc%moyselect%par%lon(:)
    ! ------- recherche de l'épicentre dans ce rayon             --------    . rayon correspond à 3 sigma des coldruns
    do j=1,nbseismes
      val1 = 3.00_wr*dc%ecselect%par%lat(j)*(pi*rT)/180.0_wr                 ! en km
      val2 = 3.00_wr*dc%ecselect%par%lon(j)*pi*(cos(dc%moyselect%par%lat(j)/180.0_wr*pi)*6371.0_wr)/180.0_wr ! en km
      p%Rayon(j) = max(val1,val2,5.0_wr)                                     ! diamètre minimum : 10,0 km
    enddo
    !
    ! -------------------------------------------------------------------    . si certains paramètres fixes ...
    call paramfixe(p)
    ! -------------------------------------------------------------------    .


    ! --------------------------------------------------------------------- .
    ! Carriere
    ! --------------------------------------------------------------------- .
    !do i=1,nbseismes
    !  midi20=D(i)%datatps(1)%tpsR
    !  midi20%date%hour=11                                                   ! midi TU, pour Vannes en hiver
    !  midi20%date%min=20
    !  call basetime(midi20)
    !  call difftime(deltaP,deltaS,midi20,D(i)%datatps(1)%tpsR)
    !  if ((abs(deltaP).lt.1.00_wr).or.(i==45)) then
    !    write(*,*)'seisme ',i,' : tire de carriere ?'
    !    write(*,*)'Zhypo max = 6 km'
    !    p%maxi%Zhypo(i) = 6.d0
    !    p%ecartype%Zhypo(i) =  p%ecartype%Zhypo(i)/4.0_wr
    !  endif
    !  ! ------------------------------------------------------------------ .
    !enddo
    ! --------------------------------------------------------------------- .
    ! --------------------------------------------------------------------- .

    ! production de fichier
    ! parce que dans la norme MPI-Fortran 2008,
    ! faire un call MPI_BCAST avec des types dérivés, c'est pas prévu !

    ! --------------------  un fichier par hotruns  ---------------------    .
    do k=1,nbChaineMVhot
      write(numberchaine(1:5),'(i5)')k
      open(unit=100,file="OUTPUT/files/Hot/In_"//trim(adjustl(numberchaine))//".bin", &
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
        write(100)dc
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
  do i=1,nbseismes
    deallocate (D(i)%datatps,pEpis(i)%pEpi)
  enddo
  deallocate(param_init)
  ! ---------------------------------------------------------------------    .

end program che2013_hotruns



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .

