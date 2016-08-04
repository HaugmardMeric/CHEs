! programme principal II ter
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
! Synthèse des Hotruns                                                       !
!                                                                            !
program che2013_hotruns
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use affiche
  use sub_param
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none
  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(parametres), dimension(:), allocatable :: param_init,param_best       ! paramètres d'inv.
  type(fcout), dimension(:), allocatable :: misfit                           ! fonction coût
  type(accept), dimension(:), allocatable :: acceptance                      ! acceptance
  type(parametresinv) :: p                                                   ! paramètres d'inv.
  type(coldmoy) :: dc                                                        ! modèles du coldrun
  type(residus), dimension(:), allocatable :: R                              ! résidus (si FLAGresSTA=.true.)
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  ! -------                                                      --------    .
  real(KIND=wr), dimension(:), allocatable :: vec
  real(KIND=wr) :: moy,ec
  real(KIND=wr) :: xmin(nbseismes), xmax(nbseismes)                          ! cercles pond.
  ! -------                                                      --------    .
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: i,j,k,ok
  integer(KIND=wi) :: nbChaineMVhot,nbChaineMVcold                           ! nombre chaînes
  integer(KIND=wi) :: maxiterhot,maxitercold                                 ! nombre d'itérations
  integer(KIND=wi) :: nbsta,nbtps(nbseismes)                                 ! nombre station possible et reel nombre de données de temps
  integer(KIND=wi), dimension(:), allocatable :: nmod                        ! nombre modèles sélectionnés par chaine
  ! -------                                                      --------    .
  ! -------                                                      --------    .
  character (LEN=5) :: numberchaine
  ! ---------------------------------------------------------------------    .

  ! ---------------------------------------------------------------------    .
  ! -------  nb de chaines de Markov et d'iterations par chaine  --------    .
  call lectparam(nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold,chut=.true.)
  allocate(param_best(nbChaineMVhot),param_init(nbChaineMVhot))
  allocate(misfit(nbChaineMVhot),acceptance(nbChaineMVhot))
  allocate(nmod(nbChaineMVhot))

  ok=0
  do k=1,nbChaineMVhot
    write(numberchaine(1:5),'(i5)')k
    ! ----------------  un fichier OUTPUT par cold runs  ----------------    .
    open(unit=100,file="OUTPUT/files/Hot/Out_"//trim(adjustl(numberchaine))//".bin", &
        status="old",form="unformatted",access="sequential",iostat = ok)
      read(100)nbsta,nbtps
      read(100)p
      do i=1,nbseismes
        if (.not.(allocated(D(i)%datatps))) allocate(D(i)%datatps(nbtps(i)))
        do j=1,nbtps(i)
          read(100)D(i)%datatps(j)
        enddo
      enddo
      read(100)misfit(k)
      read(100)acceptance(k)
      read(100)dc
      read(100)nmod(k)
      read(100)param_init(k)
      read(100)param_best(k)
      read(100)xmin,xmax
      do i=1,nbseismes
        read(100)pEpis(i)%nb
        if (.not.(allocated(pEpis(i)%pEpi))) allocate(pEpis(i)%pEpi(pEpis(i)%nb))
        do j=1,pEpis(i)%nb
          read(100)pEpis(i)%pEpi(j)
        enddo
      enddo
      read(100)mb
      read(100)acentroid
    close(100)
  enddo
  if (ok.ne.0) then
    write(*,*)"problème dans che_hotruns : le fichier OUTPUT/files/Hot/Out_"//trim(adjustl(numberchaine))//".bin n''existe pas "
    stop
  endif


  ! ----------------------------  ecriture ------------------------------    .
  allocate(vec(nbChaineMVhot))
  do j=1,nbChaineMVhot
    vec(j)=acceptance(j)%val
   enddo
  call moy_ec (vec,nbChaineMVhot,nbChaineMVhot,moy,ec)
  write(*,1111)' acceptance (%)                 : ',moy,' +ou- ',ec
  ! -------                                                      --------    .
  do j=1,nbChaineMVhot
    vec(j)=misfit(j)%best
  enddo
  call moy_ec (vec,nbChaineMVhot,nbChaineMVhot,moy,ec)
  deallocate(vec)
  write(*,1111)' fonction coût minimale         : ',moy,' +ou- ',ec


  call print_mess_3bis

  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .
  ! ------- fin MCMc                                             --------    . HOTRUNS
  ! ---------------------------------------------------------------------    .
  ! ---------------------------------------------------------------------    .

  ! ---------------------  sauve pour che_plot --------------------------    .
  open(unit=204,file="OUTPUT/files/passHot2Plot.bin",status="replace",form="unformatted",access="sequential")
    do i=1,nbseismes
      write(204)D(i)%datatps
    enddo
    write(204)nmod
    write(204)p
    write(204)misfit
    write(204)param_best
    write(204)xmin
    write(204)xmax
    write(204)acceptance
    write(204)dc
    do i=1,nbseismes
      write(204)pEpis(i)%nb
      do j=1,pEpis(i)%nb
        write(204)pEpis(i)%pEpi(j)
      enddo
    enddo
    write(204)mb
    write(204)acentroid
  close(204)
  ! ---------------------  fin du programme -----------------------------    .
  do i=1,nbseismes
    deallocate (D(i)%datatps,pEpis(i)%pEpi)
  enddo
  deallocate(param_best,param_init,misfit,acceptance,nmod)
  if(allocated(R)) deallocate(R)
  ! ---------------------------------------------------------------------    .
1111 format(a,f8.3,a,f6.2)
  ! ---------------------------------------------------------------------    .

end program che2013_hotruns



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .

