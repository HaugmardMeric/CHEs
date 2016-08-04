! programme principal I ter
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
! -----------------------------------------------------------------------    .
! Synthèse des Coldruns                                                      !   
!                                                                            !
program che2013_coldruns_init
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use affiche
  use datalecture
  use sub_param
  use tri
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none
  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(parametres), dimension(:), allocatable :: param_init, param_best      ! paramètres d'inv.
  type(fcout), dimension(:), allocatable :: misfit                           ! fonction coût
  type(accept), dimension(:), allocatable :: acceptance                      ! acceptance
  type(parametresinv) :: p                                                   ! paramètres d'inv.
  type(coldmoy) :: dc                                                        ! moyennes et écarts-types des modèles du coldrun
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  ! -------                                                      --------    .
  real(KIND=wr), dimension(:), allocatable :: vec
  real(KIND=wr) :: moy,ec
  real(KIND=wr) :: xmin(nbseismes), xmax(nbseismes)                          ! cercles pond.
  ! -------                                                      --------    .
  integer(KIND=wi) :: i,j,k,ok
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: nbChaineMVhot,nbChaineMVcold                           ! nombre chaînes
  integer(KIND=wi) :: maxiterhot,maxitercold                                 ! nombre d'itérations
  integer(KIND=wi) :: nbsta, nbtps(nbseismes)                                ! nombre station et nombre de données de temps par séismes]
  ! -------                                                      --------    .
  character (LEN=5) :: chaine
  ! ---------------------------------------------------------------------    .
  call lectparam(nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold,chut=.true.)

  allocate(param_best(nbChaineMVcold),param_init(nbChaineMVcold))
  allocate(misfit(nbChaineMVcold),acceptance(nbChaineMVcold))

  do k=1,nbChaineMVcold
    write(chaine(1:5),'(i5)')k
    open(unit=100,file="OUTPUT/files/Cold/Out_"//trim(adjustl(chaine))//".bin", &
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
  if (ok .ne. 0) then
    write(*,*)"problème dans che2013_coldruns_syn : le fichier ", &
      "OUTPUT/files/Cold/In_"//trim(adjustl(chaine))//".bin n''existe pas "
    stop
  endif
  ! ------- tri des colds runs                                   --------    .
  call triparam(nbChaineMVcold,misfit,param_best)
  ! ------- calcul des moyennes pour les coldruns                --------    .
  call moycoldruns(nbChaineMVcold,param_best,misfit,nbChaineMVhot,dc)

  ! ----------------------------  ecriture ------------------------------    .
  allocate(vec(nbChaineMVcold))
  do j=1,nbChaineMVcold
   vec(j)=acceptance(j)%val
  enddo
  call moy_ec (vec,nbChaineMVcold,nbChaineMVcold,moy,ec)
  deallocate(vec)
  write(*,1111)' acceptance (%)                 : ',moy,' +ou- ',ec
  write(*,1111)' fonction coût minimale TOTALE  : ',dc%moytot%mis,' +ou- ',dc%ectot%mis
  write(*,1111)' fonction coût minimale SELECT  : ',dc%moyselect%mis,' +ou- ',dc%ecselect%mis

  ! ---------------------  sauve pour che_hotruns -----------------------    .
  open(unit=100,file="OUTPUT/files/passCold2Hot.bin",status="replace",form="unformatted",access="sequential")
    write(100)nbsta,nbtps
    do i=1,nbseismes
      do j=1,nbtps(i)
        write(100)D(i)%datatps(j)
      enddo
    enddo
    write(100)param_best
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

  call print_mess_3
  ! ---------------------  fin du programme -----------------------------    .
  do i=1,nbseismes
    deallocate (D(i)%datatps,pEpis(i)%pEpi)
  enddo
  deallocate(param_best,param_init,misfit,acceptance)
  ! ---------------------------------------------------------------------    .
1111 format(a,f8.3,a,f6.2)
  ! ---------------------------------------------------------------------    .

end program che2013_coldruns_init



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .


