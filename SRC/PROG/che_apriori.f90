! programme principal IV
! ***********************************************************************    .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .
! -------                                                        --------    .
! ------- CHE2013_apriori version 1.5                            --------    .
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
!                                                                            !
!                                                                            !
! -----------------------------------------------------------------------    .
! période d'exploitation des distributions a priori                          !
!                                                                            !
program che2013_apriori
  ! ---------------------------------------------------------------------    .
  ! ------- modules :                                            --------    .
  use modparam
  use typetemps
  use sub_param
  use mt19937
  use affiche
  ! ---------------------------------------------------------------------    .
  ! ------- déclaration :                                        --------    .
  implicit none

    include 'mpif.h'

  ! -------                                                      --------    .
  type(dataall) :: D(nbseismes)                                              ! données de temps
  type(densityplot) :: dp
  type(priorEPI) :: pEpis(nbseismes)                                         ! prior
  type(amoho_centroid) :: acentroid                                          ! si moho non tabulaire
  ! -------                                                      --------    .
  integer(KIND=wi) :: i,j,k,ok
  integer(KIND=wi) :: mb                                                     ! prior
  integer(KIND=wi) :: nbm_ap                                                 ! nombre de modèles a priori générés
  integer(KIND=wi) :: nbtps(nbseismes)                                       ! nombre station et nombre de données de temps
  ! -------                                                      --------    .
  character (LEN=5) :: nbdeseismes
  ! ---------------------------------------------------------------------    .
  ! -------                                                      --------    . MPI :
  integer :: nb_procs,rang,code,err
  integer :: seed
  integer, dimension(:), allocatable :: allseed
  ! ---------------------------------------------------------------------    . MPI_BEGIN
  call MPI_INIT(code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  ! ---------------------------------------------------------------------    . MPI_BEGIN
  ! ------- initialisation de la graine                          --------    .
  if (rang==0) then
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

  ! ---------------------------------------------------------------------    . relecture des données
  ok=0
  write(nbdeseismes(1:5),'(i5)')rang+1
  open(unit=500+rang,file="OUTPUT/files/Plot/apriori_"//trim(adjustl(nbdeseismes))//".bin", &
    status="old",form="unformatted",access="sequential",iostat = ok)
  read(500+rang)k,nbtps,dp%temps_ref,mb
  do i=1,nbseismes
    allocate(D(i)%datatps(nbtps(i)))
    read(500+rang)D(i)%datatps
  enddo
  do i=1,nbseismes
    read(500+rang)pEpis(i)%nb
    allocate(pEpis(i)%pEpi(pEpis(i)%nb))
    do j=1,pEpis(i)%nb
      read(500+rang)pEpis(i)%pEpi(j)
    enddo
  enddo
  read(500+rang)acentroid
  close(500+rang)
  if (ok.ne.0) then
    write(*,*)"problème dans che_apriori : le fichier OUTPUT/files/Plot/", &
      "apriori_"//trim(adjustl(nbdeseismes))//".bin n''existe pas "
    call MPI_ABORT(MPI_COMM_WORLD,err,code)
  endif
  ! ---------------------------------------------------------------------    .
  ! -------  géneration des densité à priori                     --------    .
  ! ---------------------------------------------------------------------    .
  nbm_ap=100000
  call dist_apriori(k,rang,nbtps,D,nbm_ap,dp%temps_ref,pEpis,mb,acentroid)

  ! ---------------------------------------------------------------------    .
  ! ---------------------  fin du programme -----------------------------    .
  do i=1,nbseismes
    deallocate (D(i)%datatps,pEpis(i)%pEpi)
  enddo
  ! ---------------------------------------------------------------------    .
  if (rang==0) call print_mess_fin
  ! ---------------------------------------------------------------------    . MPI_END
  call MPI_FINALIZE(code)
  ! ---------------------------------------------------------------------    . MPI_END

  ! FIN

end program che2013_apriori



! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    .



