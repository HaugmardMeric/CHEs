! programme principal sac_verifmediatrice
! ***********************************************************************   .mh
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! -------                                                        --------   .
! -------                                                        --------   .
! ------- avril 2015                                             --------   .
! -------                                                        --------   .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr           --------   .
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
!  This program is distributed for research purposes and in the hope        !
!  that it will be useful however without any warranties.                   !
!  Use it on your own risk. Please, report bugs/improvements/... .          !
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
!                                                                           !
!
program verifmediatrice
  ! --------------------------------------------------------------------    .
  ! carte de recherches epicentrales a priori
  ! --------------------------------------------------------------------    .
  use modparam
  use typetemps
  use datalecture
  use rechercheepi
  ! -------                                                      --------   .
  implicit none
  ! -------                                                      --------   .
  type(dataall) :: D(nbseismes)                                             ! données de temps
  type(priorEPI) :: pEpis(nbseismes)                                        ! prior
  integer(KIND=wi) :: i, mb, nbsta, nbtps(nbseismes)
  ! --------------------------------------------------------------------    .
  ! existe OUTPUT/GMT/ et  OUTPUT/figures/
  ! existe DATA/sta.d et DATA/sta.d et DATA/seismes.d DATA/xxxx.xx.xx.xx.xx.d
  ! -------                                                      --------   .
  ! ------- lecture des données                                 --------    . phases and stations list
  call lectnbdata(nbsta,nbtps)                                              ! nombre de données
  do i=1,nbseismes
    allocate(D(i)%datatps(nbtps(i)))                                        ! alloue par séisme
  enddo
  call lectdata(nbsta,nbtps,D)
  ! -------                                                      --------   .
  print*,'initialisation du prior pour les paramètres épicentraux'
  call zoneRecherche(nbtps,D,pEpis,mb)
  ! -------                                                      --------   .
  write(*,*)'chmod +x OUTPUT/GMT/script*.sh ; ./OUTPUT/GMT/script*.sh '     ! chmod +x ; plot it
  ! -------                                                      --------   .
  ! --------------------------------------------------------------------    .
end program verifmediatrice

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   .
! ***********************************************************************   .

