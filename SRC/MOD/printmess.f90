! Librairie de subroutines permettant l'impression de messages à l'écran
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE affiche

    use modparam

    implicit none

    private

    public  :: print_mess_1,print_mess_2,print_mess_2bis,print_mess_3, &
               print_mess_3bis,print_mess_4,print_mess_5
    public  :: print_mess_fin
    public  :: print_mess_finchainemin
    public  :: print_mess_finchainemax
    public  :: print_line
    public  :: print_messchaine


CONTAINS

! ---------------------------------------------------------------------    .


  subroutine print_mess_1
    ! -------                                                  --------    .mh
    implicit none                                                          ! début du programme che_coldruns_init
    write(*,*)
    write(*,*)'---------------------------------------------------'
    write(*,*)'--------------- CHE2016 version 1.6 ---------------'
    write(*,*)'---------------------------------------------------'
    write(*,*)'                      meric.haugmard@univ-nantes.fr'
    write(*,*)'                                     méric Haugmard'
    write(*,*)'                                          2013-2016'
    write(*,*)
    write(*,*)'------------  initialisation coldruns  ------------' ;
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_1

! ---------------------------------------------------------------------    .

  subroutine print_mess_fin
    ! -------                                                  --------    .mh
    implicit none                                                          ! fin du programme che_plot
    write(*,*)
    write(*,*)'---------------------------------------------------'
    write(*,*)'--------------------  fin prog  -------------------'
    write(*,*)'---------------------------------------------------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_fin

! ---------------------------------------------------------------------    .

  subroutine print_mess_2
    ! -------                                                  --------    .mh
    implicit none                                                          ! début McMC
    write(*,*)
    write(*,*)'-----------------  début coldruns  ----------------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_2

! ---------------------------------------------------------------------    .

  subroutine print_mess_2bis
    ! -------                                                  --------    .mh
    implicit none                                                          ! début McMC
    write(*,*)
    write(*,*)'-----------------  début hotruns  -----------------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_2bis

! ---------------------------------------------------------------------    .

  subroutine  print_mess_finchainemin(a,b)
    ! -------                                                  --------    .mh
    ! fin de chaque chaîne _ première étape
    implicit none
    real(KIND=wr) :: a,b
    ! -----------------------------------------------------------------    .
    write(*,'(a46,f7.2)')' fonction coût minimale         :             ',a
    write(*,'(a43,f7.2,a2)')' acceptance                     :           ',b,' %'
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_finchainemin

! ---------------------------------------------------------------------    .

  subroutine  print_mess_finchainemax(a,b,c)
    ! -------                                                  --------    .mh
    ! fin de chaque chaîne _ seconde étape
    implicit none
    real(KIND=wr) :: a,b
    integer(KIND=wi) :: c
    ! -----------------------------------------------------------------    .
    write(*,'(a40,i15)')' nombre de modèles sélectionnés :    ',c
    write(*,'(a46,f7.2)')' fonction coût minimale         :             ',a
    write(*,'(a43,f7.2,a2)')' acceptance                     :           ',b,' %'
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_finchainemax

! ---------------------------------------------------------------------    .

  subroutine print_mess_3
    ! -------                                                  --------    .mh
    implicit none                                                          ! fin McMC
    write(*,*)
    write(*,*)'------------------  fin coldruns  -----------------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_3

! ---------------------------------------------------------------------    .

  subroutine print_mess_3bis
    ! -------                                                  --------    .mh
    implicit none                                                          ! fin McMC
    write(*,*)
    write(*,*)'------------------  fin hotruns  ------------------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_3bis

! ---------------------------------------------------------------------    .

  subroutine print_mess_4
    ! -------                                                  --------    .mh
    implicit none                                                          ! calcul des moyennes
    write(*,*)
    write(*,*)'--------------- calcul a posteriori ---------------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_4

! ---------------------------------------------------------------------    .

  subroutine print_mess_5
    ! -------                                                  --------    .mh
    implicit none                                                          ! production des script pour les figures
    write(*,*)
    write(*,*)'----- production des script pour les figures ------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_mess_5

! ---------------------------------------------------------------------    .

  subroutine print_line
    ! -------                                                  --------    .mh
    implicit none                                                          ! entre chaque chaîne
    write(*,*)
    write(*,*)'---------------------------------------------------'
    write(*,*)
    ! -----------------------------------------------------------------    .
  end subroutine print_line

! ---------------------------------------------------------------------    .

  subroutine print_messchaine(i,n)
    ! -------                                                  --------    .mh
    ! en début de chaîne
    implicit none
    integer(KIND=wi), intent(in) :: i,n
    ! -----------------------------------------------------------------    .
    write(*,*)
    write(*,*)'-- chaîne numéro :   ',i,'/',n,' --'
    ! -----------------------------------------------------------------    .
  end subroutine print_messchaine

END MODULE affiche



! *********************************************************************    .
! *********************************************************************    .


