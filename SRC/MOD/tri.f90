! Librairie de subroutines avec Algorithmes_de_tri/Tri_rapide (ou pas ...)
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE tri

    use modparam

    implicit none

    private

    public  :: tri_bulle
    public  :: tridata
    public  :: triparam
    public  :: melangetab


    ! -------                                                  --------    .
    interface tri_bulle
      module procedure tri_bulle_reel, tri_bulleRstring
    end interface tri_bulle
    ! -------                                                  --------    .
    interface triparam
      module procedure triparams, triparam1
    end interface triparam
    ! -------                                                  --------    .
    interface melangetab
      module procedure melangetab1, melangetab3
    end interface melangetab
    ! -------                                                  --------    .


CONTAINS

    ! -----------------------------------------------------------------    .


    subroutine tri_bulle_reel(InList,n,OutList)
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent(in) :: n
    real(KIND=wr), dimension(:,:), intent(in) :: InList
    real(KIND=wr), dimension(:,:), intent(out) :: OutList
    ! -----------------------------------------------------------------    .
    if (size(InList,2)==2) then
      call tri_bulle2col(InList,n,OutList)
    else if (size(InList,2)==3) then
      call tri_bulle3col(InList,n,OutList)
    else
      write(*,*)'problème dans tri_bulle, tableau non conforme ',size(InList,2)
      stop
    endif
    ! -----------------------------------------------------------------    .

    contains

      subroutine tri_bulle2col(InList,n,OutList)
        ! -------                                              --------    .
        ! permet le tri croissant d'un vecteur de 2 dimensions 
        ! en fonction de la seconde colonne
        ! -------                                              --------    .
        implicit none
        ! -------                                              --------    .
        integer(KIND=wi), intent(in) :: n
        real(KIND=wr), dimension(n,2), intent(in) :: InList
        real(KIND=wr), dimension(n,2), intent(out) :: OutList
        ! -------                                              --------    .
        real(KIND=wr) :: pass(2)
        integer(KIND=wi) :: i,j
        real(KIND=wr) :: delta
        ! -------------------------------------------------------------    .
        OutList=InList
        do j=2,n
          do i=j,2,-1
            delta=OutList(i,2)-OutList(i-1,2)
            if (delta.lt.0.0_wr) then
              pass(:) = OutList(i-1,:)
              OutList(i-1,:) = OutList(i,:)
              OutList(i,:) = pass(:)
            endif
          enddo
        enddo
        ! -------------------------------------------------------------    .
      end subroutine tri_bulle2col

        ! -------------------------------------------------------------    .

      subroutine tri_bulle3col(InList,n,OutList)
        ! -------                                              --------    .
        ! permet le tri croissant d'un vecteur de 3 dimensions
        ! en fonction de la troisième colonne
        ! -------                                              --------    .
        implicit none
        ! -------                                              --------    .
        integer(KIND=wi), intent(in) :: n
        real(KIND=wr), dimension(n,3), intent(in) :: InList
        real(KIND=wr), dimension(n,3), intent(out) :: OutList
        ! -------                                              --------    .
        real(KIND=wr) :: pass(3)
        integer(KIND=wi) :: i,j
        real(KIND=wr) :: delta
        ! -------------------------------------------------------------    .
        OutList=InList
        do j=2,n
          do i=j,2,-1
            delta=OutList(i,3)-OutList(i-1,3)
            if (delta.lt.0.0_wr) then
              pass(:) = OutList(i-1,:)
              OutList(i-1,:) = OutList(i,:)
              OutList(i,:) = pass(:)
            endif
          enddo
        enddo
        ! -------------------------------------------------------------    .
      end subroutine tri_bulle3col

  end subroutine tri_bulle_reel

! ---------------------------------------------------------------------    .

  subroutine tri_bulleRstring(InOut_reel,InOut_string,n)
    ! -------                                                --------    .
    ! permet le tri Dé-croissant d'un vecteur de 2 dimensions
    ! en fonction de la seconde colonne
    ! -------                                                --------    .
    implicit none
    ! -------                                                --------    .
    integer(KIND=wi), intent(in) :: n
    real(KIND=wr), dimension(n), intent(inout) :: InOut_reel
    character(LEN=4), dimension(n), intent(inout) :: InOut_string
    ! -------                                                --------    .
    character(LEN=4) :: passS
    real(KIND=wr) :: passR
    integer(KIND=wi) :: i,j
    real(KIND=wr) :: delta
    ! ---------------------------------------------------------------    .
    do j=2,n
      do i=j,2,-1
        delta=InOut_reel(i)-InOut_reel(i-1)
        if (delta.gt.0.0_wr) then
          passR = InOut_reel(i-1)
          passS = InOut_string(i-1)
          InOut_reel(i-1) = InOut_reel(i)
          InOut_string(i-1) = InOut_string(i)
          InOut_reel(i) = passR
          InOut_string(i) = passS
        endif
      enddo
    enddo
    ! ---------------------------------------------------------------    .
  end subroutine tri_bulleRstring

! ---------------------------------------------------------------------    .

  subroutine tridata(nbtps,datatps)
    ! -------                                                  --------    .
    ! permet le tri croissant d'un vecteur de type dataone en fonction des temps d'arrivés des ondes P
    ! -------                                                  --------    .
    use typetemps, only : dataone
    use time, only : difftime
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi),intent (in) :: nbtps
    type(dataone),intent (inout) :: datatps(nbtps)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j
    type(dataone) :: datapass
    real(KIND=wr) :: deltaP,deltaS
    ! -----------------------------------------------------------------    .
    do j=2,nbtps
      do i=j,2,-1
        call difftime(deltaP,deltaS,datatps(i)%tpsR,datatps(i-1)%tpsR)
        if (deltaP.lt.0.0_wr) then
          datapass = datatps(i-1)
          datatps(i-1) = datatps(i)
          datatps(i) = datapass
        endif
      enddo
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine tridata

! ---------------------------------------------------------------------    .

  subroutine melangetab1(n,tab)
    ! -------                                                  --------    .
    ! permet de mélanger un tableau par le mélange de Fisher-Yates (ou de Knuth) :
    ! un algorithme générant une permutation aléatoire d'un ensemble fini
    ! -------                                                  --------    .
    use mt19937
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi),intent (in) :: n
    real(KIND=wr),intent (inout) :: tab(n)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,pif
    real(KIND=wr) :: pass
    ! -----------------------------------------------------------------    .
    do j=1,3                                                               ! 3 passages
      do i=1,n
        pif=1+int(genrand_real1()*real(n,wr),wi)                           ! equi-aléatoire de 1 à n
        pass=tab(i)
        tab(i)=tab(pif)
        tab(pif)=pass
      enddo
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine melangetab1

! ---------------------------------------------------------------------    .

  subroutine melangetab3(n,tab)
    ! -------                                                  --------    .
    ! permet de mélanger un tableau par le mélange de Fisher-Yates (ou de Knuth) :
    ! un algorithme générant une permutation aléatoire d'un ensemble fini
    ! -------                                                  --------    .
    use mt19937
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi),intent (in) :: n
    real(KIND=wr),intent (inout) :: tab(3,n)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,pif
    real(KIND=wr) :: pass(3)
    ! -----------------------------------------------------------------    .
    do j=1,3                                                               ! 3 passages
      do i=1,n
        pif=1+int(genrand_real1()*real(n,wr),wi)                           ! equi-aléatoire de 1 à n
        pass(:)=tab(:,i)
        tab(:,i)=tab(:,pif)
        tab(:,pif)=pass(:)
      enddo
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine melangetab3

! ---------------------------------------------------------------------    .

  subroutine triparams(nb,misfit,param_best)
    ! -------                                                  --------    .
    ! permet le tri croissant d'un jeu de paramètres en fonction de la fonction coût
    ! pour tous les séismes
    ! -----------------------------------------------------------------    .
    use typetemps, only : parametres, fcout
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: nb
    type(parametres), intent (inout) :: param_best(nb)
    type(fcout), intent (inout) :: misfit(nb)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j
    type(parametres) :: par
    type(fcout) :: mis
    ! -----------------------------------------------------------------    .
    do j=2,nb
      do i=j,2,-1
        if (misfit(i)%best.lt.misfit(i-1)%best) then
          mis = misfit(i-1)
          misfit(i-1) = misfit(i)
          misfit(i) = mis
          par = param_best(i-1)
          param_best(i-1) = param_best(i)
          param_best(i) = par
        endif
      enddo
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine triparams

! ---------------------------------------------------------------------    .

  subroutine triparam1(nb,misfit,param_best)
    ! -------                                                  --------    .
    ! permet le tri croissant d'un jeu de paramètre en fonction de la fonction coût
    ! pour 1 séisme 
    ! -----------------------------------------------------------------    .
    use typetemps, only : parametre
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: nb
    type(parametre), intent (inout) :: param_best(nb)
    real(KIND=wr), intent (inout)  :: misfit(nb)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j
    type(parametre) :: par
    real(KIND=wr) :: mis
    ! -----------------------------------------------------------------    .
    do j=2,nb
      do i=j,2,-1
        if (misfit(i).lt.misfit(i-1)) then
          mis = misfit(i-1)
          misfit(i-1) = misfit(i)
          misfit(i) = mis
          par = param_best(i-1)
          param_best(i-1) = param_best(i)
          param_best(i) = par
        endif
      enddo
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine triparam1

! ---------------------------------------------------------------------    .

END MODULE tri



! *********************************************************************    .
! *********************************************************************    .


