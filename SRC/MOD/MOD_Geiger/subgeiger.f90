! septembre 2013 - fevrier 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .
!
! résolution d'un problème inverse : 
!
! recherche des parametres hypocentraux en admettant un modèle de 
! Terre, par la méthode Geiger (1910, 1912), méthode itérative par 
! moindres carrées, pour un ou plusierus seismes

MODULE invGEIGER

    use modparam
    use typetemps, only : date_sec, dataone, dataall, parametre, parametres
    use typetemps, only : mvPall_2_P1, mvP1_2_Pall, amoho_centroid
    use mt19937, only : genrand_real3, normal
    use time, only : basetime, difftime


    implicit none

    private

    public  :: dogeiger                                                    ! inversion par la méthode Geiger avec un modèle de terre et hypocentral connu
    public  :: dogeigerone

    ! -----------------------------------------------------------------    .
    ! mod (character (LEN=1)):
    ! -------                                                  --------    .
    ! 'I' pour un modèle simple                                            ! 2 couches -> pb direct CHEs !
    ! 'A' pour les modèles de Terre Arroucau, phD 2006                     ! 9 couches
    ! 'S' pour le modèles de Terre Si-Hex type Haslach                     ! 3 couches
    ! 'F' pour un modèles de Terre test et fun                             ! x couches
    ! 'C' pour un modèles de Terre CEA                                     ! 3 couches
    ! -----------------------------------------------------------------    .


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine dogeiger(n,D,params,acentroid,mod,chut,amisfit)
    ! -------                                                  --------    .mh
    ! dogeiger pour nbseismes seismes
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent (in) :: n(nbseismes)                          ! nombre de données de temps
    type(dataall), intent (inout) :: D(nbseismes)                          ! données
    type(parametres), intent (inout) :: params                             ! paramètres
    type (amoho_centroid), intent (in) :: acentroid
    character (LEN=1), intent(in) :: mod
    logical, intent(in), optional :: chut
    real(KIND=wr), intent(out), optional :: amisfit(nbseismes)
    ! -------                                                  --------    .
    type(parametre) :: param_init                                          ! paramètres
    integer(KIND=wi) :: j
    real(KIND=wr) :: onemisfit
    logical :: con    
    ! -----------------------------------------------------------------    .
    do j=1,nbseismes
      call  mvPall_2_P1(param_init,params,j)
      if(present(chut)) then
        call dogeigerone(j,n(j),D(j)%datatps,param_init,acentroid,mod,onemisfit,chut=chut,con=con)
      else
        call dogeigerone(j,n(j),D(j)%datatps,param_init,acentroid,mod,onemisfit,con=con)
      endif
      call mvP1_2_Pall(params,param_init,j)
      if(present(amisfit)) amisfit(j)=onemisfit
    enddo
    ! -----------------------------------------------------------------    .
    end subroutine dogeiger

! ---------------------------------------------------------------------    .

  subroutine dogeigerone(j,nbtps,datatps,param_one,acentroid,mod,onemisfit,chut,con)
    ! -------                                                  --------    .mh
    ! applique la méthode Geiger, en admettant un modèle de Terre fixe (param_one), pour un unique séisme.
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent (in) :: j
    integer(KIND=wi), intent (in) :: nbtps                                 ! nb de données de temps
    type(dataone), intent (inout) :: datatps(nbtps)                        ! données
    type(parametre), intent (inout) :: param_one                           ! paramètres
    type (amoho_centroid), intent (in) :: acentroid
    character (LEN=1), intent(in) :: mod
    real(KIND=wr), intent(out) :: onemisfit
    logical, intent(in), optional :: chut
    logical, intent(out), optional :: con                                  ! "true" si convergence
    ! -------                                                  --------    .
    integer(KIND=wi) i
    ! -----------------------------------------------------------------    .
    i=0
    if(present(chut)) then
        call geiger(nbtps,datatps,param_one,con,i,acentroid,mod,onemisfit,chut)
    else
        call geiger(nbtps,datatps,param_one,con,i,acentroid,mod,onemisfit)
        if (.not.con) write(*,*)'problème dans dogeiger : données insufiantes, méthode Geiger fixe divergente',j
    endif
    ! -----------------------------------------------------------------    .
  end subroutine dogeigerone

! ---------------------------------------------------------------------    .

  subroutine geiger(nbtps,datatps,param_init,con,nk,acentroid,mod,onemisfit,chut)
    ! -------                                                  --------    .mh
    ! application de la méthode itérative de Geiger après initialisation des paramètres
    ! con = .true. : si convergence
    ! con = .false. : si divergence
    ! -------                                                  --------    .
    use sub_misfit, only : compute_misfitone
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    TYPE elementmatrice
      real(KIND=wr) :: df_Dto,df_Dhypo_P,df_Dhypo_S,df_Dlon_P,df_Dlon_S,df_Dlat_P,df_Dlat_S
    END TYPE elementmatrice
    ! -------                                                  --------    .
    integer(KIND=wi), intent(in) :: nbtps                                  ! nombre de données de temps P et S confondues
    type(dataone), intent(inout) :: datatps(nbtps)                            ! données
    type(parametre), intent(inout) :: param_init                           ! paramètres
    logical, intent(out) :: con                                            ! "true" si convergence
    integer(KIND=wi), intent(in) :: nk                                     ! nombre de divergences précendantes
    type (amoho_centroid), intent (in) :: acentroid
    character (LEN=1), intent(in) :: mod
    real(KIND=wr), intent(out) :: onemisfit
    logical, intent(in), optional :: chut
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k
    integer(KIND=wi) :: nbdonnees                                          ! nombre de données de temps P et S non confondues
    ! -------                                                  --------    . 
    real(KIND=wr), dimension(:,:), allocatable :: mat_A,mat_A_t,invAt_A_At ! matrices cf : SRC/SUB_Geiger/DOC/DOC_Geiger.pdf
    real(KIND=wr), dimension(4,4) :: At_A,invAt_A                          ! matrices cf : SRC/SUB_Geiger/DOC/DOC_Geiger.pdf
    real(KIND=wr), dimension(4) :: x_h, delta_x                            ! matrices cf : SRC/SUB_Geiger/DOC/DOC_Geiger.pdf
    real(KIND=wr), dimension(:), allocatable :: gamma                      ! matrice différence de temps théoriques vs. réel
    ! -------                                                  --------    . 
    type(elementmatrice) :: mat_el                                         ! matrice avec les dérivées partielles
    ! -------                                                  --------    .
    real(KIND=wr) :: sumdelta_x                                            ! précision à atteindre
    type(dataone) :: one_data
    integer(KIND=wi) :: dim                                                ! dimension, ici dim = 4
    character(LEN=76) :: prog                                              ! chaîne pour affichage
    logical :: critique, rechut
    real(KIND=wr) :: xmin=1000.0_wr, xmax=1500.0_wr                        ! diamètres de pondération, mais ici on s'en fiche ...
    real(KIND=wr) :: pdfmoho
    ! -----------------------------------------------------------------    .
    if (mod=='I') then
      call tempsTheoDirectone(param_init,datatps(1),critique,acentroid)  ! calcul le temps théorique et la difference théoriques et réels
      pdfmoho=param_init%Zmoho
    elseif ((mod=='A').or.(mod=='S').or.(mod=='C').or.(mod=='F')) then
      call tempsTheoDirectone_AUTRE(param_init,datatps(1),pdfmoho,critique,mod) ! calcul le temps théorique et la difference théoriques et réels
    else
      write(*,*)'problème dans geiger : modèle inexistant !'
      stop
    endif
    ! -----------------------------------------------------------------    .
    con=.true.                                                             ! par défaut
    ! -----------------------------------------------------------------    .
    nbdonnees=0                                                            ! taille de la matrice A (nombre de données directes et réfractés, à la fois P et S)
    do i=1,nbtps
      if (datatps(i)%coefP.ne.4) then                                      ! ne prend pas en compte les plus mauvaises données
        nbdonnees=nbdonnees+1
      endif
      if (datatps(i)%andS.eq.'S') then                                     ! si ondes S existe
        if (datatps(i)%coefS.ne.4) then                                    ! ne prend pas en compte les plus mauvaises données
          nbdonnees=nbdonnees+1
        endif
      endif
    enddo
    ! -------                                                  --------    .
    allocate(mat_A(nbdonnees,4),mat_A_t(4,nbdonnees),invAt_A_At(4,nbdonnees),gamma(nbdonnees))
    ! -----------------------------------------------------------------    .
    ! ------- solution apriori initiale                        --------    .
    x_h(1) = param_init%Lon
    x_h(2) = param_init%Lat
    if ((param_init%Zhypo.gt.0.0_wr).and.( param_init%Zhypo.lt.(pdfmoho-0.1_wr))) then
      x_h(3) = param_init%Zhypo                                            ! séisme entre la surface (ici 0.0 masl) et le moho
    else
      x_h(3) = pdfmoho/2.0_wr
    endif
    x_h(4) = param_init%Tzero%sec
    ! -------                                                  --------    .
    ! initialisation
    sumdelta_x=1.0e9_wr
    k=0
    ! -------                                                  --------    .
    ! processus itératif
    do while(con.and.(sumdelta_x.gt.1.e-4_wr))                             ! itération en cours
      gamma=0.0_wr
      k = k + 1
      rechut=.true.                                                        ! écriture du nombre d'itération
      if(present(chut)) then
        if (chut) rechut=.false.
      endif
      if (rechut) then
        if(nk.eq.1) then                                                   ! premier essai
          write(prog,'(a51,i3)')" Geiger méthode - nb itération  :                ",k
        elseif(nk.lt.4) then                                               ! nk essai
          write(prog,'(a51,i3,a2,i2.2,a)')" Geiger méthode - nb itération  :                ", &
          k," (",nk-1," fois divergente)"
        else
          write(prog,'(a51,i3,a)')" Geiger méthode - nb itération  :                ", &
          k," (non convergente)"
        endif
        write(*,'(a,a,$)') prog(1:76),char(13)                             ! affichage dynamique
      endif
      ! ---------------------------------------------------------------    .
      ! ------- remplissage de la matrice A                    --------    . matrice des dérivées partielles
      nbdonnees=0
      do i=1,nbtps
      ! -------                                                --------    .
        param_init%Tzero%sec=x_h(4)
        call basetime(param_init%Tzero)                                    ! reste en base 60/12/365 ...
        call derivpart(mat_el,datatps(i),param_init,acentroid,mod)         ! calcul les dérivées partielles
        if (datatps(i)%coefP.ne.4) then                                    ! ne prend pas en compte les plus mauvaises données
          nbdonnees=nbdonnees+1
          mat_A(nbdonnees,1)=mat_el%df_Dlon_P
          mat_A(nbdonnees,2)=mat_el%df_Dlat_P
          mat_A(nbdonnees,3)=mat_el%df_Dhypo_P
          mat_A(nbdonnees,4)=mat_el%df_Dto
        endif
        if (datatps(i)%andS.eq.'S') then                                   ! si ondes S existe
          if (datatps(i)%coefS.ne.4) then                                  ! ne prend pas en compte les plus mauvaises données
            nbdonnees=nbdonnees+1
            mat_A(nbdonnees,1)=mat_el%df_Dlon_S
            mat_A(nbdonnees,2)=mat_el%df_Dlat_S
            mat_A(nbdonnees,3)=mat_el%df_Dhypo_S
            mat_A(nbdonnees,4)=mat_el%df_Dto
          endif
        endif
      enddo
      ! ---------------------------------------------------------------    .
      ! ------- remplissage de la matrice A_t, transposée de A --------    .
      do j=1,nbdonnees
        do i=1,4
          mat_A_t(i,j) = mat_A(j,i)
        enddo
      enddo
      ! ---------------------------------------------------------------    .
      ! ------- remplissage de la matrice At_A, produit de At et A ----    .
      At_A = matmul(mat_A_t,mat_A)
      ! ---------------------------------------------------------------    .
      ! ------- inversion de la matrice At_A                   --------    .
      dim=4
      call inverse(At_A,invAt_A,dim)
      ! ---------------------------------------------------------------    .
      ! -------  remplissage de la matrice produit de invAt_A et A_t --    .
      invAt_A_At = matmul(invAt_A,mat_A_t)
      ! ---------------------------------------------------------------    .
      ! ------- remplissage de la matrice gamma                --------    . difference entre temps théoriques et réels
      nbdonnees=0
      do i=1,nbtps
        ! -------                                              --------    .
        one_data=datatps(i)
        ! -------                                              --------    .
        if (mod=='I') then
          call tempsTheoDirectone(param_init,one_data,critique,acentroid)  ! calcul le temps théorique et la difference théoriques et réels
        elseif ((mod=='A').or.(mod=='S').or.(mod=='F').or.(mod=='C')) then
          call tempsTheoDirectone_AUTRE(param_init,one_data,pdfmoho,critique,mod) ! calcul le temps théorique et la difference théoriques et réels
        else
          write(*,*)'problème dans geiger : modèle inexistant !'
          stop
        endif
        ! -------                                              --------    .
        if (datatps(i)%coefP.ne.4) then                                    ! ne prend pas en compte les plus mauvaises données
          nbdonnees=nbdonnees+1
          gamma(nbdonnees)=one_data%dTP
        endif
        ! -------                                              --------    .
        if (datatps(i)%andS.eq.'S') then                                   ! si ondes S existe
          if (datatps(i)%coefS.ne.4) then                                  ! ne prend pas en compte les plus mauvaises données
            nbdonnees=nbdonnees+1
            gamma(nbdonnees)=one_data%dTS
          endif
        endif
        ! -------                                              --------    .
        datatps(i)=one_data
        ! -------                                              --------    .
      enddo
      ! ---------------------------------------------------------------    .
      call compute_misfitone (nbtps,datatps,onemisfit,xmin,xmax,CorH='C')
      onemisfit=onemisfit/real(nbtps,wr)  ! or : onemisfit=sumdelta_x
      ! ------- remplissage de la matrice delta_x              --------    .
      delta_x = matmul(invAt_A_At,gamma)
      ! ------- itération suivante                             --------    .
      x_h=x_h+delta_x
      param_init%Lon = x_h(1)
      param_init%Lat = x_h(2)
      ! -------                                                --------    . séisme entre la surface et le moho
      if ((x_h(3).gt.0.0_wr).and.(x_h(3).lt.(pdfmoho-0.1_wr))) then
        param_init%Zhypo = x_h(3)
      end if
      ! -------                                                --------    .
      param_init%Tzero%sec = x_h(4)
      call basetime(param_init%Tzero)                                      ! reste en base 60/12/365 ...
      ! sumdelta_x = abs(delta_x(1))+abs(delta_x(2))+abs(delta_x(3))+abs(delta_x(4)) ! résidus
      sumdelta_x = abs(delta_x(1))+abs(delta_x(2))+abs(delta_x(4)) ! résidus sans pdf ...
      ! ------- divergence                                     --------    .
      if ((k.gt.25).or.(IsNaN(sumdelta_x)).or.(sumdelta_x.gt.1.e5_wr) ) then
        con=.false.                                                        ! c'est pas grave, on repart pour un autre modele de terre ...
        !write(*,*)'divergence après',k,'itérations'
        sumdelta_x=1.0e9_wr
        onemisfit=1.e9_wr
      endif
      ! -------                                                --------    .
    enddo                                                                  ! fin processus itératif
    ! -------                                                  --------    .
    deallocate (mat_A,mat_A_t,invAt_A_At,gamma)
    ! -----------------------------------------------------------------    .


    ! -----------------------------------------------------------------    .

      CONTAINS !

      subroutine inverse(a_ori,c,n)
        ! -------                                              --------    .
        ! Inverse matrix, Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        ! -------------------------------------------------------------    .
        ! source : http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
        ! -------------------------------------------------------------    .
        implicit none
        ! -------                                              --------    .
        real(KIND=wr), intent (in) :: a_ori(n,n) ! array of coefficients for matrix A
        real(KIND=wr) :: a(n,n)
        real(KIND=wr), intent (out) :: c(n,n)    ! inverse matrix of A
        integer(KIND=wi), intent (in) :: n       ! dimension
        real(KIND=wr) :: L(n,n), U(n,n), b(n), d(n), x(n)
        real(KIND=wr) :: coeff
        integer(KIND=wi) i, j, k
        ! -------                                              --------    .
        a=a_ori
        ! -------                                              --------    .
        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0_wr
        U=0.0_wr
        b=0.0_wr
        do k=1, n-1                                                        ! step 1: forward elimination
          do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
              a(i,j) = a(i,j)-coeff*a(k,j)
            end do
          end do
        end do
        ! Step 2: prepare L and U matrices
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0_wr
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do
        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0_wr
          d(1) = b(1)
          ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
          ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
          x(i) = x(i)/u(i,i)
          end do
          ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            c(i,k) = x(i)
          end do
          b(k)=0.0_wr
        end do
        ! -------------------------------------------------------------    .
      end subroutine inverse

        ! -------------------------------------------------------------    .

      subroutine derivpart(mat_el,one_data,p_i,acentroid,mod)
        ! -------                                              --------    .mh
        ! calcul élements de la matrice A en dérivant les lois de temps théoriques des arrivées des ondes
        ! dérivations non analytiques : dérivée par différences finies centrées d’ordre 2
        ! -------                                              --------    .
        implicit none
        ! -------                                              --------    .
        type(parametre), intent (in) :: p_i                                ! paramètres de l'inversion
        type(dataone), intent (in) :: one_data                             ! pour une donnée
        type(elementmatrice), intent (out) :: mat_el                       ! élements de la matrice (dérivées partielles)
        type (amoho_centroid), intent (in) :: acentroid
        character (LEN=1), intent(in) :: mod
        ! -------                                              --------    .
        type(dataone) :: donnees                                           ! pour cette donnée
        type(parametre) :: param                                           ! pour ce jeu de paramètres
        ! -------                                              --------    .
        real(KIND=wr), parameter :: h=0.00001_wr                           ! précision
        real(KIND=wr) :: d1p,d2p,d1s,d2s
        real(KIND=wr) :: pdfmoho
        logical :: critique
        ! -------------------------------------------------------------    .
        donnees=one_data
        param=p_i
        ! -------------------------------------------------------------    . Tzéro (temps derivé par le temps) :
        mat_el%df_Dto = 1.0_wr
        ! -------                                              --------    . Zhypo :
        param%Zhypo = param%Zhypo - h
        if (mod=='I') then
          call tempsTheoDirectone(param,donnees,critique,acentroid)
        elseif ((mod=='A').or.(mod=='S').or.(mod=='F').or.(mod=='C')) then
          call tempsTheoDirectone_AUTRE(param,donnees,pdfmoho,critique,mod)
        else
          write(*,*)'problème dans derivpart : modèle inexistant !'
          stop
        endif
        d1p = donnees%tpsTh%secP
        d1s = donnees%tpsTh%secS
        param%Zhypo = param%Zhypo + 2.0_wr*h
        if (mod=='I') then
          call tempsTheoDirectone(param,donnees,critique,acentroid)
        else
          call tempsTheoDirectone_AUTRE(param,donnees,pdfmoho,critique,mod)
        endif
        d2p = donnees%tpsTh%secP
        d2s = donnees%tpsTh%secS
        mat_el%df_Dhypo_P = (d2p - d1p)/(2.0_wr*h)
        mat_el%df_Dhypo_S = (d2s - d1s)/(2.0_wr*h)
        ! -------                                              --------    . Lon :
        param=p_i
        param%lon = param%lon - h
        if (mod=='I') then
          call tempsTheoDirectone(param,donnees,critique,acentroid)
        else
          call tempsTheoDirectone_AUTRE(param,donnees,pdfmoho,critique,mod)
        endif
        d1p = donnees%tpsTh%secP
        d1s = donnees%tpsTh%secS
        param%lon = param%lon + 2.0_wr*h
        if (mod=='I') then
          call tempsTheoDirectone(param,donnees,critique,acentroid)
        else
          call tempsTheoDirectone_AUTRE(param,donnees,pdfmoho,critique,mod)
        endif
        d2p = donnees%tpsTh%secP
        d2s = donnees%tpsTh%secS
        mat_el%df_Dlon_P = (d2p - d1p)/(2.0_wr*h)
        mat_el%df_Dlon_S = (d2s - d1s)/(2.0_wr*h)
        ! -------                                              --------    . Lat :
        param=p_i
        param%lat = param%lat - h
        if (mod=='I') then
          call tempsTheoDirectone(param,donnees,critique,acentroid)
        else
          call tempsTheoDirectone_AUTRE(param,donnees,pdfmoho,critique,mod)
        endif
        d1p = donnees%tpsTh%secP
        d1s = donnees%tpsTh%secS
        param%lat = param%lat + 2.0_wr*h
        if (mod=='I') then
          call tempsTheoDirectone(param,donnees,critique,acentroid)
        else
          call tempsTheoDirectone_AUTRE(param,donnees,pdfmoho,critique,mod)
        endif
        d2p = donnees%tpsTh%secP
        d2s = donnees%tpsTh%secS
        mat_el%df_Dlat_P = (d2p - d1p)/(2.0_wr*h)
        mat_el%df_Dlat_S = (d2s - d1s)/(2.0_wr*h)
        ! -------------------------------------------------------------    .
      end subroutine derivpart

  end subroutine geiger

END MODULE invGEIGER



! *********************************************************************    .
! *********************************************************************    .


