! Librairie de subroutines concernant le calcul de la fonction coût
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE sub_misfit

    use modparam

    implicit none

    private

    public  :: ponderation
    public  :: compute_misfitone, compute_misfit

CONTAINS

    ! -----------------------------------------------------------------    .

  subroutine ponderation (nbtps,datatps,xmin,xmax,w)
    ! -------                                                  --------    .mh
    ! calcul de la ponderation des données (distance et qualité)
    ! -------                                                  --------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    type(dataone), intent (inout) :: datatps(nbtps)                        ! données de temps
    integer(KIND=wi), intent (in) :: nbtps                                 ! nb de données de temps P et S confondues
    real(KIND=wr), intent (in) :: xmin,xmax                                ! diamètres cercles pondération
    type(pond), intent (out) :: w                                          ! pondération
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    real(KIND=wr) :: one_coef_P, one_coef_S, two_coef                      ! coeficients
    real(KIND=wr) :: a, b                                                  ! droite de régression
    ! -----------------------------------------------------------------    . initialistaion
    w%S_Pg = 0.0_wr                                                        ! pondération
    w%S_Sg = 0.0_wr
    w%S_Pn = 0.0_wr
    w%S_Sn = 0.0_wr
    w%nPg = 0                                                              ! nb de données Pg
    w%nPn = 0                                                              ! nb de données Pn
    w%nSg = 0                                                              ! nb de données Sg
    w%nSn = 0                                                              ! nb de données Sn
    ! -------                                                  --------    .
    do i=1,nbtps
      ! ---------------------------------------------------------------    .
      ! ------- coeficient 1 : distance épicentrale            --------    .
      if (FLAGcercles) then
        if (datatps(i)%depi.lt.xmin/2.0_wr) then
          two_coef = 1.0_wr                                                ! distance faible -> gros coef [1]
        elseif(datatps(i)%depi.gt.xmax/2.0_wr) then                        ! distance forte -> petit coef [0.1]
          two_coef = 0.000001_wr                                           ! très faible mais non nulle -> sinon fonction cout -> NaN (parfois)
        else
          b = (1.0_wr-0.1_wr*(xmin/xmax))/(1.0_wr-(xmin/xmax))
          a = (0.1_wr-b)/(xmax/2.0_wr)
          two_coef = a * datatps(i)%depi + b                               ! entre deux [0.1(loin);1(près)]
        endif
      else
          two_coef = 1.0_wr
      endif
      ! ---------------------------------------------------------------    .
      ! ------- coeficient 2 : qualité des données             --------    . transforme [0(bon);4(mauvais)] -> [0(mauvais);1(bon)]
      one_coef_P = -0.25_wr*real(datatps(i)%coefP,wr) + 1.0_wr
      one_coef_S = -0.25_wr*real(datatps(i)%coefS,wr) + 1.0_wr
      ! ---------------------------------------------------------------    .
      ! -------                                                --------    . ONDE P
      if(datatps(i)%typeonde.eq.'G') then                                  ! si onde directe compressive
        w%nPg = w%nPg+1
        datatps(i)%wp = (one_coef_P * two_coef)
        w%S_Pg = w%S_Pg + datatps(i)%wp                                    ! somme de la pondération pour Pg
      elseif(datatps(i)%typeonde.eq.'N') then                              ! si onde réfractée compressive
        w%nPn = w%nPn+1
        datatps(i)%wp = (one_coef_P * two_coef)
        w%S_Pn = w%S_Pn + datatps(i)%wp                                    ! somme de la pondération pour Pn
      else
        write(*,*)'problème 1 dans ponderation : onde compressive ni directe ni réfractée '
        stop
      endif
      ! -------                                                --------    . ONDE S
      if(datatps(i)%andS.eq.'S') then                                      ! si S existe
        if(datatps(i)%typeonde.eq.'G') then                                ! si onde directe cisaillante
          w%nSg = w%nSg+1
          datatps(i)%ws = (one_coef_S * two_coef)
          w%S_Sg = w%S_Sg + datatps(i)%ws                                  ! somme de la pondération pour Sg
        elseif(datatps(i)%typeonde.eq.'N') then                            ! si onde réfractée cisaillante
          w%nSn = w%nSn+1
          datatps(i)%ws = (one_coef_S * two_coef)
          w%S_Sn = w%S_Sn + datatps(i)%ws                                  ! somme de la pondérationpour Sn
        else
          write(*,*)'problème 2 dans ponderation : onde cisaillante ni directe ni réfractée '
          stop
        endif
      endif
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine ponderation

    ! -----------------------------------------------------------------    .

  subroutine compute_misfit (nbtps,D,misfit,xmin,xmax,CorH,div)
    ! -------                                                  --------    .mh
    ! calcul de la fonction coût, pour tous les séismes
    ! -------                                                  --------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    type(dataall), intent (inout) :: D(nbseismes)                          ! données de temps
    integer(KIND=wi), intent (in) :: nbtps(nbseismes)                      ! nb de données de temps P et S confondues
    real(KIND=wr), intent (inout) :: xmin(nbseismes),xmax(nbseismes)       ! diamètres cercles pondération
    real(KIND=wr), intent (out) :: misfit                                  ! valeur de la fonction coût
    character (LEN=1), intent (in) :: CorH                                 ! cold or Hot runs
    logical, intent (inout), optional :: div
    ! -------                                                  --------    .
    real(KIND=wr) :: misfit_one
    integer(KIND=wi) :: i,N
    ! -----------------------------------------------------------------    .
    misfit = 0.0_wr
    do i=1, nbseismes
      if (present(div)) then
        ! div : chaine divergente ? -> stop en coldruns
        call compute_misfitone (nbtps(i),D(i)%datatps,misfit_one,xmin(i),xmax(i),CorH,div)
      else
        ! div : chaine divergente ? -> pas de stop en hotruns
        call compute_misfitone (nbtps(i),D(i)%datatps,misfit_one,xmin(i),xmax(i),CorH)
      endif
      misfit = misfit + misfit_one
    enddo
    ! -------                                                  --------    . nombre de données
    N=0
    do i=1, nbseismes
      N = N + nbtps(i)
    enddo
    misfit = misfit/real(N,wr)
    ! -----------------------------------------------------------------    .
  end subroutine compute_misfit

    ! -----------------------------------------------------------------    .

  subroutine compute_misfitone (nbtps,datatps,misfit,xmin,xmax,CorH,div)
    ! -------                                                  --------    .mh
    ! calcul de la fonction coût, par séisme
    ! -------                                                  --------    .
    use typetemps
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real(KIND=wr), parameter :: alpha=3.0_wr
    ! -------                                                  --------    .
    type(dataone), intent (inout) :: datatps(nbtps)                        ! données de temps
    integer(KIND=wi), intent (in) :: nbtps                                 ! nb de données de temps P et S confondues
    real(KIND=wr), intent (inout) :: xmin,xmax                             ! diamètres cercles pondération
    real(KIND=wr), intent (out) :: misfit                                  ! valeur de la fonction coût
    character (LEN=1), intent (in) :: CorH                                 ! cold or Hot runs
    ! -------                                                  --------    .
    logical, intent (inout), optional :: div
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    real(KIND=wr) :: mis_Pg, mis_Sg, mis_Pn, mis_Sn
    real(KIND=wr) :: moy
    type(pond) :: w                                                        ! coeficients
    logical :: FLAGnormale
    ! -------                                                  --------    .
    integer(KIND=wi), save :: nbdoublement=0
    ! -----------------------------------------------------------------    . si divergent :
123 continue
    ! -----------------------------------------------------------------    . initialistaion
    mis_Pg = 0.0_wr                                                        ! misfits différents par ondes
    mis_Pn = 0.0_wr
    mis_Sg = 0.0_wr
    mis_Sn = 0.0_wr
    ! -------                                                  --------    .
    call ponderation (nbtps,datatps,xmin,xmax,w)                           ! calcul de la pondération (qualité des données / distance)
    ! -------                                                  --------    .
    do i=1,nbtps
      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .
      FLAGnormale=.true.
      ! ---------------------------------------------------------------    .
      if(FLAGnormale) then
        ! -------------------------------------------------------------    .
        ! les données ont des distributions d'incertitudes de types normales symetriques
        ! norme L2
        ! pondérées
        ! -------------------------------------------------------------    .
        ! -------                                              --------    . ONDE P
        if(datatps(i)%typeonde.eq.'G') then                                ! si onde directe compressive

          mis_Pg = mis_Pg + datatps(i)%wp/2.0_wr * (datatps(i)%dTP/datatps(i)%sigP)**2.0_wr

        elseif(datatps(i)%typeonde.eq.'N') then                            ! si onde réfractée compressive

          mis_Pn = mis_Pn + datatps(i)%wp/2.0_wr * (datatps(i)%dTP/datatps(i)%sigP)**2.0_wr

        else
          write(*,*)'problème dans compute_misfit : onde compressive ni directe ni réfractée '
          stop
        endif
        if(IsNaN(mis_Pg).or.IsNaN(mis_Pn)) then
          write(*,*)'problème dans compute_misfit : misfit_P = NaN'
          stop
        endif
        ! -------------------------------------------------------------    .
        ! -------                                              --------    . ONDE S
        if(datatps(i)%andS.eq.'S') then                                    ! si S existe
          if(datatps(i)%typeonde.eq.'G') then                              ! si onde directe cisaillante

            mis_Sg = mis_Sg + datatps(i)%ws/2.0_wr * (datatps(i)%dTS/datatps(i)%sigS)**2.0_wr

          elseif(datatps(i)%typeonde.eq.'N') then                          ! si onde réfractée cisaillante

            mis_Sn = mis_Sn + datatps(i)%ws/2.0_wr * (datatps(i)%dTS/datatps(i)%sigS)**2.0_wr

          else
            write(*,*)'problème dans compute_misfit : onde cisaillante ni directe ni réfractée '
            stop
          endif
          if(IsNaN(mis_Sg).or.IsNaN(mis_Sn)) then
            write(*,*)'problème dans compute_misfit : misfit_S = NaN'
            stop
          endif
        endif
        ! -------------------------------------------------------------    .
        ! -------------------------------------------------------------    .
      else
        ! -------------------------------------------------------------    .
        ! les données ont des distributions d'incertitudes de types normales assymetriques
        ! norme L2
        ! pondérées
        ! -------------------------------------------------------------    .
        write(*,*)'problème dans compute_misfit : incertitudes de types LOG-normales non codée'
        stop
        ! -------------------------------------------------------------    .
        moy=100.0_wr
        ! -------                                              --------    . ONDE P
        if(datatps(i)%typeonde.eq.'G') then                                ! si onde directe compressive

          mis_Pg = mis_Pg + datatps(i)%wp*( 1.0_wr/2.0_wr * (datatps(i)%dTP/datatps(i)%sigP)**2.0_wr + &
            log((1._wr+erf((alpha*datatps(i)%dTP)/(datatps(i)%sigP*sqrt(2.0_wr))))/(datatps(i)%sigP*sqrt(2.0_wr*pi))))

        elseif(datatps(i)%typeonde.eq.'N') then                            ! si onde réfractée compressive

          mis_Pn = mis_Pn + datatps(i)%wp*( 1.0_wr/2.0_wr * (datatps(i)%dTP/datatps(i)%sigP)**2.0_wr + &
            log((1._wr+erf((alpha*datatps(i)%dTP)/(datatps(i)%sigP*sqrt(2.0_wr))))/(datatps(i)%sigP*sqrt(2.0_wr*pi))))

        else
          write(*,*)'problème dans compute_misfit : onde compressive ni directe ni réfractée '
          stop
        endif
        if(IsNaN(mis_Pg).or.IsNaN(mis_Pn)) then
          write(*,*)'problème dans compute_misfit : misfit_P = NaN'
          stop
        endif
        ! -------------------------------------------------------------    .
        ! -------                                              --------    . ONDE S
        if(datatps(i)%andS.eq.'S') then                                    ! si S existe
          if(datatps(i)%typeonde.eq.'G') then                              ! si onde directe cisaillante

            mis_Sg = mis_Sg + datatps(i)%ws*( 1.0_wr/2.0_wr * (datatps(i)%dTS/datatps(i)%sigS)**2.0_wr + &
              log((1._wr+erf((alpha*datatps(i)%dTS)/(datatps(i)%sigS*sqrt(2.0_wr))))/(datatps(i)%sigS*sqrt(2.0_wr*pi))))

          elseif(datatps(i)%typeonde.eq.'N') then                          ! si onde réfractée cisaillante

            mis_Sn = mis_Sn +  datatps(i)%ws*( 1.0_wr/2.0_wr * (datatps(i)%dTS/datatps(i)%sigS)**2.0_wr + &
              log((1._wr+erf((alpha*datatps(i)%dTS)/(datatps(i)%sigS*sqrt(2.0_wr))))/(datatps(i)%sigs*sqrt(2.0_wr*pi))))

          else
            write(*,*)'problème dans compute_misfit : onde cisaillante ni directe ni réfractée '
            stop
          endif
          if(IsNaN(mis_Sg).or.IsNaN(mis_Sn)) then
            write(*,*)'problème dans compute_misfit : misfit_S = NaN'
            stop
          endif
        endif
        ! -------------------------------------------------------------    .
        ! -------------------------------------------------------------    .
      endif

      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .
    enddo
    ! -------                                                  --------    . calcul de la fonction coût
    misfit = mis_Pg + mis_Pn + mis_Sg + mis_Sn
    ! -----------------------------------------------------------------    .

    ! -----------------------------------------------------------------    .
    ! si divergent, séisme loin de son épicentre : w tend vers 0 et misfit vers NaN
    ! sauf si w, jamais nul car : two_coef = 0.000001_wr
    ! -------                                                  --------    .
    if(IsNaN(misfit).or.(misfit.gt.1.e99_wr)) then                            ! w tends vers 0 (car la loc. est divergente pour ce séisme)
      ! -------                                                --------    . COLDRUNS
      if (CorH.eq."C") then
        if (present(div)) then                                             ! STOP LA CHAINE
          div=.true.
        else
          write(*,*)'problème dans compute_misfit avec "C" : misfit = NaN '
          stop
      endif
      ! -------                                                --------    . HOTRUNS
      elseif (CorH.eq."H") then
        if (present(div)) then                                             ! STOP LA CHAINE
          div=.true.
        else
          if(nbdoublement.lt.nbseismes) then
            write(*,*)'attention dans compute_misfit avec "H" : doublement des cercles de pondération (xmax et xmin) '
            nbdoublement=nbdoublement+1
            xmin=xmin*2._wr
            xmax=xmax*2._wr
          else
            write(*,*)'PB :'
            write(*,*)xmin,xmax
            write(*,*)'problème dans compute_misfit avec "H" : misfit = NaN, Pg ->',mis_Pg,w%S_Pg
            write(*,*)'problème dans compute_misfit avec "H" : misfit = NaN, Pn ->',mis_Pn,w%S_Pn
            write(*,*)'problème dans compute_misfit avec "H" : misfit = NaN, Sg ->',mis_Sg,w%S_Sg
            write(*,*)'problème dans compute_misfit avec "H" : misfit = NaN, Sn ->',mis_Sn,w%S_Sn
            do i=1,nbtps
              write(*,*)datatps(i)
            enddo
            write(*,*)w
            stop
          endif
          go to 123
        endif
      endif
    endif
    ! -------                                                  --------    .
    if(xmax.lt.xmin) then
      write(*,*)'problème dans compute_misfit : xmax < xmin :: ', xmax," ? < ? ",xmin
      stop
    endif
    ! -----------------------------------------------------------------    .
  end subroutine compute_misfitone

END MODULE sub_misfit



! *********************************************************************    .
! *********************************************************************    .


