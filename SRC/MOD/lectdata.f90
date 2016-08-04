! Librairie de subroutines permettant la lecture des données
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE datalecture

    use modparam

    implicit none

    private

    public  :: lectnbdata
    public  :: lectdata
    public  :: mksynth,mksynthallsta
    public  :: cerclespond
    public  :: catalogue
    public  :: initR

CONTAINS

! ---------------------------------------------------------------------    .

  subroutine lectnbdata(nbsta,nbtps)
    ! -------                                                  --------    .mh
    ! lecture du nombre de données (par séisme) et du nombre de stations connues
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (out) :: nbsta,nbtps(nbseismes)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,ok
    character(len=35) :: seismeNom(nbseismes)
    ! -----------------------------------------------------------------    . nombre de station connues
    nbsta = 0
    ok = 0
    open(500, FILE = 'DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      open(500, FILE = 'sta.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans lectnbdata : le fichier DATA/sta.d n''existe pas '
        write(*,*)'problème dans lectnbdata : le fichier sta.d n''existe pas '
        stop
      endif
    endif
    do while(ok .eq. 0)                                                    ! boucle pour compter le nombre de lignes du fichier
      read(500,*,iostat = ok)
      if (ok .eq. 0) nbsta = nbsta + 1
    end do
    close(500)
    ! -----------------------------------------------------------------    . nombre de séismes
    call nomseismefichier(seismeNom)                                       ! lit les différents noms de phases list
    ! -------                                                  --------    . nombre de données par séismes
    do i=1,nbseismes
      nbtps(i) = 0
      ok = 0
      open(501, FILE = 'DATA/'//seismeNom(i),status='old',iostat = ok)
      if (ok .ne. 0) then
        open(501, FILE = seismeNom(i),status='old',iostat = ok)
        if (ok .ne. 0) then
          write(*,*)'problème dans lectnbdata : le fichier DATA/',seismeNom(i),' n''existe pas '
          write(*,*)'problème dans lectnbdata : le fichier ',seismeNom(i),' n''existe pas '
          stop
        endif
      endif
      read(501,*,iostat = ok)
      do while(ok .eq. 0)                                                  ! boucle pour compter le nombre de lignes du fichier
        read(501,*,iostat = ok)
        if (ok .eq. 0) nbtps(i) = nbtps(i) + 1
      end do
      close(501)
      if (nbtps(i).lt.3) then
        write(*,*)'problème dans lectnbdata : le fichier DATA/',seismeNom(i),' contient trop peu de données : nb lignes < 4 '
        stop
      endif
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine lectnbdata

    ! -----------------------------------------------------------------    .

  subroutine lectdata(nbsta,nbtps,D)
    ! -------                                                  --------    .mh
    ! lecture des données dans DATA/xxxx.xx.xx.xx.xx.d et DATA/sta-x.d
    ! modification des données, si résidus (station correction) aux stations disponibles dans DATA/sta.d
    ! -------                                                  --------    .
    use typetemps, only : dataall, stations
    use time, only : JDATE, GDATE, basetime
    use tri, only : tridata
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: nbsta,nbtps(nbseismes)
    type(dataall), intent (inout) :: D(nbseismes)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k,ok,count
    type(stations) :: datasta(nbsta)
    character(len=35) :: seismeNom(nbseismes)
    ! -----------------------------------------------------------------    . lecture des données stations
    open(503, FILE = 'DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans lectdata : le fichier DATA/sta.d n''existe pas '
      stop
    endif
    do i=1,nbsta
      read(503,*,iostat = ok) datasta(i) ! nouveau -> (i)
    enddo
    close(503)
    ! -----------------------------------------------------------------    . nom des fichier par séismes
    call nomseismefichier(seismeNom)                                       ! lit les différents noms de phases list
    ! -------                                                  --------    . données par séismes
    a_event : do k=1,nbseismes
      !write(*,*)'DATA/'//seismeNom(k)
      open(504, FILE = 'DATA/'//seismeNom(k),status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans lectdata : le fichier DATA',seismeNom(i),' n''existe pas '
        stop
      endif
      read(504,*,iostat = ok)
      a_sta : do i=1,nbtps(k)                                              ! lecture des données temps
        read(504,1000,iostat = ok)D(k)%datatps(i)%sta%staname,D(k)%datatps(i)%typeonde, &
        D(k)%datatps(i)%coefP,D(k)%datatps(i)%tpsR%date%year,D(k)%datatps(i)%tpsR%date%month, &
        D(k)%datatps(i)%tpsR%date%day,D(k)%datatps(i)%tpsR%date%hour,D(k)%datatps(i)%tpsR%date%min, &
        D(k)%datatps(i)%tpsR%secP,D(k)%datatps(i)%tpsR%secS,D(k)%datatps(i)%andS,D(k)%datatps(i)%coefS, &
        D(k)%datatps(i)%sigP,D(k)%datatps(i)%sigS

        ! -------                                              --------    . incertitudes en secondes sur les données
        if((D(k)%datatps(i)%sigP.lt.0.0001_wr).and.(D(k)%datatps(i)%sigP.gt.10._wr).and.(IsNaN(D(k)%datatps(i)%sigP))) then
          write(*,*)'problème dans lectdata : les incertitudes sur les données P n''existent pas ',D(k)%datatps(i)%sigP
          stop
        endif
        if ((D(k)%datatps(i)%andS.ne."S").and.(D(k)%datatps(i)%sigS.gt.10._wr) &
          .and.(D(k)%datatps(i)%sigS.lt.0.0001_wr).and.(IsNaN(D(k)%datatps(i)%sigS))) then
          write(*,*)'problème dans lectdata : les incertitudes sur les données S n''existent pas ',D(k)%datatps(i)%sigS
          stop
        endif
        ! -------                                              --------    . si pas d'ondes S
        if(D(k)%datatps(i)%andS.ne."S") then
          D(k)%datatps(i)%andS="X"
          D(k)%datatps(i)%tpsR%secS=0.0_wr
          D(k)%datatps(i)%coefS=4
        endif
        ! -------                                              --------    .
        D(k)%datatps(i)%tpsR%date%year=D(k)%datatps(i)%tpsR%date%year+2000
        ! -------                                              --------    . respect du decoupage des années en mois et jours avec prise en compte des années bisextiles
        call basetime(D(k)%datatps(i)%tpsR)
        call JDATE(D(k)%datatps(i)%tpsR%date%Jday,D(k)%datatps(i)%tpsR%date%year, &
          D(k)%datatps(i)%tpsR%date%month,D(k)%datatps(i)%tpsR%date%day)
        call GDATE (D(k)%datatps(i)%tpsR%date%Jday,D(k)%datatps(i)%tpsR%date%year, &
          D(k)%datatps(i)%tpsR%date%month,D(k)%datatps(i)%tpsR%date%day)
        call basetime(D(k)%datatps(i)%tpsR)                                ! respect du decoupage temps dans la base composite 60/12/365 ...
        ! -------                                              --------    .
        count=0
        do j=1,nbsta
          if(D(k)%datatps(i)%sta%staname.eq.datasta(j)%staname) then
            D(k)%datatps(i)%sta=datasta(j)                                 ! attribution d'une station
            count=count+1
          endif
        enddo
        if(count.gt.1) then                                                ! vérification de l'absence de doublons
          write(*,*)'problème dans lectdata : station ',D(k)%datatps(i)%sta%staname,' en double in file : DATA/sta.d'
          stop
        endif
        if(count.eq.0) then
          write(*,*)'problème dans lectdata : station ',D(k)%datatps(i)%sta%staname,' non répertoriée'
          stop
        endif
        ! -------                                              --------    .
        ! ajout des résidus aux station (DATA/sta.d)
        ! les données sont donc modifié en amont
        ! -------                                              --------    .
        if (D(k)%datatps(i)%typeonde.eq.'G') then                          ! ondes directes
          D(k)%datatps(i)%tpsR%secP=D(k)%datatps(i)%tpsR%secP-D(k)%datatps(i)%sta%res_Pg
          if(D(k)%datatps(i)%andS.eq.'S') D(k)%datatps(i)%tpsR%secS=D(k)%datatps(i)%tpsR%secS-D(k)%datatps(i)%sta%res_Sg
        elseif (D(k)%datatps(i)%typeonde.eq.'N') then                      ! ondes réfractées
          ! -------                                            --------    . ondes réfractées moins pondérées !!!!!!!!!!!!!!!
          !D(k)%datatps(i)%coefP= min(D(k)%datatps(i)%coefP+2,4)
          !D(k)%datatps(i)%coefS= min(D(k)%datatps(i)%coefS+2,4)
          ! -------                                            --------    . ondes réfractées moins pondérées !!!!!!!!!!!!!!!
          D(k)%datatps(i)%tpsR%secP=D(k)%datatps(i)%tpsR%secP-D(k)%datatps(i)%sta%res_Pn
          if(D(k)%datatps(i)%andS.eq.'S') D(k)%datatps(i)%tpsR%secS=D(k)%datatps(i)%tpsR%secS-D(k)%datatps(i)%sta%res_Sn
        else
          write(*,*)'problème dans lectdata : onde ni directe ni réfractée ... ? '
          write(*,*)D(k)%datatps(i)%sta%staname
          write(*,*)
          write(*,*)D(k)%datatps(i)
          write(*,*)
          stop
        endif
        call basetime(D(k)%datatps(i)%tpsR)
        call modifPY42(D(k)%datatps(i))                                    ! gestion de cas particuliers ...
        ! -------                                              --------    .
      enddo a_sta
      close(504)
      ! -------                                                --------    .
      call tridata(nbtps(k),D(k)%datatps)                                  ! tri des donnée selon un temps d'arrivée des ondes P croissant
    enddo a_event
    ! -------                                                  --------    . format de lecture des données
    1000 format(a4,1x,a1,1x,i1,1x,5i2.2,f6.3,5x,f7.3,1x,a1,i1.1,3x,f6.3,4x,f6.3)
    ! -----------------------------------------------------------------    .

    CONTAINS

      ! ---------------------------------------------------------------    .
      subroutine modifPY42(adata)
      ! ---------------------------------------------------------------    . mh
      ! modification des données, gestion de cas particuliers ...
      ! ici modification de localisation de la station PY42 de PyrOPE
      ! PY42A et PY42B
      ! ---------------------------------------------------------------    .
      use typetemps, only : dataone,date_secPS
      use time, only : difftime
      ! ---------------------------------------------------------------    .
      implicit none
      type(dataone), intent (inout) :: adata
      ! ---------------------------------------------------------------    .
      type(date_secPS) :: datespe
      real(KIND=wr) :: diff1,diff2
      real(KIND=wr) :: deltatime
      ! ---------------------------------------------------------------    .
      if (adata%sta%staname=='PY42') then
        ! -------                                              --------    .
        datespe%date%year=2012
        datespe%date%month=10
        datespe%date%day=16
        datespe%date%hour=12
        datespe%date%min=30
        call JDATE (datespe%date%Jday,datespe%date%year,datespe%date%month,datespe%date%day)
        call basetime(datespe)
        datespe%secP=0.0_wr
        datespe%secS=0.0_wr
        ! -------                                              --------    .
        call difftime(diff1,diff2,adata%tpsR,datespe)
        deltatime=real(60*60*24*2,wr)                                      ! 2 jours
        ! -------                                              --------    .
        if (diff1.gt.deltatime) then                                       ! PY42-B
          adata%sta%lon=-1.20256_wr
          adata%sta%lat=46.39462_wr
          adata%sta%alti=53.0_wr
        elseif (diff1.lt.(-deltatime)) then                                ! PY42-A
          adata%sta%lon=-1.20075_wr
          adata%sta%lat=46.41012_wr
          adata%sta%alti=50.0_wr
        else                                                               ! n'existe pas !
          adata%sta%lon=0.0_wr
          adata%sta%lat=0.0_wr
          adata%sta%alti=0.0_wr
          write(*,*)'problème dans modifPY42 (lectdata) : avec PY42 ',adata
          stop
        endif
        ! -------                                              --------    .
      endif
      ! ---------------------------------------------------------------    .
      end subroutine modifPY42
      ! ---------------------------------------------------------------    .

  end subroutine lectdata

    ! -----------------------------------------------------------------    .

  subroutine initR(D,R,maxiter,nbtps,nbsta,nbofSta)
    ! -----------------------------------------------------------------    .mh
    ! initialise le calcul des résidus (si FLAGresSTA=.true.), en allouant R
    ! pour chaque station, chaque itération on sauve le résidus 
    ! afin de définir des tendenses -> station correction
    ! sub inR et outR -> dans MOD/subparam.f90
    ! -----------------------------------------------------------------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    type(dataall), intent (in) :: D(nbseismes)                             ! données
    integer(KIND=wi), intent (in) :: maxiter, nbtps(nbseismes),nbsta
    integer(KIND=wi), intent (out) :: nbofSta                              ! nb de phases
    ! -------                                                  --------    .
    type(residus), dimension(:), allocatable, intent (out) :: R
    character (LEN=4) :: nomDesta(2*nbsta+2)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k,l
    logical :: test
    ! -----------------------------------------------------------------    . nb de station différente pour ondes directes et refractées
    do i=1,2*nbsta+2
      nomDesta(i)='1234'
    enddo
    ! -------                                                  --------    . stations ?
    nbofSta=0
    do i=1,nbseismes
      do j=1,nbtps(i)
        test=.true.
        do k=1,nbofSta+2
          if (D(i)%datatps(j)%sta%staname==nomDesta(k)) then               ! déjà présente ?
            test=.false.
          endif
        enddo
        if (test) then
          nbofSta=nbofSta+1
          nomDesta(nbofSta)=D(i)%datatps(j)%sta%staname
        endif
      enddo
    enddo
    ! -------                                                  --------    .
    if (nbsta.lt.nbofSta) then
      write(*,*)'problème dans initR : ',nbsta,'<',nbofSta
      stop
    endif
    ! -----------------------------------------------------------------    .
    allocate(R(nbofSta))
    ! -------                                                  --------    .
    do i=1,nbofSta
      R(i)%nbPg=-1000
      R(i)%nbSg=-1000
      R(i)%nbPn=-1000
      R(i)%nbSn=-1000
      R(i)%nbPgT=-1000
      R(i)%nbSgT=-1000
      R(i)%nbPnT=-1000
      R(i)%nbSnT=-1000
      R(i)%staname=nomDesta(i)
    enddo
    ! -------                                                  --------    . existe-il des Pg, Pn, Sg et Sn pour ces stations ?
    do k=1,nbofSta                                                         ! pour chaque station
      ! -------                                                --------    .
      l=0
      do i=1,nbseismes                                                     ! Pg ?
        do j=1,nbtps(i)
          if ((D(i)%datatps(j)%sta%staname==R(k)%staname).and.(D(i)%datatps(j)%typeonde=='G')) l=l+1
        enddo
      enddo
      if (l.gt.0) R(k)%nbPg=0
      if (l.gt.0) R(k)%nbPgT=0
      if ((l.gt.0).and.(.not.allocated(R(k)%resPg))) then
        allocate(R(k)%resPg(maxiter*l,2))
        R(k)%resPg(:,1)=0.0_wr ! résidu
        R(k)%resPg(:,2)=1.e9_wr ! fonction coût
      endif
      if (l.gt.nbseismes) then
        write(*,*)'problème dans initR Pg : trop de ok ',l,nbseismes
        stop
      endif
      ! -------                                                --------    .
      l=0
      do i=1,nbseismes                                                     ! Pn ?
        do j=1,nbtps(i)
          if ((D(i)%datatps(j)%sta%staname==R(k)%staname).and.(D(i)%datatps(j)%typeonde=='N')) l=l+1
        enddo
      enddo
      if (l.gt.0) R(k)%nbPn=0
      if (l.gt.0) R(k)%nbPnT=0
      if ((l.gt.0).and.(.not.allocated(R(k)%resPn))) then
        allocate(R(k)%resPn(maxiter*l,2))
        R(k)%resPn(:,1)=0.0_wr ! résidu
        R(k)%resPn(:,2)=1.e9_wr ! fonction coût
      endif
      if (l.gt.nbseismes) then
        write(*,*)'problème dans initR Pn : trop de ok ',l,nbseismes
        stop
      endif
      ! -------                                                --------    .
      l=0
      do i=1,nbseismes                                                     ! Sg ?
        do j=1,nbtps(i)
          if ((D(i)%datatps(j)%sta%staname==R(k)%staname).and.(D(i)%datatps(j)%typeonde=='G') &
            .and.(D(i)%datatps(j)%andS=='S')) l=l+1
        enddo
      enddo
      if (l.gt.0) R(k)%nbSg=0
      if (l.gt.0) R(k)%nbSgT=0
      if ((l.gt.0).and.(.not.allocated(R(k)%resSg))) then
        allocate(R(k)%resSg(maxiter*l,2))
        R(k)%resSg(:,1)=0.0_wr ! résidu
        R(k)%resSg(:,2)=1.e9_wr ! fonction coût
      endif
      if (l.gt.nbseismes) then
        write(*,*)'problème dans initR Sg : trop de ok ',l,nbseismes
        stop
      endif
      ! -------                                                --------    .
      l=0
      do i=1,nbseismes                                                     ! Sn ?
        do j=1,nbtps(i)
          if ((D(i)%datatps(j)%sta%staname==R(k)%staname).and.(D(i)%datatps(j)%typeonde=='N') &
            .and.(D(i)%datatps(j)%andS=='S')) l=l+1
        enddo
      enddo
      if (l.gt.0) R(k)%nbSn=0
      if (l.gt.0) R(k)%nbSnT=0
      if ((l.gt.0).and.(.not.allocated(R(k)%resSn))) then
        allocate(R(k)%resSn(maxiter*l,2))
        R(k)%resSn(:,1)=0.0_wr ! résidu
        R(k)%resSn(:,2)=1.e9_wr ! fonction coût
        endif
      if (l.gt.nbseismes) then
        write(*,*)'problème dans initR Sn : trop de ok ',l,nbseismes
        stop
      endif
      ! -------                                                --------    .
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine initR

    ! -----------------------------------------------------------------    .

  subroutine mksynth(nbtps,D,acentroid)
    ! -------                                                  --------    .mh
    ! calcul données synthétiques bruitées ou non
    ! -----------------------------------------------------------------    .
    use typetemps
    use time, only : basetime
    use mt19937, only : genrand_real1, normal
    use pb_direct
    use sub_param
    use sub_misfit
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: nbtps(nbseismes)                      ! nombre de données de temps
    type(dataall), intent (in) :: D(nbseismes)                             ! données
    type (amoho_centroid), intent (in) :: acentroid
    ! -------                                                  --------    .
    type(parametre) :: pmod                                                ! paramètres d'inv. modifiés
    type(dataone), dimension(:), allocatable :: datatempsmod               ! données modifiés
    integer(KIND=wi) :: i, j, mil_1, mil_2, ok
    real(KIND=wr) :: val_1, val_2, aleatoire, bestval, xmin, xmax
    logical :: critique
    type(parametres) :: paramall
    ! -----------------------------------------------------------------    .
    open(505, FILE = 'DATA/newtemps.d',status='replace')
    ! -------                                                  --------    .
    nb_seismes : do j=1,nbseismes
      allocate(datatempsmod(nbtps(j)))
      datatempsmod=D(j)%datatps
      ! -------                                                --------    .
      ! un MODÈLE À DÉFINIR :
      pmod%VC=6.0_wr
      pmod%VM=8.0_wr
      pmod%Zmoho=30._wr
      pmod%VpVs=1.710_wr
      pmod%Tzero%date=D(j)%datatps(1)%tpsR%date
      pmod%Tzero%date%year=2012
      pmod%Tzero%date%month=12
      pmod%Tzero%date%day=15
      pmod%Tzero%date%hour=12
      pmod%Tzero%date%min=0
      pmod%Tzero%sec=30.00_wr
      pmod%Lat=48.25_wr
      pmod%Lon=-2.25_wr
      pmod%Zhypo=15.00_wr
      ! -------                                                --------    .
      write(505,*)'FICHIER synthétique : Lat=',pmod%Lat,'Lon=',pmod%Lon,'Zhypo=',pmod%Zhypo
      write(505,*)
      do i=1,nbtps(j)
        call tempsTheoDirectone(pmod,datatempsmod(i),critique,acentroid)   ! pour le jeu de paramètre tiré et chaque donnée
        if((datatempsmod(i)%typeonde.eq.'G').or.((datatempsmod(i)%typeonde.eq.'N').and.(.not.critique))) then ! distance critique
          if(datatempsmod(i)%andS.eq.'S') then
            write(505,1001)datatempsmod(i)%sta%staname,datatempsmod(i)%typeonde,datatempsmod(i)%coefP, &
            datatempsmod(i)%tpsTh%date%year-2000,datatempsmod(i)%tpsTh%date%month,datatempsmod(i)%tpsTh%date%day, &
            datatempsmod(i)%tpsTh%date%hour,datatempsmod(i)%tpsTh%date%min,datatempsmod(i)%tpsTh%secP, &
            datatempsmod(i)%tpsTh%secS,datatempsmod(i)%andS,datatempsmod(i)%coefS, &
            datatempsmod(i)%sigP,datatempsmod(i)%sigS
          else
            write(505,1002)datatempsmod(i)%sta%staname,datatempsmod(i)%typeonde,datatempsmod(i)%coefP, &
            datatempsmod(i)%tpsTh%date%year-2000,datatempsmod(i)%tpsTh%date%month,datatempsmod(i)%tpsTh%date%day, &
            datatempsmod(i)%tpsTh%date%hour,datatempsmod(i)%tpsTh%date%min,datatempsmod(i)%tpsTh%secP,datatempsmod(i)%sigP
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(505,*)
      write(505,*)' ! -------                                  --------    .'
      write(505,*)'synthétique bruité : '
      write(505,*)
      do i=1,nbtps(j)
        ! -------                                              --------    .
        aleatoire = genrand_real1()
        mil_1 = int(aleatoire*999.99999999999_wr)
        aleatoire = genrand_real1()
        mil_2 = int(aleatoire*999.99999999999_wr)
        select case(mil_1)
        case(0:555)
          val_1=0.025_wr                                                   ! 55,5 % d'être bruité avec une gaussienne de 0.025 s d'écart-type
        case(556:835)
          val_1=0.05_wr                                                    ! 28 % d'être bruité avec une gaussienne de 0.05 s d'écart-type
        case(836:905)
          val_1=0.10_wr                                                    ! 7 % d'être bruité avec une gaussienne de 0.1 s d'écart-type
        case(906:910)
          val_1=0.25_wr                                                    ! 0,5 % d'être bruité avec une gaussienne de 0.25 s d'écart-type
        case(911:1000)
          val_1=0.3_wr                                                     ! 9 % d'être bruité avec une gaussienne de 0.3 s d'écart-type
        case default                                                       ! pour les pct -> Golle, 2013 et Haugmard, 2013 (rapport M2, Annexes, p I.)
          write(*,*)'problème dans mksynth : mil_1 ',mil_1
        end select
        ! -------                                              --------    .
        select case(mil_2)
        case(0:555)
          val_2=0.025_wr
        case(556:835)
          val_2=0.05_wr
        case(836:905)
          val_2=0.10_wr
        case(906:910)
          val_2=0.25_wr
        case(911:1000)
          val_2=0.3_wr
        case default
          write(*,*)'problème dans mksynth : mil_2 ',mil_2
        end select
        ! -------                                              --------    .
        datatempsmod(i)%tpsR = datatempsmod(i)%tpsTh
        aleatoire=val_1 ! aleatoire=normal(0.0_wr,val_1)
        datatempsmod(i)%tpsR%secP = datatempsmod(i)%tpsTh%secP + aleatoire
        datatempsmod(i)%sigP=abs(aleatoire)
        call basetime(datatempsmod(i)%tpsR)                                ! respect du decoupage temps dans la base composite 60/12/365
        if((datatempsmod(i)%typeonde.eq.'G').or.((datatempsmod(i)%typeonde.eq.'N').and.(.not.critique))) then ! distance critique
          if(datatempsmod(i)%andS.eq.'S') then
            aleatoire=val_2 ! aleatoire=normal(0.0_wr,val_2)
            datatempsmod(i)%tpsR%secS = datatempsmod(i)%tpsTh%secS + aleatoire
            datatempsmod(i)%sigS=abs(aleatoire)
            call basetime(datatempsmod(i)%tpsR)                            ! respect du decoupage temps dans la base composite 60/12/365
            write(505,1001)datatempsmod(i)%sta%staname,datatempsmod(i)%typeonde,datatempsmod(i)%coefP, &
            datatempsmod(i)%tpsR%date%year-2000,datatempsmod(i)%tpsR%date%month,datatempsmod(i)%tpsR%date%day, &
            datatempsmod(i)%tpsR%date%hour,datatempsmod(i)%tpsR%date%min,datatempsmod(i)%tpsR%secP, &
            datatempsmod(i)%tpsR%secS,datatempsmod(i)%andS,datatempsmod(i)%coefS, &
            datatempsmod(i)%sigP,datatempsmod(i)%sigS
          else
            write(505,1002)datatempsmod(i)%sta%staname,datatempsmod(i)%typeonde,datatempsmod(i)%coefP, &
            datatempsmod(i)%tpsR%date%year-2000,datatempsmod(i)%tpsR%date%month,datatempsmod(i)%tpsR%date%day, &
            datatempsmod(i)%tpsR%date%hour,datatempsmod(i)%tpsR%date%min,datatempsmod(i)%tpsR%secP, datatempsmod(i)%sigP
          endif
        endif
      enddo
      ! -------                                                --------    . calcul du misfit
      ok=0
      open(507, FILE = 'PARAM/cerclesponderation.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans mksynth : le fichier PARAM/cerclesponderation.d n''existe pas '
        stop
      endif
      i=0
      do while (i.ne.j)
        i=i+1
        read(507,*,iostat = ok)xmin,xmax
      enddo
      close(507)
      call cerclespondone(nbtps(j),datatempsmod,pmod,xmin,xmax,acentroid,chut=.true.)
      ! -------                                                --------    .
      do i=1,nbtps(j)
        call tempsTheoDirectone(pmod,datatempsmod(i),critique,acentroid)
      enddo
      call mvP1_2_Pall(paramall,pmod,j)
      call compute_misfitone(nbtps(j),datatempsmod,bestval,xmin,xmax,'H')
      call mvPall_2_P1(pmod,paramall,j)
      ! -------                                                --------    .
      write(505,*)
      write(505,*)'meilleur misfit :', bestval/real(nbtps(j),wr)           ! valeur de la fonction coût pour ces donnés bruités et les paramètres fixés (minimum) pour ce séisme
      ! -------                                                --------    .
      deallocate(datatempsmod)
    enddo nb_seismes
    ! -------                                                  --------    .
    close(505)
    ! -------                                                  --------    .
    1001 format(a4,1x,a1,1x,i1,1x,5i2.2,f6.3,5x,f7.3,1x,a1,i1.1,3x,f6.3,4x,f6.3)
    1002 format(a4,1x,a1,1x,i1,1x,5i2.2,f6.3,5x,13x,f6.3)
    ! -----------------------------------------------------------------    .
  end subroutine mksynth

    ! -----------------------------------------------------------------    .

  subroutine mksynthallsta(acentroid,dp)
    ! -------                                                  --------    .mh
    ! calcul données synthétiques pour toutes les station connues
    ! -----------------------------------------------------------------    .
    use typetemps, only : amoho_centroid,parametres,stations,parametre,dataone,densityplot,mvPall_2_P1
    use pb_direct, only : tempsTheoDirectone
    use time, only : tempszero, basetime
    ! -------                                                  --------    .
    implicit none
    type (amoho_centroid), intent (in) :: acentroid
    type(densityplot), intent (in) :: dp
    ! -------                                                  --------    .
    type(parametres) :: param_best
    type(stations), dimension(:), allocatable :: datasta
    type(parametre) :: pmod                                                ! paramètres d'inv. modifiés
    type(dataone) :: datatempsmod
    integer(KIND=wi) :: nbsta,i, j, ok
    logical :: critique
    character (LEN=5) :: numberchaine
    ! -----------------------------------------------------------------    .
    open(509, FILE = 'DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans mksynthallsta : le fichier DATA/sta.d n''existe pas '
      stop
    endif
    ! -------                                                  --------    .
    nbsta=0
    do while(ok .eq. 0)                                                    ! boucle pour compter le nombre de lignes du fichier
      read(509,*,iostat = ok)
      if (ok .eq. 0) nbsta = nbsta + 1
    end do
    ! -------                                                  --------    .
    rewind(509)
    allocate(datasta(nbsta))
    ! -------                                                  --------    .
    do i=1,nbsta
      read(509,*,iostat = ok) datasta(i) ! nouveau -> (i)
    enddo
    close(509)
    ! -----------------------------------------------------------------    .
    param_best%VC=dp%VC%moy_100
    param_best%VM=dp%VM%moy_100
    param_best%Zmoho=dp%Zmoho%moy_100
    param_best%VpVs=dp%VpVs%moy_100
    do i=1,nbseismes
      param_best%Zhypo(i)=dp%Zhypo(i)%moy_100
      param_best%lon(i)=dp%lon(i)%moy_100
      param_best%lat(i)=dp%lat(i)%moy_100
      param_best%Tzero(i) = dp%temps_ref(i)
      param_best%Tzero(i)%sec = dp%Tzero(i)%moy_100
      call basetime(param_best%Tzero(i))
    enddo
    ! -----------------------------------------------------------------    .
    nb_seismes : do j=1,nbseismes
      write(numberchaine(1:5),'(i5)')j
      ! -------                                                --------    .
      open(508, FILE = 'OUTPUT/files/tempsTheoOUT_'//trim(adjustl(numberchaine))//'.d',status='replace')
      ! -------                                                --------    .
      call mvPall_2_P1(pmod,param_best,j)
      ! -------                                                --------    .
      write(508,*)'FICHIER synthétique : ',j
      do i=1,nbsta
        ! -------                                              --------    .
        datatempsmod%sta=datasta(i)
        call tempszero(datatempsmod%tpsR%date)
        datatempsmod%tpsR%secP=0.0_wr
        datatempsmod%tpsR%secS=0.0_wr
        call tempszero(datatempsmod%tpsTh%date)
        datatempsmod%tpsTh%secP=0.0_wr
        datatempsmod%tpsTh%secS=0.0_wr
        datatempsmod%typeonde='G'                                          ! Pas de N, pour le moment
        datatempsmod%coefP=0
        datatempsmod%coefS=0
        datatempsmod%andS='S'
        datatempsmod%dTP=0.0_wr
        datatempsmod%dTS=0.0_wr
        datatempsmod%ws=1.0_wr
        datatempsmod%wp=1.0_wr
        datatempsmod%tpsparcP=0.0_wr
        datatempsmod%tpsparcS=0.0_wr
        datatempsmod%depi=0.0_wr
        datatempsmod%dhypo=0.0_wr
        datatempsmod%dcritiqueH=0.0_wr
        datatempsmod%baz=0.0_wr
        datatempsmod%sigP=0.2_wr
        datatempsmod%sigS=0.5_wr
        ! -------                                              --------    .
        call tempsTheoDirectone(pmod,datatempsmod,critique,acentroid)      ! pour le jeu de paramètre tiré et chaque donnée
        ! -------                                              --------    .
        if((datatempsmod%typeonde.eq.'G').or.((datatempsmod%typeonde.eq.'N').and.(.not.critique))) then ! distance critique
          if(datatempsmod%andS.eq.'S') then
            write(508,1003)datatempsmod%sta%staname,datatempsmod%typeonde,datatempsmod%coefP, &
              datatempsmod%tpsTh%date%year-2000,datatempsmod%tpsTh%date%month,datatempsmod%tpsTh%date%day, &
              datatempsmod%tpsTh%date%hour,datatempsmod%tpsTh%date%min,datatempsmod%tpsTh%secP, &
              datatempsmod%tpsTh%secS,datatempsmod%andS,datatempsmod%coefS, &
              datatempsmod%sigP,datatempsmod%sigS
          else
            write(508,1004)datatempsmod%sta%staname,datatempsmod%typeonde,datatempsmod%coefP, &
              datatempsmod%tpsTh%date%year-2000,datatempsmod%tpsTh%date%month,datatempsmod%tpsTh%date%day, &
              datatempsmod%tpsTh%date%hour,datatempsmod%tpsTh%date%min,datatempsmod%tpsTh%secP,datatempsmod%sigP
          endif
        endif
      enddo
      ! -------                                                --------    .
      close(508)
    enddo nb_seismes
    ! -----------------------------------------------------------------    .
    1003 format(a4,1x,a1,1x,i1,1x,5i2.2,f6.3,5x,f7.3,1x,a1,i1.1,3x,f6.3,4x,f6.3)
    1004 format(a4,1x,a1,1x,i1,1x,5i2.2,f6.3,5x,13x,f6.3)
    ! -----------------------------------------------------------------    .
  end subroutine mksynthallsta

    ! -----------------------------------------------------------------    .

  subroutine cerclespond(nbtps,D,param_init,xmin,xmax,acentroid)
    ! -------                                                  --------    .mh
    ! calcul xmin et xmax, diamètres cercles de pondération, pour tous les séismes
    ! lecture dans un fichier PARAM/cerclesponderation.d ou determiné par
    ! cacul, si val=-1
    ! -------                                                  --------    .
    use typetemps, only : dataall, parametre, parametres, amoho_centroid
    use typetemps, only : mvPall_2_P1, mvP1_2_Pall
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(in) :: nbtps(nbseismes)
    type(dataall), intent(inout) :: D(nbseismes)
    type(parametres), intent(inout) :: param_init
    real(KIND=wr), intent(out) :: xmin(nbseismes),xmax(nbseismes)
    type (amoho_centroid), intent (in) :: acentroid
    ! -------                                                  --------    .
    integer(KIND=wi) :: i, ok
    type(parametre) :: param
    ! -----------------------------------------------------------------    .
    ok=0
    open(507, FILE = 'PARAM/cerclesponderation.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans cerclespond : le fichier PARAM/cerclesponderation.d n''existe pas '
      stop
    endif
    do i=1,nbseismes
      read(507,*,iostat = ok)xmin(i),xmax(i)
      call mvPall_2_P1(param,param_init,i)
      call cerclespondone(nbtps(i),D(i)%datatps,param,xmin(i),xmax(i),acentroid)
      call mvP1_2_Pall(param_init,param,i)
    enddo
    close(507)
    ! -----------------------------------------------------------------    .
   end subroutine cerclespond

    ! -----------------------------------------------------------------    .

  subroutine cerclespondone(nbtps,datatps,param_init,xmin,xmax,acentroid,chut)
    ! -------                                                  --------    .mh
    ! calcul xmin et xmax, diamètres cercles de pondération pour un séisme
    ! -------                                                  --------    .
    use typetemps, only : dataone, parametre, amoho_centroid
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(in) :: nbtps
    type(dataone), intent(inout) :: datatps(nbtps)
    type(parametre), intent(in) :: param_init
    real(KIND=wr), intent(out) :: xmin,xmax
    type (amoho_centroid), intent (in) :: acentroid
    logical, intent(in), optional :: chut
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,n
    real(KIND=wr) :: valmax
    logical :: critique
    ! -----------------------------------------------------------------    .
    if (xmax.lt.xmin) then
      valmax=xmax
      xmax=xmin
      xmin=valmax
      write(*,*)'problème dans cerclespondone : xmax < xmin !'
    endif
      ! -------                                                --------    .
    if ((xmin.le.5.0_wr).or.(xmax.le.10.0_wr)) then
      xmin = 0.0_wr
      xmax = 0.0_wr
      valmax = -1.0_wr
      ! -------                                                --------    .
      do i=1,nbtps
        call tempsTheoDirectone(param_init,datatps(i),critique,acentroid)  ! permet le calcul des distances épicentrales
        if (datatps(i)%depi.gt.valmax) valmax = datatps(i)%depi            ! sauve la plus lointaine station
      enddo
      ! -------                                                --------    .
      n = 0
      do while (n.lt.(nbtps/2))
        xmin = xmin + 50.0_wr                                              ! on élargit le cercle tant qu'au mois la moitié des stations s'y trouvent
        n = 0                                                              ! nombre des stations dans le petit cercle
        do i=1,nbtps
          if (datatps(i)%depi.lt.xmin) n=n+1
        enddo
      enddo
      ! -------                                                --------    .
      if (valmax.lt.100.0_wr) valmax = 100.0_wr                            ! limite inférieur
      if (valmax.gt.1000.0_wr) valmax = 1000.0_wr                          ! limite supérieur
      do while (xmax.lt.valmax)
        xmax = xmax + 50.0_wr                                              ! toutes les stations sont dans le grand cercle
      enddo
      ! -------                                                --------    .
      xmax = xmax*2.0_wr                                                   ! ce sont des diamètres
      xmin = xmin*2.0_wr                                                   ! ce sont des diamètres
      ! -------                                                --------    .
      if (xmax.le.xmin) then
        xmax = xmax + 0.1_wr*xmax
      endif
      if(.not.present(chut)) then
        write(*,'(a34,f7.2,1x,f7.2)')" cercles de pondération libres  : ",xmin,xmax
      else
        if(.not.chut) write(*,'(a34,f7.2,1x,f7.2)')" cercles de pondération libres  : ",xmin,xmax
      endif
    else
      if(.not.present(chut)) then
        write(*,'(a34,f7.2,1x,f7.2)')" cercles de pondération fixes   : ",xmin,xmax
      else
        if(.not.chut) write(*,'(a34,f7.2,1x,f7.2)')" cercles de pondération fixes   : ",xmin,xmax
      endif
    endif
    ! -----------------------------------------------------------------    .
  end subroutine cerclespondone

    ! -----------------------------------------------------------------    .

  subroutine catalogue(param,theseisme,find)
    ! -------                                                  --------    .mh
    ! permet la lecture d'un catalogue et de trouver la référence ReNaSS et/ou LDG du séisme
    ! -------                                                  --------    .
    use typetemps, only : parametre, seismes
    use distance_epi
    use time, only :  JDATE, GDATE, difftime, tempszero
    ! -------                                                  --------    .
    implicit none
    type(parametre), intent(in) :: param                                   ! paramètres issus de l'inversion McMC
    type(seismes), intent(out) :: theseisme(2)                             ! paramètres du catalogue
    integer(KIND=wi), intent(out) :: find                                  ! ce séisme est présent dans le catalogue, 0-1-2 fois
    ! -------                                                  --------    .
    type(seismes) :: S                                                     ! pour chaque événement du catalogue
    integer(KIND=wi) :: i,ok
    ! -----------------------------------------------------------------    .
    ok = 0
    i = 0
    find = 0
    ! -------                                                  --------    . initialise
    theseisme(1)%lon=0.0_wr
    theseisme(1)%lat=0.0_wr
    theseisme(1)%pfd=0.0_wr
    call tempszero(theseisme(1)%tps_init%date)
    theseisme(1)%tps_init%sec=0.0_wr
    theseisme(1)%mag=0.0_wr
    theseisme(1)%d_t=0.0_wr
    theseisme(1)%d_epi=0.0_wr
    theseisme(1)%d_p=0.0_wr
    theseisme(1)%name="     "
    ! -------                                                  --------    .
    theseisme(2)%lon=0.0_wr
    theseisme(2)%lat=0.0_wr
    theseisme(2)%pfd=0.0_wr
    call tempszero(theseisme(2)%tps_init%date)
    theseisme(2)%tps_init%sec=0.0_wr
    theseisme(2)%mag=0.0_wr
    theseisme(2)%d_t=0.0_wr
    theseisme(2)%d_epi=0.0_wr
    theseisme(2)%d_p=0.0_wr
    theseisme(2)%name="     "
    ! -------                                                  --------    . lecture catalogue
    open(506, FILE = 'DATA/files/catalogue.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      open(506, FILE = 'files/catalogue.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans catalogue : le fichier DATA/files/catalogue.d n''existe pas '
        write(*,*)'problème dans catalogue : le fichier files/catalogue.d n''existe pas '
        stop
      endif
    endif
    do while(ok .eq. 0)
      ! -------                                                --------    .
      read(506,*,iostat = ok)S%tps_init%date%year,S%tps_init%date%month,S%tps_init%date%day, &
        S%tps_init%date%hour,S%tps_init%date%min,S%tps_init%sec,S%lat,S%lon,S%mag,S%pfd,S%name
      if (ok==0) then
        ! -------                                              --------    .
        call JDATE(S%tps_init%date%Jday,S%tps_init%date%year,S%tps_init%date%month,S%tps_init%date%day)
        call GDATE(S%tps_init%date%Jday,S%tps_init%date%year,S%tps_init%date%month,S%tps_init%date%day)
        call difftime(S%d_t,param%Tzero,S%tps_init)
        ! -------                                              --------    .
        if (abs(S%d_t).lt.60.00_wr) then                                   ! il existe un événement proche (à plus ou moins 1 min)
          i=i+1
          call dellipsgc(param%lat,param%lon,S%Lat,S%Lon,S%d_epi)          ! distance entre les deux epicentres (km)
          S%d_p = param%Zhypo - S%pfd                                      ! différence entre les deux hypocentres (km)
          if (i.le.2) theseisme(i)=S
          find = i
        endif
      endif
    enddo
    if (i.gt.2) then                                                       ! deux séismes max !
      write(*,*)'problème dans catalogue : + de 2 séismes au catalogue'
      find = 2
    endif
    close(506)
    ! -----------------------------------------------------------------    .
  end subroutine catalogue

    ! -----------------------------------------------------------------    .

  subroutine nomseismefichier(Nom)
    ! -------                                                  --------    .mh
    ! nom des nbseismes fichier (xxxx.xx.xx.xx.xx.d) de données dans DATA/
    ! sélectionne nbseisme séismes aléatoirement dans ce fichier
    ! -------                                                  --------    .
    use mt19937, only : genrand_real1
    ! -------                                                  --------    .
    implicit none
    integer (KIND=wi) :: i,j,l,ok
    integer (KIND=wi) :: deja(nbseismes)
    integer (KIND=wi) :: nbfile
    character(len=35) :: Nom(nbseismes)
    logical :: conv
    !logical, save :: present=.false.
    ! -----------------------------------------------------------------    . nombre de lignes
    !if (present) then
    !present=.true.
    open(501, FILE = 'DATA/seismes.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      open(501, FILE = 'seismes.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans nomseismefichier : le fichier DATA/seismes.d n''existe pas '
        stop
      endif
    endif
    nbfile=0
    do while(ok .eq. 0)                                                    ! boucle pour compter le nombre de lignes du fichier
      read(501,*,iostat = ok)
      if (ok .eq. 0) nbfile = nbfile + 1
    end do
    ! -------                                                  --------    .
    rewind(501)
    ! -------                                                  --------    . lecture de nbseisme séismes, sans doublons
    if (nbfile.eq.nbseismes) then
      do i=1,nbseismes
        read(501,*)Nom(i)
      enddo
      close(501)
    ! -------                                                  --------    .
    elseif(nbfile.lt.nbseismes) then
      write(*,*)'problème dans nomseismefichier : nbfile < nb_seismes, manque fichier .dat avec les phases ?'
      stop
    ! -------                                                  --------    . trop de fichier
    else
      ! -------                                                --------    . tire au sort nbseisme sans doublons
      do i=1,nbseismes
        l=int(genrand_real1()*real(nbfile,wr)+1.0_wr)                      ! aléatoire de 1 à nbfile
        if (i==1) then
          deja(1)=l
        else
          conv=.true.
          do while (conv)
            conv=.false.
            do j=1,i-1
              if (deja(j)==l) then
                l=int(genrand_real1()*real(nbfile,wr)+1.0_wr)              ! aléatoire de 1 à nbfile
                conv=.true.
              endif
            enddo
          enddo
        endif
        deja(i)=l
      enddo
      ! -------                                                --------    . lecture des nom de fichier
      do i=1,nbseismes
        do j=1,deja(i)
          read(501,*)Nom(i)
        enddo
        rewind(501)
      enddo
      close(501)
      ! -------                                                --------    . réécriture pour après nbfile=nbseismes
      open(502, FILE = 'DATA/seismes.d',status='replace',iostat = ok)
      if (ok .ne. 0) then
        open(502, FILE = 'seismes.d',status='replace')
      endif
      do i=1,nbseismes
        write(502,'(a)')trim(adjustl(Nom(i)))
      enddo
      close(502)
    endif
    !endif
    ! -----------------------------------------------------------------    .
  end subroutine nomseismefichier

END MODULE datalecture



! *********************************************************************    .
! *********************************************************************    .


