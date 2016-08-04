! Subroutines permettant la réduction du prior pour lon/lat
! juillet 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE rechercheepi

    use modparam

    implicit none

    private

    public  :: zoneRecherche, initparam

CONTAINS

! ---------------------------------------------------------------------    .


  subroutine initparam(nbtps,D,param_init,pEpis,nb)
    ! -------                                                  --------    .mh
    ! initialise les paramertres hypocentraux et le modèle de Terre
    ! -------                                                  --------    .
    ! VC entre min et max
    ! VM entre min et max, avec VM > VC + 0.1
    ! Zmoho entre min et max
    ! vPvS entre min et max
    ! -------                                                  --------    .
    ! lon,lat -> méthodes des hémisphères, tiré dans un prior restreint par zoneRecherche
    ! Zhypo entre min et max
    ! Tzero correspond à lon,lat,VC et Zhypo
    ! -------                                                  --------    .
    use typetemps
    use time, only : basetime, difftime
    use mt19937
    use distance_epi
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent (in) :: nbtps(nbseismes)                      ! nombre station et nombre de données de temps par séismes
    type(dataall), intent (inout) :: D(nbseismes)                          ! données de temps
    type(parametres), intent (out) :: param_init                           ! paramètres
    type(priorEPI), intent (inout) :: pEpis(nbseismes)
    integer(KIND=wi), intent (in) :: nb                                    ! au carré, nb de cases possible
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j
    integer(KIND=wi), save :: ok=999
    real(KIND=wr), save :: mini_Vc,maxi_Vc,mini_Vm,maxi_Vm,mini_Zmoho
    real(KIND=wr), save :: maxi_Zmoho,mini_VpVs,maxi_VpVs,mini_Zhypo,maxi_Zhypo
    real(KIND=wr) :: ec,moy,val,depi,dhypo
    type(date_sec) :: a_ref_temps,a_time
    real(KIND=wr) :: sommecoef,a_coef
    ! -----------------------------------------------------------------    . lus dans le prior (PARAM/priorIn.d)
    if(ok==999)then                                                        ! lecture unique. -> save
      ok=0
      open(975, FILE = 'PARAM/priorIn_COLD.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans initparam : le fichier PARAM/priorIn_COLD.d n''existe pas '
        stop
      endif
      ! -------                                                --------    .
      read(975,*)mini_Vc,maxi_Vc,ec
      read(975,*)mini_Vm,maxi_Vm,ec
      read(975,*)mini_Zmoho,maxi_Zmoho,ec
      read(975,*)mini_VpVs,maxi_VpVs,ec
      read(975,*)mini_Zhypo,maxi_Zhypo,ec
      close(975)
    endif
    ! -----------------------------------------------------------------    .
    ! ------- parametres du modèle de terre                        ----    .
    ! vitesse dans la croûte des ondes P
    param_init%VC = (maxi_Vc-mini_Vc) * genrand_real3() + mini_Vc
    ! vitesse dans le manteau des ondes P
    val=mini_Vm
    if(mini_Vm.lt.(param_init%VC+0.1_wr)) val = param_init%VC+0.1_wr       ! force VC < VM (réfracction)
    param_init%VM = (maxi_Vm-val) * genrand_real3() + val
    ! profondeur du moho
    param_init%Zmoho = (maxi_Zmoho-mini_Zmoho) * genrand_real3() + mini_Zmoho
    ! ratio de vitesse
    param_init%VpVs = (maxi_VpVs-mini_VpVs) * genrand_real3() + mini_VpVs
    do i=1,nbseismes
      ! ---------------------------------------------------------------    .
      ! ------- parametres hypocentraux                            ----    .
      ! profondeur du séisme
      val=maxi_Zhypo
      if(maxi_Zhypo.gt.(param_init%Zmoho-0.1_wr)) val = param_init%Zmoho-0.1_wr  ! force Zhypo < Zmoho
      param_init%Zhypo(i) = (val-mini_Zhypo) * genrand_real3() + mini_Zhypo
      ! -------                                                --------    .
      ! épicentre
      ! prend une maille parmis les pEpis(i)%nb

      if (pEpis(i)%nb.lt.(nb*nb/6)) then                                 ! si réduction d'un moins un quart -> sinon gros gap azimutal
        ok = int(genrand_real1()*real(pEpis(i)%nb,wr)+1.0_wr)            ! aléatoire de 1 à pEpis(i)%k
        val = genrand_real1()*2.0_wr - 1.0_wr                            ! entre 1.0 et -1.0
        val = val * pEpis(i)%pEpi(1)%distcarre * 360.0_wr / ( 2.0_wr * pi * rT)
        param_init%lat(i) = pEpis(i)%pEpi(ok)%lat + val                  ! ajoute un delta + ou - pEpis(i)%pEpi(1)%distcarre au noeud de la maille choisie
        val = genrand_real1()*2.0_wr - 1.0_wr                            ! entre 1.0 et -1.0
        val = val * pEpis(i)%pEpi(1)%distcarre* 360.0_wr / ( 2.0_wr * pi * rT * cos(param_init%lat(i)*pi/180.0_wr))
        param_init%lon(i) = pEpis(i)%pEpi(ok)%lon + val                  ! ajoute un delta + ou - pEpis(i)%pEpi(1)%distcarre au noeud de la maille choisie
      else
        if (nbtps(i).ge.2) then ! si au moins deux données (ce qui dervrait être le cas, sinon on va pas très loin !)
          ! barycentre : longitudes et latitudes des stations pondérées par les temps d'arrivées des ondes P directes
          ! geometrie simple et bien mieux que l'habitude prise de prendre lon et lat de la premiere station
          do j=1,nbtps(i)

            a_ref_temps%date = D(i)%datatps(nbtps(i))%tpsR%date          ! temps de référence, le plus vieux
            a_ref_temps%sec = D(i)%datatps(nbtps(i))%tpsR%secP+30.0_wr
            call basetime(a_ref_temps)                                   ! reste en base 60/12/365 ...
            sommecoef = 0.0_wr
            param_init%Lon(i) = 0.0_wr
            param_init%Lat(i) = 0.0_wr
            if (D(i)%datatps(j)%typeonde=="G") then ! pour simplifier, on ne travaille qu'avec les ondes directes
              a_time%date=D(i)%datatps(j)%tpsR%date
              a_time%sec=D(i)%datatps(j)%tpsR%secP
              call difftime(a_coef,a_ref_temps,a_time)
              if (a_coef.gt.0.0_wr) then ! au cas ou le plus vieux, n'est pas le plus vieux
                sommecoef = sommecoef + a_coef*a_coef
                param_init%Lon(i) = param_init%Lon(i) + D(i)%datatps(j)%sta%lon * a_coef*a_coef
                param_init%Lat(i) = param_init%Lat(i) + D(i)%datatps(j)%sta%lat * a_coef*a_coef
              endif
            endif
          enddo
          moy=0.01_wr ; ec=0.1_wr
          param_init%Lon(i) = param_init%Lon(i) / sommecoef + normal(moy,ec)
          param_init%Lat(i) = param_init%Lat(i) / sommecoef + normal(moy,ec)
        else
          moy=0.01_wr ; ec=0.1_wr
          param_init%Lon(i) = D(i)%datatps(1)%sta%lon + normal(moy,ec)
          param_init%lat(i) = D(i)%datatps(1)%sta%lat + normal(moy,ec)
        endif
      endif
      ! -------                                                --------    .
      ! temps initial
      ok=-1
      do j=1,nbtps(i)
        if ((ok==-1).and.(D(i)%datatps(j)%typeonde=='G').and.(D(i)%datatps(j)%coefP.lt.3)) then ! première onde P directe
          ok=j
        endif
      enddo
      ! -------                                                --------    .
      if (ok.ne.-1) then                                                   ! il existe des ondes directes
        ! distance épicentrale avec la premiere station
        call dellipsgc(D(i)%datatps(ok)%sta%lat,D(i)%datatps(ok)%sta%lon,param_init%lat(i),param_init%lon(i),depi)
        ! distance hyponcentrale avec la premiere station
        dhypo=sqrt(depi*depi+param_init%Zhypo(i)*param_init%Zhypo(i))
        val = dhypo / param_init%VC
      else                                                                 ! aléatoire ...
        ok = 1
        moy=1.0_wr
        ec=10.0_wr
        val = abs(normal(moy,ec))
      endif

      param_init%Tzero(i)%date = D(i)%datatps(ok)%tpsR%date                ! secondes avant la première arrivée PG
      param_init%Tzero(i)%sec = D(i)%datatps(ok)%tpsR%secP - val
      call basetime(param_init%Tzero(i))                                   ! reste en base 60/12/365 ...
      ! -------                                                --------    .
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine initparam

    ! -----------------------------------------------------------------    .

  subroutine zoneRecherche(nbtps,D,pEpis,nb)
    ! -------                                                  --------    .mh
    ! recherche, sans a priori, la zone de l'épicentre par les données de temps d'arrivée
    ! chaque couple de station défini deux hémisphères (délimité par leur médiatrice)
    ! l'hémishère contenant la station ou l'arrivée de l'onde est la plus tardive est
    ! -------                                                  --------    .
    use typetemps, only : dataall, date_secPS, priorEPI
    use time, only : difftime
    use distance_epi
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: nbtps(nbseismes)
    type(dataall), intent (inout) :: D(nbseismes)
    type(priorEPI), intent (inout) :: pEpis(nbseismes)
    integer(KIND=wi), intent (out) :: nb
    ! -------                                                  --------    .
    type(date_secPS):: atime
    ! -------                                                  --------    .
    integer(KIND=wi) :: nbmax
    integer(KIND=wi), dimension(:,:), allocatable :: grille
    integer(KIND=wi) :: i, j, k, l, m
    integer(KIND=wi) :: Amin, Amax
    ! -------                                                  --------    .
    real(KIND=wr) :: deltag, gsize, dfP, dfS, d1, d2
    real(KIND=wr) :: alon, alat, alon1, alat1, alon0, alat0, alon2, alat2
    character (LEN=5) :: numberchaine
    logical :: done
    ! -----------------------------------------------------------------    .
    deltag=2.5_wr                                                          ! taille du maillage de la grille (km)
    gsize=300.0_wr                                                         ! taille de la grille (km) ; distance de recherche depuis la plus proche station
    nb=int(gsize/deltag)
    ! -------                                                  --------    . définition de la grille de recherche
    ! pour grille(i,j) -> lat = alat + j*deltag* 360.0_wr / ( 2.0_wr * pi * rT)
    ! pour grille(i,j) -> lon = alon + i*deltag* 360.0_wr / sin((90.0_wr - lat)/180.0_wr*pi)
    ! -------                                                  --------    .
    ! si i,j = o -> ~ centré sur la plus proche station
    ! -------                                                  --------    .
    allocate(grille(-nb:nb,-nb:nb))
    ! -----------------------------------------------------------------    .
    open(511, FILE = 'OUTPUT/GMT/script0.sh',status='replace')
    write(511,'(a)')"gmtset LABEL_FONT_SIZE 15"                            ! nouvelles options GMT
    write(511,'(a)')'gmtset HEADER_FONT_SIZE 15'
    write(511,'(a)')"gmtset ANNOT_FONT_PRIMARY Times-Roman"
    write(511,'(a)')"gmtset ANNOT_FONT_SECONDARY Times-Roman"
    write(511,'(a)')"gmtset PAPER_MEDIA A3"
    write(511,'(a)')"gmtset TIME_LANGUAGE FR"
    write(511,'(a)')"gmtset TRANSPARENCY 50"
    write(511,'(a)')"gmtset PLOT_DEGREE_FORMAT dddmm"
    write(511,'(a)')"gmtset BASEMAP_TYPE fancy"
    write(511,'(a)')"gmtset CHAR_ENCODING ISOLatin1+"
    ! -----------------------------------------------------------------    .
    do i=1,nbseismes                                                       ! pour chaque séisme
      nbmax=min(nbtps(i),50)                                               ! max. 50 stations ...
      done=.false.
      do while (.not.done)
        write(numberchaine(1:5),'(i5)')i
        open(508, FILE = 'OUTPUT/GMT/zoneRecherche-'//trim(adjustl(numberchaine))//'.d',status='replace')
        open(509, FILE = 'OUTPUT/GMT/zRsta-'//trim(adjustl(numberchaine))//'.d',status='replace')
        open(510, FILE = 'OUTPUT/GMT/zRepi-'//trim(adjustl(numberchaine))//'.d',status='replace')
        ! -------                                              --------    . initialisation
        grille=0
        ! -------                                              --------    . la plus proche station
        atime=D(i)%datatps(1)%tpsR
        do j=1,nbmax
          call difftime(dfP,dfS,atime,D(i)%datatps(j)%tpsR)
          if (dfP.ge.0.0_wr) then
            alon=D(i)%datatps(j)%sta%lon+0.1_wr                            ! plus ~ 100 m au Nord et à l'Est
            alat=D(i)%datatps(j)%sta%lat+0.1_wr
          endif
        enddo
        ! -------------------------------------------------------------    . on compare que les P, Pn et Pg séparément
        do j=1,nbmax
          do k=j+1,nbmax
            ! -------                                          --------    . pour chaque couple (avec coef > 3)
            if ((D(i)%datatps(k)%coefP.lt.3).and.(D(i)%datatps(j)%coefP.lt.3).and. &
              (D(i)%datatps(k)%typeonde==D(i)%datatps(j)%typeonde)) then
              ! -------                                        --------    . la station la plus proche : 1
              call difftime(dfP,dfS,D(i)%datatps(k)%tpsR,D(i)%datatps(j)%tpsR)
              if (abs(dfP).ge.Tminsec) then                                ! écart minimal de temps entre deux stations
                if (dfP.ge.0.0_wr) then
                  alon1=D(i)%datatps(j)%sta%lon
                  alat1=D(i)%datatps(j)%sta%lat
                  alon2=D(i)%datatps(k)%sta%lon
                  alat2=D(i)%datatps(k)%sta%lat
                else
                  alon1=D(i)%datatps(k)%sta%lon
                  alat1=D(i)%datatps(k)%sta%lat
                  alon2=D(i)%datatps(j)%sta%lon
                  alat2=D(i)%datatps(j)%sta%lat
                endif
                ! -------                                      --------    . données S
                if ((D(i)%datatps(k)%andS=='S').and.(D(i)%datatps(j)%andS=='S').and. &
                  (D(i)%datatps(k)%coefS.lt.3).and.(D(i)%datatps(k)%coefS.lt.3)) then
                  ! -------                                    --------    . données non cohérentes -> coef = 4
                  if ((dfP*dfS).lt.0.0_wr) then
                    write(*,*)'séisme ',i,' pondération +1 (onde S',D(i)%datatps(j)%typeonde,') pour ', &
                      D(i)%datatps(j)%sta%staname,' et ',D(i)%datatps(k)%sta%staname
                      D(i)%datatps(j)%coefS=min(D(i)%datatps(j)%coefS+1,4)
                      D(i)%datatps(k)%coefS=min(D(i)%datatps(k)%coefS+1,4)
                  endif
                endif
                ! -------                                      --------    . pour chaque point de la grille
                do l=-nb,nb
                  do m=-nb,nb
                    ! -------                                  --------    .
                    alat0=alat + real(m,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT)
                    alon0=alon + real(l,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT * cos(alat0*pi/180.0_wr))
                    ! -------                                  --------    . distance entre point de la grille et une station
                    if ((l==0).and.(m==0)) write(509,*)alon2,alat2
                    if ((l==0).and.(m==0)) write(509,*)alon1,alat1
                    if ((alon0.ne.alon1).and.(alat1.ne.alat0)) then
                      call  dellipsgc(alat0,alon0,alat1,alon1,d1)
                    else
                      d1=0.0_wr
                    endif
                    if ((alon0.ne.alon2).and.(alat2.ne.alat0)) then
                      call  dellipsgc(alat0,alon0,alat2,alon2,d2)
                    else
                      d2=0.0_wr
                    endif
                    ! -------                                  --------    . pour l'hemisphere le plus loin grille -= 1
                    if (d2.lt.d1) then
                      grille(l,m)=grille(l,m)-1
                    endif
                    ! -------                                  --------    .
                  enddo
                enddo

              endif
              ! -------                                        --------    .
            endif
            ! -------                                          --------    .
          enddo
        enddo
        ! -------------------------------------------------------------    .
        do l=-nb,nb
          do m=-nb,nb
            alat0=alat + real(m,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT)
            alon0=alon + real(l,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT * cos(alat0*pi/180.0_wr))
            write(508,*)alon0,alat0,grille(l,m)
          enddo
        enddo
        ! -------------------------------------------------------------    .
        ! -------                                              --------    . reste 0 dans la grille de recherche ?
        Amin=10000
        Amax=-10000
        do l=-nb,nb
          do m=-nb,nb
            if(Amin.gt.grille(l,m)) Amin=grille(l,m)
            if(Amax.lt.grille(l,m)) Amax=grille(l,m)
          enddo
        enddo
        ! -------------------------------------------------------------    .
        ! -------                                              --------    . recherche des données P abérentes
        if (Amax.ne.0) then
          do l=-nb,nb
            do m=-nb,nb
              if(grille(l,m)==Amax) then
                alat0=alat + real(m,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT)
                alon0=alon + real(l,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT * cos(alat0*pi/180.0_wr))
                ! -------                                      --------    . quelles stations ?
                do j=1,nbmax
                  do k=j+1,nbmax
                    ! -------                                  --------    . pour chaque couple (avec coef > 3)
                    if ((D(i)%datatps(k)%coefP.lt.3).and.(D(i)%datatps(j)%coefP.lt.3).and. &
                       (D(i)%datatps(k)%typeonde==D(i)%datatps(j)%typeonde)) then
                      ! -------                                --------    . la station la plus proche : 1
                      call difftime(dfP,dfS,D(i)%datatps(k)%tpsR,D(i)%datatps(j)%tpsR)
                      if (abs(dfP).ge.Tminsec) then                        ! écart minimal de temps entre deux stations
                        if (dfP.ge.0.0_wr) then
                          alon1=D(i)%datatps(j)%sta%lon
                          alat1=D(i)%datatps(j)%sta%lat
                          alon2=D(i)%datatps(k)%sta%lon
                          alat2=D(i)%datatps(k)%sta%lat
                        else
                          alon1=D(i)%datatps(k)%sta%lon
                          alat1=D(i)%datatps(k)%sta%lat
                          alon2=D(i)%datatps(j)%sta%lon
                          alat2=D(i)%datatps(j)%sta%lat
                        endif
                        if ((alon0.ne.alon1).and.(alat1.ne.alat0)) then
                          call  dellipsgc(alat0,alon0,alat1,alon1,d1)
                        else
                          d1=0.0_wr
                        endif
                        if ((alon0.ne.alon2).and.(alat2.ne.alat0)) then
                          call  dellipsgc(alat0,alon0,alat2,alon2,d2)
                        else
                          d2=0.0_wr
                        endif
                        ! -------                              --------    . donnée abérente pour l'une des deux stations :
                        if (d2.lt.d1) then
                          write(*,*)'séisme ',i,' pondération +1 (onde P',D(i)%datatps(j)%typeonde,') pour ', &
                            D(i)%datatps(j)%sta%staname,' et ',D(i)%datatps(k)%sta%staname
                          D(i)%datatps(j)%coefP=min(D(i)%datatps(j)%coefP+1,4)
                          D(i)%datatps(k)%coefP=min(D(i)%datatps(k)%coefP+1,4)
                        endif
                      endif
                    endif
                  enddo
                enddo
              endif
            enddo
          enddo
          close(508)
          close(509)
          close(510)
          done=.false.
        else
          done=.true.
        endif
      enddo
      ! ---------------------------------------------------------------    .
      ! -------                                                --------    . taille de la grille de recherche
      k=0
      do l=-nb,nb
        do m=-nb,nb
          if(grille(l,m)==0) then
            alat0=alat + real(m,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT)
            alon0=alon + real(l,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT * cos(alat0*pi/180.0_wr))
            k=k+1
            write(510,*)alon0,alat0
          endif
        enddo
      enddo
      ! -------                                                --------    . tirage aléatoire dans la grille de recherche
      allocate(pEpis(i)%pEpi(k))
      k=0
      do l=-nb,nb
        do m=-nb,nb
          if(grille(l,m)==0) then
            alat0=alat + real(m,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT)
            alon0=alon + real(l,wr)*deltag* 360.0_wr / ( 2.0_wr * pi * rT * cos(alat0*pi/180.0_wr))
            k=k+1
            pEpis(i)%nb=k
            pEpis(i)%pEpi(k)%lon=alon0
            pEpis(i)%pEpi(k)%lat=alat0
            pEpis(i)%pEpi(k)%distcarre=deltag/2.0_wr
          endif
        enddo
      enddo
      ! ---------------------------------------------------------------    .
      ! -------                                                --------    . plot
      d1= deltag / ( 2.0_wr * pi * rT) * 360.0_wr
      d2= deltag / ( 2.0_wr * pi * rT * cos(alat0*pi/180.0_wr)) * 360.0_wr

      write(511,*)"# echo 'execution du script GMT rechercheepi - "//trim(adjustl(numberchaine))//"'"
      write(511,*)"BEFORE=$SECONDS"
      write(511,*)'file=OUTPUT/GMT/init-'//trim(adjustl(numberchaine))//'.ps'
      write(511,*)'geoproj=-JM5i'
      write(511,*)'minmax -I0.001 OUTPUT/GMT/zoneRecherche-'//trim(adjustl(numberchaine))//'.d > toto.txt'
      write(511,*)'read geozone < toto.txt'
      write(511,*)'rm -rf toto.txt'
      write(511,'(a,E13.7,a,E13.7,2a)')'makecpt -Chot -I -T',real(Amin-1,wr),'/',real(Amax,wr),'/1.0', &
      '> OUTPUT/GMT/colorpalinit.cpt'
      write(numberchaine(1:5),'(f5.1)')deltag
      write(511,'(a,i9,a)')'psbasemap $geozone $geoproj -Ba1:."Prior \072 ',int(real(k,wr)*deltag*deltag), &
        ' km \262 (maille '//trim(adjustl(numberchaine))//' km)":SnWe -K -Xc -Yc >  $file'
      write(numberchaine(1:5),'(i5)')i
      write(511,'(a,E13.7,a1,E13.7,2a)')'nearneighbor $geozone -I',d1/5.0_wr,'/',d2/5.0_wr,' OUTPUT/GMT/zoneRecherche', &
        '-'//trim(adjustl(numberchaine))//'.d -F -N4 -S6K -GOUTPUT/GMT/topo.grd'
      write(511,'(2a)')'grdimage $geozone $geoproj OUTPUT/GMT/topo.grd -Qnan ', &
        '-COUTPUT/GMT/colorpalinit.cpt -B0 -O -Sn -K -N >> $file'
      !write(511,*)'sort OUTPUT/GMT/zRsta-'//trim(adjustl(numberchaine))//'.d | uniq | sphtriangulate -Qv -T > OUTPUT/GMT/voronoi.d'
      !write(511,*)'psxy $geoproj $geozone -m -K -O OUTPUT/GMT/voronoi.d -W1,-- >>  $file'
      write(511,*)'pscoast $geozone $geoproj -Df+ -W1 -O -K >> $file'
      write(511,*)'sort OUTPUT/GMT/zRsta-'//trim(adjustl(numberchaine))//'.d | uniq | ', &
        'psxy $geozone $geoproj -St0.25 -W6 -Gred -K -O >> $file'
      !write(511,*)'psxy $geozone $geoproj OUTPUT/GMT/zRepi-'//trim(adjustl(numberchaine))//'.d -Sa0.1 -W1 -Ggreen -O -K >> $file'
      write(511,'(a,f5.1,a1,f5.1,2a)')"awk '{ print $1, $2,"" 0.0,",deltag,",",deltag,"""}' OUTPUT/GMT/zRepi", &
        "-"//trim(adjustl(numberchaine))//".d | psxy $geozone $geoproj -SJ -W1,white -O -K >> $file"
      write(511,'(2a)')'psxy $geozone $geoproj OUTPUT/GMT/ellipse-'//trim(adjustl(numberchaine))//'.txt ', &
        '-Sa0.25 -W1,gray -Gblue -O >> $file'
      write(511,*)"ps2raster OUTPUT/GMT/init-"//trim(adjustl(numberchaine))//".ps -Tf -A -P"
      write(511,'(2a)')"mv OUTPUT/GMT/init-"//trim(adjustl(numberchaine))//".pdf ", &
        "OUTPUT/figures/init-"//trim(adjustl(numberchaine))//".pdf"
      write(511,*)"ELAPSED=$(($SECONDS-$BEFORE))"
      write(511,*)"# echo $ELAPSED secondes"
      ! -------                                                --------    .
    enddo
    close(511)
    deallocate(grille)
    ! -----------------------------------------------------------------    .
  end subroutine zoneRecherche

END MODULE rechercheepi



! *********************************************************************    .
! *********************************************************************    .


