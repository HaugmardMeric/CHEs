! permet la création des scripts GMT pour une carte de la région avec les stations
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTmap

    use modparam

    implicit none

    private

    public  :: GMT_map


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_map(l,xmin,xmax,dp,nbtps,datatps)
    ! -------                                                  --------    .mh
    ! production du script GMT produisant une carte des stations
    ! -------                                                  --------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: l
    real(KIND=wr), intent(in) :: xmin,xmax
    type(densityplot), intent(in) :: dp
    integer(KIND=wi), intent(in) :: nbtps
    type(dataone), intent(in) :: datatps(nbtps)
    ! -------                                                  --------    .
    type(stations) :: datasta
    real(KIND=wr) :: v1, v2
    real(KIND=wr) :: X, Y, tl
    real(KIND=wr) :: lon1,lon2,lat1,lat2, azim(36), max
    integer(KIND=wi) :: i,j,k,ok,Noldtime, Nnewtime, ratetime
    character (LEN=5) :: numberfile
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture du script GMT_map "
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT map'"
    write(600,*)"#########################################"
    write(600,*)"###########    carte GMT       ##########"
    write(600,*)"#########################################"
    write(numberfile(1:5),'(i5)')l
    write(600,*)"file=OUTPUT/GMT/map-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"gmtset BASEMAP_TYPE plain"
    write(600,*)"labasemap1=-Bpa2g1.f.5/a1g1.f.25WeSn"
    write(600,*)"grdfile=SRC/FILES/bath1.bin"
    ! -------                                                  --------    .
    v1 = 2.0_wr * pi * rT / 360.0_wr                                       ! km / degree en lon
    v2 = 2.0_wr * pi * rT * sin((90.0_wr-dp%lat(l)%vec10000(1,1))/180.0_wr*pi) /360.0_wr ! km / degree en lat
    ! -------                                                  --------    .
    lon1 = dp%lon(l)%vec10000(1,1) - (xmax / v2 * 1.125_wr) / 2.0_wr
    lon2 = dp%lon(l)%vec10000(1,1) + (xmax / v2 * 1.125_wr) / 2.0_wr
    lat1 = dp%lat(l)%vec10000(1,1) - (xmax / v1 * 1.125_wr) / 2.0_wr
    lat2 = dp%lat(l)%vec10000(1,1) + (xmax / v1 * 1.125_wr) / 2.0_wr
    ! -------                                                  --------    .
    write(600,'(a10,E13.7,a1,E13.7,a1,E13.7,a1,E13.7)')"geozone=-R",lon1,"/",lon2,"/",lat1,"/",lat2
    write(600,'(a11,E13.7,a1,E13.7,a3)')"geoproj=-JC",dp%lon(l)%vec10000(1,1),"/",dp%lat(l)%vec10000(1,1),"/7i"
    ! -------                                                  --------    .
    write(600,*)"bluef=""0/0/100"" "
    write(600,*)"grdgradient $grdfile -A0/270 -GOUTPUT/GMT/grd.gradients -Ne0.6 "
    write(600,*)"grdimage $geozone $geoproj $grdfile -CSRC/FILES/mytopo.cpt -Sc/1 $labasemap1 \"
    write(600,*)"OUTPUT/GMT/grd.gradients -Sn -Xc -X5.5i -Yc -K  > $file"
    write(600,*)"pscoast $geozone $geoproj -Df+ -Ia/$bluef -W1 -O -K >>  $file"
    write(600,*)"grdcontour $grdfile $geozone $geoproj -Ba0 -C1000 -L-10000/-10 -W2 -O -K >> $file"
    write(600,*)"grdcontour $grdfile $geozone $geoproj -Ba0 -C1000 -L10/1000 -W2 -O -K >> $file"
    ! -------                                                  --------    .
    write(600,*)"######### cercles de pondération ########"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"0",xmin,xmin, &
      " | psxy $geozone $geoproj -SE -W13,white -O -K >> $file"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"0",xmin,xmin, &
      " | psxy $geozone $geoproj -SE -W10 -O -K >> $file"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"0",xmax,xmax, &
      " | psxy $geozone $geoproj -SE -W13,white -O -K >> $file"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"0",xmax,xmax, &
      " | psxy $geozone $geoproj -SE -W10  -O -K >> $file"
    ! -------                                                  --------    .
    write(600,*)"######### affiche les stations ##########"
    ok = 0                                                                 ! toutes les stations
    open(800, FILE = 'DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans GMT_map : le fichier data/sta.d n''existe pas '
      stop
    endif
    do while(ok .eq. 0)
      read(800,*,iostat = ok)datasta
      if ((datasta%lon.gt.lon1).and.(datasta%lon.lt.lon2).and.(datasta%lat.gt.lat1).and.(datasta%lat.lt.lat2)) then
        write(600,*)" echo ",datasta%lon,datasta%lat," | psxy $geoproj $geozone -St0.05i -Ggrey -Lk -Wthinnest -O -K >> $file"
      endif
    enddo
    close(800)
    ! ------- relie les stations utilisées                     --------    .
    do i=1,nbtps
      k=0                                                                  ! nb de pointé pur cette station
      do j=1,nbtps
        if (datatps(i)%sta%staname.eq.datatps(j)%sta%staname)   then
          k =  k + 1
          if(datatps(i)%tpsR%secP.eq.datatps(j)%tpsR%secP) then
            if (datatps(j)%andS.eq."S") k =  k + 1
          endif
          if (datatps(i)%tpsR%secP.ne.datatps(j)%tpsR%secP) then
            if (datatps(j)%andS.eq."S") k =  k + 1
          endif
        endif
      enddo
      if ((k.lt.0).or.(k.gt.4)) write(*,*)'problème dans GMT_map : nombre de données par station incorrecte !'
      if (k .eq. 1) then
        write(600,*)" echo -e "" ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"\n", &
        datatps(i)%sta%lon,datatps(i)%sta%lat," "" | psxy $geozone $geoproj -L -K -O -W2,green >> $file"
      endif
      if (k .eq. 2) then
        write(600,*)" echo -e "" ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"\n", &
        datatps(i)%sta%lon,datatps(i)%sta%lat," "" | psxy $geozone $geoproj -L -K -O -W4,yellow >> $file"
      endif
      if (k .eq. 3) then
        write(600,*)" echo -e "" ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"\n", &
        datatps(i)%sta%lon,datatps(i)%sta%lat," "" | psxy $geozone $geoproj -L -K -O -W6,orange >> $file"
      endif
      if (k .eq. 4) then
        write(600,*)" echo -e "" ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1),"\n", &
        datatps(i)%sta%lon,datatps(i)%sta%lat," "" | psxy $geozone $geoproj -L -K -O -W8,red >> $file"
      endif
    enddo
    ! -------                                                  --------    .
    do i=1,nbtps                                                           ! affiche les stations utilisées
      write(600,*)" echo ",datatps(i)%sta%lon,datatps(i)%sta%lat, &
      " | psxy $geoproj $geozone -St0.1i -Gblue -Lk -Wthinnest -O -K >> $file"
    enddo
    do i=1,nbtps                                                           ! affiche les noms des stations utilisées
      write(600,*)" echo ",datatps(i)%sta%lon+0.4_wr,datatps(i)%sta%lat,"8 0 5 6 "//datatps(i)%sta%staname, &
      " | pstext $geozone $geoproj -O -K >> $file"
    enddo
    ! -------                                                  --------    .
    write(600,*)"###### Limites du massif Armoricain #####"
    write(600,*)"psxy $geozone $geoproj -A -W4,gray -O SRC/FILES/limitesMA -K -M >> $file"
    ! -------                                                  --------    .
    write(600,*)"############### Epicentre ###############"
    write(600,*)"echo",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1)," \"
    write(600,*)" | psxy $geozone $geoproj -L -K -O -Wthinnest -Ggreen -Sa0.20i >> $file"
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    ! ------- légende :                                            --------    .
    write(600,*)"################ légende ################"
    write(600,*)"geozone=-R0/1/0/1.25"
    write(600,*)"geoproj=-JX2i/1i"
    X=.1_wr
    Y=.25_wr
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -X-2.5i -O -K >>  $file"
    write(600,*)"echo ",X+.5_wr,Y," 15 0 5 6 ""1 donn\351es""  | pstext $geozone $geoproj -O -K >> $file"
    write(600,*)" echo -e "" ",X,Y,"\n",X+.25_wr,Y," "" | psxy $geozone $geoproj -L -K -O -W2,green >> $file"
    Y=.5_wr
    write(600,*)"echo ",X+.5_wr,Y," 15 0 5 6 ""2 donn\351es"" | pstext $geozone $geoproj -O -K >> $file"
    write(600,*)" echo -e "" ",X,Y,"\n",X+.25_wr,Y," "" | psxy $geozone $geoproj -L -K -O -W4,yellow >> $file"
    Y=.75_wr
    write(600,*)"echo ",X+.5_wr,Y," 15 0 5 6 ""3 donn\351es"" | pstext $geozone $geoproj -O -K >> $file"
    write(600,*)" echo -e "" ",X,Y,"\n",X+.25_wr,Y," "" | psxy $geozone $geoproj -L -K -O -W6,orange >> $file"
    Y=1._wr
    write(600,*)"echo ",X+.5_wr,Y," 15 0 5 6 ""4 donn\351es"" | pstext $geozone $geoproj -O -K >> $file"
    write(600,*)" echo -e "" ",X,Y,"\n",X+.25_wr,Y," "" \ | psxy $geozone $geoproj -L -K -O -W8,red >> $file"
    ! -------                                                  --------    .
    write(600,*)"################## rose #################"
    write(600,*)"geozone=-R0/1/0/1.25"
    write(600,*)"geoproj=-JX2i/2i"
    ! -------                                                  --------    . initialisation compteur azimutal (10 degrés)
    do j=1,36
      azim(j) = 0.0_wr
    enddo
    ! -------                                                  --------    . compteur azimutal (10 degrés) pour toutes données
    open(801, FILE = 'OUTPUT/GMT/baz'//"-"//trim(adjustl(numberfile))//'.d',status='replace',iostat = ok)
    do i=1,nbtps
      write(801,*)datatps(i)%baz, datatps(i)%wp
      do j=1,36
        if((datatps(i)%baz.ge.((real(j,wr)*10.0_wr)-10.00_wr)).and.(datatps(i)%baz.lt.(real(j,wr)*10.0_wr))) then
          azim(j) = azim(j) + datatps(i)%wp
        endif
      enddo
      if (datatps(i)%andS .eq. "S") then
        write(801,*)datatps(i)%baz, datatps(i)%ws
        do j=1,36
          if((datatps(i)%baz.ge.((real(j,wr)*10.0_wr)-10.00_wr)).and.(datatps(i)%baz.lt.(real(j,wr)*10.0_wr))) then
            azim(j) = azim(j) + datatps(i)%ws
          endif
        enddo
      endif
    enddo
    close(801)
    ! -------                                                  --------    . mode du compteur
    max=-1.0_wr
    do j=1,36
      if(max.lt.azim(j)) max = azim(j)
    enddo
    ! -------                                                  --------    . pour données directes
    open(802, FILE = 'OUTPUT/GMT/baz1_'//"-"//trim(adjustl(numberfile))//'.d',status='replace',iostat = ok)
    do i=1,nbtps
      if(datatps(i)%typeonde .eq. "G") then
        write(802,*)datatps(i)%baz, datatps(i)%wp/max
        if (datatps(i)%andS .eq. "S") write(802,*)datatps(i)%baz, datatps(i)%ws/max
      endif
    enddo
    close(802)
    ! -------                                                  --------    . pour données réfractées
    open(803, FILE = 'OUTPUT/GMT/baz2_'//"-"//trim(adjustl(numberfile))//'.d',status='replace',iostat = ok)
    do i=1,nbtps
      if(datatps(i)%typeonde .eq. "N") then
        write(803,*)datatps(i)%baz, datatps(i)%wp/max
        if (datatps(i)%andS .eq. "S") write(803,*)datatps(i)%baz, datatps(i)%ws/max
      endif
    enddo
    close(803)
    ! -------                                                  --------    . légende 1
    write(600,*)"echo ""0 1 15 0 4 LM couverture azimutale "" |", & 
      "pstext $geozone $geoproj -Y5i -O -K >> $file"
    write(600,*)"echo ""0 0.85 15 0 4 LM pond\351r\351e : "" |", & 
      "pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    . légende 2
    write(600,*)"echo ""0 0.6 15 0 4 LM - ondes @;blue;directes"" |", & 
      "pstext $geozone $geoproj -O -K >> $file"
    write(600,*)"echo ""0 0.45 15 0 4 LM - ondes @;olivedrab4;r\351fract\351es"" |", & 
      "pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    . rose données directes
    write(600,*)"psrose ./OUTPUT/GMT/baz"//"-"//trim(adjustl(numberfile))//".d -: -A10 -S.75in", &
      " -Ggray -R0/1/0/360 -F -L'@'/'@'/'@'/'Nord' -Y-1.5i -B.25g0.25/30g30 -O -K >> $file"
    write(600,*)"psrose ./OUTPUT/GMT/baz1_"//"-"//trim(adjustl(numberfile))//".d -: -A10 -S.75i", &
      " -Gblue -R0/1/0/360 -W2 -F -L'@'/'@'/'@'/'Nord' -O -K >> $file "
    ! -------                                                  --------    . rose données réfractées
    write(600,*)"psrose ./OUTPUT/GMT/baz"//"-"//trim(adjustl(numberfile))//".d -: -A10 -S.75in", &
      " -Ggray -R0/1/0/360 -F -L'@'/'@'/'@'/'Nord' -Y-2i -B.25g0.25/30g30 -O -K >> $file "
    write(600,*)"psrose ./OUTPUT/GMT/baz2_"//"-"//trim(adjustl(numberfile))//".d -: -A10 -S.75i", &
      " -Golivedrab4 -R0/1/0/360 -W2 -F -L'@'/'@'/'@'/'Nord' -O >> $file "
    ! -------                                                  --------    . fin
    write(600,*)"ps2raster OUTPUT/GMT/map-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/map-"//trim(adjustl(numberfile))//".pdf ", &
      "OUTPUT/figures/map-"//trim(adjustl(numberfile))//".pdf"
    ! -------                                                  --------    .
    write(600,*)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)"ELAPSED=$(($SECONDS-$BEFORE))"
    write(600,*)" echo $ELAPSED secondes"
    call system_clock(Nnewtime,ratetime)
    tl=(real(Nnewtime,wr)-real(Noldtime,wr))/real(ratetime,wr)
    write(*,'(a9,i2.2,'':'',i2.2,'':'',f9.2)')' temps : ',int(tl/3600.0_wr,wi), &
      int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi),(tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr)
    ! -----------------------------------------------------------------    .
  end subroutine GMT_map

END MODULE figure_GMTmap



! *********************************************************************    .
! *********************************************************************    .


