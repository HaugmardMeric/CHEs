! extention programme avec un moho incliné 
! FLAG_non_tabulaire=.true. dans SRC/MOD/modparam.f90
! fevrier 2015
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTmoho_inc

    use modparam

    implicit none

    private

    public  :: GMT_moho


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_moho(acentroid,l,nbtps,xmax,datatps,param)
    ! -------                                                  --------    .mh
    ! Calcul les regressions sur les hodochrones et affiche l'hodochrone
    ! -------                                                  --------    .
    use typetemps
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: l
    integer(KIND=wi), intent(in) :: nbtps
    real(KIND=wr), intent(in) :: xmax
    type(dataone), intent(in) :: datatps(nbtps)
    type(parametre), intent(in)  :: param
    type (amoho_centroid), intent(in) :: acentroid
    ! -------                                                  --------    .
    integer (KIND=wi) :: Noldtime, Nnewtime, ratetime
    integer (KIND=wi) :: j,k,ok,nmax
    real (KIND=wr) :: tl
    real(kind=wr) :: dishypo,Tps,Tpp,pfdsousseisme
    real(KIND=wr) :: v1, v2, lon1,lon2,lat1,lat2
    type(parametre) :: param_best
    character (LEN=5) :: numberfile
    ! -----------------------------------------------------------------    .
    param_best=param
    v1 = 2.0_wr * pi * rT / 360.0_wr                                       ! km / degree en lon
    v2 = 2.0_wr * pi * rT * sin((90.0_wr-param_best%lat)/180.0_wr*pi) /360.0_wr ! km / degree en lat
    ! -------                                                  --------    .
    lon1 = param_best%lon - (xmax / v2 * 1.125_wr) / 2.0_wr
    lon2 = param_best%lon + (xmax / v2 * 1.125_wr) / 2.0_wr
    lat1 = param_best%lat - (xmax / v1 * 1.125_wr) / 2.0_wr
    lat2 = param_best%lat + (xmax / v1 * 1.125_wr) / 2.0_wr
    ! -------                                                  --------    .
    ok=0
    write(numberfile(1:5),'(i5)')l
    open(101, FILE = "OUTPUT/GMT/moho-"//trim(adjustl(numberfile))//".txt",status='replace',iostat = ok)
    nmax=50
    do j=1,nmax
       do k=1,nmax
         param_best%Lon = lon1 - 1.5_wr + abs(lon2-lon1+3.0_wr)/real(nmax,wr)*real(j,wr)
         param_best%Lat = lat1 - 1.5_wr + abs(lat2-lat1+3.0_wr)/real(nmax,wr)*real(k,wr)
         ! -------                                             --------    .
         ! datatps(1)%sta%Lon = lon1 - lon1*0.05_wr + abs(lon1-lon2)/real(nmax,wr)*1.1_wr*real(j,wr)
         ! datatps(1)%sta%Lat = lat1 - lat1*0.05_wr + abs(lat1-lat2)/real(nmax,wr)*1.1_wr*real(k,wr)
         ! -------                                             --------    .
         call refracte_mohovar(acentroid,param_best,datatps(1)%sta,dishypo,Tps,Tpp,pfd=pfdsousseisme)
         write(101,*)param_best%Lon,param_best%Lat,abs(pfdsousseisme)
      enddo
    enddo
    close(101)
    ! -------                                                  --------    .
    !open(106, FILE = "OUTPUT/GMT/sta-"//trim(adjustl(numberfile))//".txt",status='replace',iostat = ok)
    !do i=1,nbtps
    !  write(106,*)datatps(i)%sta%Lon,datatps(i)%sta%Lat,datatps(i)%dTP,datatps(i)%dTS
    !enddo
    !close(106)
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture du script GMT moho"
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"#########################################"
    write(600,*)"################## moho #################"
    write(600,*)"#########################################"
    write(600,*)"echo 'execution du script GMT moho'"
    write(600,*)"file=OUTPUT/GMT/topo_moho-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"gmtset BASEMAP_TYPE plain"
    write(600,'(a10,E13.7,a1,E13.7,a1,E13.7,a1,E13.7)')"geozone=-R",lon1,"/",lon2,"/",lat1,"/",lat2
    write(600,'(a11,E13.7,a1,E13.7,a3)')"geoproj=-JC",param_best%lon,"/",param_best%lat,"/7i"
    write(600,*)"labasemap=-Bpa2.f.5/a1.f.25/a0WeSnz+"
    ! -------                                                  --------    .
    write(600,*)"nearneighbor $geozone -GOUTPUT/GMT/geoid.grd OUTPUT/GMT/moho-"//trim(adjustl(numberfile))//".txt -F -I50k -S50K"
    write(600,*)"grdsample OUTPUT/GMT/geoid.grd -N500 -GOUTPUT/GMT/geoid2.grd"
    write(600,*)"grd2cpt OUTPUT/GMT/geoid2.grd -Csealand > OUTPUT/GMT/colorpal7.cpt "
    write(600,'(2a)')"grdimage $geozone $geoproj OUTPUT/GMT/geoid2.grd -COUTPUT/GMT/colorpal7.cpt ", &
    "$labasemap -Qs -Xc -Yc -K -N0 -S -P > $file"
    ! -------                                                  --------    .
    !write(600,'(2a)')"psxy $geozone $geoproj OUTPUT/GMT/sta-"//trim(adjustl(numberfile))//".txt ", &
    !  "-Sc0.1i -COUTPUT/GMT/colorpal7.cpt -Wthinnest -O -K -P >> $file"
    ! -------                                                  --------    .
    write(600,*)"pscoast -Z0 $geozone $geoproj -Df+ -Ia/blue -W2 -O -K -P >> $file"
    ! -------                                                  --------    .
    write(600,*)"######### cercles de pondération ########"
    write(600,*)"echo ",param%Lon,param%Lat,"0",xmax,xmax, &
      " | psxy $geozone $geoproj -SE -W13,white -O -K >> $file"
    write(600,*)"echo ",param%Lon,param%Lat,"0",xmax,xmax, &
      " | psxy $geozone $geoproj -SE -W10  -O -K >> $file"
    write(600,*)"############### Epicentre ###############"
    write(600,*)"echo",param%Lon,param%Lat," \"
    write(600,*)" | psxy $geozone $geoproj -L -K -O -Wthinnest -Ggreen -Sa0.20i >> $file"
    ! -------                                                  --------    .
    write(600,*)"###### Limites du massif Armoricain #####"
    write(600,*)"psxy $geozone $geoproj -A -W4,gray -O SRC/FILES/limitesMA -K -M >> $file"
    write(600,*)"grdcontour $geozone $geoproj OUTPUT/GMT/geoid2.grd -C2.0 -A -W1 -O -K -P >> $file"
    write(600,'(2a)')"psscale -D-1.5i/4i/5i/.35i -COUTPUT/GMT/colorpal7.cpt -I -O ", &
        "-K -B1:""profondeur du moho \050km\051"": -P -Aa >> $file"
    ! -------                                                  --------    .
    write(600,*)"############### centroide ###############"
    write(600,*)"echo",acentroid%LonC,acentroid%LatC," \"
    write(600,*)" | psxy $geozone $geoproj -L -K -O -Wthinnest -Ggreen -S+0.5i >> $file"
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone $geoproj -JZ-2.5i -Ba0 -O >> $file"
    write(600,*)"ps2raster -Tf -A $file"
    write(600,'(2a)')"mv OUTPUT/GMT/topo_moho-"//trim(adjustl(numberfile))//".pdf ", &
      "OUTPUT/figures/topo_moho-"//trim(adjustl(numberfile))//".pdf"
    write(600,*)"#########################################"
    write(600,*)"ELAPSED=$(($SECONDS-$BEFORE))"
    write(600,*)" echo $ELAPSED secondes"
    call system_clock(Nnewtime,ratetime)
    tl=(real(Nnewtime,wr)-real(Noldtime,wr))/real(ratetime,wr)
    write(*,'(a9,i2.2,'':'',i2.2,'':'',f9.2)')' temps : ',int(tl/3600.0_wr,wi), &
    int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi),(tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr)
    ! -----------------------------------------------------------------    .
  end subroutine GMT_moho

END MODULE figure_GMTmoho_inc



! *********************************************************************    .
! *********************************************************************    .


