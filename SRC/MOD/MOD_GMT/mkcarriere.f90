! permet la création des scripts GMT pour une carte de la région avec les carrières
! octobre 2015
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .
! merci à Pascal Guterman pascal.guterman@cnrs.fr pour les données         !
! de carrieres (Travail de Frechet & Thouvenot 2012)                       !
! adaptation en ForTran du module python/seiscomp de Pascal Guterman       !
! ---------------------------------------------------------------------    .

MODULE figure_GMTcarriere

    use modparam
    use typetemps, only : date_sec
    use time
    use distance_epi

    implicit none

    private

    public  :: GMT_carriere

  ! -------                                                     --------    .
    TYPE catalogueSiHex
        real(KIND=wr) :: sec,lat,lon,pfd,mag,annee,distcarriere,heure360
        integer(KIND=wi) :: number,mm,jj,aaaa,hh,min,jday
        character(LEN=8) :: orga
        character(LEN=2) :: type
    END TYPE catalogueSiHex
  ! -------                                                     --------    .

CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_carriere(l,seislon,seislat,seistps)
    ! -------                                                  --------    .mh
    ! production du script GMT produisant une carte 
    ! et des diagrammes polaires en vue de discriminer les tirs de carrière
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: l
    real(KIND=wr), intent (in) :: seislon,seislat                          ! lon lat du seisme
    type(date_sec), intent (inout) :: seistps
    ! -------                                                  --------    .
    type(catalogueSiHex) :: ev,refold,refnew
    ! -------                                                  --------    .
    real(KIND=wr) :: v1, v2
    real(KIND=wr) :: clon,clat,cdist
    real(KIND=wr) :: tl
    real(KIND=wr) :: lon1,lon2,lat1,lat2
    integer(KIND=wi) :: i,ok,ok2,Noldtime, Nnewtime, ratetime
    character (LEN=5) :: numberfile
    logical :: existe1
    ! -----------------------------------------------------------------    .
    real(KIND=wr), parameter :: xmax=25.0_wr                               ! distance à l'épicentre
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture du script GMT_carrieres "
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT carrieres'"
    write(600,*)"#########################################"
    write(600,*)"###########    carte GMT       ##########"
    write(600,*)"#########################################"
    write(numberfile(1:5),'(i5)')l
    ! -----------------------------------------------------------------    . 


    inquire (FILE="DATA/files/catalogue_SiHEX_ke_FR.d",exist=existe1)
    if (existe1) then
      ! -------                                                --------    . fichiers
      ok = 0
      open(95, FILE = 'DATA/files/carriere.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème de fichier dans GMT_carriere 1'
        stop
      endif
      open(96, FILE = 'DATA/files/catalogue.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème de fichier dans GMT_carriere 1'
        stop
      endif
      open(97, FILE = 'DATA/files/catalogue_non_tecto.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème de fichier dans GMT_carriere 2'
      stop
      endif
      open(98, FILE = 'DATA/files/catalogue_SiHEX_all_Ma.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème de fichier dans GMT_carriere 1'
        stop
      endif
      open(99, FILE = 'DATA/files/catalogue_SiHEX_ke_FR.d',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème de fichier dans GMT_carriere 2'
        stop
      endif
      ! -------                                                --------    .
      open(101, FILE = "OUTPUT/GMT/cata_all-"//trim(adjustl(numberfile))//".d",status='replace',iostat = ok)
      open(103, FILE = "OUTPUT/GMT/SiHEXP1-"//trim(adjustl(numberfile))//".d",status='replace',iostat = ok)
      open(104, FILE = "OUTPUT/GMT/SiHEXP2-"//trim(adjustl(numberfile))//".d",status='replace',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème de fichier dans GMT_carriere'
        stop
      endif
      ! -------                                                --------    .
      refold%aaaa=1962
      refold%mm=1
      refold%jj=1
      call JDATE (refold%jday,refold%aaaa,refold%mm,refold%jj)

      refnew%aaaa=2015
      refnew%mm=1
      refnew%jj=1
      call JDATE (refnew%jday,refnew%aaaa,refnew%mm,refnew%jj)
      ! ---------------------------------------------------------------    .
      ok2=0
      refold%distcarriere=9999.9_wr
      do while(ok2.eq.0)
        read(95,*,iostat=ok2)clon,clat
        call dellipsgc(clat,clon,seislat,seislon,cdist)
        if (cdist.le.refold%distcarriere) refold%distcarriere=cdist
      end do
      ! ---------------------------------------------------------------    . catalogue_SiHEX_ke_FR.d
      ok=0
      do while(ok.eq.0)
        read(99,12345,iostat=ok)ev%number,ev%jj,ev%mm,ev%aaaa,ev%hh,ev%min,ev%sec,ev%lat,ev%lon,ev%pfd,ev%orga,ev%type,ev%mag
        ! -------                                              --------    .
        call JDATE (ev%jday,ev%aaaa,ev%mm,ev%jj)
        ev%heure360=real(ev%hh,wr)+real(ev%min,wr)/60._wr+ev%sec/3600._wr
        call polar2time(ev%heure360)
        ! -------                                              --------    .
        ok2=0
        ev%distcarriere=9999.9_wr
        clat=seislat
        clon=seislon
        call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
        if(cdist.lt.xmax) then
          rewind(95)
          do while(ok2.eq.0)
            read(95,*,iostat=ok2)clon,clat
            call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
            if (cdist.le.ev%distcarriere) ev%distcarriere=cdist
          end do
          ! -------                                            --------    .
          if ((ok.eq.0).and.(ev%distcarriere.le.xmax)) then
            write(101,*)ev%lon,ev%lat,1.0,'catalogue_SiHEX_ke_FR.d'        ! cata_all
            call JDATE (i,ev%aaaa,1,1)
            lat2=real(ev%aaaa,wr)+real(ev%jday-i,wr)/365.25_wr-real(refold%aaaa,wr)
            write(103,*)ev%heure360,lat2,1.0,'catalogue_SiHEX_ke_FR.d'     ! Plot polaire 1, SiHEXP1
            write(104,*)ev%heure360,ev%distcarriere,1.0,'catalogue_SiHEX_ke_FR.d' ! Plot polaire 2, SiHEXP2
          endif
        endif
        ! -------                                              --------    .
      end do
      ! ---------------------------------------------------------------    . catalogue_SiHEX_all_Ma.d
      ok=0
      do while(ok.eq.0)
        read(98,12346,iostat=ok)ev%number,ev%aaaa,ev%mm,ev%jj,ev%hh,ev%min,ev%sec,ev%lat,ev%lon,ev%pfd,ev%orga,ev%type,ev%mag
        ! -------                                              --------    .
        call JDATE (ev%jday,ev%aaaa,ev%mm,ev%jj)
        ev%heure360=real(ev%hh,wr)+real(ev%min,wr)/60._wr+ev%sec/3600._wr
        call polar2time(ev%heure360)
        ! -------                                              --------    .
        ok2=0
        ev%distcarriere=9999.9_wr
        clat=seislat
        clon=seislon
        call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
        if(cdist.lt.xmax) then
          rewind(95)
          do while(ok2.eq.0)
            read(95,*,iostat=ok2)clon,clat
            call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
            if (cdist.le.ev%distcarriere) ev%distcarriere=cdist
          end do
          ! -------                                            --------    .
          if ((ok.eq.0).and.(ev%distcarriere.le.xmax)) then
            write(101,*)ev%lon,ev%lat,-1.0,'catalogue_SiHEX_all_Ma.d'      ! cata_all
            call JDATE (i,ev%aaaa,1,1)
            lat2=real(ev%aaaa,wr)+real(ev%jday-i,wr)/365.25_wr-real(refold%aaaa,wr)
            write(103,*)ev%heure360,lat2,-1.0,'catalogue_SiHEX_all_Ma.d'   ! Plot polaire 1, SiHEXP1
            write(104,*)ev%heure360,ev%distcarriere,-1.0,'catalogue_SiHEX_all_Ma.d' ! Plot polaire 2, SiHEXP2
          endif
        endif
        ! -------                                              --------    .
      end do
      ! ---------------------------------------------------------------    . catalogue.d
      ok=0
      do while(ok.eq.0)
        read(96,*,iostat=ok)ev%aaaa,ev%mm,ev%jj,ev%hh,ev%min,ev%sec,ev%lat,ev%lon,ev%mag,ev%pfd,ev%orga
        ev%number=0
        ev%type='ke'
        ! -------                                              --------    .
        call JDATE (ev%jday,ev%aaaa,ev%mm,ev%jj)
        ev%heure360=real(ev%hh,wr)+real(ev%min,wr)/60._wr+ev%sec/3600._wr
        call polar2time(ev%heure360)
        ! -------                                              --------    .
        ok2=0
        ev%distcarriere=9999.9_wr
        clat=seislat
        clon=seislon
        call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
        if(cdist.lt.xmax) then
          rewind(95)
          do while(ok2.eq.0)
            read(95,*,iostat=ok2)clon,clat
            call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
            if (cdist.le.ev%distcarriere) ev%distcarriere=cdist
          end do
          ! -------                                            --------    .
          if ((ok.eq.0).and.(ev%distcarriere.le.xmax)) then
            write(101,*)ev%lon,ev%lat,2.0,'catalogue.d'                    ! cata_all
            call JDATE (i,ev%aaaa,1,1)
            lat2=real(ev%aaaa,wr)+real(ev%jday-i,wr)/365.25_wr-real(refold%aaaa,wr)
            write(103,*)ev%heure360,lat2,2.0,'catalogue.d'                 ! Plot polaire 1, SiHEXP1
            write(104,*)ev%heure360,ev%distcarriere,2.0,'catalogue.d'      ! Plot polaire 2, SiHEXP2
          endif
        endif
        ! -------                                              --------    .
      end do
      ! ---------------------------------------------------------------    . catalogue_non_tecto.d
      ok=0
      do while(ok.eq.0)
        read(97,12347,iostat=ok)ev%jj,ev%mm,ev%aaaa,ev%hh,ev%min,i,ev%lat,ev%lon,ev%mag,ev%type
        ev%sec=real(i,wr)
        ev%pfd=0.0_wr
        ev%orga='tir'
        ev%number=0
        ! -------                                              --------    .
        call JDATE (ev%jday,ev%aaaa,ev%mm,ev%jj)
        ev%heure360=real(ev%hh,wr)+real(ev%min,wr)/60._wr+ev%sec/3600._wr
        call polar2time(ev%heure360)
        ! -------                                              --------    .
        ok2=0
        ev%distcarriere=9999.9_wr
        clat=seislat
        clon=seislon
        call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
        if(cdist.lt.xmax) then
          rewind(95)
          do while(ok2.eq.0)
            read(95,*,iostat=ok2)clon,clat
            call dellipsgc(clat,clon,ev%lat,ev%lon,cdist)
            if (cdist.le.ev%distcarriere) ev%distcarriere=cdist
          end do
          ! -------                                            --------    .
          if ((ok.eq.0).and.(ev%distcarriere.le.xmax)) then
            write(101,*)ev%lon,ev%lat,-2.0,'catalogue_non_tecto.d'         ! cata_all
            call JDATE (i,ev%aaaa,1,1)
            lat2=real(ev%aaaa,wr)+real(ev%jday-i,wr)/365.25_wr-real(refold%aaaa,wr)
            write(103,*)ev%heure360,lat2,-2.0,'catalogue_non_tecto.d'      ! Plot polaire 1, SiHEXP1
            write(104,*)ev%heure360,ev%distcarriere,-2.0,'catalogue_non_tecto.d' ! Plot polaire 2, SiHEXP2
          endif
        endif
        ! -------                                              --------    .
      end do
      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .
      close(95)
      close(96)
      close(97)
      close(98)
      close(99)
      close(101)
      close(103)
      close(104)
      ! -------                                                --------    .
      12345 format(2x,i6,4x,i2.2,1x,i2.2,1x,i4.4,1x,i2.2,1x,i2.2,1x,f4.1,2x,f7.2,4x,f7.2,8x,f4.1,4x,a8,a2,4x,f3.1) ! catalogue_SiHEX_ke_FR.d
      12346 format(i6,1x,i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,f4.1,1x,f8.5,f8.5,f6.2,a4,1x,a2,1x,f4.2) ! catalogue_SiHEX_all_Ma.d
      12347 format(i2.2,1x,i2.2,1x,i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,f10.7,1x,f10.7,1x,f10.7,1x,a2) ! catalogue_non_tecto.d
      ! -------                                                --------    .
    endif
    ! -----------------------------------------------------------------    . Plot polaire 1 (heure / année)
    write(600,*)"file=OUTPUT/GMT/carrieresP1-"//trim(adjustl(numberfile))//".ps"
    write(600,'(a,i4.4)')"geozone=-R0/360/0/",refnew%aaaa+1-refold%aaaa
    write(600,'(a)')"geoproj='-JP3i'"
    write(600,'(2a)')"psxy $geozone $geoproj OUTPUT/GMT/SiHEXP1-"//trim(adjustl(numberfile))//".d ", &
      " -Sc0.05i -Wthinnest -COUTPUT/GMT/colortir.cpt -Xc -Yc -Bpa45f7.5g15/a10f1g5wens -K > $file"
    lat1=real(seistps%date%hour,wr)+real(seistps%date%min,wr)/60._wr+seistps%sec/3600._wr
    call polar2time(lat1)
    call JDATE (i,seistps%date%year,1,1)
    call JDATE (seistps%date%Jday,seistps%date%year,seistps%date%month,seistps%date%day)
    lat2=real(seistps%date%year,wr)+real(seistps%date%Jday-i,wr)/365.25_wr-real(refold%aaaa,wr)
    ! -------                                                  --------    . cadran
    cdist=0.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 00h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=3.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 03h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=6.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 06h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=9.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 09h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=12.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 12h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=15.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 15h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=18.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 18h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=21.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*1.16_wr," 12 0 1 CM 21h'", &
      " | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    ! -------                                                  --------    .
    cdist=19.5
    call polar2time(cdist)
    write(600,*)"echo '",cdist,0," 12 67.5 1 CM ",refold%aaaa," ' ", &
    " | pstext $geoproj $geozone -G150/150/150 -O -K -N >> $file"
    write(600,*)"echo '",cdist,real(refnew%aaaa-refold%aaaa,wr)*0.5_wr," 12 67.5 1 CM ", &
    int(refold%aaaa+(refnew%aaaa-refold%aaaa)/2),"' | pstext $geoproj $geozone -G150/150/150  -O -K -N >> $file"
    write(600,*)"echo '",cdist,refnew%aaaa-refold%aaaa," 12 67.5 1 CM ",refnew%aaaa, &
    "' | pstext $geoproj $geozone -G150/150/150  -O -K -N >> $file"
    ! -------                                                  --------    .
    write(600,'(a,2(f10.5,1x),a)')"echo '",lat1,lat2, "'| psxy $geozone $geoproj -Sc0.05i -Wthinnest -Gblack -O >> $file"
    write(600,*)"ps2raster OUTPUT/GMT/carrieresP1-"//trim(adjustl(numberfile))//".ps -Tf -A -P"
    write(600,'(2a)')"mv OUTPUT/GMT/carrieresP1-"//trim(adjustl(numberfile))//".pdf ", &
      "OUTPUT/figures/carrieresP1-"//trim(adjustl(numberfile))//".pdf"

    ! -----------------------------------------------------------------    . Plot polaire 2 (heure / distance)
    write(600,*)"file=OUTPUT/GMT/carrieresP2-"//trim(adjustl(numberfile))//".ps"
    write(600,'(a,E13.7)')"geozone=-R0/360/0/",xmax
    write(600,'(a)')"geoproj='-JP3i'"
    write(600,'(2a)')"psxy $geozone $geoproj OUTPUT/GMT/SiHEXP2-"//trim(adjustl(numberfile))//".d ", &
      " -Sc0.05i -Wthinnest -COUTPUT/GMT/colortir.cpt -Xc -Yc -Bpa45f7.5g15/a5f1g5wens -K > $file"
    lat2=refold%distcarriere
    ! -------                                                  --------    . cadran
    cdist=0.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 00h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=3.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 03h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=6.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 06h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=9.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 09h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=12.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 12h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=15.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 15h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=18.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 18h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    cdist=21.0
    call polar2time(cdist)
    write(600,*)"echo '",cdist,xmax*1.16," 12 0 1 CM 21h' | pstext $geoproj $geozone -O -K -N -Wwhite >> $file"
    ! -------                                                  --------    .
    cdist=19.5
    call polar2time(cdist)
    ! write(600,*)"echo '",cdist,0," 12 67.5 1 CM 0 km'| pstext $geoproj $geozone -O -K -N >> $file"
    write(600,*)"echo '",cdist,xmax*0.5_wr," 12 67.5 1 CM ",int(xmax*0.5_wr), &
      " km' | pstext $geoproj $geozone -G150/150/150 -O -K -N >> $file"
    write(600,*)"echo '",cdist,xmax," 12 67.5 1 CM ",int(xmax), &
      " km'| pstext $geoproj $geozone -G150/150/150 -O -K -N >> $file"
    ! -------                                                  --------    .
    write(600,'(a,2f10.5,a)')"echo '",lat1,lat2, "'| psxy $geozone $geoproj -Sc0.05i -Wthinnest -Gblack -O >> $file"
    write(600,*)"ps2raster OUTPUT/GMT/carrieresP2-"//trim(adjustl(numberfile))//".ps -Tf -A -P"
    write(600,'(2a)')"mv OUTPUT/GMT/carrieresP2-"//trim(adjustl(numberfile))//".pdf ", &
    "OUTPUT/figures/carrieresP2-"//trim(adjustl(numberfile))//".pdf"

    ! -----------------------------------------------------------------    . Plot carte
    write(600,*)"file=OUTPUT/GMT/carrieres-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"gmtset BASEMAP_TYPE plain"
    ! -------                                                  --------    .
    v1 = 2.0_wr * pi * rT / 360.0_wr                                       ! km / degree en lon
    v2 = 2.0_wr * pi * rT * sin((90.0_wr-seislat)/180.0_wr*pi) /360.0_wr   ! km / degree en lat
    ! -------                                                  --------    .
    lon1 = seislon - (xmax / v2 * 1.1_wr)
    lon2 = seislon + (xmax / v2 * 1.1_wr)
    lat1 = seislat - (xmax / v1 * 1.1_wr)
    lat2 = seislat + (xmax / v1 * 1.1_wr)
    ! -------                                                  --------    .
    write(600,'(a10,E13.7,a1,E13.7,a1,E13.7,a1,E13.7)')"geozone=-R",lon1,"/",lon2,"/",lat1,"/",lat2
    write(600,'(a11,E13.7,a1,E13.7,a3)')"geoproj=-JC",seislon,"/",seislat,"/7i"
    ! -------                                                  --------    .
    write(600,*)"bluef=""0/0/100"" "
    write(600,'(2a)')"pscoast $geozone $geoproj -Df+ -S240/255/255 -G180/238/180 -Ia/$bluef -W1 ", &
        "-Xc -Yc -K -Bpa.1f.005/a.1f.005WeSn > $file"
    ! -------                                                  --------    . cercles et croix
    write(600,*)"######### cercles ########"
    do i=1,int(xmax)
      if (mod(i,5).ne.0) then
        write(600,*)"echo '",seislon,seislat,"0",i*2,i*2,"'| psxy $geozone $geoproj -SE -W1,- -O -K >> $file"
      else
        write(600,*)"echo '",seislon,seislat,"0",i*2,i*2,"'| psxy $geozone $geoproj -SE -W2 -O -K >> $file"
        write(600,*)"echo '",seislon+(real(2*i,wr)/v2)/2.0_wr,seislat, &
            "7 0 1 15 ",i,"km ' | pstext $geoproj $geozone -O -K >> $file"
      endif
    enddo
    ! -------                                                  --------    . croix
    write(600,*)"echo -e '",seislon+(real(2*int(xmax,4)+1,wr)/v2)/2.0_wr,seislat,"\n", &
      seislon-(real(2*int(xmax,4)+1,wr)/v2)/2.0_wr,seislat, &
      "' | psxy $geozone $geoproj -W2 -O -K >> $file"
    write(600,*)"echo -e '",seislon,seislat+(real(2*int(xmax,4)+1,wr)/v1)/2.0_wr,"\n", &
      seislon,seislat-(real(2*int(xmax,4)+1,wr)/v1)/2.0_wr, &
      "' | psxy $geozone $geoproj -W2 -O -K >> $file"
    ! -------                                                  --------    . carrières
    write(600,*)"psxy $geozone DATA/files/carriere.d $geoproj -St.05i -Wthinnest,purple -Gblue -O -K >> $file"
    ! -------                                                  --------    . sismicité tectonique
    write(600,'(2a)')"psxy $geozone OUTPUT/GMT/cata_all-"//trim(adjustl(numberfile))//".d ", &
      " $geoproj -Sc0.05i -COUTPUT/GMT/colortir.cpt -Wthinnest -O -K >> $file"
    ! -------                                                  --------    .
    write(600,*)"echo 'N 3' > OUTPUT/GMT/atirlegend.d"
    write(600,*)"echo 'S 0.i c 0.075i 000/000/000 0.25p 0.3i s\351isme' >> OUTPUT/GMT/atirlegend.d"
    write(600,*)"echo 'S 0.i c 0.075i 170/000/000 0.25p 0.3i me+km catalogue' >> OUTPUT/GMT/atirlegend.d"
    write(600,*)"echo 'S 0.i c 0.075i 255/028/000 0.25p 0.3i me+km Si-Hex' >> OUTPUT/GMT/atirlegend.d"
    write(600,*)"echo 'S 0.i c 0.075i 255/255/000 0.25p 0.3i ke Si-Hex' >> OUTPUT/GMT/atirlegend.d"
    write(600,*)"echo 'S 0.i c 0.075i 255/255/170 0.25p 0.3i ke catalogue' >> OUTPUT/GMT/atirlegend.d"
    write(600,*)"echo 'S 0.i t 0.075i 000/000/255 0.25p,purple 0.3i carri\350res' >> OUTPUT/GMT/atirlegend.d"
    ! -------                                                  --------    .
    write(600,*)"pslegend -Dx3.5i/-0.4i/7i/1.i/TC $geozone $geoproj -O -K OUTPUT/GMT/atirlegend.d >> $file"
    ! -------                                                  --------    .
    write(600,*)"echo '",seislon,seislat,"' | psxy $geozone $geoproj -Sc0.05i -Gblack -Wthinnest -O >> $file"
    ! -------                                                  --------    . fin
    write(600,*)"ps2raster OUTPUT/GMT/carrieres-"//trim(adjustl(numberfile))//".ps -Tf -A -P"
    write(600,'(2a)')"mv OUTPUT/GMT/carrieres-"//trim(adjustl(numberfile))//".pdf ", &
      "OUTPUT/figures/carrieres-"//trim(adjustl(numberfile))//".pdf"
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
  end subroutine GMT_carriere

  ! -------------------------------------------------------------------    .

  subroutine polar2time(une_heure)
    ! -------                                                  --------    .mh
    ! transforme coordonnées orloge en polaire (sens inverse)
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr), intent (inout) :: une_heure
    ! -----------------------------------------------------------------    .
    une_heure=((une_heure*(-1.0_wr))+24.0_wr)
    if (une_heure.ge.24.0_wr) une_heure=une_heure-24.0_wr
    une_heure=une_heure/24.0_wr*360.0_wr+90.0_wr
    ! -----------------------------------------------------------------    .
  end subroutine polar2time

  ! -------------------------------------------------------------------    .

END MODULE figure_GMTcarriere



! *********************************************************************    .
! *********************************************************************    .


