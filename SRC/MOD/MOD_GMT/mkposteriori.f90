! étude a posteriori des paramètres
! janvier 2016
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_posteriori

    use modparam
    use typetemps

    implicit none

    private

    public  :: GMT_posteriori_lonlat
    public  :: GMT_posteriori


CONTAINS


! ---------------------------------------------------------------------    .

  subroutine GMT_posteriori_lonlat(chaine,dp)
    ! -------                                                  --------    .mh
    ! étude a posteriori des paramètres Lon/Lat
    ! -----------------------------------------------------------------    .
    use distance_epi
    use tri, only : melangetab
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    type(densityplot), intent (in) :: dp                                   ! modèles retenus par McMC
    integer(KIND=wi), intent (in) :: chaine
    ! -------                                                  --------    .
    integer(KIND=wi), parameter :: abin=16
    integer(KIND=wi), parameter :: nbptmaxbin=500
    ! -------                                                  --------    .
    real(KIND=wr) :: tl
    real(KIND=wr) :: dmax,dlon,dlat,d,daz
    real(KIND=wr), dimension(:,:), allocatable :: vect1,vect2,vect3,vect4
    integer(KIND=wi) :: i,j,k,m,Noldtime, Nnewtime, ratetime
    integer(KIND=wi) :: minbin
    integer(KIND=wi) :: ok
    integer(KIND=wi) :: bin(abin)
    character (LEN=5) :: numberchaine
    ! -----------------------------------------------------------------    .
    write(numberchaine(1:5),'(i5)')chaine
    ! -----------------------------------------------------------------    .
    ! selection des modèles autour de l'épientre et dont la distance au meilleur modèle est grande
    ! répartition des modèles en fonction de leur azimuth
    ! -----------------------------------------------------------------    .
    open(unit=19,file="OUTPUT/GMT/post_lonlat_lonlat-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=20,file="OUTPUT/GMT/post_lonlat_topdf-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=21,file="OUTPUT/GMT/post_lonlat_VCVpVs-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=22,file="OUTPUT/GMT/post_lonlat_VMZ_moho-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    ! -------                                                  --------    . le modèle le plus loin
    dmax=0.0_wr
    do i=1,dp%nbparam
      call dellipsgc(dp%lat(chaine)%vec(i),dp%lon(chaine)%vec(i), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),d,daz)
      if(d.gt.dmax)dmax=d
    enddo
    ! -------                                                  --------    . binnage
    bin(:)=0
    do i=1,dp%nbparam
      call dellipsgc(dp%lat(chaine)%vec(i),dp%lon(chaine)%vec(i), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),d,daz)
      call doselect(dmax,d,ok)
      if(ok==0) call dobinnage(daz,bin)
    enddo
    ! -------                                                  --------    . binnage max
    minbin=99999999
    do i=1,abin
      if(bin(i).lt.minbin) minbin=bin(i)
    enddo
    ! -------                                                  --------    . au moins !
    if (minbin.lt.10) minbin=10
    ! -------                                                  --------    . au plus !
    if (minbin.gt.nbptmaxbin) minbin=nbptmaxbin
    ! -------                                                  --------    .
    allocate(vect1(3,abin*minbin))
    allocate(vect2(3,abin*minbin))
    allocate(vect3(3,abin*minbin))
    allocate(vect4(3,abin*minbin))
    ! -------                                                  --------    .
    vect1=0.0_wr
    vect2=0.0_wr
    vect3=0.0_wr
    vect4=0.0_wr
    ! -------                                                  --------    .
    j=0
    k=0
    bin(:)=0
    do i=1,dp%nbparam
      ! -------                                                --------    .
      call dellipsgc(dp%lat(chaine)%vec(i),dp%lon(chaine)%vec10000(1,1), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),dlat)    ! calcul de dlat
      call dellipsgc(dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec(i), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),dlon)    ! calcul de dlon
      call dellipsgc(dp%lat(chaine)%vec(i),dp%lon(chaine)%vec(i), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),d,daz)   ! calcul de daz
      ! ---------------------------------------------------------------    .
      call doselect(dmax,d,ok)
      ! ---------------------------------------------------------------    .
      if(ok==0) then
        call dobinnage(daz,bin,amax=minbin,test=ok)
      endif
      ! ---------------------------------------------------------------    .
      if (ok==0) then
        ! -------                                              --------    .
        j=j+1
        if ((daz.ge.0.0_wr).and.(daz.le.360.0_wr)) then
          if (daz.le.90.0_wr) then
            vect1(1,j)=dlon
            vect1(2,j)=dlat
            vect1(3,j)=daz
          elseif(daz.le.180.0_wr) then
            vect1(1,j)=dlon
            vect1(2,j)=-dlat
            vect1(3,j)=daz
          elseif(daz.le.270.0_wr) then
            vect1(1,j)=-dlon
            vect1(2,j)=-dlat
            vect1(3,j)=daz
          elseif(daz.le.360.0_wr) then
            vect1(1,j)=-dlon
            vect1(2,j)=dlat
            vect1(3,j)=daz
          endif
        else
          write(*,*)'problème 1 dans GMT_posteriori_lonlat : daz = ',daz
          stop
        endif
        ! -------                                              --------    .
        vect2(1,j)=dp%Zhypo(chaine)%vec(i)
        vect2(2,j)=dp%Tzero(chaine)%vec(i)
        vect2(3,j)=daz
        ! -------                                              --------    .
        vect3(1,j)=dp%VC%vec(i)
        vect3(2,j)=dp%VpVs%vec(i)
        vect3(3,j)=daz
        ! -------                                              --------    .
        vect4(1,j)=dp%VM%vec(i)
        vect4(2,j)=dp%Zmoho%vec(i)
        vect4(3,j)=daz
        ! -------                                              --------    .
      else
        k=k+1
        ! -------                                              --------    .
        if ((daz.ge.0.0_wr).and.(daz.le.360.0_wr)) then
          if (daz.le.90.0_wr) then
            write(19,REC=k)real(dlon,8),real(dlat,8),real(ok,8)*500._8
          elseif(daz.le.180.0_wr) then
            write(19,REC=k)real(dlon,8),-real(dlat,8),real(ok,8)*500._8
          elseif(daz.le.270.0_wr) then
            write(19,REC=k)-real(dlon,8),-real(dlat,8),real(ok,8)*500._8
          elseif(daz.le.360.0_wr) then
            write(19,REC=k)-real(dlon,8),real(dlat,8),real(ok,8)*500._8
          endif
        else
          write(*,*)'problème 2 dans GMT_posteriori_lonlat : daz = ',daz
          stop
        endif
        ! -------                                              --------    .
        write(20,REC=k)real(dp%Zhypo(chaine)%vec(i),8),real(dp%Tzero(chaine)%vec(i),8),real(ok,8)*500._8
        write(21,REC=k)real(dp%VC%vec(i),8),real(dp%VpVs%vec(i),8),real(ok,8)*500._8
        write(22,REC=k)real(dp%VM%vec(i),8),real(dp%Zmoho%vec(i),8),real(ok,8)*500._8
        ! -------                                              --------    .
      endif
      ! ---------------------------------------------------------------    .
    enddo
    ! -------                                                  --------    .
    close(19)
    close(20)
    close(21)
    close(22)
    ! -----------------------------------------------------------------    .
    call melangetab(j,vect1)                                               ! randomize le tableau
    call melangetab(j,vect2)                                               ! randomize le tableau
    call melangetab(j,vect3)                                               ! randomize le tableau
    call melangetab(j,vect4)                                               ! randomize le tableau
    ! -------                                                  --------    .
    open(unit=15,file="OUTPUT/GMT/post_lonlatok_lonlat-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=16,file="OUTPUT/GMT/post_lonlatok_topdf-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=17,file="OUTPUT/GMT/post_lonlatok_VCVpVs-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=18,file="OUTPUT/GMT/post_lonlatok_VMZ_moho-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    ! -------                                                  --------    .
    m=0
    do i=1,j
      if ((abs(vect1(1,i))+abs(vect1(2,i))+abs(vect1(3,i))).gt.0.00000001_wr) then
        m=m+1
        write(15,REC=m)real(vect1(1,i),8),real(vect1(2,i),8),real(vect1(3,i),8)
        write(16,REC=m)real(vect2(1,i),8),real(vect2(2,i),8),real(vect2(3,i),8)
        write(17,REC=m)real(vect3(1,i),8),real(vect3(2,i),8),real(vect3(3,i),8)
        write(18,REC=m)real(vect4(1,i),8),real(vect4(2,i),8),real(vect4(3,i),8)
      endif
    enddo
    ! -------                                                  --------    .
    close(15)
    close(16)
    close(17)
    close(18)
    ! -------                                                  --------    .
    deallocate(vect1,vect2,vect3,vect4)
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture des script GMT_post LonLat"
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT post LonLat'"
    write(600,*)"#########################################"
    write(600,*)"###########     autocorr       ##########"
    write(600,*)"#########################################"
    ! -----------------------------------------------------------------    .
    write(600,'(a)')"makecpt -Cwysiwyg -T0.0/380.0/5.0 -Z -N > OUTPUT/GMT/colorpostlonlat.cpt"
    write(600,'(a)')"echo 'B 150 150 150 ' >> OUTPUT/GMT/colorpostlonlat.cpt"
    write(600,'(a)')"echo 'F 220 220 220 ' >> OUTPUT/GMT/colorpostlonlat.cpt"
    write(600,'(a)')"echo 'N 200 200 200 ' >> OUTPUT/GMT/colorpostlonlat.cpt"
    ! -----------------------------------------------------------------    .
    write(600,*)"geoproj=-JX4i"
    write(600,*)"file=OUTPUT/GMT/postLonLat-"//trim(adjustl(numberchaine))//".ps"
    ! -----------------------------------------------------------------    .
    ! pour lon et lat :
    write(600,*)"minmax OUTPUT/GMT/post_lonlatok_lonlat-"//trim(adjustl(numberchaine))//".bin -bi3 -I1.5 ", &
      "> OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba0:'"//dp%Lon(chaine)%char//"':/a0:'", &
      dp%Lat(chaine)%char//"':SWen -K -Xc -Yc >  $file"
    write(600,*)"psbasemap $geozone $geoproj -Ba0:'longitude (km)':/a0:'latitude (km)':SWen -K >  $file"
    ! -----------------------------------------------------------------    .
    write(600,*)"gmtset BASEMAP_FRAME_RGB gray"
    write(600,*)"gmtset GRID_PEN_PRIMARY=black GRID_PEN_SECONDARY=black TICK_PEN=black"
    ! -----------------------------------------------------------------    .
    write(600,*)"psbasemap $geozone $geoproj -Ba1f.25g100SWen -K -O >>  $file"
    ! -----------------------------------------------------------------    .
    write(600,*)"gmtset BASEMAP_FRAME_RGB black"
    write(600,*)"gmtset GRID_PEN_PRIMARY=black GRID_PEN_SECONDARY=black TICK_PEN=black"
    ! -----------------------------------------------------------------    .
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlat_lonlat-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlatok_lonlat-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.05 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    ! -------                                                  --------    . rose données
    write(600,*)"psrose ./OUTPUT/GMT/baz"//"-"//trim(adjustl(numberchaine))//".d -: -A10 -S.5in", &
    " -Ggreen -R0/1/0/360 -W1 -F -L'@'/'@'/'@'/'donn\351es' -Y2.5i -X.15i -B.25g0.25/30g30 -O -K -D >> $file"
    write(600,*)"psrose ./OUTPUT/GMT/baz1_"//"-"//trim(adjustl(numberchaine))//".d -: -A10 -S.5i", &
    " -Gblue -R0/1/0/360 -W1 -F -O -K -D >> $file "
    ! -----------------------------------------------------------------    .
    ! pour to et pfd :
    write(600,*)"minmax OUTPUT/GMT/post_lonlatok_topdf-"//trim(adjustl(numberchaine))//".bin -bi3 -I5/1", &
      " > OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba5f1:'"//dp%Zhypo(chaine)%char//"':/a1f.25:'", &
      dp%Tzero(chaine)%char//"':SEwn -K -O -X4.35i -Y-2.5i >>  $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlat_topdf-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlatok_topdf-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.05 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    ! -----------------------------------------------------------------    .
    ! pour VC et VpVs :
    write(600,*)"minmax OUTPUT/GMT/post_lonlatok_VCVpVs-"//trim(adjustl(numberchaine))//".bin -bi3 -I.25/.05", &
      " > OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba.2f.05:'"//dp%VC%char//"':/a.03f.01:'", &
      dp%VpVs%char//"':NWes -K -O -Y4.5i -X-4.5i >>  $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlat_VCVpVs-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlatok_VCVpVs-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.05 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    ! -----------------------------------------------------------------    .
    ! pour VM et Z_moho :
    write(600,*)"geoproj=-JX4i/-4i"
    write(600,*)"minmax OUTPUT/GMT/post_lonlatok_VMZ_moho-"//trim(adjustl(numberchaine))//".bin -bi3 -I.25/5", &
      " > OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba.2f.05:'"//dp%VM%char//"':/a5f1:'", &
      dp%Zmoho%char//"':NsEw -K -O -X4.5i  >>  $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlat_VMZ_moho-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_lonlatok_VMZ_moho-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.05 -COUTPUT/GMT/colorpostlonlat.cpt -K -O >> $file"
    ! -----------------------------------------------------------------    .
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"ps2raster OUTPUT/GMT/postLonLat-"//trim(adjustl(numberchaine))//".ps -Tg -A -P"
    write(600,'(2a)')"mv OUTPUT/GMT/postLonLat-"//trim(adjustl(numberchaine))//".png OUTPUT/figures/postLonLat"// &
      "-"//trim(adjustl(numberchaine))//".png"
    ! -----------------------------------------------------------------    .
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
    CONTAINS
    ! -----------------------------------------------------------------    .

    subroutine doselect(dmax,d,ok)
      ! -------                                                --------    .mh
      ! selection des modèles à colorier !
      ! ---------------------------------------------------------------    .
      implicit none
      ! -------                                                --------    .
      real(KIND=wr), intent (in)  :: dmax,d
      integer, intent (out)  :: ok
      ! -------                                                --------    .
      real(KIND=wr) :: minpct,maxpct
      ! ---------------------------------------------------------------    .
      minpct=0.05_wr
      maxpct=0.5_wr
      ! ---------------------------------------------------------------    .
      if (d.lt.(maxpct*dmax)) then
        if (d.gt.(minpct*dmax)) then
          ok=0
        else
          ok=-1
        endif
      else
        ok=1
      endif
      ! ---------------------------------------------------------------    .
    end subroutine doselect

      ! ---------------------------------------------------------------    .

    subroutine dobinnage(az,bin,amax,test)
      ! -------                                                --------    .mh
      ! répartition azimutale homogène des modèles à colorier !
      ! ---------------------------------------------------------------    .
      implicit none
      ! -------                                                --------    .
      real(KIND=wr), intent (in)  :: az
      integer(KIND=wi), intent (inout)  :: bin(abin)
      integer(KIND=wi), intent (in), optional :: amax
      integer(KIND=wi), intent (out), optional :: test
      ! -------                                                --------    .
      integer(KIND=wi) :: i
      ! ---------------------------------------------------------------    .
      i=int(az/(360._wr/real(abin,wr))+1.0_wr,wi)
      ! -------                                                --------    .
      if (present(amax)) then
        if (bin(i).lt.amax) then
          bin(i)=bin(i)+1
          if (present(amax)) test=0
        else
          if (present(amax)) test=1
        endif
      else
        bin(i)=bin(i)+1
      endif
      ! ---------------------------------------------------------------    .
    end subroutine dobinnage

    ! -----------------------------------------------------------------    .
  end subroutine GMT_posteriori_lonlat


  ! -------------------------------------------------------------------    .

  subroutine GMT_posteriori(chaine,dp,undp)
    ! -------                                                  --------    .mh
    ! étude a posteriori d'un paramètre
    ! -----------------------------------------------------------------    .
    use distance_epi
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    type(densityplot), intent (in) :: dp                                   ! modèles retenus par McMC
    type(densityplot_one), intent (in) :: undp
    integer(KIND=wi), intent (in) :: chaine
    ! -------                                                  --------    .
    real(KIND=wr) :: tl
    real(KIND=wr) :: dlon,dlat,d,daz
    real(KIND=wr) ::dpmin,dpmax,a,b,coef
    integer(KIND=wi) :: i,Noldtime, Nnewtime, ratetime
    character (LEN=5) :: numberchaine
    ! -----------------------------------------------------------------    .
    write(numberchaine(1:5),'(i5)')chaine
    ! -----------------------------------------------------------------    .
    dpmin=1.e9_wr
    dpmax=-1.e9_wr
    do i=1,dp%nbparam
      if (dpmin.gt.undp%vec(i)) dpmin=undp%vec(i)
      if (dpmax.lt.undp%vec(i)) dpmax=undp%vec(i)
    enddo
    ! -----------------------------------------------------------------    .
    open(unit=29,file="OUTPUT/GMT/post_"//undp%name//"_lonlat-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=30,file="OUTPUT/GMT/post_"//undp%name//"_topdf-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=31,file="OUTPUT/GMT/post_"//undp%name//"_VCVpVs-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    open(unit=32,file="OUTPUT/GMT/post_"//undp%name//"_VMZ_moho-"//trim(adjustl(numberchaine))//".bin", &
      STATUS="replace",access='direct',RECL=24)
    ! -----------------------------------------------------------------    .
    do i=1,dp%nbparam
      b= (1.0_wr)/(1.0_wr-(dpmax/dpmin))
      a = -b/dpmin
      coef = a * undp%vec(i) + b                                          ! redimessionne entre 0 et 1
      ! -------                                                --------    .
      call dellipsgc(dp%lat(chaine)%vec(i),dp%lon(chaine)%vec10000(1,1), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),dlat)    ! calcul de dlat
      call dellipsgc(dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec(i), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),dlon)    ! calcul de dlon
      call dellipsgc(dp%lat(chaine)%vec(i),dp%lon(chaine)%vec(i), &
        dp%lat(chaine)%vec10000(1,1),dp%lon(chaine)%vec10000(1,1),d,daz)   ! calcul de daz
      ! -------                                                --------    .
      if ((daz.ge.0.0_wr).and.(daz.le.360.0_wr)) then
        if (daz.le.90.0_wr) then
          write(29,REC=i)real(dlon,8),real(dlat,8),real(coef,8)
        elseif(daz.le.180.0_wr) then
          write(29,REC=i)real(dlon,8),-real(dlat,8),real(coef,8)
        elseif(daz.le.270.0_wr) then
          write(29,REC=i)-real(dlon,8),-real(dlat,8),real(coef,8)
        elseif(daz.le.360.0_wr) then
          write(29,REC=i)-real(dlon,8),real(dlat,8),real(coef,8)
        endif
      else
        write(*,*)'problème 2 dans GMT_posteriori_VC : daz = ',daz
        stop
      endif
      ! -------                                                --------    .
      write(30,REC=i)real(dp%Zhypo(chaine)%vec(i),8),real(dp%Tzero(chaine)%vec(i),8),real(coef,8)
      write(31,REC=i)real(dp%VC%vec(i),8),real(dp%VpVs%vec(i),8),real(coef,8)
      write(32,REC=i)real(dp%VM%vec(i),8),real(dp%Zmoho%vec(i),8),real(coef,8)
      ! -------                                                --------    .
    enddo
    ! -----------------------------------------------------------------    .
    close(29)
    close(30)
    close(31)
    close(32)
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture des script GMT_post "//undp%name
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT post "//undp%name//" '"
    write(600,*)"#########################################"
    write(600,*)"###########     autocorr       ##########"
    write(600,*)"#########################################"
    ! -----------------------------------------------------------------    .
    write(600,'(a)')"makecpt -Cwysiwyg -T0.0/1.0/0.005 -Z -N > OUTPUT/GMT/colorpostvc.cpt"
    write(600,'(a)')"echo 'B 150 150 150 ' >> OUTPUT/GMT/colorpostvc.cpt"
    write(600,'(a)')"echo 'F 220 220 220 ' >> OUTPUT/GMT/colorpostvc.cpt"
    write(600,'(a)')"echo 'N 200 200 200 ' >> OUTPUT/GMT/colorpostvc.cpt"
    ! -----------------------------------------------------------------    .
    write(600,*)"geoproj=-JX4i"
    write(600,*)"file=OUTPUT/GMT/post_"//undp%name//"-"//trim(adjustl(numberchaine))//".ps"
    ! -----------------------------------------------------------------    .
    ! pour lon et lat :
    write(600,*)"minmax OUTPUT/GMT/post_"//undp%name//"_lonlat-"//trim(adjustl(numberchaine))//".bin", &
      " -bi3 -I1.5 > OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba1f.25g100:'longitude (km)':/a1f.25g100:'latitude (km)':SWen -K -Xc  >  $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_"//undp%name//"_lonlat-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostvc.cpt -K -O >> $file"
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    ! -----------------------------------------------------------------    .
    ! pour to et pfd :
    write(600,*)"minmax OUTPUT/GMT/post_"//undp%name//"_topdf-"//trim(adjustl(numberchaine))//".bin -bi3 -I5/1 ", &
      "> OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba5f1:'"//dp%Zhypo(chaine)%char//"':/a1f.25:'", &
      dp%Tzero(chaine)%char//"':SEwn -K -O -X4.5i >>  $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_"//undp%name//"_topdf-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostvc.cpt -K -O >> $file"
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    ! -----------------------------------------------------------------    .
    ! pour VC et VpVs :
    write(600,*)"minmax OUTPUT/GMT/post_"//undp%name//"_VCVpVs-"//trim(adjustl(numberchaine))//".bin -bi3 -I.25/.05 ", &
      "> OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba.2f.05:'"//dp%VC%char//"':/a.03f.01:'", &
      dp%VpVs%char//"':NWes -K -O -Y4.5i -X-4.5i >>  $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_"//undp%name//"_VCVpVs-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostvc.cpt -K -O >> $file"
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    ! -----------------------------------------------------------------    .
    ! pour VM et Z_moho :
    write(600,*)"geoproj=-JX4i/-4i"
    write(600,*)"minmax OUTPUT/GMT/post_"//undp%name//"_VMZ_moho-"//trim(adjustl(numberchaine))//".bin -bi3 -I.25/5 ", &
      "> OUTPUT/GMT/geozone.d "
    write(600,*)"read geozone < OUTPUT/GMT/geozone.d "
    write(600,*)"psbasemap $geozone $geoproj -Ba.2f.05:'"//dp%VM%char//"':/a5f1:'", &
      dp%Zmoho%char//"':NsEw -K -O -X4.5i  >>  $file"
    write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/post_"//undp%name//"_VMZ_moho-"//trim(adjustl(numberchaine))//".bin ", &
      "-bi3 -Sa0.075 -COUTPUT/GMT/colorpostvc.cpt -K -O >> $file"
    ! -----------------------------------------------------------------    .
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"ps2raster OUTPUT/GMT/post_"//undp%name//"-"//trim(adjustl(numberchaine))//".ps -Tg -A -P"
    write(600,'(2a)')"mv OUTPUT/GMT/post_"//undp%name//"-"//trim(adjustl(numberchaine))//".png ", &
      "OUTPUT/figures/post_"//undp%name//"-"//trim(adjustl(numberchaine))//".png"
    ! -----------------------------------------------------------------    .
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
  end subroutine GMT_posteriori

  ! -------------------------------------------------------------------    .

END MODULE figure_posteriori



! *********************************************************************    .
! *********************************************************************    .


