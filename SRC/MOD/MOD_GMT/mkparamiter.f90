! permet la création des scripts GMT pour les figures sur chaque parametre vs itération et sur les fonctions d’autocorrélation
! mars-octobre 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTpar

    use modparam

    implicit none

    private

    public  :: GMT_param
    public  :: RVB


CONTAINS

! ---------------------------------------------------------------------    .

subroutine GMT_param(dp,nbChaineMV,nmod)
    ! -------                                                  --------    .mh
    ! figures des chaines après échantillonage   
    ! -----------------------------------------------------------------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    type(densityplot), intent (in) :: dp                                   ! modèles retenus par McMC
    integer(KIND=wi), intent (in) :: nbChaineMV
    integer(KIND=wi), intent (in) :: nmod(nbChaineMV)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k,l,m,maxiter,val,Noldtime, Nnewtime, ratetime
    real(KIND=wr) :: tl
    character(LEN=5) :: char
    character(LEN=11) :: color, char1, char2
    character (LEN=5) :: numberfile
    ! -----------------------------------------------------------------    . ecriture des paramètre de terre, par chaîne, avec rééchantillonnage
    k=0
    maxiter=0
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      open(unit=1200+8*i,file="OUTPUT/GMT/theVC"//char//".bin",STATUS="replace",access='direct',RECL=16)
      open(unit=1201+8*i,file="OUTPUT/GMT/theVM"//char//".bin",STATUS="replace",access='direct',RECL=16)
      open(unit=1202+8*i,file="OUTPUT/GMT/theVpVs"//char//".bin",STATUS="replace",access='direct',RECL=16)
      open(unit=1203+8*i,file="OUTPUT/GMT/theZmoho"//char//".bin",STATUS="replace",access='direct',RECL=16)
      l=0
      do j=1,nmod(i)
        k=k+1
        val = max(100,(nmod(i)/2000))
        if (mod(k,val)==0) then
          l=l+1
          write(1200+8*i,rec=l)real(dp%VC%vec(k),8),real(j,8)
          write(1201+8*i,rec=l)real(dp%VM%vec(k),8),real(j,8)
          write(1202+8*i,rec=l)real(dp%VpVs%vec(k),8),real(j,8)
          write(1203+8*i,rec=l)real(dp%Zmoho%vec(k),8),real(j,8)
        endif
      if(maxiter.lt.j)maxiter=j
    enddo
    close(1200+8*i)
    close(1201+8*i)
    close(1202+8*i)
    close(1203+8*i)
    enddo
    ! -----------------------------------------------------------------    . ecriture des paramètre hypocentraux, par chaîne, avec rééchantillonnage
    do m=1,nbseismes
      write(numberfile(1:5),'(i5)')m
      k=0
      maxiter=0
      do i=1,nbChaineMV
        write(char,'(i5)')10000+i
        open(unit=1204+8*i,file="OUTPUT/GMT/theLon"//char//"-"//trim(adjustl(numberfile))//".bin", &
        STATUS="replace",access='direct',RECL=16)
        open(unit=1205+8*i,file="OUTPUT/GMT/theLat"//char//"-"//trim(adjustl(numberfile))//".bin", &
        STATUS="replace",access='direct',RECL=16)
        open(unit=1206+8*i,file="OUTPUT/GMT/theZhypo"//char//"-"//trim(adjustl(numberfile))//".bin", &
        STATUS="replace",access='direct',RECL=16)
        open(unit=1207+8*i,file="OUTPUT/GMT/theTzero"//char//"-"//trim(adjustl(numberfile))//".bin", &
        STATUS="replace",access='direct',RECL=16)
        l=0
        do j=1,nmod(i)
          k=k+1
          val = max(100,(nmod(i)/2000))
          if (mod(k,val)==0) then
            l=l+1
            write(1204+8*i,rec=l)real(dp%Lon(m)%vec(k),8),real(j,8)
            write(1205+8*i,rec=l)real(dp%Lat(m)%vec(k),8),real(j,8)
            write(1206+8*i,rec=l)real(dp%Zhypo(m)%vec(k),8),real(j,8)
            write(1207+8*i,rec=l)real(dp%Tzero(m)%vec(k),8),real(j,8)
          endif
          if(maxiter.lt.j)maxiter=j
        enddo
      close(1204+8*i)
      close(1205+8*i)
      close(1206+8*i)
      close(1207+8*i)
      enddo
    enddo
    m=1
    ! -----------------------------------------------------------------    .
    ! plot les différentes réalisations pour chaque parametre
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture des script GMT_param "
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT param vs iter'"
    write(600,*)"#########################################"
    write(600,*)"###########  param vs iter     ##########"
    write(600,*)"#########################################"
    ! -----------------------------------------------------------------    .
    write(600,*)"geoproj=-JX13.5i/4.5i"
    write(600,*)"file=OUTPUT/GMT/VCVM_histo.ps"
    write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%VC%themin,"/",dp%VC%themax
    ! -----------------------------------------------------------------    .
    write(char1(1:11),'(i11)')(maxiter/1000)*250
    write(char2(1:11),'(i11)')val
    ! -----------------------------------------------------------------    .    
    write(600,'(3a,E13.7,2a)')"psbasemap $geozone $geoproj ", &
      "-Ba"//trim(adjustl(char1))//":'mod\350les \050\351chantillonnage 1/"//trim(adjustl(char2))//"\0", &
      "51':/a",real(int(((dp%VC%themax-dp%VC%themin)/5.0_wr)*10.0_wr,wi),wr)/10._wr,":'"//dp%VC%char//"':nSeW", &
      " -K -Y6.3i -Xc >  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theVC"//char//".bin -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%VM%themin,"/",dp%VM%themax
    write(600,'(3a,E13.7,2a)')"psbasemap $geozone $geoproj ", &
      " -Ba"//trim(adjustl(char1))//":'mod\350les \050\351chantillonnage 1/"//trim(adjustl(char2))//"\0", &
      "51':/a",real(int(((dp%VM%themax-dp%VM%themin)/5.0_wr)*10.0_wr),wr)/10.0_wr,":'"//dp%VM%char//"':nSeW", &
      " -K -O -Y-5.5i >>  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theVM"//char//".bin -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"ps2raster OUTPUT/GMT/VCVM_histo.ps -Tf -A"
    write(600,'(a)')"mv OUTPUT/GMT/VCVM_histo.pdf OUTPUT/figures/VCVM_histo.pdf"
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(600,*)"file=OUTPUT/GMT/VpVsZmoho_histo.ps"
    write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%VpVs%themin,"/",dp%VpVs%themax
    write(600,'(2a)')"psbasemap $geozone $geoproj -Ba"//trim(adjustl(char1))//":'mod\350les ", &
      "\050\351chantillonnage 1/"//trim(adjustl(char2))//"\051':/a0.05:'"//dp%VpVs%char//"':nSeW -K -Y6.3i -Xc >  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theVpVs"//char//".bin -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%Zmoho%themin,"/",dp%Zmoho%themax
    write(600,'(3a,E13.7,2a)')"psbasemap $geozone $geoproj ", &
      " -Ba"//trim(adjustl(char1))//":'mod\350les \050\351chantillonnage 1/"//trim(adjustl(char2))//"\0", &
      "51':/a",real(int(((dp%Zmoho%themax-dp%Zmoho%themin)/5.0_wr)*5.0_wr),wr)/10._wr,":'"//dp%Zmoho%char//"':nSeW", &
      " -K -O -Y-5.5i >>  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theZmoho"//char//".bin -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"ps2raster OUTPUT/GMT/VpVsZmoho_histo.ps -Tf -A"
    write(600,'(a)')"mv OUTPUT/GMT/VpVsZmoho_histo.pdf OUTPUT/figures/VpVsZmoho_histo.pdf"
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    do m=1,nbseismes
      write(numberfile(1:5),'(i5)')m
      write(600,*)"file=OUTPUT/GMT/LatLon_histo"//"-"//trim(adjustl(numberfile))//".ps"
      write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%Lon(m)%themin,"/",dp%Lon(m)%themax
      write(600,'(3a,E13.7,2a)')"psbasemap $geozone $geoproj ", &
        " -Ba"//trim(adjustl(char1))//":'mod\350les \050\351chantillonnage 1/"//trim(adjustl(char2))//"\0", &
        "51':/a",real(int(((dp%Lon(m)%themax-dp%Lon(m)%themin)/5.0_wr)*200.0_wr),wr)/200.0_wr,":'"//dp%Lon(m)%char//"':nSeW", &
        " -K -Y6.3i -Xc >  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')10000+i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theLon"//char//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
      enddo
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
      write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%Lat(m)%themin,"/",dp%Lat(m)%themax
      write(600,'(3a,E13.7,2a)')"psbasemap $geozone $geoproj ", &
        " -Ba"//trim(adjustl(char1))//":'mod\350les \050\351chantillonnage 1/"//trim(adjustl(char2))//"\0", &
        "51':/a",real(int(((dp%Lat(m)%themax-dp%Lat(m)%themin)/5.0_wr)*200.0_wr),wr)/200.0_wr,":'"//dp%Lat(m)%char//"':nSeW", &
        " -K -O -Y-5.5i >>  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')10000+i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theLat"//char//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
      enddo
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
      write(600,*)"ps2raster OUTPUT/GMT/LatLon_histo"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/LatLon_histo"//"-"//trim(adjustl(numberfile))//".pdf ", &
        " OUTPUT/figures/LatLon_histo"//"-"//trim(adjustl(numberfile))//".pdf"
    enddo
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    do m=1,nbseismes
      write(numberfile(1:5),'(i5)')m
      write(600,*)"file=OUTPUT/GMT/ZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".ps"
      write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%Tzero(m)%themin,"/",dp%Tzero(m)%themax
      write(600,'(3a,E13.7,2a)')"psbasemap $geozone $geoproj ", &
        " -Ba"//trim(adjustl(char1))//":'mod\350les \050\351chantillonnage 1/"//trim(adjustl(char2))//"\0", &
        "51':/a",real(int(((dp%Tzero(m)%themax-dp%Tzero(m)%themin)/5.0_wr)*10.0_wr),wr)/10.0_wr,":'"//dp%Tzero(m)%char//"':nSeW", &
        " -K -Y6.3i -Xc >  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')10000+i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theTzero"//char//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
      enddo
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
      write(600,'(a12,i9.9,a1,E13.7,a1,E13.7)')"geozone=-R0/",maxiter,"/",dp%Zhypo(m)%themin,"/",dp%Zhypo(m)%themax
      write(600,'(3a,E13.7,2a)')"psbasemap $geozone $geoproj ", &
        " -Ba"//trim(adjustl(char1))//":'mod\350les \050\351chantillonnage 1/"//trim(adjustl(char2))//"\0", &
        "51':/a",real(int(((dp%Zhypo(m)%themax-dp%Zhypo(m)%themin)/5.0_wr)*10.0_wr),wr)/10.0_wr,":'"//dp%Zhypo(m)%char//"':nSeW", &
        " -K -O -Y-5.5i >>  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')10000+i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/theZhypo"//char//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
      enddo
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
      write(600,*)"ps2raster OUTPUT/GMT/ZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/ZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".pdf ", &
        "OUTPUT/figures/ZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".pdf"
    enddo
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


    ! -----------------------------------------------------------------    .
    ! plot figures fonctions d’autocorrélation
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture des script GMT_autocorr "
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT autocorr'"
    write(600,*)"#########################################"
    write(600,*)"###########     autocorr       ##########"
    write(600,*)"#########################################"
    ! -----------------------------------------------------------------    .
    write(600,*)"geoproj=-JX13.5i/4.5i"
    ! -----------------------------------------------------------------    .
    write(600,'(a,i9.9,a)')"geozone=-R-100/",autocorr,"/-1/1"
    write(600,*)"file=OUTPUT/GMT/autoVCVM_histo.ps"
    ! -----------------------------------------------------------------    .
    write(600,'(a)')"psbasemap $geozone $geoproj -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur V@-C @-':nSeW -K -Y6.3i -Xc >  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/autovar_VC"//trim(adjustl(char))//".txt ", &
        "-Wthinnest,"//color//" -O -K >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    write(600,'(a)')"psbasemap $geozone $geoproj -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur V@-M @-':nSeW -K -O -Y-5.5i >>  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/autovar_VM"//trim(adjustl(char))//".txt ", &
        "-Wthinnest,"//color//" -O -K >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"ps2raster OUTPUT/GMT/autoVCVM_histo.ps -Tf -A"
    write(600,'(a)')"mv OUTPUT/GMT/autoVCVM_histo.pdf OUTPUT/figures/autoVCVM_histo.pdf"
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(600,*)"file=OUTPUT/GMT/autoVpVsZmoho_histo.ps"
    write(600,'(2a)')"psbasemap $geozone $geoproj -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur V@-P @- / V@-S @- ':nSeW ", &
      " -K -Y6.3i -Xc >  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/autovar_VpVs"//trim(adjustl(char))//".txt ", &
        "-Wthinnest,"//color//" -O -K >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
    write(600,'(2a)')"psbasemap $geozone $geoproj -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur la profondeur du moho ':nSeW ", &
      "-K -O -Y-5.5i >>  $file"
    do i=1,nbChaineMV
      write(char,'(i5)')i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/autovar_Zmoho"//trim(adjustl(char))//".txt ", &
        "-Wthinnest,"//color//" -O -K >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"ps2raster OUTPUT/GMT/autoVpVsZmoho_histo.ps -Tf -A"
    write(600,'(a)')"mv OUTPUT/GMT/autoVpVsZmoho_histo.pdf OUTPUT/figures/autoVpVsZmoho_histo.pdf"
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    do m=1,nbseismes
      write(numberfile(1:5),'(i5)')m
      write(600,*)"file=OUTPUT/GMT/autoLatLon_histo"//"-"//trim(adjustl(numberfile))//".ps"
      write(600,'(2a)')"psbasemap $geozone $geoproj -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur la longitude':nSeW", &
        " -K -Y6.3i -Xc >  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj ", &
          "OUTPUT/GMT/autovar_lon_"//trim(adjustl(numberfile))//"_"//trim(adjustl(char))//".txt ", &
          "-Wthinnest,"//color//" -O -K >>  $file"
      enddo
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
      write(600,'(2a)')"psbasemap $geozone $geoproj -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur la latitude':nSeW ", &
        "-K -O -Y-5.5i >>  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj ", &
          "OUTPUT/GMT/autovar_lat_"//trim(adjustl(numberfile))//"_"//trim(adjustl(char))//".txt ", &
          "-Wthinnest,"//color//" -O -K >>  $file"
      enddo
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
      write(600,*)"ps2raster OUTPUT/GMT/autoLatLon_histo"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/autoLatLon_histo"//"-"//trim(adjustl(numberfile))//".pdf ", &
        "OUTPUT/figures/autoLatLon_histo"//"-"//trim(adjustl(numberfile))//".pdf"
    enddo
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    do m=1,nbseismes
      write(numberfile(1:5),'(i5)')m
      write(600,*)"file=OUTPUT/GMT/autoZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".ps"

      write(600,'(2a)')"psbasemap $geozone $geoproj -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur le temps initial ':nSeW ", &
        "-K -Y6.3i -Xc >  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj ", &
          "OUTPUT/GMT/autovar_Tzero_"//trim(adjustl(numberfile))//"_"//trim(adjustl(char))//".txt ", &
          " -Wthinnest,"//color//" -O -K >>  $file"
      enddo
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
      write(600,'(2a)')"psbasemap $geozone $geoproj ", &
        " -Ba1000f500g10000:k:/a0.5g5:'r@-k@- sur la profondeur de l\234hypocentre':nSeW -K -O -Y-5.5i >>  $file"
      do i=1,nbChaineMV
        write(char,'(i5)')i
        call RVB(i,nbChaineMV,color)
        write(600,*)"psxy $geozone $geoproj ", &
          "OUTPUT/GMT/autovar_Zhypo_"//trim(adjustl(numberfile))//"_"//trim(adjustl(char))//".txt ", &
          "-Wthinnest,"//color//" -O -K >>  $file"
      enddo

      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
      write(600,*)"ps2raster OUTPUT/GMT/autoZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/autoZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".pdf ",&
        "OUTPUT/figures/autoZhypoTzero_histo"//"-"//trim(adjustl(numberfile))//".pdf"
    enddo
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

  end subroutine GMT_param

    ! -----------------------------------------------------------------    .

  subroutine RVB(i,n,color)
    ! -------                                                  --------    .mh
    ! color en RVB sur vecteur de i sur n
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent (in) :: i,n
    character(LEN=11), intent (out) :: color
    ! -------                                                  --------    .
    real(KIND=wr) :: val,r,v,b
    ! -----------------------------------------------------------------    .
    val=15.0_wr+real(int(real(n-i,wr)/real(n,wr)*250.0_wr,wi),wr)*3.0_wr   ! décale pour le fun (+15)
    ! -------                                                  --------    . rouge
    if (val.gt.255.0_wr) then
      r=0.0_wr
    endif
    if (val.lt.(255.0_wr/2.0_wr)) then
      r=real(int(val),wr)*2.0_wr
    endif
    if ((val.gt.(255.0_wr/2.0_wr)).and.(val.lt.255.0_wr)) then
      r=255.0_wr*2.0_wr-real(int(val),wr)*2.0_wr
    endif
    ! -------                                                  --------    . vert
    if (val.gt.(255.0_wr*2.0_wr)) then
      v=0.0_wr
    endif
    if (val.lt.255.0_wr) then
      v=0.0_wr
    endif
    if ((val.gt.255.0_wr).and.(val.lt.(255.0_wr+255.0_wr/2.0_wr))) then
      v=real(int(val),wr)*2.0_wr-2.0_wr*255.0_wr
    endif
    if ((val.gt.(255.0_wr+255.0_wr/2.0_wr)).and.(val.lt.(255.0_wr*2.0_wr))) then
      v=255.0_wr*4.0_wr-real(int(val),wr)*2.0_wr
    endif
    ! -------                                                  --------    . bleu
    if (val.lt.(255.0_wr*2.0_wr)) then
      b=0.0_wr
    endif
    if ((val.lt.(255.0_wr*2.0_wr+255.0_wr/2.0_wr)).and.(val.gt.(255.0_wr*2.0_wr))) then
      b=real(int(val),wr)*2.0_wr-4.0_wr*255.0_wr
    endif
    if (val.gt.(255.0_wr*2.0_wr+255.0_wr/2.0_wr)) then
      b=255.0_wr*6.0_wr-real(int(val),wr)*2.0_wr
    endif
    write(color,'(i3.3,a1,i3.3,a1,i3.3)')int(r,wi),"/",int(v,wi),"/",int(b,wi)
    ! -----------------------------------------------------------------    .
  end subroutine RVB

END MODULE figure_GMTpar



! *********************************************************************    .
! *********************************************************************    .


