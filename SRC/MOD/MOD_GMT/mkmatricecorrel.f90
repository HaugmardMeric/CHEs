! permet la création des scripts GMT pour les figures
! octobre 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTmCorr

    use modparam

    implicit none

    private

    public  :: GMT_mCorr


CONTAINS

! ---------------------------------------------------------------------    .

subroutine GMT_mCorr(dp)
    ! -------                                                  --------    .mh
    ! figures de la matrice de corrélation
    ! -----------------------------------------------------------------    .
    use typetemps
    use statistiques
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    type(densityplot), intent (in) :: dp                                   ! modèles retenus par McMC
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,nbpar,Noldtime, Nnewtime, ratetime
    real(KIND=wr) :: rp, tl
    real(KIND=wr), dimension(:,:), allocatable :: tabCorrel
    character(LEN=15), dimension(:), allocatable :: namepar
    character (LEN=5) :: numberchaine
    ! -----------------------------------------------------------------    .
    nbpar=4+nbseismes*4
    ! -------                                                  --------    .
    allocate(tabCorrel(nbpar,dp%nbparam))
    allocate(namepar(nbpar))
    tabCorrel(1,:)=dp%VC%vec ; namepar(1)='VC'
    tabCorrel(2,:)=dp%VM%vec ; namepar(2)='VM'
    tabCorrel(3,:)=dp%Zmoho%vec ; namepar(3)='Zmoho'
    tabCorrel(4,:)=dp%VpVs%vec ; namepar(4)='VpVs'
    j=1
    do i=5,nbpar,4
      write(numberchaine(1:5),'(i5)')j
      tabCorrel(i,:)=dp%Lat(j)%vec ; namepar(i)='Lat_('//trim(adjustl(numberchaine))//")"
      tabCorrel(i+1,:)=dp%Lon(j)%vec ; namepar(i+1)='Lon_('//trim(adjustl(numberchaine))//")"
      tabCorrel(i+2,:)=dp%Zhypo(j)%vec ; namepar(i+2)='Zhypo_('//trim(adjustl(numberchaine))//")"
      tabCorrel(i+3,:)=dp%Tzero(j)%vec ; namepar(i+3)='Tzero_('//trim(adjustl(numberchaine))//")"
      j=j+1
    enddo
    ! -------                                                  --------    .
    open(unit=100,file="OUTPUT/GMT/matrice.corr",STATUS="replace")
    do i=1,nbpar
      do j=1,nbpar
        call Rpcalc(tabCorrel(i,:),tabCorrel(j,:),dp%nbparam,Rp)
        write(100,'(f6.1,1x,f6.1,1x,f9.6,1x,a,1x,a)') real(i,wr)-0.5_wr, real(j,wr)-0.5_wr ,Rp, &
          trim(adjustl(namepar(i))), trim(adjustl(namepar(j)))
      enddo
    enddo
    deallocate(tabCorrel)
    close(100)
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture des script GMT_matCorr "
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT matCorr'"
    write(600,*)"#########################################"
    write(600,*)"###########     autocorr       ##########"
    write(600,*)"#########################################"
    ! -----------------------------------------------------------------    .
    write(600,*)"geoproj=-JX5i"
    write(600,'(a,i9.9,a,i9.9)')" geozone=-R-0/",nbpar,"/0/",nbpar
    write(600,*)"file=OUTPUT/GMT/matCorr.ps"
    ! -------                                                  --------    .
    write(600,'(a)')"echo '-1.00 000 000 255 -0.25 255 255 255' > OUTPUT/GMT/colorpal6.cpt"
    write(600,'(a)')"echo '-0.25 255 255 255 0.250 255 255 255' >> OUTPUT/GMT/colorpal6.cpt"
    write(600,'(a)')"echo '0.250 255 255 255 1.000 255 000 000' >> OUTPUT/GMT/colorpal6.cpt"
    write(600,*)"cat OUTPUT/GMT/matrice.corr | awk '{print $1, $2, $3}' > OUTPUT/GMT/toto.d"
    write(600,'(a,i9.9,2a)')" head -",nbpar," OUTPUT/GMT/matrice.corr | awk ", &
        "'{print -1.0, $2, 4, 0, 4, ""RM"", $5}' > OUTPUT/GMT/toto2.d"
    write(600,'(a,i9.9,2a)')" head -",nbpar," OUTPUT/GMT/matrice.corr | awk ", &
        "'{print $2, -1.0, 4, -90, 4, ""LM"", $5}' > OUTPUT/GMT/toto3.d"
    ! -------                                                  --------    .
    write(600,*)"xyz2grd -I.5 $geozone -GOUTPUT/GMT/toto1.grd OUTPUT/GMT/toto.d"
    write(600,*)"grdsample OUTPUT/GMT/toto1.grd -I1 -GOUTPUT/GMT/toto.grd -F"
    write(600,*)"grdimage OUTPUT/GMT/toto.grd $geozone $geoproj -COUTPUT/GMT/colorpal6.cpt -K -Sn -Xc -Yc > $file"
    write(600,*)"pstext OUTPUT/GMT/toto2.d $geozone $geoproj -O -K -N -Ba0 >> $file"
    write(600,*)"pstext OUTPUT/GMT/toto3.d $geozone $geoproj -O -K -N -Ba0 >> $file"
    do i=0,nbpar,4
      write(600,*)"echo -e '",i,0,"\n",i,nbpar,"' | psxy $geozone $geoproj -O -K -W3 >> $file"
      write(600,*)"echo -e '",0,i,"\n",nbpar,i,"' | psxy $geozone $geoproj -O -K -W3 >> $file"
    enddo
    write(600,*)"psscale -D5.5i/2.5i/5i/0.2i -O -K -COUTPUT/GMT/colorpal6.cpt -B0.2:'corr\351lation': >> $file"
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone $geoproj -Ba0g1 -O >>  $file"
    write(600,*)"ps2raster OUTPUT/GMT/matCorr.ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/matCorr.pdf OUTPUT/figures/matCorr.pdf"
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

  end subroutine GMT_mCorr

    ! -----------------------------------------------------------------    .

END MODULE figure_GMTmCorr



! *********************************************************************    .
! *********************************************************************    .


