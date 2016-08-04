! permet la création des scripts GMT pour les figures sur le fonction coût
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTfc

    use modparam

    implicit none

    private

    public  :: GMT_fc


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_fc(dp,nbChaineMV,nmod)
    ! -------                                                  --------    .mh
    use typetemps, only : densityplot
    use figure_GMTpar, only : RVB
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    type(densityplot), intent (in) :: dp                                       ! modèles retenus par McMC
    integer(KIND=wi), intent (in) :: nbChaineMV
    integer(KIND=wi), intent (in) :: nmod(nbChaineMV)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k,maxiter,Noldtime, Nnewtime, ratetime
    character(LEN=5) :: char
    real(KIND=wr) :: tl
    character(LEN=11) :: color
    ! -----------------------------------------------------------------    . ecriture des misfit par chaîne
    k=0
    maxiter=0
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      open(unit=850+i,file="OUTPUT/GMT/themis"//char//".bin",STATUS="replace",access='direct',RECL=16)
      do j=1,nmod(i)
        k=k+1
        write(850+i,rec=j)real(dp%mis%vec(k),8),real(j,8)
        if(maxiter.lt.j)maxiter=j
      enddo
      close(850+i)
    enddo
    ! -----------------------------------------------------------------    . script GMT
    write(*,*)"ecriture des script GMT_fc "
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT fonction coût'"
    write(600,*)"#########################################"
    write(600,*)"###########  fonction coût     ##########"
    write(600,*)"#########################################"
    ! -----------------------------------------------------------------    .
    write(600,'(a,E13.7,a1,E13.7)')"geozone=-R0/7.5/",dp%mis%themin-2.0_wr,"/",dp%mis%themax+2.0_wr
    write(600,*)"geoproj=-JX2i/4.5i"
    write(600,*)"file=OUTPUT/GMT/mishisto.ps"
    ! -----------------------------------------------------------------    . histogramme de la fonction cout
    write(600,*)"psbasemap $geozone $geoproj -Ba5f1:'effectif (%)':/a10Snew -K -X10i >  $file"
    write(600,'(a,E13.7,a1,E13.7,a)')"geozone=-R",dp%mis%themin-2.0_wr,"/",dp%mis%themax+2.0_wr,"/0/10"
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/themis"//char//".bin -bi1d -W0.2", &
      " -Ggray -Z1 -O -K -A >>  $file"
    enddo
    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/themis"//char//".bin -bi2d", &
      " -W0.2 -S -L0/0 -Z1 -O -K -A >>  $file"
    enddo
    write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/lon_lat-1_tot.bin -bi3d -T2 -W0.2", &
    " -Gblue -Z1 -K -O -A  >>  $file"
    write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/lon_lat-1_tot.bin -bi3d -T2 -W0.2 -S", &
    "  -L0/0 -Z1 -O -K -A >>  $file"
    ! -----------------------------------------------------------------    .
    write(600,'(a,i9.9,a,E13.7,a,E13.7)')"geozone=-R0/",maxiter,"/",dp%mis%themin-2.0_wr,"/",dp%mis%themax+2.0_wr
    write(600,*)"geoproj=-JX8i/4.5i"
    write(600,'(a31,i9.9,a)')"psbasemap $geozone $geoproj -Ba",(maxiter/1000)*100, &
    ":'mod\350les':/a10:'fonction co\373t':nSeW -K -O -X-8.5i >>  $file"

    do i=1,nbChaineMV
      write(char,'(i5)')10000+i
      call RVB(i,nbChaineMV,color)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/themis"//char//".bin -bi2d -Wthinnest,"//color//" -O -K -: >>  $file"
    enddo
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    ! -----------------------------------------------------------------    .
    write(600,*)"ps2raster OUTPUT/GMT/mishisto.ps -Tf -A"
    write(600,'(a)')"mv OUTPUT/GMT/mishisto.pdf OUTPUT/figures/mishisto.pdf"
    ! -------                                                  --------    . fin script GMT
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
  end subroutine GMT_fc

END MODULE figure_GMTfc



! *********************************************************************    .
! *********************************************************************    .


