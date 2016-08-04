! permet la création des scripts GMT pour le diagramme de Wadati
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTwada

    use modparam

    implicit none

    private

    public  :: GMT_wadati


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_wadati(nbtps,D,param,dp)
    ! -------                                                  --------    .mh
    ! production d'une partie du script GMT pour le Wadatiplot
    ! -------                                                  --------    .
    use typetemps
    use cpt_temps
    use time
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(in) :: nbtps(nbseismes)
    type(dataall), intent(in) :: D(nbseismes)
    type(parametres), intent(in) :: param
    type(densityplot), intent (in) :: dp
    ! -------                                                  --------    .
    real (KIND=wr) :: aREF, R2REF                                          ! coeficient directeur de la régression
    real (KIND=wr) :: aDIR, R2DIR                                          ! Chi2
    real (KIND=wr) :: X, Y, tl
    real (KIND=wr) :: Xmaxi, Ymaxi
    ! -------                                                  --------    .
    integer (KIND=wi), parameter :: taille=2500
    type(date_sec) :: one_tps_1,one_tps_2
    integer(KIND=wi) :: i, j, n, nREF, nDIR, Noldtime, Nnewtime, ratetime
    real (KIND=wr) :: XY(taille,3),XYREF(taille,3),XYDIR(taille,3)
    real (KIND=wr) :: sDIR(taille,2),sREF(taille,2)
    real (KIND=wr) :: x1, y1
    real (KIND=wr) :: a, R2                                                ! coeficient directeur et Chi2 de la régression
    character(LEN=7) :: char1, char2
    character (LEN=5) :: char_nbseismes
    ! -----------------------------------------------------------------    . diagramme de wadati pour toutes les ondes
    call Wadatiplot(nbtps,D,param,a=a,R2=R2,XY=XY,nb=n)
    ! -------                                                  --------    . diagramme de wadati pour ondes réfractées
    call Wadatiplot(nbtps,D,param,a=aREF,R2=R2REF,XY=XYREF,nb=nREF,sig=sREF,atype='N')
    ! -------                                                  --------    . diagramme de wadati pour ondes directes
    call Wadatiplot(nbtps,D,param,a=aDIR,R2=R2DIR,XY=XYDIR,nb=nDIR,sig=sDIR,atype='G')
    ! -----------------------------------------------------------------    . minmax
    Xmaxi=1.0_wr
    Ymaxi=1.0_wr
    do i=1,n
      if (Xmaxi.lt.XY(i,1)) Xmaxi=XY(i,1)
      if (Ymaxi.lt.XY(i,2)) Ymaxi=XY(i,2)
    enddo
    Xmaxi=1.1_wr*Xmaxi
    Ymaxi=1.1_wr*Ymaxi
    ! -----------------------------------------------------------------    .
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"echo 'execution du script GMT Wadatiplot'"
    write(*,*)"ecriture du script GMT Wadatiplot "
    ! -----------------------------------------------------------------    .
    do j=1,nbseismes
      write(char_nbseismes(1:5),'(i5)')j
      write(600,*)"#########################################"
      write(600,*)"############### Wadatiplot ##############"
      ! ---------------------------------------------------------------    .
      write(600,*)"geoproj=-JX13i/8i"                                      ! système de projection
      write(600,'(a12,E13.7,a3,E13.7)')"geozone=-R0/",Xmaxi,"/0/",Ymaxi
      write(600,*)"file=OUTPUT/GMT/wadatiplot"//trim(adjustl(char_nbseismes))//".ps"
      if (Xmaxi.gt.60.0_wr) then
        write(600,*)"psbasemap $geozone $geoproj -Ba10f5:""T@-P@--T@-0@- (s)"":",&
        "/a10f5:""T@-S@--T@-P@- (s)"":WenS -Xc -Yc -K >  $file"
      else
        write(600,*)"psbasemap $geozone $geoproj -Ba2f.5:""T@-P@--T@-0@- (s)"":",&
        "/a2f.5:""T@-S@--T@-P@- (s)"":WenS -Xc -Yc -K >  $file"
      endif
      ! -------                                                --------    .
      do i=1,int(Xmaxi,wi)+1
        ! -------                                              --------    . traits théoriquea min et max
        write(600,*)"echo -e '",i-1,(dp%VpVs%themax-1.0_wr)*real(i-1,wr),"\n",i,(dp%VpVs%themax-1.0_wr)*real(i,wr), &
            "' | psxy $geozone $geoproj -W1,red,- -O -K >>  $file"
        write(600,*)"echo -e '",i-1,(dp%VpVs%themin-1.0_wr)*real(i-1,wr),"\n",i,(dp%VpVs%themin-1.0_wr)*real(i,wr), &
            "' | psxy $geozone $geoproj -W1,red,- -O -K >>  $file"
        ! -------                                              --------    . traits théoriquea (droite de regression) TOTAL
        write(600,*)"echo -e "" ",i-1,a*real(i-1,wr),"\n",i,a*real(i,wr), " \"
        write(600,*)" "" | psxy $geozone $geoproj -W5 -O -K >>  $file"
        ! -------                                              --------    . traits théoriquea (droite de regression) DIR
        write(600,*)"echo -e "" ",i-1,aDIR*real(i-1,wr),"\n",i,aDIR*real(i,wr), " \"
        write(600,*)" "" | psxy $geozone $geoproj -W5,gray,-.-. -O -K >>  $file"
        ! -------                                              --------    . traits théoriquea (droite de regression) REF
        write(600,*)"echo -e "" ",i-1,aREF*real(i-1,wr),"\n",i,aREF*real(i,wr), " \"
        write(600,*)" "" | psxy $geozone $geoproj -W5,gray,-- -O -K >>  $file"
        ! -------                                              --------    .
      enddo
      ! -------                                                --------    .
      do i=1,nDIR                                                          ! points réels des ondes directes
                                                                           !  modèle de référence : 1000 best modèle
        write(600,*)"echo",XYDIR(i,1),XYDIR(i,2),XYDIR(i,3), &
        sqrt(sDIR(i,1)**2.0_wr+dp%Tzero(j)%ec_1000**2.0_wr),sqrt(sDIR(i,1)**2.0_wr+sDIR(i,2)**2.0_wr), " \"
        write(600,*)" | psxy $geozone $geoproj -O -K -St0.1i -Wthinnest -Exy -COUTPUT/GMT/colorpal3.cpt >>  $file"
      enddo
      ! -------                                                --------    .
      do i=1,nREF                                                          ! points réels des ondes réfractées
                                                                           ! modèle de référence : 1000 best modèle
        write(600,*)"echo",XYREF(i,1),XYREF(i,2),XYREF(i,3), &
        sqrt(sREF(i,1)**2.0_wr+dp%Tzero(j)%ec_1000**2.0_wr),sqrt(sREF(i,1)**2.0_wr+sREF(i,2)**2.0_wr), " \"
        write(600,*)" | psxy $geozone $geoproj -O -K -Si0.1i -Wthinnest -Exy -COUTPUT/GMT/colorpal3.cpt >>  $file"
      enddo
      ! -------                                                --------    . POUR CE SEISME
      if (nbseismes.gt.1) then
        do i=1,nbtps(j)
          if(D(j)%datatps(i)%andS.eq.'S') then
            one_tps_1%date = D(j)%datatps(i)%tpsR%date
            one_tps_1%sec = D(j)%datatps(i)%tpsR%secP
            call basetime(one_tps_1)
            one_tps_2 = dp%temps_ref(j)
            one_tps_2%sec = dp%Tzero(j)%moy_1000                           !  modèle de référence : 1000 best modèle
            ! one_tps_2%sec = dp%Tzero(j)%best                             !  modèle de référence : best modèle
            call basetime(one_tps_2)
            call difftime(x1,one_tps_1,one_tps_2)
            one_tps_2%date = D(j)%datatps(i)%tpsR%date
            one_tps_2%sec = D(j)%datatps(i)%tpsR%secS
            call basetime(one_tps_2)
            call difftime(y1,one_tps_2,one_tps_1)
            write(600,*)"echo",x1,y1," | psxy $geozone $geoproj -O -K -Sc0.3i -Wthinnest >>  $file"
          endif
        enddo
      endif
      ! -------                                                --------    . légende
      write(600,*)"#########################################"
      X = Xmaxi * 0.1_wr
      Y = Ymaxi * 0.8_wr
      if(XYREF(3,1).gt.0.0_wr) then
        ! si il existe des ondes réfractées :
        write(600,*)"echo",X,Y," | psxy $geozone $geoproj -O -K -St0.1i -Wthinnest -Gyellow >>  $file"
        write(600,*)"echo -e """,X-0.25_wr*X,Y," \n", X-0.5_wr*X,Y,""" | psxy $geozone $geoproj -O -K -W5,gray,-.-. >>  $file"
        write(char1,'(f7.4)')1.0_wr+aDIR
        write(char2,'(f7.4)')R2DIR
        write(600,*)"echo """,X+0.1_wr*X,Y,"15 0 4 LM ondes directes : V@-P@-/V@-S@- =",char1," \"
        write(600,*)" (@~\143@~@-2@- =",char2,")"," \"
        write(600,*)""" | pstext $geozone $geoproj -O -K >> $file"
        Y = Ymaxi * 0.75_wr
        write(char1,'(f7.4)')1.0_wr+aREF
        write(char2,'(f7.4)')R2REF
        write(600,*)"echo",X,Y," | psxy $geozone $geoproj -O -K -Si0.1i -Wthinnest -Gyellow >>  $file"
        write(600,*)"echo -e """,X-0.25_wr*X,Y," \n", X-0.5_wr*X,Y,""" | psxy $geozone $geoproj -O -K -W5,gray,-- >>  $file"
        write(600,*)"echo """,X+0.1_wr*X,Y,"15 0 4 LM ondes r\351fract\351es : V@-P@-/V@-S@- =",char1," \"
        write(600,*)" (@~\143@~@-2@- =",char2,")"" | pstext $geozone $geoproj -O -K >> $file"
        Y = Ymaxi * 0.70_wr
        write(char1,'(f7.4)')1.0_wr+a
        write(char2,'(f7.4)')R2
        write(600,*)"echo -e """,X-0.25_wr*X,Y," \n", X-0.5_wr*X,Y,""" | psxy $geozone $geoproj -O -K -W5 >>  $file"
        write(600,*)"echo """,X+0.1_wr*X,Y,"15 0 4 LM ensembles : V@-P@-/V@-S@- ="," \"
        write(600,*)char1," (@~\143@~@-2@- =",char2,")"" | pstext $geozone $geoproj -O -K >> $file"
      else
        Y = Ymaxi * 0.75_wr
        write(char1,'(f7.4)')1.0_wr+a
        write(char2,'(f7.4)')R2
        write(600,*)"echo -e """,X-0.25_wr*X,Y," \n", X-0.5_wr*X,Y,""" | psxy $geozone $geoproj -O -K -W5 >>  $file"
        write(600,*)"echo """,X,Y,"15 0 4 LM ondes directes : V@-P@-/V@-S@- =",char1," \"
        write(600,*)" (@~\143@~@-2@- =",char2,")"" | pstext $geozone $geoproj -O -K >> $file"
      endif
      ! -------                                                  --------    .
      write(600,*)"psscale -D1/-1/0.50E+01/0.25ch -B.25:""pond\351ration"": -S -I -COUTPUT/GMT/colorpal3.cpt -O -K >> $file"
      ! -------                                                  --------    .
      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
      write(600,*)"#########################################"
      write(600,*)"ps2raster OUTPUT/GMT/wadatiplot"//trim(adjustl(char_nbseismes))//".ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/wadatiplot"//trim(adjustl(char_nbseismes))//".pdf ", &
        "OUTPUT/figures/wadatiplot"//trim(adjustl(char_nbseismes))//".pdf"
      write(600,*)"#########################################"
      write(600,*)"ELAPSED=$(($SECONDS-$BEFORE))"
      write(600,*)" echo $ELAPSED secondes"
      call system_clock(Nnewtime,ratetime)
      tl=(real(Nnewtime,wr)-real(Noldtime,wr))/real(ratetime,wr)
      write(*,'(a9,i2.2,'':'',i2.2,'':'',f9.2)')' temps : ',int(tl/3600.0_wr,wi), &
      int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi),(tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr)
      ! ---------------------------------------------------------------    .
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine GMT_wadati

END MODULE figure_GMTwada



! *********************************************************************    .
! *********************************************************************    .


