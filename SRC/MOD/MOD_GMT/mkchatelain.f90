! permet la création des scripts GMT pour le diagramme de Wadati modifié (châtelain, 1978)
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTchat

    use modparam

    implicit none

    private

    public  :: GMT_chatelain


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_chatelain(nbtps,nbsta,D,dp)
    ! -------                                                  --------    .mh
    ! production d'une partie du script GMT pour le diagramme de Châtelain (Châtelain, 1978)
    ! le diagramme de Châtelain (ou Wadati modifié) est indépendant des parametres
    ! et peux ainsi identifier une erreur de pointé sur les ondes          !
    ! -------                                                  --------    .
    use typetemps
    use cpt_temps
    use time
    use pb_direct
    use tri
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent(in) :: nbtps(nbseismes),nbsta
    type(dataall), intent(in) :: D(nbseismes)
    type(densityplot), intent (in) :: dp
    ! -------                                                  --------    .
    real (KIND=wr) :: aREF, R2REF                                          ! coeficient directeur de la régression
    real (KIND=wr) :: aDIR, R2DIR                                          ! Chi2
    real (KIND=wr) :: X, Y, tl
    real (KIND=wr) :: Xmaxi, Ymaxi
    ! -------                                                  --------    .
    integer (KIND=wi), parameter :: taille=5000
    ! -------                                                  --------    .
    type(date_sec) :: one_tps_1,one_tps_2
    integer(KIND=wi) :: i, j, k, l
    integer(KIND=wi) :: n, nREF, nDIR
    integer(KIND=wi) :: Noldtime, Nnewtime, ratetime
    integer(KIND=wi) :: test1,test2
    real (KIND=wr) :: XY(taille,3),XYREF(taille,3),XYDIR(taille,3)
    real (KIND=wr) :: sDIR(taille,2),sREF(taille,2)
    real (KIND=wr) :: x1, y1, sta_dist(nbsta+2)
    real (KIND=wr) :: a, R2                                                ! coeficient directeur et Chi2 de la régression
    real (KIND=wr) :: nsta(nbsta+2)
    character(len=4) :: nomstaDIR(taille,2),nomstaREF(taille,2)
    character(LEN=4) :: sta(nbsta+2)
    character(LEN=7) :: char1, char2
    character (LEN=5) :: char_nbseismes
    ! -----------------------------------------------------------------    . diagramme de chatelainplot pour toutes les ondes
    call chatelainplot(nbtps,D,a=a,R2=R2,XY=XY,nb=n)
    ! -------                                                  --------    . diagramme de chatelainplot pour ondes réfractées
    call chatelainplot(nbtps,D,a=aREF,R2=R2REF,XY=XYREF,nb=nREF,sig=sREF,atype='N',nom_sta=nomstaREF)
    ! -------                                                  --------    . diagramme de chatelainplot pour ondes directes
    call chatelainplot(nbtps,D,a=aDIR,R2=R2DIR,XY=XYDIR,nb=nDIR,sig=sDIR,atype='G',nom_sta=nomstaDIR)
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
    ! -----------------------------------------------------------------    .
    ! ecart des stations à la courbe théorique
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    . DIR :
    ! -----------------------------------------------------------------    .
    sta_dist(:)=0.0_wr
    sta(:)='    '
    nsta(:)=0.0_wr
    l=1
    ! -------                                                  --------    .
    do i=1,n
      if (nomstaDIR(i,1).ne.nomstaDIR(i,2)) then
        ! -------                                              --------    .
        test1=0
        test2=0
        ! -------                                              --------    .
        do k=1,l
          ! -------                                            --------    . si la station existe déjà
          if (sta(k)==nomstaDIR(i,1)) then
            call distancePlot(aDIR,XYDIR(i,1),XYDIR(i,2),sta_dist(k))
            nsta(k)=nsta(k)+1.0_wr
            test1=test1+1
          endif
          ! -------                                            --------    .
          if (sta(k)==nomstaDIR(i,2)) then
            call distancePlot(aDIR,XYDIR(i,1),XYDIR(i,2),sta_dist(k))
            nsta(k)=nsta(k)+1.0_wr
            test2=test2+1
          endif
          ! -------                                            --------    .
        enddo
        ! -------                                              --------    . si la station n'existe pas encore
        if (test1==0) then
          sta(l)=nomstaDIR(i,1)
          call distancePlot(aDIR,XYDIR(i,1),XYDIR(i,2),sta_dist(l))
          nsta(l)=1.0_wr
          l=l+1
        elseif(test1.gt.1) then
          write(*,*)'problème dans GMT_chatelain 1 : test1 =',test1
          stop
        endif
        ! -------                                              --------    .
        if (test2==0) then
          sta(l)=nomstaDIR(i,2)
          call distancePlot(aDIR,XYDIR(i,1),XYDIR(i,2),sta_dist(l))
          nsta(l)=1.0_wr
          l=l+1
        elseif(test1.gt.1) then
          write(*,*)'problème dans GMT_chatelain 1 : test2 =',test2
          stop
        endif
        ! -------                                              --------    .
        if (l.gt.nbsta+2) then
            write(*,*)'problème dans GMT_chatelain 1 : l > nbsta : ',l,nbsta
            stop
        endif
        ! -------                                              --------    .
      endif
    enddo
    ! -------                                                  --------    . moyenne
    do i=1,nbsta
      if (nsta(i).gt.0.0_wr) sta_dist(i)=sta_dist(i)/nsta(i)
    enddo
    ! -------                                                  --------    . tri
    call tri_bulle(sta_dist,sta,nbsta+2)
    ! -------                                                  --------    .
    X = Xmaxi * 0.75_wr
    Y = Ymaxi * 0.075_wr
    open(unit=33,file="OUTPUT/GMT/chat-stationsDIR.d",status='replace')
    do i=1,nbsta
      Y = Y + Ymaxi * 0.05_wr
      if ((sta_dist(i).gt.0.5_wr).and.(i.le.5)) then
        write(33,'(2f10.5,a,f8.2,a)')X,Y," 15 0 5 6 @~\104@~ "//sta(i)//" g = ",sta_dist(i)," s"
      else
        write(33,'(a,f6.2,a)')"-1000 -1000 15 0 5 6 @~\104@~ "//sta(i)//" g = ",sta_dist(i)," s"
      endif
    enddo
    close(33)
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    . REF :
    ! -----------------------------------------------------------------    .
    sta_dist(:)=0.0_wr
    sta(:)='    '
    nsta(:)=0.0_wr
    l=1
    ! -------                                                  --------    .
    do i=1,n
      if (nomstaREF(i,1).ne.nomstaREF(i,2)) then
        ! -------                                              --------    .
        test1=0
        test2=0
        ! -------                                              --------    .
        do k=1,l
          ! -------                                            --------    . si la station existe déjà
          if (sta(k)==nomstaREF(i,1)) then
            call distancePlot(aREF,XYREF(i,1),XYREF(i,2),sta_dist(k))
            nsta(k)=nsta(k)+1.0_wr
            test1=test1+1
          endif
          ! -------                                            --------    .
          if (sta(k)==nomstaREF(i,2)) then
            call distancePlot(aREF,XYREF(i,1),XYREF(i,2),sta_dist(k))
             nsta(k)=nsta(k)+1.0_wr
            test2=test2+1
          endif
          ! -------                                            --------    .
        enddo
        ! -------                                              --------    . si la station n'existe pas encore
        if (test1==0) then
          sta(l)=nomstaREF(i,1)
          call distancePlot(aREF,XYREF(i,1),XYREF(i,2),sta_dist(l))
          nsta(l)=1.0_wr
          l=l+1
        elseif(test1.gt.1) then
          write(*,*)'problème dans GMT_chatelain 2 : test1 =',test1
          stop
        endif
        ! -------                                              --------    .
        if (test2==0) then
          sta(l)=nomstaREF(i,2)
          call distancePlot(aREF,XYREF(i,1),XYREF(i,2),sta_dist(l))
          nsta(l)=1.0_wr
          l=l+1
        elseif(test1.gt.1) then
          write(*,*)'problème dans GMT_chatelain 2 : test2 =',test2
          stop
        endif
        ! -------                                              --------    .
        if (l.gt.nbsta+2) then
          write(*,*)'problème dans GMT_chatelain 2 : l > nbsta : ',l,nbsta
          stop
        endif
        ! -------                                              --------    .
      endif
    enddo
    ! -------                                                  --------    . moyenne
    do i=1,nbsta
      if (nsta(i).gt.0.0_wr) sta_dist(i)=sta_dist(i)/nsta(i)
    enddo
    ! -------                                                  --------    . tri
    call tri_bulle(sta_dist,sta,nbsta+2)
    ! -------                                                  --------    .
    X = Xmaxi * 0.9_wr
    Y = Ymaxi * 0.075_wr
    open(unit=33,file="OUTPUT/GMT/chat-stationsREF.d",status='replace')
    do i=1,nbsta
      Y = Y + Ymaxi * 0.05_wr
      if ((sta_dist(i).gt.0.5_wr).and.(i.le.5)) then
        write(33,'(2f10.5,a,f8.2,a)')X,Y," 15 0 5 6 @~\104@~ "//sta(i)//" n = ",sta_dist(i)," s"
      else
        write(33,'(a,f6.2,a)')"-1000 -1000 15 0 5 6 @~\104@~ "//sta(i)//" n = ",sta_dist(i)," s"
      endif
    enddo
    close(33)
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    . SCRIPT GMT :
    ! -----------------------------------------------------------------    .
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"echo 'execution du script GMT chatelainplot'"
    write(*,*)"ecriture du script GMT chatelainplot "
    ! -----------------------------------------------------------------    .
    do j=1,nbseismes
      write(char_nbseismes(1:5),'(i5)')j
      write(600,*)"#########################################"
      write(600,*)"############## chatelainplot #############"
      ! ---------------------------------------------------------------    .
      write(600,*)"geoproj=-JX13i/8i"                                      ! système de projection
      write(600,'(a12,E13.7,a3,E13.7)')"geozone=-R0/",Xmaxi,"/0/",Ymaxi
      write(600,*)"file=OUTPUT/GMT/chatelainplot"//trim(adjustl(char_nbseismes))//".ps"
      if (Xmaxi.gt.60.0_wr) then
        write(600,*)"psbasemap $geozone $geoproj -Ba10f5:""T@-P1@--T@-P2@- (s)"":",&
        "/a10f5:""T@-S1@--T@-S2@- (s)"":WenS -Xc -Yc -K >  $file"
      else
        write(600,*)"psbasemap $geozone $geoproj -Ba2f.5:""T@-P1@--T@-P2@- (s)"":",&
        "/a2f.5:""T@-S1@--T@-S2@- (s)"":WenS -Xc -Yc -K >  $file"
      endif
      ! -------                                                --------    .
      do i=1,int(Xmaxi,wi)+1
        ! -------                                              --------    . traits théoriquea min et max
        write(600,*)"echo -e '",i-1,(dp%VpVs%themax)*real(i-1,wr),"\n",i,(dp%VpVs%themax)*real(i,wr), &
            "' | psxy $geozone $geoproj -W1,red,- -O -K >>  $file"
        write(600,*)"echo -e '",i-1,(dp%VpVs%themin)*real(i-1,wr),"\n",i,(dp%VpVs%themin)*real(i,wr), &
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
      open(unit=30,file="OUTPUT/GMT/chat-dir"//trim(adjustl(char_nbseismes))//".txt",status='replace')
      do i=1,nDIR                                                          ! points réels des ondes directes
        write(30,*)XYDIR(i,1),XYDIR(i,2),XYDIR(i,3),sDIR(i,1),sDIR(i,2)
      enddo
      close(30)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/chat-dir"//trim(adjustl(char_nbseismes))//".txt -O -K -St0.1i ", &
        "-Wthinnest -Exy -COUTPUT/GMT/colorpal3.cpt >>  $file"
      ! -------                                                --------    .
      open(unit=31,file="OUTPUT/GMT/chat-ref"//trim(adjustl(char_nbseismes))//".txt",status='replace')
      do i=1,nREF                                                          ! points réels des ondes réfractées
        write(31,*)XYREF(i,1),XYREF(i,2),XYREF(i,3),sREF(i,1),sREF(i,2)
      enddo
      close(31)
      write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/chat-ref"//trim(adjustl(char_nbseismes))//".txt -O -K -Si0.1i ", &
        "-Wthinnest -Exy -COUTPUT/GMT/colorpal3.cpt >>  $file"
      ! -------                                                --------    . POUR CE SEISME
      open(unit=32,file="OUTPUT/GMT/chat-all"//trim(adjustl(char_nbseismes))//".txt",status='replace')
      if (nbseismes.gt.1) then
        do i=1,nbtps(j)
          do k=1,nbtps(j)
            if ((D(j)%datatps(i)%andS=='S').and.(D(j)%datatps(k)%andS=='S') &
              .and.(D(j)%datatps(i)%typeonde==D(j)%datatps(k)%typeonde)) then
              one_tps_1%date = D(j)%datatps(i)%tpsR%date
              one_tps_1%sec = D(j)%datatps(i)%tpsR%secP
              one_tps_2%date = D(j)%datatps(k)%tpsR%date
              one_tps_2%sec = D(j)%datatps(k)%tpsR%secP
              call difftime(x1,one_tps_1,one_tps_2)
              one_tps_1%date = D(j)%datatps(i)%tpsR%date
              one_tps_1%sec = D(j)%datatps(i)%tpsR%secS
              one_tps_2%date = D(j)%datatps(k)%tpsR%date
              one_tps_2%sec = D(j)%datatps(k)%tpsR%secS
              call difftime(y1,one_tps_1,one_tps_2)
              write(32,*)abs(x1),abs(y1)
            endif
          enddo
        enddo
        write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/chat-all"//trim(adjustl(char_nbseismes))//".txt ", &
          " -O -K -Sc0.3i -Wthinnest >>  $file"
        close(32)
      endif
      ! ---------------------------------------------------------------    .
      write(600,*)"pstext OUTPUT/GMT/chat-stationsDIR.d $geozone $geoproj -O -K -Gorange >> $file"
      write(600,*)"pstext OUTPUT/GMT/chat-stationsREF.d $geozone $geoproj -O -K -Gorange >> $file"
      ! ---------------------------------------------------------------    . légende
      write(600,*)"#########################################"
      X = Xmaxi * 0.1_wr
      Y = Ymaxi * 0.8_wr
      if(XYREF(3,1).gt.0.0_wr) then
        ! si il existe des ondes réfractées :
        write(600,*)"echo",X,Y," | psxy $geozone $geoproj -O -K -St0.1i -Wthinnest -Gyellow >>  $file"
        write(600,*)"echo -e """,X-0.25_wr*X,Y," \n", X-0.5_wr*X,Y,""" | psxy $geozone $geoproj -O -K -W5,gray,-.-. >>  $file"
        write(char1,'(f7.4)')aDIR
        write(char2,'(f7.4)')R2DIR
        write(600,*)"echo """,X+0.1_wr*X,Y,"15 0 4 LM ondes directes : V@-P@-/V@-S@- =",char1," \"
        write(600,*)" (@~\143@~@-2@- =",char2,")"," \"
        write(600,*)""" | pstext $geozone $geoproj -O -K >> $file"
        Y = Ymaxi * 0.75_wr
        write(char1,'(f7.4)')aREF
        write(char2,'(f7.4)')R2REF
        write(600,*)"echo",X,Y," | psxy $geozone $geoproj -O -K -Si0.1i -Wthinnest -Gyellow >>  $file"
        write(600,*)"echo -e """,X-0.25_wr*X,Y," \n", X-0.5_wr*X,Y,""" | psxy $geozone $geoproj -O -K -W5,gray,-- >>  $file"
        write(600,*)"echo """,X+0.1_wr*X,Y,"15 0 4 LM ondes r\351fract\351es : V@-P@-/V@-S@- =",char1," \"
        write(600,*)" (@~\143@~@-2@- =",char2,")"" | pstext $geozone $geoproj -O -K >> $file"
        Y = Ymaxi * 0.70_wr
        write(char1,'(f7.4)')a
        write(char2,'(f7.4)')R2
        write(600,*)"echo -e """,X-0.25_wr*X,Y," \n", X-0.5_wr*X,Y,""" | psxy $geozone $geoproj -O -K -W5 >>  $file"
        write(600,*)"echo """,X+0.1_wr*X,Y,"15 0 4 LM ensembles : V@-P@-/V@-S@- ="," \"
        write(600,*)char1," (@~\143@~@-2@- =",char2,")"" | pstext $geozone $geoproj -O -K >> $file"
      else
        Y = Ymaxi * 0.75_wr
        write(char1,'(f7.4)')a
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
      write(600,*)"ps2raster OUTPUT/GMT/chatelainplot"//trim(adjustl(char_nbseismes))//".ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/chatelainplot"//trim(adjustl(char_nbseismes))//".pdf ", &
      "OUTPUT/figures/chatelainplot"//trim(adjustl(char_nbseismes))//".pdf"
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
  end subroutine GMT_chatelain

  ! -------------------------------------------------------------------    .

  subroutine distancePlot(a,x,y,d)
    ! -------                                                  --------    .mh
    ! 1) calcule de d, distance la plus courte entre un point (x,y)
    ! et la droite de coef directeur a
    ! ou
    ! 2) calcul de d, distance horizontale ou verticale maximale
    ! ou
    ! 3) calcul de d, distance
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    real (KIND=wr), intent(in) :: a                                        ! coeficient directeur
    real (KIND=wr), intent(in) :: x,y                                      ! point de coordonnées (x,y)
    real (KIND=wr), intent(inout) :: d
    ! -------                                                  --------    .
    real (KIND=wr) :: deltaX,deltaY
    ! -------                                                  --------    .
    deltaX=abs(y/a-x)
    deltaY=abs(y-a*x)
    ! -------                                                  --------    .
    ! d=d+deltaY*sin(atan(deltaX/deltaY))                                  ! 1)
    ! d=d+max(deltaX,deltaY)                                               ! 2)
    !d=d+sqrt(deltaX*deltaX+deltaY*deltaY)                                 ! 3)
    ! -------                                                  --------    .
    d=d+max(deltaX,deltaY)
    ! -----------------------------------------------------------------    .
  end subroutine distancePlot

END MODULE figure_GMTchat



! *********************************************************************    .
! *********************************************************************    .


