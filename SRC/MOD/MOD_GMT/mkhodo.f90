! permet la création des scripts GMT pour l'hodochrone
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMThodo

    use modparam

    implicit none


    private

    public  :: GMT_Hodochrone


    ! -----------------------------------------------------------------    .
    ! on défini _FILE_DIR_ en fonction du compilateur
    ! variable permettant de tester l'existance d'un dossier

#ifdef __INTEL_COMPILER
#define _FILE_DIR_ DIRECTORY
#elif __GFORTRAN__
#define _FILE_DIR_ FILE
#endif

    ! -----------------------------------------------------------------    .


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_Hodochrone(j,nbtps,datatps,param_best,xmaxcercle,acentroid)
    ! -------                                                  --------    .mh
    ! Calcul les regressions sur les hodochrones et affiche l'hodochrone
    ! -------                                                  --------    .
    use typetemps
    use time
    use statistiques
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: j
    integer(KIND=wi), intent(in) :: nbtps
    type(dataone), intent(in) :: datatps(nbtps)
    type(parametre), intent(inout) :: param_best
    real(KIND=wr), intent (in) :: xmaxcercle
    type(amoho_centroid), intent (in) :: acentroid
    ! -------                                                  --------    .
    type(date_sec) :: one_tps
    type(stations) :: a_sta
    integer (KIND=wi) :: i,k,l
    integer (KIND=wi) :: nPg, nPn ,nSg, nSn, Noldtime, Nnewtime, ratetime
    real (KIND=wr), dimension(:,:), allocatable :: pt_Pg, pt_Pn ,pt_Sg, pt_Sn
    real (KIND=wr) :: a_Pg,R2_Pg,a_Pn,b_Pn,R2_Pn,a_Sg,R2_Sg,a_Sn,b_Sn,R2_Sn
    real (KIND=wr) :: min, max_0, max, discritiqueH
    real (KIND=wr) :: X, Y, tl, val
    real (KIND=wr) :: Tpsmin,Tppmin,Tpsmax,Tppmax
    character (LEN=5) :: numberfile
    character (LEN=4) :: nomstadoublets(nbtps+1)
    logical :: deja,existe1
    ! -----------------------------------------------------------------    .
    nPg = 0                                                                ! nombre de données
    nPn = 0
    nSg = 0
    nSn = 0
    do i=1,nbtps
      if(datatps(i)%typeonde.eq.'N') then
        nPn = nPn + 1
      elseif(datatps(i)%typeonde.eq.'G') then
        nPg = nPg + 1
      else
        write(*,*)'problème dans GMT_Hodochrone : onde P ni N ni G'
        stop
      endif
      if(datatps(i)%andS.eq.'S') then
        if(datatps(i)%typeonde.eq.'N') then
          nSn = nSn + 1
        elseif(datatps(i)%typeonde.eq.'G') then
          nSg = nSg + 1
        else
          write(*,*)'problème dans GMT_Hodochrone : onde S ni N ni G'
          stop
        endif
      endif
    enddo
    ! -------                                                  --------    . pour chaque vecteur
    allocate (pt_Pg(nPg,4))
    allocate (pt_Pn(nPn,4))
    allocate (pt_Sg(nSg,4))
    allocate (pt_Sn(nSn,4))
    nPg = 0                                                                ! nombre de données
    nPn = 0
    nSg = 0
    nSn = 0
    max = 0.0_wr
    discritiqueH = 0.0_wr
    do i=1,nbtps
      one_tps%date = datatps(i)%tpsR%date
      one_tps%sec = datatps(i)%tpsR%secP
      if(datatps(i)%typeonde.eq.'N') then
        nPn = nPn + 1
        call difftime(pt_Pn(nPn,2),one_tps,param_best%Tzero)
        pt_Pn(nPn,1) = datatps(i)%dhypo
        pt_Pn(nPn,3) = datatps(i)%wp
        if((pt_Pn(nPn,2).gt.max).and.(datatps(i)%depi.lt.xmaxcercle/2.0_wr)) max =pt_Pn(nPn,2)
        discritiqueH = discritiqueH + datatps(i)%dcritiqueH
        pt_Pn(nPn,4) = datatps(i)%sigP
      elseif(datatps(i)%typeonde.eq.'G') then
        nPg = nPg + 1
        call difftime(pt_Pg(nPg,2),one_tps,param_best%Tzero)
        pt_Pg(nPg,1) = datatps(i)%dhypo
        pt_Pg(nPg,3) = datatps(i)%wp
        if((pt_Pg(nPg,2).gt.max).and.(datatps(i)%depi.lt.xmaxcercle/2.0_wr)) max = pt_Pg(nPg,2)
        pt_Pg(nPg,4) = datatps(i)%sigP
      else
        write(*,*)'problème dans GMT_Hodochrone : onde P ni N ni G'
        stop
      endif
      if(datatps(i)%andS.eq.'S') then
        if(datatps(i)%typeonde.eq.'N') then
          nSn = nSn + 1
          one_tps%date = datatps(i)%tpsR%date
          one_tps%sec = datatps(i)%tpsR%secS
          call difftime(pt_Sn(nSn,2),one_tps,param_best%Tzero)
          pt_Sn(nSn,1) = datatps(i)%dhypo
          pt_Sn(nSn,3) = datatps(i)%ws
          if((pt_Sn(nSn,2).gt.max).and.(datatps(i)%depi.lt.xmaxcercle/2.0_wr)) max =pt_Sn(nSn,2)
          discritiqueH = discritiqueH + datatps(i)%dcritiqueH
          pt_Sn(nSn,4) = datatps(i)%sigS
        elseif(datatps(i)%typeonde.eq.'G') then
          nSg = nSg + 1
          one_tps%date = datatps(i)%tpsR%date
          one_tps%sec = datatps(i)%tpsR%secS
          call difftime(pt_Sg(nSg,2),one_tps,param_best%Tzero)
          pt_Sg(nSg,1) = datatps(i)%dhypo
          pt_Sg(nSg,3) = datatps(i)%ws
          pt_Sg(nSg,4) = datatps(i)%sigS
          if((pt_Sg(nSg,2).gt.max).and.(datatps(i)%depi.lt.xmaxcercle/2.0_wr)) max = pt_Sg(nSg,2)
        else
          write(*,*)'problème dans GMT_Hodochrone : onde S ni N ni G'
          stop
        endif
      endif
    enddo
    discritiqueH = discritiqueH / real(nSn+nPn,wr)
    ! -------                                                  --------    .
    call correlationaffpond(a_Pg,R2_Pg,nPg,pt_Pg)
    call correlationpond(a_Pn,b_Pn,R2_Pn,nPn,pt_Pn)
    call correlationaffpond(a_Sg,R2_Sg,nSg,pt_Sg)
    call correlationpond(a_Sn,b_Sn,R2_Sn,nSn,pt_Sn)
    ! -------                                                  --------    .
    max_0 = max
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture du script GMT Hodochrone"
    write(600,*)"thegray1=200/200/200"
    write(600,*)"thegray2=175/175/175"
    write(600,*)"thegray3=150/150/150"
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#########################################"
    write(600,*)"############### Hodochrone ##############"
    write(600,*)"echo 'execution du script GMT Hodochrone'"
    write(numberfile(1:5),'(i5)')j
    write(600,*)"file=OUTPUT/GMT/hodochrone"//"-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"geoproj=-JX13i/8i"                                        ! système de projection
    if (xmaxcercle/2.0_wr.lt.datatps(nbtps)%dhypo) then
      val=xmaxcercle/2.0_wr
    else
      val=datatps(nbtps)%dhypo
    endif
    write(600,'(a12,E13.7,a3,E13.7)')"geozone=-R0/",val*1.2_wr,"/0/",max_0*1.1_wr
    if (max_0.gt.60.0_wr) then
      write(600,*)"psbasemap $geozone $geoproj -Ba50f25:""distance hypocentrale (km)"":",&
      "/a15f5g60:""temps d'arriv\351es des ondes (s)"":WenS -Xc -Yc -K >  $file"
    else
      write(600,*)"psbasemap $geozone $geoproj -Ba20f5:""distance hypocentrale (km)"":",&
      "/a5f1g60:""temps d'arriv\351es des ondes (s)"":WenS -Xc -Yc -K >  $file"
    endif
    ! -------                                                  --------    . plot distance critique
    write(600,*)"echo -e "" ",discritiqueH,1,"\n",discritiqueH,max_0*1.09_wr," \"
    write(600,*)" "" | psxy $geozone $geoproj -W0.01i,gray -O -K >>  $file"

    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    inquire (_FILE_DIR_="DATA/sac-"//trim(adjustl(numberfile)),exist=existe1) ! option différente selon compilo !
    if ((existe1).and.(tracessac)) then
      ! ---------------------------------------------------------------    .
      ! double si Pg et Pn (ou Sg et Sn), mais pas tres grave ....
      ! ---------------------------------------------------------------    . plot traces si existes
      do i=1,nbtps
        ! lecture un peu archaïque, mais permet un peu de souplesse dans le non des station et des fichiers sac
        ! -------                                              --------    . COMPOSANTE Z
        if (datatps(i)%sta%staname(4:4)=='0') then                         ! nom station en trois caractere + "0"
          write(600,'(9a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
            " 2>/dev/null | uniq | while read nom "
        else
          write(600,'(5a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//".*Z.*sa* ", &
            " 2>/dev/null | uniq | while read nom "
        endif
        write(600,*)" do "
        write(600,*)" nombis=${nom/'DATA'/'OUTPUT'}"
        write(600,*)" nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
        write(600,*)" echo '",param_best%Tzero,datatps(i)%dhypo,val/15.0_wr,"' > toto2.txt"
        write(600,*)" ./BIN/sac_bin2txt.exe $nom $nomter.txt < toto2.txt "
        write(600,*)" psxy $geozone $geoproj -W2,$thegray1 -O -K $nomter.txt -: >>  $file"
        write(600,*)" done "
        write(600,*)"rm -rf toto2.txt"
        ! -------                                              --------    . COMPOSANTE E
        if (datatps(i)%sta%staname(4:4)=='0') then                         ! nom station en trois caractere + "0"
          write(600,'(9a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*E.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*E.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*E.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*E.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
            " 2>/dev/null | uniq | while read nom "
        else
          write(600,'(5a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*E.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*E.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//".*Z.*sa* ", &
            " 2>/dev/null | uniq | while read nom "
        endif
        write(600,*)" do "
        write(600,*)" nombis=${nom/'DATA'/'OUTPUT'}"
        write(600,*)" nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
        write(600,*)" echo '",param_best%Tzero,datatps(i)%dhypo,val/15.0_wr,"' > toto2.txt"
        write(600,*)" ./BIN/sac_bin2txt.exe $nom $nomter.txt < toto2.txt "
        write(600,*)" psxy $geozone $geoproj -W2,$thegray2 -O -K $nomter.txt -: >>  $file"
        write(600,*)" done "
        write(600,*)"rm -rf toto2.txt"
        ! -------                                              --------    . COMPOSANTE N
        if (datatps(i)%sta%staname(4:4)=='0') then                         ! nom station en trois caractere + "0"
          write(600,'(9a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*N.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*N.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*N.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*N.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
            " 2>/dev/null | uniq | while read nom "
        else
          write(600,'(5a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*N.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*N.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*sa* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//".*Z.*sa* ", &
            " 2>/dev/null | uniq | while read nom "
        endif
        write(600,*)" do "
        write(600,*)" nombis=${nom/'DATA'/'OUTPUT'}"
        write(600,*)" nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
        write(600,*)" echo '",param_best%Tzero,datatps(i)%dhypo,val/15.0_wr,"' > toto2.txt"
        write(600,*)" ./BIN/sac_bin2txt.exe $nom $nomter.txt < toto2.txt "
        write(600,*)" psxy $geozone $geoproj -W2,$thegray3 -O -K $nomter.txt -: >>  $file"
        write(600,*)" done "
        write(600,*)"rm -rf toto2.txt"
      enddo
      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .
    endif
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .


    ! prendre en considération un modo incliné !  ... à revoir !


    ! -------                                                  --------    . traits réels (lois appliquées au jeu de paramètres)
    ! -------                                                  --------    . ONDES REFLECHIES puis REFRACTEES
    max = val*1.09_wr
    min = 0.0_wr
    ! -------                                                  --------    . pPn par modèle
    !call pPn_sSn(param_best,max,Tpsmax,Tppmax)
    !call pPn_sSn(param_best,min,Tpsmin,Tppmin)
    !write(600,*)"echo -e "" ",min,Tppmin,"\n",max,Tppmax," \"
    !write(600,*)" "" | psxy $geozone $geoproj -W0.01i,LIGHTGREEN -O -K >>  $file"
    ! -------                                                  --------    . sSn par modèle
    !write(600,*)"echo -e "" ",min,Tpsmin,"\n",max,Tpsmax," \"
    !write(600,*)" "" | psxy $geozone $geoproj -W0.01i,LIGHTGREEN -O -K >>  $file"
    ! -------                                                  --------    . ONDES REFLECHIES
    !do i=int(param_best%Zhypo+2.0_wr),int(val*1.09_wr+0.5_wr),2
      !min = real(i-1,wr)
      !max = real(i+1,wr)
      ! -------                                                --------    . PmP par modèle
      !call reflechie(param_best,min,Tpsmin,Tppmin)
      !call reflechie(param_best,max,Tpsmax,Tppmax)
      !write(600,*)"echo -e "" ",min,Tppmin,"\n",max,Tppmax," \"
      !write(600,*)" "" | psxy $geozone $geoproj -W0.01i,LIGHTORANGE -O -K >>  $file"
      ! -------                                                --------    . SmS par modèle
      !write(600,*)"echo -e "" ",min,Tpsmin,"\n",max,Tpsmax," \"
      !write(600,*)" "" | psxy $geozone $geoproj -W0.01i,LIGHTORANGE -O -K >>  $file"
    !enddo
    ! -------                                                  --------    . ONDES REFLECHIES 2
    !do i=int(param_best%Zhypo+2.0_wr),int(val*1.09_wr+0.5_wr),2
      !min = real(i-1,wr)
      !max = real(i+1,wr)
      ! -------                                                --------    . 2PmP par modèle
      !call reflechie2(param_best,min,Tpsmin,Tppmin)
      !call reflechie2(param_best,max,Tpsmax,Tppmax)
      !write(600,*)"echo -e "" ",min,Tppmin,"\n",max,Tppmax," \"
      !write(600,*)" "" | psxy $geozone $geoproj -W0.01i,LIGHTORANGE -O -K >>  $file"
      ! -------                                                --------    . 2SmS par modèle
      !write(600,*)"echo -e "" ",min,Tpsmin,"\n",max,Tpsmax," \"
      !write(600,*)" "" | psxy $geozone $geoproj -W0.01i,LIGHTORANGE -O -K >>  $file"
    !enddo
    ! -------                                                  --------    . ONDES DIRECTES
    min = 0.0_wr
    max = val*1.09_wr
    ! -------                                                  --------    . Pg par modèle
    call directe(param_best,min,Tpsmin,Tppmin)
    call directe(param_best,max,Tpsmax,Tppmax)
    write(600,*)"echo -e "" ",min,Tppmin,"\n",max,Tppmax," \"
    write(600,*)" "" | psxy $geozone $geoproj -W0.05i,$pp -O -K >>  $file"
    ! -------                                                  --------    . Sg par modèle
    write(600,*)"echo -e "" ",min,Tpsmin,"\n",max,Tpsmax," \"
    write(600,*)" "" | psxy $geozone $geoproj -W0.05i,$ss -O -K >>  $file"
    ! -------                                                  --------    . ONDES REFRACTEES
    min = discritiqueH
    ! -------                                                  --------    . Pn par modèle
    if(FLAG_non_tabulaire) then                                            ! moho incliné
      a_sta%lon=param_best%lon
      a_sta%lat=param_best%lat
      a_sta%alti=0.0_wr
      call refracte_mohovar(acentroid,param_best,a_sta,min,Tpsmin,Tppmin)
      call refracte_mohovar(acentroid,param_best,a_sta,max,Tpsmax,Tppmax)
    else
      call refracte(param_best,min,Tpsmin,Tppmin)
      call refracte(param_best,max,Tpsmax,Tppmax)
    endif
    write(600,*)"echo -e "" ",min,Tppmin,"\n",max,Tppmax," \"
    write(600,*)" "" | psxy $geozone $geoproj -W0.05i,$pp,- -O -K >>  $file"
    ! -------                                                  --------    . Sn par modèle
    write(600,*)"echo -e "" ",min,Tpsmin,"\n",max,Tpsmax," \"
    write(600,*)" "" | psxy $geozone $geoproj -W0.05i,$ss,- -O -K >>  $file"
    ! -------                                                  --------    . Pn par regression
    write(600,*)"echo -e "" ",min,a_Pn*min+b_Pn,"\n",max,a_Pn*max+b_Pn," \"
    write(600,*)" ""| psxy $geozone $geoproj -W0.001i -O -K >> $file"
    ! -------                                                  --------    . Sn par regression
    write(600,*)"echo -e "" ",min,a_Sn*min+b_Sn,"\n",max,a_Sn*max+b_Sn," \"
    write(600,*)" ""| psxy $geozone $geoproj -W0.001i -O -K >> $file"
    ! -------                                                  --------    . Pg par regression
    min = 0.0_wr
    write(600,*)"echo -e "" ",min,a_Pg*min,"\n",max,a_Pg*max," \"
    write(600,*)" ""| psxy $geozone $geoproj -W0.001i -O -K >>  $file"
    ! -------                                                  --------    . Sg par regression
    write(600,*)"echo -e "" ",min,a_Sg*min,"\n",max,a_Sg*max," \"
    write(600,*)" ""| psxy $geozone $geoproj -W0.001i -O -K >>  $file"
    ! -------                                                  --------    .
    do i=1,nPg                                                             ! points réels
      write(600,*)"echo",pt_Pg(i,1),pt_Pg(i,2),pt_Pg(i,3),pt_Pg(i,4)," \"
      write(600,*)" | psxy $geozone $geoproj -O -K -St0.1i -Wthinnest -Ey -COUTPUT/GMT/colorpal3.cpt >>  $file"
    enddo
    ! -------                                                  --------    .
    do i=1,nPn                                                             ! points réels
      write(600,*)"echo",pt_Pn(i,1),pt_Pn(i,2),pt_Pn(i,3),pt_Pn(i,4)," \"
      write(600,*)" | psxy $geozone $geoproj -O -K -Si0.1i -Wthinnest -Ey -COUTPUT/GMT/colorpal3.cpt >>  $file"
    enddo
    ! -------                                                  --------    .
    do i=1,nSg                                                             ! points réels
      write(600,*)"echo",pt_Sg(i,1),pt_Sg(i,2),pt_Sg(i,3),pt_Sg(i,4)," \"
      write(600,*)" | psxy $geozone $geoproj -O -K -Ss0.1i -Wthinnest -Ey -COUTPUT/GMT/colorpal3.cpt >>  $file"
    enddo
    ! -------                                                  --------    .
    do i=1,nSn                                                             ! points réels
      write(600,*)"echo",pt_Sn(i,1),pt_Sn(i,2),pt_Sn(i,3),pt_Sn(i,4)," \"
      write(600,*)" | psxy $geozone $geoproj -O -K -Sd0.1i -Wthinnest -Ey -COUTPUT/GMT/colorpal3.cpt >>  $file"
    enddo
    ! -------                                                  --------    . nom sta
    l=1
    do i=1,nbtps+1
      nomstadoublets(i)='xxxx'
    enddo
    do i=1,nbtps
      ! -------                                                --------    . nom déja affiché ?
      deja =.true.
      do k=1,l
        if (nomstadoublets(k).eq.datatps(i)%sta%staname) deja=.false.
      enddo
      ! -------                                                --------    . affiche
      if (deja) then
        if (mod(i,3)==0) then
          write(600,*)"echo '",datatps(i)%dhypo,max_0*0.9_wr, &
          " 7 90 1 5 "//datatps(i)%sta%staname//"' | pstext $geoproj $geozone -O -K -C2 >> $file"
          nomstadoublets(l)=datatps(i)%sta%staname
          l=l+1
        elseif (mod(i,2)==0) then
          write(600,*)"echo '",datatps(i)%dhypo,max_0*0.8_wr, &
          " 7 90 1 5 "//datatps(i)%sta%staname//"' | pstext $geoproj $geozone -O -K -C2 >> $file"
          nomstadoublets(l)=datatps(i)%sta%staname
          l=l+1
        else
          write(600,*)"echo '",datatps(i)%dhypo,max_0*0.7_wr, &
          " 7 90 1 5 "//datatps(i)%sta%staname//"' | pstext $geoproj $geozone -O -K -C2 >> $file"
          nomstadoublets(l)=datatps(i)%sta%staname
          l=l+1
        endif
      endif
    enddo
    ! -------                                                  --------    .
    write(600,*)"#########################################"                ! Légende
    ! -------                                                  --------    . figurés Pg
    X = val * 1.2_wr * 0.085_wr
    Y = max_0 * 1.1_wr  * 0.745_wr
    write(600,*)"echo",X,Y," \"
    write(600,*)" | psxy $geozone $geoproj -O -K -St0.1i -Wthinnest -Gyellow >>  $file"
    X = val * 1.2_wr  * 0.1_wr
    Y = max_0 * 1.1_wr  * 0.75_wr
    write(600,*)"echo """,X,Y," \"
    write(600,*)"15 0 4 LM ondes compressives directes"" | pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    . figurés Pn
    X = val * 1.2_wr * 0.085_wr
    Y = max_0 * 1.1_wr  * 0.695_wr
    write(600,*)"echo",X,Y," \"
    write(600,*)" | psxy $geozone $geoproj -O -K -Ss0.1i -Wthinnest -Gyellow >>  $file"
    X = val * 1.2_wr  * 0.1_wr
    Y = max_0 * 1.1_wr  * 0.70_wr
    write(600,*)"echo """,X,Y," \"
    write(600,*)"15 0 4 LM ondes cisaillantes directes"" | pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    . figurés Sg
    X = val * 1.2_wr * 0.085_wr
    Y = max_0 * 1.1_wr  * 0.645_wr
    write(600,*)"echo",X,Y," \"
    write(600,*)" | psxy $geozone $geoproj -O -K -Si0.1i -Wthinnest -Gyellow >>  $file"
    X = val * 1.2_wr  * 0.1_wr
    Y = max_0 * 1.1_wr  * 0.65_wr
    write(600,*)"echo """,X,Y," \"
    write(600,*)"15 0 4 LM ondes compressives r\351fract\351es"" | pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    . figurés Sn
    X = val * 1.2_wr * 0.085_wr
    Y = max_0 * 1.1_wr  * 0.595_wr
    write(600,*)"echo",X,Y," \"
    write(600,*)" | psxy $geozone $geoproj -O -K -Sd0.1i -Wthinnest -Gyellow >>  $file"
    X = val * 1.2_wr  * 0.1_wr
    Y = max_0 * 1.1_wr  * 0.60_wr
    write(600,*)"echo """,X,Y," \"
    write(600,*)"15 0 4 LM ondes cisaillantes r\351fract\351es"" | pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    . doites des modèles
    X = val * 1.2_wr  * 0.1_wr
    Y = max_0 * 1.1_wr  * 0.80_wr
    write(600,*)"echo """,X,Y," \"
    write(600,*)"15 0 4 LM selon le mod\350le de terre"" | pstext $geozone $geoproj -O -K >> $file"
    X = val * 1.2_wr * 0.085_wr
    write(600,*)"echo -e "" ",X-0.075_wr*X,Y,"\n",X-0.01_wr*X,Y," \"
    write(600,*)" "" | psxy $geozone $geoproj -W0.05i,$ss -O -K >>  $file"
    write(600,*)"echo -e "" ",X+0.01_wr*X,Y,"\n",X+0.075_wr*X,Y," \"
    write(600,*)" "" | psxy $geozone $geoproj -W0.05i,$pp -O -K >>  $file"
    ! -------                                                  --------    . doites de régressions
    X = val * 1.2_wr  * 0.1_wr
    Y = max_0 * 1.1_wr  * 0.85_wr
    write(600,*)"echo """,X,Y," \"
    write(600,*)"15 0 4 LM r\351gressions lin\351aires"" | pstext $geozone $geoproj -O -K >> $file"
    X = val * 1.2_wr * 0.085_wr
    write(600,*)"echo -e "" ",X-0.05_wr*X,Y,"\n",X+0.05_wr*X,Y," \"
    write(600,*)" ""| psxy $geozone $geoproj -W0.001i -O -K >> $file"
    ! -------                                                  --------    .
    write(600,*)"psscale -D1/-1/0.50E+01/0.25ch -B.25:""pond\351ration"": -S -I -COUTPUT/GMT/colorpal3.cpt -O -K >> $file"
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"#########################################"
    write(600,*)"ps2raster OUTPUT/GMT/hodochrone"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/hodochrone"//"-"//trim(adjustl(numberfile))//".pdf ", &
    "OUTPUT/figures/hodochrone"//"-"//trim(adjustl(numberfile))//".pdf "
    write(600,*)"#########################################"
    write(600,*)"ELAPSED=$(($SECONDS-$BEFORE))"
    write(600,*)" echo $ELAPSED secondes"
    call system_clock(Nnewtime,ratetime)
    tl=(real(Nnewtime,wr)-real(Noldtime,wr))/real(ratetime,wr)
    write(*,'(a9,i2.2,'':'',i2.2,'':'',f9.2)')' temps : ',int(tl/3600.0_wr,wi), &
      int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi),(tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr)
    ! -----------------------------------------------------------------    .
    deallocate (pt_Pg,pt_Pn,pt_Sg,pt_Sn)
    ! -----------------------------------------------------------------    .
  end subroutine GMT_Hodochrone

END MODULE figure_GMThodo



! *********************************************************************    .
! *********************************************************************    .


