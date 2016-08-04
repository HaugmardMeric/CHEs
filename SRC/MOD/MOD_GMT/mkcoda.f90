! permet la création des scripts GMT pour l'hodochrone et le calcul de la magnitude Md
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTcoda

    use modparam

    implicit none


    private

    public  :: GMT_coda

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

  subroutine GMT_coda(j,nbtps,datatps,param_best,xmaxcercle,lon,lat,acentroid)
    ! -------                                                  --------    .mh
    ! Calcul les regressions sur les hodochrones
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
    real(KIND=wr), intent (in) :: lon,lat
    type(amoho_centroid), intent (in) :: acentroid
    ! -------                                                  --------    .
    type(date_sec) :: one_tps
    type(stations) :: a_sta
    integer (KIND=wi) :: i,k,l,ok
    integer (KIND=wi) :: nPg, nPn ,nSg, nSn, Noldtime, Nnewtime, ratetime
    real (KIND=wr), dimension(:,:), allocatable :: pt_Pg, pt_Pn ,pt_Sg, pt_Sn
    real (KIND=wr) :: a_Pg,R2_Pg,a_Pn,b_Pn,R2_Pn,a_Sg,R2_Sg,a_Sn,b_Sn,R2_Sn
    real (KIND=wr) :: min, max_0, max, discritiqueH
    real (KIND=wr) :: tl, val,dx,dy,coefa1,coefb1,coefa2,coefb2
    real (KIND=wr) :: duree,duree1,duree2,ml,depi,depi2
    real (KIND=wr) :: Tpsmin,Tppmin,Tpsmax,Tppmax
    real(KIND=wr) :: lon1,lon2,lat1,lat2,v1,v2
    character (LEN=5) :: numberfile
    character (LEN=4) :: nomstadoublets(nbtps+1),kstnm
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
        write(*,*)'problème dans GMT_coda : onde P ni N ni G'
        stop
      endif
      if(datatps(i)%andS.eq.'S') then
        if(datatps(i)%typeonde.eq.'N') then
          nSn = nSn + 1
        elseif(datatps(i)%typeonde.eq.'G') then
          nSg = nSg + 1
        else
          write(*,*)'problème dans GMT_coda : onde S ni N ni G'
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
        write(*,*)'problème dans GMT_coda : onde P ni N ni G'
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
          write(*,*)'problème dans GMT_coda : onde S ni N ni G'
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
    write(*,*)"ecriture du script GMT Coda"
    write(600,*)"thegray1=240/240/240"
    write(600,*)"thegray2=230/230/230"
    write(600,*)"thegray3=220/220/220"
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#########################################"
    write(600,*)"################## coda #################"
    write(600,*)"echo 'execution du script GMT coda Md'"
    write(numberfile(1:5),'(i5)')j
    write(600,*)"file=OUTPUT/GMT/coda"//"-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"geoproj=-JX13i/8i"                                        ! système de projection
    if (xmaxcercle/2.0_wr.lt.datatps(nbtps)%dhypo) then
      val=xmaxcercle/2.0_wr
    else
      val=datatps(nbtps)%dhypo
    endif
    write(600,'(a,E13.7,a,E13.7)')"geozone=-R0.0/",val*1.1_wr,"/-30.0/",max_0+500.0_wr
    if (max_0.gt.60.0_wr) then
      write(600,*)"psbasemap $geozone $geoproj -Ba50f25:""distance hypocentrale (km)"":",&
      "/a100f25:""temps d'arriv\351es des ondes (s)"":WenS -Xc -Yc -K >  $file"
    else
      write(600,*)"psbasemap $geozone $geoproj -Ba20f5:""distance hypocentrale (km)"":",&
      "/a100f25:""temps d'arriv\351es des ondes (s)"":WenS -Xc -Yc -K >  $file"
    endif
    ! -----------------------------------------------------------------    .
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
    if (IsNaN(discritiqueH)) min=max-5.0_wr ! arbitraire ...
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
    ! -----------------------------------------------------------------    .
    inquire (_FILE_DIR_="DATA/sac-"//trim(adjustl(numberfile)),exist=existe1) ! option différente selon compilo !
    if ((existe1).and.(tracessac)) then
      ! ---------------------------------------------------------------    . plot traces si existes
      do i=1,nbtps
        ! lecture un peu archaïque, mais permet un peu de souplesse dans le non des station et des fichiers sac
        ! -------                                              --------    . COMPOSANTE Z
        if (datatps(i)%sta%staname(4:4)=='0') then                         ! nom station en trois caractere + "0"
          write(600,'(5a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*Z.*SA* ", &
            " 2>/dev/null | uniq | while read nom "
        else
          write(600,'(3a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//".*Z.*SA* ", &
            " 2>/dev/null | uniq | while read nom "
        endif
        write(600,*)" do "
        write(600,*)" nombis=${nom/'DATA'/'OUTPUT'}"
        write(600,*)" nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
        write(600,*)" psxy $geozone $geoproj -W2,$thegray2 -O -K $nomter.txt -: >>  $file"
        write(600,*)" done "
        write(600,*)"rm -rf toto2.txt"
        ! -------                                              --------    . COMPOSANTE E
        if (datatps(i)%sta%staname(4:4)=='0') then                         ! nom station en trois caractere + "0"
        write(600,'(5a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*E.*SA* ", &
          "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*E.*SA* ", &
          "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*E.*SA* ", &
          "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*E.*SA* ", &
            " 2>/dev/null | uniq | while read nom "
        else
          write(600,'(3a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*E.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//".*Z.*SA* ", &
            " 2>/dev/null | uniq | while read nom "
        endif
        write(600,*)" do "
        write(600,*)" nombis=${nom/'DATA'/'OUTPUT'}"
        write(600,*)" nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
        write(600,*)" psxy $geozone $geoproj -W2,$thegray2 -O -K $nomter.txt -: >>  $file"
        write(600,*)" done "
        write(600,*)"rm -rf toto2.txt"
        ! -------                                              --------    . COMPOSANTE N
        if (datatps(i)%sta%staname(4:4)=='0') then                         ! nom station en trois caractere + "0"
        write(600,'(5a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*N.*SA* ", &
          "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*N.*SA* ", &
          "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*N.*SA* ", &
          "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*N.*SA* ", &
            " 2>/dev/null | uniq | while read nom "
        else
          write(600,'(3a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*SA* ", &
            "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//".*N.*SA* ", &
            " 2>/dev/null | uniq | while read nom "
        endif
        write(600,*)" do "
        write(600,*)" nombis=${nom/'DATA'/'OUTPUT'}"
        write(600,*)" nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
        write(600,*)" psxy $geozone $geoproj -W2,$thegray3 -O -K $nomter.txt -: >>  $file"
        write(600,*)" done "
        write(600,*)"rm -rf toto2.txt"
        ! -------                                              --------    .
      enddo
      ! ---------------------------------------------------------------    .
      ! calcul de la Magnitude Md                                          !
      ! ---------------------------------------------------------------    .
      ok=0
      open(599-j, FILE = 'OUTPUT/GMT/scriptmag'//trim(adjustl(numberfile))//'.sh',status='replace',iostat = ok)
      if(ok.ne.0) then
        write(*,*)'problème dans GMT_coda : OUTPUT/GMT/scriptmag'//trim(adjustl(numberfile))//'.sh n''existe pas '
        stop
      endif
      ! -------                                                --------    . distance Pn > Pg
      call directe(param_best,min,Tpsmin,Tppmin)                           ! coef dir. droite Pg
      call directe(param_best,max,Tpsmax,Tppmax)
      coefa1=(Tppmax-Tppmin)/(max-min)
      coefb1=-coefa1*max+Tppmax
      if(FLAG_non_tabulaire) then                                          ! moho incliné
        a_sta%lon=param_best%lon
        a_sta%lat=param_best%lat
        a_sta%alti=0.0_wr
        call refracte_mohovar(acentroid,param_best,a_sta,min,Tpsmin,Tppmin)
        call refracte_mohovar(acentroid,param_best,a_sta,max,Tpsmax,Tppmax)
      else
        call refracte(param_best,min,Tpsmin,Tppmin)                        ! coef dir. droite Pg
        call refracte(param_best,max,Tpsmax,Tppmax)
      endif
      coefa2=(Tppmax-Tppmin)/(max-min)
      coefb2=-coefa2*max+Tppmax
      call deuxdroites(coefa1,coefb1,coefa2,coefb2,dx,dy)                  ! intersection des droites Pg et Pn et dx
      if (IsNaN(dx)) then
        write(*,*)'problème dans GMT_coda : dx = NaN'
        stop
      endif
      if (dx.lt.0.0_wr) then
        write(*,*)'problème dans GMT_coda : dx < 0, bizarre !'
        stop
      endif
      ! -------                                                --------    .
      do i=1,nbtps
        if ((((datatps(i)%typeonde=='G').and.(datatps(i)%dhypo.le.dx)).or. &
            ((datatps(i)%typeonde=='N').and.(datatps(i)%dhypo.ge.dx))).and. &
            (datatps(i)%depi.gt.10.0_wr)) then ! premieres arrivée P (Pg ou Pn), distance épi > 10 km sinon saturation
          ! lecture un peu archaïque, mais permet un peu de souplesse dans le non des station et des fichiers sac
          ! -------                                            --------    .
          if (datatps(i)%sta%staname(4:4)=='0') then                       ! nom station en trois caractères + "0"
            write(599-j,'(9a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:3)//".*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname(1:4)//".*Z.*sa* ", &
              " 2>/dev/null | uniq | while read nom "
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
            write(599-j,'(9a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*sa* ", &
              " 2>/dev/null | uniq | while read nom "
            write(600,'(9a)')"ls DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*SA* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/"//datatps(i)%sta%staname//"*Z.*sa* ", &
              "DATA/sac-"//trim(adjustl(numberfile))//"/*"//datatps(i)%sta%staname//"*Z.*sa* ", &
              " 2>/dev/null | uniq | while read nom "
          endif
          ! -------                                            --------    .
          write(599-j,'(a)')"do"
          write(600,'(a)')"do"
          write(599-j,'(a)')"deux=${nom/HZ/HE}"
          write(599-j,'(a)')"trois=${nom/HZ/HN}"
          ! -------                                            --------    . script sac
          write(599-j,'(a)')"sac << EOF >/dev/null"
          write(599-j,'(a)')"r $nom"
          write(599-j,'(a)')"mulf $nom"
          write(599-j,'(a)')"w 1.sac"
          write(599-j,'(a)')"r $deux"
          write(599-j,'(a)')"mulf $deux"
          write(599-j,'(a)')"w 2.sac"
          write(599-j,'(a)')"r $trois"
          write(599-j,'(a)')"mulf $trois"
          write(599-j,'(a)')"w 3.sac"
          write(599-j,'(a)')"r 1.sac"
          write(599-j,'(a)')"addf 2.sac"
          write(599-j,'(a)')"addf 3.sac"
          ! -------                                            --------    .
          ! new.sac = env(Z)**2 + env(N)**2 + env(E)**2 !
          ! -------                                            --------    .
          write(599-j,'(a)')"w new.SAC"
          write(599-j,'(a)')"quit"
          write(599-j,'(a)')"EOF"
          ! -------                                            --------    .
          write(599-j,'(a)')"nombis=${nom/'DATA'/'OUTPUT'}"
          write(599-j,'(a)')"nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
          write(600,'(a)')"nombis=${nom/'DATA'/'OUTPUT'}"
          write(600,'(a)')"nomter=${nombis/sac-"//trim(adjustl(numberfile))//"/GMT}"
          write(599-j,*)" echo '",j,datatps(i)%tpsTh,param_best%Tzero,datatps(i)%depi, &
            datatps(i)%dhypo,val/15.0_wr,"'> toto1.txt"
          write(599-j,'(a)')" ./BIN/sac_coda.exe new.SAC $nomter-coda-"//trim(adjustl(numberfile))//".txt < toto1.txt "
          write(599-j,'(a)')"rm -rf toto1.txt "
          write(599-j,'(a)')"done "
          ! -------                                            --------    .
          write(600,*)"if test -f $nomter-coda-"//trim(adjustl(numberfile))//".txt ; then"
          write(600,*)" psxy $geozone $geoproj -W4,red -O -K $nomter-coda-"//trim(adjustl(numberfile))//".txt >>  $file"
          write(600,*)"fi"
          write(600,'(a)')"done "
          ! -------                                            --------    .
          write(599-j,'(a)')"rm -rf 1.sac 2.sac 3.sac new.SAC"
          ! -------                                            --------    .
        else
          write(599-j,'(a)')'### no mag : '//datatps(i)%sta%staname//" - "//trim(adjustl(numberfile))
        endif
      enddo
      Close(599-j)
      ! -------                                                --------    . pour gfortran :
      !call execute_command_line ("chmod +x OUTPUT/GMT/scriptmag"//trim(adjustl(numberfile))//".sh", wait=.true.)
      !call execute_command_line ("./OUTPUT/GMT/scriptmag"//trim(adjustl(numberfile))//".sh", wait=.true.)
      ! -------                                                --------    . pour gfortran & ifort :
      call system ("chmod +x OUTPUT/GMT/scriptmag"//trim(adjustl(numberfile))//".sh")
      call system ("./OUTPUT/GMT/scriptmag"//trim(adjustl(numberfile))//".sh")
      ! -------                                                --------    .
      ! ---------------------------------------------------------------    .
    endif
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
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


    ! -------                                                  --------    . trace les magnitudes théoriques :
    do i=1,5 ! Ml
      ml=real(i,wr)
      write(numberfile(1:5),'(i5)')i
      depi=0.0_wr
      do while (depi.lt.(val*1.1_wr))
        duree1=10.0_wr**((ml+ 0.87_wr-0.0035_wr*depi)/2.0_wr)+depi/param_best%VC
        depi2=depi+2._wr
        duree2=10.0_wr**((ml+ 0.87_wr-0.0035_wr*depi2)/2.0_wr)+depi2/param_best%VC
        write(600,*)"echo -e "" ",depi,duree1,"\n",depi2,duree2," \"
        write(600,*)" "" | psxy $geozone $geoproj -W0.01i,blue -O -K >> $file"
        depi=depi2
      enddo
      depi=10.0_wr
      duree1=10.0_wr**((ml+ 0.87_wr-0.0035_wr*depi)/2.0_wr)!depi/param_best%VC
      write(600,*)"echo '",depi,duree1," 15 0 1 5 Md ="//trim(adjustl(numberfile))//"'",&
      " | pstext $geoproj $geozone -O -K -C2 >> $file"
    enddo
    write(numberfile(1:5),'(i5)')j
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone $geoproj -Ba0 -O >>  $file"
    write(600,*)"#########################################"
    write(600,*)"ps2raster OUTPUT/GMT/coda"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/coda"//"-"//trim(adjustl(numberfile))//".pdf ", &
    "OUTPUT/figures/coda"//"-"//trim(adjustl(numberfile))//".pdf "
    ! -----------------------------------------------------------------    .
    deallocate (pt_Pg,pt_Pn,pt_Sg,pt_Sn)
    ! -----------------------------------------------------------------    .
    ! plot carte avec des cercles à chaque station, focntion de la magnitude
    ! -----------------------------------------------------------------    .
    write(600,*)"#########################################"
    write(*,*)"ecriture du script GMT_coda_map "
    write(600,*)"#########################################"
    write(numberfile(1:5),'(i5)')j
    write(600,*)"file=OUTPUT/GMT/coda_map-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"gmtset BASEMAP_TYPE plain"
    write(600,*)"labasemap1=-Bpa2g1.f.5/a1g1.f.25WeSn"
    ! -------                                                  --------    .
    v1 = 2.0_wr * pi * rT / 360.0_wr                                       ! km / degree en lon
    v2 = 2.0_wr * pi * rT * sin((90.0_wr-lat)/180.0_wr*pi) /360.0_wr       ! km / degree en lat
    ! -------                                                  --------    .
    lon1 = lon - (xmaxcercle / v2 * 1.125_wr) / 2.0_wr
    lon2 = lon + (xmaxcercle / v2 * 1.125_wr) / 2.0_wr
    lat1 = lat - (xmaxcercle / v1 * 1.125_wr) / 2.0_wr
    lat2 = lat + (xmaxcercle / v1 * 1.125_wr) / 2.0_wr
    ! -------                                                  --------    .
    write(600,'(a10,E13.7,a1,E13.7,a1,E13.7,a1,E13.7)')"geozone=-R",lon1,"/",lon2,"/",lat1,"/",lat2
    write(600,'(a11,E13.7,a1,E13.7,a3)')"geoproj=-JC",lon,"/",lat,"/7i"
    ! -------                                                  --------    .
    write(600,*)"bluef=""0/0/100"" "
    write(600,*)"makecpt -Cseis -I -T1.5/4.5/0.01 -Z > OUTPUT/GMT/neis.cpt"
    write(600,*)"pscoast $geozone $geoproj -Df+ -Ia/$bluef -S240/255/255 -G180/238/180 -W1 -K -Xc -X5.5i -Yc $labasemap1 >  $file"
    write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -m >> $file"
    ! -------                                                  --------    .
    inquire (FILE="OUTPUT/files/mag-"//trim(adjustl(numberfile))//".d",exist=existe1) ! option différente selon compilo !
    ok=0
    ! -------                                                  --------    .
    do while (.not.existe1)
      ok=ok+1
      call sleep(5)
      inquire (FILE="OUTPUT/files/mag-"//trim(adjustl(numberfile))//".d",exist=existe1) ! option différente selon compilo !
      if (ok.gt.5) then
        ! write(*,*)'problème dans GMT_coda : le fichier OUTPUT/files/mag-'//trim(adjustl(numberfile))//'.d n''existe pas, ok=',ok
        existe1=.true.
      endif
    enddo
    ! -------                                                  --------    .
    ok=0
    ! -------                                                  --------    .
    open(111, FILE ="OUTPUT/files/mag-"//trim(adjustl(numberfile))//".d",status='old',iostat = ok)
    if (ok .ne. 0) then
      write(600,*)'psxy $geozone $geoproj OUTPUT/GMT/ellipse-'//trim(adjustl(numberfile))//'.txt ', &
        '-Sa0.5 -W1,gray -Gblue -O >> $file'
    else
      do while(ok .eq. 0)
        ! -------                                              --------    .
        read(111,*,iostat = ok)kstnm,ml,duree,depi
        if (kstnm(4:4)==' ')kstnm(4:4)='0'
        if (ok .eq. 0) then
          do i=1,nbtps
            if (datatps(i)%sta%staname==kstnm) then ! si bonne station
               write(600,*)"echo ",datatps(i)%sta%lon,datatps(i)%sta%lat,ml,ml*0.08_wr-0.04_wr, &
                 " | psxy $geozone $geoproj -Sci -Wthinnest -O -K -Ba0 -COUTPUT/GMT/neis.cpt >>  $file"
            endif
          enddo
        endif
        ! -------                                              --------    .
      end do
      close(111)
      ! -----------------------------------------------------------------    .
      write(600,*)"echo '",lon,lat,"0 300 300' | psxy $geozone $geoproj -SE -W5,-- -O -K -N >> $file"
      ! -------                                                --------    .
      write(600,*)'psxy $geozone $geoproj OUTPUT/GMT/ellipse-'//trim(adjustl(numberfile))//'.txt ', &
        '-Sa0.5 -W1,gray -Gblue -O -K >> $file'
      write(600,*)"echo '>' > OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'N4' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'G0.3i' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.08i  0/0/205 0.5p 0.525i  M@-d@- : 1.5' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.12i  0/160/183 0.5p 0.525i  M@-d@- : 2.0' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.16i  090/255/030 0.5p 0.525i  M@-d@- : 2.5' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.20i  255/255/0 0.5p 0.525i  M@-d@- : 3.0' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'G0.15i' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.24i  255/170/0 0.5p 0.525i  M@-d@- : 3.5' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.28i  255/042/0 0.5p 0.525i  M@-d@- : 4.0' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.32i  173/0/0 0.5p 0.525i  M@-d@- : 4.5' >> OUTPUT/GMT/mag.legend"
      write(600,*)"echo 'S 0.28i c 0.36i  0/0/0 0.5p .525i M@-d@- : 5.0' >> OUTPUT/GMT/mag.legend"
      write(600,*)"pslegend -Dx4.5i/-0.4i/7i/1.i/TC $geozone $geoproj -O OUTPUT/GMT/mag.legend >> $file"
    endif
    ! -------                                                  --------    .
    write(600,*)"ps2raster OUTPUT/GMT/coda_map-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/coda_map-"//trim(adjustl(numberfile))//".pdf ", &
      "OUTPUT/figures/coda_map-"//trim(adjustl(numberfile))//".pdf"
    ! -----------------------------------------------------------------    .
    write(600,*)"#########################################"
    write(600,*)"ELAPSED=$(($SECONDS-$BEFORE))"
    write(600,*)" echo $ELAPSED secondes"
    call system_clock(Nnewtime,ratetime)
    tl=(real(Nnewtime,wr)-real(Noldtime,wr))/real(ratetime,wr)
    write(*,'(a9,i2.2,'':'',i2.2,'':'',f9.2)')' temps : ',int(tl/3600.0_wr,wi), &
    int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi),(tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr)
    ! -----------------------------------------------------------------    .
  end subroutine GMT_coda

    ! -----------------------------------------------------------------    .

  subroutine deuxdroites(a1,b1,a2,b2,x,y)
    ! -------                                                  --------    .mh
    ! point d'intesection (x,y) de deux droites 1 et 2, 
    ! de coef. dir. a et ordonnée à l'origine b
    ! -------                                                  --------    .
    implicit none
    real(KIND=wr), intent (in) :: a1,b1,a2,b2
    real(KIND=wr), intent (out) :: x,y
    ! -----------------------------------------------------------------    .
    if ((a1-a2).ne.0.0_wr) then
      x=(b2-b1)/(a1-a2)
      y=a1*x+b1
    else
      x=0.0_wr
      y=a1*x+b1
    endif
    ! -----------------------------------------------------------------    .
  end subroutine deuxdroites

END MODULE figure_GMTcoda



! *********************************************************************    .
! *********************************************************************    .


