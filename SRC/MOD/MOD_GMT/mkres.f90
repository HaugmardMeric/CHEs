! permet la création des scripts GMT pour une carte des résidus par stations
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMTres

    use modparam

    implicit none

    private

    public  :: GMT_res, GMT_resSTA


    ! -----------------------------------------------------------------    .

    TYPE sta_tps
        character(LEN=4) :: staname                                        ! nom de la station
        real(KIND=wr) :: lonSTA, latSTA                                    ! coordonnées de la station et résidus
        real(KIND=wr) :: Ttot, TPn, TPg, TSn, TSg                          ! résidus (somme absolue, puis par type d'onde)
        real(KIND=wr) :: TpsPn, TpsPg, TpsSn, TpsSg                        ! temps (absolue)
        real(KIND=wr) :: pdsPn, pdsPg, pdsSn, pdsSg                        ! pondérations par type d'onde
        real(KIND=wr) :: disthypo
    END TYPE sta_tps

CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMT_resSTA(nbsta,nbtps,D,nomsta)
    ! -------                                                  --------    .mh
    ! plot residus des tous les séismes à chaque station pour un modèle (ensemble des 1000 meilleurs)
    ! -------                                                  --------    .
    use typetemps
    use time
    use sub_param
    use statistiques
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent(in) :: nbtps(nbseismes),nbsta                 ! nombre de données de temps et de station possible
    type(dataall), intent(in) :: D(nbseismes)                              ! données
    character(LEN=4), dimension(:), allocatable, intent(out)  :: nomsta
    ! -------                                                  --------    .
    ! pour "biner" les résidus et definir si la distribution est gaussienne
    ! delta = 0.05 sec ;  borne inf = numero * delta - delta/2
    real(KIND=wr) :: binage(-5000:5000), binageR(-5000:5000)
    real(KIND=wr) :: a,b,R2,XY(10001,3),Rp
    real(KIND=wr), dimension(:), allocatable :: XYbis(:,:)
    ! -------                                                  --------    .
    real(KIND=wr) :: val,X,Y,ymay
    real(KIND=wr) :: lon,lat,alti
    real(KIND=wr) :: moy(4),ec(4),med(4),gauss(4)
    real(KIND=wr), dimension(:), allocatable  :: vec
    integer(KIND=wi) :: i,j,k,ok,sta,noctet
    integer(KIND=wi) :: iPg,iPn,iSg,iSn, nb(4)
    logical :: test, existe(4)
    character(LEN=4) :: aname(nbsta+1)
    integer(KIND=wi), dimension(:,:), allocatable :: anbname
    ! -----------------------------------------------------------------    . nb de stations différentes
    ! -------                                                  --------    . initialistaion
    aname(1)=D(1)%datatps(1)%sta%staname
    do i=2,nbsta+1
      aname(i)='123_'
    enddo
    sta=0
    ! -------                                                  --------    .
    do i=1,nbseismes
      do j=1,nbtps(i)
        test=.true.
        k=1
        do while (aname(k).ne.'123_')
          if (aname(k)==D(i)%datatps(j)%sta%staname) test=.false.
          k=k+1
        enddo
        if(test) then
          sta=sta+1
          aname(sta)=D(i)%datatps(j)%sta%staname
        endif
      enddo
    enddo
    allocate(nomsta(sta))
    ! -----------------------------------------------------------------    . un fichier de résidus par station
    do k=1,sta
      open(950, FILE ='OUTPUT/GMT/res-PG-'//aname(k)//'.d',status='replace')
      open(951, FILE ='OUTPUT/GMT/res-PN-'//aname(k)//'.d',status='replace')
      open(952, FILE ='OUTPUT/GMT/res-SG-'//aname(k)//'.d',status='replace')
      open(953, FILE ='OUTPUT/GMT/res-SN-'//aname(k)//'.d',status='replace')
      nomsta(k)=aname(k)
      iPg=0
      iPn=0
      iSg=0
      iSn=0
      val=0.0_wr
      do i=1,nbseismes
        do j=1,nbtps(i)
          if (aname(k)==D(i)%datatps(j)%sta%staname) then
            if (D(i)%datatps(j)%typeonde=='G') then
              write(950,*)D(i)%datatps(j)%dTP
              iPg=iPg+1
              if(D(i)%datatps(j)%andS=='S') then
                write(952,*)D(i)%datatps(j)%dTS
                iSg=iSg+1
              endif
            else if (D(i)%datatps(j)%typeonde=='N') then
              write(951,*)D(i)%datatps(j)%dTP
              iPn=iPn+1
              if(D(i)%datatps(j)%andS=='S') then
                write(953,*)D(i)%datatps(j)%dTS
                iSn=iSn+1
              endif
            else
              write(*,*)'problème dans GMT_resSTA : onde ni directe ni réfractée'
            endif
            if (abs(D(i)%datatps(j)%dTP).gt.val) val=abs(D(i)%datatps(j)%dTP)
            if (abs(D(i)%datatps(j)%dTS).gt.val) val=abs(D(i)%datatps(j)%dTS)
          endif
        enddo
      enddo
      val= max((real(int(val+0.55_wr,wi),wr)*10.0_wr)/10._wr,1.1_wr)       ! bornes du graph
      close(950)
      close(951)
      close(952)
      close(953)
      X = -val + 0.3_wr*(2.0_wr*val)
      Y = 95.0_wr
      ! ---------------------------------------------------------------    . script GMT
      write(*,*)"ecriture des script GMT_res STA "//aname(k)
      write(600,*)"#***************************************#"
      write(600,*)"#***************************************#"
      write(600,*)
      write(600,*)"echo 'execution du script GMT res STA'"//aname(k)
      write(600,*)"#########################################"
      write(600,*)"###########      histo         ##########"
      write(600,*)"#########################################"
      ! ---------------------------------------------------------------    .
      write(600,'(a10,E13.7,a1,E13.7,a6)')"geozone=-R",-val,"/",val,"/0/100"
      write(600,*)"geoproj=-JX4.25i"
      write(600,*)"file=OUTPUT/GMT/resSTA"//"-"//aname(k)//".ps"
      ! -------                                                --------    .
      if (iPg.gt.0) then
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:", &
        "'r\351sidus ondes compressives directes':/a10:'effectif (%)':nSeW ", &
        " -X2.5i -K > $file "
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Z1 OUTPUT/GMT/res-PG-"//aname(k)//".d -O -K >> $file"
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iPg,""" | pstext $geozone $geoproj -O -K >> $file"
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$pp -Z1 OUTPUT/GMT/res-PG-"//aname(k)//".d -L0/0  -O -K >> $file"
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iPg,""" | pstext $geozone $geoproj -O -K >> $file"
        write(600,*)"psbasemap $geozone $geoproj -Bg1/a0 -O -K >> $file "

      else
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:", &
        "'r\351sidus ondes compressives directes':/a10:'effectif (%)':nSeW ", &
        " -X2.5i -K > $file "
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : 0"" | pstext $geozone $geoproj -O -K >> $file"
      endif
      ! -------                                                --------    .
      if (iPn.gt.0) then
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:'r\351sidus ondes compressives r\351fract\351es'", &
        ":/a10:'effectif (%)':nSew -O -K -X5.5i >>  $file"
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Ba0 ", &
        " -Z1 -K OUTPUT/GMT/res-PN-"//aname(k)//".d -O -K >> $file "
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$pp -L0/0 -Ba0 ", &
        " -Z1 -K OUTPUT/GMT/res-PN-"//aname(k)//".d -O -K >> $file "
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iPn,""" | pstext $geozone $geoproj -O -K >> $file"
        write(600,*)"psbasemap $geozone $geoproj -Bg1/a0 -O -K >> $file "
      else
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:'r\351sidus ondes compressives r\351fract\351es'", &
        ":/a10:'effectif (%)':nSew -O -K -X5.5i >>  $file"
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : 0"" | pstext $geozone $geoproj -O -K >> $file"
      endif
      ! -------                                                --------    .
      if (iSg.gt.0) then
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:", &
        "'r\351sidus ondes cisaillantes directes':/a10:'effectif (%)':nSew -O  -K -Y5.5i >>  $file"
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Ba0", &
        " -Z1 -K OUTPUT/GMT/res-SG-"//aname(k)//".d -O -K >> $file "
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$ss -L0/0 -Ba0", &
        " -Z1 -K OUTPUT/GMT/res-SG-"//aname(k)//".d -O -K >> $file "
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iSg,""" | pstext $geozone $geoproj -O -K >> $file"
        write(600,*)"psbasemap $geozone $geoproj -Bg1/a0 -O -K >> $file "
      else
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:", &
        "'r\351sidus ondes cisaillantes directes':/a10:'effectif (%)':nSew -O  -K -Y5.5i >>  $file"
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : 0"" | pstext $geozone $geoproj -O -K >> $file"
      endif
      ! -------                                                --------    .
      if (iSn.gt.0) then
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:'r\351sidus ondes cisaillantes r\351fract\351es'", &
        ":/a10:'effectif (%)':nSeW -O -K -X-5.5i >>  $file"
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Ba0 ", &
        " -Z1 OUTPUT/GMT/res-SN-"//aname(k)//".d -O -K >> $file "
        write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$ss -L0/0 -Ba0 ", &
        " -Z1 OUTPUT/GMT/res-SN-"//aname(k)//".d -O -K >> $file "
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iSn,""" | pstext $geozone $geoproj -O -K >> $file"
        write(600,*)"psbasemap $geozone $geoproj -Bg1/a0 -O >> $file "
      else
        write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:'r\351sidus ondes cisaillantes r\351fract\351es'", &
        ":/a10:'effectif (%)':nSeW -O -K  -X-5.5i >>  $file"
        write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : 0"" | pstext $geozone $geoproj -O >> $file"
      endif
      ! -------                                                --------    .
      write(600,*)"ps2raster OUTPUT/GMT/resSTA-"//aname(k)//".ps -Tf -A"
      write(600,'(a)')"mv OUTPUT/GMT/resSTA-"//aname(k)//".pdf OUTPUT/figures/resSTA-"//aname(k)//".pdf"
      write(600,*)
    ! -----------------------------------------------------------------    .
    enddo
    ! -----------------------------------------------------------------    .

    write(600,*)"#########################################"
    write(600,*)"###########    all histo       ##########"
    write(600,*)"#########################################"

    write(600,*)"geoproj=-JX5i"
    write(600,*)"file=OUTPUT/GMT/Allres.ps"
    write(600,*)"cat OUTPUT/GMT/residus*pg OUTPUT/GMT/residus*pn OUTPUT/GMT/residus*sg ", &
      "OUTPUT/GMT/residus*sn > OUTPUT/GMT/allres.sn"
    write(600,*)"cat OUTPUT/GMT/residus*pg OUTPUT/GMT/residus*pn OUTPUT/GMT/residus*sg > OUTPUT/GMT/allres.sg"
    write(600,*)"cat OUTPUT/GMT/residus*pg OUTPUT/GMT/residus*pn > OUTPUT/GMT/allres.pn"
    write(600,*)"cat OUTPUT/GMT/residus*pg > OUTPUT/GMT/allres.pg"

    write(600,*)"pshistogram -W0.02  -Gred  OUTPUT/GMT/allres.sn -IO > OUTPUT/GMT/allres.all 2>/dev/null "
    write(600,*)"minmax -I1.5 OUTPUT/GMT/allres.all > OUTPUT/GMT/geozonehistoall 2>/dev/null "
    write(600,*)"read geozone < OUTPUT/GMT/geozonehistoall"

    write(600,*)"pshistogram  $geozone $geoproj -W0.02  -Gred  OUTPUT/GMT/allres.sn -K ", &
      "-Ba.5f.25g10:'r\351sidus (s)':/a5f1:'effectif':nSeW > $file"

    write(600,*)"pshistogram  $geozone $geoproj -W0.02  -G$ss  OUTPUT/GMT/allres.sg -K -O >> $file"
    write(600,*)"pshistogram  $geozone $geoproj -W0.02  -Gblue OUTPUT/GMT/allres.pn -K -O >> $file"
    write(600,*)"pshistogram  $geozone $geoproj -W0.02  -G$pp  OUTPUT/GMT/allres.pg -O >> $file"

    write(600,*)"ps2raster $file -Tf -A"
    write(600,'(a)')"mv OUTPUT/GMT/Allres.pdf OUTPUT/figures/Allres.pdf "
    write(600,*)
    ! -----------------------------------------------------------------    .
    val=0.0_wr
    if (FLAGresSTA) then
      nb(:)=-1
      ! -------                                                --------    .
      ! plot residus des tous les séismes à chaque station pour tous les modèles sélectionnés
      ! -------                                                --------    .
      inquire (iolength = noctet)val
      ! ---------------------------------------------------------------    .
      open(956, FILE ='OUTPUT/input/sta_new.d',status='replace')
      open(957, FILE ='OUTPUT/GMT/sta_RES_TOT2latex.txt',status='replace')
      open(958, FILE ='OUTPUT/GMT/sta_RES_TOT.txt',status='replace')
      ! ---------------------------------------------------------------    .
      ! -------  relecture des fichier et calcul de la moyenne et médiane  .
      do k=1,sta
        moy=0.0_wr
        med=0.0_wr
        ec=0.0_wr
        ! -------                                              --------    . PG
        inquire(file='OUTPUT/files/STA/'//aname(k)//'-PG.bin',exist=existe(1)) ! on teste l'existence du fichier
        ok=0
        if (existe(1)) then
          open(111, FILE ='OUTPUT/files/STA/'//aname(k)//'-PG.bin',status='old',access='direct',RECL=(noctet),iostat = ok)
          if (ok .ne. 0) then
            write(*,*)'problème dans GMT_resSTA : le fichier OUTPUT/files/STA/'//aname(k)//'-PG.bin n''existe pas'
            stop
          endif
          nb(1)=0
          do while(ok.eq.0)                                                ! boucle pour compter le nombre de lignes du fichier
            nb(1)=nb(1)+1
            read(111,rec=nb(1),iostat = ok)
          end do
          nb(1)=nb(1)-1
          allocate(vec(nb(1)))
          do j=1,nb(1)
            read(111,rec=j)vec(j)
            if (abs(vec(j)).gt.val)val=abs(vec(j)) ! garde le max pour affichage
          enddo
          close(111)
          call moy_ec(vec,nb(1),nb(1),moy(1),ec(1))
          call mediane(2.0_wr,vec,nb(1),med(1))
          ! -----------------------------------------------------------    .
          ! calcul de la droite de Henry
          ! méthode graphique pour verifier si une distribution est gaussienne
          ! -----------------------------------------------------------    .
          binage(:)=0.0_wr                                                 ! initialisation
          XY(:,:)=0.0_wr
          do j=1,nb(1)
            i=int((vec(j)+0.0005_wr)/0.001_wr,wi)
            if (i.gt.5000) i=5000
            if (i.lt.-5000) i=-5000
            binage(i)=binage(i)+1.0_wr                                     ! binage
          enddo
          binage(-5000)=binage(-5000)/real(nb(1),wr)
          do j=-4999,5000
            binage(j)=binage(j)/real(nb(1),wr)+binage(j-1)                 ! normalise et transforme en distribution cumulée
          enddo
          i=0
          do j=-5000,5000
            if ((binage(j).gt.0.15_wr).and.(binage(j).lt.0.85_wr)) then
              call inv_normal_cumulative_distrib_func(binage(j),binageR(j))! computing the inverse normal cumulative distribution function
              i=i+1
              XY(i,1)=real(j,wr)*0.001_wr-0.0005_wr
              XY(i,2)=binageR(j)
              XY(i,3)=1.0_wr                                               ! pondération constante
            endif
          enddo
          allocate (XYbis(i,3))
          do j=1,i
            XYbis(j,1)=XY(j,1)
            XYbis(j,2)=XY(j,2)
            XYbis(j,3)=XY(j,3)
          enddo
          call correlationpond(a,b,R2,i,XYbis)                             ! regession linéaire sur binageR
          do j=1,i
            XYbis(j,1)=a*XYbis(j,1)+b
          enddo
          if (i.gt.2) then
            Call Rpcalc(XYbis,i,i,Rp)                                      ! coeficient de corrélation linéaire de Bravais-Pearson
          else
            Rp=0.0_wr
          endif
          deallocate (XYbis)
          gauss(1)=Rp
          ! -----------------------------------------------------------    .
          deallocate(vec)
        else
          nb(1)=0
          moy(1)=0.0_wr
          ec(1)=0.0_wr
          gauss(1)=0.0_wr
        endif
        ! -------                                              --------    . PN
        inquire(file='OUTPUT/files/STA/'//aname(k)//'-PN.bin',exist=existe(2)) ! on teste l'existence du fichier
        ok=0
        if (existe(2)) then
          open(111, FILE ='OUTPUT/files/STA/'//aname(k)//'-PN.bin',status='old',access='direct',RECL=(noctet),iostat = ok)
          if (ok .ne. 0) then
            write(*,*)'problème dans GMT_resSTA : le fichier OUTPUT/files/STA/'//aname(k)//'-PN.bin n''existe pas'
            stop
          endif
          nb(2)=0
          do while(ok.eq.0)                                                ! boucle pour compter le nombre de lignes du fichier
            nb(2)=nb(2)+1
            read(111,rec=nb(2),iostat = ok)
          end do
          nb(2)=nb(2)-1
          allocate(vec(nb(2)))
          do j=1,nb(2)
            read(111,rec=j)vec(j)
            if (abs(vec(j)).gt.val)val=abs(vec(j)) ! garde le max pour affichage
            enddo
          close(111)
          call moy_ec(vec,nb(2),nb(2),moy(2),ec(2))
          call mediane(2.0_wr,vec,nb(2),med(2))
          ! -----------------------------------------------------------    .
          ! calcul de la droite de Henry
          ! méthode graphique pour verifier si une distribution est gaussienne
          ! -----------------------------------------------------------    .
          binage(:)=0.0_wr                                                 ! initialisation
          XY(:,:)=0.0_wr
          do j=1,nb(2)
            i=int((vec(j)+0.0005_wr)/0.001_wr,wi)
            if (i.gt.5000) i=5000
            if (i.lt.-5000) i=-5000
            binage(i)=binage(i)+1.0_wr                                     ! binage
          enddo
          binage(-5000)=binage(-5000)/real(nb(2),wr)
          do j=-4999,5000
            binage(j)=binage(j)/real(nb(2),wr)+binage(j-1)                 ! normalise et transforme en distribution cumulée
          enddo
          i=0
          do j=-5000,5000
            if ((binage(j).gt.0.15_wr).and.(binage(j).lt.0.85_wr)) then
              call inv_normal_cumulative_distrib_func(binage(j),binageR(j))! computing the inverse normal cumulative distribution function
              i=i+1
              XY(i,1)=real(j,wr)*0.001_wr-0.0005_wr
              XY(i,2)=binageR(j)
              XY(i,3)=1.0_wr                                               ! pondération constante
            endif
          enddo
          allocate (XYbis(i,3))
          do j=1,i
            XYbis(j,1)=XY(j,1)
            XYbis(j,2)=XY(j,2)
            XYbis(j,3)=XY(j,3)
          enddo
          call correlationpond(a,b,R2,i,XYbis)                             ! regession linéaire sur binageR
          do j=1,i
            XYbis(j,1)=a*XYbis(j,1)+b
          enddo
          if (i.gt.2) then
            Call Rpcalc(XYbis,i,i,Rp)                                      ! coeficient de corrélation linéaire de Bravais-Pearson
          else
            Rp=0.0_wr
          endif
          deallocate (XYbis)
          gauss(2)=Rp
          ! -----------------------------------------------------------    .
          deallocate(vec)
        else
          nb(2)=0
          moy(2)=0.0_wr
          ec(2)=0.0_wr
          gauss(2)=0.0_wr
        endif
        ! -------                                              --------    . SG
        inquire(file='OUTPUT/files/STA/'//aname(k)//'-SG.bin',exist=existe(3)) ! on teste l'existence du fichier
        ok=0
        if (existe(3)) then
          open(111, FILE ='OUTPUT/files/STA/'//aname(k)//'-SG.bin',status='old',access='direct',RECL=(noctet),iostat = ok)
          if (ok .ne. 0) then
            write(*,*)'problème dans GMT_resSTA : le fichier OUTPUT/files/STA/'//aname(k)//'-SG.bin n''existe pas'
            stop
          endif
          nb(3)=0
          do while(ok.eq.0)                                                ! boucle pour compter le nombre de lignes du fichier
            nb(3)=nb(3)+1
            read(111,rec=nb(3),iostat = ok)
          end do
          nb(3)=nb(3)-1
          allocate(vec(nb(3)))
          do j=1,nb(3)
            read(111,rec=j)vec(j)
            if (abs(vec(j)).gt.val)val=abs(vec(j)) ! garde le max pour affichage
          enddo
          close(111)
          call moy_ec(vec,nb(3),nb(3),moy(3),ec(3))
          call mediane(2.0_wr,vec,nb(3),med(3))
          ! -----------------------------------------------------------    .
          ! calcul de la droite de Henry
          ! méthode graphique pour verifier si une distribution est gaussienne
          ! -----------------------------------------------------------    .
          binage(:)=0.0_wr                                                 ! initialisation
          XY(:,:)=0.0_wr
          do j=1,nb(3)
            i=int((vec(j)+0.0005_wr)/0.001_wr,wi)
            if (i.gt.5000) i=5000
            if (i.lt.-5000) i=-5000
            binage(i)=binage(i)+1.0_wr                                     ! binage
          enddo
            binage(-5000)=binage(-5000)/real(nb(3),wr)
          do j=-4999,5000
            binage(j)=binage(j)/real(nb(3),wr)+binage(j-1)                 ! normalise et transforme en distribution cumulée
          enddo
          i=0
          do j=-5000,5000
            if ((binage(j).gt.0.15_wr).and.(binage(j).lt.0.85_wr)) then
              call inv_normal_cumulative_distrib_func(binage(j),binageR(j))! computing the inverse normal cumulative distribution function
              i=i+1
              XY(i,1)=real(j,wr)*0.001_wr-0.0005_wr
              XY(i,2)=binageR(j)
              XY(i,3)=1.0_wr                                               ! pondération constante
            endif
          enddo
          allocate (XYbis(i,3))
          do j=1,i
            XYbis(j,1)=XY(j,1)
            XYbis(j,2)=XY(j,2)
            XYbis(j,3)=XY(j,3)
          enddo
          call correlationpond(a,b,R2,i,XYbis)                             ! regession linéaire sur binageR
          do j=1,i
            XYbis(j,1)=a*XYbis(j,1)+b
          enddo
          if (i.gt.2) then
            Call Rpcalc(XYbis,i,i,Rp)                                      ! coeficient de corrélation linéaire de Bravais-Pearson
          else
            Rp=0.0_wr
          endif
          deallocate (XYbis)
          gauss(3)=Rp
          ! -----------------------------------------------------------    .
          deallocate(vec)
        else
          nb(3)=0
          moy(3)=0.0_wr
          ec(3)=0.0_wr
          gauss(3)=0.0_wr
        endif
        ! -------                                              --------    . SN
        inquire(file='OUTPUT/files/STA/'//aname(k)//'-SN.bin',exist=existe(4)) ! on teste l'existence du fichier
        ok=0
        if (existe(4)) then
          open(111, FILE ='OUTPUT/files/STA/'//aname(k)//'-SN.bin',status='old',access='direct',RECL=(noctet),iostat = ok)
          if (ok .ne. 0) then
            write(*,*)'problème dans GMT_resSTA : le fichier OUTPUT/files/STA/'//aname(k)//'-SN.bin n''existe pas'
            stop
          endif
          nb(4)=0
          do while(ok.eq.0)                                                ! boucle pour compter le nombre de lignes du fichier
            nb(4)=nb(4)+1
            read(111,rec=nb(4),iostat = ok)
          end do
          nb(4)=nb(4)-1
          allocate(vec(nb(4)))
          do j=1,nb(4)
            read(111,rec=j)vec(j)
            if (abs(vec(j)).gt.val)val=abs(vec(j)) ! garde le max pour affichage
          enddo
          close(111)
          call moy_ec(vec,nb(4),nb(4),moy(4),ec(4))
          call mediane(2.0_wr,vec,nb(4),med(4))
          ! -----------------------------------------------------------    .
          ! calcul de la droite de Henry
          ! méthode graphique pour verifier si une distribution est gaussienne
          ! -----------------------------------------------------------    .
          binage(:)=0.0_wr                                                 ! initialisation
          XY(:,:)=0.0_wr
          do j=1,nb(4)
            i=int((vec(j)+0.0005_wr)/0.001_wr,wi)
            if (i.gt.5000) i=5000
            if (i.lt.-5000) i=-5000
            binage(i)=binage(i)+1.0_wr                                     ! binage
          enddo
          binage(-5000)=binage(-5000)/real(nb(4),wr)
          do j=-4999,5000
            binage(j)=binage(j)/real(nb(4),wr)+binage(j-1)                 ! normalise et transforme en distribution cumulée
          enddo
          i=0
          do j=-5000,5000
            if ((binage(j).gt.0.15_wr).and.(binage(j).lt.0.85_wr)) then
              call inv_normal_cumulative_distrib_func(binage(j),binageR(j))! computing the inverse normal cumulative distribution function
              i=i+1
              XY(i,1)=real(j,wr)*0.001_wr-0.0005_wr
              XY(i,2)=binageR(j)
              XY(i,3)=1.0_wr                                               ! pondération constante
            endif
          enddo
          allocate (XYbis(i,3))
          do j=1,i
            XYbis(j,1)=XY(j,1)
            XYbis(j,2)=XY(j,2)
            XYbis(j,3)=XY(j,3)
          enddo
          call correlationpond(a,b,R2,i,XYbis)                             ! regession linéaire sur binageR
          do j=1,i
            XYbis(j,1)=a*XYbis(j,1)+b
          enddo
          if (i.gt.2) then
            Call Rpcalc(XYbis,i,i,Rp)                                      ! coeficient de corrélation linéaire de Bravais-Pearson
          else
            Rp=0.0_wr
          endif
          deallocate (XYbis)
          gauss(4)=Rp
          ! -----------------------------------------------------------    .
          deallocate(vec)
        else
          nb(4)=0
          moy(4)=0.0_wr
          ec(4)=0.0_wr
          gauss(4)=0.0_wr
        endif
        ymay = 1.0_wr
        val = min(val,2.0_wr)
        X = -val + 0.1_wr*val
        Y = 0.9_wr*ymay
        ! -------------------------------------------------------------    . print figure
        write(*,*)"ecriture des script GMT_res STA TOT "//aname(k)
        write(600,*)"#***************************************#"
        write(600,*)"#***************************************#"
        write(600,*)
        write(600,*)"echo 'execution du script GMT res STA TOT ",aname(k),"'"
        write(600,*)"#########################################"
        write(600,*)"###########      histo         ##########"
        write(600,*)"#########################################"
        ! -------------------------------------------------------------    .
        write(600,'(a,E13.7,a,E13.7,a,E13.7)')"geozone=-R",-val,"/",val,"/0/",ymay
        write(600,*)"geoproj=-JX4.25i"
        write(600,*)"file=OUTPUT/GMT/resSTA_TOT"//"-"//aname(k)//".ps"
        ! -------                                              --------    . PG
        if (nb(1).gt.0) then
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5g1/a0 -Bsf0.1:", &
          "'r\351sidus ondes compressives directes':/a.15:'effectif (%)':nSeW ", &
          " -X2.5i -K > $file "
          write(600,*)"echo -e ' ",moy(1)," 0.01 \n ",moy(1),ymay*.99_wr,"'| psxy $geozone $geoproj -W2,green -O -K >> $file"
          write(600,*)"echo -e ' ",moy(1)+ec(1)," 0.01 \n ",moy(1)+ec(1),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",moy(1)-ec(1)," 0.01 \n ",moy(1)-ec(1),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",med(1)," 0.01 \n ",med(1),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,blue -O -K >> $file"
          write(600,*)"pshistogram -F $geozone $geoproj -W0.001 -G$pp -Z1 OUTPUT/files/STA/"//aname(k)//"-PG.bin ", &
          "-L0/0 -bi1d -O -K >> $file"
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : ",nb(1),""" | pstext $geozone $geoproj -O -K >> $file"
          Y = 0.8_wr*ymay
          write(600,'(a,f13.5,1x,f13.5,a,f13.5,a)')"echo """,X,Y," 15 0 5 LT Rp : ",gauss(1), &
            """ | pstext $geozone $geoproj -O -K >> $file"
          Y = 0.9_wr*ymay
        else
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:", &
          "'r\351sidus ondes compressives directes':/a.15:'effectif (%)':nSeW ", &
          " -X2.5i -K > $file "
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : 0"" | pstext $geozone $geoproj -O -K >> $file"
        endif
        ! -------                                              --------    . PN
        if (nb(2).gt.0) then
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5g1/a0 -Bsf0.1:'r\351sidus ondes compressives r\351fract\351es'", &
          ":/a.15nSew -O -K -X5.5i >>  $file"
          write(600,*)"echo -e ' ",moy(2)," 0.01 \n ",moy(2),ymay*.99_wr,"'| psxy $geozone $geoproj -W2,green -O -K >> $file"
          write(600,*)"echo -e ' ",moy(2)+ec(2)," 0.01 \n ",moy(2)+ec(2),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",moy(2)-ec(2)," 0.01 \n ",moy(2)-ec(2),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",med(2)," 0.01 \n ",med(2),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,blue -O -K >> $file"
          write(600,*)"pshistogram -F $geozone $geoproj -W0.001 -G$pp -L0/0 -Ba0 ", &
          " -Z1 -K OUTPUT/files/STA/"//aname(k)//"-PN.bin -bi1d -O >> $file "
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : ",nb(2),""" | pstext $geozone $geoproj -O -K >> $file"
          Y = 0.8_wr*ymay
          write(600,'(a,f13.5,1x,f13.5,a,f13.5,a)')"echo """,X,Y," 15 0 5 LT Rp : ",gauss(2), &
            """ | pstext $geozone $geoproj -O -K >> $file"
          Y = 0.9_wr*ymay
        else
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:'r\351sidus ondes compressives r\351fract\351es'", &
          ":/a.15:'effectif (%)':nSew -O -K -X5.5i >>  $file"
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : 0"" | pstext $geozone $geoproj -O -K >> $file"
        endif
        ! -------                                              --------    . SG
        if (nb(3).gt.0) then
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5g1/a0 -Bsf0.1:", &
          "'r\351sidus ondes cisaillantes directes':/a.15:'effectif (%)':nSew -O  -K -Y5.5i >>  $file"
          write(600,*)"echo -e ' ",moy(3)," 0.01 \n ",moy(3),ymay*.99_wr,"'| psxy $geozone $geoproj -W2,green -O -K >> $file"
          write(600,*)"echo -e ' ",moy(3)+ec(3)," 0.01 \n ",moy(3)+ec(3),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",moy(3)-ec(3)," 0.01 \n ",moy(3)-ec(3),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",med(3)," 0.01 \n ",med(3),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,blue -O -K >> $file"
          write(600,*)"pshistogram -F $geozone $geoproj -W0.001 -G$ss -L0/0 -Ba0 -bi1d ", &
          " -Z1 -K OUTPUT/files/STA/"//aname(k)//"-SG.bin -O >> $file "
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : ",nb(3),""" | pstext $geozone $geoproj -O -K >> $file"
          Y = 0.8_wr*ymay
          write(600,'(a,f13.5,1x,f13.5,a,f13.5,a)')"echo """,X,Y," 15 0 5 LT Rp : ",gauss(3), &
            """ | pstext $geozone $geoproj -O -K >> $file"
          Y = 0.9_wr*ymay
        else
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:", &
          "'r\351sidus ondes cisaillantes directes':/a.15:'effectif (%)':nSew -O  -K -Y5.5i >>  $file"
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : 0"" | pstext $geozone $geoproj -O -K >> $file"
        endif
        ! -------                                              --------    . SN
        if (nb(4).gt.0) then
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5g1/a0 -Bsf0.1:'r\351sidus ondes cisaillantes r\351fract\351es'", &
          ":/a.15:'effectif (%)':nSeW -O -K -X-5.5i >>  $file"
          write(600,*)"echo -e ' ",moy(4)," 0.01 \n ",moy(4),ymay*.99_wr,"'| psxy $geozone $geoproj -W2,green -O -K >> $file"
          write(600,*)"echo -e ' ",moy(4)+ec(4)," 0.01 \n ",moy(4)+ec(4),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",moy(4)-ec(4)," 0.01 \n ",moy(4)-ec(4),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,green,-- -O -K >> $file"
          write(600,*)"echo -e ' ",med(4)," 0.01 \n ",med(4),ymay*.99_wr, &
          "'| psxy $geozone $geoproj -W2,blue -O -K >> $file"
          write(600,*)"pshistogram -F $geozone $geoproj -W0.001 -G$ss -L0/0 -Ba0 ", &
          " -Z1 OUTPUT/files/STA/"//aname(k)//"-SN.bin -O -K -bi1d >> $file "
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : ",nb(4),""" | pstext $geozone $geoproj -O -K >> $file"
          Y = 0.8_wr*ymay
          write(600,'(a,f13.5,1x,f13.5,a,f13.5,a)')"echo """,X,Y," 15 0 5 LT Rp : ",gauss(4), &
            """ | pstext $geozone $geoproj -O >> $file"
          Y = 0.9_wr*ymay
        else
          write(600,*)"psbasemap $geozone $geoproj -Bpa0.5/a0 -Bsf0.1:'r\351sidus ondes cisaillantes r\351fract\351es'", &
          ":/a.15:'effectif (%)':nSeW -O -K  -X-5.5i >>  $file"
          write(600,*)"echo """,X,Y," 15 0 5 LT donn\351es : 0"" | pstext $geozone $geoproj -O >> $file"
        endif
        ! -------                                              --------    .
        write(600,*)"ps2raster OUTPUT/GMT/resSTA_TOT-"//aname(k)//".ps -Tf -A"
        write(600,'(a)')"mv OUTPUT/GMT/resSTA_TOT-"//aname(k)//".pdf OUTPUT/figures/resSTA_TOT-"//aname(k)//".pdf"
        write(600,*)
        ! ---------------------------------------------------------------    .
        allocate(anbname(sta,4))
        anbname=0

        do i=1,nbseismes
          do j=1,nbtps(i)
            if (aname(k)==D(i)%datatps(j)%sta%staname) then
              lon=D(i)%datatps(j)%sta%lon
              lat=D(i)%datatps(j)%sta%lat
              alti=D(i)%datatps(j)%sta%alti

              if (D(i)%datatps(j)%typeonde=='G') then
                anbname(k,1)=anbname(k,1)+1
                if (D(i)%datatps(j)%andS=='S') then
                  anbname(k,3)=anbname(k,3)+1
                endif
              elseif(D(i)%datatps(j)%typeonde=='N') then
                anbname(k,2)=anbname(k,2)+1
                if (D(i)%datatps(j)%andS=='S') then
                  anbname(k,4)=anbname(k,4)+1
                endif
              else
                write(*,*)'problème dans GMT_resSTA : onde ni G ni N'
                stop
              endif


            endif
          enddo
        enddo

        ! ---------------------------------------------------------------    .

        write(957,1122)aname(k)," & \np{",moy(1),"} {\small ($\pm$ \np{",2.0_wr*ec(1),"}}) & \np{",med(1),"} \\", &
          " & \np{",moy(2),"} {\small ($\pm$ \np{",2.0_wr*ec(2),"}}) & \np{",med(2),"} \\", &
          " & \np{",moy(3),"} {\small ($\pm$ \np{",2.0_wr*ec(3),"}}) & \np{",med(3),"} \\", &
          " & \np{",moy(4),"} {\small ($\pm$ \np{",2.0_wr*ec(4),"}}) & \np{",med(4),"} \\"

        ! -------------------------------------------------------------    .
        !  au moins trois stations ...
        ! -------------------------------------------------------------    .
        if (anbname(k,1).lt.2) moy(1)=0.0_wr ! PG
        if (anbname(k,2).lt.2) moy(2)=0.0_wr ! PN
        if (anbname(k,3).lt.2) moy(3)=0.0_wr ! SG
        if (anbname(k,4).lt.2) moy(4)=0.0_wr ! SN

        if (anbname(k,1).lt.2) med(1)=0.0_wr ! PG
        if (anbname(k,2).lt.2) med(2)=0.0_wr ! PN
        if (anbname(k,3).lt.2) med(3)=0.0_wr ! SG
        if (anbname(k,4).lt.2) med(4)=0.0_wr ! SN

        ! -------------------------------------------------------------    .
        write(958,*)aname(k),lat,lon,alti,moy(1),2.0_wr*ec(1),med(1),(gauss(1)*1000.0_wr-990.0_wr)/40.0_wr, &
                                          moy(2),2.0_wr*ec(2),med(2),(gauss(2)*1000.0_wr-990.0_wr)/40.0_wr, &
                                          moy(3),2.0_wr*ec(3),med(3),(gauss(3)*1000.0_wr-990.0_wr)/40.0_wr, &
                                          moy(4),2.0_wr*ec(4),med(4),(gauss(4)*1000.0_wr-990.0_wr)/40.0_wr

        ! -------------------------------------------------------------    .
        write(956,*)aname(k),lat,lon,alti,med(1),med(2),med(3),med(4)

        deallocate(anbname)
      enddo
      close(956)
      close(957)
      close(958)
      ! ---------------------------------------------------------------    .
    endif
    ! -----------------------------------------------------------------    .
    1122 format (2a,f10.3,a,f10.2,a,f10.3,2a,f10.3,a,f10.2,a,f10.3,2a,f10.3,a,f10.2,a,f10.3,2a,f10.3,a,f10.2,a,f10.3,a)

    ! -----------------------------------------------------------------    .
  end subroutine GMT_resSTA

  ! -------------------------------------------------------------------    .


  subroutine GMT_res(l,xmax,dp,nbtps,datatps)
    ! -------                                                  --------    .mh
    ! production du script GMT produisant une carte des stations
    ! -------                                                  --------    .
    use typetemps
    use time
    implicit none
    ! -------                                                  --------    .
    real(KIND=wr), intent(in) :: xmax
    type(densityplot), intent(in) :: dp
    integer(KIND=wi), intent(in) :: nbtps, l
    type(dataone), intent(in) :: datatps(nbtps)
    ! -------                                                  --------    .
    type(sta_tps), dimension(:), allocatable :: statps
    integer(KIND=wi) :: iPn, iPg, iSn, iSg
    type(stations) :: datasta
    real(KIND=wr) :: v1, v2, amax, val, size, tl, X, Y
    real(KIND=wr) :: lat,lon,alti,moy(4),ec(4),med(4),gauss(4)
    type(date_sec) :: tps_ref, a_time
    real(KIND=wr) :: lon1,lon2,lat1,lat2
    integer(KIND=wi) :: i,j,k,n,n2,ok,Noldtime, Nnewtime, ratetime, nbsta
    character (LEN=13) :: sizename
    character (LEN=5) :: numberfile
    character (LEN=4) :: aname
    ! -----------------------------------------------------------------    .
    write(numberfile(1:5),'(i5)')l
    ! -----------------------------------------------------------------    . nombre de stations
    n=0
    n2=0
    do i=1,nbtps
      k=0
      do j=1,nbtps
        if (datatps(i)%sta%staname.eq.datatps(j)%sta%staname) then
          k=k+1
        endif
      enddo
      if (k.eq.1) n=n+1
      if (k.eq.2) n2=n2+1
      if (k.eq.3) then
        write(*,*)'problème dans GMT_res 1 : station présentes 3 fois'
        stop
      endif
    enddo
    n=n+n2/2
    allocate(statps(n))
    ! -----------------------------------------------------------------    . initialisation
    do i=1,n
      statps(i)%staname = "0000"
      statps(i)%lonSTA = 0.0_wr
      statps(i)%latSTA = 0.0_wr
      statps(i)%TPn = 0.0_wr
      statps(i)%TSn = 0.0_wr
      statps(i)%TPg = 0.0_wr
      statps(i)%TSg = 0.0_wr
      statps(i)%Ttot = 0.0_wr
      statps(i)%pdsPn = 0.0_wr
      statps(i)%pdsSn = 0.0_wr
      statps(i)%pdsPg = 0.0_wr
      statps(i)%pdsSg = 0.0_wr
      statps(i)%TpsPn = 0.0_wr
      statps(i)%TpsPg = 0.0_wr
      statps(i)%TpsSn = 0.0_wr
      statps(i)%TpsSg = 0.0_wr
      statps(i)%disthypo = 0.0_wr
    enddo
    iPn=0
    iPg=0
    iSn=0
    iSg=0
    ! -----------------------------------------------------------------    .
    i=0
    tps_ref=dp%temps_ref(l)
    tps_ref%sec=dp%Tzero(l)%vec(1)
    call basetime(tps_ref)
    donnees : do j=1,nbtps
      n2 = 1
      do k=1,n
        if (datatps(j)%sta%staname.eq.statps(k)%staname) n2 = n2 +1        ! station déjà présente ?
      enddo
      if (n2.eq.1) then                                                    ! non
        i = i +1
        statps(i)%staname = datatps(j)%sta%staname
        statps(i)%disthypo = datatps(j)%dhypo
        statps(i)%lonSTA = datatps(j)%sta%lon
        statps(i)%latSTA = datatps(j)%sta%lat
        if (datatps(j)%typeonde.eq."N") then
          iPn = iPn +1
          statps(i)%TPn = datatps(j)%dTP
          statps(i)%pdsPn = datatps(j)%wp
          a_time%date=datatps(j)%tpsR%date
          a_time%sec=datatps(j)%tpsR%secP
          call basetime(a_time)
          call difftime(statps(i)%TpsPn,a_time,tps_ref)
          if (datatps(j)%andS.eq."S") then
            iSn = iSn +1
            statps(i)%TSn = datatps(j)%dTS
            statps(i)%pdsSn = datatps(j)%ws
            a_time%date=datatps(j)%tpsR%date
            a_time%sec=datatps(j)%tpsR%secS
            call basetime(a_time)
            call difftime(statps(i)%TpsSn,a_time,tps_ref)
          endif
        elseif (datatps(j)%typeonde.eq."G") then
          iPg = iPg +1
          statps(i)%TPg = datatps(j)%dTP
          statps(i)%pdsPg = datatps(j)%wp
          a_time%date=datatps(j)%tpsR%date
          a_time%sec=datatps(j)%tpsR%secP
          call basetime(a_time)
          call difftime(statps(i)%TpsPg,a_time,tps_ref)
          if (datatps(j)%andS.eq."S") then
            iSg = iSg +1
            statps(i)%TSg = datatps(j)%dTS
            statps(i)%pdsSg = datatps(j)%ws
            a_time%date=datatps(j)%tpsR%date
            a_time%sec=datatps(j)%tpsR%secS
            call basetime(a_time)
            call difftime(statps(i)%TpsSg,a_time,tps_ref)
           endif
        else
          write(*,*)'problème dans GMT_res 1 : onde ni réfractée ni directe'
          stop
        endif
    ! -----------------------------------------------------------------    .
      elseif (n2.eq.2) then                                                ! oui
        do k=1,n
          if (datatps(j)%sta%staname.eq.statps(k)%staname) n2 = k
        enddo
        if (datatps(j)%typeonde.eq."N") then
          iPn = iPn +1
          statps(n2)%TPn = datatps(j)%dTP
          statps(n2)%pdsPn = datatps(j)%wp
          a_time%date=datatps(j)%tpsR%date
          a_time%sec=datatps(j)%tpsR%secP
          call basetime(a_time)
          call difftime(statps(n2)%TpsPn,a_time,tps_ref)
          if (datatps(j)%andS.eq."S") then
            iSn = iSn +1
            statps(n2)%TSn = datatps(j)%dTS
            statps(n2)%pdsSn = datatps(j)%ws
            a_time%date=datatps(j)%tpsR%date
            a_time%sec=datatps(j)%tpsR%secS
            call basetime(a_time)
            call difftime(statps(n2)%TpsSn,a_time,tps_ref)
          endif
        elseif (datatps(j)%typeonde.eq."G") then
          iPg = iPg +1
          statps(n2)%TPg = datatps(j)%dTP
          statps(n2)%pdsPg = datatps(j)%wp
          a_time%date=datatps(j)%tpsR%date
          a_time%sec=datatps(j)%tpsR%secP
          call basetime(a_time)
          call difftime(statps(n2)%TpsPg,a_time,tps_ref)
          if (datatps(j)%andS.eq."S") then
            iSg = iSg +1
            statps(n2)%TSg = datatps(j)%dTS
            statps(n2)%pdsSg = datatps(j)%ws
            a_time%date=datatps(j)%tpsR%date
            a_time%sec=datatps(j)%tpsR%secS
            call basetime(a_time)
            call difftime(statps(n2)%TpsSg,a_time,tps_ref)
          endif
        else
          write(*,*)'problème dans GMT_res 2 : onde ni réfractée ni directe'
          stop
        endif
      else                                                                 ! problème
        write(*,*)'problème dans GMT_res 2 : station présentes 3 fois'
        stop
      endif
    enddo donnees
    ! -----------------------------------------------------------------    .
    amax=-1.0_wr
    do i=1,n
      statps(i)%Ttot = abs(statps(i)%TPg) + abs(statps(i)%TPn) + &
      abs(statps(i)%TSg) + abs(statps(i)%TSn)
      if (amax.lt.statps(i)%TPg) amax = statps(i)%TPg
      if (amax.lt.statps(i)%TSg) amax = statps(i)%TSg
      if (amax.lt.statps(i)%TPn) amax = statps(i)%TPn
      if (amax.lt.statps(i)%TSn) amax = statps(i)%TSn
    enddo
    val= max(real(int(real(int(amax*10._wr),wr)/10._wr+1._wr),wr),1.1_wr)  ! bornes du graph
    ! -----------------------------------------------------------------    .
    open(950, FILE ='OUTPUT/residus'//"-"//trim(adjustl(numberfile))//'.d',status='replace')
    do i=1,n

        if(statps(i)%TPg.ne.0.0_wr) write(950,1001)"PG ",statps(i)%staname,statps(i)%lonSTA, &
          statps(i)%latSTA,statps(i)%TPg,statps(i)%pdsPg,abs(statps(i)%TPg*100._wr/statps(i)%TpsPg), &
          statps(i)%disthypo
        if(statps(i)%TPn.ne.0.0_wr) write(950,1001)"PN ",statps(i)%staname,statps(i)%lonSTA, &
          statps(i)%latSTA,statps(i)%TPn,statps(i)%pdsPn,abs(statps(i)%TPn*100._wr/statps(i)%TpsPn), &
          statps(i)%disthypo
        if(statps(i)%TSg.ne.0.0_wr) write(950,1001)"SG ",statps(i)%staname,statps(i)%lonSTA, &
          statps(i)%latSTA,statps(i)%TSg,statps(i)%pdsSg,abs(statps(i)%TSg*100._wr/statps(i)%TpsSg), &
          statps(i)%disthypo
        if(statps(i)%TSn.ne.0.0_wr) write(950,1001)"SN ",statps(i)%staname,statps(i)%lonSTA, &
          statps(i)%latSTA,statps(i)%TSn,statps(i)%pdsSn,abs(statps(i)%TSn*100._wr/statps(i)%TpsSn), &
          statps(i)%disthypo
    enddo
    1001 format (a,a,f10.4,f10.4,f10.4,f10.4,1x,f15.10,1x,f15.10)
    close(950)
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(*,*)"ecriture des script GMT_res "
    write(600,*)"BEFORE=$SECONDS"
    call system_clock(Noldtime)
    write(600,*)"#***************************************#"
    write(600,*)"#***************************************#"
    write(600,*)
    write(600,*)"echo 'execution du script GMT res'"
    write(600,*)"#########################################"
    write(600,*)"###########      histo         ##########"
    write(600,*)"#########################################"
    ! -----------------------------------------------------------------    .
    write(600,'(a10,E13.7,a1,E13.7,a6)')"geozone=-R",-val,"/",val,"/0/100"
    write(600,*)"geoproj=-JX4.25i"
    write(600,*)"file=OUTPUT/GMT/reshisto"//"-"//trim(adjustl(numberfile))//".ps"
    ok=0
    ! -------                                                  --------    .
    open(955, FILE = "OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot",status='replace',iostat = ok)
    do i=1,n
      if(statps(i)%TPg.ne.0.0_wr) write(955,*)statps(i)%TPg
      if(statps(i)%TPn.ne.0.0_wr) write(955,*)statps(i)%TPn
      if(statps(i)%TSg.ne.0.0_wr) write(955,*)statps(i)%TSg
      if(statps(i)%TSn.ne.0.0_wr) write(955,*)statps(i)%TSn
    enddo
    close(955)
    ! -------                                                  --------    .
    open(951, FILE = "OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".pg",status='replace',iostat = ok)
    j=0
    do i=1,n
      if(statps(i)%TPg.ne.0.0_wr) then
        write(951,*)statps(i)%TPg
        j=j+1
      endif
    enddo
    close(951)

    if(j.gt.0) then
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Ba0.0 -Z1 -K OUTPUT/GMT/", &
      "residus"//"-"//trim(adjustl(numberfile))//".pg -X2.5i > $file "
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes compressives directes':/a10:'effectif (%)':nSeW ", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -O >> $file "
    else
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes compressives directes':/a10:'effectif (%)':nSeW ", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -X2.5i > $file "
    endif
    X = -val + 0.3_wr*(2.0_wr*val)
    Y = 95.0_wr
    if(j.gt.0) write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$pp -L0/0 -Ba0 -Z1 -O -K", &
    " OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".pg >> $file "
    write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iPg,""" | pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    .
    open(952, FILE = "OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".pn",status='replace',iostat = ok)
    j=0
    do i=1,n
      if(statps(i)%TPn.ne.0.0_wr) then
        write(952,*)statps(i)%TPn
        j=j+1
      endif
    enddo
    close(952)
    if(j.gt.0) then
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Ba0.0 -Z1 -K OUTPUT/GMT/", &
      "residus"//"-"//trim(adjustl(numberfile))//".pn -O -X5.5i >> $file "
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes compressives r\351fract\351es':/a10:'effectif (%)':nSew ", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -O >> $file "
    else
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes compressives r\351fract\351es':/a10:'effectif (%)':nSew ", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -O -X5.5i >> $file "
    endif
    if(j.gt.0) write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$pp -L0/0 -Ba0 -Z1 -O -K", &
    " OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".pn >> $file "
    write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iPn,""" | pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    .
    open(953, FILE = "OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".sg",status='replace',iostat = ok)
    j=0
    do i=1,n
      if(statps(i)%TSg.ne.0.0_wr) then
        write(953,*)statps(i)%TSg
        j=j+1
      endif
    enddo
    close(953)

    if(j.gt.0) then
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Ba0.0 -Z1 -K OUTPUT/GMT/", &
      "residus"//"-"//trim(adjustl(numberfile))//".sg -O -Y5.5i >> $file "
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes cisaillantes directes':/a10:'effectif (%)':nSew", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -O >> $file "
    else
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes cisaillantes directes':/a10:'effectif (%)':nSew", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -O -Y5.5i >> $file "
    endif
    if(j.gt.0) write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$ss -L0/0 -Ba0 -Z1 -O -K", &
    " OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".sg >> $file "
    write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iSg,""" | pstext $geozone $geoproj -O -K >> $file"
    ! -------                                                  --------    .
    open(954, FILE = "OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".sn",status='replace',iostat = ok)
    j=0
      do i=1,n
    if(statps(i)%TSn.ne.0.0_wr) then
      write(954,*)statps(i)%TSn
      j=j+1
    endif
    enddo
    close(954)
    if(j.gt.0) then
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G225 -Q -Ba0.0 -Z1 -K OUTPUT/GMT/", &
      "residus"//"-"//trim(adjustl(numberfile))//".sn -O -X-5.5i  >> $file "
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes cisaillantes r\351fract\351es':/a10:'effectif (%)':nSeW ", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -O >> $file "
    else
      write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -Ggray -L0/0 -Bpa1.0g10/a0 -Bsf0.1:", &
      "'r\351sidus ondes cisaillantes r\351fract\351es':/a10:'effectif (%)':nSeW ", &
      " -Z1 -K OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".tot -O -X-5.5i >> $file "
    endif
    if(j.gt.0) write(600,*)"pshistogram -F $geozone $geoproj -W0.05 -G$ss -L0/0 -Ba0 -Z1 -O -K", &
    " OUTPUT/GMT/residus"//"-"//trim(adjustl(numberfile))//".sn  >> $file "
    write(600,*)"echo """,X,Y," 15 0 5 6 donn\351es : ",iSn,""" | pstext $geozone $geoproj -O >> $file"
    ! -------                                                  --------    .
    write(600,*)"ps2raster OUTPUT/GMT/reshisto"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/reshisto"//"-"//trim(adjustl(numberfile))//".pdf ", &
    " OUTPUT/figures/reshisto"//"-"//trim(adjustl(numberfile))//".pdf"
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(600,*)"gmtset BASEMAP_TYPE plain"
    write(600,*)"grdfile=SRC/FILES/bath1.bin"
    write(600,*)"bluef=""0/0/100"" "
    ! -------                                                  --------    .
    v1 = 2.0_wr * pi * rT / 360.0_wr                                       ! km / degree en lon
    v2 = 2.0_wr * pi * rT * sin((90.0_wr-dp%lat(l)%vec10000(1,1))/180.0_wr*pi) /360.0_wr ! km / degree en lat
    ! -------                                                  --------    .
    lon1 = dp%lon(l)%vec10000(1,1) - (xmax / v2 * .99_wr) / 2.0_wr
    lon2 = dp%lon(l)%vec10000(1,1) + (xmax / v2 * .99_wr) / 2.0_wr
    lat1 = dp%lat(l)%vec10000(1,1) - (xmax / v1 * .99_wr) / 2.0_wr
    lat2 = dp%lat(l)%vec10000(1,1) + (xmax / v1 * .99_wr) / 2.0_wr
    if (lon1.lt.-6.0_wr) lon1=-6.0_wr
    if (lat2.gt.52.0_wr) lat2=52.0_wr
    ! -------                                                  --------    .
    write(600,'(a10,E13.7,a1,E13.7,a1,E13.7,a1,E13.7)')"geozone=-R",lon1,"/",lon2,"/",lat1,"/",lat2
    write(600,'(a12,E13.7,a1,E13.7,a1,E13.7,a1,E13.7,a3,E13.7)')"geozone3d=-R",lon1,"/",lon2,"/",lat1,"/",lat2,"/0/",val
    write(600,'(a11,E13.7,a1,E13.7,a3)')"geoproj=-JC",dp%lon(l)%vec10000(1,1),"/",dp%lat(l)%vec10000(1,1),"/20.i"
    ! -----------------------------------------------------------------    .
    ! pour Pg
    ! -----------------------------------------------------------------    .
    write(600,*)"file=OUTPUT/GMT/resPg"//"-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
    " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
    ! -------                                                  --------    .
    write(600,*)"######### affiche les stations ##########"
    ok = 0                                                                 ! toutes les stations
    open(960, FILE ='DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans GMT_map : le fichier data/sta.d n''existe pas'
      stop
    endif
    do while(ok .eq. 0)
      read(960,*,iostat = ok)datasta
      if ((datasta%lon.gt.lon1).and.(datasta%lon.lt.lon2).and.(datasta%lat.gt.lat1).and.(datasta%lat.lt.lat2)) then
        write(600,*)" echo ",datasta%lon,datasta%lat, &
        " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Ggrey -Lk -Wthinnest -O -K -E200/50 >> $file"
      endif
    end do
    close(960)
    ! -------                                                  --------    .
    do i=1,nbtps                                                           ! affiche les stations utilisées
      write(600,*)" echo ",datatps(i)%sta%lon,datatps(i)%sta%lat, &
      " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Gblue -Lk -Wthinnest -O -K -E200/50 >> $file"
    enddo
    ! -------                                                  --------    .
    write(600,*)"###### Limites du massif Armoricain #####"
    write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
    ! -------                                                  --------    .
    write(600,*)"############### Epicentre ###############"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1)," \"
    write(600,*)" 0.0 | psxyz $geozone3d $geoproj -JZ2.5i -L -K -O -Wthinnest -Ggreen -Sa0.5i -E200/50 >> $file"
    ! -------                                                  --------    .
    write(600,*)"################  barres  ###############"
    do i=1,n
      size = 0.05_wr + 0.15_wr*statps(i)%pdsPg
      val = abs(statps(i)%TPg)
      write(sizename,'(E13.7)')size
      if (statps(i)%TPg.gt.0.0_wr) then
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -G$pp -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      else
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -Gblue -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      endif
    enddo
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
    ! ------- légende :                                        --------    .
    write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
    " 20 0 4 LM ondes compressives directes | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
    ! -------                                                  --------    . fin
    write(600,*)"ps2raster OUTPUT/GMT/resPg"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/resPg"//"-"//trim(adjustl(numberfile))//".pdf ", &
    "OUTPUT/figures/resPg"//"-"//trim(adjustl(numberfile))//".pdf"
    ! -----------------------------------------------------------------    .
    ! pour Pn :
    ! -----------------------------------------------------------------    .
    write(600,*)"file=OUTPUT/GMT/resPn"//"-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
    " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
    ! -------                                                  --------    .
    write(600,*)"######### affiche les stations ##########"
    ok = 0                                                                 ! toutes les stations
    open(961, FILE ='DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans GMT_map : le fichier data/sta.d n''existe pas'
      stop
    endif
    do while(ok .eq. 0)
      read(961,*,iostat = ok)datasta
      if ((datasta%lon.gt.lon1).and.(datasta%lon.lt.lon2).and.(datasta%lat.gt.lat1).and.(datasta%lat.lt.lat2)) then
        write(600,*)" echo ",datasta%lon,datasta%lat, &
        " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Ggrey -Lk -Wthinnest -O -K -E200/50 >> $file"
      endif
    end do
    close(961)
    ! -------                                                  --------    .
    do i=1,nbtps                                                           ! affiche les stations utilisées
      write(600,*)" echo ",datatps(i)%sta%lon,datatps(i)%sta%lat, &
      " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Gblue -Lk -Wthinnest -O -K -E200/50 >> $file"
    enddo
    ! -------                                                  --------    .
    write(600,*)"###### Limites du massif Armoricain #####"
    write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
    ! -------                                                  --------    .
    write(600,*)"############### Epicentre ###############"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1)," \"
    write(600,*)" 0.0 | psxyz $geozone3d $geoproj -JZ2.5i -L -K -O -Wthinnest -Ggreen -Sa0.5i -E200/50 >> $file"
    ! -------                                                  --------    .
    write(600,*)"################  barres  ###############"
    do i=1,n
      size = 0.05_wr + 0.15_wr*statps(i)%pdsPn
      val = abs(statps(i)%TPn)
      write(sizename,'(E13.7)')size
      if (statps(i)%TPn.gt.0.0_wr) then
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -G$pp -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      else
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -Gblue -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      endif
    enddo
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
    ! ------- légende :                                        --------    .
    write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
    " 20 0 4 LM 'ondes compressives r\351fract\351es' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
    ! -------                                                  --------    . fin
    write(600,*)"ps2raster OUTPUT/GMT/resPn"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/resPn"//"-"//trim(adjustl(numberfile))//".pdf ", &
    "OUTPUT/figures/resPn"//"-"//trim(adjustl(numberfile))//".pdf"
    ! -----------------------------------------------------------------    .
    ! pour Sg
    ! -----------------------------------------------------------------    .
    write(600,*)"file=OUTPUT/GMT/resSg"//"-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
    " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
    ! -------                                                  --------    .
    write(600,*)"######### affiche les stations ##########"
    ok = 0                                                                 ! toutes les stations
    open(962, FILE ='DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans GMT_map : le fichier data/sta.d n''existe pas'
      stop
    endif
    do while(ok .eq. 0)
      read(962,*,iostat = ok)datasta
      if ((datasta%lon.gt.lon1).and.(datasta%lon.lt.lon2).and.(datasta%lat.gt.lat1).and.(datasta%lat.lt.lat2)) then
        write(600,*)" echo ",datasta%lon,datasta%lat, &
        " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Ggrey -Lk -Wthinnest -O -K -E200/50 >> $file"
      endif
    end do
    close(962)
    ! -------                                                  --------    .
    do i=1,nbtps                                                           ! affiche les stations utilisées
      write(600,*)" echo ",datatps(i)%sta%lon,datatps(i)%sta%lat, &
      " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Gblue -Lk -Wthinnest -O -K -E200/50 >> $file"
    enddo
    ! -------                                                  --------    .
    write(600,*)"###### Limites du massif Armoricain #####"
    write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
    ! -------                                                  --------    .
    write(600,*)"############### Epicentre ###############"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1)," \"
    write(600,*)" 0.0 | psxyz $geozone3d $geoproj -JZ2.5i -L -K -O -Wthinnest -Ggreen -Sa0.5i -E200/50 >> $file"
    ! -------                                                  --------    .
    write(600,*)"################  barres  ###############"
    do i=1,n
      size = 0.05_wr + 0.15_wr*statps(i)%pdsSg
      val = abs(statps(i)%TSg)
      write(sizename,'(E13.7)')size
      if (statps(i)%TSg.gt.0.0_wr) then
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -G$ss -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      else
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -Gred -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      endif
    enddo
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
    ! ------- légende :                                        --------    .
    write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
    " 20 0 4 LM ondes cisaillantes directes | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
    ! -------                                                  --------    . fin
    write(600,*)"ps2raster OUTPUT/GMT/resSg"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/resSg"//"-"//trim(adjustl(numberfile))//".pdf ", &
    " OUTPUT/figures/resSg"//"-"//trim(adjustl(numberfile))//".pdf"
    ! -----------------------------------------------------------------    .
    ! pour Sn :
    ! -----------------------------------------------------------------    .
    write(600,*)"file=OUTPUT/GMT/resSn"//"-"//trim(adjustl(numberfile))//".ps"
    write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
    " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
    ! -------                                                  --------    .
    write(600,*)"######### affiche les stations ##########"
    ok = 0                                                                 ! toutes les stations
    open(963, FILE ='DATA/sta.d',status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans GMT_map : le fichier data/sta.d n''existe pas'
      stop
    endif
    do while(ok .eq. 0)
      read(963,*,iostat = ok)datasta
      if ((datasta%lon.gt.lon1).and.(datasta%lon.lt.lon2).and.(datasta%lat.gt.lat1).and.(datasta%lat.lt.lat2)) then
        write(600,*)" echo ",datasta%lon,datasta%lat, &
        " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Ggrey -Lk -Wthinnest -O -K -E200/50 >> $file"
      endif
    end do
    close(963)
    ! -------                                                  --------    .
    do i=1,nbtps                                                           ! affiche les stations utilisées
      write(600,*)" echo ",datatps(i)%sta%lon,datatps(i)%sta%lat, &
      " 0.0 | psxyz $geoproj $geozone3d -JZ2.5i -St0.1i -Gblue -Lk -Wthinnest -O -K -E200/50 >> $file"
    enddo
    ! -------                                                  --------    .
    write(600,*)"###### Limites du massif Armoricain #####"
    write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
    ! -------                                                  --------    .
    write(600,*)"############### Epicentre ###############"
    write(600,*)"echo ",dp%lon(l)%vec10000(1,1),dp%lat(l)%vec10000(1,1)," \"
    write(600,*)" 0.0 | psxyz $geozone3d $geoproj -JZ2.5i -L -K -O -Wthinnest -Ggreen -Sa0.5i -E200/50 >> $file"
    ! -------                                                  --------    .
    write(600,*)"################  barres  ###############"
    do i=1,n
      size = 0.05_wr + 0.15_wr*statps(i)%pdsSn
      val = abs(statps(i)%TSn)
      write(sizename,'(E13.7)')size
      if (statps(i)%TSn.gt.0.0_wr) then
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -G$ss -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      else
        if (val.ne.0.0_wr) write(600,*)"echo ",statps(i)%lonSTA,statps(i)%latSTA," ",val, &
        " | psxyz $geoproj $geozone3d -JZ2.5i -So"//sizename//"ib0.0 -Gred -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
      endif
    enddo
    ! -------                                                  --------    .
    write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
    ! ------- légende :                                        --------    .
    write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
    " 20 0 4 LM 'ondes cisaillantes r\351fract\351es' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
    ! -------                                                  --------    . fin
    write(600,*)"ps2raster OUTPUT/GMT/resSn"//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
    write(600,'(2a)')"mv OUTPUT/GMT/resSn"//"-"//trim(adjustl(numberfile))//".pdf ", &
    "OUTPUT/figures/resSn"//"-"//trim(adjustl(numberfile))//".pdf"
    ! -----------------------------------------------------------------    .
    deallocate(statps)
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    if ((FLAGresSTA).and.(l==nbseismes)) then
      write(600,'(a12,E13.7,a1,E13.7,a1,E13.7,a1,E13.7,a)')"geozone3d=-R",lon1,"/",lon2,"/",lat1,"/",lat2,"/0/3.0"
      nbsta=0
      ok=0
      open(959, FILE ='OUTPUT/GMT/sta_RES_TOT.txt',status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans GMT_res : le fichier OUTPUT/GMT/sta_RES_TOT.txt n''existe pas'
        stop
      endif
      do while(ok .eq. 0)                                                  ! boucle pour compter le nombre de lignes du fichier
        read(959,*,iostat = ok)
        if (ok .eq. 0) nbsta = nbsta + 1
      end do
      rewind(959)
      ok=0
      ! ---------------------------------------------------------------    .
      ! Pg
      ! ---------------------------------------------------------------    . MOY
      write(600,*)"echo 'execution du script GMT res Pg - TOT'"
      write(600,*)"file=OUTPUT/GMT/resTOTPgmoy.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
        " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)

        val = abs(moy(1))
        if (gauss(1).ne.0.0_wr) then
          if (moy(1).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(1), &
              "ib0.0 -G$pp -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(1)-ec(1)).gt.0.0_wr) then

              write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val-ec(1), &
                " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(1), &
                "ib0.0 -G$pp -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            endif
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val+ec(1), &
              " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(1), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (moy(1).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val," | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(1), &
              "ib0.0 -Gblue -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(1)-ec(1)).gt.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val-ec(1)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(1), &
              "ib0.0 -Gblue -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            if (moy(1).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val+ec(1)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(1), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                        --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
      " 20 0 4 LM 'ondes compressives directes' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTPgmoy.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTPgmoy.pdf OUTPUT/figures/resTOTPgmoy.pdf"
      ! -------                                                --------    .
      rewind(959)
      ! -------                                                --------    . MODE
      write(600,*)"file=OUTPUT/GMT/resTOTPgmed.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
        " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)
        val = abs(med(1))
        if (gauss(1).ne.0.0_wr) then
          if (med(1).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(1), &
              "ib0.0 -G$pp -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (med(1).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
            " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(1), &
            "ib0.0 -Gblue -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                        --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
        " 20 0 4 LM 'ondes compressives directes' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTPgmed.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTPgmed.pdf OUTPUT/figures/resTOTPgmed.pdf"
      ! -------                                                --------    .
      rewind(959)
      ! ---------------------------------------------------------------    .
      ! Sg
      ! ---------------------------------------------------------------    . MOY
      write(600,*)"echo 'execution du script GMT res Sg - TOT'"
      write(600,*)"file=OUTPUT/GMT/resTOTSgmoy.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
        " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)
        val = abs(moy(3))
        if (gauss(3).ne.0.0_wr) then
          if (moy(3).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(3), &
              "ib0.0 -G$ss -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(3)-ec(3)).gt.0.0_wr) then
              write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val-ec(3), &
                " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(3), &
                "ib0.0 -G$ss -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            endif
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val+ec(3), &
              " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(3), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (moy(3).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val," | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(3), &
              "ib0.0 -Gred -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(3)-ec(3)).gt.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val-ec(3)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(3), &
              "ib0.0 -Gred -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            if (moy(3).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val+ec(3)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(3), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                        --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
        " 20 0 4 LM 'ondes cisaillantes directes' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTSgmoy.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTSgmoy.pdf OUTPUT/figures/resTOTSgmoy.pdf"
      ! -------                                                --------    .
      rewind(959)
      ! -------                                                --------    . MODE
      write(600,*)"file=OUTPUT/GMT/resTOTSgmed.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
      " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)
        val = abs(med(3))
        if (gauss(3).ne.0.0_wr) then
          if (med(3).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(3), &
              "ib0.0 -G$ss -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (med(3).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(3), &
              "ib0.0 -Gred -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                        --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
      " 20 0 4 LM 'ondes cisaillantes directes' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTSgmed.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTSgmed.pdf OUTPUT/figures/resTOTSgmed.pdf"
      ! -------                                                --------    .
      rewind(959)
      ! ---------------------------------------------------------------    .
      ! Pn
      ! ---------------------------------------------------------------    . MOY
      write(600,*)"echo 'execution du script GMT res Pn - TOT'"
      write(600,*)"file=OUTPUT/GMT/resTOTPnmoy.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
        " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)
        val = abs(moy(2))
        if (gauss(2).ne.0.0_wr) then
          if (moy(2).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(2), &
              "ib0.0 -G$pp -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(2)-ec(2)).gt.0.0_wr) then
              write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val-ec(2), &
                " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(2), &
                "ib0.0 -G$pp -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            endif
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val+ec(2), &
              " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(2), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (moy(2).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val," | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(2), &
              "ib0.0 -Gblue -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(2)-ec(2)).gt.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val-ec(2)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(2), &
              "ib0.0 -Gblue -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            if (moy(2).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val+ec(2)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(2), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                        --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
        " 20 0 4 LM 'ondes compressives r\351fract\351es' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTPnmoy.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTPnmoy.pdf OUTPUT/figures/resTOTPnmoy.pdf"
      ! -------                                                --------    .
      rewind(959)
      ! -------                                                --------    . MODE
      write(600,*)"file=OUTPUT/GMT/resTOTPnmed.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
        " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)
        val = abs(med(2))
        if (gauss(2).ne.0.0_wr) then
          if (med(2).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(2), &
              "ib0.0 -G$pp -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (med(2).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(2), &
              "ib0.0 -Gblue -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                        --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
        " 20 0 4 LM 'ondes compressives r\351fract\351es' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTPnmed.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTPnmed.pdf OUTPUT/figures/resTOTPnmed.pdf"
      ! -------                                                --------    .
      rewind(959)
      ! ---------------------------------------------------------------    .
      ! Sn
      ! ---------------------------------------------------------------    . MOY
      write(600,*)"echo 'execution du script GMT res Sn - TOT'"
      write(600,*)"file=OUTPUT/GMT/resTOTSnmoy.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
        " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)
        val = abs(moy(4))
        if (gauss(4).ne.0.0_wr) then
          if (moy(4).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(4), &
              "ib0.0 -G$ss -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(4)-ec(4)).gt.0.0_wr) then
              write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val-ec(4), &
                " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(4), &
                "ib0.0 -G$ss -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            endif
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val+ec(4), &
              " | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(4), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (moy(4).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val," | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(4), &
              "ib0.0 -Gred -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
            if ((moy(4)-ec(4)).gt.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val-ec(4)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(4), &
              "ib0.0 -Gred -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
            if (moy(4).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ", &
              lon,lat," ",val+ec(4)," | psxyz $geoproj $geozone3d -JZ2.5i -SO",gauss(4), &
              "ib0.0 -Wthinner,black -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                      --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
        " 20 0 4 LM 'ondes cisaillantes r\351fract\351es' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTSnmoy.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTSnmoy.pdf OUTPUT/figures/resTOTSnmoy.pdf"
      ! -------                                                --------    .
      rewind(959)
      ! -------                                                --------    . MODE
      write(600,*)"file=OUTPUT/GMT/resTOTSnmed.ps"
      write(600,*)"pscoast $geozone3d $geoproj -JZ2.5i -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -K ", &
        " -E200/50 -Ba2f1g1/a2f1g1/a1:secondes:WSneZ -Xc -Yc >  $file"
      write(600,*)"###### Limites du massif Armoricain #####"
      write(600,*)"psxyz $geozone $geoproj -W4,gray -O SRC/FILES/limitesMA -K -E200/50 -M >> $file"
      ! -------                                                --------    .
      write(600,*)"################  barres  ###############"
      do i=1,nbsta
        read(959,*)aname,lat,lon,alti,moy(1),ec(1),med(1),gauss(1), &
                                      moy(2),ec(2),med(2),gauss(2), &
                                      moy(3),ec(3),med(3),gauss(3), &
                                      moy(4),ec(4),med(4),gauss(4)
        val = abs(med(4))
        if (gauss(4).ne.0.0_wr) then
          if (med(4).gt.0.0_wr) then
            write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(4), &
              "ib0.0 -G$ss -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          else
            if (med(4).ne.0.0_wr) write(600,'(a,f13.5,1x,f13.5,a,f13.5,a,E13.7,a)')"echo ",lon,lat," ",val, &
              " | psxyz $geoproj $geozone3d -JZ2.5i -So",gauss(4), &
              "ib0.0 -Gred -Wthinner -O -K -Ba0 -E200/50 -N >> $file"
          endif
        endif
      enddo
      ! -------                                                --------    .
      write(600,*)"psbasemap $geozone3d $geoproj -JZ2.5i -Ba0 -O -K -E200/50 >>  $file"
      ! ------- légende :                                      --------    .
      write(600,*)"echo ",lon1+(lon2-lon1)*0.05_wr,lat1+(lat2-lat1)*1.05_wr, &
        " 20 0 4 LM 'ondes cisaillantes r\351fract\351es' | pstext $geozone3d -JZ2.5i $geoproj -O -N -E200/50 >> $file"
      ! -------                                                --------    . fin
      write(600,*)"ps2raster OUTPUT/GMT/resTOTSnmed.ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/resTOTSnmed.pdf OUTPUT/figures/resTOTSnmed.pdf"
      close(959)
    endif
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
  end subroutine GMT_res

END MODULE figure_GMTres



! *********************************************************************    .
! *********************************************************************    .


