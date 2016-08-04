! permet la création des fichiers et autres scripts GMT en vue de la production des figures
! mars 2014
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE figure_GMT

    use modparam

    implicit none

    private

    public  :: GMTfull
    public  :: affiche_temps_ref

    ! -----------------------------------------------------------------    .
    ! on défini _FILE_DIR_ en fonction du compilateur -> test existance dossier

#ifdef __INTEL_COMPILER
#define _FILE_DIR_ DIRECTORY
#elif __GFORTRAN__
#define _FILE_DIR_ FILE
#endif

    ! -----------------------------------------------------------------    .

CONTAINS

! ---------------------------------------------------------------------    .

  subroutine GMTfull(dp,nmod,nbChaineMV,xmin,xmax,nbtps,nbsta,D,Ellips,nomsta,acentroid)
    ! -------                                              --------    .mh
    ! production d'un large script bash (souvent > 10 000 lignes)
    ! pour produire l'ensemble des figures sous G.M.T.
    ! -----------------------------------------------------------------    .
    use typetemps
    use time
    use pb_direct
    use sub_param
    use distance_epi
    use sub_misfit
    ! -------                                                  --------    .
    use figure_GMTwada
    use figure_GMTchat
    use figure_GMTpar
    use figure_GMTmCorr
    use figure_GMTfc
    use figure_GMTres
    use figure_GMTmap
    use figure_GMThodo
    use figure_GMTcoda
    use figure_GMTmoho_inc
    use figure_GMTcarriere
    use figure_posteriori
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent(in) :: nbtps(nbseismes),nbsta                 ! nombre de données de temps et de station
    integer(KIND=wi), intent (in) :: nbChaineMV                            ! nombre de chaîne de Markov
    integer(KIND=wi), intent (in) :: nmod(nbChaineMV)                      ! nombre de modèles retenus par chaîne de Markov
    type(densityplot), intent (inout) :: dp                                ! modèles retenus par McMC
    type(dataall), intent(inout) :: D(nbseismes)                           ! données
    real(KIND=wr), intent(in) :: xmin(nbseismes),xmax(nbseismes)           ! cercles de pondération
    type(ellip), intent(out) :: Ellips(nbseismes)                          ! ellipse
    character(LEN=4), dimension(:), allocatable, intent(out) :: nomsta
    type(amoho_centroid), intent (in) :: acentroid
    ! -------                                                  --------    .
    type(parametres) :: param_best
    type(parametre) :: aparambest
    type(pond) :: w
    integer(KIND=wi) :: i,j
    logical :: critique
    real(KIND=wr) :: delta
    integer(KIND=wi) :: triseismes(nbseismes)
    type(date_sec) :: datetriseismes(nbseismes)
    character (LEN=5) :: numberfile
    logical :: existe1
    logical :: plot=.true.
    ! -----------------------------------------------------------------    .
    open(unit=600,file="OUTPUT/GMT/script.sh",STATUS="replace")
    open(unit=601,file="OUTPUT/GMT/files.txt",STATUS="replace")
    write(600,'(a)')"echo ' '"
    write(600,'(a)')"echo '---------------------------------------------------'"
    write(600,'(a)')"echo '-------------------  script GMT -------------------'"
    write(600,'(a)')"echo '---------------------------------------------------'"
    write(600,'(a)')"echo ' '"
    write(600,'(a)')"gmtset LABEL_FONT_SIZE 15"                            ! nouvelles options GMT
    write(600,'(a)')'gmtset HEADER_FONT_SIZE 15'
    write(600,'(a)')"gmtset ANNOT_FONT_PRIMARY Times-Roman"
    write(600,'(a)')"gmtset ANNOT_FONT_SECONDARY Times-Roman"
    write(600,'(a)')"gmtset PAPER_MEDIA A3"
    write(600,'(a)')"gmtset TIME_LANGUAGE FR"
    write(600,'(a)')"gmtset TRANSPARENCY 50"
    write(600,'(a)')"gmtset PLOT_DEGREE_FORMAT dddmm"
    write(600,'(a)')"gmtset BASEMAP_TYPE fancy"
    write(600,'(a)')"gmtset CHAR_ENCODING ISOLatin1+"
    ! -----------------------------------------------------------------    . tables couleurs pour GMT
    write(600,'(a)')"echo '0	0	0	0	5	255	0	255' > OUTPUT/GMT/colorpal1.cpt"
    !write(600,'(a)')"echo '0	0	0	0	.5	1	254	255' >> OUTPUT/GMT/colorpal1.cpt"
    !write(600,'(a)')"makecpt -Ccool -T.5/5/0.01 > OUTPUT/GMT/colorpal1.cpt"  ! pour la fonction coût
    write(600,'(a)')"makecpt -Crainbow -T5/33/.01 -N -Z >> OUTPUT/GMT/colorpal1.cpt"
    write(600,'(a)')"echo '33	0	1	1	100	0	1	1' >> OUTPUT/GMT/colorpal1.cpt"
    ! -------                                                  --------    . pour le diagramme de densité
    write(600,'(a)')"makecpt -Cno_green -T1.0/97.5/.1 -N > OUTPUT/GMT/colorpal2.cpt"
    write(600,'(a)')"echo '97.5 255 0 0 100 255 0 255' >> OUTPUT/GMT/colorpal2.cpt"
    write(600,'(a)')"echo 'B 0 0 255' >> OUTPUT/GMT/colorpal2.cpt"
    write(600,'(a)')"echo 'F 255 0 0' >> OUTPUT/GMT/colorpal2.cpt"
    ! -------                                                  --------    . pour la pondération
    write(600,'(a)')"makecpt -Cno_green -T0/1/.01 -N > OUTPUT/GMT/colorpal3.cpt"
    write(600,'(a)')"pp=0/171/235"                                         ! couleur cyan
    write(600,'(a)')"ss=255/0/255"                                         ! couleur fushia
    ! -------                                                  --------    . pour le diagramme de densité du gif
    write(600,'(a)')"makecpt -Cno_green -T1/95/.1 -N > OUTPUT/GMT/colorpal5.cpt"
    write(600,'(a)')"echo 'N 255 255 255' >> OUTPUT/GMT/colorpal5.cpt"
    write(600,'(a)')"echo 'B 255 255 255' >> OUTPUT/GMT/colorpal5.cpt"
    write(600,'(a)')"echo 'F 255 0 0' >> OUTPUT/GMT/colorpal5.cpt"
    ! -----------------------------------------------------------------    . pour les tirs de carrières
    write(600,'(a)')"makecpt -Chot -T-3/3/1 > OUTPUT/GMT/colortir.cpt"
    ! -----------------------------------------------------------------    .
    ! tri des séismes dans l'ordre :
    do i=1,nbseismes
      datetriseismes(i)=dp%temps_ref(i)
    enddo
    do j=1,nbseismes
      triseismes(j)=1
      do i=1,nbseismes
        call difftime(delta,datetriseismes(i),datetriseismes(j))
        if (delta.lt.0.0_wr) triseismes(j)=triseismes(j)+1
      enddo
    enddo
    j=-1
    do i=1,nbseismes
      if (triseismes(i)==1)j=i
    enddo
    ! -----------------------------------------------------------------    . modèle de référence : best modèle
    !param_best%VC=dp%VC%best
    !param_best%VM=dp%VM%best
    !param_best%Zmoho=dp%Zmoho%best
    !param_best%VpVs=dp%VpVs%best
    !do i=1,nbseismes
    !  param_best%Zhypo(i)=dp%Zhypo(i)%best
    !  param_best%lon(i)=dp%lon(i)%best
    !  param_best%lat(i)=dp%lat(i)%best
    !  param_best%Tzero(i) = dp%temps_ref(i)
    !  param_best%Tzero(i)%sec = dp%Tzero(i)%best
    !  call basetime(param_best%Tzero(i))
    !enddo
    ! (à modifier aussi dans mkwada.f90 et mklatex.f90)
    ! -------                                                  --------    . modèle de référence : 100 best modèle
    param_best%VC=dp%VC%moy_100
    param_best%VM=dp%VM%moy_100
    param_best%Zmoho=dp%Zmoho%moy_100
    param_best%VpVs=dp%VpVs%moy_100
    do i=1,nbseismes
      param_best%Zhypo(i)=dp%Zhypo(i)%moy_100
      param_best%lon(i)=dp%lon(i)%moy_100
      param_best%lat(i)=dp%lat(i)%moy_100
      param_best%Tzero(i) = dp%temps_ref(i)
      param_best%Tzero(i)%sec = dp%Tzero(i)%moy_100
      call basetime(param_best%Tzero(i))
    enddo
    ! -------                                                  --------    . modèle de référence : un modèle fixe
    !param_best%VC=6.0_wr
    !param_best%VM=8.0_wr
    !param_best%Zmoho=30.0_wr
    !param_best%VpVs=1.71_wr
    !do i=1,nbseismes
    !  param_best%Zhypo(i)=15.0_wr
    !  param_best%lon(i)=-2.25_wr
    !  param_best%lat(i)=48.25_wr
    !  param_best%Tzero(i) = dp%temps_ref(i)
    !  param_best%Tzero(i)%sec = 30.0_wr
    !  call basetime(param_best%Tzero(i))
    !enddo ! (à modifier aussi dans mkwada.f90 et mklatex.f90)
    ! -----------------------------------------------------------------    .
    !write(*,*)param_best
    ! -----------------------------------------------------------------    .
    call tempsTheoDirect(nbtps,param_best,D,critique,acentroid)
    call mvPall_2_P1(aparambest,param_best,1)
    ! -----------------------------------------------------------------    . production des scripts pour chaque couple de parametres
    i=0
    ! -----------------------------------------------------------------    . pour VC versus VpVs
    if(plot) call GMT_2paramplot(i,dp%VC,dp%VpVs,dp%mis,dp%deltaxy,dp%nbparam,nbChaineMV,dp%temps_ref,nmod,aparambest)
    ! -----------------------------------------------------------------    . pour VM versus Zmoho
    if(plot) call GMT_2paramplot(i,dp%VM,dp%Zmoho,dp%mis,dp%deltaxy,dp%nbparam,nbChaineMV,dp%temps_ref,nmod,aparambest)
    ! -----------------------------------------------------------------    .
    do i=1,nbseismes
      call  mvPall_2_P1(aparambest,param_best,i)
      ! ---------------------------------------------------------------    . pour lon versus lat (format carte)
      if(plot) call GMT_2paramplot(i,dp%lon(i),dp%lat(i),dp%mis,dp%deltaxy,dp%nbparam, &
        nbChaineMV,dp%temps_ref,nmod,aparambest,t=triseismes,E=Ellips(i))
      ! ---------------------------------------------------------------    . pour Zhypo versus Tzero
      if(plot) call GMT_2paramplot(i,dp%Zhypo(i),dp%Tzero(i),dp%mis,dp%deltaxy,dp%nbparam,nbChaineMV,dp%temps_ref,nmod,aparambest)
      ! ---------------------------------------------------------------    .
      !if(plot) call GMT_2paramplot(i,dp%lon(i),dp%Zhypo(i),dp%mis,dp%deltaxy,dp%nbparam,nbChaineMV,dp%temps_ref,nmod,aparambest)
      !if(plot) call GMT_2paramplot(i,dp%Zhypo(i),dp%lat(i),dp%mis,dp%deltaxy,dp%nbparam,nbChaineMV,dp%temps_ref,nmod,aparambest)
      !if(plot) call GMT_2paramplot(i,dp%Zhypo(i),dp%Zmoho,dp%mis,dp%deltaxy,dp%nbparam,nbChaineMV,dp%temps_ref,nmod,aparambest)
      !if(plot) call GMT_2paramplot(i,dp%Tzero(i),dp%Zmoho,dp%mis,dp%deltaxy,dp%nbparam,nbChaineMV,dp%temps_ref,nmod,aparambest)
      ! ---------------------------------------------------------------    .
    enddo
    ! -----------------------------------------------------------------    . production d'autres scripts
    if((plot).and.(FLAGresSTA)) call GMT_resSTA(nbsta,nbtps,D,nomsta)     
    ! -------                                                  --------    .
    do i=1,nbseismes
      write(numberfile(1:5),'(i5)')i
      call mvPall_2_P1(aparambest,param_best,i)
      call ponderation(nbtps(i),D(i)%datatps,xmin(i),xmax(i),w)
      ! -------                                                --------    .
      if((plot).and.(plotposteriori)) then
        call GMT_posteriori_lonlat(i,dp)
        call GMT_posteriori(i,dp,dp%VC)
        call GMT_posteriori(i,dp,dp%VM)
        call GMT_posteriori(i,dp,dp%VpVs)
        call GMT_posteriori(i,dp,dp%Zmoho)
        call GMT_posteriori(i,dp,dp%Zhypo(i))
        call GMT_posteriori(i,dp,dp%Tzero(i))
      endif
      ! -------                                                --------    .
      ! if((plot).and.(plotposteriori)) call GMT_posteriori ... atres plots
      ! -------                                                --------    .
      if(plot) call GMT_map(i,xmin(i),xmax(i),dp,nbtps(i),D(i)%datatps)
      ! -------                                                --------    .
      if((plot).and.(FLAG_non_tabulaire)) call GMT_moho(acentroid,i,nbtps(i),xmax(i),D(i)%datatps,aparambest)
      ! -------                                                --------    .
      if(plot) call GMT_res(i,xmax(i),dp,nbtps(i),D(i)%datatps)
      ! -------                                                --------    .
      if(plot) call GMT_Hodochrone(i,nbtps(i),D(i)%datatps,aparambest,xmax(i),acentroid)
      ! -------                                                --------    .
      if(plot) then
        inquire (_FILE_DIR_="DATA/sac-"//trim(adjustl(numberfile)),exist=existe1) ! option différente selon compilo !
        if ((existe1).and.(tracessac)) call GMT_coda(i,nbtps(i),D(i)%datatps, &
        aparambest,xmax(i),dp%lon(i)%vec10000(1,1),dp%lat(i)%vec10000(1,1),acentroid)
      endif
      ! -------                                                --------    .
      if(plot) call GMT_carriere(i,dp%lon(i)%vec10000(1,1),dp%lat(i)%vec10000(1,1),dp%temps_ref(i))
      ! -------                                                --------    .
    enddo
    ! -------                                                  --------    .
    if(plot) call GMT_chatelain(nbtps,nbsta,D,dp)
    ! -------                                                  --------    .
    if(plot) call GMT_wadati(nbtps,D,param_best,dp)
    ! -------                                                  --------    .
    if((plot).and.(plotgraph)) call GMT_fc(dp,nbChaineMV,nmod)
    ! -------                                                  --------    .
    if((plot).and.(plotgraph)) call GMT_param(dp,nbChaineMV,nmod)
    ! -------                                                  --------    .
    if(plot) call GMT_mCorr(dp)
    ! -----------------------------------------------------------------    .
    write(600,'(a8)')"echo ' '"
    write(600,'(a58)')"echo '---------------------------------------------------'"
    write(600,'(a58)')"echo '-----------------  fin script GMT -----------------'"
    write(600,'(a58)')"echo '---------------------------------------------------'"
    write(600,'(a8)')"echo ' '"
    write(600,'(a)')"rm -rf OUTPUT/GMT/*.ps OUTPUT/GMT/*.eps"
    close(600)
    close(601)
    ! -----------------------------------------------------------------    .

    CONTAINS

      ! ---------------------------------------------------------------    .

    subroutine GMT_2paramplot(mm,param1,param2,mis,deltaxy,nbparam,nbChaineMV,temps_ref,nmod,param_best,t,E)
      ! -------                                                --------    .mh
      ! production d'une partie du script GMT pour le diagramme de densité,
      ! les distributions de probabilités marginalles et la représentation 2D de la fonction coût,
      ! d'un couple de paramètres (param1,param2)
      ! script différent si param1 = lon et param2 = lat -> carte
      ! ---------------------------------------------------------------    .
      use typetemps
      use time
      use datalecture
      ! -------                                                --------    .
      implicit none
      ! -------                                                --------    .
      integer(KIND=wi), intent (in) :: mm
      type(densityplot_one), intent (inout) :: param1, param2, mis         ! les deux paramètres
      type(date_sec), intent (in) :: temps_ref(nbseismes)                  ! temps de réference (si param1 = Tzéro ou param2 = Tzéro)
      integer(KIND=wi), intent (in) :: deltaxy                             ! nombre de discrétisations pour le mode et le diagramme de densité
      integer(KIND=wi), intent (in) :: nbparam                             ! nombre de modèles
      integer(KIND=wi), intent (in) :: nbChaineMV                          ! nombre de chaînes
      integer(KIND=wi), intent (in) :: nmod(nbChaineMV)                    ! nombre de modèles par chaîne
      ! -------                                                --------    .
      type(ellip), intent(out), optional :: E                              ! ellipse
      integer(KIND=wi), intent (in), optional :: t(nbseismes)
      ! -------                                                --------    .
      integer(KIND=wi) :: i,j,k,l,m,ok
      ! -------                                                --------    .
      integer(KIND=wi), parameter :: deltaxymis = 1000                     ! nombre de discrétisations pour la représentation 2D de la fonction coût
      type(stations) :: datasta                                            ! propriétés d'une une station
      ! -------                                                --------    .
      real(KIND=wr) :: diff1,diff2,d_diff                                  ! pour les échelles
      real(KIND=wr) :: themax,themin,minmax1,minmax2,val,val1,val2,X,Y
      real(KIND=wr) :: delta_1,delta_2                                     ! pas de discrétisation pour la représentation 2D de la fonction coût
      ! -------                                                --------    . quelques chaînes de caractères
      real(KIND=wr) :: moyBAZ,sumBAZ,baz,p_a,p_b,p_c,dist,bazV(360),bazVbis(360)
      real(KIND=wr) :: tl
      ! -------                                                --------    .
      integer(KIND=wi) :: Noldtime, Nnewtime, ratetime
      integer(KIND=wi) :: find                                             ! séisme trouvé dans le catalogue (find=1 ou 2)
      ! -------                                                --------    .
      character(LEN=30) :: char_0
      character(LEN=13) :: char_1,char_2,char_3,char_4
      character(LEN=10) :: char_map(4)
      character(LEN=5) :: char_5
      character(LEN=2) :: char_6
      character(LEN=7) :: filename                                         ! base du nom des fichiers de sorties
      type(parametre) :: param_best
      logical :: existe1, existe2, existe3, existe4, findtest              ! a priori dispo (existe), séisme trouvé dans le catalogue
      type(seismes) :: refseisme(2)
      character (LEN=5) :: numberfile
      ! ---------------------------------------------------------------    . initialisation
      findtest=.false.
      ! ---------------------------------------------------------------    .
      m=mm
      filename=param1%name//"_"//param2%name                               ! nom des fichiers de sorties
      if (m==0) then
      m=1
        write(numberfile(1:5),'(i5)')m
        do i=1,nbseismes
          write(601,'(i6,1x,a)')i,filename//"-"//trim(adjustl(numberfile))//".pdf"
        enddo
      else
      write(numberfile(1:5),'(i5)')m
      write(601,'(i6,1x,a)')m,filename//"-"//trim(adjustl(numberfile))//".pdf"
      endif
      ! ---------------------------------------------------------------    .
      call mkdensityplot(m,param1,param2,deltaxy,nbparam,filename)         ! créer la grille de densité
      ! ---------------------------------------------------------------    .
      ! identifie les points permettant la représentation de la fonction coût
      call mkfcoutplot(m,param1,param2,mis,deltaxymis,nbparam,filename,delta_1,delta_2)
      ! ---------------------------------------------------------------    .
      ! ------- ecriture dans un fichier du vecteur 1, 2 et misfit  ---    .
      open(unit=603,file="OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_v12.bin", &
           STATUS="replace",access='direct',RECL=24)
      do i=1,10000
        write(603,rec=i)real(param1%vec10000(i,1),8),real(param2%vec10000(i,1),8),real(param1%vec10000(i,2),8)
      enddo
      close(603)
      ! ------- ecriture dans un fichier du vecteur 1, 2 et misfit  ---    .
      open(unit=604,file="OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_tot.bin", &
           STATUS="replace",access='direct',RECL=24)
      do i=1,nbparam
        write(604,rec=i)real(param1%vec(i),8),real(param2%vec(i),8),real(mis%vec(i),8)
      enddo
      close(604)
      ! ------- ecriture du script GMT                         --------    .
      write(*,'(2a)')" ecriture du script GMT pour ",filename
      write(600,'(a)')"BEFORE=$SECONDS"
      call system_clock(Noldtime)
      write(600,'(a)')"#***************************************#"
      write(600,'(a)')"#***************************************#"
      write(600,*)
      write(600,'(2a)')"echo 'execution du script GMT pour '",filename
      write(600,'(a)')"#########################################"
      write(600,'(3a)')"#####     density plot : ",filename,"     ####"
      write(600,'(a)')"#########################################"
      write(numberfile(1:5),'(i5)')m
      write(600,'(a)')"file=OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".ps"

      ! ---------------------------------------------------------------    .      
      ! ---------------------------------------------------------------    .
      ! ------- production du script si différent d'une carte  --------    .
      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .

      if (.not.((param1%name.eq."lon").and.(param2%name.eq."lat"))) then
        write(600,'(a)')"geoproj=-JX5i/5i"                                 ! système de projection
        write(600,'(a10,E13.7,a1,E13.7,a1,E13.7,a1,E13.7)')"geozone=-R",param1%themin,"/",param1%themax,"/", &
        param2%themin,"/",param2%themax                                    ! bornes minimales et maximales
        write(600,'(a)')"############## xyz -> grid #############"
        ! ------- grille pour le diagramme de densité          --------    .
        if ((param1%name.eq."lon").and.(param2%name.eq."_zh")) then
          write(600,'(a19,E13.7,a1,E13.7,a2)')"xyz2grd $geozone -I",param1%delta*1.01_wr,"/",param2%delta*1.01_wr," \"
          write(600,'(a)')" OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".bin -bi3d -Nnan \"
          write(600,'(a)')" -GOUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd "
        elseif ((param1%name.eq."_zh").and.(param2%name.eq."lat")) then
          write(600,'(a19,E13.7,a1,E13.7,a2)')"xyz2grd $geozone -I",param1%delta*1.01_wr,"/",param2%delta*1.01_wr," \"
          write(600,'(a)')" OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".bin -bi3d -Nnan \"
          write(600,'(a)')" -GOUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd "
        else
          write(600,'(a19,E13.7,a1,E13.7,a2)')"xyz2grd $geozone -I",param1%delta*1.01_wr,"/",param2%delta*1.01_wr," \"
          write(600,'(a)')" OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".bin -bi3d -Nnan \"
          write(600,'(a)')" -GOUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd "
        endif
        ! ------- grille pour la fonction coût                 --------    .
        write(600,'(a19,E13.7,a1,E13.7,a2)')"xyz2grd $geozone -I",delta_1*1.5_wr,"/",delta_2*1.5_wr," \"
        write(600,'(a)')" OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"mis.bin -bi3d -Nnan \"
        write(600,'(a)')" -GOUTPUT/GMT/topo2_"//trim(adjustl(numberfile))//".grd "
        ! -------------------------------------------------------------    .
        ! ------- choix des échelles et incréments             --------    .
        d_diff=(param1%themax-param1%themin)/4.5_wr
        if (d_diff.gt.0.0001_wr) diff1=real(int(d_diff*100000.0_wr,wi),wr)/100000.0_wr
        if (d_diff.gt.0.001_wr) diff1=real(int(d_diff*10000.0_wr,wi),wr)/10000.0_wr
        if (d_diff.gt.0.001_wr) diff1=real(int(d_diff*1000.0_wr,wi),wr)/1000.0_wr
        if (d_diff.gt.0.01_wr) diff1=real(int(d_diff*100.0_wr,wi),wr)/100.0_wr
        if (d_diff.gt.0.1_wr) diff1=real(int(d_diff*10.0_wr,wi),wr)/10.0_wr
        if (d_diff.gt.1.0_wr) diff1=real(int(d_diff*1.0_wr,wi),wr)/1.0_wr
        if (d_diff.gt.10.0_wr) diff1=real(int(d_diff*0.1_wr,wi),wr)/0.1_wr
        if (d_diff.gt.100.0_wr) diff1=real(int(d_diff*0.01_wr,wi),wr)/0.01_wr
        write(char_1,'(E13.7)')diff1
        !write(*,*)d_diff,dp%st_1,diff1
        ! -------                                              --------    .
        d_diff=(param2%themax-param2%themin)/4.5_wr
        if (d_diff.gt.0.0001_wr) diff2=real(int(d_diff*100000.0_wr,wi),wr)/100000.0_wr
        if (d_diff.gt.0.001_wr) diff2=real(int(d_diff*10000.0_wr,wi),wr)/10000.0_wr
        if (d_diff.gt.0.001_wr) diff2=real(int(d_diff*1000.0_wr,wi),wr)/1000.0_wr
        if (d_diff.gt.0.01_wr) diff2=real(int(d_diff*100.0_wr,wi),wr)/100.0_wr
        if (d_diff.gt.0.1_wr) diff2=real(int(d_diff*10.0_wr,wi),wr)/10.0_wr
        if (d_diff.gt.1.0_wr) diff2=real(int(d_diff*1.0_wr,wi),wr)/1.0_wr
        if (d_diff.gt.10.0_wr) diff2=real(int(d_diff*0.1_wr,wi),wr)/0.1_wr
        if (d_diff.gt.100.0_wr) diff2=real(int(d_diff*0.01_wr,wi),wr)/0.01_wr
        write(char_2,'(E13.7)')diff2
        ! -------------------------------------------------------------    .
        write(char_3,'(E13.7)')param1%delta
        write(char_4,'(E13.7)')param2%delta
        ! -------------------------------------------------------------    .
        write(600,'(a)')"#########################################"
        write(600,'(a)')"############ diag. fct coût #############"
        write(600,'(a)')"#########################################"
        write(600,'(2a)')"psbasemap $geozone $geoproj -Ba"//char_1//":"""//param1%char//""":/a"//char_2//":"""// &
        param2%char//""":WenS -K -X2.5i -Yc >  $file"
        write(600,'(a)')"########### affiche les points ##########"
        write(600,'(2a)')"psxy $geozone $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"mis.bin -bi3d", &
        " -Sc0.015i -COUTPUT/GMT/colorpal1.cpt -O -K >> $file"
        write(600,'(a)')"############ contour densité ############"
        ! -------------------------------------------------------------    .
        if ((param1%name.eq."lon").and.(param2%name.eq."_zh")) then
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd $geozone ", &
          "$geoproj -Ba0 -C75 -L74/75 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd $geozone ", &
          "$geoproj -Ba0 -C50 -L49/50 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd $geozone ", &
          "$geoproj -Ba0 -C25 -L24/25 -A+s15 -W2 -O -K >> $file"
        elseif ((param1%name.eq."_zh").and.(param2%name.eq."lat")) then
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd $geozone ", &
          "$geoproj -Ba0 -C75 -L74/75 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd $geozone ", &
          "$geoproj -Ba0 -C50 -L49/50 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd $geozone ", &
          "$geoproj -Ba0 -C25 -L24/25 -A+s15 -W2 -O -K >> $file"
        else
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd $geozone ", &
          "$geoproj -Ba0 -C75 -L74/75 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd $geozone ", &
          "$geoproj -Ba0 -C50 -L49/50 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd $geozone ", &
          "$geoproj -Ba0 -C25 -L24/25 -A+s15 -W2 -O -K >> $file"
        endif
        ! -------------------------------------------------------------    .
        write(600,'(a)')"###### moy modèles des chaînes  1 #######"
        ! ------- affiche la moyennes de l'ensemble du meilleur modèle de chaque chaîne
        minmax1 = param1%moy_bestchaine - param1%ec_bestchaine
        minmax2 = param1%moy_bestchaine + param1%ec_bestchaine
        themin = param2%themin  + 0.035_wr * (param2%themax - param2%themin)
        themax = param2%themax  - 0.035_wr * (param2%themax - param2%themin)
        write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%moy_bestchaine,themin," \n ",param1%moy_bestchaine,themax, &
        """ | psxy $geozone $geoproj -W3,grey -O -K >> $file"
        minmax1 = param2%moy_bestchaine - param2%ec_bestchaine
        minmax2 = param2%moy_bestchaine + param2%ec_bestchaine
        themin = param1%themin  + 0.035_wr * (param1%themax - param1%themin)
        themax = param1%themax  - 0.035_wr * (param1%themax - param1%themin)
        write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,themin,param2%moy_bestchaine," \n ",themax,param2%moy_bestchaine, &
        """ | psxy $geozone $geoproj -W3,grey -O -K >> $file"
        write(600,'(a)')"######### 10000 meilleurs modèles #######"
        ! ------- affiche les 10000 meilleurs modèles          --------    .
        ! write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_v12.bin -Sc0.001i -O -K -bi3d >> $file"
        write(600,'(a)')"########### barre de couleur ############"
        ! ------- affiche la barre de couleur                  --------    .
        write(600,'(2a)')"psscale -D1/-1/0.50E+01/0.25ch -B25.0:""fonction co\373t"": -S -I ", &
        "-COUTPUT/GMT/colorpal1.cpt -O -K >> $file"
        write(600,'(a)')"######### 10 meilleurs modèles ##########"
        ! ------- affiche les 10 meilleurs modèles (tous différents) sous la forme d'étoiles jaunes.
        i=1
        l=1
        write(600,'(a,2f15.5,a)')"echo ",param1%vec10000(i,1),param2%vec10000(i,1), &
        "| psxy $geozone $geoproj  -Sa0.1i -Gyellow -Wthinnest,black -O -K >> $file"
        do while(l.lt.10)
          i=i+1
          if (mis%vec10000(i,1).ne.mis%vec10000(i-1,1)) then
            l=l+1
            write(600,'(a,2f15.5,a)')"echo ",param1%vec10000(i,1),param2%vec10000(i,1), &
            "| psxy $geozone $geoproj  -Sa0.1i -Gyellow -Wthinnest,black -O -K >> $file"
          endif
        enddo
        ! -------                                              --------    .
        write(600,'(a)')"psbasemap $geozone $geoproj -Ba0 -K -O >>  $file"
        write(600,'(a)')"###### moy des meilleurs modèles ########"
        ! ------- affiche la moyenne +ou- un ecart-type des 10000, 1000 et 100 meilleurs modèles, avec des flêches sur les côtés
        minmax1 = 0.015_wr * (param1%themax - param1%themin)
        minmax2 = 0.015_wr * (param2%themax - param2%themin)
        ! -------                                              --------    .
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_10000+param1%ec_10000, &
        param2%themin+minmax2,param1%moy_10000-param1%ec_10000,param2%themin+minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themin+minmax1, &
        param2%moy_10000+param2%ec_10000,param1%themin+minmax1,param2%moy_10000-param2%ec_10000, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_1000+param1%ec_1000, &
        param2%themin+minmax2,param1%moy_1000-param1%ec_1000,param2%themin+minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themin+minmax1, &
        param2%moy_1000+param2%ec_1000,param1%themin+minmax1,param2%moy_1000-param2%ec_1000, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_100+param1%ec_100, &
        param2%themin+minmax2,param1%moy_100-param1%ec_100,param2%themin+minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themin+minmax1, &
        param2%moy_100+param2%ec_100,param1%themin+minmax1,param2%moy_100-param2%ec_100, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        ! -------                                              --------    .
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_10000+param1%ec_10000, &
        param2%themax-minmax2,param1%moy_10000-param1%ec_10000,param2%themax-minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themax-minmax1, &
        param2%moy_10000+param2%ec_10000,param1%themax-minmax1, param2%moy_10000-param2%ec_10000, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_1000+param1%ec_1000, &
        param2%themax-minmax2,param1%moy_1000-param1%ec_1000,param2%themax-minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themax-minmax1, &
        param2%moy_1000+param2%ec_1000,param1%themax-minmax1,param2%moy_1000-param2%ec_1000, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_100+param1%ec_100, &
        param2%themax-minmax2,param1%moy_100-param1%ec_100,param2%themax-minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themax-minmax1, &
        param2%moy_100+param2%ec_100,param1%themax-minmax1,param2%moy_100-param2%ec_100, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        ! ------- affiche la moyennes de l'ensemble du meilleur modèle de chaque chaîne
        write(600,'(a)')"###### moy modèles des chaînes 2 #######"
        minmax1 = param1%moy_bestchaine - param1%ec_bestchaine
        minmax2 = param1%moy_bestchaine + param1%ec_bestchaine
        themin = param2%themin  + 0.035_wr * (param2%themax - param2%themin)
        themax = param2%themax  - 0.035_wr * (param2%themax - param2%themin)
        write(600,'(a,2f15.5,a)')"echo ",minmax1,themin," 0 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",minmax2,themin," 0 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",param1%moy_bestchaine,themin," 0 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.05i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",minmax1,themax," 180 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",minmax2,themax," 180 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",param1%moy_bestchaine,themax," 180 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        minmax1 = param2%moy_bestchaine - param2%ec_bestchaine
        minmax2 = param2%moy_bestchaine + param2%ec_bestchaine
        themin = param1%themin  + 0.035_wr * (param1%themax - param1%themin)
        themax = param1%themax  - 0.035_wr * (param1%themax - param1%themin)
        write(600,'(a,2f15.5,a)')"echo ",themin,minmax1," 90 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",themin,minmax2," 90 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",themin,param2%moy_bestchaine," 90 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.05i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",themax,minmax1," 270 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",themax,minmax2," 270 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,'(a,2f15.5,a)')"echo ",themax,param2%moy_bestchaine," 270 0.1i | psxy $geozone $geoproj \"
        write(600,'(a)')"-SV0.04i/0.06i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        ! ------- affiche le temps de réference                --------    .
        if((param1%name.eq."_to").or.(param2%name.eq."_to")) then
          write(600,'(a)')"############ temps réference ############"
          call affiche_temps_ref(temps_ref(m),char_0,-1)
          X = param1%themin + 0.3_wr*(param1%themax-param1%themin)
          Y = param2%themin + 0.05_wr*(param2%themax-param2%themin)
          write(600,'(a,2f15.5,a)')"echo """,X,Y," 15 0 5 6 "//char_0//""" | pstext $geozone $geoproj -O -K >> $file"
        endif
        write(600,'(a)')"#########################################"
        write(600,'(a)')"############# diag. densité #############"
        write(600,'(a)')"#########################################"
        write(600,'(2a)')"psbasemap $geozone $geoproj -Ba"//char_1//":"""//param1%char//""":/a"//char_2//":"""// &
        param2%char//""":WenS -K -O -X6.25i >>  $file"
        write(600,'(a)')"############### grid image ##############"
        ! -------                                              --------    .
        if ((param1%name.eq."lon").and.(param2%name.eq."_zh")) then
          write(600,'(2a)')"grdimage $geozone $geoproj OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd ",&
          "-Qnan -COUTPUT/GMT/colorpal2.cpt -B0 -O -K -Sn >> $file"
          write(600,'(a)')"############## grid contour #############"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd ", &
          "$geozone $geoproj -Ba0 -C75 -L74/75 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd ", &
          "$geozone $geoproj -Ba0 -C50 -L49/50 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".lon_zh.grd ", &
          "$geozone $geoproj -Ba0 -C25 -L24/25 -A+s15 -W2 -O -K >> $file"
        elseif ((param1%name.eq."_zh").and.(param2%name.eq."lat")) then
          write(600,'(2a)')"grdimage $geozone $geoproj OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd ",&
          "-Qnan -COUTPUT/GMT/colorpal2.cpt -B0 -O -K -Sn >> $file"
          write(600,'(2a)')"############## grid contour #############"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd ", &
          "$geozone $geoproj -Ba0 -C75 -L74/75 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd ", &
          "$geozone $geoproj -Ba0 -C50 -L49/50 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".zh_lat.grd ", &
          "$geozone $geoproj -Ba0 -C25 -L24/25 -A+s15 -W2 -O -K >> $file"
        else
          write(600,'(2a)')"grdimage $geozone $geoproj OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd ",&
          "-Qnan -COUTPUT/GMT/colorpal2.cpt -B0 -O -K -Sn >> $file"
          write(600,'(a)')"############## grid contour #############"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd ", &
          "$geozone $geoproj -Ba0 -C75 -L74/75 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd ", &
          "$geozone $geoproj -Ba0 -C50 -L49/50 -W2 -O -K >> $file"
          write(600,'(2a)')"grdcontour OUTPUT/GMT/topo1_"//trim(adjustl(numberfile))//".grd ", &
          "$geozone $geoproj -Ba0 -C25 -L24/25 -A+s15 -W2 -O -K >> $file"
        endif
        ! -------                                              --------    .
        write(600,'(a)')"psbasemap $geozone $geoproj -Ba0 -K -O >>  $file"
        write(600,'(a)')"psscale -D1/-1/0.50E+01/0.25ch -B25.0:""densit\351"": -S -I -COUTPUT/GMT/colorpal2.cpt -O -K >> $file"
        call catalogue(param_best,refseisme,find)
        ! -------                                              --------    . 1 séisme
        if ((find==1).or.(find==2)) then
          findtest=.false.
          write(600,'(a)')"################## catalogue 1 ################"
          if (param1%name.eq."lon") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,refseisme(1)%lon,param2%themin+0.001_wr," \n ", &
            refseisme(1)%lon,param2%themax-0.001_wr," "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."lon") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,refseisme(1)%lon," \n ", &
            param1%themax-0.001_wr,refseisme(1)%lon, " "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param1%name.eq."lat") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,refseisme(1)%lat,param2%themin+0.001_wr," \n ",refseisme(1)%lat, &
            param2%themax-0.001_wr, " "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."lat") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,refseisme(1)%lat," \n ",param1%themax-0.001_wr, &
            refseisme(1)%lat, " "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param1%name.eq."_zh") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,refseisme(1)%pfd,param2%themin+0.001_wr," \n ",refseisme(1)%pfd, &
            param2%themax-0.001_wr, " "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."_zh")then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,refseisme(1)%pfd," \n ",param1%themax-0.001_wr, &
            refseisme(1)%pfd, " "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param1%name.eq."_to") then
            call difftime(val,refseisme(1)%tps_init,temps_ref(m))
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,val,param2%themin+0.001_wr," \n ",val,param2%themax-0.001_wr, &
            " "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."_to") then
            call difftime(val,refseisme(1)%tps_init,temps_ref(m))
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,val," \n ",param1%themax-0.001_wr,val, &
            " "" | psxy $geozone $geoproj -W4,orange,- -O -K -N >> $file"
            findtest=.true.
          endif
          write(600,'(a)')"#########################################"
        endif
        ! -------                                              --------    . 2 séismes
        if (find==2) then
          findtest=.false.
          write(600,'(a)')"################## catalogue 2 ################"
          if (param1%name.eq."lon") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,refseisme(2)%lon,param2%themin+0.001_wr," \n ", &
            refseisme(2)%lon,param2%themax-0.001_wr, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."lon") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,refseisme(2)%lon," \n ", &
            param1%themax-0.001_wr,refseisme(2)%lon, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param1%name.eq."lat") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,refseisme(2)%lat,param2%themin+0.001_wr," \n ", &
            refseisme(2)%lat,param2%themax-0.001_wr, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."lat") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,refseisme(2)%lat," \n ", &
            param1%themax-0.001_wr,refseisme(2)%lat, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param1%name.eq."_zh") then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,refseisme(2)%pfd,param2%themin+0.001_wr," \n ", &
            refseisme(2)%pfd,param2%themax-0.001_wr, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."_zh")then
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,refseisme(2)%pfd," \n ", &
            param1%themax-0.001_wr,refseisme(2)%pfd, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param1%name.eq."_to") then
            call difftime(val,refseisme(2)%tps_init,temps_ref(m))
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,val,param2%themin+0.001_wr," \n ",val, &
            param2%themax-0.001_wr, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          if(param2%name.eq."_to") then
            call difftime(val,refseisme(2)%tps_init,temps_ref(m))
            write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%themin+0.001_wr,val," \n ", &
            param1%themax-0.001_wr,val, " "" | psxy $geozone $geoproj -W4,orange,-.- -O -K -N >> $file"
            findtest=.true.
          endif
          write(600,'(a)')"#########################################"
        endif
      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .

      else

      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .
      ! ------- carte pour les paramètres lon et lat           --------    .
      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .

        ! -------------------------------------------------------------    .
        ! ----- calcul de l'ellipse des 1000 meilleurs modèles --------    .
        ! -------------------------------------------------------------    . baz moyen [0:360]
        bazV(:)=0.0_wr
        open(unit=101,file="OUTPUT/GMT/ellipse-"//trim(adjustl(numberfile))//".txt")
        do i=1,1000
          write(101,*)param1%vec10000(i,1),param2%vec10000(i,1),param1%vec10000(i,2)
          call dellipsgc(param2%moy_1000,param1%moy_1000,param2%vec10000(i,1),param1%vec10000(i,1),dist,baz)
            do j=1,360
                if ((baz.gt.real(j-1,wr)).and.(baz.le.real(j,wr))) bazV(j)=bazV(j)+dist**2
            enddo
        enddo
        close(101)
        ! -------------------------------------------------------------    . lissage baz
        do i=1,100
          ! -----------------------------------------------------------    .
          bazVbis(1)=bazV(1)+(bazV(360)+bazV(2))/2.0_wr
          do j=2,359
            bazVbis(j)=bazV(j)+(bazV(j-1)+bazV(j+1))/2.0_wr
          enddo
          bazVbis(360)=bazV(360)+(bazV(359)+bazV(1))/2.0_wr
          ! -----------------------------------------------------------    .
          bazV(1)=bazVbis(1)+(bazVbis(360)+bazVbis(2))/2.0_wr
          do j=2,359
            bazV(j)=bazVbis(j)+(bazVbis(j-1)+bazVbis(j+1))/2.0_wr
          enddo
          bazV(360)=bazVbis(360)+(bazVbis(359)+bazVbis(1))/2.0_wr
          ! -----------------------------------------------------------    .
        enddo
        ! -------------------------------------------------------------    . sélection baz
        sumBAZ=-1.0_wr
        do j=1,360
          if (bazV(j).gt.sumBAZ) then
            moyBAZ=real(j,wr)-0.5_wr
            sumBAZ=bazV(j)
          endif
        enddo
        ! -------------------------------------------------------------    .
        ! -------------------------------------------------------------    . axes a et b, demi axes [km], plus o moins 1 sigma
        p_a=0.0_wr
        p_b=0.0_wr
        do i=1,1000
          call dellipsgc(param2%moy_1000,param1%moy_1000,param2%vec10000(i,1),param1%vec10000(i,1),dist,baz)
          p_c = min( abs(mod(moyBAZ-BAZ,90.0_wr)), &
                     abs(mod(mod(moyBAZ,180.0_wr)-BAZ,90.0_wr)), &
                     abs(mod(moyBAZ-mod(BAZ,180.0_wr),90.0_wr)), &
                     abs(mod(mod(moyBAZ,180.0_wr)-mod(BAZ,180.0_wr),90.0_wr)))
          p_a = p_a + (cos(p_c/180.0_wr*pi)*dist)**2.0_wr                  ! ecartype de l'axe a, moyenne nulle [km]
          p_b = p_b + (sin(p_c/180.0_wr*pi)*dist)**2.0_wr                  ! ecartype de l'axe b, moyenne nulle [km]
        enddo
        p_a=sqrt(p_a/real(1000,wr))
        p_b=sqrt(p_b/real(1000,wr))
        ! -------------------------------------------------------------    . axe a > b
        if (p_b.gt.p_a) then
          p_c=p_b
          p_b=p_a
          p_a=p_c
        endif
        ! -------------------------------------------------------------    . sauve baz, axes a et b
        E%ang=moyBAZ
        E%axeA=p_a
        E%axeB=p_b
        ! -------------------------------------------------------------    . futures gif plots
        open(unit=100,file="OUTPUT/GMT/doc-"//trim(adjustl(numberfile))//".txt")
        write(100,*) t(m)
        write(100,*) temps_ref(m)
        write(100,*) param1%moy_1000,param2%moy_1000,E%ang,E%axeA*2.0_wr,E%axeB*2.0_wr
        write(100,*) param1%moy_1000,param2%moy_1000,E%ang,E%axeA*4.0_wr,E%axeB*4.0_wr
        write(100,*) param1%moy_1000,param2%moy_1000,E%ang,E%axeA*6.0_wr,E%axeB*6.0_wr
        close(100)
        ! -------------------------------------------------------------    .
        write(600,'(a)')"geoproj=-JQ5i"
        ! ------- rééquilible min /max -> carte de au moins 5 x 5 km --    .
        ! ------- km / degree en latitude                      --------    .
        val = 2.0_wr * pi * rT * sin((90.0_wr-param2%vec10000(1,1))/180.0_wr*pi) /360.0_wr
        do while(((param1%themax-param1%themin)*val).lt.8.0_wr)
          param1%themax = param1%themax + (param1%themax-param1%themin)/100.0_wr
          param1%themin = param1%themin - (param1%themax-param1%themin)/100.0_wr
        enddo
        ! ------- km / degree en longitude                     --------    .
        val = 2.0_wr * pi * rT / 360.0_wr
        do while(((param2%themax-param2%themin)*val).lt.8.0_wr)
          param2%themax = param2%themax + (param2%themax-param2%themin)/100.0_wr
          param2%themin = param2%themin - (param2%themax-param2%themin)/100.0_wr
        enddo
        ! ------- rééquilible min /max -> map carrée           --------    .
        if ((param1%themax-param1%themin).gt.(param2%themax-param2%themin)) then
          val= (param2%themax-param2%themin)/2.0_wr
          param2%themax = param2%themax - val + (param1%themax-param1%themin)/2.0_wr
          param2%themin = param2%themin + val - (param1%themax-param1%themin)/2.0_wr
        else
          val= (param1%themax-param1%themin)/2.0_wr
          param1%themax = param1%themax - val + (param2%themax-param2%themin)/2.0_wr
          param1%themin = param1%themin + val - (param2%themax-param2%themin)/2.0_wr
        endif
        write(600,'(a10,E13.7,a1,E13.7,a1,E13.7,a1,E13.7)')"geozone=-R",param1%themin, &
        "/",param1%themax,"/",param2%themin,"/",param2%themax
        write(600,'(a)')"############## xyz -> grid #############"
        write(600,'(a19,E13.7,a1,E13.7,2a)')"xyz2grd $geozone -I",param1%delta*1.01_wr,"/",param2%delta*1.01_wr, &
        " OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".bin -bi3d -Nnan ",&
        " -GOUTPUT/GMT/topo0_"//trim(adjustl(numberfile))//".grd "
        write(600,'(a19,E13.7,a1,E13.7,2a)')"xyz2grd $geozone -I",delta_1*1.5_wr,"/",delta_2*1.5_wr, &
        " OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"mis.bin -bi3d -Nnan ", &
        " -GOUTPUT/GMT/topo2_"//trim(adjustl(numberfile))//".grd "
        ! -------------------------------------------------------------    .
        ! ------- choix des échelles et incréments             --------    .
        d_diff=(param1%themax-param1%themin)/4.0_wr
        if (d_diff.gt.0.0001_wr) diff1=real(int(d_diff*100000.0_wr,wi),wr)/100000.0_wr
        if (d_diff.gt.0.001_wr) diff1=real(int(d_diff*10000.0_wr,wi),wr)/10000.0_wr
        if (d_diff.gt.0.001_wr) diff1=real(int(d_diff*1000.0_wr,wi),wr)/1000.0_wr
        if (d_diff.gt.0.01_wr) diff1=real(int(d_diff*100.0_wr,wi),wr)/100.0_wr
        if (d_diff.gt.0.1_wr) diff1=real(int(d_diff*10.0_wr,wi),wr)/10.0_wr
        if (d_diff.gt.1.0_wr) diff1=real(int(d_diff*1.0_wr,wi),wr)/1.0_wr
        if (d_diff.gt.10.0_wr) diff1=real(int(d_diff*0.1_wr,wi),wr)/0.1_wr
        if (d_diff.gt.100.0_wr) diff1=real(int(d_diff*0.01_wr,wi),wr)/0.01_wr
        write(char_1,'(E13.7)')diff1
        ! -------                                              --------    .
        d_diff=(param2%themax-param2%themin)/4.0_wr
        if (d_diff.gt.0.0001_wr) diff2=real(int(d_diff*100000.0_wr,wi),wr)/100000.0_wr
        if (d_diff.gt.0.001_wr) diff2=real(int(d_diff*10000.0_wr,wi),wr)/10000.0_wr
        if (d_diff.gt.0.001_wr) diff2=real(int(d_diff*1000.0_wr,wi),wr)/1000.0_wr
        if (d_diff.gt.0.01_wr) diff2=real(int(d_diff*100.0_wr,wi),wr)/100.0_wr
        if (d_diff.gt.0.1_wr) diff2=real(int(d_diff*10.0_wr,wi),wr)/10.0_wr
        if (d_diff.gt.1.0_wr) diff2=real(int(d_diff*1.0_wr,wi),wr)/1.0_wr
        if (d_diff.gt.10.0_wr) diff2=real(int(d_diff*0.1_wr,wi),wr)/0.1_wr
        if (d_diff.gt.100.0_wr) diff2=real(int(d_diff*0.01_wr,wi),wr)/0.01_wr
        write(char_2,'(E13.7)')diff2
        ! -------------------------------------------------------------    .
        write(char_3,'(E13.7)')param1%delta
        write(char_4,'(E13.7)')param2%delta
        ! -------------------------------------------------------------    .
        write(600,'(a)')"#########################################"
        write(600,'(a)')"############ diag. fct coût #############"
        write(600,'(a)')"#########################################"
        write(600,'(2a)')"psbasemap $geozone $geoproj -Ba"//char_1//":"""//param1%char//""":/a"//char_2//":"""// &
        param2%char//""":WenS -K -X2.5i -Yc >  $file"
        write(600,'(a)')"pscoast $geozone $geoproj -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -O -K >>  $file"
        write(600,'(a)')"########### affiche les points ##########"
        write(600,'(2a)')"psxy $geozone $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"mis.bin -bi3d", &
        " -Sc0.015i -COUTPUT/GMT/colorpal1.cpt -O -K >> $file"
        write(600,'(a)')"############ contour densité ############"
        write(600,'(2a)')"grdcontour OUTPUT/GMT/topo0_"//trim(adjustl(numberfile))//".grd ", &
        "$geozone $geoproj -Ba0 -C75 -L74/75 -W2 -O -K >> $file"
        write(600,'(2a)')"grdcontour OUTPUT/GMT/topo0_"//trim(adjustl(numberfile))//".grd ", &
        "$geozone $geoproj -Ba0 -C50 -L49/50 -W2 -O -K >> $file"
        write(600,'(a)')"###### moy modèles des chaînes  1 #######"
        ! ------- affiche la moyennes de l'ensemble du meilleur modèle de chaque chaîne
        minmax1 = param1%moy_bestchaine - param1%ec_bestchaine
        minmax2 = param1%moy_bestchaine + param1%ec_bestchaine
        themin = param2%themin  + 0.035_wr * (param2%themax - param2%themin)
        themax = param2%themax  - 0.035_wr * (param2%themax - param2%themin)
        write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,param1%moy_bestchaine,themin," \n ",param1%moy_bestchaine,themax, &
        " "" | psxy $geozone $geoproj -W3,grey -O -K >> $file"
        minmax1 = param2%moy_bestchaine - param2%ec_bestchaine
        minmax2 = param2%moy_bestchaine + param2%ec_bestchaine
        themin = param1%themin  + 0.035_wr * (param1%themax - param1%themin)
        themax = param1%themax  - 0.035_wr * (param1%themax - param1%themin)
        write(600,'(a,2f15.5,a,2f15.5,a)')"echo -e """,themin,param2%moy_bestchaine," \n ",themax,param2%moy_bestchaine, &
        " "" | psxy $geozone $geoproj -W3,grey -O -K >> $file"
        write(600,'(a)')"######### 10000 meilleurs modèles #######"
        ! ------- affiche les 10000 meilleurs modèles          --------    .
        ! write(600,*)"psxy $geozone $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_v12.bin -Sc0.001i -O -K -bi3d >> $file"
        write(600,'(a)')"########### barre de couleur ############"
        ! ------- affiche la barre de couleur                  --------    .
        write(600,'(2a)')"psscale -D1/-1/0.50E+01/0.25ch -B25.0:""fonction co\373t"": -S -I  ", &
        "-COUTPUT/GMT/colorpal1.cpt -O -K >> $file"
        write(600,'(a)')"######### 10 meilleurs modèles ##########"
        ! ------- affiche les 10 meilleurs modèles (tous différents) sous la forme d'étoiles jaunes
        i=1
        l=1
        write(600,'(a,2f15.5,a)')"echo ",param1%vec10000(i,1),param2%vec10000(i,1), &
        "| psxy $geozone $geoproj  -Sa0.1i -Gyellow -Wthinnest,black -O -K >> $file"
        do while(l.lt.10)
          i=i+1
          if (mis%vec10000(i,1).ne.mis%vec10000(i-1,1)) then
            l=l+1
            write(600,'(a,2f15.5,a)')"echo ",param1%vec10000(i,1),param2%vec10000(i,1), &
            "| psxy $geozone $geoproj  -Sa0.1i -Gyellow -Wthinnest,black -O -K >> $file"
          endif
        enddo
        ! -------                                              --------    .
        write(600,'(a)')"psbasemap $geozone $geoproj -Ba0 -K -O >>  $file"
        write(600,'(a)')"###### moy des meilleurs modèles ########"
        ! ------- affiche la moyenne +ou- un ecart-type des 10000, 1000 et 100 meilleurs modèles, avec des flêches sur les côtés
        minmax1 = 0.015_wr * (param1%themax - param1%themin)
        minmax2 = 0.015_wr * (param2%themax - param2%themin)
        ! -------                                              --------    .
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_10000+param1%ec_10000, &
        param2%themin+minmax2,param1%moy_10000-param1%ec_10000,param2%themin+minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themin+minmax1,param2%moy_10000+param2%ec_10000, &
        param1%themin+minmax1,param2%moy_10000-param2%ec_10000, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_1000+param1%ec_1000,param2%themin+minmax2, &
        param1%moy_1000-param1%ec_1000,param2%themin+minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themin+minmax1,param2%moy_1000+param2%ec_1000, &
        param1%themin+minmax1,param2%moy_1000-param2%ec_1000, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_100+param1%ec_100,param2%themin+minmax2, &
        param1%moy_100-param1%ec_100,param2%themin+minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i  -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themin+minmax1,param2%moy_100+param2%ec_100, &
        param1%themin+minmax1,param2%moy_100-param2%ec_100, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        ! -------                                                    --------    .
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_10000+param1%ec_10000,param2%themax-minmax2, &
        param1%moy_10000-param1%ec_10000,param2%themax-minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themax-minmax1,param2%moy_10000+param2%ec_10000, &
        param1%themax-minmax1,param2%moy_10000-param2%ec_10000, &
        " "" | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_1000+param1%ec_1000,param2%themax-minmax2, &
        param1%moy_1000-param1%ec_1000,param2%themax-minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themax-minmax1,param2%moy_1000+param2%ec_1000, &
        param1%themax-minmax1,param2%moy_1000-param2%ec_1000, &
        " "" | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i  -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%moy_100+param1%ec_100,param2%themax-minmax2, &
        param1%moy_100-param1%ec_100,param2%themax-minmax2, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i  -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        write(600,'(a,4f15.5,a)')"echo -e """,param1%themax-minmax1,param2%moy_100+param2%ec_100, &
        param1%themax-minmax1,param2%moy_100-param2%ec_100, &
        " "" | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
        ! ------- affiche la moyennes de l'ensemble du meilleur modèle de chaque chaîne
        write(600,'(a)')"###### moy modèles des chaînes 2 #######"
        minmax1 = param1%moy_bestchaine - param1%ec_bestchaine
        minmax2 = param1%moy_bestchaine + param1%ec_bestchaine
        themin = param2%themin + 0.035_wr * (param2%themax - param2%themin)
        themax = param2%themax - 0.035_wr * (param2%themax - param2%themin)
        write(600,*)"echo ",minmax1,themin," 0 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",minmax2,themin," 0 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",param1%moy_bestchaine,themin," 0 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.05i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",minmax1,themax," 180 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",minmax2,themax," 180 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",param1%moy_bestchaine,themax," 180 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        minmax1 = param2%moy_bestchaine - param2%ec_bestchaine
        minmax2 = param2%moy_bestchaine + param2%ec_bestchaine
        themin = param1%themin  + 0.035_wr * (param1%themax - param1%themin)
        themax = param1%themax  - 0.035_wr * (param1%themax - param1%themin)
        write(600,*)"echo ",themin,minmax1," 90 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",themin,minmax2," 90 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",themin,param2%moy_bestchaine," 90 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.05i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",themax,minmax1," 270 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",themax,minmax2," 270 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gorange -Wthinnest,black -O -K >> $file"
        write(600,*)"echo ",themax,param2%moy_bestchaine," 270 0.1i | psxy $geozone $geoproj ", "\"
        write(600,*)"-SV0.04i/0.06i/0.04i -Gred -Wthinnest,black -O -K >> $file"
        ! -------------------------------------------------------------    .
        write(600,'(a)')"######### affiche les stations ##########"
        ok = 0
        open(605, FILE = 'DATA/sta.d',status='old',iostat = ok)
        if (ok .ne. 0) then
          write(*,*)'problème dans GMT_2paramplot : le fichier data/sta.d n''existe pas '
          stop
        endif
        do while(ok .eq. 0)
          read(605,*,iostat = ok)datasta
          if ((datasta%lon.gt.param1%themin).and.(datasta%lon.lt.param1%themax) &
                .and.(datasta%lat.gt.param2%themin).and.(datasta%lat.lt.param2%themax)) then
            write(600,*)" echo ",datasta%lon,datasta%lat," | psxy $geoproj $geozone -St0.1i -Gred  ", &
                "-Lk -Wthinnest -O -K >> $file"
          endif
        enddo
        close(605)
        ! -------------------------------------------------------------    .
        ! ------- valeurs pour placer la rose des vents et l'échelle --    .
        val = param2%themax - (param2%themax-param2%themin)*0.175_wr       ! LAT max rose
        if ((val.ge.-90.0_wr).and.(val.le.-10.0_wr)) then
          write(char_map(1),'(f10.6)')  val
        endif
        if ((val.gt.-10.0_wr).and.(val.lt.0.0_wr)) then
          write(char_map(1),'(f10.7)')  val
        endif
        if ((val.ge.0.0_wr).and.(val.lt.10.0_wr)) then
          write(char_map(1),'(f10.8)')  val
        endif
        if ((val.ge.10.0_wr).and.(val.le.90.0_wr)) then
          write(char_map(1),'(f10.7)')  val
        endif
        val = param2%themin + (param2%themax-param2%themin)*0.075_wr       ! LAT min barre
        if ((val.ge.-90.0_wr).and.(val.le.-10.0_wr)) then
          write(char_map(3),'(f10.6)')  val
        endif
        if ((val.gt.-10.0_wr).and.(val.lt.0.0_wr)) then
          write(char_map(3),'(f10.7)')  val
        endif
        if ((val.ge.0.0_wr).and.(val.lt.10.0_wr)) then
          write(char_map(3),'(f10.8)')  val
        endif
        if ((val.ge.10.0_wr).and.(val.le.90.0_wr)) then
          write(char_map(3),'(f10.7)')  val
        endif
        val = param1%themin + (param1%themax-param1%themin)*0.2_wr         ! LON min barre
        if ((val.ge.-180.0_wr).and.(val.le.-100.0_wr)) then
          write(char_map(2),'(f10.5)')  val
        endif
        if ((val.gt.-100.0_wr).and.(val.le.-10.0_wr)) then
          write(char_map(2),'(f10.6)')  val
        endif
        if ((val.gt.-10.0_wr).and.(val.lt.0.0_wr)) then
          write(char_map(2),'(f10.7)')  val
        endif
        if ((val.ge.0.0_wr).and.(val.lt.10.0_wr)) then
          write(char_map(2),'(f10.8)')  val
        endif
        if ((val.ge.10.0_wr).and.(val.lt.100.0_wr)) then
          write(char_map(2),'(f10.7)')  val
        endif
        if ((val.ge.100.0_wr).and.(val.le.180.0_wr)) then
          write(char_map(2),'(f10.6)')  val
        endif
          val = param1%themin + (param1%themax-param1%themin)*0.125_wr     ! LON min rose
        if ((val.ge.-180.0_wr).and.(val.le.-100.0_wr)) then
          write(char_map(4),'(f10.5)')  val
        endif
        if ((val.gt.-100.0_wr).and.(val.le.-10.0_wr)) then
          write(char_map(4),'(f10.6)')  val
        endif
        if ((val.gt.-10.0_wr).and.(val.lt.0.0_wr)) then
          write(char_map(4),'(f10.7)')  val
        endif
        if ((val.ge.0.0_wr).and.(val.lt.10.0_wr)) then
          write(char_map(4),'(f10.8)')  val
        endif
        if ((val.ge.10.0_wr).and.(val.lt.100.0_wr)) then
          write(char_map(4),'(f10.7)')  val
        endif
        if ((val.ge.100.0_wr).and.(val.le.180.0_wr)) then
          write(char_map(4),'(f10.6)')  val
        endif
        ! -------------------------------------------------------------    .
        ! ------- taille de l'échelle (en km)                  --------    .
        i=int((param1%themax-param1%themin)/8.0_wr/360._wr*2._wr*pi*rT)
        write(char_6,'(i2.2)')  i
        write(600,'(a)')"######### rose vents et échelle #########"
        write(600,*) "psbasemap $geozone $geoproj -Ba0g0.06", &
        " -Lf"//char_map(2)//"/"//char_map(3)//"/"//char_map(3)//"/"//char_6//"k+l+jl", &
        " -Tf"//char_map(4)//"/"//char_map(1)//"/0.75i/3:@,@,@,@-N@-:", &
        " -O -K >>  $file"
        write(600,'(a)')"#########################################"
        write(600,'(a)')"############# diag. densité #############"
        write(600,'(a)')"#########################################"
        write(600,*)"psbasemap $geozone $geoproj -Ba"//char_1//":"""//param1%char//""":/a"//char_2//":"""// &
        param2%char//""":WenS -K -O -X6.25i >>  $file"
        write(600,*)"pscoast $geozone $geoproj -S240/255/255 -G180/238/180 -N1 -Df+ -Ia/blue -W1 -O -K >>  $file"
        write(600,'(a)')"############### grid image ##############"
        write(600,*)"grdimage $geozone $geoproj OUTPUT/GMT/topo0_"//trim(adjustl(numberfile))//".grd  ", &
        "-QNaN -COUTPUT/GMT/colorpal2.cpt -B0 -O -K -Sn >> $file"
        write(600,'(a)')"############## grid contour #############"
        write(600,*)"grdcontour OUTPUT/GMT/topo0_"//trim(adjustl(numberfile))//".grd $geozone $geoproj  ", &
        "-Ba0 -C75 -L74/75 -W2 -O -K >> $file"
        write(600,*)"grdcontour OUTPUT/GMT/topo0_"//trim(adjustl(numberfile))//".grd $geozone $geoproj  ", &
        "-Ba0 -C50 -L49/50 -W2 -O -K >> $file"
        write(600,'(a)')"######### cercles concentriques #########"
        do i=2,100,2
          write(600,*)"echo "" ",param1%vec10000(1,1),param2%vec10000(1,1),"0",i,i," "" \"
          write(600,*)" | psxy $geozone $geoproj -SE -W3,grey,-- -O -K >> $file"
        enddo
        ! -------------------------------------------------------------    . ellipses
        write(numberfile(1:5),'(i5)')m
        write(600,*)"echo ",param1%moy_1000,param2%moy_1000,E%ang,E%axeA*2.0_wr,E%axeB*2.0_wr, &
            " > OUTPUT/GMT/ellipse-"//trim(adjustl(numberfile))//".txt"
        write(600,*)"echo ",param1%moy_1000,param2%moy_1000,E%ang,E%axeA*4.0_wr,E%axeB*4.0_wr, &
            " >> OUTPUT/GMT/ellipse-"//trim(adjustl(numberfile))//".txt"
        write(600,*)"echo ",param1%moy_1000,param2%moy_1000,E%ang,E%axeA*6.0_wr,E%axeB*6.0_wr, &
            " >> OUTPUT/GMT/ellipse-"//trim(adjustl(numberfile))//".txt"
        write(600,*)"psxy $geozone $geoproj -SE -W2,red,-- OUTPUT/GMT/ellipse-"//trim(adjustl(numberfile))//".txt -O -K -N >> $file"
        ! -------------------------------------------------------------    .
        write(600,'(a)')"######### rose vents et échelle #########"
        write(600,*) "psbasemap $geozone $geoproj -Ba0g0.06", &
        " -Lf"//char_map(2)//"/"//char_map(3)//"/"//char_map(3)//"/"//char_6//"k+l+jl", &
        " -Tf"//char_map(4)//"/"//char_map(1)//"/0.75i/3:@,@,@,@-N@-:", &
        " -O -K >>  $file"
        write(600,*)"psscale -D1/-1/0.50E+01/0.25ch -B25.0:""densit\351"": -S -I -COUTPUT/GMT/colorpal2.cpt -O -K >> $file"
        write(600,'(a)')"################## catalogue 1 ##################"
        call catalogue(param_best,refseisme,find)
        if ((find==1).or.(find==2)) then
          findtest=.true.
          write(600,*)"echo -e """,refseisme(1)%lon,param2%themin+0.001_wr," \n ",refseisme(1)%lon,param2%themax-0.001_wr," \"
          write(600,*)" "" | psxy $geozone $geoproj -W4,orange,- -O -K >> $file"
          write(600,*)"echo -e """,param1%themin+0.001_wr,refseisme(1)%lat," \n ",param1%themax-0.001_wr,refseisme(1)%lat," \"
          write(600,*)" "" | psxy $geozone $geoproj -W4,orange,- -O -K >> $file"
          write(600,'(a)')"#########################################"
        endif
        if (find==2) then
          findtest=.true.
          write(600,*)"echo -e """,refseisme(2)%lon,param2%themin+0.001_wr," \n ",refseisme(2)%lon,param2%themax-0.001_wr," \"
          write(600,*)" "" | psxy $geozone $geoproj -W4,orange,-.- -O -K >> $file"
          write(600,*)"echo -e """,param1%themin+0.001_wr,refseisme(2)%lat," \n ",param1%themax-0.001_wr,refseisme(2)%lat," \"
          write(600,*)" "" | psxy $geozone $geoproj -W4,orange,-.- -O -K >> $file"
          write(600,'(a)')"#########################################"
        endif

        ! -------------------------------------------------------------    . geiger avec Modèle Terre Arroucau

        write(600,'(2a)')"psxy $geoproj $geozone -Sa0.15i -Gred -Wthinnest ", & ! étoile rouge
          " OUTPUT/GMT/ArrALL-"//trim(adjustl(numberfile))//".txt -O -K -N >> $file"
        write(600,'(2a)')"head -n 1 OUTPUT/GMT/ArrALLt-"//trim(adjustl(numberfile))//".txt ", &
          " | pstext $geoproj $geozone -D.5/.5v1 -O -K -N >> $file"
        write(600,'(2a)')"tail -n 1 OUTPUT/GMT/ArrALLt-"//trim(adjustl(numberfile))//".txt ", &
          " | pstext $geoproj $geozone -D.5/.5v1 -O -K -N >> $file"

        ! -------------------------------------------------------------    . geiger avec Modèle Terre Si-Hex

        write(600,'(2a)')"psxy $geoproj $geozone -Sa0.15i -Gblue -Wthinnest ", & ! étoile bleue
          "  OUTPUT/GMT/SiHexALL-"//trim(adjustl(numberfile))//".txt -O -K -N >> $file"
        write(600,'(2a)')"head -n 1 OUTPUT/GMT/SiHexALLt-"//trim(adjustl(numberfile))//".txt ", &
          " | pstext $geoproj $geozone -D.5/.5v1 -O -K -N >> $file"
        write(600,'(2a)')"tail -n 1 OUTPUT/GMT/SiHexALLt-"//trim(adjustl(numberfile))//".txt ", &
          " | pstext $geoproj $geozone -D.5/.5v1 -O -K -N >> $file"

        ! -------------------------------------------------------------    . geiger avec Modèle Terre CÉA

        write(600,'(2a)')"psxy $geoproj $geozone -Sa0.15i -Ggreen -Wthinnest ", & ! étoile green
          "  OUTPUT/GMT/CEAALL-"//trim(adjustl(numberfile))//".txt -O -K -N >> $file"
        write(600,'(2a)')"head -n 1 OUTPUT/GMT/CEAALLt-"//trim(adjustl(numberfile))//".txt ", &
          " | pstext $geoproj $geozone -D.5/.5v1 -O -K -N >> $file"
        write(600,'(2a)')"tail -n 1 OUTPUT/GMT/CEAALLt-"//trim(adjustl(numberfile))//".txt ", &
          " | pstext $geoproj $geozone -D.5/.5v1 -O -K -N >> $file"

        ! -------------------------------------------------------------    .
        write(600,'(a)')"############## légende n0 ###############"
        val = param2%themax-(param2%themax-param2%themin)*0.035_wr
        val1 = param1%themin+(param1%themax-param1%themin)*0.275_wr
        val2 = param1%themin+(param1%themax-param1%themin)*0.325_wr
        write(600,*)"echo  -e '",val1,val," \n ",val2,val,"' | psxy $geozone $geoproj -W2,red,-- -O -K >> $file "
        val1 = param1%themin+(param1%themax-param1%themin)*0.35_wr
        write(600,*)"echo '",val1,val," 12 0 4 LM ellipse @:10:\050\2611, 2 et 3@~\163@~\051@:: ", &
        "des 1000 meilleurs mod\350les ' | pstext $geozone $geoproj -Gred -O -K >> $file"
        ! -------------------------------------------------------------    .
        ! plots                                                       -    .
        ! -------------------------------------------------------------    .

      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .

      endif ! (param lon-lat ou autres ?)

      ! ---------------------------------------------------------------    .
      ! ---------------------------------------------------------------    .

      write(600,'(a)')"#########################################"
      write(600,'(a)')"################ histo 1 ################"
      write(600,'(a)')"#########################################"
      write(numberfile(1:5),'(i5)')m
      write(600,*)"geoproj=-JX5i/1.5i"
      write(600,'(a10,E13.7,a1,E13.7,a5)')"geozone=-R",param1%themin,"/",param1%themax,"/0/15"
      write(600,*)"psbasemap $geozone $geoproj -Ba"//char_1//"/a10:""probabilit\351 (%)"":sWen -K -O -Y5.25i  >>  $file"
      write(600,'(a)')"################# mode #################"
      write(600,*)"echo -e """,param1%mode,"0 \n",param1%mode," 100 "" | psxy $geozone $geoproj -W2,red -O -K >> $file"
      write(600,*)"echo -e """,param1%mediane,"0 \n",param1%mediane," 100 "" | psxy $geozone $geoproj -W2,blue -O -K >> $file"
      write(600,*)"echo -e """,param1%moy_tot,"0 \n",param1%moy_tot," 100 "" | psxy $geozone $geoproj -Wgreen -O -K >> $file"
      write(600,*)"echo -e """,param1%moy_tot+param1%ec_tot,"0 \n",param1%moy_tot+param1%ec_tot, &
      " 100 "" | psxy $geozone $geoproj -Wgreen,:-: -O -K >> $file"
      write(600,*)"echo -e """,param1%moy_tot-param1%ec_tot,"0 \n",param1%moy_tot-param1%ec_tot, &
      " 100 "" | psxy $geozone $geoproj -Wgreen,:-: -O -K >> $file"
      write(600,'(a)')"############ chaque chaîne 1 ############"
      ! ------- distribution de probabilités marginales a posteriori pour chaque chaîne
      k=0
      do i=1,nbChaineMV                                                    ! fond gris, puis barres noires
        write(char_5,'(i5)')10000+i
        open(unit=650+i,file="OUTPUT/GMT/"//filename//char_5//"-"//trim(adjustl(numberfile))//".bin" &
        ,STATUS="replace",access='direct',RECL=16)
        do j=1,nmod(i)
          k=k+1
          write(650+i,rec=j)real(param1%vec(k),8),real(param2%vec(k),8)
        enddo
        close(650+i)
        if (nmod(i).gt.100) then
          write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/"//filename//char_5//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -T0 -W"//char_3//" -Ggray -Ba0/a0 -O -K -Z1 >>  $file"
        endif
      enddo
      write(600,'(a)')"################# geiger ################"
      existe3=.true.
      if(param1%name.eq."_vc") existe3=.false.
      if(param1%name.eq."_vm") existe3=.false.
      if(param1%name.eq."_zm") existe3=.false.
      if(param1%name.eq."vps") existe3=.false.
      if (FLAGhypofixe) then
        if(param1%name.eq."lon") existe3=.false.
        if(param1%name.eq."lat") existe3=.false.
        if(param1%name.eq."_zh") existe3=.false.
        if(param1%name.eq."_to") existe3=.false.
      endif
      if(existe3) then
        ! ------- distribution de probabilités marginales a priori ------    .
        write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/geiger"//param1%name//"-"//trim(adjustl(numberfile))//".bin", &
        " -bi1d -W"//char_3//" -Gred -Ba0/a0 -O -K -Z1 >>  $file"
        write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/geiger"//param1%name//"-"//trim(adjustl(numberfile))//".bin", &
        " -bi1d -W"//char_3//" -S -L0/0 -Ba0/a0 -O -K -Z1 >>  $file"
      endif
      write(600,'(a)')"################# apriori ################"
      existe1=.true.
      if (FLAGterrefixe) then
        if(param1%name.eq."_vc") existe1=.false.
        if(param1%name.eq."_vm") existe1=.false.
        if(param1%name.eq."_zm") existe1=.false.
        if(param1%name.eq."vps") existe1=.false.
      endif
      if (FLAGhypofixe) then
        if(param1%name.eq."lon") existe1=.false.
        if(param1%name.eq."lat") existe1=.false.
        if(param1%name.eq."_zh") existe1=.false.
        if(param1%name.eq."_to") existe1=.false.
      endif
      if(existe1) then
      ! ------- distribution de probabilités marginales a priori ------    .
        if((param1%name.eq."_vc").or.(param1%name.eq."_vm").or.(param1%name.eq."_zm").or.(param1%name.eq."vps"))then
          write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/aprio"//param1%name//"-1.bin", &
            " -bi1d -W"//char_3//" -Gorange -Ba0/a0 -O -K -Z1 >>  $file"
          write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/aprio"//param1%name//"-1.bin", &
            " -bi1d -W"//char_3//" -S -L0/0 -Ba0/a0 -O -K -Z1 >>  $file"
        else
          write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/aprio"//param1%name//"-"//trim(adjustl(numberfile))//".bin", &
            " -bi1d -W"//char_3//" -Gorange -Ba0/a0 -O -K -Z1 >>  $file"
          write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/aprio"//param1%name//"-"//trim(adjustl(numberfile))//".bin", &
            " -bi1d -W"//char_3//" -S -L0/0 -Ba0/a0 -O -K -Z1 >>  $file"
        endif
      endif
      write(600,'(a)')"################ totale #################"
      ! ------- distribution de probabilités marginales a posteriori TOTALE  !!! -Gp300/1:BblueF-
      write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_tot.bin", &
      " -bi3d -T0  -W"//char_3//" -Gblue -Ba0/a0 -O -K -Z1 >>  $file" !
      write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_tot.bin", &
      " -bi3d -T0  -W"//char_3//" -L1/0 -S -Ba0/a0 -O -K -Z1 >>  $file" !
      write(600,'(a)')"############ chaque chaîne 2 ############"
      do i=1,nbChaineMV
        write(char_5,'(i5)')10000+i
        if (nmod(i).gt.100) then
          write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/"//filename//char_5//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -T0  -W"//char_3//" -L0/0 -Ba0/a0 -O -K -Z1 -S >>  $file"
        endif
      enddo

      if(existe3) then
        ! ------- distribution de probabilités marginales a priori ----    .
        write(600,*)"pshistogram $geozone $geoproj OUTPUT/GMT/geiger"//param1%name//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi1d -W"//char_3//" -Sred -L1/1 -Ba0/a0 -O -K -Z1 >>  $file"
      endif

      write(600,*)"psbasemap $geozone $geoproj -Ba0 -K -O >>  $file"
      ! ---------------------------------------------------------------    .
      write(600,'(a)')"#########################################"
      write(600,'(a)')"################ histo 2 ################"
      write(600,'(a)')"#########################################"
      write(600,*)"geoproj=-JX1.5i/5i"
      write(600,'(a15,E13.7,a1,E13.7)')"geozone=-R0/15/",param2%themin,"/",param2%themax
      write(600,'(a13,E13.7,a1,E13.7,a5)')"geozonebis=-R",param2%themin,"/",param2%themax,"/0/15"
      write(600,*)"psbasemap $geozone $geoproj -Ba10:""probabilit\351 (%)"":/a"//char_2//"Swne -K -O -Y-5.25i -X5.25i  >>  $file"
      write(600,'(a)')"################# mode #################"
      write(600,*)"echo -e "" 0",param2%mode,"\n 100 ",param2%mode,""" | psxy $geozone $geoproj -W2,red -O -K >> $file"
      write(600,*)"echo -e "" 0",param2%mediane,"\n 100 ",param2%mediane,""" | psxy $geozone $geoproj -W2,blue -O -K >> $file"
      write(600,*)"echo -e "" 0",param2%moy_tot," \n 100 ",param2%moy_tot,""" | psxy $geozone $geoproj -Wgreen -O -K >> $file"
      write(600,*)"echo -e "" 0",param2%moy_tot+param2%ec_tot," \n 100 ",param2%moy_tot+param2%ec_tot, &
      """ | psxy $geozone $geoproj -Wgreen,:-: -O -K >> $file"
      write(600,*)"echo -e "" 0",param2%moy_tot-param2%ec_tot," \n 100 ",param2%moy_tot-param2%ec_tot, &
      """ | psxy $geozone $geoproj -Wgreen,:-: -O -K >> $file"
      write(600,'(a)')"############ chaque chaîne 1 ############"
      ! ------- distribution de probabilités marginales a posteriori pour chaque chaîne
      do i=1,nbChaineMV                                                    ! fond gris, puis barres noires
        write(char_5,'(i5)')10000+i
        if (nmod(i).gt.100) then
          write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/"//filename//char_5//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -T1  -W"//char_4//" -Ggray -Ba0/a0 -K -O -Z1 -A0 >>  $file" ! -bi3d
        endif
      enddo
      write(600,'(a)')"################# geiger ################"
      existe4=.true.
      if(param2%name.eq."_vc") existe4=.false.
      if(param2%name.eq."_vm") existe4=.false.
      if(param2%name.eq."_zm") existe4=.false.
      if(param2%name.eq."vps") existe4=.false.
      if (FLAGhypofixe) then
        if(param2%name.eq."lon") existe4=.false.
        if(param2%name.eq."lat") existe4=.false.
        if(param2%name.eq."_zh") existe4=.false.
        if(param2%name.eq."_to") existe4=.false.
      endif
      if(existe4) then
        ! ------- distribution de probabilités marginales a priori ------    .
        write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/geiger"//param2%name//"-"//trim(adjustl(numberfile))//".bin", &
        " -bi1d -W"//char_4//" -Gred -Ba0/a0 -O -K -Z1 -A0 >>  $file"
        write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/geiger"//param2%name//"-"//trim(adjustl(numberfile))//".bin", &
        " -bi1d -W"//char_4//" -S -L0/0 -Ba0/a0 -O -K -Z1 -A0 >>  $file"
      endif
      write(600,'(a)')"################# apriori ################"
      existe2=.true.
      if (FLAGterrefixe) then
        if(param2%name.eq."_vc") existe2=.false.
        if(param2%name.eq."_vm") existe2=.false.
        if(param2%name.eq."_zm") existe2=.false.
        if(param2%name.eq."vps") existe2=.false.
      endif
        if (FLAGhypofixe) then
        if(param2%name.eq."lon") existe2=.false.
        if(param2%name.eq."lat") existe2=.false.
        if(param2%name.eq."_zh") existe2=.false.
        if(param2%name.eq."_to") existe2=.false.
      endif
      if(existe2) then
        ! ------- distribution de probabilités marginales a priori ----    .
        if((param2%name.eq."_vc").or.(param2%name.eq."_vm").or.(param2%name.eq."_zm").or.(param2%name.eq."vps"))then
          write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/aprio"//param2%name//"-1.bin", &
            " -bi1d -W"//char_4//" -Gorange -Ba0/a0 -O -K -Z1 -A0 >>  $file"
          write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/aprio"//param2%name//"-1.bin", &
            " -bi1d -W"//char_4//" -S -L0/0 -Ba0/a0 -O -K -Z1 -A0 >>  $file"
        else
          write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/aprio"//param2%name//"-"//trim(adjustl(numberfile))//".bin", &
            " -bi1d -W"//char_4//" -Gorange -Ba0/a0 -O -K -Z1 -A0 >>  $file"
          write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/aprio"//param2%name//"-"//trim(adjustl(numberfile))//".bin", &
            " -bi1d -W"//char_4//" -S -L0/0 -Ba0/a0 -O -K -Z1 -A0 >>  $file"
        endif
      endif
      write(600,'(a)')"################ totale #################"
      ! ------- distribution de probabilités marginales a posteriori TOTALE  !!! -Gp300/1:BblueF-
      write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_tot.bin -T1", &
      " -bi3d -W"//char_4//" -Gblue -Ba0/a0 -O -K -Z1 -A0 >>  $file"
      write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"_tot.bin -T1", &
      " -bi3d -W"//char_4//" -L1/0 -S -Ba0/a0 -O -K -Z1 -A0 >>  $file"
      write(600,'(a)')"############ chaque chaîne 2 ############"
      do i=1,nbChaineMV
        write(char_5,'(i5)')10000+i
        if (nmod(i).gt.100) then
          write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/"//filename//char_5//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi2d -T1  -W"//char_4//" -L0/0 -Ba0/a0 -K -O -Z1 -A0 -S >>  $file" ! -bi3d
        endif
      enddo
      if(existe4) then
        ! ------- distribution de probabilités marginales a priori ----    .
        write(600,*)"pshistogram $geozonebis $geoproj OUTPUT/GMT/geiger"//param2%name//"-"//trim(adjustl(numberfile))//".bin", &
          " -bi1d -W"//char_4//" -Sred -L1/1 -Ba0/a0 -O -K -Z1 -A0 >>  $file"
      endif

      write(600,*)"psbasemap $geozone $geoproj -Ba0 -O -K >>  $file"
      ! ---------------------------------------------------------------    .
      write(600,'(a)')"#########################################"
      write(600,'(a)')"############## légende n1 ###############"
      write(600,'(a)')"#########################################"
      write(600,*)"geoproj=-JX1.5i/1.5i"
      write(600,*)"geozone=-R0/1/0/1"
      write(600,*)"echo -e '0.0 .10 \n 0.125 .10' | psxy $geozone $geoproj -Wgreen -O -K -Y5.25i >> $file"
      write(600,*)"echo '0.15 .10 15 0 4 LM moyenne @:10: \050\2611@~\163@~\051 @::' | pstext $geozone $geoproj -O -K >> $file"
      write(600,*)"echo -e '0.0 .12 \n 0.125 .12' | psxy $geozone $geoproj -Wgreen,:-: -O -K >> $file"
      write(600,*)"echo -e '0.0 .08 \n 0.125 .08' | psxy $geozone $geoproj -Wgreen,:-: -O -K >> $file"
      write(600,*)"echo -e '0.0 .25 \n 0.125 .25' | psxy $geozone $geoproj -W2,red  -O -K  >> $file"
      write(600,*)"echo -e '0.5 .25 \n 0.625 .25' | psxy $geozone $geoproj -W2,blue  -O -K  >> $file"
      write(600,*)"echo '0.15 .25 15 0 4 LM mode' | pstext $geozone $geoproj -O -K >> $file"
      write(600,*)"echo '0.65 .25 15 0 4 LM m\351diane' | pstext $geozone $geoproj -O -K -N >> $file"

       if (findtest) then
          write(600,*)"echo -e '0.01 .4 \n 0.125 .4' | psxy $geozone $geoproj -W4,orange,- -O -K  >> $file"
          write(600,*)"echo '0.15 .4 15 0 4 LM catalogue' | pstext $geozone $geoproj -O -K >> $file"
      endif

      write(600,*)"echo '0.02 .875 0.35 0.25' | psxy $geozone $geoproj -Ggray -Sr -O -K -L -N >> $file"
      write(600,*)"echo '0.02 .875 0.35 0.25' | psxy $geozone $geoproj -W1,black -Sr -O -K -L -N >> $file"
      write(600,*)"echo '0.02 .735 0.35 0.25' | psxy $geozone $geoproj -Gblue -Sr -O -K -L -N >> $file"
      write(600,*)"echo '0.02 .735 0.35 0.25' | psxy $geozone $geoproj -W1,black -Sr -O -K -L -N >> $file"

      if (existe1.or.existe2) write(600,*)"echo '0.02 .58 0.35 0.25' | psxy $geozone $geoproj -Gorange -Sr -O -K -L -N >> $file"
      if (existe1.or.existe2) write(600,*)"echo '0.02 .58 0.35 0.25' | psxy $geozone $geoproj -W1,black -Sr -O -K -L -N >> $file"
      if (existe3.or.existe4) write(600,*)"echo '0.635 .58 0.35 0.25' | psxy $geozone $geoproj -Gred -Sr -O -K -L -N >> $file"
      if (existe3.or.existe4) write(600,*)"echo '0.635 .58 0.35 0.25' | psxy $geozone $geoproj -W1,black -Sr -O -K -L -N >> $file"
      if (existe1.or.existe2) write(600,*)"echo '0.15 .60 15 0 4 LM init.' | pstext $geozone $geoproj -O -K >> $file"
      if (existe3.or.existe4) write(600,*)"echo '0.70 .60 15 0 4 LM Geiger' | pstext $geozone $geoproj -O -K -N >> $file"

      write(600,*)"echo '0.15 .90 15 0 4 LM \050par cha\356ne\051' | pstext $geozone $geoproj -O -K >> $file"
      write(600,*)"echo '0.15 .75 15 0 4 LM a posteriori' | pstext $geozone $geoproj -O -K >> $file"
      ! ---------------------------------------------------------------    .
      write(600,'(a)')"#########################################"
      write(600,'(a)')"############## légende n2 ###############"
      write(600,'(a)')"#########################################"
      write(600,*)"geoproj=-JX5i/1.5i"
      write(600,*)"geozone=-R0/1/0/1"
      write(600,'(a)')"################ figure #################"
      write(600,*)"echo '0.075 .05 0 0.1i' | psxy $geozone $geoproj", &
        " -SV0.04i/0.05i/0.04i -Gred -Wthinnest,black  -O -K -X-11.5i -UTL/0/1.5i  >> $file"
      write(600,*)"echo '0.100 .05 0 0.1i' | psxy $geozone $geoproj", &
        " -SV0.04i/0.05i/0.04i -Gorange -Wthinnest,black  -O -K >> $file"
      write(600,*)"echo '0.050 .05 0 0.1i' | psxy $geozone $geoproj", &
        " -SV0.04i/0.05i/0.04i -Gorange -Wthinnest,black  -O -K >> $file"
      write(600,*)"echo '0.01 .25 0.15 .25' | psxy $geozone $geoproj -SVS0.08i/0.12i/0.08i", &
        " -Gp300/73:Bgray45F- -Wthinnest,black -O -K >> $file"
      write(600,*)"echo '0.01 .40 0.15 .40' | psxy $geozone $geoproj -SVS0.06i/0.10i/0.06i", &
        " -Gp300/73:Bgray80F- -Wthinnest,black -O -K >> $file"
      write(600,*)"echo '0.01 .55 0.15 .55' | psxy $geozone $geoproj -SVS0.04i/0.08i/0.04i", &
        " -Gp300/73:Bgray100F- -Wthinnest,black -O -K >> $file"
      write(600,*)"echo '0.075 .7' | psxy $geozone $geoproj -Sa0.1i -Gyellow -Wthinnest,black -O -K >> $file"
      write(600,'(a)')"################# texte #################"
      write(600,*)"echo '0.175 .10 15 0 4 LM moyenne des meilleurs mod\350les de chaque cha\356ne'", &
        " | pstext $geozone $geoproj -O -K >> $file"
      write(600,*)"echo '0.175 .25 15 0 4 LM moyenne @:10: \050\2611@~\163@~\051 @:: des 10000 meilleurs mod\350les'", &
        " | pstext $geozone $geoproj -O -K >> $file"
      write(600,*)"echo '0.175 .40 15 0 4 LM moyenne @:10: \050\2611@~\163@~\051 @:: des 1000 meilleurs mod\350les'", &
        " | pstext $geozone $geoproj -O -K >> $file"
      write(600,*)"echo '0.175 .55 15 0 4 LM moyenne @:10: \050\2611@~\163@~\051 @:: des 100 meilleurs mod\350les'", &
        " | pstext $geozone $geoproj -O -K >> $file"
      write(600,*)"echo '0.174 .70  15 0 4 LM 10 meilleurs mod\350les diff\351rents' | pstext $geozone $geoproj -O >> $file"
      ! ---------------------------------------------------------------    .
      write(600,'(a)')"#########################################"
      write(600,'(a)')"############# conversions ###############"
      write(600,'(a)')"#########################################"
      write(600,*)"ps2raster OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".ps -Tf -A"
      write(600,'(2a)')"mv OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".pdf ", &
        "OUTPUT/figures/"//filename//"-"//trim(adjustl(numberfile))//".pdf"
      write(600,'(a)')"#########################################"
      write(600,*)
      write(600,*)"#***************************************#"
      write(600,*)"#***************************************#"
      write(600,*)"ELAPSED=$(($SECONDS-$BEFORE))"
      write(600,*)" echo $ELAPSED secondes"
      call system_clock(Nnewtime,ratetime)
      tl=(real(Nnewtime,wr)-real(Noldtime,wr))/real(ratetime,wr)
      write(*,'(a9,i2.2,'':'',i2.2,'':'',f9.2)')' temps : ',int(tl/3600.0_wr,wi), &
      int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi),(tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr)
      ! ---------------------------------------------------------------    .

    end subroutine GMT_2paramplot
    ! -----------------------------------------------------------------    .
  end subroutine GMTfull

  ! -------------------------------------------------------------------    .

  subroutine affiche_temps_GMT(temps_ref,char_1)
    ! -------                                                  --------    .mh
    ! conversion du temps de référence -> GMT, pas utilisée mais peut être utile
    ! -----------------------------------------------------------------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    character(LEN=30), intent (out) :: char_1
    type(date_sec), intent (in) :: temps_ref
    ! -----------------------------------------------------------------    .
    write(char_1,"(i4.4,4(1a,i2.2))") temps_ref%date%year,"-",temps_ref%date%month, &
    "-",temps_ref%date%day,"T",temps_ref%date%hour,":",temps_ref%date%min
    ! -----------------------------------------------------------------    .
  end subroutine affiche_temps_GMT

  ! -------------------------------------------------------------------    .

  subroutine affiche_temps_ref(temps_ref,char_1,ok)
    ! -------                                                  --------    .mh
    ! affiche le temps de référence sur les figures
    ! -------                                                  --------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    character(LEN=30), intent (out) :: char_1
    type(date_sec), intent (in) :: temps_ref
    integer(KIND=wi), intent (in)  :: ok
    ! -------                                                  --------    .
    character(LEN=30) :: lemois
    character(LEN=2) :: char_2
    ! -----------------------------------------------------------------    .
    if (ok == 1) then                                                      ! style LaTeX pour les accents
      if(temps_ref%date%month.eq.1) lemois="janvier"
      if(temps_ref%date%month.eq.2) lemois="f\'evrier"
      if(temps_ref%date%month.eq.3) lemois="mars"
      if(temps_ref%date%month.eq.4) lemois="avril"
      if(temps_ref%date%month.eq.5) lemois="mai"
      if(temps_ref%date%month.eq.6) lemois="juin"
      if(temps_ref%date%month.eq.7) lemois="juillet"
      if(temps_ref%date%month.eq.8) lemois="ao\^ut"
      if(temps_ref%date%month.eq.9) lemois="septembre"
      if(temps_ref%date%month.eq.10) lemois="octobre"
      if(temps_ref%date%month.eq.11) lemois="novembre"
      if(temps_ref%date%month.eq.12) lemois="d\'ecembre"
    else                                                                   ! style GMT pour les accents
      if(temps_ref%date%month.eq.1) lemois="janvier"
      if(temps_ref%date%month.eq.2) lemois="f\351vrier"
      if(temps_ref%date%month.eq.3) lemois="mars"
      if(temps_ref%date%month.eq.4) lemois="avril"
      if(temps_ref%date%month.eq.5) lemois="mai"
      if(temps_ref%date%month.eq.6) lemois="juin"
      if(temps_ref%date%month.eq.7) lemois="juillet"
      if(temps_ref%date%month.eq.8) lemois="ao\373t"
      if(temps_ref%date%month.eq.9) lemois="septembre"
      if(temps_ref%date%month.eq.10) lemois="octobre"
      if(temps_ref%date%month.eq.11) lemois="novembre"
      if(temps_ref%date%month.eq.12) lemois="d\351cembre"
    endif
    ! -------                                                  --------    .
    if ((temps_ref%date%day.gt.0).and.(temps_ref%date%day.lt.32)) then
      write(char_2,'(i2.2)')len(trim(lemois))
      write(char_1,"(i2.2,1x,a"//char_2//",1x,i4.4,1x,i2.2,a1,i2.2)")temps_ref%date%day,lemois,temps_ref%date%year, &
      temps_ref%date%hour,":",temps_ref%date%min
    else
      write(char_1,'(a30)')"problème in affiche_temps_ref"
    endif
    ! -----------------------------------------------------------------    .
  end subroutine affiche_temps_ref

  ! -------------------------------------------------------------------    .

  subroutine mkdensityplot(l,param1,param2,deltaxy,nbparam,filename)
    ! -------                                                  --------    .mh
    ! Calcul du diagramme de densité pour deux paramètres (param1,param2) et
    ! stock une grille xyz en .bin de deltaxy x deltaxy arguments
    ! -------                                                  --------    .
    use typetemps
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent (in) :: l
    type(densityplot_one), intent (inout) :: param1, param2                ! les deux paramètres
    integer(KIND=wi), intent (in) :: deltaxy                               ! pas de discrétisation pour le mode et le diagramme de densité
    integer(KIND=wi), intent (in) :: nbparam                               ! nombre de modèles
    character(LEN=7), intent (in) :: filename                              ! base du nom des fichiers de sorties
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k
    ! -------                                                  --------    .
    real(KIND=wr) :: themax
    real(KIND=wr), dimension(:,:), allocatable :: tab                      ! grille pour le diagramme de densité
    character (LEN=5) :: numberfile
    ! -----------------------------------------------------------------    .
    write(numberfile(1:5),'(i5)')l
    allocate(tab(deltaxy,deltaxy))
    do i=1,deltaxy
      do j=1,deltaxy
        tab(i,j)=0.0_wr                                                    ! initialisation
      enddo
    enddo
    ! ------- delta de pas                                     --------    .
    param1%delta=(param1%themax-param1%themin)/real(deltaxy,wr)
    param2%delta=(param2%themax-param2%themin)/real(deltaxy,wr)
    ! -----------------------------------------------------------------    . calcul de densité
    do k=1,nbparam
      i=int((param1%vec(k)-param1%themin)/param1%delta)+1
      j=int((param2%vec(k)-param2%themin)/param2%delta)+1
      if ((i.ge.1).and.(i.le.deltaxy)) then
        if ((j.ge.1).and.(j.le.deltaxy)) then
          tab(i,j)=tab(i,j)+1.0_wr
        else
          write(*,*)'problème dans mkdensityplot : calcul de densité, incrément incorrect, boucle 1',i,i,k
        endif
      else
        write(*,*)'problème dans mkdensityplot : calcul de densité, incrément incorrect, boucle 2',i,i,k
      endif
    enddo
    ! ------- normalisation (min max entre 0 et 100%)          --------    .
    themax=-1.0_wr
    do i=1,deltaxy
      do j=1,deltaxy
        if (themax.lt.tab(i,j)) themax=tab(i,j)                            ! garde le maximum
      enddo
    enddo
    do i=1,deltaxy
      do j=1,deltaxy
        tab(i,j)= 100.0_wr/(themax)*tab(i,j)
      enddo
    enddo
    ! -----------------------------------------------------------------    .
    ! ------- ecriture dans un fichier de la grille de densité --------    .
    open(unit=606,file="OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//".bin",STATUS="replace",access='direct',RECL=24)
    k=1
    do i=1,deltaxy
      do j=1,deltaxy
        if (tab(i,j).gt.0.5_wr) then                                       ! au moins 0.5 %
          write(606,rec=k)real(param1%themin+(real(i,wr)-.5_wr)*param1%delta,8) , &
          real(param2%themin +(real(j,wr)-.5_wr)*param2%delta,8),real(tab(i,j),8)
          k=k+1
        endif
      enddo
    enddo
    close(606)
    ! -----------------------------------------------------------------    .
    deallocate(tab)
    ! -----------------------------------------------------------------    .
  end subroutine mkdensityplot

  ! -------------------------------------------------------------------    .

  subroutine mkfcoutplot(l,param1,param2,mis,deltaxymis,nbparam,filename,delta_1,delta_2)
    ! -------                                                  --------    .mh
    ! Identifie les points représentatifs pour la représentation de la fonction coût
    ! Un grille échantillonant les plus petits misfits pour les deux paramètres (param1,param2)
    ! puis stock les points triés en xyz en .bin
    ! -> permet de faire des figures avec moins de recouvrement de point (donc plus légères)
    ! -> permet le tri des modèles selon la fonction coût
    ! -------                                                  --------    .
    use typetemps
    use tri
    ! -------                                                  --------    .
    implicit none
    ! -------                                                  --------    .
    integer(KIND=wi), intent (in) :: l
    type(densityplot_one), intent (inout) :: param1, param2, mis           ! les deux paramètres
    integer(KIND=wi), intent (in) :: deltaxymis                            ! nombre de discrétisation (entier)
    integer(KIND=wi), intent (in) :: nbparam                               ! nombre de modèles
    character(LEN=7), intent (in) :: filename                              ! base du nom des fichiers de sorties
    real(KIND=wr), intent (out) :: delta_1, delta_2                        ! pas de disrétisation (réels)
    ! -------                                                  --------    .
    integer(KIND=wi) :: i,j,k,n
    ! -------                                                  --------    .
    real(KIND=wr) :: themax,themin
    real(KIND=wr), dimension(:,:,:), allocatable :: tabmis                 ! grille pour la représentation 2D de la fonction coût
    real(KIND=wr), dimension(:,:), allocatable :: tabmistri,tabmistribis   ! points sélectionnés et triés
    character (LEN=5) :: numberfile
    ! -----------------------------------------------------------------    .
    write(numberfile(1:5),'(i5)')l
    allocate(tabmis(deltaxymis,deltaxymis,3))
    do i=1,deltaxymis
      do j=1,deltaxymis
        tabmis(i,j,3)=100000.0_wr                                          ! initialisation de grille de la fonction coût
      enddo
    enddo
    ! ------- delta de pas                                     --------    .
    delta_1=(param1%themax-param1%themin)/real(deltaxymis,wr)
    delta_2=(param2%themax-param2%themin)/real(deltaxymis,wr)
    ! -----------------------------------------------------------------    .
    ! ------- garde le meilleur modèle dans chaque case        --------    .
    do k=1,nbparam
      i=int((param1%vec(k)-param1%themin)/delta_1)+1
      j=int((param2%vec(k)-param2%themin)/delta_2)+1
      if ((i.ge.1).and.(i.le.deltaxymis)) then
        if ((j.ge.1).and.(j.le.deltaxymis)) then
          if (tabmis(i,j,3).gt.mis%vec(k)) then
            tabmis(i,j,1)=param1%vec(k)
            tabmis(i,j,2)=param2%vec(k)
            tabmis(i,j,3)=mis%vec(k)
          endif
        else
          write(*,*)'problème dans mkfcoutplot : calcul de densité, incrément incorrect, boucle 1',i,i,k
        endif
      else
        write(*,*)'problème dans mkfcoutplot : calcul de densité, incrément incorrect, boucle 2',i,i,k
      endif
    enddo
    ! ------- nombre de cases retenues                          --------    .
    n=0
    do i=1,deltaxymis
      do j=1,deltaxymis
        if ((tabmis(i,j,3).lt.10000.0_wr).and.(tabmis(i,j,3).gt.0.0_wr)) n=n+1
      enddo
    enddo
    ! ------- transforme la grille en points                   --------    .
    allocate(tabmistri(n,3))                                               ! non triés
    allocate(tabmistribis(n,3))                                            ! triés
    k=0
    do i=1,deltaxymis
      do j=1,deltaxymis
        if ((tabmis(i,j,3).lt.10000.0_wr).and.(tabmis(i,j,3).gt.0.0_wr)) then
          k=k+1
          tabmistri(k,1)=tabmis(i,j,1)
          tabmistri(k,2)=tabmis(i,j,2)
          tabmistri(k,3)=tabmis(i,j,3)
        endif
      enddo
    enddo
    ! ------- normalisation (min max entre 0 et 100%)          --------    .
    themax=-1.0_wr
    themin=100000._wr
    do i=1,n
      if (themax.lt.tabmistri(i,3)) themax=tabmistri(i,3)
      if (themin.gt.tabmistri(i,3)) themin=tabmistri(i,3)
    enddo
    do i=1,n
      tabmistri(i,3)= 100.0_wr/(themax-themin)*tabmistri(i,3) - 100.0_wr/(themax-themin)*themin
    enddo
    ! ------- tri                                              --------    .
    call tri_bulle(tabmistri,n,tabmistribis)
    ! ------- ecriture dans un fichier des points              --------    .
    open(unit=607,file="OUTPUT/GMT/"//filename//"-"//trim(adjustl(numberfile))//"mis.bin",STATUS="replace",access='direct',RECL=24)
    k=0
    do i=n,1,-1                                                            ! tri décroissant
      k=k+1
      write(607,rec=k)real(tabmistribis(i,1),8),real(tabmistribis(i,2),8),real(tabmistribis(i,3),8)
    enddo
    close(607)
    ! -----------------------------------------------------------------    .
    deallocate(tabmis,tabmistri,tabmistribis)
    ! -----------------------------------------------------------------    .

  end subroutine mkfcoutplot

END MODULE figure_GMT



! *********************************************************************    .
! *********************************************************************    .


