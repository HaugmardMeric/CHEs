! subroutine permettant la création d'un script LaTeX
! octobre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE latexscript

    use modparam
    use typetemps
    use time
    use datalecture, only : catalogue
    use sub_param, only : moy_ec
    use figure_GMT, only : affiche_temps_ref

    implicit none

    private

    public  :: latexfull


    ! -----------------------------------------------------------------    .
    ! on défini _FILE_DIR_ en fonction du compilateur
    ! variable permettant de tester l'existance du dossier

#ifdef __INTEL_COMPILER
#define _FILE_DIR_ DIRECTORY
#elif __GFORTRAN__
#define _FILE_DIR_ FILE
#endif

    ! -----------------------------------------------------------------    .

CONTAINS

! ---------------------------------------------------------------------    .
  subroutine latexfull(dc,dp,xmin,xmax,nbChaineMV,acceptance,pbest,misfit,E,nomsta,acentroid,nbtps,D)
    ! -------                                                  --------    .mh
    ! création de n scripts LaTeX pour la sortie des tableaux et des figures ...
    ! -------                                                  --------    .
    use sub_param, only : lectparam
    ! -------                                                  --------    .
    implicit none
    type(densityplot), intent (in) :: dp                                   ! modèles retenus par McMC
    type(coldmoy), intent (in) :: dc                                       ! modèles du coldrun
    real(KIND=wr), intent (inout) :: xmin(nbseismes), xmax(nbseismes)      ! cercles pond.
    integer(KIND=wi), intent (in) :: nbChaineMV
    type(accept), intent (in) :: acceptance(nbChaineMV)                    ! acceptance
    type(parametres), intent (inout) :: pbest(nbChaineMV)
    type(fcout), intent (in) :: misfit (nbChaineMV)                        ! fonction coût
    type(ellip), intent(in) :: E(nbseismes)
    character(LEN=4), dimension(:), intent(in) :: nomsta
    type(amoho_centroid), intent(in)  :: acentroid                         ! si moho non tabulaire
    integer(KIND=wi), intent(in) :: nbtps(nbseismes)                       ! nombre de données de temps
    type(dataall), intent(inout) :: D(nbseismes)                           ! données
    ! -------                                                  --------    .
    integer(KIND=wi) :: i
    integer(KIND=wi) :: nbChaineMVCold, nbChaineMVhot, maxiterhot, maxitercold
    ! -----------------------------------------------------------------    .
    call lectparam(nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold,chut=.true.)
    ! -------                                                  --------    .
    do i=1, nbseismes
      call latexone(i,dc,dp,xmin,xmax,nbChaineMVCold,nbChaineMVhot,maxiterhot,maxitercold, &
        acceptance,pbest,misfit,E(i),nbtps,D,acentroid)
    enddo
    if (FLAGresSTA) call mk_sta(nomsta)                                    ! nécessite GMT_resSTA dans mkGMT.f90
    ! -----------------------------------------------------------------    .
  end subroutine latexfull

! ---------------------------------------------------------------------    .

  subroutine mk_sta(nomsta)
    ! -------                                                  --------    .mh
    ! création d'un script LaTeX pour la sortie des résidus aux stations ...
    ! -------                                                  --------    .
    implicit none
    character(LEN=4), dimension(:), intent(in) :: nomsta
    ! -------                                                  --------    .
    integer(KIND=wi) :: i, ok
    character (LEN=500) :: a_char
    ! -----------------------------------------------------------------    .
    open(unit=905,file="OUTPUT/LATEX/stations.tex",STATUS="replace")
    write(905,*)"%latex plan.tex ; pdflatex plan.tex"
    write(905,*)"\documentclass[11pt,a4paper]{article}"
    write(905,*)"\usepackage[T1]{fontenc} \usepackage[frenchb]{babel}"
    write(905,*)"\usepackage{geometry}  \geometry{top=20mm, bottom=20mm, left=20mm, right=20mm}"
    write(905,*)"\usepackage{lscape, tabularx,array,multirow,pdflscape}"
    write(905,*)"\newcommand{\kms}{km$\cdotp$s{\footnotesize $^{-1}$}}"
    write(905,*)"\usepackage{gensymb} % degree"
    write(905,*)"\usepackage[pdftex]{graphicx} \DeclareGraphicsExtensions{.jpg,.pdf,.png,.bmp,.jpeg,.ps,.eps}"
    write(905,*)"\usepackage[clockwise]{rotating}"
    write(905,*)"\usepackage{longtable}"
    write(905,*)"\usepackage[np]{numprint} \npdecimalsign{,}"
    write(905,*)"\usepackage{datetime}"
    write(905,*)"\usepackage{lastpage}"
    write(905,*)"\usepackage{fancyhdr}"
    write(905,*)"\usepackage{color}"
    write(905,*)"\definecolor{dkblue}{rgb}{0,0,0.55}"
    write(905,*)"\pagestyle{fancy}"
    write(905,*)"\fancyhead[R]{\today ~ -- ~\currenttime}"
    write(905,*)"\fancyhead[L]{M\'eric Haugmard (\href{mailto:meric.haugmard@univ-nantes.fr}{meric.haugmard@univ-nantes.fr})}"
    write(905,*)"\fancyfoot[L] {\thepage /\pageref{LastPage}}"
    write(905,*)"\fancyfoot[C] {}"
    write(905,*)"\fancyfoot[R] {\href{http://lpgnantes.fr/haugmard-m}{http://lpgnantes.fr/haugmard-m} ; tel : 02 51 12 54 31 }"
    write(905,*)"\usepackage{hyperref}   % pdf interactif"
    write(905,*)"\hypersetup"
    write(905,*)"{"
    write(905,*)"pdftitle={Th\`ese de Doctorat},"
    write(905,*)"pdfauthor={M\'eric Haugmard},"
    write(905,*)"pdfsubject={sismicit\'e du Massif armoricain}"
    write(905,*)"pdfKeywords={sismicit\'e, Massif armoricain}"
    write(905,*)"pdfproducer={Laboratoire de Plan\'etologie et G\'eodynamique de Nantes}"
    write(905,*)"pdftoolbar=true, %barre d'outils non visible"
    write(905,*)"pdfmenubar=true, %barre de menu visible"
    write(905,*)"pdfpagelayout=TwoColumnLeft,"
    write(905,*)"pdfpagemode=UseThumbs,"
    write(905,*)"colorlinks=true,"
    write(905,*)"linkcolor= dkblue,"
    write(905,*)"filecolor= dkblue,"
    write(905,*)"urlcolor= dkblue"
    write(905,*)"}"
    write(905,*)"\usepackage[labelformat=empty]{caption}"
    write(905,*)"\newcommand{\PN}{$P_{n}$}"
    write(905,*)"\newcommand{\PG}{$P_{g}$}"
    write(905,*)"\newcommand{\SN}{$S_{n}$}"
    write(905,*)"\newcommand{\SG}{$S_{g}$}"
    ! -----------------------------------------------------------------    .
    write(905,*)"\begin{document}"
    do i=1,size(nomsta)
      write(905,*)"\centering \begin{landscape} \centering"
      write(905,*)"\vfill"
      write(905,*)"\begin{figure}[!ht]   	%%%%%%%%%%%% FIGURE reshisto %%%%%%%%%%"
      write(905,*)"\centerline{\includegraphics[width=0.65\linewidth,angle=-90]{../figures/resSTA-"//nomsta(i)//".pdf}}"
       write(905,*)"\caption{\large r\'esidus pour l'ensemble des \np{1000} meilleurs mo\`eles, par station : "//nomsta(i)//"}"
      write(905,*)"\end{figure}"
      write(905,*)"\vfill"
      write(905,*)"\end{landscape}"
    enddo
    ! -----------------------------------------------------------------    . 
    if (FLAGresSTA) then
      do i=1,size(nomsta)
        write(905,*)"\centering \begin{landscape} \centering"
        write(905,*)"\vfill"
        write(905,*)"\begin{figure}[!ht]   	%%%%%%%%%%%% FIGURE reshisto %%%%%%%%%%"
        write(905,*)"\centerline{\includegraphics[width=0.65\linewidth,angle=-90]{../figures/resSTA_TOT-"//nomsta(i)//".pdf}}"
        write(905,*)"\caption{\large r\'esidus pour toutes les it\'erations : "//nomsta(i)//"}"
        write(905,*)"\end{figure}"
        write(905,*)"\vfill"
        write(905,*)"\end{landscape}"
      enddo
      ! ---------------------------------------------------------------    .
      write(905,'(a)')"\begin{landscape} \centering"
      write(905,'(a)')"\renewcommand{\arraystretch}{1.75}"
      write(905,'(a)')"\begin{longtable}[!ht]"
      write(905,'(3a)')"{|>{\centering}p{2.5cm}|>{\centering}p{2.9cm}|>{\centering}p{1.5cm}|>{\centering}", &
        "p{2.9cm}|>{\centering}p{1.5cm}|>{\centering}p{2.9cm}|>{\centering}p{1.5cm}|", &
        ">{\centering}p{2.9cm}|p{1.5cm}<{\centering}|}"
      write(905,'(a)')"% Entete de la première page"
      write(905,'(2a)')"\hline"
      write(905,'(2a)')"\multirow{2}{*}{station} & \multicolumn{2}{c|}{\PG} & \multicolumn{2}{c|}{\PN} & ", &
        " \multicolumn{2}{c|}{\SG} & \multicolumn{2}{c|}{\SN} \\ "
      write(905,'(a)')"\cline{2-9} "
      write(905,'(2a)')" & moy. {\small($\pm2\sigma$)} & m\`ed. & moy. {\small($\pm2\sigma$)} & m\`ed. & ", &
        " moy. {\small($\pm2\sigma$)} & m\`ed. & moy. {\small($\pm2\sigma$)} & m\`ed. \\ "
      write(905,'(a)')"\hline"
      write(905,'(a)')"\bf & & & & & & & & \\"
      write(905,'(a)')"\endfirsthead"
      write(905,'(a)')"% Entête de toutes les pages"
      write(905,'(a)')"\hline"
      write(905,'(2a)')"\multirow{2}{*}{station} & \multicolumn{2}{c|}{\PG} & \multicolumn{2}{c|}{\PN} & ", &
        " \multicolumn{2}{c|}{\SG} & \multicolumn{2}{c|}{\SN} \\ "
      write(905,'(a)')"\cline{2-9} "
      write(905,'(2a)')" & moy. {\small($\pm2\sigma$)} & m\`ed. & moy. {\small($\pm2\sigma$)} & m\`ed. & ", &
        " moy. {\small($\pm2\sigma$)} & m\`ed. & moy. {\small($\pm2\sigma$)} & m\`ed. \\ "
      write(905,'(a)')"\hline"
      write(905,'(a)')"\bf & & & & & & & & \\"
      write(905,'(a)')"\bf ... & ... & ... & ... & ... & ... & ... & ... & ... \\"
      write(905,'(a)')"\bf & & & & & & & & \\"
      write(905,'(a)')"\endhead"
      write(905,'(a)')"% Bas de toutes les pages"
      write(905,'(a)')"\bf ... & ... & ... & ... & ... & ... & ... & ... & ... \\"
      write(905,'(a)')"\bf & & & & & & & & \\"
      write(905,'(a)')"\hline"
      write(905,'(a)')"\endfoot"
      write(905,'(a)')"\endlastfoot"
      write(905,'(a)')"% Contenu du tableau"
      ! ---------------------------------------------------------------    .
      ok=0
      open(902, FILE = "OUTPUT/GMT/sta_RES_TOT2latex.txt",status='old',iostat = ok)
      if (ok .ne. 0) then
        write(*,*)"problème dans mk_sta : le fichier OUTPUT/GMT/sta_RES_TOT2latex.txt n''existe pas "
        stop
      endif
      do while(ok .eq. 0)
        read(902,'(a500)',iostat = ok)a_char
        if(ok .eq. 0) write(905,'(a)')a_char
      enddo
      write(905,*)"\hline "
      write(905,*)"\end{longtable}"
      write(905,*)"\end{landscape}"
      close(902)
      ! ---------------------------------------------------------------    . affiche la carte des résidus aux stations, moy.
      write(905,*)"\begin{landscape} \centering"
      write(905,*)"\begin{figure}[!ht]   	%%%%%%%%%%%% FIGURE residus %%%%%%%%%%"
      write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
        "{../figures/resTOTPgmoy.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\hfill"
       write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
        "{../figures/resTOTSgmoy.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\vfill"
      write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
        "{../figures/resTOTPnmoy.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\hfill"
      write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
        "{../figures/resTOTSnmoy.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\caption{moyenne ($\pm2\sigma$) des r\'esidus pour tous les s\'eismes ", &
        "-- la taille de l'histogramme est proportionnelle au crit\`ere R$_p$ ", &
        "(R$_p$ = 1, si la distribution est gaussienne, 0 sinon)}"
      write(905,*)"\end{figure}"
      write(905,*)"\end{landscape}"
      ! ---------------------------------------------------------------    . affiche la carte des résidus aux stations, med.
      write(905,*)"\begin{landscape} \centering"
      write(905,*)"\begin{figure}[!ht]   	%%%%%%%%%%%% FIGURE residus %%%%%%%%%%"
      write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
        "{../figures/resTOTPgmed.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\hfill"
      write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
      "{../figures/resTOTSgmed.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\vfill"
      write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
      "{../figures/resTOTPnmed.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\hfill"
      write(905,*)"\begin{minipage}[t]{0.47\linewidth}"
      write(905,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
        "{../figures/resTOTSnmed.pdf}}"
      write(905,*)"\end{minipage}"
      write(905,*)"\caption{m\'ediane des r\'esidus pour tous les s\'eismes}"
      write(905,*)"\end{figure}"
      write(905,*)"\end{landscape}"

      ! ---------------------------------------------------------------    . histo des résidus aux stations
      write(905,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE matCorr %%%%%%%%%"
      write(905,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/Allres.pdf}}"
      write(905,*)"\definecolor{ss}{rgb}{1,0,1}"
      write(905,*)"\definecolor{pp}{rgb}{0,0.67,0.92}"
      write(905,*)"\caption{\large histogramme emplil\'e des r\'esidus pour tous les s\'eismes (ondes {\textcolor{pp}{\PG}}, ", &
        "{\textcolor{blue}{\PN}}, {\textcolor{ss}{\SG}} et {\textcolor{red}{\SN}}). }"
      write(905,*)"\end{figure}"

    endif
    write(905,*)"\end{document}"
    ! -----------------------------------------------------------------    .
    close(905)
    ! -----------------------------------------------------------------    .
  end subroutine mk_sta

! ---------------------------------------------------------------------    .

  subroutine latexone(i,dc,dp,xmin,xmax,nbChaineMVCold,nbChaineMV,maxiter,maxitercold, &
                      acceptance,pbest,misfit,E,nbtps,D,acentroid)
    ! -------                                                  --------    .mh
    ! création d'un script LaTeX pour la sortie des tableaux et des figures (pour un unique séisme)
    ! -------                                                  --------    .
    use invGEIGER
    use sub_misfit
    use pb_direct
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: i
    type(densityplot), intent (in) :: dp                                   ! modèles retenus par McMC
    type(coldmoy), intent (in) :: dc                                       ! modèles du coldrun
    real(KIND=wr), intent (inout) :: xmin(nbseismes), xmax(nbseismes)         ! cercles pond.
    integer(KIND=wi), intent (in) :: nbChaineMVCold,nbChaineMV,maxiter,maxitercold
    type(accept), intent (in) :: acceptance(nbChaineMV)                    ! acceptance
    type(parametres), intent (inout) :: pbest(nbChaineMV)
    type(fcout), intent (in) :: misfit (nbChaineMV)                        ! fonction coût
    type(ellip), intent(in) :: E
    type(amoho_centroid), intent(in)  :: acentroid                         ! si moho non tabulaire
    integer(KIND=wi), intent(in) :: nbtps(nbseismes)                       ! nombre de données de temps
    type(dataall), intent(inout) :: D(nbseismes)                           ! données
    ! -------                                                  --------    .
    type(seismes) :: refseisme(2)
    type(parametre) :: onep
    type(parametres) :: param_best
    type(date_sec) :: tpsref
    logical :: existe1, chut, critique
    character(LEN=30) :: char_0
    character(LEN=20) :: charname
    character(LEN=4) :: staname
    character(LEN=2) :: ondetype
    character(LEN=33) :: nomfile
    integer(KIND=wi) :: j,k,l,ok,x1,x2,Noldtime,Nnewtime,ratetime,find
    real(KIND=wr) :: moy,val,val2,ec,vec(nbChaineMV),a_ponderation,residu,atime
    real(KIND=wr) :: pct,atimemoy,atimeec,tl,dist,valec,valec2,tps,depi
    real(KIND=wr) :: misfitval
    character (LEN=5) :: numberfile
    character (LEN=1) :: mod
    ! -----------------------------------------------------------------    .
    write(numberfile(1:5),'(i5)')i
    write(*,*)"ecriture du script LaTeX"
    call system_clock(Noldtime)
    ! -------                                                  --------    .
    write(nomfile(1:13),'(a13)')"OUTPUT/LATEX/"
    write(nomfile(14:18),'(i4.4)')dp%temps_ref(i)%date%year
    write(nomfile(18:19),'(a1)')"."
    write(nomfile(19:21),'(i2.2)')dp%temps_ref(i)%date%month
    write(nomfile(21:22),'(a1)')"."
    write(nomfile(22:24),'(i2.2)')dp%temps_ref(i)%date%day
    write(nomfile(24:25),'(a1)')"T"
    write(nomfile(25:27),'(i2.2)')dp%temps_ref(i)%date%hour
    write(nomfile(27:28),'(a1)')"h"
    write(nomfile(28:30),'(i2.2)')dp%temps_ref(i)%date%min
    write(nomfile(30:31),'(a1)')"."
    write(nomfile(31:33),'(i2.2)')int(dp%Tzero(i)%best,wi)
    open(unit=900,file=trim(adjustl(nomfile))//"-"//trim(adjustl(numberfile))//".tex",STATUS="replace")
    write(900,*)"%latex plan.tex ; pdflatex plan.tex"
    write(900,*)"%",dp%temps_ref(i)
    write(900,*)"\documentclass[11pt,a4paper]{article}"
    write(900,*)"\usepackage[T1]{fontenc} \usepackage[frenchb]{babel}"
    write(900,*)"\usepackage{geometry}  \geometry{top=20mm, bottom=20mm, left=20mm, right=20mm}"
    write(900,*)"\usepackage{lscape, tabularx,array,multirow,pdflscape}"
    write(900,*)"\newcommand{\kms}{km$\cdotp$s{\footnotesize $^{-1}$}}"
    write(900,*)"\usepackage{gensymb} % degree"
    write(900,*)"\usepackage[pdftex]{graphicx} \DeclareGraphicsExtensions{.jpg,.pdf,.png,.bmp,.jpeg,.ps,.eps}"
    write(900,*)"\usepackage[clockwise]{rotating}"
    write(900,*)"\usepackage[labelformat=empty]{caption}"
    write(900,*)"\usepackage{longtable}"
    write(900,*)"\usepackage[np]{numprint} \npdecimalsign{,}"
    write(900,*)"\usepackage{datetime}"
    write(900,*)"\usepackage{lastpage}"
    write(900,*)"\usepackage{fancyhdr}"
    write(900,*)"\usepackage{color}"
    write(900,*)"\definecolor{dkblue}{rgb}{0,0,0.55}"
    write(900,*)"\pagestyle{fancy}"
    write(900,*)"\fancyhead[R]{\today ~ -- ~\currenttime}"
    write(900,*)"\fancyhead[L]{M\'eric Haugmard (\href{mailto:meric.haugmard@univ-nantes.fr}{meric.haugmard@univ-nantes.fr})}"
    write(900,*)"\fancyfoot[L] {\thepage /\pageref{LastPage}}"
    write(900,*)"\fancyfoot[C] {}"
    write(900,*)"\fancyfoot[R] {\href{http://lpgnantes.fr/haugmard-m}{http://lpgnantes.fr/haugmard-m} ; tel : 02 51 12 54 31 }"
    write(900,*)"\usepackage{hyperref}   % pdf interactif"
    write(900,*)"\hypersetup"
    write(900,*)"{"
    write(900,*)"pdftitle={},"
    write(900,*)"pdfauthor={M\'eric Haugmard},"
    write(900,*)"pdfsubject={Th\`ese - sismicit\'e du Massif armoricain }"
    write(900,*)"pdfKeywords={sismicit\'e, Vannes, Massif armoricain}"
    write(900,*)"pdfproducer={Laboratoire de Plan\'etologie et G\'eodynamique de Nantes}"
    write(900,*)"pdftoolbar=true, %barre d'outils non visible"
    write(900,*)"pdfmenubar=true, %barre de menu visible"
    write(900,*)"pdfpagelayout=TwoColumnLeft,"
    write(900,*)"pdfpagelayout=TwoColumnLeft,"
    write(900,*)"pdfpagemode=UseThumbs,"
    write(900,*)"linkcolor= dkblue,"
    write(900,*)"filecolor= dkblue,"
    write(900,*)"urlcolor= dkblue"
    write(900,*)"}"
    write(900,*)"\newcommand{\PN}{$P_{n}$}"
    write(900,*)"\newcommand{\PG}{$P_{g}$}"
    write(900,*)"\newcommand{\SN}{$S_{n}$}"
    write(900,*)"\newcommand{\SG}{$S_{g}$}"
    ! -----------------------------------------------------------------    .
    write(900,*)"\begin{document}"
    ! -----------------------------------------------------------------    .
    write(900,*)
    call affiche_temps_ref(dp%temps_ref(i),char_0,1)
    write(900,*)"\centering{\NoAutoSpaceBeforeFDP \Huge S\'eisme du ",char_0,"} \\"
    write(900,*)" "
    write(900,*)"\vspace{1.0cm}"
    write(900,*)" "
    write(900,*)"\sloppy \raggedright"
    ! -----------------------------------------------------------------    .
    do j=1,nbChaineMV
      vec(j)=acceptance(j)%val
    enddo
    call moy_ec (vec,nbChaineMV,nbChaineMV,moy,ec)
    write(900,*)"Param\`etres de l'inversion:"
    write(900,'(a,f9.2,a,f5.2,a)')"moyenne des acceptances : \np{",moy,"} $\pm$ \np{",ec,"} \\"
    write(900,'(a,i12,a)')"nombre de cha\^ines de Markov (coldrun) : \np{",nbChaineMVCold,"}\\"
    write(900,'(a,i12,a)')"nombre d'it\'erations par cha\^ine (coldrun) : \np{",maxitercold,"}\\"
    write(900,'(a,i12,a)')"nombre de cha\^ines de Markov (hotrun) : \np{",nbChaineMV,"}\\"
    write(900,'(a,i12,a)')"nombre d'it\'erations par cha\^ine (hotrun) : \np{",maxiter,"}\\"
    write(900,'(a,i12,a)')"nombre mod\`eles test\'es : \np{",maxiter*nbChaineMV+nbChaineMVCold*maxitercold,"}\\"
    write(900,'(a,i12,a)')"nombre de mod\`eles retenus : \np{",dp%nbparam,"}\\"
    write(900,'(a,i6,a)')"discr\'etisation pour le diagramme de densit\'e : \np{",dp%deltaxy,"}\\"
    x1=int(xmin(i),wi)
    x2=int(xmax(i),wi)
    write(900,'(a,i5,a,i5,a)')"cercles de pond\'erations (km) : \np{",x1,"} et \np{",x2,"}\\"
    ! -----------------------------------------------------------------    . ellipse
    write(900,'(a)')"~\\"
    write(900,'(a)')"ellipse (1$\sigma$) des 1000 meilleurs mod\`eles"
    write(900,'(a,f8.2,a)')" azimuth : \np[\degree]{",E%ang,"}\\"
    write(900,'(a,f10.2,a)')" demi axe a : \np[m]{",E%axeA*1000.0_wr,"}\\"
    write(900,'(a,f10.2,a)')" demi axe b : \np[m]{",E%axeB*1000.0_wr,"}\\"
    write(900,'(a,f10.2,a,f13.2,a)')" Aire : \np[km^2]{",E%axeB*E%axeA*pi,"} {\small (\np[ha]{",E%axeB*E%axeA*pi/0.01_wr,"})} \\"
    write(900,'(a)')"~\\"
    ! -----------------------------------------------------------------    . recherche meilleur modèle
    moy=100000.0_wr
    do j=1,nbChaineMV
      if (misfit(j)%best.lt.moy) then
      moy = misfit(j)%best
      k = j
      endif
    enddo
    call mvPall_2_P1(onep,pbest(k),i)
    call catalogue(onep,refseisme,find)                                    ! comparaison avec le catalogue
    ! -------                                                  --------    .
    if (find==1) then
      write(900,*)"{ \small S\'eisme pr\'esent dans le catalogue : "//refseisme(1)%name//"\\"
      write(900,*)"\begin{itemize}"
      write(900,'(a,f9.2,a3)')"\item magnitude $M_l$ : \np{",refseisme(1)%mag,"}\\"
      write(900,'(a55,f9.2,a3)')"\item longitude : \np[\degree]{",refseisme(1)%lon,"}\\"
      write(900,'(a55,f9.2,a3)')"\item latitude : \np[\degree]{",refseisme(1)%lat,"}\\"
      write(900,'(a55,f9.2,a3)')"\item profondeur hypocentre : \np[km]{",refseisme(1)%pfd,"}\\"
      call affiche_temps_ref(refseisme(1)%tps_init,char_0,1)
      write(900,'(a45,a30,a8,f4.1,a13)')"\item temps initial : \NoAutoSpaceBeforeFDP",char_0," et \np{", &
        refseisme(1)%tps_init%sec,"} secondes \\"
      write(900,'(a65,f9.3,a4)')"\item diff\'erence de temps avec le meilleur mod\`ele : \np[s]{",refseisme(1)%d_t,"} \\"
      write(900,'(a75,f9.2,a4)')"\item diff\'erence de profondeur avec le meilleur mod\`ele : \np[km]{",refseisme(1)%d_p,"} \\"
      write(900,'(a65,f12.3,a4)')"\item distance \'epicentrale : \np[m]{",refseisme(1)%d_epi*1000.0_wr,"} \\"
      write(900,*)"\end{itemize}}"

    elseif (find==2) then
      write(900,*)"{ \small S\'eisme pr\'esent dans le catalogue 2 fois : \\"
      write(900,*)"catalogue 1 : "//refseisme(1)%name//"\\"
      write(900,*)"\begin{itemize}"
      write(900,'(a,f9.2,a3)')"\item magnitude $M_l$ : \np{",refseisme(1)%mag,"}\\"
      write(900,'(a55,f9.2,a3)')"\item longitude : \np[\degree]{",refseisme(1)%lon,"}\\"
      write(900,'(a55,f9.2,a3)')"\item latitude : \np[\degree]{",refseisme(1)%lat,"}\\"
      write(900,'(a55,f9.2,a3)')"\item profondeur hypocentre : \np[km]{",refseisme(1)%pfd,"}\\"
      call affiche_temps_ref(refseisme(1)%tps_init,char_0,1)
      write(900,'(a45,a30,a8,f4.1,a13)')"\item temps initial : \NoAutoSpaceBeforeFDP",char_0," et \np{", &
        refseisme(1)%tps_init%sec,"} secondes \\"
      write(900,'(a65,f9.3,a4)')"\item diff\'erence de temps avec le meilleur mod\`ele : \np[s]{",refseisme(1)%d_t,"} \\"
      write(900,'(a75,f9.2,a4)')"\item diff\'erence de profondeur avec le meilleur mod\`ele : \np[km]{",refseisme(1)%d_p,"} \\"
      write(900,'(a65,f12.3,a4)')"\item distance \'epicentrale : \np[m]{",refseisme(1)%d_epi*1000.0_wr,"} \\"
      write(900,*)"\end{itemize}"
      write(900,*)"~\\"
      write(900,*)"catalogue 2 : "//refseisme(2)%name//"\\"
      write(900,*)"\begin{itemize}"
      write(900,'(a,f9.2,a3)')"\item magnitude $M_l$ : \np{",refseisme(2)%mag,"}\\"
      write(900,'(a55,f9.2,a3)')"\item longitude : \np[\degree]{",refseisme(2)%lon,"}\\"
      write(900,'(a55,f9.2,a3)')"\item latitude : \np[\degree]{",refseisme(2)%lat,"}\\"
      write(900,'(a55,f9.2,a3)')"\item profondeur hypocentre : \np[km]{",refseisme(2)%pfd,"}\\"
      call affiche_temps_ref(refseisme(2)%tps_init,char_0,1)
      write(900,'(a45,a30,a8,f4.1,a13)')"\item temps initial : \NoAutoSpaceBeforeFDP",char_0," et \np{", &
        refseisme(2)%tps_init%sec,"} secondes \\"
      write(900,'(a65,f9.3,a4)')"\item diff\'erence de temps avec le meilleur mod\`ele : \np[s]{",refseisme(2)%d_t,"} \\"
      write(900,'(a75,f9.2,a4)')"\item diff\'erence de profondeur avec le meilleur mod\`ele : \np[km]{",refseisme(2)%d_p,"} \\"
      write(900,'(a65,f12.3,a4)')"\item distance \'epicentrale : \np[m]{",refseisme(2)%d_epi*1000.0_wr,"} \\"
      write(900,*)"\end{itemize}}"
    else
      write(900,*)"S\'eisme absent du catalogue"
    endif
    write(900,*)
    ! -----------------------------------------------------------------    .
    write(900,*)"\begin{figure}[!ht] \centering   	%%%%%%%%% FIGURE MAP init %%%%%%%%"
    write(900,*)"\centerline{\includegraphics[width=0.99\linewidth]", &
      "{../figures/init-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\caption{\large prior sur l'\'epicentre, m\'ethode des arriv\'ees les plus proches}"
    write(900,*)"\end{figure}"
    write(900,*)
    ! -----------------------------------------------------------------    . affiche le diagramme de chatelain
    write(900,*)"\begin{landscape} \centering"
    write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE chatelain %%%%%%%%"
    write(900,*)"\centerline{\includegraphics[width=.6\linewidth,angle=-90]", &
      "{../figures/chatelainplot"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\caption{\large diagramme de Ch\^atelain }"
    write(900,*)"\end{figure}"
    write(900,*)"\end{landscape}"
    write(900,*)
    ! -----------------------------------------------------------------    .
    write(900,*)"\begin{landscape} \centering"
    write(900,*)"\begin{figure}[!ht] \centering   	%%%%%%%%%%%% FIGURE MAP %%%%%%%%%%"
    write(900,*)"\centerline{\includegraphics[width=0.99\textwidth,angle=-90]", &
      "{../figures/map"//"-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{figure}"
    write(900,*)"\end{landscape}"
    write(900,*)
    ! -----------------------------------------------------------------    .
    write(900,*)"\begin{landscape} \centering   	%%%%%%%%%% Table données %%%%%%%%%"

    write(900,*)"\begin{table}[!ht] \scriptsize \centering  \renewcommand{\arraystretch}{1.5}"
    write(900,*)"\begin{tabular}{|>{\centering}m{2.5cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
    write(900,*)"|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
    write(900,*)"|>{\centering}m{2.2cm}|>{\centering}m{2.2cm}|m{2.1cm}<{\centering}|}"
    write(900,*)"\hline "
    write(900,*)"{\bf \large coldruns} & {fonction co\^ut} & {V$_C$ (\kms)} & {V$_M$ (\kms)} & ", &
      " {Z$_{moho}$ (km)} & {V$_{P}$V$_{S}$} & {Z$_{hypo}$ (km)} & {longitude (\degree)} & ", &
      " {latitude (\degree)} & {T$_{z\acute{e}ro}$ (s)} \\"
    write(900,*)"\hline "
    call difftime(atime,dc%tempsrefcold(i),dp%temps_ref(i)) !!!!!!!!!
    atimeec=dc%ectot%par%Tzero(i)%sec
    atimemoy=dc%moytot%par%Tzero(i)%sec + atime
    write(900,4997)" {moyenne ($\pm2\sigma$) des meilleurs mod\`eles de toutes les cha\^ines} & \np{",dc%moytot%mis, &
      "} $\pm$ \np{",2.0_wr*dc%ectot%mis,"} & \np{",dc%moytot%par%VC,"} $\pm$ \np{",2.0_wr*dc%ectot%par%VC, &
      "} & \np{",dc%moytot%par%VM,"} $\pm$ \np{",2.0_wr*dc%ectot%par%VM,"} & \np{",dc%moytot%par%Zmoho, &
      "} $\pm$ \np{",2.0_wr*dc%ectot%par%Zmoho,"} & \np{",dc%moytot%par%VpVs,"} $\pm$ \np{",2.0_wr*dc%ectot%par%VpVs, &
      "} & \np{",dc%moytot%par%Zhypo(i),"} $\pm$ \np{",2.0_wr*dc%ectot%par%Zhypo(i),"} & \np{",dc%moytot%par%lon(i), &
      "} $\pm$ \np{",2.0_wr*dc%ectot%par%lon(i),"} & \np{",dc%moytot%par%lat(i),"} $\pm$ \np{",2.0_wr*dc%ectot%par%lat(i), &
      "} & \np{",atimemoy,"} $\pm$ \np{",2.0_wr*atimeec,"} \\ "
    atimeec=dc%ecselect%par%Tzero(i)%sec
    atimemoy=dc%moyselect%par%Tzero(i)%sec + atime
    write(900,*)"\hline "
    write(900,4998)" {moyenne ($\pm2\sigma$) des meilleurs mod\`eles de chaque cha\^ine s\'electionn\'ee} & \np{", &
      dc%moyselect%mis,"} $\pm$ \np{",2.0_wr*dc%ecselect%mis,"} & \np{",dc%moyselect%par%VC,"} $\pm$ \np{", &
      2.0_wr*dc%ecselect%par%VC,"} & \np{",dc%moyselect%par%VM,"} $\pm$ \np{",2.0_wr*dc%ecselect%par%VM,"} & \np{", &
      dc%moyselect%par%Zmoho,"} $\pm$ \np{",2.0_wr*dc%ecselect%par%Zmoho,"} & \np{",dc%moyselect%par%VpVs, &
      "} $\pm$ \np{",2.0_wr*dc%ecselect%par%VpVs,"} & \np{",dc%moyselect%par%Zhypo(i),"} $\pm$ \np{", &
      2.0_wr*dc%ecselect%par%Zhypo(i),"} & \np{",dc%moyselect%par%lon(i),"} $\pm$ \np{",2.0_wr*dc%ecselect%par%lon(i), &
      "} & \np{",dc%moyselect%par%lat(i),"} $\pm$ \np{",2.0_wr*dc%ecselect%par%lat(i),"} & \np{",atimemoy, &
      "} $\pm$ \np{",2.0_wr*atimeec,"} \\ "
    4997 format (a128,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.1,a8,f9.3,a12,f9.2,a8,f9.2,a12,f9.1,a8, &
    f10.4,a12,f10.3,a8,f10.4,a12,f10.3,a8,f10.2,a12,f10.2,a4)
    4998 format (a140,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.1,a8,f9.3,a12,f9.2,a8,f9.2,a12,f9.1,a8, &
    f7.4,a12,f9.3,a8,f7.4,a12,f9.3,a8,f9.2,a12,f9.2,a4)
    write(900,*)"\hline "
    write(900,*)"\end{tabular}"
    write(900,*)"\end{table}"

    write(900,*)"\begin{table}[!ht] \scriptsize \centering  \renewcommand{\arraystretch}{1.5}"
    write(900,*)"\begin{tabular}{|>{\centering}m{2.5cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
    write(900,*)"|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
    write(900,*)"|>{\centering}m{2.2cm}|>{\centering}m{2.2cm}|m{2.1cm}<{\centering}|}"
    write(900,*)"\hline "
    write(900,*)" {\bf \large hotruns} & {fonction co\^ut} & {V$_C$ (\kms)} & {V$_M$ (\kms)} & {Z$_{moho}$ (km)} ", &
      " & {V$_{P}$V$_{S}$} & {Z$_{hypo}$ (km)} & {longitude (\degree)} & {latitude (\degree)} & {T$_{z\acute{e}ro}$ (s)} \\"
    write(900,*)"\hline "
    write(900,4999)" {mode} & \np{",dp%mis%mode,"} & \np{",dp%VC%mode,"} & \np{",dp%VM%mode, &
      "} & \np{",dp%Zmoho%mode,"} & \np{",dp%VpVs%mode,"} & \np{",dp%Zhypo(i)%mode,"} & \np{",dp%lon(i)%mode, &
      "} & \np{",dp%lat(i)%mode,"} & \np{",dp%Tzero(i)%mode,"} \\ "
    write(900,*)"\hline "
    write(900,5000)" {m\'ediane} & \np{",dp%mis%mediane,"} & \np{",dp%VC%mediane,"} & \np{",dp%VM%mediane, &
      "} & \np{",dp%Zmoho%mediane,"} & \np{",dp%VpVs%mediane,"} & \np{",dp%Zhypo(i)%mediane,"} & \np{",dp%lon(i)%mediane, &
      "} & \np{",dp%lat(i)%mediane,"} & \np{",dp%Tzero(i)%mediane,"} \\ "
    write(900,*)"\hline "
    write(900,5001)" {meilleur mod\`ele} & \np{",dp%mis%best,"} & \np{",dp%VC%best,"} & \np{",dp%VM%best, &
      "} & \np{",dp%Zmoho%best,"} & \np{",dp%VpVs%best,"} & \np{",dp%Zhypo(i)%best,"} & \np{",dp%lon(i)%best, &
      "} & \np{",dp%lat(i)%best,"} & \np{",dp%Tzero(i)%best,"} \\ "
    write(900,*)"\hline "
    write(900,5003)" {moyenne ($\pm2\sigma$) des 100 meilleurs mod\`eles} & \np{",dp%mis%moy_100,"} $\pm$ \np{", &
      2.0_wr*dp%mis%ec_100,"} & \np{",dp%VC%moy_100,"} $\pm$ \np{",2.0_wr*dp%VC%ec_100,"} & \np{",dp%VM%moy_100, &
      "} $\pm$ \np{",2.0_wr*dp%VM%ec_100,"} & \np{",dp%Zmoho%moy_100,"} $\pm$ \np{",2.0_wr*dp%Zmoho%ec_100, &
      "} & \np{",dp%VpVs%moy_100,"} $\pm$ \np{",2.0_wr*dp%VpVs%ec_100,"} & \np{",dp%Zhypo(i)%moy_100,"} $\pm$ \np{", &
      2.0_wr*dp%Zhypo(i)%ec_100,"} & \np{",dp%lon(i)%moy_100,"} $\pm$ \np{",2.0_wr*dp%lon(i)%ec_100,"} & \np{",dp%lat(i)%moy_100, &
      "} $\pm$ \np{",2.0_wr*dp%lat(i)%ec_100,"} & \np{",dp%Tzero(i)%moy_100,"} $\pm$ \np{",2.0_wr*dp%Tzero(i)%ec_100,"} \\ "
    write(900,*)"\hline "
    write(900,5004)" {moyenne ($\pm2\sigma$) des 1000 meilleurs mod\`eles} & \np{",dp%mis%moy_1000,"} $\pm$ \np{", &
      2.0_wr*dp%mis%ec_1000,"} & \np{",dp%VC%moy_1000,"} $\pm$ \np{",2.0_wr*dp%VC%ec_1000,"} & \np{",dp%VM%moy_1000, &
      "} $\pm$ \np{",2.0_wr*dp%VM%ec_1000,"} & \np{",dp%Zmoho%moy_1000,"} $\pm$ \np{",2.0_wr*dp%Zmoho%ec_1000,"} & \np{", &
      dp%VpVs%moy_1000,"} $\pm$ \np{",2.0_wr*dp%VpVs%ec_1000,"} & \np{",dp%Zhypo(i)%moy_1000,"} $\pm$ \np{", &
      2.0_wr*dp%Zhypo(i)%ec_1000,"} & \np{",dp%lon(i)%moy_1000,"} $\pm$ \np{",2.0_wr*dp%lon(i)%ec_1000,"} & \np{", &
      dp%lat(i)%moy_1000,"} $\pm$ \np{",2.0_wr*dp%lat(i)%ec_1000,"} & \np{",dp%Tzero(i)%moy_1000,"} $\pm$ \np{", &
      2.0_wr*dp%Tzero(i)%ec_1000,"} \\ "
    write(900,*)"\hline "
    write(900,5005)" {moyenne ($\pm2\sigma$) des 10000 meilleurs mod\`eles} & \np{",dp%mis%moy_10000,"} $\pm$ \np{", &
      2.0_wr*dp%mis%ec_10000,"} & \np{",dp%VC%moy_10000,"} $\pm$ \np{",2.0_wr*dp%VC%ec_10000,"} & \np{",dp%VM%moy_10000, &
      "} $\pm$ \np{",2.0_wr*dp%VM%ec_10000,"} & \np{",dp%Zmoho%moy_10000,"} $\pm$ \np{",2.0_wr*dp%Zmoho%ec_10000, &
      "} & \np{",dp%VpVs%moy_10000,"} $\pm$ \np{",2.0_wr*dp%VpVs%ec_10000,"} & \np{",dp%Zhypo(i)%moy_10000,"} $\pm$ \np{", &
      2.0_wr*dp%Zhypo(i)%ec_10000,"} & \np{",dp%lon(i)%moy_10000,"} $\pm$ \np{",2.0_wr*dp%lon(i)%ec_10000,"} & \np{", &
      dp%lat(i)%moy_10000,"} $\pm$ \np{",2.0_wr*dp%lat(i)%ec_10000,"} & \np{",dp%Tzero(i)%moy_10000, &
      "} $\pm$ \np{",2.0_wr*dp%Tzero(i)%ec_10000,"} \\ "
    write(900,*)"\hline "
    write(900,5006)" {moyenne ($\pm2\sigma$) des meilleurs mod\`eles de chaque cha\^ine} & \np{",dp%mis%moy_bestchaine, &
      "} $\pm$ \np{",2.0_wr*dp%mis%ec_bestchaine,"} & \np{",dp%VC%moy_bestchaine,"} $\pm$ \np{",2.0_wr*dp%VC%ec_bestchaine, &
      "} & \np{",dp%VM%moy_bestchaine,"} $\pm$ \np{",2.0_wr*dp%VM%ec_bestchaine,"} & \np{",dp%Zmoho%moy_bestchaine, &
      "} $\pm$ \np{",2.0_wr*dp%Zmoho%ec_bestchaine,"} & \np{",dp%VpVs%moy_bestchaine,"} $\pm$ \np{", &
      2.0_wr*dp%VpVs%ec_bestchaine,"} & \np{",dp%Zhypo(i)%moy_bestchaine,"} $\pm$ \np{",2.0_wr*dp%Zhypo(i)%ec_bestchaine, &
      "} & \np{",dp%lon(i)%moy_bestchaine,"} $\pm$ \np{",2.0_wr*dp%lon(i)%ec_bestchaine,"} & \np{", &
      dp%lat(i)%moy_bestchaine,"} $\pm$ \np{",2.0_wr*dp%lat(i)%ec_bestchaine,"} & \np{",dp%Tzero(i)%moy_bestchaine, &
      "} $\pm$ \np{",2.0_wr*dp%Tzero(i)%ec_bestchaine,"} \\ "
    write(900,*)"\hline "
      write(900,5002)" {moyenne ($\pm2\sigma$) totale} & \np{",dp%mis%moy_tot,"} $\pm$ \np{",2.0_wr*dp%mis%ec_tot, &
      "} & \np{",dp%VC%moy_tot,"} $\pm$ \np{",2.0_wr*dp%VC%ec_tot,"} & \np{",dp%VM%moy_tot, &
      "} $\pm$ \np{",2.0_wr*dp%VM%ec_tot,"} & \np{",dp%Zmoho%moy_tot,"} $\pm$ \np{",2.0_wr*dp%Zmoho%ec_tot, &
      "} & \np{",dp%VpVs%moy_tot,"} $\pm$ \np{",2.0_wr*dp%VpVs%ec_tot,"} & \np{",dp%Zhypo(i)%moy_tot, &
      "} $\pm$ \np{",2.0_wr*dp%Zhypo(i)%ec_tot,"} & \np{",dp%lon(i)%moy_tot,"} $\pm$ \np{",2.0_wr*dp%lon(i)%ec_tot, &
      "} & \np{",dp%lat(i)%moy_tot,"} $\pm$ \np{",2.0_wr*dp%lat(i)%ec_tot,"} & \np{",dp%Tzero(i)%moy_tot, &
      "} $\pm$ \np{",2.0_wr*dp%Tzero(i)%ec_tot,"} \\ "
    write(900,*)"\hline "
    4999 format (a14,f9.2,a8,f9.2,a8,f9.2,a8,f9.2,a8,f9.3,a8,f9.2,a8,f7.4,a8,f7.4,a8,f9.2,a4)
    5000 format (a19,f9.2,a8,f9.2,a8,f9.2,a8,f9.2,a8,f9.3,a8,f9.2,a8,f7.4,a8,f7.4,a8,f9.2,a4)
    5001 format (a27,f9.2,a8,f9.2,a8,f9.2,a8,f9.2,a8,f9.3,a8,f9.2,a8,f7.4,a8,f7.4,a8,f9.2,a4)
    5002 format (a87,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.1,a8,f9.3,a12,f9.2,a8,f9.2,a12,f9.1,a8, &
    f7.4,a12,f9.3,a8,f7.4,a12,f9.3,a8,f9.2,a12,f9.2,a4)
    5003 format (a128,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.1,a8,f9.3,a12,f9.2,a8,f9.2,a12,f9.1,a8, &
    f7.4,a12,f9.3,a8,f7.4,a12,f9.3,a8,f9.2,a12,f9.2,a4)
    5004 format (a129,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.1,a8,f9.3,a12,f9.2,a8,f9.2,a12,f9.1,a8, &
    f7.4,a12,f9.3,a8,f7.4,a12,f9.3,a8,f9.2,a12,f9.2,a4)
    5005 format (a80,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.1,a8,f9.3,a12,f9.2,a8,f9.2,a12,f9.1,a8, &
    f7.4,a12,f9.3,a8,f7.4,a12,f9.3,a8,f9.2,a12,f9.2,a4)
    5006 format (a123,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.2,a8,f9.2,a12,f9.1,a8,f9.3,a12,f9.2,a8,f9.2,a12,f9.1,a8, &
    f7.4,a12,f9.3,a8,f7.4,a12,f9.3,a8,f9.2,a12,f9.2,a4)
    write(900,*)"\end{tabular}"
    write(900,*)"\end{table}"


    ! -----------------------------------------------------------------    .
    ! modèle Arroucau
    ! -----------------------------------------------------------------    .
    write(900,*)" "
    write(900,*)"\input{../GMT/Arroucau_all-"//trim(adjustl(numberfile))//"}"
    write(900,*)" "
    ! -----------------------------------------------------------------    .
    ! modèle Si-Hex
    ! -----------------------------------------------------------------    .
    write(900,*)" "
    write(900,*)"\input{../GMT/SiHex_all-"//trim(adjustl(numberfile))//"}"
    write(900,*)" "
    ! -----------------------------------------------------------------    .
    ! modèle CEA
    ! -----------------------------------------------------------------    .
    write(900,*)" "
    write(900,*)"\input{../GMT/CEA_all-"//trim(adjustl(numberfile))//"}"
    write(900,*)" "
    ! -----------------------------------------------------------------    .



    ! -----------------------------------------------------------------    .
    ! ------- verif avec méthode de Geiger                     --------    .
    ! -----------------------------------------------------------------    .
    param_best%VC=dp%VC%best ! moy_100
    param_best%VM=dp%VM%best ! moy_100
    param_best%Zmoho=dp%Zmoho%best ! moy_100
    param_best%VpVs=dp%VpVs%best ! moy_100
    do j=1,nbseismes
      param_best%Zhypo(j)=dp%Zhypo(j)%best ! moy_100
      param_best%lon(j)=dp%lon(j)%best ! moy_100
      param_best%lat(j)=dp%lat(j)%best ! moy_100
      param_best%Tzero(j) = dp%temps_ref(j)
      param_best%Tzero(j)%sec = dp%Tzero(j)%best ! 100
      call basetime(param_best%Tzero(j))
    enddo
    call tempsTheoDirect(nbtps,param_best,D,critique,acentroid)
    call compute_misfit(nbtps,D,misfitval,xmin,xmax,'H')
    ! -------                                                  --------    .
    write(900,*)"\begin{table}[!ht] \scriptsize \centering  \renewcommand{\arraystretch}{1.5}"
    write(900,*)"\begin{tabular}{|>{\centering}m{2.5cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
    write(900,*)"|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}|>{\centering}m{2.1cm}"
    write(900,*)"|>{\centering}m{2.2cm}|>{\centering}m{2.2cm}|m{2.1cm}<{\centering}|}"
    write(900,*)"\hline "
    write(900,*)"{\bf \large geiger} & {fonction co\^ut} & {V$_C$ (\kms)} & {V$_M$ (\kms)} & ", &
      " {Z$_{moho}$ (km)} & {V$_{P}$V$_{S}$} & {Z$_{hypo}$ (km)} & {longitude (\degree)} & ", &
      " {latitude (\degree)} & {T$_{z\acute{e}ro}$ (s)} \\"
    write(900,*)"\hline "
    ! -------                                                  --------    .
    write(900,5010)" {entr\'ee} & \np{",misfitval,"} & \np{",param_best%VC,"} & \np{",param_best%VM, &
      "} & \np{",param_best%Zmoho,"} & \np{",param_best%VpVs,"} & \np{",param_best%Zhypo(i),"} & \np{",param_best%lon(i), &
      "} & \np{",param_best%lat(i),"} & \np{",param_best%Tzero(i)%sec,"} \\ "
    ! -------                                                  --------    .
    chut=.true.
    mod='I'
    call dogeiger(nbtps,D,param_best,acentroid,mod,chut)
    call tempsTheoDirect(nbtps,param_best,D,critique,acentroid)
    call compute_misfit(nbtps,D,misfitval,xmin,xmax,'H')
    ! -------                                                  --------    .
    write(900,*)"\hline "

    if (misfitval.lt.1000.0_wr)  then
      write(900,5010)" {sortie} & \np{",misfitval,"} & \np{",param_best%VC,"} & \np{",param_best%VM, &
        "} & \np{",param_best%Zmoho,"} & \np{",param_best%VpVs,"} & \np{",param_best%Zhypo(i),"} & \np{",param_best%lon(i), &
        "} & \np{",param_best%lat(i),"} & \np{",param_best%Tzero(i)%sec,"} \\ "
    else
      write(900,5010)" {sortie} & \np{",0.0_wr,"} & \np{",param_best%VC,"} & \np{",param_best%VM, &
        "} & \np{",param_best%Zmoho,"} & \np{",param_best%VpVs,"} & \np{",0.0_wr,"} & \np{",0.0_wr, &
        "} & \np{",0.0_wr,"} & \np{",0.0_wr,"} \\ "
    endif
    5010 format (a25,f9.2,a8,f9.2,a8,f9.2,a8,f9.2,a8,f9.3,a8,f9.2,a8,f7.4,a8,f7.4,a8,f9.2,a4)
    write(900,*)"\hline "
    write(900,*)"\end{tabular}"
    write(900,*)"\end{table}"
    write(900,*)"\end{landscape}"
    write(900,*)
    ! -----------------------------------------------------------------    .
    ! -----------------------------------------------------------------    .
    write(900,*)"\begin{longtable}[!ht]"
    write(900,*)"{|>{\centering}m{2.5cm}|>{\centering}m{1.5cm}|>{\centering}m{2.5cm}|>{\centering}m{2.5cm}", &
                "|>{\centering}m{2.5cm}|m{2.5cm}<{\centering}|}"
    write(900,*)"% Entete de la première page"
    write(900,*)"\hline \multicolumn{6}{c}{\bf r\'esidus aux stations} \\ \hline"
    write(900,*)"\hline   non station & onde & r\'esidus (s) & pond\'eration & r\'esidus / temps total ", &
                "& distance hypocentrale (km) \\ \hline"
    write(900,*)"\endfirsthead"
    write(900,*)"% Entête de toutes les pages"
    write(900,*)"\hline   non station & onde & r\'esidus (s) & pond\'eration & r\'esidus / temps total ", &
                "& distance hypocentrale (km) \\ \hline"
    write(900,*)"\endhead"
    write(900,*)"% Bas de toutes les pages"
    write(900,*)"\bf ... & ... & ... & ... & ... & ... \\ \hline"
    write(900,*)"\endfoot"
    write(900,*)"\endlastfoot"
    write(900,*)"% Contenu du tableau"
    ok=0
    j=0
    val=0.0_wr
    val2=0.0_wr
    open(901, FILE = "OUTPUT/residus"//"-"//trim(adjustl(numberfile))//".d",status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)"problème dans latexfull 1 : le fichier OUTPUT/residus"//"-"//trim(adjustl(numberfile))//".d n''existe pas "
      stop
    endif
    do while(ok .eq. 0)
      read(901,*,iostat = ok)ondetype,staname,moy,moy,residu,a_ponderation,pct, dist
      if(ok .eq. 0) then
        j=j+1
        val = val + abs(residu)
        val2 = val2 + (residu)
        if (a_ponderation.gt.0.001_wr) write(900,'(a4,a4,a2,a7,f8.4,a8,f8.4,a8,f8.2,a,f10.3,a)')staname," & \", &
        ondetype," & \np{",residu,"} & \np{",a_ponderation,"} & \np{",pct,"}\% & \np{",dist,"} \\"
      endif
    end do
    close(901)
    val = val / real(j,wr)
    val2 = val2 / real(j,wr)
    ok=0
    valec=0.0_wr
    valec2=0.0_wr
    open(904, FILE = "OUTPUT/residus"//"-"//trim(adjustl(numberfile))//".d",status='old',iostat = ok)
    if (ok .ne. 0) then
      write(*,*)"problème dans latexfull 2 : le fichier OUTPUT/residus"//"-"//trim(adjustl(numberfile))//".d n''existe pas "
      stop
    endif
    do k=1,j
      read(904,*,iostat = ok)ondetype,staname,moy,moy,residu,a_ponderation,pct,dist
      valec = valec + (abs(residu)-val)**2.0_wr
      valec2 = valec2 + (residu-val2)**2.0_wr
    end do
    close(904)
    valec = sqrt(valec/real(j,wr))
    valec2 = sqrt(valec2/real(j,wr))
    write(900,*)"\hline "
    write(900,*)"\end{longtable}"
    write(900,'(a,f18.5,a,f18.5,a)')"moyenne des r\'esidus : \np{",val2,"} $\pm$ \np{",1.0_wr*valec2,"} ($\pm1\sigma$) \\"
    write(900,'(a,f18.5,a,f18.5,a)')"moyenne des r\'esidus absolus : \np{",val,"} $\pm$ \np{",1.0_wr*valec,"} ($\pm1\sigma$) \\"

    ! -----------------------------------------------------------------    . affiche magnitude si calculée

    inquire (file="OUTPUT/files/mag-"//trim(adjustl(numberfile))//".d",exist=existe1)
    open(907, FILE = "OUTPUT/GMT/magfin-"//trim(adjustl(numberfile))//".d",status='replace')

    if ((tracessac).and.(existe1)) then

      open(906, FILE = "OUTPUT/files/mag-"//trim(adjustl(numberfile))//".d",status='old')
      write(900,*)"\vspace{1cm}"
      write(900,*)" {\Large \bf magnitude :\\}"
      write(900,*)"\vspace{1cm}"
      write(900,*)"formule de Lee et al. (1972) : "
      write(900,*)"\begin{equation}"
      write(900,*)"M_d =  -\np{0,87}  + 2 log_{10}(coda) + \np{0,0035} \Delta ,"
      write(900,*)"\end{equation}"
      write(900,'(2a)')"avec $ \Delta$, la distance \'epicentrale (en km) et {$coda$}, la dur\'ee du sigal depuis ", &
        " la premi\`ere arriv\'ee de de l'onde {\em P} jusqu'\`a la fin du signal. \\"
      write(900,*)"\vspace{1cm}"

      write(900,*)"\begin{longtable}[!ht]"
      write(900,*)"{>{\centering}m{2.5cm}>{\centering}m{2.5cm}>{\centering}m{2.5cm}m{2.5cm}<{\centering}}"
      write(900,*)"\hline"
      write(900,*)"station & M$_d$ & dur\'ee (s) & distance \'epicentrale (km)\\"
      write(900,*)"\hline"
      write(900,*)"\endfirsthead"
      write(900,*)"% Entête de toutes les pages"
      write(900,*)"\hline"
      write(900,*)"station & M$_d$ & dur\'ee (s) & distance \'epicentrale (km)\\"
      write(900,*)"\hline"
      write(900,*)"\endhead"
      write(900,*)"% Bas de toutes les pages"
      write(900,*)"\bf ... & ... & ... & ... \\ \hline"
      write(900,*)"\endfoot"
      write(900,*)"\endlastfoot"
      write(900,*)"% Contenu du tableau"
      ok=0
      j=0
      val=0.0_wr
      do while(ok .eq. 0)
        read(906,*,iostat = ok)staname,moy,tps,depi
        if (ok .eq. 0) then
          val=val+moy
          j=j+1
          write(900,'(a4,a,f10.2,a,f10.2,a,f10.2,a)')staname," & \np{",moy,"} & \np{",tps,"} & \np{",depi,"} \\"
        endif
      end do
      moy=val/real(j,wr)
      rewind(906)
      ok=0
      j=0
      ec=0.0_wr
      do while(ok .eq. 0)
        read(906,*,iostat = ok)staname,val,tps,depi
        if (ok .eq. 0) then
          ec=ec+(moy-val)**2.0_wr
          j=j+1
        endif
      end do
      close(906)
      ec = sqrt(ec/real(j,wr))
      write(900,*)"\hline"
      write(900,'(a,f10.2,a,f10.2,a)')"moyenne & \multicolumn{3}{l}{\np{",moy,"} {\small $\pm$ \np{",ec,"} (1$\sigma$)}} \\"
      ! -----------------------------------------------------------------    .
      write(907,"(i4.4,4(1a,i2.2),f10.2,f10.2,2f10.2)")dp%temps_ref(i)%date%year,"-",dp%temps_ref(i)%date%month, &
        "-",dp%temps_ref(i)%date%day,"T",dp%temps_ref(i)%date%hour,":",dp%temps_ref(i)%date%min,moy,ec, &
        refseisme(1)%mag,refseisme(2)%mag
      ! -----------------------------------------------------------------    .
      write(900,*)"\hline"
      write(900,*)"\end{longtable}"
    endif
    close(907)
    ! -----------------------------------------------------------------    . affiche les parametres de Terre si fixes
    if (FLAGterre) then
      write(900,*)"\newpage"
      write(900,*)"{\bf Param\`etres de Terre fixes} {\small (PARAM/paramTerre.d)} : \\"
      write(900,*)"\vspace{0.5cm}"
      ok=0
      open(unit=50,file="PARAM/paramTerre.d",STATUS="old",iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans latexone : le fichier PARAM/paramTerre.d n''magn pas '
        stop
      endif
      read(50,*)moy,ec
      write(900,'(a,f18.2,a,f18.3,a)')"V$_C$ (\kms) = \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
      read(50,*)moy,ec
      write(900,'(a,f18.2,a,f18.3,a)')"V$_M$ (\kms) = \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
      read(50,*)moy,ec
      write(900,'(a,f18.2,a,f18.3,a)')"Z$_{moho}$ (km) = \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
      read(50,*)moy,ec
      write(900,'(a,f18.2,a,f18.3,a)')"V$_{P}$V$_{S}$ = \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
      close(50)
    endif
    ! -----------------------------------------------------------------    . affiche les parametres hypocentaux si fixes
    if (FLAGhypo) then
      if(FLAGterre)write(900,*)"\vspace{1cm}"
      if(.not.FLAGterre) write(900,*)"\newpage"
      write(900,*)"{\bf Param\`etres hypocentaux fixes} {\small (PARAM/paramHypo.d)} : \\"
      write(900,*)"\vspace{0.5cm}"
      ok=0
      open(unit=51,file="PARAM/paramHypo.d",STATUS="old",iostat = ok)
      if (ok .ne. 0) then
        write(*,*)'problème dans latexone : le fichier PARAM/paramHypo.d n''existe pas '
        stop
      endif
      do j=1,nbseismes
        write(900,*)"s\'eisme {\tiny \#}",j,"\\"
        read(51,*)moy,ec
        write(900,'(a,f18.3,a,f18.4,a)')"longitude (\degree) = \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
        read(51,*)moy,ec
        write(900,'(a,f18.3,a,f18.4,a)')"latitude (\degree) = \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
        read(51,*)moy,ec
        write(900,'(a,f18.2,a,f18.3,a)')"Z$_{hypo}$ (km)= \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
        read(51,*)tpsref
        read(51,*)moy,ec
        call affiche_temps_ref(tpsref,char_0,1)
        write(900,'(a)')"{\NoAutoSpaceBeforeFDP"//char_0//"} \\"
        write(900,'(a,f18.2,a,f18.3,a)')"T$_{z\acute{e}ro}$ (s) = \np{",moy,"} $\pm$ \np{",2.0_wr*ec,"} ($\pm2\sigma$)\\"
        write(900,*)"\vspace{0.5cm}"
      enddo
        write(900,*)"\newpage"
      close(51)
    endif
    ! -----------------------------------------------------------------    . affiche les figures param1 VS param2
    ok = 0
    open(unit=103,file="OUTPUT/GMT/files.txt",STATUS="old",iostat = ok)
    if (ok .ne. 0) then
      write(*,*)'problème dans latexfull : le fichier OUTPUT/GMT/files.txt n''existe pas '
      stop
    endif
    do while(ok .eq. 0)
      read(103,'(i6,1x,a20)',iostat = ok)l,charname
      if (l==i) then
        if (ok .eq. 0) then
          write(900,*)"\begin{landscape} \centering "
          write(900,*)"\begin{figure}[!ht] %%%%%%%%% FIGURE _vc__vm"//"-"//trim(adjustl(numberfile))//".pdf %%%%%%%%"
          ! write(900,*)"\centerline{\includegraphics[width=1.0\linewidth]{../figures/"//charname//"}}"
            write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/"//trim(adjustl(charname))//"}}"
          write(900,*)"\end{figure}"
          write(900,*)"\end{landscape}"
        endif
      endif
    end do
    close(103)
    ! -----------------------------------------------------------------    . affiche la carte des résidus aux stations : ondes directes
    write(900,*)"\begin{landscape} \centering"
    write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%%% FIGURE residus %%%%%%%%%%"
    write(900,*)"\begin{minipage}[t]{0.45\linewidth}"
    write(900,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
      "{../figures/resPg"//"-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{minipage}"
    write(900,*)"\hfill"
    write(900,*)"\begin{minipage}[t]{0.45\linewidth}"
    write(900,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
      "{../figures/resSg"//"-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{minipage}"
    write(900,*)"\vfill"
    write(900,*)"\begin{minipage}[t]{0.45\linewidth}"
    write(900,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
      "{../figures/resPn"//"-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{minipage}"
    write(900,*)"\hfill"
    write(900,*)"\begin{minipage}[t]{0.45\linewidth}"
    write(900,*)"\centerline{\includegraphics[width=.7\linewidth,angle=-90]", &
      "{../figures/resSn"//"-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{minipage}"
    write(900,*)"\caption{\large r\'esidus aux stations (s) -- ", &
      "la taille de l'histogramme est proportionnelle aux coefficient de pond\'eration}"
    write(900,*)"\end{figure}"
    write(900,*)"\end{landscape}"
    ! -----------------------------------------------------------------    . affiche le histogrammes des résidus
    write(900,*)"\begin{landscape} \centering"
    write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%%% FIGURE reshisto %%%%%%%%%%"
    write(900,*)"\centerline{\includegraphics[width=.6\linewidth,angle=-90]", &
      "{../figures/reshisto"//"-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{figure}"
    write(900,*)"\end{landscape}"
    write(900,*)
    ! -----------------------------------------------------------------    . affiche le diagramme de Wadati
    write(900,*)"\begin{landscape} \centering"
    write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%%% FIGURE Wadati %%%%%%%%%%"
    write(900,*)"\centerline{\includegraphics[width=.6\linewidth,angle=-90]", &
      "{../figures/wadatiplot"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\caption{\large diagramme de Wadati }"
    write(900,*)"\end{figure}"
    write(900,*)"\end{landscape}"
    write(900,*)
    ! -----------------------------------------------------------------    . affiche l'hodochrone
    write(900,*)"\begin{landscape} \centering"
    write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%%%% FIGURE Hodo %%%%%%%%%%%"
    write(900,*)"\centerline{\includegraphics[width=.6\linewidth,angle=-90]", &
      "{../figures/hodochrone"//"-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\caption{\large hodochrone}"
    write(900,*)"\end{figure}"
    write(900,*)"\end{landscape}"
    write(900,*)
    ! -----------------------------------------------------------------    . affiche le plot coda/magnitude
    inquire (_FILE_DIR_="DATA/sac-"//trim(adjustl(numberfile)),exist=existe1) ! option différente selon compilo !
    if ((existe1).and.(tracessac)) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%%%% FIGURE coda %%%%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.6\linewidth,angle=-90]", &
        "{../figures/coda-"//trim(adjustl(numberfile))//".pdf}}"
      write(900,*)"\caption{\large dur\'ee de la coda et magnitude de dur\'ee, $M_d$}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
      write(900,*)
      ! ---------------------------------------------------------------    . affiche le plot coda/magnitude -> map
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%%%% FIGURE coda %%%%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.6\linewidth,angle=-90]", &
        "{../figures/coda_map-"//trim(adjustl(numberfile))//".pdf}}"
      write(900,*)"\caption{\large magnitude de dur\'ee, $M_d$, calcul\'ee \`a chaque station, pour ce s\'eisme}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
      write(900,*)
    endif
    ! -----------------------------------------------------------------    . affiche l'histogramme de la fonction coût
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE mishisto %%%%%%%%%"
      !    write(900,*)"\centerline{\includegraphics[width=1.0\linewidth]{../figures/mishisto.pdf}}"
      write(900,*)"\centerline{\includegraphics[width=.66\textwidth,angle=-90]{../figures/mishisto.pdf}}"
      write(900,*)"\caption{\large cha\^ines de Markov : fonction co\^ut des mod\`eles s\'electionn\'es}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche VC et VM en fonction des itérations
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE VC et VM %%%%%%%%%"
      !    write(900,*)"\centerline{\includegraphics[width=1.0\linewidth]{../figures/VCVM_histo.pdf}}"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/VCVM_histo.pdf}}"
      write(900,*)"\caption{\large cha\^ines de Markov : $V_c$ et $V_m$ }"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche  VpVs et Zmoho en fonction des itérations
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE VpVs Zmo %%%%%%%%%"
      !    write(900,*)"\centerline{\includegraphics[width=1.0\linewidth]{../figures/VpVsZmoho_histo.pdf}}"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/VpVsZmoho_histo.pdf}}"
      write(900,*)"\caption{\large cha\^ines de Markov : $\frac{V_P}{V_S}$ et $Z_{moho}$}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche Lat et Lon en fonction des itérations
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE lon lat  %%%%%%%%%"
      !    write(900,*)"\centerline{\includegraphics[width=1.0\linewidth]{../figures/LatLon_histo.pdf}}"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/LatLon_histo"// &
        "-"//trim(adjustl(numberfile))//".pdf}}"
      write(900,*)"\caption{\large cha\^ines de Markov : $lon$ et $lat$}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche Zhypo et Tzero en fonction des itérations
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE Zhypo To %%%%%%%%%"
      ! write(900,*)"\centerline{\includegraphics[width=1.0\linewidth]{../figures/ZhypoTzero_histo.pdf}}"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/ZhypoTzero_histo"// &
        "-"//trim(adjustl(numberfile))//".pdf}}"
      write(900,*)"\caption{\large cha\^ines de Markov : $T_{z\acute ero}$ et $Z_{hypo}$}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche fonction d'autovariance de VC et VM
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE VC et VM %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/autoVCVM_histo.pdf}}"
      write(900,*)"\caption{\large fonction d'autovariance pour $V_c$ et $V_m$ }"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche fonction d'autovariance de VpVs et Zmoho
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE VpVs Zmo %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/autoVpVsZmoho_histo.pdf}}"
      write(900,*)"\caption{\large fonction d'autovariance pour $\frac{V_P}{V_S}$ et $Z_{moho}$}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche fonction d'autovariance de Lat et Lon
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE lon lat  %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/autoLatLon_histo"// &
        "-"//trim(adjustl(numberfile))//".pdf}}"
      write(900,*)"\caption{\large fonction d'autovariance pour $lon$ et $lat$}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . affiche fonction d'autovariance de Zhypo et Tzero
    if (plotgraph) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE Zhypo To %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/autoZhypoTzero_histo"// &
        "-"//trim(adjustl(numberfile))//".pdf}}"
      write(900,*)"\caption{\large fonction d'autovariance pour $T_{z\acute ero}$ et $Z_{hypo}$}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    . matrice correlation parametres
    write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE matCorr %%%%%%%%%"
    write(900,*)"\centerline{\includegraphics[width=.9\textwidth,angle=-90]{../figures/matCorr.pdf}}"
    write(900,*)"\caption{\large matrice de corr\'elation}"
    write(900,*)"\end{figure}"
    ! -----------------------------------------------------------------    . étude a posteriori
    if (plotposteriori) then
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE post %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/postLonLat"// &
        "-"//trim(adjustl(numberfile))//".png}}"
      write(900,*)"\caption{\large \'Etude {\em a posteriori} des param\`etres Lon et Lat}"
      write(900,*)"\end{figure}"
      ! -------                                                --------    .
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE post %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/post__vc"// &
        "-"//trim(adjustl(numberfile))//".png}}"
      write(900,*)"\caption{\large \'Etude {\em a posteriori} du param\`etre $V_{c}$ }"
      write(900,*)"\end{figure}"
      ! -------                                                --------    .
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE post %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/post_vps"// &
        "-"//trim(adjustl(numberfile))//".png}}"
      write(900,*)"\caption{\large \'Etude {\em a posteriori} du param\`etre $\frac{V_P}{V_S}$ }"
      write(900,*)"\end{figure}"
      ! -------                                                --------    .
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE post %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/post__vm"// &
        "-"//trim(adjustl(numberfile))//".png}}"
      write(900,*)"\caption{\large \'Etude {\em a posteriori} du param\`etre $V_{m}$ }"
      write(900,*)"\end{figure}"
      ! -------                                                --------    .
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE post %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/post__zm"// &
        "-"//trim(adjustl(numberfile))//".png}}"
      write(900,*)"\caption{\large \'Etude {\em a posteriori} du param\`etre $Z_{moho}$ }"
      write(900,*)"\end{figure}"
      ! -------                                                --------    .
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE post %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/post__to"// &
        "-"//trim(adjustl(numberfile))//".png}}"
      write(900,*)"\caption{\large \'Etude {\em a posteriori} du param\`etre $T_{z\acute ero}$ }"
      write(900,*)"\end{figure}"
      ! -------                                                --------    .
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE post %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/post__zh"// &
        "-"//trim(adjustl(numberfile))//".png}}"
      write(900,*)"\caption{\large \'Etude {\em a posteriori} du param\`etre $Z_{hypo}$ }"
      write(900,*)"\end{figure}"
      ! -------                                                --------    .
    endif
    ! -----------------------------------------------------------------    . carrieres
    write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%% FIGURE carrieres %%%%%%%%"
    write(900,*)"\vspace{1cm}"
    write(900,*)"\centering"
    write(900,*)"\begin{minipage}{.45\textwidth}"
    write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/carrieresP1-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{minipage}"
    write(900,*)"\hfill"
    write(900,*)"\begin{minipage}{.45\textwidth}"
    write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/carrieresP2-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\end{minipage}"
    write(900,*)
    write(900,*)"\vspace{1cm}"
    write(900,*)
    write(900,*)"\centerline{\includegraphics[width=.9\textwidth]{../figures/carrieres-"//trim(adjustl(numberfile))//".pdf}}"
    write(900,*)"\caption{\large discrimination des tirs de carri\`eres}"
    write(900,*)"\end{figure}"
    write(900,*)"\vspace{1cm}"
    ! -----------------------------------------------------------------    . moho
    if (FLAG_non_tabulaire) then
      write(900,*)"\begin{landscape} \centering"
      write(900,*)"\begin{figure}[!ht]   	%%%%%%%%%%% FIGURE moho_inc %%%%%%%%%"
      write(900,*)"\centerline{\includegraphics[width=0.8\linewidth]{../figures/topo_moho"// &
        "-"//trim(adjustl(numberfile))//".pdf}}"
      write(900,*)"\caption{\large mohographie}"
      write(900,*)"\end{figure}"
      write(900,*)"\end{landscape}"
    endif
    ! -----------------------------------------------------------------    .
    write(900,*)"\end{document}"
    ! -----------------------------------------------------------------    .
    write(900,*)
    close(900)
    ! -----------------------------------------------------------------    .
    call system_clock(Nnewtime,ratetime)
    tl=(real(Nnewtime,wr)-real(Noldtime,wr))/real(ratetime,wr)
    write(*,'(a9,i2.2,'':'',i2.2,'':'',f9.2)')' temps : ',int(tl/3600.0_wr,wi), &
    int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi),(tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr)
    ! -----------------------------------------------------------------    .
  end subroutine latexone

END MODULE latexscript



! *********************************************************************    .
! *********************************************************************    .


