#####################################################################
# Méric Haugmard meric.haugmard@univ-nantes.fr
# !/bin/bash
#####################################################################
# -------                                               --------    .
# -------          CHE2013_coldruns version 2.0         --------    .
# -------          octobre 2013 - décembre 2014         --------    .
# -------                                               --------    .
#####################################################################
# mpiexec <- pour executer mpi
# mpd &
#####################################################################

clear
#####################################################################
cd SRC
make all
cd ..

#####################################################################
########## # permet la sauvegarde du dernier repertoire OUTPUT
########## # par défaut : no
########## ans="No"
########## # attend 3 secondes une réponse positive
########## read -p "Sauver dernier run [y/n]?" -t 3 ans
########## if [ $ans = y -o $ans = Y -o $ans = yes -o $ans = Yes -o $ans = YES ]
########## then
##########     # on creer un repertoire et sauve dernier run dans OLD
##########     echo yes !
##########     LA_DATE=$(date +%Y"_"%m"_"%d"_"%H"h"%M"_"%S)
##########     mkdir OLD/$LA_DATE ; mv OUTPUT OLD/$LA_DATE/OUTPUT
########## else
##########     # on supprime le repertoire OUTPUT sans sauver le dernier run
##########     echo no !
rm -rf OUTPUT
########## fi
#####################################################################
# (re)création de l'arborescence
mkdir OUTPUT
mkdir OUTPUT/figures
mkdir OUTPUT/files
mkdir OUTPUT/files/Cold
mkdir OUTPUT/files/Hot
mkdir OUTPUT/files/STA
mkdir OUTPUT/files/Plot
mkdir OUTPUT/LOG
mkdir OUTPUT/input
mkdir OUTPUT/GMT
mkdir OUTPUT/LATEX
#####################################################################
cd DATA
ls -f *.dat > seismes.d 2>/dev/null
cd ..

#####################################################################
head -1 PARAM/iteration.d > toto.d
read nbchainecold itercold < toto.d
tail -1 PARAM/iteration.d > toto.d
read nbchainehot iterhot < toto.d
rm -rf toto.d
#####################################################################
# programmes principaux
# FORT_FMT_RECL=1000 -> permet d'écrire des fichiers textes de plus de 1000 caracteres par lignes (pour ifort)
#####################################################################
# coldruns
FORT_FMT_RECL=1000 ./BIN/che_coldruns_init.exe || exit  # exécute {exit} uniquement si {/BIN/che_coldruns_init.exe} échoue
    #################################################################
read nbseisme < OUTPUT/GMT/nbseisme.d
echo '/dev/null' > cmd.exe
    #################################################################
FORT_FMT_RECL=1000 mpiexec -n $nbchainecold ./BIN/che_coldruns.exe < cmd.exe || exit
FORT_FMT_RECL=1000 ./BIN/che_coldruns_syn.exe || exit
    #################################################################
# hotruns
FORT_FMT_RECL=1000 ./BIN/che_hotruns_init.exe || exit
FORT_FMT_RECL=1000 mpiexec -n $nbchainehot ./BIN/che_hotruns.exe < cmd.exe || exit
FORT_FMT_RECL=1000 ./BIN/che_hotruns_syn.exe || exit
    #################################################################
# plots
FORT_FMT_RECL=1000 ./BIN/che_plot.exe || exit
FORT_FMT_RECL=1000 mpiexec -n $nbseisme ./BIN/che_apriori.exe < cmd.exe || exit
#####################################################################
rm -rf .gmtcommands4 .gmtdefaults4 cmd.exe
# supprime les anciennes options par défaut
chmod +x OUTPUT/GMT/script0.sh
# figure recherche_initiale
./OUTPUT/GMT/script0.sh || exit
#####################################################################
# diverses copies
cp PARAM/priorIn_HOT.d OUTPUT/input/priorIn_HOT.d
cp PARAM/priorIn_COLD.d OUTPUT/input/priorIn_COLD.d
cp PARAM/paramHypo.d OUTPUT/input/paramHypo.d 2>/dev/null
cp PARAM/paramTerre.d OUTPUT/input/paramTerre.d 2>/dev/null
cp PARAM/iteration.d OUTPUT/input/iteration.d
cp DATA/*.d OUTPUT/input/
#####################################################################
# execution des scripts GMT
chmod +x OUTPUT/GMT/script.sh
./OUTPUT/GMT/script.sh || exit
#####################################################################
# execution des scripts LaTeX
rm -rf OUTPUT/LOG/gslog.d
cd OUTPUT/LATEX/
grep '*' 2*.tex
ls 2*.tex sta*.tex 2>/dev/null | while read afile
do
  echo $afile
  pdflatex $afile >> ../LOG/afilelog.d
  pdflatex $afile >> ../LOG/afilelog.d
  afilepdf=${afile/tex/pdf}
  gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -sOutputFile=../$afilepdf $afilepdf >> ../LOG/gslog.d
done
cd ../..
#####################################################################

if [ -d DOC/SCRIPT ]
then
    #################################################################
    rm -rf ./DOC/SCRIPT/*.aux ./DOC/SCRIPT/*.pdf ./DOC/SCRIPT/*.log ./DOC/SCRIPT/*.gz
    rm -rf ./DOC/SCRIPT/prog.txt ./DOC/SCRIPT/pbashacc.txt
    #################################################################
    # execution d'autres scripts
    chmod +x DOC/SCRIPT/makesumfiles.sh
    ./DOC/SCRIPT/makesumfiles.sh
    echo EDITscripts.tex
    cd DOC/SCRIPT/
    pdflatex EDITscripts.tex > afilelog.d || exit
    pdflatex EDITscripts.tex > afilelog.d || exit
    cd ../..
    #################################################################
fi
#####################################################################


