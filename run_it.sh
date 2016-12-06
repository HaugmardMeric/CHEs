 #####################################################################
# Méric Haugmard meric.haugmard@univ-nantes.fr
# !/bin/bash
#####################################################################
# -------                                               --------    .
# -------          CHE2013_coldruns version 2.0         --------    .
# -------          octobre 2013 - décembre 2014         --------    .
# -------                                               --------    .
#####################################################################
# version 1.1 : inversion d'un séisme (jan 2014)                    !
# version 1.2 : inversion de plusieurs séismes (juillet 2014)       !
# version 1.3 : parallelisation OpenMP (septembre 2014)             !
# version 1.4 : parallelisation MPI (octobre 2014)                  !
# version 1.5 : initialisation du prior auto (novembre 2014)        !
# version 1.6 : parallelisation full MPI (décembre 2014)            !
# version 1.7 : ajout d'une notice (décembre 2014)                  !
# version 1.8 : compilation avec ifort et gfortran (janvier 2015)   !
# version 1.9 : moho non tabulaire (fevrier 2015)                   !
# version 2.0 : test 50 séismes (septembre 2015)                    !
# version 2.1 : gestion des carrières, d'après Pascal Guterman (octobre 2015)!
# version 2.2 : ajout de modèle de terre différents pour le problème directe (novembre 2015)!
# version 2.3 : calculs a posteriori (janvier 2016)                !
#####################################################################
# compile et execute le programme CHE
# permet aussi l'écriture du scripte LOG
#####################################################################
# The default process manager is called MPD,
# which is a ring of daemons on the machines
# where you will run your MPI programs.
# mpd &
#####################################################################


T="$(date +%s)"

chmod +x SRC/run.sh

./SRC/run.sh 2> >(tee stderrlog.d | tee -a alllog.d > /dev/tty ) | tee stdoutlog.d | tee -a alllog.d 

mv *log.d OUTPUT/LOG
T="$(($(date +%s)-T))"
echo "execution time (secs) ${T}"


rm -rf .gmtcommands4 .gmtdefaults4

#####################################################################
