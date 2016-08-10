#!/bin/sh
# meric haugmard, avril 2015
#############################################################
########## 			ppk 		    #########
############################################################# à modifié dans mk_ppk.sh
# pointe sur le accelerogrammes :
pointe=ACC
# sur les séisme num :
seisme=1
############################################################# pick les ondes
#  Pg (T1-2) Pn (T5-6) Sg (T3-4) Sn (T7-8)
# --- #
# T1 : temps arrivée de l'onde Pg
# T2 : incertitudes absolue sur l'onde Pg -> calcul du delta après
# --- #
# T3 : temps arrivée de l'onde Sg
# T4 : incertitudes absolue sur l'onde Sg
# --- #
# T5 : temps arrivée de l'onde Pn
# T6 : incertitudes absolue sur l'onde Pn
# --- #
# T7 : temps arrivée de l'onde Sn
# T8 : incertitudes absolue sur l'onde Sn
# --- #
#  coef : (bon) 0 1 2 3 (mauvais)
#  coef : 4 pas pris en compte
#  écrire : H
#  sauver : N
#  zoom X, retour O
#############################################################

./mk_ppk.sh

date > data.d
date > err.log.d
#############################################################
if test -f  datafile.d  ; then
    cat  datafile.d  | cut -c 1-6 | sed '1d' | sed -e 's/JSA0/JSA-/g'  | sed -e 's/DYA0/DYA-/g'  | sed -e 's/HTL0/HTL-/g'  | sed -e 's/MFF0/MFF-/g'  | sed -e 's/PG/ /g' | sed -e 's/PN/ /g' | sort | uniq | sed '/^$/d' > toto.deux

    while read sta
    do

	if [ ${sta:3:1} = '-' ]
	then 
		sta2=${sta:0:3}
	else
		sta2=$sta
	fi


        if test -f sac$seisme/*$sta2*HZ*$pointe ; then

            if test -f sac$seisme/*$sta2*HN*$pointe ; then
                ./../BIN/sac_stalta_kurtosis.exe sac$seisme/*$sta2*HZ*$pointe sac$seisme/*$sta2*HN*$pointe sac$seisme/*$sta2*HE*$pointe
            else
                ./../BIN/sac_stalta_kurtosis.exe sac$seisme/*$sta2*HZ*$pointe
            fi

sac << EOF >> log.d 2>&1 | tee -a err.log.d
qdp off
read sac$seisme/*$sta2*$pointe sac$seisme/*$sta2*kurt
rmean 
rtrend
taper
SYNCHRONIZE BEGIN ON
ch T9 -12345
wh
PLOTPK ABSOLU REFERENCE OFF MARKALL ON SAVELOCS ON BELL OFF
w h
w over
quit
EOF

    		echo sac$seisme/*$sta2*$pointe
    		./../BIN/sac_readpick.exe sac$seisme/*$sta2*HZ*$pointe >> data.d
        fi
    done < toto.deux
fi
#############################################################
rm -rf toto.deux
cat err.log.d
echo
cat data.d > datafile.d
#############################################################
