#!/bin/sh
# meric haugmard, avril 2015
#############################################################
##########             ppk             ######### 
############################################################# à modifié dans ppk.sh
# pointe sur le accelerogrammes :
pointe=ACC
# sur les séisme num :
seisme=1
#############################################################
cat *-$seisme.dat > datafile.d
############################################################# recherches sous toutes les stations

    for sta in PY91 PY92 PY94 PY95 PY96 PY97 PY41 PY42 PY43 PY44 PY45 PY46 PY47 PY48 PY49 PY4A PY40 PY33 PY26 PY34 CCA1 ELSH HMNX SWN1 RENF CHIF LRVF ROT2 ROSF DYA- HTL- JSA- MFF-
    do 

	if [ ${sta:3:1} = '-' ]
	then 
		sta2=${sta:0:3}
	else
		sta2=$sta
	fi

	rm -rf toto.d
	if test -f datafile.d ; then
	  grep $sta2 datafile.d | wc -l > toto.d
	  read avalue <  toto.d
	else
	  avalue=-999
	fi
	if [[ $avalue -eq 0 ]]; then
          echo ${sta}PG >> datafile.d
	fi
    done
    rm -rf toto.d


############################################################# temps précédemment pointés

if test -f datafile.d ; then
    cat  datafile.d | cut -c 1-4  | sed '1d'  | sed -e 's/JSA0/JSA-/g'  | sed -e 's/DYA0/DYA-/g'  | sed -e 's/HTL0/HTL-/g'  | sed -e 's/MFF0/MFF-/g' | sort | uniq | sed '/^$/d' > toto.un
    while read sta
    do
	if [ ${sta:3:1} = '-' ]
	then 
		sta2=${sta:0:3}
		sta=${sta:0:3}'0'
	else
		sta2=$sta
	fi

        if test -f sac$seisme/*$sta2*HZ*$pointe ; then
            ls sac$seisme/*$sta2*HZ*$pointe | while read sacfile
            do
                grep $sta datafile.d > picks.d
                head -n 2 picks.d > apick.d
                ./../BIN/sac_writepick.exe $sacfile apick.d
                tail -n 1 picks.d > apick.d
                ./../BIN/sac_writepick.exe $sacfile apick.d
                rm -rf picks.d apick.d
            done
        fi
    done < toto.un
fi
rm -rf toto.un

############################################################# temps de dernier meilleur modèle : TA, TO
# rm -rf ../OUTPUT/files/tempsTheoOUT_1.d

if test -f ../OUTPUT/files/tempsTheoOUT_1.d ; then
    cat ../OUTPUT/files/tempsTheoOUT_1.d | cut -c 1-4  | sed '1d'  | sed -e 's/JSA0/JSA-/g'  | sed -e 's/DYA0/DYA-/g'  | sed -e 's/HTL0/HTL-/g'  | sed -e 's/MFF0/MFF-/g'  | sort | uniq  > toto.un
    while read sta
    do


	if [ ${sta:3:1} = '-' ]
	then 
		sta2=${sta:0:3}
		sta=${sta:0:3}'0'
	else
		sta2=$sta
	fi

        if test -f sac$seisme/*$sta2*HZ*$pointe ; then
            ls sac$seisme/*$sta2*$pointe | while read sacfile
            do
                grep $sta ../OUTPUT/files/tempsTheoOUT_1.d > picks.d

                head -n 2 picks.d > apick.d
                ./../BIN/sac_writepickTheo.exe $sacfile apick.d
                tail -n 1 picks.d > apick.d
                ./../BIN/sac_writepickTheo.exe $sacfile apick.d
                rm -rf picks.d apick.d
            done
        fi
    done < toto.un 
    rm -rf toto.un
else
  cat  datafile.d | cut -c 1-6  | sed '1d' | sed -e 's/PG/ /g' | sed -e 's/PN/ /g' | sort | uniq | sed '/^$/d' > toto.un
  while read sta
  do

	if [ ${sta:3:1} = '-' ]
	then 
		sta2=${sta:0:3}
	else
		sta2=$sta
	fi

        if test -f sac$seisme/*$sta2*HZ*$pointe ; then
            ls sac$seisme/*$sta2*$pointe | while read sacfile
            do
                ./../BIN/sac_writepickCata.exe $sacfile
            done
        fi
    done < toto.un 
    rm -rf toto.un

fi    

#############################################################
