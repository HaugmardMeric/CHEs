# ------------------------------------------------------ #
#	Makefile programme CHE (cold-hot-plot) 
# ------------------------------------------------------ #
# Makefile produit à partir de makedepf90 (version 2.8.8 -- Erik.Edelmann@iki.fi)
# -> makedepf90 *.f90 '-r $(COMPIL) $(FFLAGS) -c $<' -o aout > makefile
# ------------------------------------------------------ #

# -------------------- Compilos ------------------------ # test gfortran ou ifort

ifeq ($(shell if ifort -v 2> /dev/null ; then echo "OKifort" ; fi),OKifort)
    include macros/ifort.d
else
    ifeq ($(shell if gfortran -v 2> /dev/null ; then echo "OKgfortran" ; fi),OKgfortran)
        include macros/gfortran.d
    else

stop:
		@echo "pas de compilateur fortran !"
		exit 0
    endif
endif

# --------------- Dependent o-files -------------------- #

FOBJ1=lib/modparam.o lib/mt19937ar.o lib/types.o lib/avancement.o lib/time.o lib/dellipsgc.o lib/subgeiger.o lib/McMC.o lib/intersect.o lib/tirage.o lib/tracer_rais.o lib/pbdirect.o lib/misfit.o lib/lectdata.o lib/printmess.o lib/stat.o lib/subparam.o lib/tri.o lib/rechercheinit.o

FOBJ2=lib/mkparamiter.o lib/mkmatricecorrel.o lib/mkfcout.o lib/mkhodo.o lib/mkcoda.o lib/mklatex.o lib/mkmap.o lib/mkres.o lib/mkwada.o lib/mkchatelain.o lib/mkmoho_inc.o lib/mkGMT.o lib/mkcarriere.o lib/mkposteriori.o

FOBJ3=lib/mod_sac_io.o lib/modparam.o lib/types.o lib/time.o lib/stat.o

FOBJ4=lib/modparam.o lib/mt19937ar.o lib/types.o lib/avancement.o lib/time.o lib/dellipsgc.o lib/subgeiger.o lib/tracer_rais.o lib/pbdirect.o lib/misfit.o lib/lectdata.o lib/stat.o lib/subparam.o lib/tri.o lib/rechercheinit.o

FOBJ5=lib/mod_sac_io.o lib/modparam.o lib/types.o lib/time.o lib/stat.o lib/lectdata.o lib/pbdirect.o lib/dellipsgc.o lib/mt19937ar.o lib/tracer_rais.o lib/misfit.o lib/tri.o


# --------------------- makes -------------------------- #

all: rmv1 mess1 coldruns hotruns rmv2 mess2

coldruns: che_coldruns_init che_coldruns che_coldruns_syn
hotruns: che_hotruns_init che_hotruns che_hotruns_syn
plot: che_plot
apriori: rmv1 che_apriori
others: rmv1 sac_bin2txt sac_coda sac_readpick sac_spectre sac_stalta_kurtosis verifmediatrice sac_ZNE_2_LQT
otherssac: sac_writepick sac_writepickTheo sac_writepickCata

# ----------------- main programs ---------------------- #

# ----------------- coldruns

che_coldruns_init: che_coldruns_init.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/che_coldruns_init.o $(FOBJ1) lib/mkposteriori.o
	mv che_coldruns_init ../BIN/che_coldruns_init.exe

che_coldruns: che_coldruns.o
	$(COMPILMPI) -o $@ $(OPTIONC) $(FFLAGS) lib/che_coldruns.o $(FOBJ1)
	mv che_coldruns ../BIN/che_coldruns.exe

che_coldruns_syn: che_coldruns_syn.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/che_coldruns_syn.o $(FOBJ1)
	mv che_coldruns_syn ../BIN/che_coldruns_syn.exe

# ----------------- hotruns

che_hotruns_init: che_hotruns_init.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/che_hotruns_init.o $(FOBJ1)
	mv che_hotruns_init ../BIN/che_hotruns_init.exe

che_hotruns: che_hotruns.o
	$(COMPILMPI) -o $@ $(OPTIONC) $(FFLAGS) lib/che_hotruns.o $(FOBJ1)
	mv che_hotruns ../BIN/che_hotruns.exe

che_hotruns_syn: che_hotruns_syn.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/che_hotruns_syn.o $(FOBJ1)
	mv che_hotruns_syn ../BIN/che_hotruns_syn.exe

# ----------------- plot

che_plot: che_plot.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/che_plot.o $(FOBJ1) $(FOBJ2)
	mv che_plot ../BIN/che_plot.exe

# ----------------- apriori

che_apriori: che_apriori.o
	$(COMPILMPI) -o $@ $(OPTIONC) $(FFLAGS) lib/che_apriori.o $(FOBJ1)
	mv che_apriori ../BIN/che_apriori.exe

# ----------------- others

sac_bin2txt: sac_bin2txt.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_bin2txt.o $(FOBJ3)
	mv sac_bin2txt ../BIN/sac_bin2txt.exe

sac_coda: sac_coda.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_coda.o $(FOBJ3)
	mv sac_coda ../BIN/sac_coda.exe

sac_readpick: sac_readpick.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_readpick.o $(FOBJ3)
	mv sac_readpick ../BIN/sac_readpick.exe

sac_spectre: sac_spectre.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_spectre.o $(FOBJ3)
	mv sac_spectre ../BIN/sac_spectre.exe

sac_ZNE_2_LQT: sac_ZNE_2_LQT.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_ZNE_2_LQT.o $(FOBJ3) lib/dellipsgc.o
	mv sac_ZNE_2_LQT ../BIN/sac_ZNE_2_LQT.exe

sac_stalta_kurtosis: sac_stalta_kurtosis.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_stalta_kurtosis.o $(FOBJ3)
	mv sac_stalta_kurtosis ../BIN/sac_stalta_kurtosis.exe

verifmediatrice : verifmediatrice.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/verifmediatrice.o $(FOBJ4)
	mv verifmediatrice ../BIN/verifmediatrice.exe

# ----------------- otherssac

sac_writepick: sac_writepick.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_writepick.o $(FOBJ3)
	mv sac_writepick ../BIN/sac_writepick.exe

sac_writepickTheo: sac_writepickTheo.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_writepickTheo.o $(FOBJ3)
	mv sac_writepickTheo ../BIN/sac_writepickTheo.exe

sac_writepickCata: sac_writepickCata.o
	$(COMPIL) -o $@ $(OPTIONC) $(FFLAGS) lib/sac_writepickCata.o $(FOBJ5)
	mv sac_writepickCata ../BIN/sac_writepickCata.exe

# ------------------------------------------------------ #

clean : rmv1

clear : rmv1

# --------------------- rules -------------------------- # 

rmv2 :
	clear ; rm -rf lib/*.mod lib/*.o 

rmv1 :
	clear ; rm -rf lib/*.mod lib/*.o ../BIN/*.exe

mess1 :
	cat < MES/message1.d

mess2 :
	cat < MES/message2.d

# ------------ rules dependent o-files ----------------- #

# --------------------- programmes --------------------- #

che_coldruns.o : PROG/coldruns/che_coldruns.f90 misfit.o subparam.o printmess.o time.o pbdirect.o tirage.o McMC.o avancement.o mt19937ar.o types.o modparam.o
	$(COMPILMPI) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
che_coldruns_init.o : PROG/coldruns/che_coldruns_init.f90 rechercheinit.o subparam.o printmess.o lectdata.o mt19937ar.o types.o modparam.o McMC.o intersect.o tirage.o mkposteriori.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
che_coldruns_syn.o : PROG/coldruns/che_coldruns_syn.f90 tri.o subparam.o lectdata.o printmess.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
che_hotruns.o : PROG/hotruns/che_hotruns.f90 dellipsgc.o misfit.o subparam.o printmess.o pbdirect.o lectdata.o tirage.o McMC.o avancement.o mt19937ar.o time.o types.o modparam.o
	$(COMPILMPI) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
che_hotruns_init.o : PROG/hotruns/che_hotruns_init.f90 subparam.o printmess.o time.o tirage.o mt19937ar.o types.o modparam.o McMC.o misfit.o lectdata.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
che_hotruns_syn.o : PROG/hotruns/che_hotruns_syn.f90 subparam.o printmess.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
che_plot.o : PROG/che_plot.f90 mt19937ar.o mkGMT.o subparam.o printmess.o lectdata.o mklatex.o types.o modparam.o McMC.o intersect.o tirage.o mkposteriori.o time.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
che_apriori.o : PROG/che_apriori.f90 printmess.o mt19937ar.o subparam.o types.o modparam.o McMC.o intersect.o tirage.o misfit.o lectdata.o
	$(COMPILMPI) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_bin2txt.o : PROG/sac_bin2txt.f90 time.o types.o modparam.o mod_sac_io.o stat.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_coda.o : PROG/sac_coda.f90 mod_sac_io.o stat.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_readpick.o : PROG/sac_readpick.f90 mod_sac_io.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_stalta_kurtosis.o : PROG/sac_stalta_kurtosis.f90 mod_sac_io.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
verifmediatrice.o : PROG/verifmediatrice.f90 rechercheinit.o lectdata.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_writepick.o : PROG/sac_writepick.f90 mod_sac_io.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_spectre.o : PROG/sac_spectre.f90 mod_sac_io.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_writepickTheo.o : PROG/sac_writepickTheo.f90 mod_sac_io.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_writepickCata.o : PROG/sac_writepickCata.f90 mod_sac_io.o time.o types.o modparam.o lectdata.o pbdirect.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib
sac_ZNE_2_LQT.o : PROG/sac_ZNE_2_LQT.f90 mod_sac_io.o modparam.o dellipsgc.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o lib

# ---------------    for modules    -------------------- #

McMC.o : MOD/McMC.f90 mt19937ar.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
avancement.o : MOD/avancement.f90 modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
modparam.o : MOD/modparam.f90
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
dellipsgc.o : MOD/dellipsgc.f90 modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
intersect.o : MOD/intersect.f90 modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
lectdata.o : MOD/lectdata.f90 dellipsgc.o misfit.o subparam.o pbdirect.o mt19937ar.o tri.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
misfit.o : MOD/misfit.f90 pbdirect.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkGMT.o : MOD/MOD_GMT/mkGMT.f90 tri.o lectdata.o mkmoho_inc.o mkcoda.o mkhodo.o mkmap.o mkres.o mkfcout.o mkmatricecorrel.o mkparamiter.o mkchatelain.o mkwada.o misfit.o dellipsgc.o subparam.o pbdirect.o time.o types.o modparam.o mkcarriere.o mkposteriori.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkchatelain.o : MOD/MOD_GMT/mkchatelain.f90 pbdirect.o time.o avancement.o types.o modparam.o tri.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkcoda.o : MOD/MOD_GMT/mkcoda.f90 pbdirect.o stat.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkfcout.o : MOD/MOD_GMT/mkfcout.f90 mkparamiter.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkhodo.o : MOD/MOD_GMT/mkhodo.f90 pbdirect.o stat.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mklatex.o : MOD/MOD_LaTeX/mklatex.f90 pbdirect.o misfit.o subgeiger.o mkGMT.o subparam.o lectdata.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkmap.o : MOD/MOD_GMT/mkmap.f90 types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkmatricecorrel.o : MOD/MOD_GMT/mkmatricecorrel.f90 stat.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkmoho_inc.o : MOD/MOD_GMT/mkmoho_inc.f90 pbdirect.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkparamiter.o : MOD/MOD_GMT/mkparamiter.f90 types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkres.o : MOD/MOD_GMT/mkres.f90 stat.o subparam.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkposteriori.o : MOD/MOD_GMT/mkposteriori.f90 types.o modparam.o dellipsgc.o tri.o avancement.o tirage.o pbdirect.o misfit.o rechercheinit.o subparam.o time.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkcarriere.o : MOD/MOD_GMT/mkcarriere.f90 modparam.o types.o dellipsgc.o time.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mkwada.o : MOD/MOD_GMT/mkwada.f90 pbdirect.o time.o avancement.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mod_sac_io.o : MOD/MOD_sac/mod_sac_io.f90
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
mt19937ar.o : MOD/MOD_rand/mt19937ar.f90 modparam.o
	$(COMPIL) -c -Ilib $(FFLAGS) $(OPTIOND) $< ; mv *.o *.mod lib
pbdirect.o : MOD/pbdirect.f90 stat.o time.o dellipsgc.o types.o modparam.o tracer_rais.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
printmess.o : MOD/printmess.f90 modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
rechercheinit.o : MOD/rechercheinit.f90 dellipsgc.o mt19937ar.o time.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
stat.o : MOD/stat.f90 modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
tracer_rais.o : MOD/tracer_rais.f90 modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
subgeiger.o : MOD/MOD_Geiger/subgeiger.f90 dellipsgc.o pbdirect.o time.o mt19937ar.o types.o modparam.o misfit.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
subparam.o : MOD/subparam.f90 rechercheinit.o subgeiger.o mt19937ar.o tri.o avancement.o stat.o time.o types.o modparam.o tri.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
time.o : MOD/time.f90 types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
tirage.o : MOD/tirage.f90 time.o intersect.o mt19937ar.o types.o modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
tri.o : MOD/tri.f90 time.o types.o modparam.o mt19937ar.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib
types.o : MOD/types.f90 modparam.o
	$(COMPIL) -c -Ilib $(OPTIONC) $(FFLAGS) $< ; mv *.o *.mod lib

# ------------------------------------------------------ #
