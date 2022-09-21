#
#
# makefile for gints
#
#
#
export ON_DARTER   = 0
ifeq (${ON_DARTER}, 1)
include compile.cfg.darter
else
include compile.cfg
endif
NAME                 = SCF_DERIV
OBJC                 = ${BIN}/parseshell.o    ${BIN}/shellprop.o    ${BIN}/shell.o \
							  ${BIN}/purecart.o      ${BIN}/shellsize.o    ${BIN}/cptrans.o \
							  ${BIN}/excep.o         ${BIN}/globalinfor.o  ${BIN}/xcfunc.o \
							  ${BIN}/filerw.o        ${BIN}/textread.o     ${BIN}/parameterparsing.o \
							  ${BIN}/localmemscr.o   ${BIN}/blas.o         ${BIN}/lapack.o \
							  ${BIN}/matrix.o        ${BIN}/in_place_transpose.o \
							  ${BIN}/blockmatrix.o   ${BIN}/gints2d.o \
							  ${BIN}/element.o       ${BIN}/molecule.o  \
							  ${BIN}/denmtrx.o       ${BIN}/mo.o           ${BIN}/denmtrxinfor.o \
							  ${BIN}/sigshellpairinfor.o  ${BIN}/gints4deriv.o\
							  ${BIN}/oneintformula.o ${BIN}/shellpair.o \
							  ${BIN}/gintsinfor.o    ${BIN}/integralinfor.o ${BIN}/symmjkdigest.o     \
							  ${BIN}/histdataman.o   ${BIN}/oneemtrx.o \
							  ${BIN}/digestutil.o    ${BIN}/integraljobs.o \
							  ${BIN}/scfparam.o      ${BIN}/xcintsinfor.o \
							  ${BIN}/scfderiv.o      ${BIN}/main.o

ALL: PREPARE $(OBJC) ${LIBERI_D1} ${LIBERI_D2} ${LIBTWOBODYOVERLAP} ${LIBKINETIC} ${LIBNAI} ${LIBERI}
	$(CC) -o $(NAME) $(OBJC) $(LIBRARY) ${LIBERI_D1} ${LIBERI_D2} ${LIBTWOBODYOVERLAP} ${LIBKINETIC} ${LIBNAI} ${LIBERI}

PREPARE:
	@(if [ ! -d "${BIN}" ]; then \
		${MKDIR}  ${BIN} ;\
	fi)

# liberi_d1.a
${LIBERI_D1}:
	${MAKE} -f makefiles_intengine/Makefile_eri_d1

# liberi_d2.a
${LIBERI_D2}:
	${MAKE} -f makefiles_intengine/Makefile_eri_d2

# libtwobodyoverlap.a
${LIBTWOBODYOVERLAP}:
	${MAKE} -f makefiles_intengine/Makefile_twobodyoverlap

# libkinetic.a
${LIBKINETIC}:
	${MAKE} -f makefiles_intengine/Makefile_kinetic

# libtwobodyoverlap.a
${LIBNAI}:
	${MAKE} -f makefiles_intengine/Makefile_nai

# liberi.a
${LIBERI}:
	${MAKE} -f makefiles_intengine/Makefile_eri

###########################################################
#    here below are normal binary files
###########################################################

# this folder
${BIN}/main.o: main.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# gints module
${BIN}/digestutil.o: ../../gints/src/digestutil.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gints2d.o: ../../gints/src/gints2d.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gints4d.o: ../../gints/src/gints4d.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gints4deriv.o: ../../gints/src/gints4deriv.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/sigshellpairinfor.o: ../../gints/src/sigshellpairinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/shellpair.o: ../../gints/src/shellpair.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/oneintformula.o: ../../gints/src/oneintformula.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gintsinfor.o: ../../gints/src/gintsinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/exrho.o: ../../gints/src/exrho.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/symmjkdigest.o: ../../gints/src/symmjkdigest.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/integralinfor.o: ../../gints/src/integralinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/integraljobs.o: ../../gints/src/integraljobs.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/oneemtrx.o: ../../gints/src/oneemtrx.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# shell module
${BIN}/parseshell.o: ../../shell/src/parseshell.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/shellprop.o: ../../shell/src/shellprop.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/shell.o: ../../shell/src/shell.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/purecart.o: ../../shell/src/purecart.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/cptrans.o: ../../shell/src/cptrans.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/shellsize.o: ../../shell/src/shellsize.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# general module
${BIN}/globalinfor.o: ../../general/src/globalinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/textread.o: ../../general/src/textread.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/xcfunc.o: ../../general/src/xcfunc.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/excep.o: ../../general/src/excep.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/filerw.o: ../../general/src/filerw.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/localmemscr.o: ../../general/src/localmemscr.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/parameterparsing.o: ../../general/src/parameterparsing.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/histdataman.o: ../../general/src/histdataman.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# molecule module
${BIN}/element.o: ../../geom/src/element.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/molecule.o: ../../geom/src/molecule.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# result module
${BIN}/denmtrx.o: ../../result/src/denmtrx.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/denmtrxinfor.o: ../../result/src/denmtrxinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/mo.o: ../../result/src/mo.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# math module
${BIN}/blas.o: ../../math/src/blas.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/lapack.o: ../../math/src/lapack.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/matrix.o: ../../math/src/matrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/blockmatrix.o: ../../math/src/blockmatrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/blockmatrixlist.o: ../../math/src/blockmatrixlist.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/in_place_transpose.o: ../../math/src/in_place_transpose.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# scf_deriv module
${BIN}/scfderiv.o: ../../scf_deriv/src/scfderiv.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# other modules
${BIN}/xcintsinfor.o: ../../xcints/src/xcintsinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/scfparam.o: ../../scf/src/scfparam.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@


.PHONY: clean
clean:
	@rm -rf ${BIN}
	@rm -rf test

