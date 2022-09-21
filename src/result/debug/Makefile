#
#
# makefile for scf folder
#
#
#
include compile.cfg
NAME                 = test

# now list all of OBJs
OBJ_RESULT           = ${BIN}/denmtrx.o       ${BIN}/mo.o  \
							  ${BIN}/atomdenmtrx.o
OBJ_GINTS            = ${BIN}/digest.o        ${BIN}/sigshellpairinfor.o \
							  ${BIN}/oneintformula.o ${BIN}/shellpair.o \
							  ${BIN}/kdigestinfor.o  ${BIN}/gintsinfor.o \
							  ${BIN}/gints2d.o       ${BIN}/gints4d.o    \
							  ${BIN}/kdigest.o       ${BIN}/digestkints.o  \
							  ${BIN}/jdigest.o       ${BIN}/exrho.o    \
							  ${BIN}/oneemtrx.o
OBJ_XCINTS           = ${BIN}/xcintsinfor.o   ${BIN}/xcvarinfor.o  \
							  ${BIN}/xcvar.o         ${BIN}/atomgrids.o \
							  ${BIN}/molegrids.o     ${BIN}/lebedev.o \
							  ${BIN}/batchgrid.o     ${BIN}/sigatombasis.o \
							  ${BIN}/dftbasisderiv1.o ${BIN}/dftbasisderiv2.o \
							  ${BIN}/dftbasisderiv3.o ${BIN}/dftbasisderiv4.o \
							  ${BIN}/batchbasis.o    ${BIN}/dftmatrix.o  \
							  ${BIN}/batchvar.o      ${BIN}/batchfunc.o  \
							  ${BIN}/batchxcmtrx.o   ${BIN}/xcints.o
OBJ_XCFUNC           = ${BIN}/initfuncderives.o  ${BIN}/becke88.o   \
							  ${BIN}/lyp.o              ${BIN}/pbex.o   \
							  ${BIN}/pw86x.o            ${BIN}/pw91x.o  \
							  ${BIN}/s.o                ${BIN}/vwn1rpa.o  \
							  ${BIN}/pbec.o             ${BIN}/pw86c.o \
							  ${BIN}/pw91c.o            ${BIN}/pw92c.o \
							  ${BIN}/vdwbr89.o          ${BIN}/vwn5.o  \
							  ${BIN}/becke05ex.o        ${BIN}/br94coor.o  \
							  ${BIN}/br89exhole.o
OBJ_SHELL            = ${BIN}/parseshell.o    ${BIN}/shellprop.o    \
							  ${BIN}/shell.o         ${BIN}/purecart.o   \
						  	  ${BIN}/shellsize.o     ${BIN}/cptrans.o 
OBJ_GENERAL          = ${BIN}/excep.o         ${BIN}/globalinfor.o  \
							  ${BIN}/xcfunc.o        ${BIN}/filerw.o      \
							  ${BIN}/textread.o      ${BIN}/parameterparsing.o \
							  ${BIN}/localmemscr.o   ${BIN}/histdataman.o
OBJ_GEOM             = ${BIN}/element.o       ${BIN}/molecule.o  
OBJ_MATH             = ${BIN}/blas.o          ${BIN}/lapack.o \
							  ${BIN}/matrix.o        ${BIN}/in_place_transpose.o \
							  ${BIN}/blockmatrix.o   ${BIN}/blockmatrixlist.o  
OBJ_SCF              = ${BIN}/scf.o           ${BIN}/scfparam.o  \
							  ${BIN}/fock.o          ${BIN}/scfdiis.o   \
							  ${BIN}/scfeadiis.o     ${BIN}/gslscfminimizer.o \
							  ${BIN}/scfconv.o       ${BIN}/scfintscontroller.o 
OBJ_TEST             = ${BIN}/coreguess_test.o ${BIN}/atomdenmtrx_test.o \
							  ${BIN}/main.o
OBJC                 = ${OBJ_MATH} 
OBJC                += ${OBJ_GEOM} 
OBJC                += ${OBJ_GENERAL} 
OBJC                += ${OBJ_SHELL} 
OBJC                += ${OBJ_GINTS} 
OBJC                += ${OBJ_XCINTS} 
OBJC                += ${OBJ_XCFUNC} 
OBJC                += ${OBJ_RESULT} 
OBJC                += ${OBJ_SCF}
OBJC                += ${OBJ_TEST}

ALL: PREPARE $(OBJC) ${LIBTWOBODYOVERLAP} ${LIBKINETIC} ${LIBNAI} ${LIBERI} ${LIBESP}
	$(CC) -o $(NAME) $(OBJC) $(LIBRARY) ${LIBTWOBODYOVERLAP} ${LIBKINETIC} ${LIBNAI} ${LIBERI} ${LIBESP}

PREPARE:
	@(if [ ! -d "${BIN}" ]; then \
		${MKDIR}  ${BIN} ;\
	fi)

# libtwobodyoverlap.a
${LIBTWOBODYOVERLAP}:
	${MAKE} -f Makefile_twobodyoverlap

# libkinetic.a
${LIBKINETIC}:
	${MAKE} -f Makefile_kinetic

# libtwobodyoverlap.a
${LIBNAI}:
	${MAKE} -f Makefile_nai

# liberi.a
${LIBERI}:
	${MAKE} -f Makefile_eri

# libesp.a
${LIBESP}:
	${MAKE} -f Makefile_esp

###########################################################
#    here below are normal binary files
###########################################################

# testing codes  
${BIN}/main.o: main.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/coreguess_test.o: coreguess_test.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/atomdenmtrx_test.o: atomdenmtrx_test.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# scf module
${BIN}/scf.o: ../../scf/src/scf.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/scfparam.o: ../../scf/src/scfparam.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/scfdiis.o: ../../scf/src/scfdiis.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/fock.o: ../../scf/src/fock.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gslscfminimizer.o: ../../scf/src/gslscfminimizer.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/scfconv.o: ../../scf/src/scfconv.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/scfintscontroller.o: ../../scf/src/scfintscontroller.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/scfeadiis.o: ../../scf/src/scfeadiis.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# result module
${BIN}/denmtrx.o: ../src/denmtrx.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/mo.o: ../src/mo.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/atomdenmtrx.o: ../src/atomdenmtrx.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# gints module
${BIN}/digest.o: ../../gints/src/digest.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/exrho.o: ../../gints/src/exrho.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/oneemtrx.o: ../../gints/src/oneemtrx.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gints2d.o: ../../gints/src/gints2d.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gints4d.o: ../../gints/src/gints4d.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/sigshellpairinfor.o: ../../gints/src/sigshellpairinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/shellpair.o: ../../gints/src/shellpair.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/oneintformula.o: ../../gints/src/oneintformula.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/kdigestinfor.o: ../../gints/src/kdigestinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/kdigest.o: ../../gints/src/kdigest.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/jdigest.o: ../../gints/src/jdigest.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/digestkints.o: ../../gints/src/digestkints.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/gintsinfor.o: ../../gints/src/gintsinfor.cpp
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

# xcints module
${BIN}/xcints.o: ../../xcints/src/xcints.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/xcintsinfor.o: ../../xcints/src/xcintsinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/atomgrids.o: ../../xcints/src/atomgrids.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/batchgrid.o: ../../xcints/src/batchgrid.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/dftbasisderiv1.o: ../../xcints/src/dftbasisderiv1.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/dftbasisderiv2.o: ../../xcints/src/dftbasisderiv2.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/dftbasisderiv3.o: ../../xcints/src/dftbasisderiv3.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/dftbasisderiv4.o: ../../xcints/src/dftbasisderiv4.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/molegrids.o: ../../xcints/src/molegrids.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/batchbasis.o: ../../xcints/src/batchbasis.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/sigatombasis.o: ../../xcints/src/sigatombasis.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/xcvarinfor.o: ../../xcints/src/xcvarinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/xcvar.o: ../../xcints/src/xcvar.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/dftmatrix.o: ../../xcints/src/dftmatrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/batchvar.o: ../../xcints/src/batchvar.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/batchfunc.o: ../../xcints/src/batchfunc.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/batchxcmtrx.o: ../../xcints/src/batchxcmtrx.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/lebedev.o: ../../xcints/src/lebedev.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

#xcfunc module
${BIN}/initfuncderives.o: ../../xcfunc/src/initfuncderives.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/becke88.o: ../../xcfunc/src/becke88.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/lyp.o: ../../xcfunc/src/lyp.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/pbex.o: ../../xcfunc/src/pbex.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/pw86x.o: ../../xcfunc/src/pw86x.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/pw91x.o: ../../xcfunc/src/pw91x.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/s.o: ../../xcfunc/src/s.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/vwn1rpa.o: ../../xcfunc/src/vwn1rpa.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/pbec.o: ../../xcfunc/src/pbec.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/pw86c.o: ../../xcfunc/src/pw86c.f 
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/pw91c.o: ../../xcfunc/src/pw91c.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/pw92c.o: ../../xcfunc/src/pw92c.f 
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/vdwbr89.o: ../../xcfunc/src/vdwbr89.f  
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/vwn5.o: ../../xcfunc/src/vwn5.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/becke05ex.o: ../../xcfunc/src/becke05ex.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/br94coor.o: ../../xcfunc/src/br94coor.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

${BIN}/br89exhole.o: ../../xcfunc/src/br89exhole.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

.PHONY: clean
clean:
	@rm -rf ${BIN}
	@rm -rf test

