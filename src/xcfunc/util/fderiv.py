"""
This simple python script is used to generate the higher functional
derivatives common block for Fortran code

now we rewrite it to drop the use of common block
we just make everything in Fortran as an array
and we do not need to initilization code anymore
the initilization code is directly written by fortran
see function of initfuncderives.f in the xcfunc folder
"""
import sys
import os

__author__  = "Fenglai Liu"
__date__    = "June, 2012"

# setting up basic var group
# if you have more variables want to add here, just edit this line
# order is important!!!
varGroup = ["RHO","GAMMA","TAU","LAP","EXRHO"]

# include name, up to 3rd order
fortran_include_files = ["fderiv1.inc", "fderiv2.inc", "fderiv3.inc"]

# lines indentation record in the program
fortran_indent = 0

# the maximum order for the functional derivatives
MAX_ORDER = 3

# current variable number, another place you need to modify
MAX_VAR_NUM = 11

###################################################################
#       basic function related to the varaible information        #
###################################################################
def getVars(varType):
    """
    get the group of vars by given this group name
    """
    if varType == "RHO":
        return ("RA","RB")
    elif varType == "GAMMA":
        return ("GAA","GAB","GBB")
    elif varType == "TAU":
        return ("TA","TB")
    elif varType == "LAP":
        return ("LA","LB")
    elif varType == "EXRHO":
        return ("EXA","EXB")

def getVarType(var):
    """
    get the var type by given the var name
    """
    if var == "RA" or var == "RB":
        return "RHO"
    elif var == "GAA" or var == "GAB" or var == "GBB":
        return "GAMMA"
    elif var == "TA" or var == "TB":
        return "TAU"
    elif var == "LA" or var == "LB":
        return "LAP"
    elif var == "EXA" or var == "EXB":
        return "EXRHO"

###################################################################
#                         printing functions                      #
###################################################################
def printCode(f,line):
    """
    indenting the code and print the given line
    """
    global fortran_indent

    # for fortran files, there are more things to watch
    if line[0] == "C" and line[1] == " ":  # this is comment line
        f.write(line)
        f.write("\n")
        return
    if line[0] == "$":  # lines should be connected to the above line
        f.write("     ")
        f.write(line)
        f.write("\n")
        return

    # now it's normal printing
    indentLength = 6
    indent = fortran_indent + indentLength
    for i in range(indent):
        f.write(" ")
    f.write(line)
    f.write("\n")


def printEmptyLine(f,n=1):
    """
    print out empty lines
    """
    if n > 0:
        for i in range(n):
            f.write("\n")
    else:
        print "n is not correct in printEmptyLine function\n"
        sys.exit()


def increaseIndentation():
    """
    increase the indent according to the file type
    """
    global fortran_indent
    fortran_indent = fortran_indent + 3


def decreaseIndentation():
    """
    decrease the indent for 3
    """
    global fortran_indent
    fortran_indent = fortran_indent - 3


def initializeIndentation():
    """
    initilize the indent for cpp and fortran files
    """
    global fortran_indent
    fortran_indent = 0

###################################################################
#      code generation, for 1st to 3rd functional derivatives     #
###################################################################
def varGeneration(order):
    """
    generate the variables list
    In this list we consider all possible vars arrangement
    This function should generate the same order variables as in emul program
    see the constructor of xcvar class
    """
    varNames = [ ]
    varPrefix = "ID"
    if order == 1:
        for var in varGroup:
            varlist = getVars(var)
            for v in varlist:
                name = varPrefix + "_" + v
                varNames.append(name)
    elif order == 2:
        for ivar in varGroup:
            for jvar in varGroup:
                if varGroup.index(jvar) > varGroup.index(ivar):
                    continue
                ivarlist = getVars(ivar)
                jvarlist = getVars(jvar)
                for iv in ivarlist:
                    for jv in jvarlist:
                        if ivar == jvar and jvarlist.index(jv) > ivarlist.index(iv):
                            continue
                        name = varPrefix + "_" + jv + "_" + iv
                        varNames.append(name)
    elif order == 3:
        for ivar in varGroup:
            for jvar in varGroup:
                for kvar in varGroup:
                    if varGroup.index(jvar) > varGroup.index(ivar):
                        continue
                    if varGroup.index(kvar) > varGroup.index(jvar):
                        continue
                    ivarlist = getVars(ivar)
                    jvarlist = getVars(jvar)
                    kvarlist = getVars(kvar)
                    for iv in ivarlist:
                        for jv in jvarlist:
                            for kv in kvarlist:
                                if ivar == jvar and jvarlist.index(jv) > ivarlist.index(iv):
                                    continue
                                if kvar == jvar and kvarlist.index(kv) > jvarlist.index(jv):
                                    continue
                                name = varPrefix + "_"  + kv + "_" + jv + "_" + iv
                                varNames.append(name)
    else:
        print "Order is <=3 right now in the varGeneration"
        sys.exit()

    return varNames


def printHeadFiles():
    """
    generate the head files for fortran codes and the cpp codes
    """
    # initilize the indent
    initializeIndentation()

    # now doing fortan codes
    for order in range(MAX_ORDER):
        order = order + 1

        # open the corresponding include file
        file = fortran_include_files[order-1]
        f = open(file, "w")
        line = "C Here it's the functional derivatives information for Fortran codes"
        printCode(f,line)
        line = "C The information here is about the position of variables."
        printCode(f,line)
        line = "C How can we use it? It's quite simple."
        printCode(f,line)
        line = "C For example, as for ID_RA_GAA(second functional derivative"
        printCode(f,line)
        line = "C with respect to RA and GAA), we can get the real functional"
        printCode(f,line)
        line = "C derivatives postion of POS_RA_GAA as:"
        printCode(f,line)
        line = "C POS_RA_GAA = D2VARS[ID_RA_GAA]"
        printCode(f,line)
        line = "C D2VARS is the array holds the functional deriavtives infor"
        printCode(f,line)
        line = "C see the init_func_deriv_2 function for more information"
        printCode(f,line)
        printEmptyLine(f,2)

        # generate the variable names
        varNames = varGeneration(order)
        for v in varNames:
            line = "INTEGER " + v
            printCode(f,line)
        for v in varNames:
            pos = varNames.index(v) + 1
            line = "PARAMETER(" + v + " = " + str(pos) + ")"
            printCode(f,line)

        # print out empty lines to proceed to next order
        printEmptyLine(f,2)

        # finally about the array
        # we just print the dimension of the array
        num = MAX_VAR_NUM
        arrayLen  = "N_FUNC_DERIV_1"
        if order == 2:
            num = MAX_VAR_NUM*(MAX_VAR_NUM+1)/2
            arrayLen  = "N_FUNC_DERIV_2"
        elif order == 3:
            num = MAX_VAR_NUM*(MAX_VAR_NUM+1)*(MAX_VAR_NUM+2)/6
            arrayLen  = "N_FUNC_DERIV_3"
        elif order > MAX_ORDER:
            print "Improper order given in the printHeadFiles, only 1, 2 or 3"
            sys.exit()

        # now print out array
        line = "C Here finally it's the array length"
        printCode(f,line)
        line = "INTEGER " + arrayLen
        printCode(f,line)
        line = "PARAMETER(" + arrayLen + " = " + str(num) + ")"
        printCode(f,line)

        # print out empty lines to end the file
        printEmptyLine(f,2)

        # close the fortran file
        f.close()


########################################################
#  real working part
########################################################
printHeadFiles()
