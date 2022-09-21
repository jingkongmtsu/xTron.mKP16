C Here it's the functional derivatives information for Fortran codes
C The information here is about the position of variables.
C How can we use it? It's quite simple.
C For example, as for ID_RA_GAA(second functional derivative
C with respect to RA and GAA), we can get the real functional
C derivatives postion of POS_RA_GAA as:
C POS_RA_GAA = D2VARS[ID_RA_GAA]
C D2VARS is the array holds the functional deriavtives infor
C see the init_func_deriv_2 function for more information


      INTEGER ID_RA
      INTEGER ID_RB
      INTEGER ID_GAA
      INTEGER ID_GAB
      INTEGER ID_GBB
      INTEGER ID_TA
      INTEGER ID_TB
      INTEGER ID_LA
      INTEGER ID_LB
      INTEGER ID_EXA
      INTEGER ID_EXB
      PARAMETER(ID_RA = 1)
      PARAMETER(ID_RB = 2)
      PARAMETER(ID_GAA = 3)
      PARAMETER(ID_GAB = 4)
      PARAMETER(ID_GBB = 5)
      PARAMETER(ID_TA = 6)
      PARAMETER(ID_TB = 7)
      PARAMETER(ID_LA = 8)
      PARAMETER(ID_LB = 9)
      PARAMETER(ID_EXA = 10)
      PARAMETER(ID_EXB = 11)


C Here finally it's the array length
      INTEGER N_FUNC_DERIV_1
      PARAMETER(N_FUNC_DERIV_1 = 11)


