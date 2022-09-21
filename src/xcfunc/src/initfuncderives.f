C>
C>
C> comment for the code initilizing the functional derivatives:
C> the code here is is very general, it does not care about
C> the variable sequence and variable types etc. If you have new
C> variables adding in, and you have modify the varlist.inc;
C> then here is fine
C> the only place which you need to take care is get_num_var_d1
C> if you have variables more than gamma(3), you need to revise
C> the code just below
C> \author  fenglai liu and jing kong
C>
C>

      FUNCTION get_num_var_d1(VAR)
C>  get the variable number for the given type
C>  of variables
      IMPLICIT NONE
#include "varlist.inc" 
      INTEGER VAR, N
      INTEGER get_num_var_d1

      N = -1
      IF (VAR == ID_GAMMA) THEN
         N = 3
      ELSE
         N = 2
      END IF
      get_num_var_d1 = N
      END

      FUNCTION get_num_var_d2(VAR1,VAR2)
C>  get the variable pairs number for the given two type
C>  of variables
      IMPLICIT NONE
      INTEGER VAR1, VAR2, NVAR1, NVAR2, N
      INTEGER get_num_var_d2
      INTEGER get_num_var_d1
      EXTERNAL get_num_var_d1

      NVAR1 = get_num_var_d1(VAR1)
      NVAR2 = get_num_var_d1(VAR2)
      N = -1
      IF (VAR1 == VAR2) THEN
         N = NVAR1*(NVAR1+1)/2
      ELSE
         N = NVAR1*NVAR2
      END IF
      get_num_var_d2 = N
      END

      FUNCTION get_num_var_d3(VAR1,VAR2,VAR3)
C>  get the variable pairs number for the given three type
C>  of variables
C>  here we have additional note:
C>  the var1, var2 and var3 has the following sequence:
C>  var1<=var2<=var3
C>  the < relation is defined as:
C>  rho < gamma < tau < lap < exrho
C>  therefore, if var1 == var3, then we must have
C>  var1 == var2 == var3
C>  that's the reason why we did not test var1 == var3 and 
C>  var2 != var3
C>
      IMPLICIT NONE
      INTEGER VAR1, VAR2, VAR3, NVAR1, NVAR2, NVAR3, N
      INTEGER get_num_var_d3
      INTEGER get_num_var_d1
      EXTERNAL get_num_var_d1

      NVAR1 = get_num_var_d1(VAR1)
      NVAR2 = get_num_var_d1(VAR2)
      NVAR3 = get_num_var_d1(VAR3)
      N = -1
      IF (VAR1 == VAR2 .and. VAR1 == VAR3) THEN
         N = NVAR1*(NVAR1+1)*(NVAR1+2)/6
      ELSE IF (VAR1 == VAR2 .and. VAR2 /= VAR3) THEN
         N = NVAR1*(NVAR1+1)*NVAR3/2
      ELSE IF (VAR1 /= VAR2 .and. VAR2 == VAR3) THEN
         N = NVAR1*NVAR2*(NVAR2+1)/2
      ELSE
         N = NVAR1*NVAR2*NVAR3
      END IF
      get_num_var_d3 = N
      END

      SUBROUTINE set_dvar_array(NVAR,VAL,OFFSET,DVAR)
C>  set the dvar array according to the value
C>  \param NVAR:   length of array that going to be set
C>  \param VAL:    the value used for setting  
C>  \param OFFSET: the offset for setting dvar array
C>  \param DVAR:   the functional derivatives information array
C>
      IMPLICIT NONE
      INTEGER NVAR, VAL, OFFSET, I, VAR_INDEX
      INTEGER DVAR(*)

      ! if value is < 0, then this part of contents
      ! is not touched for current functional
      ! else, the value would be the real variable 
      ! index
      IF (VAL < 0) THEN
         DO I = 1, NVAR
            DVAR(OFFSET+I) = VAL
         END DO
      ELSE
         VAR_INDEX = VAL
         DO I = 1, NVAR
            DVAR(OFFSET+I) = VAR_INDEX
            VAR_INDEX = VAR_INDEX + 1
         END DO
      END IF
      END

      SUBROUTINE init_func_deriv_3(VAR_INFOR,D3VARS)
C> we use this function to initilize the D3VARS
C> \param VAR_INFOR: input variable type information
C> \return D3VARS: third order functional derivatives 
C>                 position infor array
      IMPLICIT NONE
#include "varlist.inc" 
      INTEGER VAR_INFOR(*) 
      INTEGER D3VARS(*)
      INTEGER VT1,VT2,VT3,NVAR
      INTEGER OFFSET       ! offset of index for the D2VARS
      INTEGER VAR_INDEX    ! value for the D2VARS
      INTEGER get_num_var_d3
      EXTERNAL get_num_var_d3
      
      ! init the d3vars
      OFFSET    = 0
      VAR_INDEX = 1
      DO VT3 = 1, MAX_VAR_TYPE
         DO VT2 = 1, VT3
            DO VT1 = 1, VT2
               NVAR = get_num_var_d3(VT1,VT2,VT3)
               IF (VAR_INFOR(VT1) < 0 .or. VAR_INFOR(VT2) < 0 .or. 
     $ VAR_INFOR(VT3) < 0) THEN
                  CALL set_dvar_array(NVAR,-1,OFFSET,D3VARS)
               ELSE
                  CALL set_dvar_array(NVAR,VAR_INDEX,OFFSET,D3VARS)
                  VAR_INDEX = VAR_INDEX + NVAR
               END IF
               OFFSET = OFFSET + NVAR
            END DO
         END DO
      END DO
      END

      SUBROUTINE init_func_deriv_2(VAR_INFOR,D2VARS)
C> we use this function to initilize the D2VARS
C> \param VAR_INFOR: input variable type information
C> \return D2VARS: second order functional derivatives 
C>                 position infor array
      IMPLICIT NONE 
#include "varlist.inc" 
      INTEGER VAR_INFOR(*) 
      INTEGER D2VARS(*)
      INTEGER VT1,VT2,NVAR
      INTEGER OFFSET       ! offset of index for the D2VARS
      INTEGER VAR_INDEX    ! value for the D2VARS
      INTEGER get_num_var_d2
      EXTERNAL get_num_var_d2
      
      ! init the d2vars
      OFFSET    = 0
      VAR_INDEX = 1
      DO VT2 = 1, MAX_VAR_TYPE
         DO VT1 = 1, VT2
            NVAR = get_num_var_d2(VT1,VT2)
            IF (VAR_INFOR(VT1) < 0 .or. VAR_INFOR(VT2) < 0) THEN
               CALL set_dvar_array(NVAR,-1,OFFSET,D2VARS)
            ELSE
               CALL set_dvar_array(NVAR,VAR_INDEX,OFFSET,D2VARS)
               VAR_INDEX = VAR_INDEX + NVAR
            END IF
            OFFSET = OFFSET + NVAR
         END DO
      END DO
      END

      SUBROUTINE init_func_deriv_1(VAR_INFOR, D1VARS)
C> we use this function to initilize the D1VARS
C> \param VAR_INFOR: input variable type information
C> \return D1VARS: first order functional derivatives 
C>                 position infor array
      IMPLICIT NONE
#include "varlist.inc" 
      INTEGER VAR_INFOR(*)
      INTEGER D1VARS(*)
      INTEGER VAR_INDEX,OFFSET,V,NVAR
      INTEGER get_num_var_d1
      EXTERNAL get_num_var_d1

      ! init the d1vars
      OFFSET    = 0
      VAR_INDEX = 1
      DO V = 1, MAX_VAR_TYPE
         NVAR = get_num_var_d1(V)
         IF (VAR_INFOR(V) < 0) THEN
            CALL set_dvar_array(NVAR,-1,OFFSET,D1VARS)
         ELSE
            CALL set_dvar_array(NVAR,VAR_INDEX,OFFSET,D1VARS)
            VAR_INDEX = VAR_INDEX + NVAR
         END IF
         OFFSET = OFFSET + NVAR
      END DO
      END

