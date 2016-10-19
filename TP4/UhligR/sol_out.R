# VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
# SOL_OUT.M prints the coefficients of the decision rules,
# delivered by SOLVE.M.
# It is assumed, that VARNAMES, a matrix with m+n+k rows has
# been set, containing the names of all the variables.
# This program overwrites m_states, k_exog and n_endog.


# Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
# However, you are not allowed to sell this software or otherwise impinge
# on its free distribution.


m_states=dim(QQ)[1]
k_exog=dim(QQ)[2]
n_endog = dim(SS)[1]
k_exog = dim(SS)[2]

cat('Exogenous states z(t): \n')
print(matrix(VARNAMES[(m_states+n_endog+1):(m_states+n_endog+k_exog)]),quote=F)
cat('\n ')
cat('Endogenous states x(t): \n')
print(matrix(VARNAMES[1:m_states]),quote=F)
cat('\n ')
if (DISPLAY_ROOTS) {
   cat('All the roots are: \n')
   cat('    root         abs(root)   \n')
   print(cbind(diag(Xi_eigval[Xi_sortindex,Xi_sortindex]),
   abs(diag(Xi_eigval[Xi_sortindex,Xi_sortindex])))) 
   cat('\n The chosen roots are: \n')
   cat('    root         abs(root)    \n')
   print(cbind(diag(Lambda_mat),abs(diag(Lambda_mat))))
   cat('\n ')
}
cat('PP: Recursive equilibrium law of motion for x(t) on x(t-1): \n')
print(PP)
cat('\n QQ: Recursive equilibrium law of motion for x(t) on z(t): \n')
print(QQ)
cat('\n ')
cat('\n Other endogenous variables y(t): \n')
print(matrix(VARNAMES[(m_states+1):(m_states+n_endog)]),quote=F)
cat('\n ')
cat('\n RR: Recursive equilibrium law of motion for y(t) on x(t-1): \n')
print(RR)
cat('\n SS: Recursive equilibrium law of motion for y(t) on z(t): \n')
print(SS)
cat('\n ')