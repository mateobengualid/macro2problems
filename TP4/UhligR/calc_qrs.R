# VERSION 4.0, November 2002, COPYRIGHT H. UHLIG.
# CALC_QRS.M calculates the matrices Q, R and S as well
# as some other matrices given the matrix P. I.e.W with the property [x(t)',y(t)',z(t)']=W [x(t)',z(t)'].
# (The program uses the name PP for the matrix P, etc.)

# Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
# However, you are not allowed to sell this software or otherwise impinge
# on its free distribution.

suppressPackageStartupMessages(library(matlab))
if (l_equ == 0) {
	RR = diag(0,m_states)
	VV = matrix(cbind(kron(t(NN),FF)+kron(diag(1,k_exog),(FF%*%PP+GG)),
                             kron(t(NN),JJ)+kron(diag(1,k_exog),KK) ))
} else {
	RR = - CC_plus%*%(AA%*%PP+BB)
	VV = rbind( cbind(kron(diag(1,k_exog),AA),   kron(diag(1,k_exog),CC)),
           cbind(kron(t(NN),FF)+kron(diag(1,k_exog),(FF%*%PP+JJ%*%RR+GG)),
                             kron(t(NN),JJ)+kron(diag(1,k_exog),KK) ))
} # end if l_equ==0

if ((rk(VV) < k_exog*(m_states+n_endog)) & (IGNORE_VV_SING!=1)) {
   message = rbind('SOLVE.M: Sorry! V is not invertible.  Cannot solve for QQ and SS. You  ',
              '         can try setting IGNORE_VV_SING = 1 and wish for the best...   ')
   if (DISPLAY_IMMEDIATELY) {print(message)}
   Warnings = rbind(Warnings,message)
} else {
   if ((rk(VV) < k_exog*(m_states+n_endog))) {
     message = rbind('SOLVE.M: Warning! V is not invertible.  However, you have set          ',
                '         IGNORE_VV_SING = 1, and thus, since you have told me to       ',
                '         ignore this, I will proceed.  Keep your fingers crossed...    ')
     if (DISPLAY_IMMEDIATELY) {print(message)}
     Warnings = rbind(Warnings,message)
   }
   LLNN_plus_MM = LL%*%NN + MM
   QQSS_vec = solve(- VV) %*% matrix(c(DD,
                       LLNN_plus_MM))
   if (max(abs(QQSS_vec)) == Inf) {
      message = rbind('SOLVE.M: You probably are in trouble!  QQ or SS contain undefined      ',
                 '         entries! Most likely, the matrix VV is not invertible.        ')
      if (DISPLAY_IMMEDIATELY) {print(message)}
      Warnings = rbind(Warnings,message)
   }           
   QQ = vec(QQSS_vec[1:(m_states*k_exog)])	
   QQ = matrix(QQ,nrow=m_states,ncol=k_exog)
   SS = vec(QQSS_vec[(m_states*k_exog+1):((m_states+n_endog)*k_exog)])
   SS = matrix(SS,nrow=n_endog,ncol=k_exog)
   WW = rbind(cbind( diag(1,m_states)         , matrix(0,nrow=m_states,ncol=k_exog)),
          cbind(RR%*%ginv(PP)           , (SS-RR%*%ginv(PP)%*%QQ)), 
          cbind(matrix(0,nrow=k_exog,ncol=m_states), diag(1,k_exog))            )
}