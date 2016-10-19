# VERSION 4.1, MAY 2003, COPYRIGHT H. UHLIG.
#   May 2003: major correction: autocorrelation calculation interchanged lags and leads before,
#      due to overlooking, that for complex vectors v, say, v' is the complex-conjugate transpose,
#      not just the transpose.  This has now been corrected.
#      Thanks goes to Mathias Trabandt for pointing out the error.
#      
# MOMENTS.M calculates variances, covariances, and autocorrelation 
# tables, using frequency domain techniques with and without HP-Filtering.
# It is assumed that the endogeneous state variables x(t), the other endogenous
# variables y(t) and the exogeneous state variables z(t) obey the law of motion
# x(t) = PP x(t-1) + QQ z(t),
# y(t) = RR x(t-1) + SS z(t),
# z(t) = NN z(t-1) + epsilon(t), VAR(epsilon(t)) = Sigma,
# It is assumed that the matrix WW with the property
# [x(t)',y(t)',z(t)']' = WW [x(t)', z(t)']' has been calculated with e.g. SOLVE.M.
# The program concentrates on 
# v' = [ x(t)' y(t)' z(t)' ]'(HP_SELECT),
# where HP_SELECT can be chosen to restrict attention to desired variables
# or to reorder them so that e.g. GNP is at the top.
#
# The following options need be chosen beforehand:
# PERIOD: number of periods per year (e.g. = 12 for monthly, = 4 for quarterly)
# HP_SELECT: A vector selecting the variables for which the computations should be
#     performed.  HP_SELECT = 1: (m+n+k) selects all variables.
# GNP_INDEX:  The autocorrelation table is calculated with respect to v(GNP_INDEX).
#     If GNP_INDEX is the index of GNP in vector, the autocorrelation table will provide
#     information about the correlations of the variables with leads and lags of GNP.
#     Note that GNP_INDEX must be the index of GNP in the vector v, which just
#     contains a subset of all variables, if HP_SELECT is not equal to 1: (m+n+k).
# 
# The program provides the following results, both for the unfiltered, ``raw'' series 
# (extension: _raw) as well as for the HP-filtered series (extension: _fil).  For example,
# covmat_raw is the contemporaneous covariance matrix for the unfiltered series, whereas
# covmat_fil is the contemporaneous covariance matrix for the filtered series.  The computed
# variables (without extension) are:
#
#  autcov:  matrix-autocovariance function, autcov(j+1,:) = E[v(t) * v(t-j)'], columnwise vectorized.
#  covmat: contemporaneous covariance matrix
#  cor: Vector of instantaneous correlations of v(t) with v(t,GNP_INDEX);
#  autcor: autcor_fil(N_LEADS_LAGS + j,:) is the correlation of v(t+j) with v(t,GNP_INDEX),
#        where v is HP-filtered (autcor_raw similar for unfiltered v).  Leaders have negative j.
#        The last row indicates the lag j. This is the autocorrelation table often used in the RBC literature.
#  au1mat: = E[v(t) * v(t-1)'], i.e. the matrix of first-order autocovariances.  Useful for calculating regressions.
#  varvec: vector of the variances, i.e. the diagonal of covmat
#
# Other variables that are computed, thus possibly overwriting variables in use with the same
# name are:
# m_states, n_endog, k_exog, n_select, WW_sel, freqs, svv_raw, svv_fil,  hp, grdpt, im, z, zi, ss, ssd 
# gnp_select, diag_select, cov0_raw, cov0_fil, autcor1_raw, autcor1_fil,
# autcor2_raw, autcor2_fil, leadslags


# Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
# However, you are not allowed to sell this software or otherwise impinge
# on its free distribution.



cat('MOMENTS: Please be patient for a few seconds...\n')
cat('MOMENTS: Performing frequency-domain calculations...\n')
suppressPackageStartupMessages(library(matlab))
m_states=dim(QQ)[1]
k_exog=dim(QQ)[2]
n_endog=dim(SS)[1]
n_select = max(dim(matrix(HP_SELECT)))

#It is assumed that the following matrix WW has been calculated:
#WW = [ eye(m_states)         , zeros(m_states,k_exog)
#       RR*pinv(PP)           , (SS-RR*pinv(PP)*QQ) 
#       zeros(k_exog,m_states), eye(k_exog)            ];
WW_sel = WW[HP_SELECT,]


#CALCULATIONS:
freqs = seq(0,(2*pi*(1 - .5/N_GRIDPOINTS)),((2*pi)/N_GRIDPOINTS))
# frequencies run from 0 , 2*pi/N_GRIDPOINTS , ... , 2*pi - 2*pi/N_GRIDPOINTS
hp = 4*HP_LAMBDA*(1 - cos(freqs))^2 / (1 + 4*HP_LAMBDA*(1 - cos(freqs))^2)
# Transfer function for the HP-filter
svv_raw <- NULL
svv_fil <- NULL
im <- 0+1i
PP_may_be_singular = ( min(abs(eigen(PP)$values)) < MOM_TOL )

for (grdpt in 1 : N_GRIDPOINTS) {
   z   = exp( im*freqs[grdpt])
   zi  = exp(-im*freqs[grdpt])
   ss  = 1.0/(2*pi) * rbind( solve(diag(1,m_states)-PP*zi)%*%QQ , diag(k_exog) ) %*% 
                      solve(diag(1,k_exog)-NN*zi) %*% Sigma %*%
                      solve(diag(1,k_exog)-t(NN)*z) %*%
                      cbind( t(QQ)%*%solve(diag(1,m_states)-t(PP)*z) , diag(1,k_exog) )
   if (PP_may_be_singular) {
     ssd = rbind(diag(1,m_states,m_states+k_exog) ,cbind(RR*zi,SS ),cbind(diag(0,k_exog,m_states),diag(1,k_exog)))%*%
        ss%*%cbind(t(diag(1,m_states,m_states+k_exog)),rbind(t(RR)*z,t(SS)),rbind(diag(0,m_states,k_exog),diag(1,k_exog)))
     ssd = ssd[HP_SELECT,HP_SELECT]
   } else {
      ssd = WW_sel %*% ss %*% t(WW_sel)
   } # end if
   # OLD VERSION, BEFORE MAY 2003:
   # svv_raw = [ svv_raw ; (ssd(:))' ];   # Columnwise vectorization of ssd, add to svv_raw
   # ssd = hp(grdpt)^2* ssd;    # Filtered version
   # svv_fil = [ svv_fil ; (ssd(:))' ];   # Columnwise vectorization of ssd, add to svv_fil
   #  ... with the intention of storing the column ssd(:) as a row-vector.  But I overlooked
   #      that (ssd(:))' will not just take the transpose but also calculate the complex-conjugate.
   #      Since we only need the transpose, we take the complex-conjugate once more.
   # Thus, NEW VERSION, AFTER MAY 2003:
   svv_raw = rbind( svv_raw , Conj(t(matrix(ssd))))   # Columnwise vectorization of ssd, add to svv_raw
   ssd = hp[grdpt]^2* ssd;    # Filtered version
   svv_fil = rbind( svv_fil , Conj(t(matrix(ssd))))   # Columnwise vectorization of ssd, add to svv_fil

} # end for

# disp('MOMENTS: Done stacking spectral density matrix...');

# For the unfiltered, "raw" series:
# ======================
autcov_raw = Re(mvfft(svv_raw,inverse=T)/nrow(svv_raw)) * (2 * pi);  # matrix-autocovariance function, raw series. 
                                              # autcov(j+1,:) = E[v(t) * v(t-j)'], columnwise vectorized.
covmat_raw = autcov_raw[1,] # covariance matrix at lag zero
covmat_raw = matrix(covmat_raw,n_select,n_select); # makes it into quadratic matrix,
                                        # inverse to columnwise vectorization
gnp_select = n_select*(0:(n_select-1)) + GNP_INDEX;  # Selecting the row A(GNP_INDEX,:) from the 
                                                                                    # columnwise vectorization
diag_select = (n_select+1)*(0:(n_select-1)) + 1; # Selects the diagonal
cov00_raw = autcov_raw[1,((GNP_INDEX-1)*n_select+GNP_INDEX)]  # Variance of v(t,GNP_INDEX)
cov0_raw = sqrt(cov00_raw)*matrix(1,N_LEADS_LAGS,1)%*%sqrt(autcov_raw[1,diag_select])

# weird thing here: autcor2 and autcor1 swapped positions

autcor2_raw = autcov_raw[1:N_LEADS_LAGS,gnp_select] / cov0_raw
# autcor1_raw(j+1) = autocorr. function of v(t,GNP_INDEX) with v(t-j).  If this is positive, 
# v(t) is leading as compared to v(t,GNP_INDEX).
gnp_select = ((GNP_INDEX - 1)*n_select + 1) : (GNP_INDEX*n_select) # Selecting the column A(:,GNP_INDEX) from the columnwise vectorization
diag_select = (n_select+1)*(0:(n_select-1)) + 1
autcor1_raw = autcov_raw[1:N_LEADS_LAGS,gnp_select] / cov0_raw

# autcor2_raw(j+1) = autocorr. function of v(t) with v(t-j,GNP_INDEX).  If this is positive,
# v(t) is lagging as compared to v(t,GNP_INDEX).
cor_raw = matrix(autcor1_raw[1,])  # Vector of instantaneous correlations of v(t) with v(t,GNP_INDEX);
leadlags = (1 - N_LEADS_LAGS) :  (N_LEADS_LAGS-1)
autcor_raw = t(cbind(rbind(flipud(autcor1_raw) , autcor2_raw[2:N_LEADS_LAGS,]),matrix(leadlags)))
# autcor_raw(N_LEADS_LAGS + j,:) is the corr of v(t+j) with v(t,GNP_INDEX).  Leaders have negative j.
# This is the "autocorrelation table" often showing up in the RBC literature.
# The last row indicates the lag.
au1mat_raw = autcov_raw[2,] # covariance matrix at lag one
au1mat_raw = matrix(au1mat_raw,n_select,n_select,byrow=T) # makes it into quadratic matrix,
                                # inverse to columnwise vectorization. au1mat_raw = E[v(t) * v(t-1)'],
varvec_raw = matrix(diag(covmat_raw))  # just the variances

# For the HP-filtered series:
# ======================
autcov_fil = Re(mvfft(svv_fil,inv=T)/nrow(svv_fil)) * (2 * pi);  # matrix-autocovariance function, raw series. 
                                              # autcov(j+1,:) = E[v(t) * v(t-j)'], columnwise vectorized.
covmat_fil = autcov_fil[1,] # covariance matrix at lag zero
covmat_fil = matrix(covmat_fil,n_select,n_select) # makes it into quadratic matrix,
                                        # inverse to columnwise vectorization
gnp_select = n_select*(0:(n_select-1)) + GNP_INDEX  # Selecting the row A(GNP_INDEX,:) from the 
                                                                                    # columnwise vectorization
diag_select = (n_select+1)*(0:(n_select-1)) + 1 # Selects the diagonal
cov00_fil = autcov_fil[1,((GNP_INDEX-1)*n_select+GNP_INDEX)]  # Variance of v(t,GNP_INDEX)
cov0_fil = sqrt(cov00_fil)*matrix(1,N_LEADS_LAGS,1)%*%sqrt(autcov_fil[1,diag_select])
#same weird thing happens here with autcor2 and autcor1 switching viz matlab
autcor2_fil = autcov_fil[1:N_LEADS_LAGS,gnp_select] / cov0_fil
# autcor1_fil(j+1) = autocorr. function of v(t,GNP_INDEX) with v(t-j).  If this is positive, 
# v(t-j) is leading as compared to v(t,GNP_INDEX).
gnp_select = ((GNP_INDEX - 1)*n_select + 1) : (GNP_INDEX*n_select) # Selecting the column A(:,GNP_INDEX) from the columnwise vectorization
diag_select = (n_select+1)*(0:(n_select-1)) + 1
autcor1_fil = autcov_fil[1:N_LEADS_LAGS,gnp_select] / cov0_fil
# autcor2_fil(j+1) = autocorr. function of v(t) with v(t-j,GNP_INDEX).  If this is positive,
# v(t+j) is lagging as compared to v(t,GNP_INDEX).
cor_fil = matrix(autcor1_fil[1,])  # Vector of instantaneous correlations of v(t) with v(t,GNP_INDEX);
leadlags = (1 - N_LEADS_LAGS) :  (N_LEADS_LAGS-1)
autcor_fil = t(cbind( rbind( flipud(autcor1_fil) , autcor2_fil[2:N_LEADS_LAGS,] ) , matrix(leadlags) ))
# autcor_fil(N_LEADS_LAGS + j,:) is the corr of v(t+j) with v(t,GNP_INDEX).  Leaders have negative j.
# This is the "autocorrelation table" often showing up in the RBC literature.
# The last row indicates the lag.
au1mat_fil = autcov_fil[2,] # covariance matrix at lag one
au1mat_fil = matrix(au1mat_fil,n_select,n_select,byrow=T) # makes it into quadratic matrix,
                                # inverse to columnwise vectorization. au1mat_fil = E[v(t) * v(t-1)'],
varvec_fil = matrix(diag(covmat_fil))  # just the variances

cat('MOMENTS: Done.\n')
