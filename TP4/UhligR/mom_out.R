# VERSION 4.1, May 2003, COPYRIGHT H. UHLIG.
# Changes: Some improvements concerning graphics.
# Plus no "HP-filtered" in title of tables, if no HP Filtering was done.
# 
# MOM_OUT produces output from the calculations done with MOMENTS.M,
# which is assumed to have been run just before.
# This program should be modified to suit tastes and needs.  Some options
# are given in the first few lines of this program.
# It is assumed that
# VARNAMES, a matrix with (m+n+k) rows, containing the variable names, has been set.
# The program overwrites freq1, diag_select, hndl, var_index


# Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
# However, you are not allowed to sell this software or otherwise impinge
# on its free distribution.

op <- par(no.readonly=T)
if (DO_HP_GRAPH) {
    
#-----------------------------------MOM JOINT------------------------------------------------------------------------    
   if (MOM_JOINT) {
       freq1 = floor(N_GRIDPOINTS/24)
       diag_select = (n_select+1)*(0:(n_select-1)) + 1 # Selects the diagonal
      if (MOM_PLOT_RAW) {
          title='Spectral densities, unfiltered'
      } else {
          title='Spectral densities, Hodrick-Prescott filtered'
      }
          if (MOM_PLOT_RAW) {
                  matplot(freqs[freq1:(N_GRIDPOINTS/2)], Re(svv_raw[freq1:(N_GRIDPOINTS/2),diag_select]),
			type='l',lty=1,lwd=2,main=title,xlab='Frequency')
                    for (var_index in 1:n_select) {
                    text(freqs[MOM_TXT_MARKER], Re(svv_raw[MOM_TXT_MARKER,diag_select[var_index]]),
                    VARNAMES[HP_SELECT[var_index],])
                    } # end for
              } else { #MOM_PLOT_RAW
                  matplot(freqs[1:(N_GRIDPOINTS/2)], Re(svv_fil[1:(N_GRIDPOINTS/2),diag_select]),
			type='l',lty=1,lwd=2,main=title,xlab='Frequency')
                     for (var_index in 1:n_select) {
                     text(freqs[MOM_TXT_MARKER], Re(svv_fil[MOM_TXT_MARKER,diag_select[var_index]]),
                     VARNAMES[HP_SELECT[var_index],])
                     }
           } # end if MOM_PLOT_RAW    
        grid()
      dummy<-readline('Inspect figure. Hit ENTER when ready...')
      
  } # end if MOM_JOINT
#----------------------------------------------END MOM JOINT------------------------------------------------------------------------


#-----------------------------------MOM SUBPLOT-------------------------------------------------------------------------------------   
   if (MOM_SUBPLOT) {
       
       freq1 = floor(N_GRIDPOINTS/24)
       diag_select = t(matrix((n_select+1)*(0:(n_select-1)) + 1)) # Selects the diagonal
       diag_length=length(diag_select)
       
        size_plot_x=ceil(sqrt(n_select))
        size_plot_y=round(sqrt(n_select))
        par(mfrow=c(size_plot_y,size_plot_x))
          
       for (mom_index in 1:diag_length) {
      if (MOM_PLOT_RAW) {
          title='Spectrum, unfiltered'
      } else {
          title='Spectrum,HP-filtered'
      }
          if (MOM_PLOT_RAW) {      
                  matplot(freqs[freq1:(N_GRIDPOINTS/2)], Re(svv_raw[freq1:(N_GRIDPOINTS/2),diag_select[1,mom_index]]),
			type='l',lty=1,xlab='Frequency',main=title,ylab="")
                  text((freqs[MOM_TXT_MARKER]*4), Re(svv_raw[MOM_TXT_MARKER,diag_select[1,mom_index]]),
                  VARNAMES[HP_SELECT[mom_index],])
            } else { #MOM_PLOT_RAW
                  matplot(freqs[1:(N_GRIDPOINTS/2)], Re(svv_fil[1:(N_GRIDPOINTS/2),diag_select[1,mom_index]]),
			type='l',lty=1,xlab='Frequency',main=title,ylab="")
                     text((freqs[MOM_TXT_MARKER]*4), Re(svv_fil[MOM_TXT_MARKER,diag_select[1,mom_index]]),
                     VARNAMES[HP_SELECT[mom_index],])
             } # end MOM_PLOT_RAW    
         #grid;
   } # end for mom_index
      
      dummy <-readline('Inspect figure. Hit ENTER when ready...')
	par(op)  
  } #if MOM_SUBPLOT
#----------------------------------------------END MOM SUBPLOT-------------------------------------------------------------------

   

#----------------------------------------------MOM SINGLE------------------------------------------------------------------------    
   if (MOM_SINGLE) {
       freq1 = floor(N_GRIDPOINTS/24)
       diag_select = t(matrix((n_select+1)*(0:(n_select-1)) + 1)) # Selects the diagonal
       diag_length=length(diag_select)
       for (mom_index in 1:diag_length) {
          if (MOM_PLOT_RAW) {
                  matplot(freqs[freq1:(N_GRIDPOINTS/2)], Re(svv_raw[freq1:(N_GRIDPOINTS/2),diag_select[1,mom_index]]),
			type='l', lty=1,lwd=2)
                    text(freqs[MOM_TXT_MARKER], Re(svv_raw[MOM_TXT_MARKER,diag_select[1,mom_index]]),
                    VARNAMES[HP_SELECT[mom_index],])
            } else { #MOM_PLOT_RAW      
                  matplot(freqs[1:(N_GRIDPOINTS/2)], Re(svv_fil[1:(N_GRIDPOINTS/2),diag_select[1,mom_index]]),
			type='l',lty=1,lwd=2)
                     text(freqs[MOM_TXT_MARKER], Re(svv_fil[MOM_TXT_MARKER,diag_select[1,mom_index]]),
                     VARNAMES[HP_SELECT[mom_index],])
             } # end MOM_PLOT_RAW    
           
        grid()
      if (MOM_PLOT_RAW) {
          title(main='Spectral densities, unfiltered',xlab='Frequency')
      } else { 
          title(main='Spectral densities, Hodrick-Prescott filtered',
		    xlab='Frequency')
      }
      dummy<-readline('Inspect figure. Hit ENTER when ready...')
  } # end for mom_index
  } # end if MOM_SINGLE
#----------------------------------------------END MOM SINGLE------------------------------------------------------------------------
} # end if DO_HP_GRAPH


if (DO_DISP1 | DO_DISP2 | DO_DISP3) {
   cat('\n')
   cat('MOM_OUT.M: Frequency-domain method based calculation of moments\n')
   cat('The variables are:\n')
   print(matrix(VARNAMES[HP_SELECT,]),quote=F)
   cat('\n')
}
if (DO_HP_FILTER) {
   if (DO_DISP1) {
      cat('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j\n')
      for (var_index in 1:n_select) {
         cat(sprintf('  %5.2f',autcor_fil[var_index,]))
	   cat('\n')
      } 
      cat(sprintf('  %5.0f',autcor_fil[n_select+1,]))
      cat('\n')
   }
   if (DO_DISP2) {
      cat('Variance-Covariance Matrix, HP-filtered series:\n')
      for (var_index in 1:n_select) {
         cat(sprintf(' %6.3f',covmat_fil[var_index,]))
	   cat('\n')
      }
      cat('\n')
   }
   if (DO_DISP3) {
      cat('Standard deviations, HP-filtered series:\n')
     print(sqrt(varvec_fil))
   }
} else {
   if (DO_DISP1) {
       cat('Autocorrelation Table, corr(v(t+j),GNP(t)).  Last row shows j\n')
      for (var_index in 1:n_select) {
         cat(sprintf('  %5.2f',autcor_raw[var_index,]))
	   cat('\n')
      }
      cat(sprintf('  %5.0f',autcor_raw[n_select+1,]))
      cat('\n')
   }
  if (DO_DISP2) {
      cat('Variance-Covariance Matrix:\n')
      for (var_index in 1:n_select) {
         cat(sprintf(' %6.3f',covmat_raw[var_index,]))
	   cat('\n')
      }
      cat('\n')
   }
   if (DO_DISP3) {
      cat('Standard deviations:\n')
     print(sqrt(varvec_raw))
   }
}  