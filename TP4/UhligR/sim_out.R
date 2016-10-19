# VERSION 4.0, November 2002, COPYRIGHT H. UHLIG.
# SIM_OUT.M produces output for SIMUL.M.  It is controlled 
# by several options to be described
# below.  All these options are set in OPTIONS.M,
# see that file for further details.
#
# The series will be plotted
# if SIM_GRAPH is set to 1.  In that case, further
# modifications can be done in particular with:
#   DO_ENLARGE : = 1, if you want large font sizes for the text on your plots.
#                     Good for slides.
#   SIM_JOINT  : = 1, if you want all series on the same graph, else = 0.
#   PRINT_FIG  : = 1, if you want plots to be printed on your printer
#   SAVE_FIG   : = 1, if you want plots to be saved as encaps. postscript. 
#   SAVE_FIG_JPG : = 1, if you want plots to be saved as jpg-files.
#                     Set PRINT_FIG = 0 also. The filenames are sim_ser1.eps, ...
#                     if SIM_JOINT = 0, and sim_data.eps is SIM_JOINT = 1.
#   SIM_PLOT_RAW : = 1, if you want a plot of the raw, unfiltered series, 
#                       even though you have chosen to HP-filter, DO_HP_FILTER = 1.                       
#                       Note, that if you have chosen to save figures, then
#                       the previously saved figures will be overwritten.
#                       This option is useful, if you want to look at the plot
#                       of the raw simulations, after having already seen the filtered
#                       ones: simply type
#                       SIM_PLOT_RAW = 1;
#                       sim_out;
#
# For printing the numbers of the autocorrelation table, 
# the following options are helpful:
#   
#  SIM_DO_DISP1: Set to = 1 to see printout of the autocorrelation matrix. 
#  SIM_DO_DISP2: Set to = 1 to see printout of the variance-covariance matrix.
#  SIM_DO_DISP3: Set to = 1 to see printout of the vector of variances.


# Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
# However, you are not allowed to sell this software or otherwise impinge
# on its free distribution.

n_select = max(dim(matrix(SIM_SELECT)))
op <- par(no.readonly=T)
if (SIM_GRAPH) {
   time_axis = (0:(SIM_LENGTH-1))/PERIOD + SIM_DATE0
   
#---------------------------------------------JOINT PLOTS-----------------------------------------------------------------
   if (SIM_JOINT) {
      if (DO_HP_FILTER & !SIM_PLOT_RAW) {
          title = 'Simulated data (HP-filtered)'
      } else {
          title = 'Simulated data'
      } # end if
      if (SIM_PLOT_RAW) {
	   matplot(matrix(c(time_axis[1:min(SIM_LENGTH,SIM_MAX)])),
		     t(sim_raw[SIM_SELECT,1:min(SIM_LENGTH,SIM_MAX)]),
		     type='l',lwd=2,lty=1,main=title,xlab='Year',
		     ylab='Percent deviation from steady state')
    	   mx = apply(abs(t(sim_raw[SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX))])),2,max)
	   pos = apply(abs(t(sim_raw[SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX))])),2,
			   which.max)
         for (var_index in (1:n_select)) {
#            text(time_axis[SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
#              VARNAMES(SIM_SELECT(var_index),:));
            text(time_axis[pos[var_index]], sim_raw[var_index,pos[var_index]],
              VARNAMES[SIM_SELECT[var_index],])
         }
      } else {
         matplot(matrix(c(time_axis[1:min(SIM_LENGTH,SIM_MAX)])),
		     t(sim_xyz[SIM_SELECT,1:min(SIM_LENGTH,SIM_MAX)]),
		     type='l', lwd=2,main=title,xlab='Year',lty=1,
		     ylab='Percent deviation from steady state')
    	   mx = apply(abs(t(sim_xyz[SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX))])),2,max)
	   pos = apply(abs(t(sim_xyz[SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX))])),2,
			   which.max)
         for (var_index in 1:n_select) {
#            text(time_axis[SIM_TXT_MARKER], sim_xyz[var_index,SIM_TXT_MARKER],
#              VARNAMES[SIM_SELECT[var_index],])
            text(time_axis[pos[var_index]], sim_xyz[var_index,pos[var_index]],
              VARNAMES[SIM_SELECT[var_index],])
         } # end for
      } # end if
      grid()
	dummy <- readline("Inspect figure. Hit ENTER when ready...")      
      } # end if SIM_JOINT; 
#--------------------------------------------END JOINT PLOTS-------------------------------------------------------------

#--------------------------------------------SIM SUBPLOT-----------------------------------------------------------------
if (SIM_SUBPLOT) {
        size_plot_x=ceiling(sqrt(n_select))
        size_plot_y=round(sqrt(n_select))
	  par(mfrow=c(size_plot_x,size_plot_y))

      for (var_index in 1:n_select) {
         if (DO_HP_FILTER & !SIM_PLOT_RAW) {
            title=paste('Sim. data(HP-filt.):',VARNAMES[SIM_SELECT[var_index],])
         } else {
            title=paste('Sim. data:',VARNAMES[SIM_SELECT[var_index],])
         } # end title                    
         if (SIM_PLOT_RAW) {
		matplot(time_axis[1:min(SIM_LENGTH,SIM_MAX)],
			  sim_raw[SIM_SELECT[var_index],1:min(SIM_LENGTH,SIM_MAX)],
			  type='l',lwd=1,main=title,xlab='Year',
			  ylab='Percent dev. from SS')
            # text(time_axis(SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
            #    VARNAMES(SIM_SELECT(var_index),:));
         } else {
		matplot(time_axis[1:min(SIM_LENGTH,SIM_MAX)],
			  sim_xyz[SIM_SELECT[var_index],1:min(SIM_LENGTH,SIM_MAX)],
			  type='l',lwd=1,main=title,xlab='Year',
			  ylab='Percent dev. from SS')

            # text(time_axis(SIM_TXT_MARKER), sim_xyz(var_index,SIM_TXT_MARKER),...
            #    VARNAMES(SIM_SELECT(var_index),:));
         }
         #grid;
      } # end for varindex
         dummy<-readline('Inspect figure. Hit ENTER when ready...')
         par(op)
  } # end if SIM_SUBPLOT

#---------------------------------------------END SIM SUBPLOT------------------------------------------------------------

#---------------------------------------------SINGLE PLOTS---------------------------------------------------------------
      
   if (SIM_SINGLE) {
      for (var_index in 1:n_select) {
         if (DO_HP_FILTER & !SIM_PLOT_RAW) {
            title=paste('Simulated data (HP-filtered):',VARNAMES[SIM_SELECT[var_index],])
         } else {
            title=paste('Simulated data:',VARNAMES[SIM_SELECT[var_index],])
         } # end title
         if (SIM_PLOT_RAW) {
		matplot(time_axis[1:min(SIM_LENGTH,SIM_MAX)],
		        sim_raw[SIM_SELECT[var_index],1:min(SIM_LENGTH,SIM_MAX)],
			  type='l',lwd=2,lty=1,xlab='Year',main=title,
			  ylab='Percent deviation from steady state')
            # text(time_axis(SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
            #    VARNAMES(SIM_SELECT(var_index),:));
         } else {
		matplot(time_axis[1:min(SIM_LENGTH,SIM_MAX)],
		        sim_xyz[SIM_SELECT[var_index],1:min(SIM_LENGTH,SIM_MAX)],
			  type='l',lwd=2,lty=1,xlab='Year',main=title,
			  ylab='Percent deviation from steady state')
            # text(time_axis(SIM_TXT_MARKER), sim_xyz(var_index,SIM_TXT_MARKER),...
            #    VARNAMES(SIM_SELECT(var_index),:));
         }
         grid()
            dummy<-readline('Inspect figure. Hit ENTER when ready...')            
        } # end for varindex
   } # end if sim single
   
   #--------------------------------------------------END SIM SINGLE----------------------------------------------------------
} #end If SIM_GRAPH

cat('\n')
if (SIM_DO_DISP1 | SIM_DO_DISP2 | SIM_DO_DISP3) {
   cat(        'SIM_OUT.M: Simulation-based calculation of moments\n')
   cat(sprintf('           Simulation length = %10d \n',SIM_LENGTH))
   if (SIM_MODE == 2) {
      cat(sprintf('           repeated %10d times.\n',SIM_N_SERIES))
   }      
   cat('The variables are:\n')
   print(matrix(VARNAMES[SIM_SELECT,]),quote=F);
   cat('\n ')
}
if (SIM_DO_DISP1) {
   if (DO_HP_FILTER) {
      cat('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j\n')
      cat('(Simulation-based calculations)\n')
   } else {
      cat('Autocorrelation Table, corr(v(t+j),GNP(t)).  Last row shows j\n')
      cat('(Simulation-based calculations)\n')
   }
   for (var_index in 1:n_select) {
      cat(sprintf('  %5.2f',autcor_sim[var_index,]))
	cat('\n')
   }
   cat(sprintf('  %5.0f',autcor_sim[n_select+1,]))
   cat('\n');
   if ((SIM_N_SERIES > 3) & (SIM_MODE == 2)) {
      cat('Small sample standard errors for the Autocorrelation Table:\n')
      cat('(Simulation-based calculations)\n')
      for (var_index in 1 : n_select) {
         cat(sprintf('  %5.2f',autcor_std[var_index,]))
	   cat('\n')
      } 
      cat(sprintf('  %5.0f',autcor_sim[n_select+1,]))
      cat('\n')
   }
}
if (SIM_DO_DISP2) {
   if (DO_HP_FILTER) {
      cat('Variance-Covariance Matrix (HP-filtered series):\n')
      cat('(Simulation-based calculations)\n')
   } else {
      cat('Variance-Covariance Matrix:\n')
      cat('(Simulation-based calculations)\n')
   }
   for (var_index in 1 : n_select) {
      cat(sprintf(' %6.3f',covmat_sim[var_index,]))
	cat('\n')
   }
   cat('\n')
   if ((SIM_N_SERIES > 3) & (SIM_MODE == 2)) {
      cat('Small sample standard errors for the Variance-Covariance Matrix:\n')
      cat('(Simulation-based calculations)\n')
      for (var_index in 1 : n_select) {
         cat(sprintf(' %6.3f',covmat_std[var_index,]))
	   cat('\n')
      }
      cat('\n')
   }
} # end do disp2
if (SIM_DO_DISP3) {
   if (DO_HP_FILTER) {
      cat('Standard deviations (HP-filtered series):\n')
      cat('(Simulation-based calculations)\n')
   } else {
      cat('Standard deviations:\n')
      cat('(Simulation-based calculations)\n')
   }
   print(matrix(stdvec_sim))
   if ((SIM_N_SERIES > 3) & (SIM_MODE == 2)) {
      cat('Small sample standard errors for the Standard deviations:\n')
      cat('(Simulation-based calculations)\n')
      cat(matrix(stdvec_std))
   }
} # end do disp3
