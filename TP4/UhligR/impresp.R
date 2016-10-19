# VERSION 4.0, November 2002, COPYRIGHT H. UHLIG.
# IMPRESP.M calculates and plots impulse responses to shocks one percent
# in size to each one of the exogenous variables z(t), keeping all
# others at zero in the first period.  
# It also calculates impulse responses to deviations of state
# variables from steady state of one percent in size.
# It is assumed,
# that SOLVE.M has been executed before, so that the matrices
# NN, PP, QQ, RR and SS are available, describing the law of motion
#   x(t) = PP x(t-1) + QQ z(t)
#   y(t) = RR x(t-1) + SS z(t)
#   z(t) = NN z(t-1) + epsilon(t)
# The following options should have been chosen beforehand:
#   default values are declared in OPTIONS.M:
#   INIT_DATE  : number of periods prior to shock.  Default = 0.
#                For INIT_DATE = 1, one sees the inital value for the state variable.
#                For INIT_DATE > 0, variables are initially zero, before responding.
#   HORIZON    : how far out should the impulse responses be calculated
#   PERIOD     : number of periods per year, i.e. 12 for monthly, 4 for quarterly
#   DO_PLOTS   : = 1, if plots should be made, = 0, if not.
#   IMP_SELECT : a vector containing the indices of the variables to be
#             plotted.  To plot all variables, set IMP_SELECT = 1 : (m+n+k),
#             where m=dim(x), n=dim(y) and k=dim(z);
#   VARNAMES   : an array with (m+n+k) rows, containing the variable names.
#   DO_SHOCK_RESP: = 1, if responses to shocks should be plotted
#   DO_STATE_RESP: = 1, if responses to state variable deviations should be plotted.
#   The following options are set in OPTIONS.M and can be changed there:
#   TXT_MARKER : a number indicating where the label for the responses should be written.
#   DO_ENLARGE : = 1, if you want large font sizes for the text on your plots.  Good for slides.
#   PRINT_FIG  : = 1, if you want plots to be printed on your printer
#   SAVE_FIG   : = 1, if you want plots to be saved as encaps. postscript.  Set PRINT_FIG = 0 also.
#   SAVE_FIG_JPG : = 1, if you want plots to be saved as .jpg-files.  Set PRINT_FIG = 0 also.
# The program calculates:
# Response: the response of all variables x(t), y(t), z(t) to each shock.
#   Response is of size (m_states+n_endog+k_exog)*HORIZON.
#   Since Response is overwritten, each time a new shock is analyzed, 
#   the results are collected in 
# Resp_mat = [ Response to first shock
#              Response to second shock
#              ...                     ]
#
# The program also defines IMP_SUB_SELECT,
# Time_axis,m_states,n_endog,k_exog,II_contemp,II_lag,hndl,init_zero_response,
# thus overwriting variables with these names that might have been used before.
# 

# Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
# However, you are not allowed to sell this software or otherwise impinge
# on its free distribution.

# Most options are now set in OPTIONS.M

# The default option for plots is IMP_JOINT=1 

 
 

 
    

# Calculations
m_states = dim(QQ)[1]
n_endog = dim(SS)[1]
k_exog = dim(SS)[2]

Time_axis = matrix((-INIT_DATE : (HORIZON-1))/PERIOD)
Resp_mat <- NULL
op <- par(no.readonly=T)

#------------------------------------SHOCK RESP----------------------------------------------------------------
if (DO_SHOCK_RESP) {
   cat('Impulse responses to shocks... \n')
   for (shock_counter in SELECT_SHOCKS) { # 1 : k_exog,
      Response = matrix(0,nrow=m_states+n_endog+k_exog,ncol=HORIZON)
      Response[m_states+n_endog+shock_counter,1] = 1
      II_lag = rbind(cbind(PP, matrix(0,nrow=m_states,ncol=n_endog),matrix(0,m_states,k_exog)),
                 cbind(RR, matrix(0,n_endog, n_endog),matrix(0,n_endog, k_exog)),
                 cbind(matrix(0,k_exog,(m_states+n_endog)), NN)                )
      II_contemp = diag(1,m_states+n_endog+k_exog) + 
           rbind(cbind( matrix(0,m_states,(m_states+n_endog)), QQ),
             cbind(matrix(0,n_endog, (m_states+n_endog)), SS),
             cbind(matrix(0,k_exog,  (m_states+n_endog)), matrix(0,k_exog,k_exog) ))
      # describing [x(t)',y(t)',z(t)']'= II_contemp*II_lag*[x(t-1)',y(t-1)',z(t-1)']';
      Response[,1] = II_contemp%*%Response[,1]
      for (time_counter in 2 : HORIZON) {
         Response[,time_counter] = II_contemp%*%II_lag%*%Response[,(time_counter-1)]
      }
      Resp_mat = rbind( Resp_mat, 
                   Response )
      init_zero_response = matrix(0,m_states+n_endog+k_exog,INIT_DATE)
      
      IMP_SUB_SELECT <- NULL
      if (DO_PLOTS) {
         if (DO_NO_ZERO_RESPONSE) {
             for (varindex in IMP_SELECT) {
                 if (max(abs(matrix(c(init_zero_response[varindex,],Response[varindex,])))) > ZERO_RESPONSE_LEVEL) {
                     IMP_SUB_SELECT = matrix(t(c( IMP_SUB_SELECT, varindex )))
                 }
             }
         } else {
             IMP_SUB_SELECT = t(matrix(IMP_SELECT))
         }
                 
         
         if (is.null(IMP_SUB_SELECT)) {
             IMP_SUB_SELECT = matrix(c(GNP_INDEX))
         }  
         
#----------------------------------------JOINT PLOTS------------------------------------------------------------------
     if (IMP_JOINT) {
         #hndl = plot(Time_axis,0*Time_axis, ...
         #            Time_axis,[init_zero_response(IMP_SUB_SELECT,:),Response(IMP_SUB_SELECT,:)])
         plotmat=cbind(init_zero_response[IMP_SUB_SELECT,],Response[IMP_SUB_SELECT,])
         title=paste('Impulse responses to a shock in',
             VARNAMES[m_states+n_endog+shock_counter,])
         matplot(Time_axis,t(plotmat),type='l',
		xlab='Years after shock',
		ylab='Percent deviation from steady state',lwd=2,lty=1,
		main=title)
	   abline(h=0)
         if (DO_SCALE_AXIS) {
             if (DO_GIVEN_AXIS_SCALE) {
                 matplot(Time_axis,t(plotmat),type='l',
	           xlab='Years after shock',
		     ylab='Percent deviation from steady state',
		     xlim=AXIS_SCALE[1:2],
		     ylim=AXIS_SCALE[3:4],lwd=2,lty=1,main=title)
	           abline(h=0)
             } else {
                 matplot(Time_axis,t(plotmat),type='l',
	           xlab='Years after shock',
		     ylab='Percent deviation from steady state',
		     xlim=c(-INIT_DATE/PERIOD,HORIZON/PERIOD),lwd=2,lty=1,
		     ylim=c(min(0,min(min(Response[AXIS_SELECT,]))),max(1,max(max(Response[AXIS_SELECT,])))),
		     main=title)
	           abline(h=0)
             } # end if do given axis scale
         } # end if do scale axis
         grid()
         for (varindex in IMP_SUB_SELECT) {
             if (DO_FIXED_TEXT_POS) {
                 text(Time_axis[min(c((TXT_MARKER+INIT_DATE),(HORIZON/2+INIT_DATE)))], 
                     Response[varindex,min(c(TXT_MARKER,HORIZON/2))],
                     VARNAMES[varindex,])
             } else {
                 max_resp = max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
		     max_ind = which.max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
                 text(Time_axis[INIT_DATE+max_ind], 
                     Response[varindex,max_ind],
                     VARNAMES[varindex,])  
             } # end do fitted text pos           
         } # end for varindex

#            if (PRINT_FIG) {
#            print(rbind('IMPRESP.M: Printing response to a shock in ',VARNAMES[m_states+n_endog+shock_counter,],'...'))
#            if (DO_COLOR_PRINT) {
#               print -dwinc
#            } else {
#               print
#            }       
#        elseif (SAVE_FIG | SAVE_FIG_JPG )
#            if SAVE_FIG,
#                imp_filename=['impresp_joint_',FILENAMEPIECE,sprintf('#1.0f',shock_counter),'.eps'];
#                disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
#                disp(['         as encapsulated postscript file. Filename is ',imp_filename]);
#                 if COLOR_FIG 
#                    print('-depsc',imp_filename);
#                else
#                    print('-deps',imp_filename);
#                end;
#            end;
#            if SAVE_FIG_JPG,
#                imp_filename=['impresp_joint_',FILENAMEPIECE,sprintf('#1.0f',shock_counter),'.jpg'];
#                disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
#                disp(['         as jpg file. Filename is ',imp_filename]);
#                print('-djpeg90',imp_filename);
#            end;
#       end;                   
            dummy <- readline('Inspect figure. Hit ENTER when ready...');
            } # end if IMP_JOINT 
           
         
 #----------------------------------------END JOINT PLOTS-------------------------------------------------------------
 
 #----------------------------------------IMP SUBPLOT-----------------------------------------------------------------
                
      if (IMP_SUBPLOT) {
          
          #IMP_SUB_SELECT =[];#It is necessary to sort out zero respones because they will lead to an error in axis([..]) below
          #   for varindex = IMP_SELECT,
          #       if max(abs([init_zero_response(varindex,:),Response(varindex,:)])) > ZERO_RESPONSE_LEVEL,
          #           IMP_SUB_SELECT = [ IMP_SUB_SELECT, varindex ];
          #       end;
          #    end;
          size_plot_x=ceiling(sqrt(length(IMP_SUB_SELECT)))
          size_plot_y=round(sqrt(length(IMP_SUB_SELECT)))
          par(mfrow=c(size_plot_x,size_plot_y))

           for (step in 1:length(IMP_SUB_SELECT)) {
             axis=matrix(c(-INIT_DATE/PERIOD,HORIZON/PERIOD,min(0,min(min(Response[IMP_SUB_SELECT[1,step],]))),
                     max(0.1,max(max(c(init_zero_response[IMP_SUB_SELECT[1,step],],Response[IMP_SUB_SELECT[1,step],]))))))
             plotmat=matrix(c(init_zero_response[IMP_SUB_SELECT[1,step],],Response[IMP_SUB_SELECT[1,step],]))
             title=paste('Response to ', 
             VARNAMES[m_states+n_endog+shock_counter,])
             matplot(Time_axis,plotmat,type='l',
		         xlab='Years',
		         ylab='% deviation from SS',lwd=1,
			   xlim=axis[1:2],
                     ylim=axis[3:4],
		         main=title)
	       abline(h=0)
             
             if (DO_GIVEN_AXIS_SCALE) {
             matplot(Time_axis,plotmat,type='l',
		         xlab='Years',
		         ylab='% deviation from SS',lwd=1,
			   xlim=AXIS_SCALE[1:2],
                     ylim=AXIS_SCALE[3:4],
		         main=title)
	       abline(h=0)            
		 } # end if do given axis scale          
         
         # set(gca,'FontSize',IMP_SUB_FONT);   
         #grid;
         
# change from Uhlig code here to 'max', to fit text in first argument of text()    
         for (varindex in IMP_SUB_SELECT) {
                     text(Time_axis[max(c((TXT_MARKER+INIT_DATE),(HORIZON/2+INIT_DATE)))],
                     Response[IMP_SUB_SELECT[1,step],min(c(TXT_MARKER,HORIZON/2))],
                     VARNAMES[IMP_SUB_SELECT[1,step],])
                 
         } # end for varindex
     } # end for step  

#      if PRINT_FIG,
#            disp(['IMPRESP.M: Printing response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:),'...']);
#            if DO_COLOR_PRINT,
#               print -dwinc
#            else
#               print;
#            end;           
#        elseif (SAVE_FIG | SAVE_FIG_JPG )
#            if SAVE_FIG,
#                imp_filename=['impresp_subplot_',FILENAMEPIECE,sprintf('#1.0f',shock_counter),'.eps'];
#                disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
#                disp(['         as encapsulated postscript file. Filename is ',imp_filename]);
#                if COLOR_FIG 
#                    print('-depsc',imp_filename);
#                else
#                    print('-deps',imp_filename);
#                end;
#            end;
#            if SAVE_FIG_JPG,
#                imp_filename=['stateresp_subplot_',FILENAMEPIECE,sprintf('#1.0f',shock_counter),'.jpg'];
#                disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
#                disp(['         as jpg file. Filename is ',imp_filename]);
#                print('-djpeg90',imp_filename);
#            end;
#        end;                   
     
       dummy<-readline('Inspect figure. Hit ENTER when ready...')
	 par(op)
    } #if IMP_SUBPLOT
#----------------------------------------------END SUBPLOT--------------------------------------------------------------

#----------------------------------------------SINGLE PLOT--------------------------------------------------------------
 if (IMP_SINGLE) {
           for (step in 1:length(IMP_SUB_SELECT)) {
             plotmat <- matrix(c(init_zero_response[IMP_SUB_SELECT[1,step],],Response[IMP_SUB_SELECT[1,step],]))
             title=paste('Response to a one percent deviation in ',
                          VARNAMES[m_states+n_endog+shock_counter,])
		 matplot(Time_axis,plotmat,type='l',xlab='Years after shock',
                     ylab='Percent deviation from steady state',lwd=2,lty=1,main=title)		
		 abline(h=0)
         if (DO_SCALE_AXIS) {
             if (DO_GIVEN_AXIS_SCALE) {
	           title=paste('Response to a one percent deviation in ',
                              VARNAMES[m_states+n_endog+shock_counter,])
                 matplot(Time_axis,plotmat,type='l',xlim=AXIS_SCALE[1:2],
			       ylim=AXIS_SCALE[3:4],lwd=2,lty=1,
                         xlab='Years after shock',ylab='Percent deviation from steady state',
				 main=title)		
		     abline(h=0)
             } else {
		     title=paste('Response to a one percent deviation in',
                              VARNAMES[m_states+n_endog+shock_counter,])
		     matplot(Time_axis,plotmat,type='l',xlim=c(-INIT_DATE/PERIOD,HORIZON/PERIOD),
				 ylim=c(min(0,min(min(Response[IMP_SUB_SELECT[1,step],]))),
				 max(1,max(max(Response[IMP_SUB_SELECT[1,step],])))),lwd=2,lty=1,
				 xlab='Years after shock',ylab='Percent deviation from steady state',
				 main=title)
		     abline(h=0)
             }
         }
         grid()
         for (varindex in IMP_SUB_SELECT[1,step]) {
             if (DO_FIXED_TEXT_POS) {
                 text(Time_axis[min(c((TXT_MARKER+INIT_DATE),(HORIZON/2+INIT_DATE)))],
                     Response[varindex,min(c(TXT_MARKER,HORIZON/2))],
                     VARNAMES[varindex,])
             } else {
       	     max_resp = max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
		     max_ind = which.max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
                 text(Time_axis[INIT_DATE+max_ind], 
                     Response[varindex,max_ind],
                     VARNAMES[varindex,])  
             } # end do fitted text       
         } # end for varindex

        dummy<-readline('Inspect figure. Hit ENTER when ready... \n')   
      
     } # end for step        
    } # end if IMP_SINGLE   
 #---------------------------------------------END SINGLE PLOT----------------------------------------------------------

} #if DO_PLOTS
} # for shock_counter
} #if DO_SHOCK_RESP


#------------------------------------------STATE RESP------------------------------------------------------------------
if (DO_STATE_RESP) {
   cat('Impulse responses to deviation of state variables from steady state... \n')
   for (state_counter in SELECT_STATES) { # 1 : m_states,
      Response = matrix(0,m_states+n_endog+k_exog,HORIZON)
      Response[state_counter,1] = 1
      II_lag = rbind( cbind(PP, matrix(0,m_states,n_endog),matrix(0,m_states,k_exog)),
                 cbind(RR, matrix(0,n_endog, n_endog),matrix(0,n_endog, k_exog)),
                 cbind(matrix(0,k_exog,(m_states+n_endog)), NN                ))
      II_contemp = diag(1,m_states+n_endog+k_exog) + 
           rbind(cbind( matrix(0,m_states,(m_states+n_endog)), QQ),
             cbind(matrix(0,n_endog, (m_states+n_endog)), SS),
             cbind(matrix(0,k_exog,  (m_states+n_endog)), matrix(0,k_exog,k_exog)))
      # describing [x(t)',y(t)',z(t)']'= II_contemp*II_lag*[x(t-1)',y(t-1)',z(t-1)']';
      Response[,1] = II_contemp%*%II_lag%*%Response[,1]
      for (time_counter in 2 : HORIZON) {
         Response[,time_counter] = II_contemp%*%II_lag%*%Response[,(time_counter-1)]
      }
      Resp_mat = rbind( Resp_mat, 
                   Response )
      init_zero_response = matrix(0,m_states+n_endog+k_exog,INIT_DATE)
      if (INIT_DATE > 0) {
          init_zero_response[state_counter,INIT_DATE] = 1;
      }
    IMP_SUB_SELECT <- NULL
      if (DO_PLOTS) {
         if (DO_NO_ZERO_RESPONSE) {
             for (varindex in IMP_SELECT) {
                 if (max(abs(matrix(c(init_zero_response[varindex,],Response[varindex,])))) > ZERO_RESPONSE_LEVEL) {
                     IMP_SUB_SELECT = matrix(t(c( IMP_SUB_SELECT, varindex )))
                 }
             }
         } else {
             IMP_SUB_SELECT = t(matrix(IMP_SELECT))
         }
                 
         
         if (is.null(IMP_SUB_SELECT)) {
             IMP_SUB_SELECT = matrix(c(GNP_INDEX))
         }  
            
 #----------------------------------------STATE JOINT PLOTS------------------------------------------------------------------
     if (IMP_JOINT) {
         plotmat=cbind(init_zero_response[IMP_SUB_SELECT,],Response[IMP_SUB_SELECT,])
         title=paste('Impulse responses to a one percent deviation in',
             VARNAMES[state_counter,])
         matplot(Time_axis,t(plotmat),type='l',
		xlab='Years after shock',
		ylab='Percent deviation from steady state',lwd=2,lty=1,
		main=title)
	   abline(h=0)
         if (DO_SCALE_AXIS) {
             if (DO_GIVEN_AXIS_SCALE) {
                 matplot(Time_axis,t(plotmat),type='l',
	           xlab='Years after shock',
		     ylab='Percent deviation from steady state',
		     xlim=AXIS_SCALE[1:2],
		     ylim=AXIS_SCALE[3:4],lwd=2,lty=1,main=title)
	           abline(h=0)
             } else {
                 matplot(Time_axis,t(plotmat),type='l',
	           xlab='Years after shock',
		     ylab='Percent deviation from steady state',
		     xlim=c(-INIT_DATE/PERIOD,HORIZON/PERIOD),lwd=2,lty=1,
		     ylim=c(min(0,min(min(Response[AXIS_SELECT,]))),max(1,max(max(Response[AXIS_SELECT,])))),
		     main=title)
	           abline(h=0)
             } # end if do given axis scale
         } # end if do scale axis
         grid()
         for (varindex in IMP_SUB_SELECT) {
             if (DO_FIXED_TEXT_POS) {
                 text(Time_axis[min(c((TXT_MARKER+INIT_DATE),(HORIZON/2+INIT_DATE)))], 
                     Response[varindex,min(c(TXT_MARKER,HORIZON/2))],
                     VARNAMES[varindex,])
             } else {
                 max_resp = max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
		     max_ind = which.max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
                 text(Time_axis[INIT_DATE+max_ind], 
                     Response[varindex,max_ind],
                     VARNAMES[varindex,])  
             } # end do fitted text pos           
         } # end for varindex
          dummy<-readline('Inspect figure. Hit ENTER when ready... \n')  
     } # end if IMP_JOINT 
           
         
 #----------------------------------------END STATE JOINT PLOTS-------------------------------------------------------------
          
 #----------------------------------------STATE SUBPLOT-----------------------------------------------------------------
  
             
      if (IMP_SUBPLOT) {
          
          #IMP_SUB_SELECT =[]; #It is necessary to sort out zero respones because they will lead to an error in axis([..]) below
          #   for varindex = IMP_SELECT,
          #       if max(abs([init_zero_response(varindex,:),Response(varindex,:)])) > ZERO_RESPONSE_LEVEL,
          #           IMP_SUB_SELECT = [ IMP_SUB_SELECT, varindex ];
          #       end;
          #   end;
             
          size_plot_x=ceiling(sqrt(length(IMP_SUB_SELECT)))
          size_plot_y=round(sqrt(length(IMP_SUB_SELECT)))
          par(mfrow=c(size_plot_x,size_plot_y))

           for (step in 1:length(IMP_SUB_SELECT)) {
             axis=matrix(c(-INIT_DATE/PERIOD,HORIZON/PERIOD,min(0,min(min(Response[IMP_SUB_SELECT[1,step],]))),
                     max(0.1,max(max(c(init_zero_response[IMP_SUB_SELECT[1,step],],Response[IMP_SUB_SELECT[1,step],]))))))
             plotmat=matrix(c(init_zero_response[IMP_SUB_SELECT[1,step],],Response[IMP_SUB_SELECT[1,step],]))
             title=paste('Resp. to 1% dev. of ',
             VARNAMES[state_counter,])
             matplot(Time_axis,plotmat,type='l',
		         xlab='Years',
		         ylab='% deviation from SS',lwd=1,
			   xlim=axis[1:2],
                     ylim=axis[3:4],
		         main=title)
	       abline(h=0)
             
             if (DO_GIVEN_AXIS_SCALE) {
             matplot(Time_axis,plotmat,type='l',
		         xlab='Years',
		         ylab='% deviation from SS',lwd=1,
			   xlim=AXIS_SCALE[1:2],
                     ylim=AXIS_SCALE[3:4],
		         main=title)
	       abline(h=0)            
		 } # end if do given axis scale          
         #grid;                  
         for (varindex in IMP_SUB_SELECT) {
                     text(Time_axis[max(c((TXT_MARKER+INIT_DATE),(HORIZON/2+INIT_DATE)))],
                     Response[IMP_SUB_SELECT[1,step],min(c(TXT_MARKER,HORIZON/2))],
                     VARNAMES[IMP_SUB_SELECT[1,step],])
                 
         } # end for varindex
     } # end for step   
        
       dummy<-readline('Inspect figure. Hit ENTER when ready... \n')  
       par(op) 
    } # end if IMP_SUBPLOT
#----------------------------------------------END STATE SUBPLOT--------------------------------------------------------------

#----------------------------------------------STATE SINGLE PLOT--------------------------------------------------------------
 if (IMP_SINGLE) {
          for (step in 1:length(IMP_SUB_SELECT)) {
             plotmat <- matrix(c(init_zero_response[IMP_SUB_SELECT[1,step],],Response[IMP_SUB_SELECT[1,step],]))
             title=paste('Impulse responses to a one percent deviation in',
                          VARNAMES[state_counter,])
		 matplot(Time_axis,plotmat,type='l',xlab='Years after shock',
                     ylab='Percent deviation from steady state',lwd=2,lty=1,main=title)		
		 abline(h=0)
        if (DO_SCALE_AXIS) {
             if (DO_GIVEN_AXIS_SCALE) {
                 matplot(Time_axis,plotmat,type='l',xlim=AXIS_SCALE[1:2],
			       ylim=AXIS_SCALE[3:4],lwd=2,lty=1,
                         xlab='Years after shock',ylab='Percent deviation from steady state',
				 main=title)		
		     abline(h=0)
             } else {
		     matplot(Time_axis,plotmat,type='l',xlim=c(-INIT_DATE/PERIOD,HORIZON/PERIOD),
				 ylim=c(min(0,min(min(Response[IMP_SUB_SELECT[1,step],]))),
				 max(1,max(max(Response[IMP_SUB_SELECT[1,step],])))),lwd=2,lty=1,
				 xlab='Years after shock',ylab='Percent deviation from steady state',
				 main=title)
		     abline(h=0)
             }
         }
         grid()
         for (varindex in IMP_SUB_SELECT[1,step]) {
             if (DO_FIXED_TEXT_POS) {
                 text(Time_axis[min(c((TXT_MARKER+INIT_DATE),(HORIZON/2+INIT_DATE)))],
                     Response[varindex,min(c(TXT_MARKER,HORIZON/2))],
                     VARNAMES[varindex,])
             } else {
       	     max_resp = max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
		     max_ind = which.max(abs(Response[varindex,1:min(MAX_TXT_MARKER,HORIZON)]))
                 text(Time_axis[INIT_DATE+max_ind], 
                     Response[varindex,max_ind],
                     VARNAMES[varindex,])  
             } # end do fitted text       
         } # end for varindex
                 
        dummy<-readline('Inspect figure. Hit ENTER when ready... \n')       
      } # end for step        
    } # end if IMP_SINGLE   
 #---------------------------------------------END SINGLE PLOT--------------------------------------------------
                       
        } # end for state_counter
   } # end if DO_PLOTS
} # end f DO_STATE_RESP