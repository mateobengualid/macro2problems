# VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
# DO_IT.M performs all the calculations and calls output-
# creating routines.  It assumes, that all required variables have been set

# Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
# However, you are not allowed to sell this software or otherwise impinge
# on its free distribution.


l_equ = dim(AA)[1]
m_states = dim(AA)[2]
n_endog=dim(CC)[1]
k_exog=dim(DD)[2]

message = '                                                                       ';
Warnings = NULL

source("options.R")
source("solve.R")
source("sol_out.R")
if (DISPLAY_LATER & (max(dim(Warnings)) > 1)) {
   cat('======================================================================= \n');
   cat('Your messages: (You can turn me off with DISPLAY_LATER = 0) \n');
   print(Warnings,quote=F);
   cat('======================================================================= \n');
}
if (DO_IMPRESP) {
   source("impresp.R")
}
if (DO_SIMUL) {
   source("simul.R")
   source("sim_out.R")
}
if (DO_MOMENTS) {
   source("moments.R")
   source("mom_out.R")
}
if (DISPLAY_AT_THE_END & (max(dim(Warnings)) > 1)) {
   cat('======================================================================= \n')
   if (DISPLAY_LATER | DISPLAY_IMMEDIATELY) {
      cat('Again, your (warning) messages: (You can turn me off with DISPLAY_AT_THE_END = 0) \n') }
   else {
      cat('Your messages: (You can turn me off with DISPLAY_AT_THE_END = 0) \n')
	}
   print(Warnings,quote=F)
   cat('======================================================================= \n')
} else
   cat('(Note: Messages will be displayed here with DISPLAY_AT_THE_END = 1) \n')
