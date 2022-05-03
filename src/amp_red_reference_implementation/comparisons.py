#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:40:10 2020

@author: drknew
"""

import rappture_reader_single as rr

import matplotlib.pyplot as plot

import pandas as pd

import numpy as np

#%%
#phase 1 plots
plot.close('all')

case = 4;

fn1 = '/home/drknew/app-adac-AttardAmpRed/rappture/amp_red_avg/testcase' + str(case) + '.xml'
fn2 = '/home/drknew/app-adac-AttardAmpRed/rappture/dac/testcase' + str(case) + '.xml'


df1 = rappture_reader_single( fn1)
df2 = rappture_reader_single( fn2)

plot.figure()
plot.plot( df1['Amp_x'], df1['Amp_y'],  df2['Amp_x'], df2['Amp_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Amplitude (nm pk)')
plot.grid()
plot.title('First harmonic amplitude')
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_Amp.png')

plot.figure()
plot.plot( df1['Phase_x'], df1['Phase_y'],  df2['Phase_x'], df2['Phase_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Phase (Deg)')
plot.grid()
plot.title('First harmonic phase')
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_Phase.png')

plot.figure()
s= 'MeanForce'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title(s)
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_MeanForce.png')

plot.figure()
s= 'FPeakRep'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title('Peak Force Repulsive')
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_PeakForceRep.png')


plot.figure()
s= 'FPeakAtt'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title('Peak Force Attractive')
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_PeakForceAtt.png')

plot.figure()
s= 'virial'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Energy')
plot.grid()
plot.title(s)
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_virial.png')


plot.figure()
s= 'Indent'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Indentation (nm)')
plot.grid()
plot.title(s)
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_dissipation.png')


plot.figure()
s= 'tcontact'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Contact time (microseconds)')
plot.grid()
plot.title('Contact time')
plot.legend([ 'New tool', 'Original VEDA'])
plot.savefig( 'Case'  + str(case) + '_contactTime.png')


plot.figure()
s= 'Fts1'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Periods (t/T)')
plot.ylabel('Tip sample force (nN)')
plot.grid()
plot.title('Time history @ Asp=0.55, tip-sample force')
plot.legend([ 'New tool', 'Original VEDA'])
plot.xlim([0 ,1])
plot.savefig( 'Case'  + str(case) + '_TimeHistoryForce.png')

#%%

#phase 2 plots

plot.close('all')

#case = 5; mesh = 70; per=0.0000036; A0 = 75; nstep = 10e3; th = 1;
#case = 6; mesh = 100; per=1/100e3; A0 = 20; nstep =10e3; th = 1;
#case = 6.1; mesh = 20; per=1/100e3; A0 = 20; nstep = 150; th = 1;
case = 7; mesh = 50; per=1/100e3; A0 = 20; nstep =10e3; th = 1;
#case = 8; mesh = 50; per=1/100e3; A0 = 50; nstep =10e3; th = 2;

fn1 = '/home/drknew/app-adac-AttardAmpRed/rappture/amp_red_avg/testcase' +  str(case)  + 'A.xml'
df1A = rr.rappture_reader_single( fn1)

fn1 = '/home/drknew/app-adac-AttardAmpRed/rappture/amp_red_avg/testcase' +  str(case)  + 'B.xml'
df1B = rr.rappture_reader_single( fn1)


fn2 = '/home/drknew/OneDrive/research/case_' +  str(case) + '_matlab.csv'
#fn2 = '/home/drknew/OneDrive/research/case_' +  str(case) + 'spatial_matlab.csv'
#fn2 = '/home/drknew/app-adac-AttardAmpRed/src/amp_red_reference_implementation/case_' + str(case) + '_octave.csv'
df2 = pd.read_csv(fn2)

df2.columns = ['arstep', 'Ar(arstep)', 'Ar_guess', 'Zguess(arstep)',  'amp(arstep)', 'dissipation', 'virial', 'Phase']


#fn3 = '/home/drknew/app-adac-AttardAmpRed/src/amp_red_reference_implementation/case' + str(case) + '_th' + str(th)  + '_meshnum' + str(mesh) + '_nstep' + str(int(nstep)) + '_octave.csv'
#fn3 = '/home/drknew/OneDrive/research/case' + str(case) + '_th1_meshnum' + str(mesh) + '_nstep10000spatial_matlab.csv'
fn3 = '/home/drknew/OneDrive/research/case' + str(case) + '_th' + str(th) + '_meshnum' + str(mesh) + '_nstep10000_matlab.csv'
df_th_tr = pd.read_csv(fn3);
df_th_tr.columns = ['time', 'force', 'gap_to_undeform', 'dpr', 'surf_coord']

pk =1

np.count_nonzero(  df1A['Fts1_y'] )


#%%
plot.close('all')


plot.figure()
plot.plot( df2['Ar(arstep)'], df2['Ar_guess'], df2['Ar(arstep)'], df1A['Amp_iter_y'][0:9] / A0, df2['Ar(arstep)'], df1B['Amp_iter_y'][0:9] / A0, '+-')
#plot.plot( df2['Ar(arstep)'], df2['Ar_guess'], df2['Ar(arstep)'], df1A['Amp_iter_y'][0:9] / A0)
plot.xlabel('Target setpoint')
plot.ylabel('Output setpoint')
plot.grid()
#plot.title('Time history @ Asp=0.55, tip-sample force')
plot.legend([ 'Reference matlab', 'VEDA (backwards compat)', 'VEDA (new)'])
plot.xlim([0 ,1])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1



plot.figure()

plot.plot( df2['Ar(arstep)'], df2['Zguess(arstep)'], df2['Ar(arstep)'], df1A['Amp_iter_x'][0:9], '+-', df2['Ar(arstep)'], df1B['Amp_iter_x'][0:9], '+-' )
plot.xlabel('Target setpoint')
plot.ylabel('Z')
plot.grid()
#plot.title('Time history @ Asp=0.55, tip-sample force')
plot.legend([ 'Reference matlab', 'VEDA (backwards compat)', 'VEDA (new)'])
plot.xlim([0 ,1])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1




plot.figure()

plot.plot( df2['Zguess(arstep)'], -df2['Ar(arstep)']+df2['Ar_guess'], '+-',  df1A['Amp_x'][0:9] ,  -df2['Ar(arstep)']+df1A['Amp_iter_y'][0:9] / A0, '+-',  df1B['Amp_x'][0:9] ,  -df2['Ar(arstep)']+df1B['Amp_iter_y'][0:9] / A0, '+-')
mz = max(df2['Zguess(arstep)'])
if (case == 7) or (case == 8) :
  tol = 0.001
else:
  tol = 0.01
  
plot.plot( [0, mz], [tol, tol], 'k--', [0, mz], [-tol, -tol], 'k--')
plot.xlabel('Z distance')
plot.ylabel('Setpoint error')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1




plot.figure()
plot.plot( df2['Zguess(arstep)'], A0*df2['Ar_guess'],  df1A['Amp_x'][0:9] , df1A['Amp_iter_y'][0:9],  df1B['Amp_x'][0:9] , df1B['Amp_iter_y'][0:9] )
plot.xlabel('Z distance')
plot.ylabel('Amplitude')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['Phase'],'--', df1A['Phase_x'] , df1A['Phase_y'], df1B['Phase_x'] , df1B['Phase_y'] )
plot.xlabel('Z distance')
plot.ylabel('Phase')
plot.grid()
#plot.title('Time history @ Asp=0.55, tip-sample force')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['dissipation'],'--', df1A['Pts_x'] , df1A['Pts_y'], df1B['Pts_x'] , df1B['Pts_y'] )
plot.xlabel('Z distance')
plot.ylabel('Dissipation')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df2['Zguess(arstep)'], -df2['virial'],'--', df1A['Pts_x'] , df1A['virial_y'], df1B['Pts_x'] , df1B['virial_y'] )
plot.xlabel('Z distance')
plot.ylabel('Virial')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1




plot.figure()
plot.plot( df_th_tr['time'] / per, df_th_tr['force']/1e-9,  df1A['Fts' + str(th) + '_x'] , df1A['Fts' + str(th) + '_y'], '--',  df1B['Fts' + str(th) + '_x'] , df1B['Fts' + str(th) + '_y'], '--')
plot.xlabel('Periods (t/T)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ' tip-sample force')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1



plot.figure()
#plot.plot( df_th_rc['gap_to_undeform']/1e-9 , df_th_rc['force']/1e-9, df_th_tr['gap_to_undeform']/1e-9 , df_th_tr['force']/1e-9,df1['gap1_y'] , df1['Fts1_y'])
plot.plot(  df_th_tr['gap_to_undeform']/1e-9 , df_th_tr['force']/1e-9,df1A['gap' +str(th) + '_y'] , df1A['Fts' + str(th) + '_y'], df1B['gap' +str(th) + '_y'] , df1B['Fts' + str(th) + '_y'])
plot.xlabel('h_0 (nm)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', tip-sample force')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1



plot.figure()
plot.plot( df_th_tr['time'], df_th_tr['gap_to_undeform']/1e-9, df1A['Fts' + str(th) + '_x']*per , df1A['gap' + str(th) + '_y'],  df1B['Fts' + str(th) + '_x']*per , df1B['gap' + str(th) + '_y'])
plot.xlabel('time')
plot.ylabel('h_0 (nm)')
plot.grid()
plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', tip-sample gap h_0')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1



plot.figure()
plot.plot( df_th_tr['time'], df_th_tr['dpr'], df1A['Fts' +str(th) + '_x']*per , df1A['Pspace' +str(th) +'_y']*1e-6, '--', df1B['Fts' +str(th) + '_x']*per , df1B['Pspace' +str(th) +'_y']*1e-6, '--')
plot.xlabel('time')
plot.ylabel('h_0_dot (m/s)')
plot.grid()
plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', tip-sample velocity h_0_dot')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df_th_tr['time'] / per, df_th_tr['surf_coord']/1e-9,  df1A['surf' +str(th) +'_x'] , df1A['surf' +str(th) +'_y'],  df1B['surf' +str(th) +'_x'] , df1B['surf' +str(th) +'_y'])
plot.xlabel('Periods (t/T)')
plot.ylabel('u (nm)')
plot.grid()
plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', u(r_1) = Deflection at r=0')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df_th_tr['time'] / per, df_th_tr['gap_to_undeform']/1e-9 - df_th_tr['surf_coord']/1e-9,  df1A['surf' +str(th) +'_x'] ,df1A['gap' +str(th) +'_y'] - df1A['surf' +str(th) +'_y'],  df1B['surf' +str(th) +'_x'] ,df1B['gap' +str(th) +'_y'] - df1B['surf' +str(th) +'_y'])
plot.xlabel('Periods (t/T)')
plot.ylabel('deformed gap h (nm)')
plot.grid()
plot.title('Time history @ ' + str(df2['Ar(arstep)'][th-1]) )
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1








plot.figure()
plot.plot(df_th_tr['time'] / per , df_th_tr['surf_coord']/1e-9,  df1A['surf' + str(th) + '_x'] , df1A['surf' + str(th) + '_y'],  df1B['surf' + str(th) + '_x'] , df1B['surf' + str(th) + '_y'])
plot.xlabel('Periods (t/T)')
plot.ylabel('u (nm)')
plot.grid()
plot.title('Time history @ Asp~' + str(df2['Ar(arstep)'][th-1])  + ', surface coord u(r_1)')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1



plot.figure()
plot.plot(  df_th_tr['surf_coord']/1e-9, df_th_tr['force']/1e-9, df1A['surf' + str(th) + '_y'],  df1A['Fts' + str(th) + '_y'],'--', df1B['surf' + str(th) + '_y'],  df1B['Fts' + str(th) + '_y'])
plot.xlabel('surf coord)')
plot.ylabel('force')
plot.grid()
plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', surface coord u(r_1)')
plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


#%%


#phase 3A plots, case 9 & 10 (new tool vs ADAC tool)
plot.close('all')

case = 10;

fn1 = '/home/drknew/app-adac-AttardAmpRedBimodal/rappture/amp_red_avg/testcase' + str(case) + '.xml'
fn2 = '/home/drknew/app-adac-AttardAmpRedBimodal/rappture/adac/testcase' + str(case) + '.xml'


df2 = rr.rappture_reader_single( fn1)
df1 = rr.rappture_reader_single( fn2)

#%%

plot.figure()
plot.plot( df1['Amp_x'], df1['Amp_y'],  df2['Amp_x'], df2['Amp_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Amplitude (nm pk)')
plot.grid()
plot.title('First harmonic amplitude')
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_Amp.png')

plot.figure()
plot.plot( df1['Phase_x'], df1['Phase_y'],  df2['Phase_x'], df2['Phase_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Phase (Deg)')
plot.grid()
plot.title('First harmonic phase')
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_Phase.png')


plot.figure()
plot.plot( df1['Amp2_x'], df1['Amp2_y'],  df2['Amp2_x'], df2['Amp2_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Amplitude 2 (nm pk)')
plot.grid()
plot.title('2nd frequency amplitude')
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_Amp2.png')

plot.figure()
plot.plot( df1['Phase2_x'], df1['Phase2_y'],  df2['Phase2_x'], df2['Phase2_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Phase 2 (Deg)')
plot.grid()
plot.title('Second frequency phase')
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_Phase2.png')



plot.figure()
s= 'MeanForce'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title(s)
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_MeanForce.png')

plot.figure()
s= 'FPeakRep'
plot.plot( df1[s + '_x'], df1[s + '_y'], '+-', df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title('Peak Force Repulsive')
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_PeakForceRep.png')


plot.figure()
s= 'FPeakAtt'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Force (nN)')
plot.grid()
plot.title('Peak Force Attractive')
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_PeakForceAtt.png')

plot.figure()
s= 'virial'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Energy')
plot.grid()
plot.title(s)
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_virial.png')


plot.figure()
s= 'Pts'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Energy')
plot.grid()
plot.title(s)
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_dissipation.png')


plot.figure()
s= 'Indent'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Indentation (nm)')
plot.grid()
plot.title(s)
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_indentation.png')


plot.figure()
s= 'tcontact'
plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
plot.xlabel('Z distance (nm)')
plot.ylabel('Contact time (microseconds)')
plot.grid()
plot.title('Contact time')
plot.legend([ 'Original VEDA','New tool'])
plot.savefig( 'Case'  + str(case) + '_contactTime.png')

#
#plot.figure()
#s= 'Fts1'
#plot.plot( df1[s + '_x'], df1[s + '_y'],  df2[s + '_x'], df2[s + '_y'], '--'  )
#plot.xlabel('Periods (t/T)')
#plot.ylabel('Tip sample force (nN)')
#plot.grid()
#plot.title('Time history @ Asp=0.55, tip-sample force')
#plot.legend([ 'Original VEDA','New tool'])
#plot.xlim([0 ,1])
#plot.savefig( 'Case'  + str(case) + '_TimeHistoryForce.png')
#

#%%



#phase 3 plots 

plot.close('all')

#case = 11; A0 = 75;per=0.0000036; Af2=3.5; th = 4
case = 12; A0 = 50;per=1/(100e3); Af2=3.5; th = 0;


fn1 = '/home/drknew/app-adac-AttardAmpRedBimodal/rappture/amp_red_avg/testcase' +  str(case)  + '.xml'
df1A = rr.rappture_reader_single( fn1)


fn2 = '/home/drknew/OneDrive/research/case_' +  str(case) + '_matlab.csv'
#fn2 = '/home/drknew/OneDrive/research/case_' +  str(case) + 'simplified_matlab.csv'
df2 = pd.read_csv(fn2)

df2.columns = ['arstep', 'Ar(arstep)', 'Ar_guess', 'Zguess(arstep)',  'Phase1', 'A2', 'Asp2_err', 'Phase2',  'dissipation2', 'virial2', 'dissipation1', 'virial1'];

pk =1
if (th>0):
  #fn3 = '/home/drknew/app-adac-AttardAmpRedBimodal/src/amp_red_reference_implementation/case' +  str(case)  + '_th' +str(th) + '.csv'
  fn3 = '/home/drknew/OneDrive/research/case' +  str(case)  + '_th' +str(th) + '.csv'
  df_th_tr = pd.read_csv(fn3);
  df_th_tr.columns = ['time', 'force', 'gap_to_undeform', 'dpr']
  
  
  
  np.count_nonzero(  df1A['Fts1_y'] )


#%%
plot.close('all')


plot.figure()
plot.plot( range(0,9), df2['Zguess(arstep)']-  df1A['Amp_iter_x'][0:9] ,'+-')
plot.xlabel('setpoint number')
plot.ylabel('difference in Z')
plot.grid()
#plot.legend([ 'Reference matlab', 'VEDA'])
#plot.xlim([0 ,1])
#plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1

#%%


plot.figure()
plot.plot( df2['Ar(arstep)'], df2['Ar_guess'], df1A['Amp_y'][0:9]/A0, df1A['Amp_iter_y'][0:9] / A0 ,'+-')
plot.xlabel('Target setpoint')
plot.ylabel('Output setpoint')
plot.grid()
#plot.title('Time history @ Asp=0.55, tip-sample force')
plot.legend([ 'Reference matlab', 'VEDA'])
#plot.xlim([0 ,1])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


#
#plot.figure()
#plot.plot( df2['Ar(arstep)'], df2['Zguess(arstep)'], df2['Ar(arstep)'], df1A['Amp_iter_x'][0:9], '+-', df2['Ar(arstep)'], df1B['Amp_iter_x'][0:9], '+-' )
#plot.xlabel('Target setpoint')
#plot.ylabel('Z')
#plot.grid()
##plot.title('Time history @ Asp=0.55, tip-sample force')
#plot.legend([ 'Reference matlab', 'VEDA (backwards compat)', 'VEDA (new)'])
#plot.xlim([0 ,1])
#plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
#pk=pk+1




plot.figure()
plot.plot( df2['Zguess(arstep)'], -df2['Ar(arstep)']+df2['Ar_guess'], '+-',  df1A['Amp_x'][0:9] ,  -df2['Ar(arstep)']+df1A['Amp_iter_y'][0:9] / A0, '+-')
mz = max(df2['Zguess(arstep)'])
if (case == 7) or (case == 8) :
  tol = 0.001
else:
  tol = 0.01
  
plot.plot( [0, mz], [tol, tol], 'k--', [0, mz], [-tol, -tol], 'k--')
plot.xlabel('Z distance')
plot.ylabel('Setpoint error')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA', 'tolerance'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1




plot.figure()
plot.plot( df2['Zguess(arstep)'], A0*df2['Ar_guess'],  df1A['Amp_x'][0:9] , df1A['Amp_iter_y'][0:9])
plot.xlabel('Z distance')
plot.ylabel('Amplitude 1')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1





plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['Asp2_err'],  df1A['Amp_x'][0:9] , df1A['Amp2_iter_y'][0:9]/Af2 -   df1A['Amp2_y'][0:9]/Af2  )
plot.plot( [0, mz], [tol, tol], 'k--', [0, mz], [-tol, -tol], 'k--')
plot.xlabel('Z distance')
plot.ylabel('Amplitude Reduction 2 error')
plot.title('Amplitude Reduction 2 error')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA', 'tolerance'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['A2'],  df1A['Amp2_x'] , df1A['Amp2_y'], '--' )
plot.xlabel('Z distance (nm)')
plot.ylabel('Amp 2 (nm)')
plot.grid()
plot.title('Second frequency Amp')
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 




plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['Phase1'],'--', df1A['Phase_x'] , df1A['Phase_y'] )
plot.xlabel('Z distance')
plot.ylabel('Phase1')
plot.grid()
#plot.title('Time history @ Asp=0.55, tip-sample force')
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1






plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['Phase2'],  df1A['Phase2_x'] , df1A['Phase2_y'], '--' )
plot.xlabel('Z distance (nm)')
plot.ylabel('Phase 2 (Deg)')
plot.grid()
plot.title('Second frequency phase')
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 




plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['dissipation1'],'--', df1A['Pts1_x'] , df1A['Pts1_y'] )
plot.xlabel('Z distance')
plot.ylabel('Dissipation 1')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df2['Zguess(arstep)'], -df2['virial1'],'--', df1A['Pts_x'] , df1A['virial1_y'] )
plot.xlabel('Z distance')
plot.ylabel('Virial 1')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1




plot.figure()
plot.plot( df2['Zguess(arstep)'], df2['dissipation2'],'--', df1A['Pts2_x'] , df1A['Pts2_y'] )
plot.xlabel('Z distance')
plot.ylabel('Dissipation 2')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1


plot.figure()
plot.plot( df2['Zguess(arstep)'], -df2['virial2'],'-+', df1A['virial2_x'] , df1A['virial2_y'], '+-' )
plot.xlabel('Z distance')
plot.ylabel('Virial 2')
plot.grid()
plot.legend([ 'Reference matlab', 'VEDA'])
plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
pk=pk+1



if (th>0):

  plot.figure()
  plot.plot( df_th_tr['time'] / per, df_th_tr['force']/1e-9,  df1A['Fts' + str(th) + '_x'] , df1A['Fts' + str(th) + '_y'], '--')
  plot.xlabel('Periods (t/T)')
  plot.ylabel('Force (nN)')
  plot.grid()
  plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ' tip-sample force')
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
  
  
  
  plot.figure()
  plot.plot(  df_th_tr['gap_to_undeform'] , df_th_tr['force']/1e-9,df1A['gap' +str(th) + '_y'] , df1A['Fts' + str(th) + '_y'])
  plot.xlabel('h_0 (nm)')
  plot.ylabel('Force (nN)')
  plot.grid()
  plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', tip-sample force')
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
  
  
  
  plot.figure()
  plot.plot( df_th_tr['time'], df_th_tr['gap_to_undeform'], df1A['Fts' + str(th) + '_x']*per , df1A['gap' + str(th) + '_y'])
  plot.xlabel('time')
  plot.ylabel('h_0 (nm)')
  plot.grid()
  plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', tip-sample gap h_0')
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
  
  
  
  plot.figure()
  plot.plot( df_th_tr['time'], df_th_tr['dpr']*1e-6, df1A['Fts' +str(th) + '_x']*per , df1A['Pspace' +str(th) +'_y']*1e-6, '--')
  #plot.plot( df_th_tr['time'], df_th_tr['dpr']*1e-6)
  #plot.plot( df1A['Fts' +str(th) + '_x']*per , df1A['Pspace' +str(th) +'_y']*1e-6, '--')
  plot.xlabel('time')
  plot.ylabel('h_0_dot (m/s)')
  plot.grid()
  plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', tip-sample velocity h_0_dot')
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
  
  
  plot.figure()
  plot.plot( df_th_tr['time'] / per, df_th_tr['surf_coord']/1e-9,  df1A['surf' +str(th) +'_x'] , df1A['surf' +str(th) +'_y'],  df1B['surf' +str(th) +'_x'] , df1B['surf' +str(th) +'_y'])
  plot.xlabel('Periods (t/T)')
  plot.ylabel('u (nm)')
  plot.grid()
  plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', u(r_1) = Deflection at r=0')
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
  
  
  plot.figure()
  plot.plot( df_th_tr['time'] / per, df_th_tr['gap_to_undeform']/1e-9 - df_th_tr['surf_coord']/1e-9,  df1A['surf' +str(th) +'_x'] ,df1A['gap' +str(th) +'_y'] - df1A['surf' +str(th) +'_y'],  df1B['surf' +str(th) +'_x'] ,df1B['gap' +str(th) +'_y'] - df1B['surf' +str(th) +'_y'])
  plot.xlabel('Periods (t/T)')
  plot.ylabel('deformed gap h (nm)')
  plot.grid()
  plot.title('Time history @ ' + str(df2['Ar(arstep)'][th-1]) )
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
  
  
  
  
  
  
  
  
  plot.figure()
  plot.plot(df_th_tr['time'] / per , df_th_tr['surf_coord']/1e-9,  df1A['surf' + str(th) + '_x'] , df1A['surf' + str(th) + '_y'],  df1B['surf' + str(th) + '_x'] , df1B['surf' + str(th) + '_y'])
  plot.xlabel('Periods (t/T)')
  plot.ylabel('u (nm)')
  plot.grid()
  plot.title('Time history @ Asp~' + str(df2['Ar(arstep)'][th-1])  + ', surface coord u(r_1)')
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
  
  
  
  plot.figure()
  plot.plot(  df_th_tr['surf_coord']/1e-9, df_th_tr['force']/1e-9, df1A['surf' + str(th) + '_y'],  df1A['Fts' + str(th) + '_y'],'--', df1B['surf' + str(th) + '_y'],  df1B['Fts' + str(th) + '_y'])
  plot.xlabel('surf coord)')
  plot.ylabel('force')
  plot.grid()
  plot.title('Time history @ Asp~=' + str(df2['Ar(arstep)'][th-1])  + ', surface coord u(r_1)')
  plot.legend([ 'Reference matlab', 'VEDA (backwards compatible)', 'VEDA (new)'])
  plot.savefig( 'Case'  + str(case) + '_' +str(pk) + '.png') 
  pk=pk+1
