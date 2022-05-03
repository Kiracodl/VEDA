#function rappture_reader_single( filename, suffix) 
#this file reads a single rappture file into python
#returns a diction object that maps curve names to data
#suffix is optional.
#example rappture_reader_single('run123456.xml')

import os
import csv
import pandas as pd
import xml.etree.ElementTree as ET


def rappture_reader_single( filename, suffix=None) :


        # for very large exponents, such as 1.234E-200, fortran likes to drop
        # the E to save one character, so it gets 1.2345-200. no programs other than fortran
        # understand this, and I don't know how to make fortran stop this, so we have this hack
        # windows machines will not typically have sed installed. If you get errors,
        # you can comment this line out, (but if you have such large exponents in your
        # file you'll have to edit them manually)
       os.system('sed -e ''s/\(\.[0-9]*\)[0-9]\([-+]\)\([0-9]*\)/\1E\2\3/g'' ' + filename  + ' > tmp.out')

       root = ET.parse('tmp.out').getroot()    
    
       curves = root.iter('curve')        

       dflist = []
    
       #numCurves = xdoc.getElementsByTagName('curve').getLength;
       for curve in curves:        
            curvename = curve.get('id')
        
            print(curvename)

            xy  = curve.find('component/xy')

            if (xy != None):
              curvedata = xy.text
        
        
              df = pd.DataFrame([x.split() for x in curvedata.split('\n')])
              df.columns = [curvename + '_x', curvename + '_y']
            
              df[curvename + '_x'] = pd.to_numeric(df[curvename + '_x'],errors='coerce')
              df[curvename + '_y'] = pd.to_numeric(df[curvename + '_y'],errors='coerce')
            
              dflist.append(df)
       
       #d = pd.concat( dflist)  
       #d = pd.join(dflist) #see if this works better
       #d = dflist[0].join( dflist[1:])
       d = pd.concat( dflist, axis=1)  
       return d