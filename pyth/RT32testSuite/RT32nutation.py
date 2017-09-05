'''
Created on Jun 13, 2016

@author: blew
'''

import numpy as np

class RT32nutation():
    '''
    test possibilityo of periodic jumps 
    '''


    def __init__(self, jdSt, jdEn,RA,DEC):
        '''
        Constructor
        '''
        sec2rad=4.848136811095e-6
        deg2rad=0.0174532925199
        cosRA=np.cos(deg2rad*RA)
        sinRA=np.sin(deg2rad*RA)
        cosDEC=np.cos(deg2rad*DEC)

        self.dRa=list()
        self.dDec=list()
        
        for DJ in np.arange(jdSt,jdEn,0.1):
            ddj = DJ-2451545.0
            SL = (280.460+0.9856474*ddj)*deg2rad
            aM = (357.528+0.9856003*ddj)*deg2rad
            Om =(47.8-0.05295376*(DJ-2453004.5))*deg2rad
            
            sinE = np.sin((23.439-4.0-7.0*ddj)*deg2rad); 
            cosE = np.sqrt(1.0-sinE*sinE);
    
            dl2 = np.fmod(218.316+481267.881*ddj/36525.0,180.0)*np.pi/90.0;

            dpsi=(-17.2*np.sin(Om) - 1.319*np.sin(SL+SL) + 0.206*np.sin(2.0*Om) + 0.143*np.sin(aM) - 0.227*np.sin(dl2))*sec2rad
            deps=(+9.203*np.cos(Om) + 0.574*np.cos(SL+SL) + 0.098*np.cos(dl2) - 0.09*np.cos(2.0*Om))*sec2rad
    
            self.dRa.append( [DJ, ((cosE + sinE*sinRA*np.tan(DEC))*dpsi - cosRA*np.tan(DEC)*deps)/deg2rad ] )
            self.dDec.append( [DJ, (sinE*cosRA*dpsi + sinRA*deps)/deg2rad ] )
            
            
        np.savetxt("dRa.JD", self.dRa)
        np.savetxt("dDec.JD", self.dDec)
        
        
        