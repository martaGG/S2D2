from clases import SFRegion2D
import auxFuncts
import numpy as np 
import argparse

#load and initialization
parser = argparse.ArgumentParser()   
parser.add_argument('--fname', required=True)
args = parser.parse_args()
inp = auxFuncts.readInput(args.fname)
spl=args.fname.split('.')
#print(spl)
ofil='.'.join(spl[:-1])+'.out'
ofilplot='.'.join(spl[:-1])+'.pdf'
print(ofil)
auxFuncts.checkInput(inp)
headerOutput='x  y  cluster'
if(inp['coord']=='Ra Dec'):
    headerOutput='ra  dec  cluster'
pts=np.genfromtxt(inp['filename'], skip_header=1, delimiter='\t')
if(inp['dim']==2):
    sfR2D=SFRegion2D(pts, inp['coord'])
    print('2D region created')
    print('Q parameter',sfR2D.q)
    if(inp['eps']=='None'):
        assert sfR2D.q <inp['Qlim'], "Check input Q limit for structure searching."        
        sfR2D.calculateRho()
        print('RHO',sfR2D.rho)
        sfR2D.calculateEps()
        print('EPS',sfR2D.eps)
        sfR2D.setSignif(inp['Signif'])
        print('SIGNIF', sfR2D.signif)
        sfR2D.calculateNmin()
        print('NMIN',sfR2D.Nmin)
        sfR2D.detectStructs()
        npts=pts.shape[0]
        salid = np.hstack((pts,sfR2D.db.labels_.reshape(npts,1)))
        np.savetxt(ofil,salid,header=headerOutput)
        sfR2D.plotStructs(ofilplot)
    elif(inp['eps']>0):
        assert inp['Nmin']>0, "please input also a minimum number of points for dbscan"
        sfR2D.calculateRho()
        sfR2D.setEps(inp['eps'])
        sfR2D.setNmin(inp['Nmin'])
        print(sfR2D.eps)
        print(sfR2D.Nmin)
        sfR2D.detectStructs()
        npts=pts.shape[0]
        salid = np.hstack((pts,sfR2D.db.labels_.reshape(npts,1)))
        np.savetxt(ofil,salid,header=headerOutput)
        sfR2D.plotStructs(ofilplot)
    else:
        print('invalid value of eps, check documentation')
else:
    print('Now, only 2 dimensions')



