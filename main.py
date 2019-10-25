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
auxFuncts.checkInput(inp)

pts=np.genfromtxt(inp['filename'], skip_header=1, delimiter='\t')
if(inp['dim']==2):
    sfR2D=SFRegion2D(pts)
    print('2D region created')
    print('Q parameter',sfR2D.q)
    if(inp['eps']=='None'):
        assert sfR2D.q <0.7, "to search for small significant structures the region must be structured"
        sfR2D.calculateRho()
        sfR2D.calculateEps()
        print(sfR2D.eps)
        sfR2D.calculateNmin()
        print(sfR2D.Nmin)
        sfR2D.detectStructs()
        npts=pts.shape[0]
        salid = np.hstack((pts,sfR2D.db.labels_.reshape(npts,1)))
        np.savetxt(ofil,salid,header='ra  dec  cluster')
else:
    print('Now, only 2 dimensions')



