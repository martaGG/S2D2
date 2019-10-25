import numpy as np
from astropy.coordinates import SkyCoord
from scipy import stats, special, integrate
import yaml

def readInput(fyaml):
    stream=open(fyaml, 'r')
    inp=yaml.safe_load(stream)
    return inp

def checkInput(inp):
    assert inp["dim"]== 2, "only can use Dimension2"
    if(inp['eps']!='None'):
        assert type(inp['eps']) is float, "Eps must be a float"
        assert type(inp['Nmin']) is int, "Must provide integer min number of points"

def distMatGcirc(points):
    npoints=len(points[:,0])
    distMat=np.zeros([npoints, npoints],dtype=float)
    for i in range(npoints):
        jotas=range(i+1, npoints)
        ci=SkyCoord(ra=points[i,0], dec=points[i,1], unit='deg')
        cj=SkyCoord(ra=points[jotas,0], dec=points[jotas,1], unit='deg')
        distMat[i,jotas]=ci.separation(cj)
        distMat[jotas,i]=distMat[i,jotas]
    return distMat

def MST(distMat):
    """
    Uses prim's algorithm to calculate the minimum spanning tree
    from a distance matrix
    """
    nedge = len(distMat[1,:])-1
    v1 = 0
    matAux = np.copy(distMat[v1,:])
    matAux[v1]=np.nan
    v2 = np.nanargmin(matAux)
    d = matAux[v2]
    tree = []
    tree.append([v1,v2,d])
    intree = [v1,v2]
    while (len(intree) <= nedge):
        notInTree=[i for i in np.arange(len(distMat[1,:])) if i not in intree]
        mndist=np.nanmin(distMat[intree,:][:,notInTree])
        wh = np.where(distMat==mndist)
        w1 = wh[0][0]
        ww1 = wh[1][0]
        if(w1 in intree):
            v1 = w1
        else:   
            v1 = ww1

        matAux = np.copy(distMat[v1,])
        matAux[intree] = np.nan
        v2 = np.nanargmin(matAux)
        d = matAux[v2]
        tree.append([v1,v2,d])
        intree.append(v2)
    #print(intree)
    return(np.array(tree))


def paramQ(distMat):
    s = np.nanmean(distMat)
    req = np.sqrt(1/np.pi)
    sbar = s/req
    npoints = len(distMat[:,1])
    ms = MST(distMat)
    mdMST = np.mean(ms[:,2])
    mbar=mdMST/(np.sqrt(1*npoints)/(npoints-1))
    q=mbar/sbar
    return q

def nearestNeighbor(distMat, n):
    """
    Calculates the distance distribution of the n-th nearest neigbour
    from a matrix distance.
    """
    indis = np.argsort(distMat, axis=0)[n,:]
    dists = distMat[np.arange(len(distMat[:,0])),indis]
    return(indis, dists)

def calcSilvermanFactor(distrib):
    """
    Calculate the Silverman factor for bandwith used in R to calculater Kernel densities. 
    Different from the one used by python gaussian_kde. 
    The gaussian_kde function's implementation (at least in 1D) multiplies by the std dev.
    So the Silv. Factor fed fo gaussian_kde must be divided by std.dev to obtain the same result.

    Input: distrib, vector sampling a distribution
    Output: std deviation of distrib
            Silverman Factor.

    """
    stdev=np.std(distrib)
    iqr=stats.iqr(distrib)
    A = min(stdev, iqr/ 1.34)
    h = 0.9 * A * len(distrib) ** (-1 / 5)
    return(stdev,h)

def densNNDim(k,dim,rho,r):
    """
    Density function (of variable r) of the nearest neighbour distance distibution from a uniform spatial point process for:
    k = k-th nearest neighbour
    dim = space of dimension dim
    rho = density of the spatial process
    r = variable
    """
    cd = np.pi**(dim/2)/special.gamma(dim/2+1)
    Vol_r = cd*r**dim
    dens = dim/r*(rho*Vol_r)**k/special.gamma(k)*np.exp(-rho*Vol_r)
    return(dens)

def integraDensNN(k,dim, rho, a, b):
    """
    Calculates the approximate integral of the densNNDim from a to b.
    densNNdim: density function fo 
    k: k-th nearest neighbour
    """
    [integ, err]=integrate.quad(lambda x: densNNDim(k,dim,rho,x), a,b)
    return(integ)


def integraRNDensNN(n,k,dim, rho, a, b):
    """
    Calculates the approximate integral of the r^n*densNNDim from a to b.
    densNNdim: density function fo 
    k: k-th nearest neighbour
    """
    [integ, err]=integrate.quad(lambda x: densNNDim(k,dim,rho,x)*x**n, a,b)
    return(integ)

