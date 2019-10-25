import numpy as np
from scipy import spatial, stats, integrate
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import stats as astrostats
from sklearn.cluster import DBSCAN
import auxFuncts
class StarFormingRegion:
    def __init__(self, puntos):
        """
        Check dimensions of matrix are appropriate, 
        initialization of attributes
        """
        if(not len(puntos.shape)==2):
            print('not valid dimension!!!!')
            return
        try:
            self.points=puntos.astype('float')
        except TypeError:
            print('not valid kind of points!!!')
            return
        self.distMat=self.calculateDist()
        self.q=self.calculateQ()
        self.bBox=bBoxParams(self.points)

    def calculateDist(self):
        """
        Calculates a matrix distance from a set of points
        """
        distMat=spatial.distance_matrix(self.points, self.points)
        return distMat

    def calculateQ(self):
        return auxFuncts.paramQ(self.distMat)


class bBoxParams:
    def __init__(self,points):
        self.minCoords=np.min(points,axis=0)
        self.maxCoords=np.max(points, axis=0)
        self.volEucl=np.prod(self.maxCoords-self.minCoords)
        

    
class SFRegion2D(StarFormingRegion):
    def __init__(self, puntos):
        super(SFRegion2D,self).__init__(puntos)
        self.distMat=self.calculateDist()
        self.q=self.calculateQ()
        self.estK=self.ripleyFun()

    def calculateDist(self):
        distMat=auxFuncts.distMatGcirc(self.points)
        return(distMat)



    def calculateRho(self):
        """
         Estimates density from the 6th nearest neighbour mean
        """
        [self.nn6ind, self.nn6dist]=auxFuncts.nearestNeighbor(self.distMat,6)
        dmean=np.mean(self.nn6dist)
        self.rho=5/(np.pi*dmean**2)
        return

    def calculateEps(self):
        """
          Calculate epsilon
        """
        [self.nn1ind, self.nn1dist]=auxFuncts.nearestNeighbor(self.distMat,1)
        [std,bwR]=auxFuncts.calcSilvermanFactor(self.nn1dist)
        self.rads=np.linspace(min(self.nn1dist), max(self.nn1dist), 5000)
        densSamp=stats.gaussian_kde(self.nn1dist,bw_method=bwR/std)
        densTheo=auxFuncts.densNNDim(1,2,self.rho,self.rads)
        opcf=densSamp(self.rads)/densTheo
        difSigLog=np.diff(np.sign(np.log(opcf)))
        indice=np.where(difSigLog<0)[0]
        #print(indice)
        self.eps=0.5*(self.rads[indice[0]]+self.rads[indice[0]+1])
        return

    def calculateNmin(self):
        """
        Calculate minimum number of points
        """

        pval=0.0015
        toler=1.e-5
        k=0
        integ=1000
        while ((integ-pval)>toler) :
            k+=1
            integ=auxFuncts.integraDensNN(k,2, self.rho, 0, self.eps)


        self.Nmin=k+1
        return
    
    def setEps(self, eps):
        """
        sets the value of the distance attribute epsilon to detect structure
        """
        try:
            self.eps=float(eps)
        except ValueError:
            print("The distance epsilon must be a float (or convertible to)")

        return

    def setNmin(self, Nmin):
        """
        Sets the value of the attribute Nmin="minimum number of points" to detect structure
        """
        try:
            self.Nmin=int(Nmin)
        except ValueError:
            print("the minimum number of points must be an integer (or convertible to)")
        return
    
    def ripleyFun(self):
        vol = self.bBox.volEucl
        xmax = self.bBox.maxCoords[0]
        ymax = self.bBox.maxCoords[1]
        xmin = self.bBox.minCoords[0]
        ymin = self.bBox.minCoords[1]
        Kest = astrostats.RipleysKEstimator(area=vol , x_max=xmax, y_max=ymax, x_min=xmin, y_min=ymin)
        return (Kest)

    def calculateHomogRadius(self,n=100):
        mod='ohser'
        obs=self.estK(self.points,self.rads,mode=mod)
        [envHi, envLo]=self.envelopeRipley(n, mod)
        envTheo=self.estK.poisson(self.rads)
       # print(envTheo)
        high2=envHi+(envHi-envTheo)
        inEnvelope=(high2-obs <0)
       # print(inEnvelope)
        self.radHomog=0
        if(sum(inEnvelope)>0):
            self.radHomog=self.rads[inEnvelope][0]
        return


    def envelopeRipley(self,n,mod):
        env=np.zeros([n,len(self.rads)])
        for i in range(n):
           # print(i)
            rPoints=np.random.uniform(low=self.bBox.minCoords, high=self.bBox.maxCoords, size=self.points.shape)
            env[i,:]=self.estK(data=rPoints,radii=self.rads,mode=mod)

        hi=np.max(env,axis=0)
        lo=np.min(env,axis=0)

        return(hi,lo)
    
    def generateCSR(self):
        rand=np.random.uniform(low=self.bBox.minCoords, high=self.bBox.maxCoords, size=self.points.shape)
        return rand

    def increaseEpsilon(self, nSd=1):
        """
        Increases epsilon by nSd times the std deviation of the Nmin nearest neighbour distance distribution
        """
        exp=auxFuncts.integraRNDensNN(1,self.Nmin,2, self.rho, 0, np.inf)
        exp2=auxFuncts.integraRNDensNN(2,self.Nmin,2, self.rho, 0, np.inf)

        sd=np.sqrt(exp2-exp**2)
        print(sd)
        self.setEps(self.eps+nSd*sd)
        return

     
    def detectStructs(self):
        """
        Uses DBSCAN to detect significant substructure
        """
        self.db = DBSCAN(eps=self.eps, min_samples=self.Nmin, metric='precomputed').fit(self.distMat)

    