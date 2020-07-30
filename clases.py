import numpy as np
from scipy import spatial, stats, integrate
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import stats as astrostats
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import auxFuncts
class StarFormingRegion:
    def __init__(self, puntos):
        """
        Check dimensions of matrix are appropriate, 
        initialization of attributes
        """
        if(not len(puntos.shape)==2):
            print('not valid dimension!!!!')
            print(puntos)
            return
        try:
            self.points=puntos.astype('float')
        except TypeError:
            print('not valid kind of points!!!')
            return
        self.distMat=spatial.distance_matrix(self.points, self.points)
        self.q=self.calculateQ()
        self.bBox=bBoxParams(self.points)
        self.signif=99.85


    def calculateQ(self):
        return auxFuncts.paramQ(self.distMat)


class bBoxParams:
    def __init__(self,points):
        self.minCoords=np.min(points,axis=0)
        self.maxCoords=np.max(points, axis=0)
        self.volEucl=np.prod(self.maxCoords-self.minCoords)
        

    
class SFRegion2D(StarFormingRegion):
    def __init__(self, puntos, coords):
        super(SFRegion2D,self).__init__(puntos)
        self.coords=coords
        print(coords)
        if (coords == 'Ra Dec'):
            print('Great circle distance')
            self.distMat=auxFuncts.distMatGcirc(self.points)
        else:
            print('Euclidean distance')
            self.distMat=spatial.distance_matrix(self.points, self.points)
        self.q=self.calculateQ()
        self.estK=self.ripleyFun()




    def calculateRho(self):
        """
         Estimates density from the 6th nearest neighbour mean
        """
        [self.nn6ind, self.nn6dist]=auxFuncts.nearestNeighbor(self.distMat,6)
        dmean=np.mean(self.nn6dist)
        self.rho=1.151082*5/(np.pi*dmean**2)
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
        #print('self.signif', self.signif)
        pval=1-self.signif/100.0#0.0015
        #print('pval',pval)
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

    
    def setSignif(self, signif):
        """
        sets the value of the distance attribute epsilon to detect structure
        """
        try:
            self.signif=float(signif)
        except ValueError:
            print("The significant level must be a valid percentage number (between 0 and 100)")

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

     
    def detectStructs(self):
        """
        Uses DBSCAN to detect significant substructure
        """
        self.db = DBSCAN(eps=self.eps, min_samples=self.Nmin, metric='precomputed').fit(self.distMat)

    def plotStructs(self, ofile):
        """
        plots a map with the structures found by dbscan which is in self.db
        """
        labels=self.db.labels_
        fig = plt.figure()
        ax = fig.add_subplot(111)
        nstruct=len(set(labels)) - (1 if -1 in labels else 0)
        cm=plt.get_cmap('viridis')
        cs=cm(np.arange(nstruct)/nstruct)
        plt.scatter(self.points[:,0],self.points[:,1], marker='o', c='gray', s=5, alpha=0.7)
        for x in range(nstruct):
            inds=(labels==x)
            plt.scatter(self.points[inds,0],self.points[inds,1], marker='o', c=cs[x].reshape(1,4), s=10, alpha=0.7)
   
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        if(self.coords=='Ra Dec'):
            plt.xlim(max(self.points[:,0]),min(self.points[:,0]))
            plt.ylabel('DEC (deg)', fontsize=14)
            plt.xlabel('RA (deg)', fontsize=14)
        else:
            plt.ylabel('X', fontsize=14)
            plt.xlabel('Y', fontsize=14)           

        plt.savefig(ofile)
        plt.close(fig)

    



    
class SFRegion3D(StarFormingRegion):
    def __init__(self, puntos, coords):
        super(SFRegion3D,self).__init__(puntos)
        self.coords=coords
        #if (coords == 'Ra Dec'):
            #print('Converting to XYZ automatically with astropy skycoord!!!!')
            #print('assume ra dec in deg and dist in pc')
            #c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, distance=770*u.kpc)
        #else:
        #    print('Euclidean distance')
            
        self.distMat=spatial.distance_matrix(self.points, self.points)
        self.q=self.calculateQ()

    def calculateDist(self):
        distMat=auxFuncts.distMatGcirc(self.points)
        return(distMat)



    def calculateRho(self):
        """
         Estimates density from the 6th nearest neighbour mean
        """
        [self.nn6ind, self.nn6dist]=auxFuncts.nearestNeighbor(self.distMat,6)
        dmean=np.mean(self.nn6dist)
        self.rho=3*5.673173/(4*np.pi*dmean**3)##3D version!!! The expression is different
        ##1.151082*5/(np.pi*dmean**2)
        return

    def calculateEps(self):
        """
          Calculate epsilon
        """
        [self.nn1ind, self.nn1dist]=auxFuncts.nearestNeighbor(self.distMat,1)
        [std,bwR]=auxFuncts.calcSilvermanFactor(self.nn1dist)
        self.rads=np.linspace(min(self.nn1dist), max(self.nn1dist), 5000)
        densSamp=stats.gaussian_kde(self.nn1dist,bw_method=bwR/std)
        densTheo=auxFuncts.densNNDim(1,3,self.rho,self.rads)
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
        #print('self.signif', self.signif)
        pval=1-self.signif/100.0#0.0015
        #print('pval',pval)
        toler=1.e-5
        k=0
        integ=1000
        while ((integ-pval)>toler) :
            k+=1
            integ=auxFuncts.integraDensNN(k,3, self.rho, 0, self.eps)


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

    
    def setSignif(self, signif):
        """
        sets the value of the distance attribute epsilon to detect structure
        """
        try:
            self.signif=float(signif)
        except ValueError:
            print("The significant level must be a valid percentage number (between 0 and 100)")

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

     
    def detectStructs(self):
        """
        Uses DBSCAN to detect significant substructure
        """
        self.db = DBSCAN(eps=self.eps, min_samples=self.Nmin, metric='precomputed').fit(self.distMat)

    def plotStructs(self):
        """
        plots a map with the structures found by dbscan which is in self.db
        """
        labels=self.db.labels_
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        nstruct=len(set(labels)) - (1 if -1 in labels else 0)
        cm=plt.get_cmap('viridis')
        cs=cm(np.arange(nstruct)/nstruct)
        inds=(labels==-1)
        ax.scatter3D(self.points[inds,0],self.points[inds,1],self.points[inds,2], marker='o', c='gray', s=3, alpha=0.5)
        for x in range(nstruct):
            inds=(labels==x)
            ax.scatter3D(self.points[inds,0],self.points[inds,1],self.points[inds,2], marker='o', c=cs[x].reshape(1,4), s=10)

        plt.show()

    
