# S2D2: Small Significant DBSCAN Detection

Python3 code for the detection of small pristine structure with user defined associated to position/velocity of stars (or young stellar objects) in starForming Regions. Developed under the SFM (StarFormMapper) EU project.

S2D2 uses DBSCAN detect the smallest significant structure in a spatial/spatial-kinematic space. The procedure proposes, in structured regions, calculation of the epsilon and Nmin parameters (as described in González et al. 2020 and references therein) for DBSCAN to retrieve the smallest structures in the region with a minimum level of significance. If the region is not structured, or the user wants it, eps and Nmin can be supplied, and Dbscan will be performed with these parameters, without, however, guaranteeing any level of significance.

Usage example:
```
python3 main.py --fname data/inputExample.yaml 
```

## Input

yaml file with a dictionary containing the input parameters. Must be provided after the call to main.py with --fname

Example call:
```
python3 main.py --fname data/inputExample.yaml 
```

Example file content:
```
{
  filename: './data/Frac1p6_1_RADEC.dat',
  dim: 2,
  coord: 'Ra Dec',
  eps: 'None',
  Nmin: 'None',
  Qlim: 0.7,
  Signif: 99.85
}
```

### Input parameters

##### 1. filename
Path for file with coordinates of stars/objects in your region.
V1: Ascii file with header and tab delimiter

##### 2. dim
Dimension of the space of search.
V1: Integer, only=2.
Future: limited options ('2D','3D','2+2D', '3+2D'...)

##### 3. coord
Coordinate frame of the input, depending on the dimension
- 2D
  - 'Ra Dec': expects input data in Right ascension, declination, in degrees to calculate the great circle distance.
  - 'X Y': expects float arbitrary input data coordinates, calculates euclidean distance.

  **if a different string is input, by default it will go into 'X Y' mode 

##### 4. eps
Scale parameter supplied to DBSCAN, associated with the size of the structures to search.
-Default: 'None' to search for the smallest scale in structured regions.
-If a float eps is supplied, Nmin must also be supplied.

##### 5. Nmin
Number of points supplied to DBSCAN, associated with the density of a neighbourhood with radius eps.
-Default: 'None' to calculate the Nmin guaranteeing the significance supplied by the user above random expectation.
-If an integer Nmin is supplied, eps must also be supplied.

##### 6. Qlim
limit of Q parameter for considering a region structured, and calculate automatically the eps and Nmin values according to the procedure, as described in González et al. 2020.

We note that the classical limit Q parameter for structured regions is 0.8, and that a conservative limit of 0.7 avoids the possibility of analyising a region withous structure (as described in Gonzalez et al. 2020)

##### 7. Signif
Significance limit above which structures will be retrieved. Ussed to calculate the minimum number of points.  
Must be an appropriate percentage value

## Output
- Console outputs some values of variables and status 

- Ascii file named as the input file with the extension .out containing the coordinates of each star in the region as in the input, and an additional column with an integer representing the number of substructure assigned. The program uses the default python convention, so value -1 represents noise stars (those not assigned to any cluster).

- pdf file named as the input file with the extension .pdf with a plot of the region where:
  - grey stars are noise.
  - stars in significant structures are overplotted in colour. Each nest will be plotted in a different colour taken from a viridis colour table with as many different shades as NESTs, so the specific colours will depend on the amount of structures retrieved. 
## Requirements
Python3 with libraries:
-	astropy
-	numpy
-	scipy
-	scikit.learn
-	yaml  (can be found in pip as pyyaml for installation)
- argparse
- matplotlib
## Description

### 2D
We refer to [González et al. 2021](https://www.aanda.org/component/article?access=doi&doi=10.1051/0004-6361/202038123) and references therein for a complete description of the procedure.


### Structured regions
We will consider that a starForming region is structured when the Q parameter (Cartwright & Witworth, 2003) is lower than the user supplied limit Q lim. In that case, and if the user has not provided eps and Nmin (default) we propose our procedure to calculate them and obtain the smallest scale significant structure in the region.

#### eps calculation: Small-scale
We calculate the length scale for DBSCAN (epsilon) using the One point correlation function, OPCF, comparing the (REF!!! Joncour et al. paper I) first nearest neighbour distance distribution of the sample with the first nearest neighbour distance distribution of a homogeneous random distribution (CSR, or complete spatial randomness) with intensity equal to the local density derived from the mean of the 6th neighbour distribution of the sample.

#### Nmin calculation: Significance
We iteratively calculate the significance of a structure of that scale and a fixed number of points k until we reach 3-sigma significance (~99.85%). The significance of a structure of size eps and k points, as described in Joncour et al 2018, is given by the the probability of having k-1 nearest neighbours in an eps neighbourhood under a homogeneous random distribution with intensity rho.

#### DBSCAN detection

We run scikit.learn’s DBSCAN for the previously calculated eps and Nmin. 


### Unstructured regions

If the Q parameter is larger than the user supplied region, and the user has not supplied eps and Nmin values, we display an error message (we cannot guarantee that the region is structured), and suggest the user to try the procedure providing eps and Nmin.

#### DBSCAN detection
We run scikit.learn’s DBSCAN for the user defined eps and Nmin.

## Acknowledging this
Please cite [González et al. 2021](https://www.aanda.org/component/article?access=doi&doi=10.1051/0004-6361/202038123) if you use this code. 
## References
To be completed
- Cartwright & Withworth, 2003
- [González et al. 2021](https://www.aanda.org/component/article?access=doi&doi=10.1051/0004-6361/202038123)
- Joncour et al. 2017.
- Joncour et al. 2018.



