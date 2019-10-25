# S2D2: Small-scale Significant DBSCAN Detection

Python3 code for the detection of small pristine structure with 3-sigma significance associated to position/velocity of stars (or young stellar objects) in starForming Regions. Developed under the SFM (StarFormMapper) EU project.

S2D2 uses DBSCAN detect the smallest significant structure in a spatial/spatial-kinematic space. The procedure proposes, in structured regions, calculation of the epsilon and Nmin parameters (as described in our paper REFS!!) for DBSCAN to retrieve the smallest structures in the region with a 3-sigma level of significance. If the region is not structured, or the user wants it, eps and Nmin can be supplied, and Dbscan will be performed with these parameters, without, however, guaranteeing any level of significance.

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
  filename: './data/Frac1p6_1_RADEC.dat',
  dim: 2,
  coord: 'Ra Dec',
  eps: 'None',
  Nmin: 'None'
}
```

### Input parameters

#### filename
Path for file with coordinates of stars/objects in your region.
V0: Ascii file with header and tab delimiter

#### dim
Dimension of the space of search.
V0: Integer, only=2.
Future: limited options ('2D','3D','2+2D', '3+2D'...)

#### coord
Coordinate frame of the input, depending on the dimension
##### 2D
V0: 'Ra, Dec'
Future: 'l,b'

#### eps
Scale parameter supplied to DBSCAN, associated with the size of the structures to search.
Default: 'None' to search for the smallest scale in structured regions.
If a float eps is supplied, Nmin must also be supplied.


#### Nmin
Number of points supplied to DBSCAN, associated with the density of a neighbourhood with radius eps.
Default: 'None' to calculate the Nmin guaranteeing 3 sigma significance.
If an integer Nmin is supplied, eps must also be supplied.

## Output
Ascii file named as the input file with the extension .out containing the coordinates of each star in the region (in RA, DEC) and an additional column with an integer representing the number of substructure assigned. The value -1 represents noise stars (those not assigned to any cluster).

## Requirements
Python3 with libraries:
	astropy
	numpy
	scipy
	scikit.learn
	yaml
    argparse
## Description

### 2D
We refer to (SFM paper, REF!!!) and references therein for a complete description of the procedure.

Our default coordinates are Ra, Dec, and the distance used is the great circle one. 

### Structured regions
We will consider that a starForming region is structured when the Q parameter (Cart&With…REF!!!) is lower than 0.7. In that case, and if the user has not provided eps and Nmin (default) we propose our procedure to calculate them and obtain the smallest scale significant structure in the region.

#### eps calculation: Small-scale
We calculate the length scale for DBSCAN (epsilon) using the One point correlation function, OPCF, comparing the (REF!!! Joncour et al. paper I) first nearest neighbour distance distribution of the sample with the first nearest neighbour distance distribution of a homogeneous random distribution (CSR, or complete spatial randomness) with intensity equal to the local density derived from the mean of the 6th neighbour distribution of the sample.

#### Nmin calculation: Significance
We iteratively calculate the significance of a structure of that scale and a fixed number of points k until we reach 3-sigma significance (~99.85%). The significance of a structure of size eps and k points, as described in REF!! Joncour et al Paper 2, is given by the the probability of having k-1 nearest neighbours in an eps neighbourhood under a homogeneous random distribution with intensity rho.

#### DBSCAN detection

We run scikit.learn’s DBSCAN for the previously calculated eps and Nmin. 


### Unstructured regions

If the Q parameter is larger than 0.7, and the user has not supplied an eps and Nmin, we display an error message (we cannot guarantee that the region is structured), and suggest the user to try the procedure providing eps and Nmin.

#### DBSCAN detection
We run scikit.learn’s DBSCAN for the user defined eps and Nmin.

## Acknowledging this
Please cite (SFM, our paper, REF!!) if you use this code. 

## References

To be completed



