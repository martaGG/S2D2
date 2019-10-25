# R2D2: Small-scale Significant DBSCAN Detection

Python3 code for the detection of small pristine structure with 3-sigma significance associated to position/velocity of stars (or young stellar objects) in starForming Regions. Developed under the SFM (StarFormMapper) EU project.

S2D2 uses DBSCAN detect the smallest significant structure in a spatial/spatial-kinematic space. The procedure proposes, in structured regions, calculation of the epsilon and Nmin parameters (as described in REFS!!) for DBSCAN to retrieve the smallest structures in the region with a 3-sigma level of significance. If the region is not structured, or the user wants it, eps and Nmin can be supplied, and Dbscan will be performed with these parameters, without, however, guaranteeing any level of significance.

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

#### coordinates
Coordinate frame of the input, dependi
V0: 'Ra, Dec'
Future: 'l,b'

