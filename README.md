# -dynamic_spatial_panel
http://myweb.ttu.edu/dolacomb/matlab.html

MATLAB Code
Dynamic Spatial Panel Data Model Code
The following code is an adaptation of Paul Ehhorst's dynamic spatial panel data code with two additional features. First, the stability condition is calculated along with its associated 95% confidence interval. Second, the effects estimates are printed out in an easy to read format suitable for cutting and pasting into papers, presentations, etc.

The file "handbook82_Lacombe_original.m" is the same file as in the original Elhorst code but with the two additions mentioned above.

Dynamic Spatial Panel Data Code with Improvements: Download Here

KNN Weight Matrix Code
The following code produces a k nearest neighbors spatial weight matrix using the Great Circle formula.

The routine generates a sparse spatial weight matrix using user supplied coordinates.

The results can be more accurate than using Euclidean distance formulas. The routine can be slower than other
routines if the data set is large.

A demonstration file is included in the zip archive.

KNN Weight Matrix Code Great Circle Distance: Download Here

Lagrange Multiplier Testing Suite
The Lagrange Multiplier Testing Suite contains the following tests for spatial dependence:

LM Lag Test
LM Error Test
LM Lag Robust Test
LM Error Robust Test
LM Combined Lag / Error Test
LM Spatial Error Components Test
LeSage and Pace Spatial Hausman Test

The Lagrange Multiplier Testing Suite is designed to be used in conjunction with Jim LeSage's Spatial Econometric Toolbox for MATLAB.

LM Testing Suite: Download Here

Lagrange Multiplier Testing Suite for Panel Data
The Lagrange Multiplier Testing Suite contains the following tests for spatial dependence:

LM Lag Panel Test
LM Error Panel Test
LM Lag Robust Panel Test
LM Error Robust Panel Test

A demonstration file is included in the zip folder. These files are designed to be used with Jim LeSage's Spatial Econometrics Toolbox for MATLAB.

LM Panel Testing Code: Download Here

Elhorst Spatial Panel Code with the LeSage and Pace Effects Estimates
The Elhorst Spatial Panel MATLAB code that has been extended to include the bias correction procedure of Lee and Yu (2010) now include the LeSage and Pace effects estimates that calculate the correct marginal effects in the case of a spatially lagged dependent variable.

These routines print out the effects estimates in a similar fashion to the LeSage and Pace spatial econometrics toolbox routines.

A demonstration file ("panel_effects_demo") is included in the zip file and these routines are designed to be used in conjunction with the LeSage and Pace Spatial Econometrics toolbox.

Elhorst Bias Corrected Spatial Panel Code with Effects Estimates

Distance Based Weight Matrix Code
The following zip file has code to calculate two different types of distance based weight matrices:

Distance Based on Cutoff Value
Inverse Distance Based on Cutoff Value
The code is currently designed to build weight matrices to ensure that all geographic entities have at least one neighbor. The code can be modified to allow for distances in kilometers or miles, or for an arbitrary distance cutoff, such as the average of all distances, or for a specific value.

Distance Based Weight Matrices: Download Here

Distance Based Weight Matrix with Bands
The following zip file contains code to calculate a distance based weight matrix with bands, where geographic entities that are withtin a certain distance band are considered neighbors.

Included in the zip file is the code and a demonstration file.

Distance Based Weight Matrix with Bands: Download Here
