# DVGardener
Matlab toolbox for analysis, manipulation and compression of Deformation Vector Fields (DVFs).
## Getting Started
We recommend using the "Example.m" script to see how DVF manipulations and analysis can be performed. This script is written in-line rather than using helper functions, but the helper functions are spun out from content of this script.

Example images are taken from The Cancer Imaging Archives (TCIA)
https://www.cancerimagingarchive.net/collection/4d-lung/

Registration is performed using Elastix
https://www.cancerimagingarchive.net/collection/4d-lung/

We recommend adding the "modules" to your matlab path (these primarily help with reading/writing medical images). Note these were written by Andy Shieh for the Image X Institute.

We also recommend adding "imshow3Dfull" to your modules for the sake of visualisation inside Matlab
https://au.mathworks.com/matlabcentral/fileexchange/47463-imshow3dfull

For viewing outside Matlab, install 3DSlicer
https://www.slicer.org/

## Citation
This work was included as a poster presentation at the 2024 AAPM annual meeting. Cite as

Dillon, O.T., Lau, B.K.F., Kertesz, H., Keall, P.J., 2024, July. Dvgardener: An Open-Source Toolbox for Manipulation, Quantification, Generation and Compression of Deformation Vector Fields. In AAPM 66th Annual Meeting & Exhibition. AAPM.
