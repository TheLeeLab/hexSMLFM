# Single Molecule Light Field Microscopy Reconstruction (hexSMLFM)
**Authors:** Ruth R. Sims, Sohaib Abdul Rehman, Ezra Bruggeman, Sam Daly, and Kevin O'Holleran
_Cambridge Advanced Imaging Centre, Downing Site, Cambridge, CB2 3DY, UK_
_Yusuf Hamied Department of Chemistry, Lensfield Road, University of Cambridge, Cambridge, CB2 1EW, UK_

This 3D reconstruction code accompanies the paper of S. Daly _et al._ "[High-density volumetric super-resolution microscopy](https://www.nature.com/articles/s41467-024-45828-5)". In brief, 2D localisation data (x,y) captured using a (square or hexagonal) fourier light field microscope are turned into 3D localisations (x,y,z). It does this by first assigning (x,y) to (x,y,u,v) space and using microscope parameters to calculate z position via parallax. For more information, see: [R. R. Sims, *et al.* Optica, 2020, 7, 1065](https://doi.org/10.1364/OPTICA.397172).

The code hosted here on GitHub is actively maintained by the authors of the preprint and is freely available for download and use. Please get in touch if you encounter any difficulties. 

Example data has been included (example_data_bcr_fixed_bcell.csv), which comprises 5,000 frames of a whole B cell membrane captured using single molecule light field microscopy.

<b>Note for Python users:</b> You may wish to try our recently released Python version of hexSMLFM produced in collaboration with "[Teledyne Photometrics](https://www.photometrics.com/)". This version of the code can be found on GitHub "[here](https://github.com/Photometrics/PySMLFM)".
