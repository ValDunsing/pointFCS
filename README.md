# pointFCS
A Matlab package to analyze point fluorescence correlation spectroscopy data, associated to "Spatial patterning of contractility by a mechanogen gradient underlies Drosophila gastrulation" by Mundhe, Dunsing-Eichenauer, Philippe et al., bioRxiv 2025, https://doi.org/10.1101/2025.04.11.648359.

The main code can be found in pFCS_analysis_czi_batch.m, which loads multiple exported tif files (each corresponding to one FCS measurement) from a single directory and executes an FCS batch analysis of these files. 
The code to measure the distance of the FCS spot from the invagination (capture in a fulle FOV image) can be found in pFCS_image_distances.m.
To pool background intensities from yw embryos, use pool__bg_intensities.m
