# Cell_medial_axis_definitions
This repository includes custom functions that can be used to draw the medial cell axis

Bivariate_medial_axis_estimation:
  The function get_medial_axis draws the central line of elongated (rod-shaped) cells from pole to pole. A set to test the medial axis estimation is also provided along with a notebook that demonstrates how the function works. The get_oned_coordinates function projects all the cell pixels onto the central line,
  in relative or absolute arch-length coordinates from pole to pole. The distance from the central line is multiplied by the sign of the cross product to get the 
  orientation of the cell pixels around the medial axis (above or below). An example for implementing the medial axis estimation is provided in the test_medial_axis.ipynb notebook.

  Cite:
    https://www.biorxiv.org/content/10.1101/2024.10.08.617237v2.full
    
    DNA/polysome phase separation and cell width confinement couple nucleoid segregation 
    to cell growth in Escherichia coli
    
    Alexandros Papagiannakis, Qiwei Yu, Sander K. Govers, Wei-Hsiang Lin,  Ned S. Wingreen, Christine Jacobs-Wagner
    
    bioRxiv, https://doi.org/10.1101/2024.10.08.617237, October 22, 2024
