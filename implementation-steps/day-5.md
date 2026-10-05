# Day 5

Copy the saved TFM tutorial .mfx into cases/tut_tfm2d/. Write src/post/read_monitors.py to load MFiX monitor CSV output into pandas, and src/post/read_vtk.py to load VTK output with PyVista or meshio. Write a script that computes (a) the time series of pressure drop across the bed and (b) the time-averaged solids fraction field over the last half of the run, then plots both. Then write a short notes/mfx_anatomy.md that annotates the tutorial .mfx file line by line, explaining what each keyword does as documented in the keyword reference (rule 1). 

Done when: both plots are reproduced from the tutorial output, and the annotated .mfx has no undocumented keywords.