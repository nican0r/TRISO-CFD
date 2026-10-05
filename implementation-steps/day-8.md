# Day 8

Build make_case.py templates for: (a) V_2d: a 2D axisymmetric (cylindrical coordinates; this is chosen because the Shallbetter thesis demonstrated that cylindrical coordinates agree better with experimental beds that were built, make a note of this in the docs) spouted bed with the cone and inlet orifice, and (b) V_3d: a 3D bed built from an STL generated in Python (cylinder + 60° cone + orifice), using MFiX's cut-cell/STL workflow. Restrict the domain height to about 3× the static bed height plus fountain headroom rather than the full 1.14 m column, and document this choice. Mesh targets: the 2D inlet orifice spans ≥ 4 cells; for 3D, use a cell size of ~2–3 mm to start, estimate the cell count, and report it. 

Initial condition: a packed bed at ε_s ≈ 0.6 up to H, with gas at rest. Inlet: mass inflow through the orifice. 
Outlet: pressure outflow at the top. Add monitors for the pressure at the inlet plane and the top of the bed, and for the solids mass flux through a horizontal plane in the spout at mid-bed height. 

Done when: both cases initialize without errors and the geometry is visually confirmed in the MFiX GUI or ParaView (screenshots in results/).