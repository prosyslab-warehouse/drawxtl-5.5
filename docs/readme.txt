 Changes in V5.5 Released March 10, 2011

New Features:
-------------

1. New POV output option for stereo images based on new features in POV-Ray version 3.7
    (which no longer require a specially modified version of POV-Ray)

 2. New command "x3d" to create VRML files in Web3D-compliant "X3D Classic VRML" encoding 
    instead of plain VRML97.  

 3. The "clear omit flag" option in the Parameters menu is now unset after successful
    application

 4. Added direct export of the structure image to a GIF file

 5. Modified the WIEN2k import filter to support file format changes introduced by WIEN2k_09

 6. MSMS-generated surfaces can now be colored according to atom type

 7. The number of dash/space segments in dashed bonds (default 5) can now be specified

 8. CIF import now supports space group declaration by IT number (space_group_IT_number, 
    _symmetry_int_tables_number)

 9. Reading input files from an unwritable directory (or medium) is now supported
    (in the sense that the file is copied to a readable location before continuing)

10. New command "mapslice" for generation of 2D sections from fourier maps
    at arbitrary angles, displaying either contours or color-coded images. A
    GUI has been added to set/alter the parameters.

11. New output option for Asymptote input files (offering 3D PDF)

12. New command "qvector" for specifying modulation wave vectors directly
    in the str file instead of importing them from CIF data.

13. New output window for measured distances and angles

14. Added reading of map files in the .xsf format used by XCrysDen (a format
    written by many "quantum chemistry" codes). 

15. Adjusted some windows so that they fit on the 1024 x 576 screen of a Netbook.

16. New dialog window for Surface-related options

17. New screen to generate "movies" from frames computed with differing phase angles
    for the modulation parameters.

Bugs Fixed:
-----------

 1. Correct conversion of stf map sections that contain negative coordinate 
 2. Correct recognition of space group Fmmm in shelx files (was misinterpreted
    as cubic due to the LATT 4 instruction)
 3. Correct handling of the "transparency" attribute in the polyhedra menu
 4. Correct a Windows-specific bug that could lead to overwriting of map files
    when the structure file is saved
 5. Correct potential loss of color information during the conversion from spheres to ellipsoids
 6. OpenGL "picking" for object deletion did not work in all cases
 7. Fixed memory management issues in the cavities code (voids mode 1)
 8. Ellipsoid orientation was incorrect for fourfold symmetric sites in
    tetragonal spacegroups
 9. Fixed an error in the ideal tetrahedral angle used in calculating angle variance.
10. Fixed a problem with the eigenvalues for some ellipsoids of revolution
    on Windows. The cube root routine failed when the input was negative.

Uncorrected BUGS:
-----------------


Possible additions (not yet implemented)
----------------------------------------

 * add keyboard command to add atoms to display list that are a
   specified radius around the graphics cursor

 * add command to clear atoms cursor list

 * "add frame" widget (not sure what it should contain - probably
   input widgets for 'offset' and a checkbox for 'copy cell and atoms')

 * label rotation

 * DIAMOND-style legend box (one labeled sphere for each type of atom)

 * "add ellipsoids..." options for other import types

 * create mng or swf animations

 * convert rhombohedral setting to hexagonal

 * convert between conventional and primitive (Niggli) cell

 * add keyboard command to position graphics cursor at given point

 * add command for viewing along a crystallographic direction

 * report current camera and light source position in the POV options menu
   and let the user modify them (to solve lighting problems without resorting
   to direct editing of the POV file, to fix camera position for a series of images
   to be converted to an animation)

 * replace the .frm files with internal arrays and combine .cns and .out to
   reduce "directory pollution"

 * provide more convenient (interactive?) option for choosing "mapslice" plane
   orientation, perhaps use x/y/z rotation angles instead of normal vector on
   mapslice card

