HISTORICAL MILESTONES
=====================

Rev#002
-------
- The program now  uses argv[1] as a filename to read  all the input files.  The
  idea is to create in bash (for i in  $(ls |grep .pos)) a txt file that will be
  read by the program.  Note  that the txt file has a limit  of 6 characters per
  line.
- Changed gsl_matrix_get_col for gsl_vector_views.
- Defined a struct to store information about all the output files.
- Now the program does not create a file for each snapshot. All the information
  is saved in the same file if it corresponds to the same mesoscopic variable.
- All the vectors and matrices are reset to zero before computation.
- Definition of the vectors and matrices are now outside the main loop.

Rev#001
-------
- Changed strcpy in stdout to accomplish with -O2 flag
- Clean up of the source code
- Splitted functions.c into mesofunctions.c and microfunctions.c
- Added PBC for the computation of the virial stress tensor
- The stress tensor is now computed only for the fluid
