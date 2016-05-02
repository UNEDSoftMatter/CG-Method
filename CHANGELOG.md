HISTORICAL MILESTONES
=====================


Rev#002
-------
- The program now uses argv[1] as a filename to read all the input files.
  The idea is to create through bash (for i in $(ls |grep .pos)) a txt
  file that will be read by the program.

Rev#001
-------
- Changed strcpy in stdout to accomplish with -O2 flag
- Clean up of the source code
- Splitted functions.c into mesofunctions.c and microfunctions.c
- Added PBC for the computation of the virial stress tensor
- The stress tensor is now computed only for the fluid
