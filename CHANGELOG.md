HISTORICAL MILESTONES
=====================

Rev#008 by @DiegoDZ
-------
- Changed parallel processing of velocities and positions in io.c. I had problems whith it. 
- Changed the number of decimals in PrintScalarWithIndex (in io.c) from 6 to 10. 


Rev#007 by @DiegoDZ
-------
- Added macrofunctions.c
- Added computation of MacroEnergy, MacroMomentum and module of the momentum. 
- Added some lines in order to compute the center of mass of the walls. 


Rev#007
-------
- Added computation of micro and meso momentum.
- Added computation of meso velocity.
- Fixed an  issue in the computation  of meso  velocity when  MesoDensity equals
  0.0.
- Added computation of internal energy.

Rev#006
-------
- The program now computes density profiles for type1 particles, type2 particles
  and the global density.
- New Parallelized regions:  `Compute_MesoSigma2`  (through force calculations),
  creation of positions and velocities snapshots, and mesoscopic profiles.
- New parallelized region: `Compute_Mean_Values`.

Rev#005
-------
- Cleaned up IO functions to make the code more readable.
- Deleted old commented lines.
- Now Sigma,  Sigma1, and Sigma2 are defined as `gsl_matrix(NNodes,9)`.  All the
  nine components are  computed at once and stored  in the following order:  00,
  01, 02, 10, 11, 12, 20, 21, 22.
- `Compute_Force_ij` now resets forces to zero each snapshot.
- The  program  now  uses  lammps  raw  trajectories  and  velocities  to create
  snapshots.  File locations **should be written in `src/params.h`**.
- The program  can be called with  or without arguments.  If you  do not provide
  an argument,  it will generate the  snapshot from lammps output.  If `argv[1]`
  is  provided,  then the  program will  use this  file list  without processing
  lammps output.
- Added information about blas, lapack and gsl.

Rev#004
-------
- Added a system call to create output directory. 
- Fixed  an  error  computing mean  values  (Sigma1_01  and  Sigma1_12  were not
  correctly closed).
- Fixed an  error computing Sigma2 (virial stress  tensor).  If a j-particle was
  over an i-particle, then zij gave a negative value.
- Code optimized. We changed gsl_matrix_column to gsl_vector_views.
- Fixed  some  memory  overflow.  Changed  GetLJParameters  from  `double  *` to
  `void`  and  added  deallocation  of   `int  *Verlet`  and  `double  *fij`  in
  Compute_Sigma2.

Rev#003
-------
- Parameters are  now in  a separated  file `params.h`.  The  file is  stored in
  `~/examples`  directory.   User  should  copy   the  file  to  `~/src`  before
  compilation.
- Fixed an error  on PBC in the  computation of the virial  stress tensor (index
  out of range for bin mu or nu = NNodes-1)
- Open issue?  There  is an error building  the linked  list on  some snapshots.
  lammps output gives negative values  for the positions,  which seems something
  weird.  The  program  now  checks  that  all  the  particles  are  inside  the
  simulation box, and it applies PBC for atoms outside the simulation box.
- There  is  a  new  second computation  step  in  which  averaged  profiles are
  obtained for mesoscopic quantities.
- Compute_Meso_Density  now computes  the  density  of  `type` particles.  Here,
  `type` is  a int value  which identifies a wall  (type 1 particle)  or a fluid
  (type 2 particle).
- Reset `gsl_vector * MesoDensity` each loop iteration.
- Print z coordinate in the microscopic files.
- Compute and print `\sigma_xy` and `\sigma_xz`.
- Added ERROR message if PBC happen.

Rev#002
-------
- The program  now uses `argv[1]`  as a filename  to read  all the  input files.
  The idea is  to create in bash `(for  i in $(ls |grep .pos))`  a txt file that
  will be read  by  the  program.  Note  that  the  txt  file  has  a limit of 6
  characters per line.
- Changed `gsl_matrix_get_col` to `gsl_vector_views`.
- Defined a struct to store information about all the output files.
- Now  the  program  does  not  create  a  file  for  each  snapshot.   All  the
  information  is  saved  in  the  same  file  if  it  corresponds  to  the same
  mesoscopic variable.
- All the vectors and matrices are reset to zero before computation.
- Definition of the vectors and matrices are now outside the main loop.

Rev#001
-------
- Changed strcpy in stdout to accomplish with `-O2` flag.
- Clean up of the source code.
- Splitted `functions.c` into `mesofunctions.c` and `microfunctions.c`.
- Added PBC for the computation of the virial stress tensor.
- The stress tensor is now computed only for the fluid.
