CG-Method
=========

*CG-Method*  obtains  mesoscopic   information  like  density  profiles,   force
interactions,  stress  tensors,  and so  on,  from  a  microscopic configuration
obtained through a MD simulation.  A MD  input script for LAMMPS can be obtained
in https://github.com/UNEDSoftMatter/Nanopore.

All the parameters  are  specified  in  the  header  file `params.h`.  Note that
**you should compile whenever the  parameters are changed**.  A good practice is
to do  a `make clean`,  followed by  `make`.  An example `params.h`  file can be
found in `~/examples/`  directory.  Note  that  **you  should  copy  the file to
~/src/ before compiling the code**. If you obtain the following error
```
~/CG-Method/src$ make
gcc  -c cg.c -o cg.c.o -lm -lblas -llapack -Wall -lgsl -lgslcblas -fopenmp -O2 -std=gnu99
In file included from cg.c:15:0:
cg.h:19:20: fatal error: params.h: No existe el fichero o el directorio
 #include "params.h"
                     ^
compilation terminated.
Makefile:13: recipe for target 'cg.c.o' failed
make: *** [cg.c.o] Error 1
```
it is because you did not copy `params.h` into the `src` directory.

Note  that  *CG-Method*  is currently  limited  to  100000  snaphots  (see `char
basename[7];` in `cg.c`).

The file  `params.h` should contains the  location (see below)  of the positions
file and the  velocities file created by LAMMPS.  In  LAMMPS,  use a custom dump
output like this one

```
dump custom positions 100 id type  x  y  z output.positions
dump custom positions 100 id type vx vy vz output.velocities

```

Please,  see LAMMPS documentation  to understand the meaning  of these commands.

*CG-Method* admits arguments.

1.  If you do not provide any argument,  *CG-Method* calls to the function `void
PrepareInputFiles(void)`   (located   in   `src/io.c`).   This   function   will
create  snapshots  for   the  positions  and  velocities   using  system  calls.
Then,  it will  create a file  list called `sim`  which will be  used to process
each snapshot.

Essentially, `PrepareInputFiles` will do the following:
```
~$ if [ ! -d data/positions ]; then mkdir -p data/positions; fi
~$ cd data/positions
~/data/positions$ ln -s ../../PositionsFileStr ./output.positions 
~/data/positions$ split -a 4 -d --lines=NParticles+9 output.positions
~/data/positions$ for i in $(ls |grep x); do cat $i | tail -n +10 | sort -n |awk '{print $2,$3,$4,$5}' > $i.pos ; rm $i ; done
~/data/positions$ rm output.positions
~/data/positions$ cd ../../ 
~$ if [ -d data/velocities ]; then mkdir -p data/velocities; fi
~$ cd data/velocities
~/data/velocities$ ln -s ../../VelocitiesFileStr ./output.velocities 
~/data/velocities$ split -a 4 -d --lines=NParticles+9 output.velocities
~/data/velocities$ for i in $(ls |grep x); do cat $i | tail -n +10 | sort -n |awk '{print $3,$4,$5}' > $i.vel ; rm $i ; done
~/data/velocities$ rm output.velocities
~/data/velocities$ cd ../../ 
```

Note  that  both  `PositionsFilesStr`  and  `VelocitiesFilesStr`  are  given  in
`src/params.h`.  **You should  specify the relative  path from where  the binary
`CG` is located**.

Once  the  snapshots  are created,  *CG-Method*  creates  a  txt  file  with the
basename of each snapshot:

```
~$ for i in $(ls data/positions/ | grep .pos); do basename $i .pos; done > sim 
```

The program uses the file `sim` to process each snapshot.

2.  If you  provide an argument to  *CG-Method*,  it will  omit the  creation of
snapshots and it will use `argv[1]` as the file list.

Note  that the  position files  have  4  columns:  the  type  of  atom,  and the
coordinates x,  y and z.  The velocitiy files have 3 columns:  the velocities in
x, y and z.

