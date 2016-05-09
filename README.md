CG-Method
=========

CG-Method  obtains  mesoscopic  information  from  a  microscopic  configuration
(positions    and     velocities     in     a     lammps     log    file    (see
https://github.com/UNEDSoftMatter/Nanopore)  for  information  about  the lammps
script files.

All the parameters  are  specified  in  the  header  file `params.h`.  Note that
**you should compile whenever the  parameters are changed**.  A good practice is
to perform a `make clean`,  followed  by `make`.  An example `params.h` file can
be found in  `~/examples/`  directory.  **You  should  copy  the  file to ~/src/
before compiling the code**.

The output  files created by lammps  need to be formatted  in snapshots.  A good
practice is, for example:

```
~$ mkdir data/positions
~$ cp output.positions data/positions 
~$ cd data/positions
~/data/positions$ split -a 4 -d --lines=NAtoms+9 output.positions
~/data/positions$ for i in $(ls |grep x0); do cat $i | tail -n +10 | sort -n |awk '{print $2,$3,$4,$5}' > $i.pos ; done
~/data/positions$ rm x???? output.positions
~/data/positions$ cd ../../ 
~$ mkdir data/velocities
~$ cp output.velocities data/velocities 
~$ cd data/velocities
~/data/velocities$ split -a 4 -d --lines=NAtoms+9 output.velocities
~/data/velocities$ for i in $(ls |grep x0); do cat $i | tail -n +10 | sort -n |awk '{print $3,$4,$5}' > $i.vel ; done
~/data/velocities$ rm x???? output.velocities
~/data/velocities$ cd ../../ 
```
Once the snapshots  are created,  we may construct  a txt file with  the name of
the position files:

```
~$ for i in $(ls data/positions/); do basename $i .pos; done > files.txt
```

The program is the called as

```
~$ ./CG files.txt
```

Note that positions file has 4  columns:  the type of atom,  and the coordinates
x, y and z.  The velocities file has 3 columns: the velocities in x, y and z.

**If  your input  files have  other format  use  *awk*  to  adapt  to  the input
format**, or change the source code.


