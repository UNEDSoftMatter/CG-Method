CG-Method
=========

*CG-Method*  obtains  mesoscopic   information  like  density  profiles,   force
interactions,  stress  tensors,  and so  on,  from  a  microscopic configuration
obtained through a MD simulation.  A MD  input script for LAMMPS can be obtained
in https://github.com/UNEDSoftMatter/Nanopore.

All the parameters  are  specified  in  the  header  file `params.h`.  Note that
**you should compile whenever the  parameters are changed**.  A good practice is
to perform a `make clean`,  followed  by `make`.  An example `params.h` file can
be found  in `~/examples/` directory.  Note that  **you should copy  the file to
~/src/ before compiling the code**.

Note  that the  program is  currently  limited  to  100000  snaphots  (see `char
basename[7];` in `cg.c`).

*CG-Method*  needs the  location (see  below)  of  the  positions  file  and the
velocities file  created by LAMMPS.  In  LAMMPS,  use a custom  dump output like
this one

```
dump custom positions 100 id type  x  y  z output.positions
dump custom positions 100 id type vx vy vz output.velocities

```

Please,  see LAMMPS documentation  to understand the meaning  of these commands.
CG-Method  uses   the   function   `void   PrepareInputFiles(void)`  (locate  in
`src/io.c`)  to create  snapshots from  the original  LAMMPS dumps  using system
calls. The following instructions are performed by `PrepareInputFiles`

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

The program uses the file `sim` to process each snapshot

Note that position files have 4 columns:  the type of atom,  and the coordinates
x, y and z.  The velocitiy files have 3 columns: the velocities in x, y and z.

