CG-Method
=========

CG-Method obtains mesoscopic information from a microscopic configuration.

All the  parameters are specified in  the header  file *cg.h*.  Note  that **you
should compile  whenever the  parameters are changed**.  A  good practice  is to
perform a *make clean*, followed by *make*.

The microscopic information enter as two  matrices stored in two different files
(see *cg.h*):

    define iFilePosStr  "data/Positions.sort.3col.dat"
    define iFileVelStr "data/Velocities.sort.3col.dat"

The positions file has 5 columns:  the label of an atom,  the type of atom,  and
the coordinates x,  y and z.  The velocities file has 4 columns: the label of an
atom and the velocities in x, y and z.

**If  your input  files have  other format  use  *awk*  to  adapt  to  the input
format**, or change the source code.


