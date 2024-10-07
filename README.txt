I have 2 fortran files: gol.f90 and gol2.f90

The gol.f90 file is the file which I use for the simple pattern and question 3 of report.
I also use timers for this file to answer the last part too.

On the other hand, I use gol2.f90 solely just for column-wise domain decomposition. I was unsure
of how to merge the parts in gol.f90 and this file together. It works very well with multiple processes, 
the grid isn't initialized in the same way.

Compile gol.f90: mpif90 gol.f90 -o gol
Run gol.f90: mpirun -np <Number of Processors> ./gol

Compile gol2.f90: mpif90 gol2.f90 -o gol2
Run gol.f90: mpirun -np <Number of Processors> ./gol2

To adjust rows and cols, adjust nrows and ncols
To adjust steps, modify N

If you want to initialize grid randomly, you can uncomment the section of my code which does that
