# NuWTun compressible flow solver

## Contents

* `src-flo `: Flow solver files
* `src-adj `: Adjoint solver files, NOT FINISHED YET
* `examples`: Example test cases
* `docs    `: Manual, license, etc.

## Compiling

You will need a fortran 90 compiler since dynamic memory allocation is used.  Compilation has been tested with gfortran. Set your compiler and flags in `makefile.in` and define a shell variable called `NUWTUN_HOME` which contains the absolute path to the root directory. You also need the archive tool "ar"; this should be available on most systems. To compile the flow solver execute the following commands at the shell prompt:

```shell
cd $NUWTUN_HOME/src-flo
make
```

## Call-tree

If you have ftnchek, then you can generate a call-tree inside the `docs/tree/src-flo` directory. Inside the directory `src-flo` type `make tree` to generate the tree. This has been tested to work only on GNU/Linux. You can change the location where the call-tree is stored by changing the `CTDIR` variable in `makefile.in` to appropriate value.

## Visualization

Nuwtun can output the solution in plot3d format. Plot3d files can be read by the following visualization softwares:

* Tecplot: Commercial software, but probably the best one.
* Paraview: Free software, based on VTK, can read only grid and solution files, not the function file.
* Mayavi: Free software, similar to Paraview and is based on VTK. Only 3-D files can be read.
* Calculix: Free software, can read both 2-d and 3-d files, only formatted, can read grid, solution and function files. It is quite good for 2-D.

## Running flow solver

Compile

```shell
cd $NUWTUN_HOME
make flo
```

Run

```shell
$NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log 2>&1 &
```

Monitor the run with tail

```shell
tail -f flo.log
```

Some output files:

* `fort.17`: residual history
* `fort.18`: cl,cd history
* `fort.19`: some global results, useful for optimization

Plot convergence history (requires gnuplot)

```shell
gnuplot $NUWTUN_HOME/src-utl/nwthist.gnu
```

This creates files clcd.eps and res.eps

## Postprocessing

The output is in plot3d format which can be viewed using Calculix cgx or Paraview.

There are some programs in src-utl to process the output.

* `extract`: This reads flo.log to extract some surface data.
* `ext2cp`: Run this after running extract
* `plot2vtk_2d`: Writes a vtk file for 2d (use Visit/Paraview to see vtk file)
* `plot2vtk_3d`: Writes a vtk file for 3d (use Visit/Paraview to see vtk file)

Run make inside this directory to compile all the programs.


## Grid deformation tools

`src-grd` has radial-basis function based deformation tools for shape and grid.  See the `README` file inside this directory.

---

* `Origin`: https://codeberg.org/cpraveen/nuwtun
* `Mirror`: https://git.sr.ht/~cpraveen/nuwtun
* `Mirror`: https://github.com/cpraveen/nuwtun
