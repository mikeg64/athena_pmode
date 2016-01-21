# athena_pmode
#Version of athena developed to model pmode oscillations

### compiling building and testing
General guidelines to do this are in the quick start tutorial at

https://trac.princeton.edu/Athena/wiki/AthenaDocsTutOT1

To build this model

- % make clean
- % configure --with-problem=solp --with-order=3
- % make all

To run the model

- % cd bin
- % athena -i ../tst/2D-mhd/athinput.solp
