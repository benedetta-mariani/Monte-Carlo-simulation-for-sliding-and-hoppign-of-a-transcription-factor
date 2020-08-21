# Monte Carlo SImulation for sliding and hopping

The code here provided aims at simulating the search of a Binding factor of its specific binding site, which happens in through a 3D diffusion (Hopping) and 
a 1D diffusion along the DNA chain (sliding). The model is reproduced following https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.022721.
All the datails of the code, the organization of the classes and how to run the simulation are provided in the file Relazione.pdf.

To run the code it is necessary the ROOT framework. Run on the terminal:
```root -l 
.L Compilemyclass;
Compylemyclass()
Simulation(1000)```
