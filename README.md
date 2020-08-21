# Monte Carlo Simulation for sliding and hopping

The code here provided aims at simulating the search of its specific binding site by a transcription factor; the search consists in a 3D diffusion (**hopping**) and in a 1D diffusion along the DNA chain (**sliding**). The model is reproduced following the article https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.022721.
All the details of the model, the organization of the C++ classes and how to execute the simulation are provided in the file *Relazione.pdf*.

To run the code it is necessary the ROOT framework. Run on the terminal:
```
root -l 
.L Compilemyclass.C
compilemyclass();
simulation(1000);
```
where here the parameter ```searches``` is set to 1000 (number of independent searches of the binding site one wants to perform, in order to obtain average values). When ```searches``` is set to 1000, good agreement with the experimental expectations starts to appear (ten minutes of execution time). However the results presented in the pdf are obtained with ```searches = 10*5```, which requires a nine hours long execution.

The parameters in the simulations are set to the one suggested in the paper, since they fit experimental data: these values of parameters (s*2 = 700) imply that hopping has very low probability to happen, and sliding together with macroscopic dissociations seems to play major contributions in the search of the site. 
