# Monte Carlo Simulation for sliding and hopping

The code here provided aims at simulating the search of its specific binding site by a transcription factor (TF); the search consists in a 3D diffusion (**hopping**) and in a 1D diffusion along the DNA chain (**sliding**). The model is reproduced following the article https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.022721.
All the details of the model, the organization of the C++ classes and how to execute the simulation are provided in the file *Relazione.pdf*.

To run the code it is necessary the ROOT framework. Run on the terminal:
```
root -l 
.L Compilemyclass.C
compilemyclass();
simulation(1000);
```
where here the parameter ```searches``` is set to 1000 (number of independent searches of the binding site one wants to perform, in order to obtain average values). When ```searches``` is set to 1000, good agreement with the experimental expectations starts to appear (ten minutes of execution time). However the results presented in the pdf are obtained with searches = 10<sup>5</sup>, which requires a nine hours long execution.

The parameters in the simulations are set to some suggested in the paper, since they fit experimental data: these values of parameters (s<sup>2</sup> = 700 bp<sup>2</sup>, rmax = 11 nm) imply that hopping has very low probability to happen, and sliding together with macroscopic dissociations seems to play major contributions in the search of the site. 

Even in the case of other values of parameters that fit experimental data (s<sup>2</sup> = 10<sup>2</sup>, rmax = 30 nm), which are considered more biologically plausible by the authors of the paper and that allow for more probable hopping events, hopping turns out not to play a significant role in the correlated search of the site, since even if the TF often dissociates from DNA to hop, it does return to the same base pair and therefore effectively stays in close proximity to
DNA during many short hopping events.
