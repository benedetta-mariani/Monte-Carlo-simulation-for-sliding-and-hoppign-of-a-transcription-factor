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

The parameters in the simulations are set to some suggested in the paper, since they fit experimental data: these values of parameters (s<sup>2</sup> = 700 bp<sup>2</sup>, rmax = 11 nm) imply that hopping has low probability to happen, and sliding together with macroscopic dissociations seems to play major contributions in the search of the site. 

Other parameters that can be tried and that fit experimental data according to the paper (Fig. 3b) are s<sup>2</sup> = 700 bp<sup>2</sup>, rmax = 20 nm (so it is necessary just to change the value of ```rmax``` in the file *Simulation.cpp* to 20.)
This value of ```rmax``` is considered more biologically plausible by the authors of the paper, and allows for more probable hopping events. Even in this case however the authors argue that hopping turns out not to play a significant role in the correlated search of the site. Indeed, even if the TF often dissociates from DNA to hop, it does often return to the same base pair, as also my histograms of hops lengths show (an example is contained in the file *hopslength.png*, where the lengths result peaked in zero).  Therefore the TF effectively stays in close proximity to DNA during many short hopping events.
