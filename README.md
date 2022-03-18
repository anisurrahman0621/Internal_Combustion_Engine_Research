# Internal_Combustion_Engine_Research

### About

I conducted research for the internal combustion engine group while in grad school. I was tasked with finding the effects of knock on engine thermodynamics and identifying the modes of efficiency loss. Engine's gain thermal efficiency under various operating conditions such as operating lean (higher air:fuel ration), igniting spark plugs earlier, and using higher compression ratios. All of these produce higher peak pressures which, on average, produces favorable thermodynamics for internal combustion engines. However, there comes a point of diminishing returns. At some point, the pressures become too high when the fuel starts "knocking". This phenomona occurs when the fuel ignites before the flame propogated by the spark plug has reached it. This causes the pressure to rapidly rise. The thermal efficiencies and fuel consumptions all noteiceably suffer. However, earlier it was stated that higher pressures were favorable. I evaluated many different thermodynamic properties to explain how this counter-intuitive result comes to be. 

### Process

The professor operated an engine at 11 degrees before top dead center (TDC) spark timing at a compression ratio of 10. Operating at these conditions induced knocking but not for every cycle. The cycles had to be seperated into 2 categories-- cycles that knocked and cycles that didn't. This was done by using an ignition delay equation for each cycle.   

### MATLAB Program

The pressure and volume were recorded every 0.2 degrees of the crank shaft for 200 cycles resulting in 720,000 data points. The data points came as one 720,000 x 1 array. Various post-processing techniques were used to analyze the large data set and calculate the relevant thermodynamic properties using equations learned in class. The final report contains many different plots and charts to meticulously explain the findings of the code. The PDF of the final report is uploaded to the main folder.

### Missing Files

The original csv files containing the pressures are not uploaded to the repository. I took this class a few years ago and had to change computers since then. I didn't fully recover all of the files I needed to. I also lost some of the functions/methods I created to simplify reptitive codes. I no longer have accesss to my courses from school so the program won't run on a different computer. 
