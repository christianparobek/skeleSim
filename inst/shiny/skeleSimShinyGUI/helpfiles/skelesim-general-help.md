This tab is required to set up some general parameters upon which all simulations depend.
 1. Title: this is used to create names for output files.  It is a good idea to use something other than the default
 2. Date: self-explanatory
 3. Quiet: print less during simulation.  Right now this does not actually change the amount of information produced
 4. 'Coalescent simulator?': this checkbox indicates whether the simulation is a coalescent reverse time simulation (checked) or a forward-time individual-based simulation (unchecked)
 5. 'Number of simulation reps': this is the number of replicates for each simulated scenario
 6. 'types of analyses requested': These checkboxes select analyses that are performed on every simulation rep
 7. 'number of permutations': Several of the analysis functions perform significance tests using permutation. This specifies the number of reps.  If you don't care about significance tests, set to 0 (default). Otherwise choose a number like 100, while remaining aware that this will slow down analyses of simulated datasets.
 8. 'Temporary subdirectory': This is the subdirectory in which simulations are executed. It is created below the 'root directory for simulations' chosen in the 'Actions' tab
 
