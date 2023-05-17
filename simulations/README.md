# Segregational drift hinders the evolution of antibiotic resistance on polyploid replicons

## File descriptions in simulation folder

### Simulation modules

* ``module.py`` Python module containing function to simulate dynamics  
* ``data.py`` Python module containing function to load data from experimental data csv-file  
* ``plot.py`` Python module containing functions for plotting simulation data   
</br>

### Simulations under neutral conditions (Figure 1D)
* ``neutral.ipynb`` Jupyter Notebook containing simulation scripts and plotting scripts to produce Figure 1D  
* ``neutral-SI.ipynb`` Jupyter Notebook for neutral simulations with different copy numbers to produce SI Figures.

### Simulations under selective conditions (Figure 3)

* ``fit-plasmid.ipynb``, ``fit-chromosome.ipynb`` Jupyter notebook used for fitting the experimental data to the mathematical model to determine the selection coefficient for the plasmid experiment and the chromosome experiment respectively.  
</br>
* ``selection.ipynb`` Jupyter Notebook containing simulation scripts and plotting scripts to produce Figure 3  
* ``selection_cluster-batch.sh`` Slurm batch script for computer cluster 
* ``selection_input/`` Folder containing input files with one simulation parameterset each created by ``selection.ipynb``
* ``selection_output/`` Folder for output files from calculations on computer cluster
* ``selection-SI...`` Jupyter Notebook, batch script, input and output files for simulations with $s=0.1$ instead of $s$-value estimated from experimental data to produce SI Figure.  
</br>

### Simulations for establishment probabilities (Figure 4)

* ``establishment.ipynb`` Jupyter Notebook containing simulation scripts and plotting scripts to produce Figure 4  
</br>

### Experimental data

* ``expdata/`` Folder containing experimental data files 
  * ``Nc_chromosome.csv``, ``Nc_plasmid.csv`` csv-file of carrying capacities for chromosome and plasmid experiment



## Installation & Cluster
The Python modules & notebooks rely on the following Python packages: 
numpy scipy pandas matplotlib openpyxl lmfit seaborn ipywidgets colour
(Note that lmfit is not part of the official Anaconda package libraries)

