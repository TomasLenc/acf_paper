Code that runs simulations for the paper: 
"**Measuring self-similarity in empirical signals to understand musical beat perception**"

Dependencies: 
* MATLAB >= 2018a
* R >= 4.3.2
* [acf_tools](https://github.com/TomasLenc/acf_tools)
* [rnb_tools](https://github.com/TomasLenc/rnb_tools)
* [letswave6](https://github.com/NOCIONS/letswave6)



The code was tested with Ubuntu 20.04 and macOS 13.6. 

First, set up the paths to all required libraries, and the base path for the project in `get_par.m`. 

Then, run `main.m`. This will run all the simulations and save the data to disk. 

Once the simulations are run, generate the plots of the results using `main_plots.m`. 

Finally, run `analysis.R` in R to generate rest of the figures and do statistical tests. 