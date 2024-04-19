Code that runs simulations and analyses for the paper: 
"**Measuring self-similarity in empirical signals to understand musical beat perception**"



## Requirements


* MATLAB >= 2018a
* R >= 4.3.2
* [acf_tools](https://github.com/TomasLenc/acf_tools)
* [rnb_tools](https://github.com/TomasLenc/rnb_tools)
* [letswave6](https://github.com/NOCIONS/letswave6)

The code was tested with Ubuntu 20.04 and macOS 13.6. 



## Data

The analyses depend on several empirical datasets:  

<br>

1. **A test-retest resting and cognitive state EEG dataset**

The data is available on [openneuro](https://openneuro.org/datasets/ds004148) (number ds004148).  

The data can be conveniently obtained with [datalad](https://www.datalad.org/). Navigate to the project directory using the command line and the run the following: 

```bash
# install the dataset in current directory
datalad install https://github.com/OpenNeuroDatasets/ds004148.git

# only get the eyes-open data from session 1: 
datalad get sub-*/ses-session1/eeg/sub-*_ses-session1_task-eyesopen*

# in order to load without error, you need to unlock the files!!!
datalad unlock sub-*/ses-session1/eeg/sub-*_ses-session1_task-eyesopen*
```


<br> 

2. **InfantRhythm dataset** 

Raw data available on [OSF](https://osf.io/9wf5u). Here, we already provide preprocessed data in a simple `letswave6` format on project [OSF](https://osf.io/5s3j7).   



<br> 

3. **Lowhigh dataset**

The data can be downloaded from project [OSF](https://osf.io/5s3j7). 





## Directory structure 

First, set up the paths to all required libraries, and the base path for the project in `get_par.m`. 

The variable `experiment_path` should point to the project directory that has the following structure (make sure you download all the required data):  


```
experiment_path
└── data
    ├── infant
    │   ├── cochlear_model
    │   │   └── Slaney_128coch_meddis_timeDomain_meanAcrossCF.mat
    │   └── eeg
    │       ├── AB high sync P001.lw6
    │       ├── AB high sync P001.mat
    │       ├── AB high sync P002.lw6
    │       ├── AB high sync P002.mat
    │       ├── AB high sync P003.lw6
    │       ├── AB high sync P003.mat
    │       ├── AB high sync P004.lw6
    │       ├── AB high sync P004.mat
    │       ├── AB high sync P005.lw6
    │		...
    ├── lowhigh
    │   ├── cochlear_model
    │   │   └── Slaney_128coch_meddis_timeDomain_meanAcrossCF.mat
    │   ├── eeg
    │   │   ├── H_syncopated.lw6
    │   │   ├── H_syncopated.mat
    │   │   ├── H_unsyncopated.lw6
    │   │   ├── H_unsyncopated.mat
    │   │   ├── L_syncopated.lw6
    │   │   ├── L_syncopated.mat
    │   │   ├── L_unsyncopated.lw6
    │   │   └── L_unsyncopated.mat
    │   └── tapping
    │       ├── H_syncopated.mat
    │       ├── H_unsyncopated.mat
    │       ├── L_syncopated.mat
    │       └── L_unsyncopated.mat
    └── resting_eeg
        └── ds004148
            ├── CHANGES
            ├── README
            ├── dataset_description.json
            ├── derivatives
            │   ├── README
            │   ├── code
            │   ├── figures
            │   ├── five-states_events.json
            │   └── preprocessed data
            ├── participants.json
            ├── participants.tsv
            ├── sub-01
            │   ├── ses-session1
            │   ├── ses-session2
            │   └── ses-session3
            ├── sub-02
            │   ├── ses-session1
            │   ├── ses-session2
            │   └── ses-session3
            ...

```


## Run analyses


Once everything is set up, run `main.m`. This will run all the simulations and save the data to disk. 

Once the simulations are run, generate the plots of the results using `main_plots.m`. 

Other plots required to assemble the figures in the paper are generated by individually running the scripts named `explain_<description>.m`. 

Finally, run `analysis.R` in R to generate the rest of the figures and do statistical tests. A .html report will be generated and figures also saved to disk ask image files.  

All the output will be saved into the `<experiment_path>/results` directory. 





