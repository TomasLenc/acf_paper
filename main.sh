#!/usr/bin/env bash

/home/tomo/.local/bin/matlab_run_terminal.py main_only_noise.m && \
/home/tomo/.local/bin/matlab_run_terminal.py main_snr_vs_nlags.m &&  \
/home/tomo/.local/bin/matlab_run_terminal.py main_syncrange_eeg_boot.m && \
/home/tomo/.local/bin/matlab_run_terminal.py main_syncrange_eeg.m && \
/home/tomo/.local/bin/matlab_run_terminal.py main_syncrange_tapping.m 
