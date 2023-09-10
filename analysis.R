rm(list=ls())

# get general parameters
source('config.R')

data_dirs <- c(
    'maxlag-halfTrial_meterUnrel-0.2', 
    'maxlag-2.4_meterUnrel-0.2', 
    'maxlag-halfTrial_meterUnrel-0.6_1.0_1.4',
    'maxfreq-5_excl5', 
    'maxfreq-30_excl5',
    'maxfreq-5_excl0.416,0.833'
)

i_data_dir <- 1

for (i_data_dir in c(1:length(data_dirs))) {
    
    data_path = file.path(experiment_path, 'data', data_dirs[i_data_dir])
    fig_path = file.path(experiment_path, 'figures', data_dirs[i_data_dir])
    
    if (! dir.exists(fig_path)) {
        dir.create(fig_path)
    }

    # Noise
    results_fname <- sprintf('report_noise')
    rmarkdown::render(sprintf('report_noise.Rmd'),
                      output_file=sprintf('%s/%s.html', data_path, results_fname))

    # Syncrange
    results_fname <- sprintf('report_syncrange')
    rmarkdown::render(sprintf('report_syncrange.Rmd'),
                      output_file=sprintf('%s/%s.html', data_path, results_fname))

    # N-lags vs. SNR
    results_fname <- sprintf('report_nalgsVsSnr')
    rmarkdown::render(sprintf('report_nlags_nfrex.Rmd'), 
                      output_file=sprintf('%s/%s.html', data_path, results_fname))

    
}


# FOOOF vs IRASA
data_path = file.path(experiment_path, 'data')
results_fname <- sprintf('report_fooofVsIrasa')
rmarkdown::render(sprintf('report_fooof_irasa.Rmd'),
                  output_file=sprintf('%s/%s.html', data_path, results_fname))



# infant
data_path = file.path(experiment_path, 'data')
results_fname <- sprintf('report_infant')
rmarkdown::render(sprintf('report_infant.Rmd'),
                  output_file=sprintf('%s/%s.html', data_path, results_fname))

# lowhigh
ap_fit_method <- 'irasa'
roi_name <- 'frontocentral'
data_path = file.path(experiment_path, 'data')
results_fname <- sprintf(sprintf('report_xp-lowhigh_apFitMethod-%s_roi-%s',
                                 ap_fit_method, roi_name))
rmarkdown::render(sprintf('report_lowhigh.Rmd'),
                  output_file=sprintf('%s/%s.html', data_path, results_fname))

# attention
ap_fit_method <- 'irasa'
roi_name <- 'all'
data_path = file.path(experiment_path, 'data')
results_fname <- sprintf(sprintf('report_xp-attention_apFitMethod-%s_roi-%s',
                                 ap_fit_method, roi_name))
rmarkdown::render(sprintf('report_attention.Rmd'),
                  output_file=sprintf('%s/%s.html', data_path, results_fname))




