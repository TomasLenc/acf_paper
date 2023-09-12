
read_features_file <- function(response, suffix=NA, area=NA, func_rois=NA){
    fname <- sprintf('response-%s', response)
    if (!is.na(area)) fname <- sprintf('%s_area-%s', fname, area)
    if (any(!is.na(func_rois))){
        if (!is.list(func_rois)) func_rois <- list(func_rois)
        for (fr in func_rois) {
            fname <- sprintf('%s_funcROI-%s', fname, fr)
        }
    }    
    if (!is.na(suffix)) fname <- sprintf('%s_%s', fname, suffix)
    fname <- sprintf('%s.tsv', fname)
    df <- fix_factors(read.csv(file.path(csv_path, fname), sep='\t'))
    return(df)
}

fix_factors <- function(df){
    # This function takes care of 
    # (1) changing char columns to factors
    # (2) renaming factor names and levels
    # (3) reordering factor levels if needed
    if ('subject' %in% names(df)){
        df$subject <- factor(df$subject)
    }
    if ('response' %in% names(df)){
        df$response <- factor(df$response)
    }
    if ('ses' %in% names(df)){
        df$ses <- factor(df$ses)
    }
    if ('rhythm' %in% names(df)){
        # recode_factor() remaps factor levels, and also changes their order 
        df$rhythm <- recode_factor(df$rhythm, !!!rhythm_label_map)
    }
    if ('area' %in% names(df)){
        df$area <- factor(df$area)
    }
    return(df)
}


make_frex_table <- function(df_frex){
    # This function takes a dataframe with frequencies of interest, and returns a formatted table 
    # that can be printed with pander. 
    df2print <- df_frex
    names(df2print) <- c('hamonic', 'frequency (Hz)', 'is_meter_related', 'is_meter_unrelated')
    df2print$is_meter_related <- ifelse(df2print$is_meter_related==1, '✓', '')
    df2print$is_meter_unrelated <- ifelse(df2print$is_meter_unrelated==1, '✓', '')
    return(df2print)
}


make_lagz_table <- function(df_acf_ind){
    # This function takes a dataframe with ACF at individual lags and returns a formatted table 
    # that can be printed with pander, containing unique meter-related and meter-unrealted lags.
    lags_meter_rel <- unique(df_acf_ind %>% filter(meter_rel == 1) %>% select(lag) %>% unlist())
    lags_meter_unrel <- unique(df_acf_ind %>% filter(meter_rel == 0) %>% select(lag) %>% unlist())
    lags <- c(lags_meter_rel, lags_meter_unrel)
    is_meter_rel <- c(rep(TRUE, length.out=length(lags_meter_rel)),
                      rep(FALSE, length.out=length(lags_meter_unrel)))
    idx = order(lags)
    df2print <- data.frame(
        lag = lags[idx],
        is_meter_related = is_meter_rel[idx],
        is_meter_unrelated = !is_meter_rel[idx]
    )
    df2print$is_meter_related <- ifelse(df2print$is_meter_related==1, '✓', '')
    df2print$is_meter_unrelated <- ifelse(df2print$is_meter_unrelated==1, '✓', '')
    return(df2print)
}



save_fig <- function(fname, plt, width, height){
    # strip file extension
    fname <- tools::file_path_sans_ext(fname)
    # # save as png
    # ggsave(paste(fname, '.png', sep=''), plt, width=width, height=height, bg='white')
    # save as svg
    ggsave(paste(fname, '.svg', sep=''), plt, width=width, height=height,  bg='white')
}


get_intercept_from_lm <- function(m_summary) {
    if ('t value' %in% colnames(m_summary$coefficients)) {
        t <- round(m_summary$coefficients['(Intercept)', 't value'], 2)
    } else {
        t <- NA
    }
    if ('df' %in% colnames(m_summary$coefficients)) {
        df <- round(m_summary$coefficients['(Intercept)', 'df'], 2)
    } else {
        df <- NA
    }
    if ('Estimate' %in% colnames(m_summary$coefficients)) {
        beta <- round(m_summary$coefficients['(Intercept)', 'Estimate'], 2)
    } else {
        beta <- NA
    }
    if ('Pr(>|t|)' %in% colnames(m_summary$coefficients)) {
        p <- round(m_summary$coefficients['(Intercept)', 'Pr(>|t|)'], 6)
    } else {
        p <- NA
    }
    
    return(c(t=t, df=df, beta=beta, p=p))
}


ses_to_date <- function(ses) {
    
    ses_split <- str_split(ses, '_')
    
    ses_split <- lapply(ses_split, function(x) x[1])
    
    ses_dates <- lapply(ses_split, function(x) {
            if (nchar(x) != 8) return(NA)
            return(as_date(x, format='%d%m%Y'))
            }
        )
    
    return(as_date(unlist(ses_dates)))
}


get_min_n_responsive <- function(df, min_over_factor=NULL){
    # Get minimum responsive electrodes per rhythm per session. If one rhythm is completely missing, 
    # return 0
    if (is.null(min_over_factor)){
        min_n_conditions <- 1
    } else {
        min_n_conditions <- length(levels(df[, min_over_factor]))
    }
    print(sprintf('putting "0" if less than %d values per session found', min_n_conditions))
    get_min_n <- function(d) {
        min_n <- min(d$n)
        if (nrow(d) < min_n_conditions) min_n <- 0
        return(data.frame(n=min_n))
    } 
    return(df %>% 
               group_by(across(any_of(c('ses', min_over_factor)))) %>% 
               tally() %>% 
               group_by(ses) %>% 
               group_modify(~get_min_n(.x))
    )
}


get_boot_summary <- function(df, grouping_val, feature_name){
    # summarize boostrapped samples to get confidence intervals
    data.frame(
        ci_low = quantile(df[, feature_name] %>% pull(), 0.025),
        ci_high = quantile(df[, feature_name] %>% pull(), 1-0.025)
    )
}


r_to_z <- function(r){
    # fisher-transform
    .5 * (log(1+r) - log(1-r))
}


# sigm <- function(x, params) {
#     {params['min_'] + (params['max_'] - params['min_']) / (1 + 10.^((params['x50_'] - x) * params['slope_']))}
# }
sigm <- function(x, min_, max_, x50_, slope_) {
    # sigmoid function 
    {min_ + (max_ - min_) / (1 + 10.^((x50_ - x) * slope_))}
}

fit_sigm <- function(x, y, init_params=NULL) {
    # fit sigmoid function to data and return fitted model object
    df <- data.frame(x=x, y=y)
    if (is_empty(init_params)){
        init_params <- c(quantile(y, 0.05), quantile(y, 0.95), NA, 1)
        if (sum(y == quantile(y, 0.5)) == 0){
            temp <- x[y == quantile(y[2:length(y)], 0.5)]
        } else {
            temp <- x[y == quantile(y, 0.5)]
        }
        init_params[3] <- temp[1]
        names(init_params) <- c('min_', 'max_', 'x50_', 'slope_')
        init_params <- as.list(init_params)        
    }
    m <- nls(y ~ I(min_ + (max_ - min_) / (1 + 10.^((x50_ - x) * slope_))), 
             data=df,
             start=init_params)
    return(m)
}
# ## TEST ##
# x <- runif(100)
# y <- sigm(x, 0, 1, 0.3, 10) + rnorm(100)*0.04
# df <- data.frame(x=x, y=y)
# plot(df)
# params <- fit_sigm(x, y)




