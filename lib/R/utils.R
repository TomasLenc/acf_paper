
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
    ggsave(paste(fname, '.png', sep=''), plt, width=width, height=height, bg='white')
    # save as svg
    ggsave(paste(fname, '.svg', sep=''), plt, width=width, height=height,  bg='white')
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


render_sigmoid_fit <- function(m){
    tmp <- as.data.frame(confint(m))
    tmp <- rownames_to_column(tmp, 'parameter')
    tmp$parameter <- str_replace(tmp$parameter, '_', '')
    tmp
}


posthoc_ttest_low_vs_high <- function(df, var_name_eeg){
    vals_L <- df %>% 
        filter(tone=='L') %>% 
        arrange(subject) %>%
        pull(.data[[var_name_eeg]])
    vals_H <- df %>% 
        filter(tone=='H') %>% 
        arrange(subject) %>%
        pull(.data[[var_name_eeg]])
    res <- t.test(vals_L, vals_H, paired=TRUE)
    data.frame(mean_L_minus_H = unname(res$estimate), 
               t = unname(res$statistic),
               df = unname(res$parameter),
               p = unname(res$p.value))
}


ttest_0 <- function(df, var_name){
    vals_eeg <- df %>% pull(.data[[var_name]])
    res <- t.test(vals_eeg, mu=0, alternative='greater')
    data.frame(t=res$statistic, p=res$p.value)
}

