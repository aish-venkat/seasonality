# Function to calculate harmonic peaks and nadirs from weekly, monthly, or daily time series
# Last Updated: April 2, 2023

# Use pacman package to load/install needed packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, imputeTS, forecast, multitaper, Rssa)

######################################################
# REFERENCES #########################################
######################################################

# https://stackoverflow.com/questions/41435777/perform-fourier-analysis-to-a-time-series-in-r

# Percival, D.B. and Walden, A.T. (1993) 
# Spectral analysis for physical applications. Cambridge University Press.

# Alsova, O. K., Loktev, V. B., & Naumova, E. N. (2019). 
# Rotavirus seasonality: An application of singular spectrum analysis and 
# polyharmonic modeling. IJERPH 16(22), 4309.

# Simpson, R. B., Zhou, B., & Naumova, E. N. (2020). Seasonal synchronization of 
# foodborne outbreaks in the United States, 1996–2017. Scientific reports, 10(1), 17500.
# https://www.nature.com/articles/s41598-020-74435-9
# Equations for negative binomial regression in Supplementary Table 3

######################################################
# SUPPORTING FUNCTIONS ###############################
######################################################

# Function to return amplitude and frequency from FFT with optional upsampling interval 
nff = function(x = NULL, n = NULL, upsample = FALSE, up = 10L){
  # Get time index of time series
  t <- 1:length(x); N <- length(x);
  #The direct transformation
  #The first frequency is DC, the rest are duplicated
  dff = fft(x)
  #The time
  t = seq(from = 1, to = length(x))
  
  if(isTRUE(upsample)){
    
    # Upsampled time
    nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
    #New spectrum
    ndff = array(data = 0, dim = c(length(nt), 1L))
    ndff[1] = dff[1] #First component = average (DC component), not meaningful
    if(n != 0){
      ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
      #The negative ones are trickier
      ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
    }
    # Add noise up to +/- 0.01 for each time period
    ndff_n <- ndff + 0.5*rnorm(N, sd=0.01)
    time = nt;
    
  } else{
    ndff_n <- dff; time = t;
  }
  
  # The inverses
  indff = fft(ndff_n/length(ndff_n), inverse = TRUE)
  #idff = fft(dff/N, inverse = TRUE) # reconstructed input time series
  
  # Amplitude and phase of inverse
  amp <- Mod(indff); pha <- Arg(indff); 
  
  # Frequency
  freq <- (1:(length(ndff_n)))*(omega)/(2*length(ndff_n))
  
  ret = data.frame(time = time, y = Mod(indff), freq = freq, amplitude = amp, phase = pha)
  return(ret)
}

# Function to convert from polar to cartesian coordinates
polar2cart <- function(r, theta) {
  data.frame(x = r * cos(theta), y = r * sin(theta))
}

# Function to convert from cartesian to polar coordinates
cart2polar <- function(x, y) {
  data.frame(r = sqrt(x^2 + y^2), theta = atan2(y, x))
}

######################################################
# MAIN FUNCTIONS ######################################
######################################################

# Function to calculate peak/nadir timings and values 
peaktimecalc <- function(mod, omega){
  
  # New data for which to predict
  preddf <- data.frame(INDEX=seq(1, omega, by=0.125))
  preddf <- preddf %>%
    mutate(SIN2PI = sin(2*pi*(1/omega)*INDEX),
           COS2PI = cos(2*pi*(1/omega)*INDEX),
           SIN4PI = sin(4*pi*(1/omega)*INDEX),
           COS4PI = cos(4*pi*(1/omega)*INDEX)) %>%
    mutate(INDEX2 = INDEX^2, INDEX3 = INDEX^3) 
  
  if(any(grepl("SIN4PI", mod$call))){
    
    # Ref: https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
    preddf <- preddf %>%
      # Develop a predictive model
      mutate(PRED = as.vector(predict(mod, newdata=preddf, type="response"))) %>%
      # Calculate derivatives
      mutate(dydx_1 = c(NA, diff(PRED)), dydx_2 = c(NA, NA, diff(diff(PRED)))) %>% 
      # Calculate peak values
      mutate(PEAKS = c(NA, (diff(sign(diff(PRED))<0)), NA )) %>% 
      mutate(MAXIMA = ifelse(PEAKS==1, 1, 0), MINIMA = ifelse(PEAKS==-1, 1, 0)) 
    
    peaktiming <- preddf %>% filter(MAXIMA==1) %>% arrange(-PRED) %>% pull(INDEX)
    peakvalue <- preddf %>% filter(MAXIMA==1) %>% arrange(-PRED) %>% pull(PRED)
    
    nadirtiming <- preddf %>% filter(MINIMA==1) %>% arrange(-PRED) %>% pull(INDEX)
    nadirvalue <- preddf %>% filter(MINIMA==1) %>% arrange(-PRED) %>% pull(PRED)
    
    fft_pred <- nff(preddf$PRED, upsample=F, n=0)
    
    pt <- data.frame(MODEL = "4PI", 
                     PEAKTIMING = peaktiming[1], PEAKVALUE = peakvalue[1],
                     NADIRTIMING = nadirtiming[1], NADIRVALUE = nadirvalue[1],
                     PEAK2TIMING = peaktiming[2], PEAK2VALUE = peakvalue[2],
                     NADIR2TIMING = nadirtiming[2], NADIR2VALUE = nadirvalue[2],
                     AMPLITUDE = max(fft_pred$amplitude) - min(fft_pred$amplitude),
                     AMP_CI = sd(fft_pred$amplitude)) 
    
  } else {
    
    M <- (omega/(2*pi))
    coefs<-as.vector(coef(mod))
    intercept<-coefs[grep("(Intercept)", names(coef(mod)))[1]] 
    coef_sin<-coefs[grep("sin", names(coef(mod)), ignore.case=TRUE)[1]] 
    coef_cos<-coefs[grep("cos", names(coef(mod)), ignore.case=TRUE)[1]] 
    
    # Amplitude
    amp <- sqrt(coef_sin^2 + coef_cos^2)
    
    # CI of amplitude
    sd_sin <- sd(timedf$SIN2PI) 
    sd_cos <- sd(timedf$COS2PI)
    sd_sincos <- sd(timedf$SIN2PI * timedf$COS2PI)
    
    var_amp <- ((coef_cos^2 + sd_cos^2) + (coef_sin^2 + sd_sin^2) + (2* sd_sincos * beta_sin * beta_cos)) / (beta_sin^2 + beta_cos^2)
    ci_amp <- 1.96 * sqrt(var_amp)
    
    # Angle
    ang <- coef_sin / coef_cos
    
    # Phi or Shift 
    shift <- atan(ang)
    
    # Variance of phi
    var_phi <- ( (coef_cos^2 + sd_sin^2) + (coef_sin^2 + sd_cos^2) - (2* sd_sincos * beta_sin * beta_cos) ) / ( (beta_sin^2 + beta_cos^2)^2 )
    
    # Peak Timing
    peaktiming <- ifelse(coef_cos < 0, ((shift + pi)*M),
                         ifelse(coef_cos > 0 & coef_sin > 0, (shift*M),
                                ((shift + (2 * pi))*M) ))
    
    # CI of peak timing
    ci_peaktiming <- 1.96 * sqrt(var_phi) * (M / (2 * pi))
    
    # Peak Value
    peakvalue = intercept + amp
    
    # Nadir Timing
    nadirtiming <- ifelse(coef_cos < 0, shift + (pi/2),
                          ifelse(coef_cos > 0 & coef_sin > 0, (shift/2), (shift + (2 * pi))/2 ))
    
    # Nadir Value
    nadirvalue = intercept - amp
    
    # Intensity 
    intensity = (intercept + amplitude) - (intercept - amplitude)
    var_intensity = 4*(var_amp)
    ci_intensity = intensity + (1.96 * sqrt((var_intensity)))
    
    pt <- data.frame(MODEL = "2PI", 
                      PEAKTIMING = peaktiming, PEAKVALUE = peakvalue,
                      NADIRTIMING = nadirtiming, NADIRVALUE = nadirvalue,
                      PEAK2TIMING = NA, PEAK2VALUE = NA,
                      NADIR2TIMING = NA, NADIR2VALUE = NA,
                      AMPLITUDE = amp, AMP_CI = amp_ci, 
                      INTENSITY = intensity, INTENSITY_CI = ci_intensity)
    
  }
  
  # Return estimated peak timings and predicted values
  return(list(pt, preddf))
  
}

# Function to calculate seasonality 
seasonalitycalc <- function(df, tfield, omega, outcome, transform_ln = TRUE){

  # Process time fields
  if(omega==12){
    omegastr = "months"
  } else if(omega==365 | omega==365.25){
    omegastr = "days"
  } else if(omega==52){
    omegastr = 52
  }
  
  timefield <- df %>% pull(tfield)
  
  if(class(timefield)=="Date"){
    
    # Round to beginning of year
    mindate <- min(timefield, na.rm=T); mindate <- floor_date(mindate, unit='year');
    maxdate <- max(timefield, na.rm=T); maxdate <- ceiling_date(maxdate, unit='year');
    
    # Create sequential series
    indexdf <- data.frame(DATE = seq.Date(mindate, maxdate, by=omegastr)) %>%
      mutate(INDEX = 1:nrow(.), INDEX2 = INDEX^2, INDEX3 = INDEX^3,
             SIN2PI = sin(2*pi*(1/omega)*INDEX), 
             COS2PI = cos(2*pi*(1/omega)*INDEX),
             SIN4PI = sin(4*pi*(1/omega)*INDEX), 
             COS4PI = cos(4*pi*(1/omega)*INDEX))
    
    # Merge sequence into time series
    timedf <- data.frame(DATE = timefield, value = df %>% pull(outcome)) %>% 
      full_join(indexdf, by="DATE") %>% 
      arrange(DATE) %>%
      distinct()
    
  } else if(class(timefiled)=="Numeric"){
    
    # Create sequential series
    
    mindate <- min(timefield, na.rm=T); maxdate <- max(timefield, na.rm=T)
    indexdf <- data.frame(DATE = seq(mindate, maxdate, by=omega),
                          value = df %>% pull(outcome)) %>%
      mutate(INDEX = 1:nrow(.)) %>%
      mutate(INDEX = 1:nrow(.), INDEX2 = INDEX^2, INDEX3 = INDEX^3,
             SIN2PI = sin(2*pi*(1/omega)*INDEX), 
             COS2PI = cos(2*pi*(1/omega)*INDEX),
             SIN4PI = sin(4*pi*(1/omega)*INDEX), 
             COS4PI = cos(4*pi*(1/omega)*INDEX))
    
    # Merge sequence into time series
    timedf <- data.frame(DATE = timefield, value = df %>% pull(outcome)) %>% 
      full_join(indexdf, by="DATE") %>% 
      arrange(DATE) %>% distinct()
    
  }
  
  # Log transform if needed
  if(isTRUE(transform_ln)){
    timedf <- timedf %>% mutate(value = log(value))
  }
    
  # 1. Impute missing values and create time series -----------
  vals <- ts(timedf$value,
             start = c(year(timedf$DATE[1]), month(timedf$DATE[1])),
             deltat = 1/omega)
  
  vals_filled <- na_seadec(vals, find_frequency = T) # in log scale
  N <- length(vals_filled)
  
  # 2. Detrend data ---------------
  
  # Detrending removes the first eigenvector which best captures complex trend.
  # This step is retained as we are interested in seasonality
  
  vals_stl <- stl(vals_filled, s.window='periodic')
  
  # Add these values back into dataset
  timedf$seasonal <- vals_stl$time.series[,1]
  timedf$trend <- vals_stl$time.series[,2]
  timedf$remainder <- vals_stl$time.series[,3]
  timedf$detrend <- vals_filled - vals_stl$time.series[,2]
  
  vals_detrend <- vals_filled - vals_stl$time.series[,2] # in log scale
  
  # 3. Spectral Analysis ----------------------
  
  # s <- spec.mtm(vals_detrend, Ftest = TRUE, jackknife = T)
  # 
  # # calculate spectrum
  # spec1 <- s$spec[1]
  # s.df <- data.frame(freq = s$freq,  spec = s$spec, spec_scaled = s$spec/spec1, 
  #                    f = s$mtm$Ftest, var_jk = s$mtm$jk$varjk,
  #                    ci_lower = s$mtm$jk$lowerCI, ci_upper = s$mtm$jk$upperCI) %>%
  #   # Multiply spectra by 2
  #   #mutate(spec = 2*spec) %>%
  #   # Remove zero value as it is meaningless
  #   dplyr::filter(freq!=0) %>%
  #   # Convert from frequency to periods 
  #   mutate(period = freq / (1/omega),
  #          time = round(period - (omega * (period %/% omega)), 0) )
  # 
  # # Which frequencies are statistically significant?
  # sig_months <- s.df %>% dplyr::filter(f <= 0.05)
  # 
  # # Are same frequencies statistically significant across cycles?
  # # ggplot(data = s.df, aes(x=period, y=spec)) + geom_line() +
  # #   geom_vline(xintercept = s.df %>% filter(f <= 0.05) %>% pull(period), lty=2, col='red') +
  # #   scale_x_continuous("Period (years)", breaks = seq(0, length(vals_detrend)*2, by=omega)) +
  # #   scale_y_log10()
  # # 
  # # yrs.period <- rev(c(1/12, 1/10, 1/8, 1/6, 1/5, 1/4, 1/3, 0.5, 1, 2))
  # # yrs.labels <- rev(c("1/12", "1/10", "1/8", "1/6", "1/5", "1/4", "1/3", "1/2", "1", "2"))
  # # yrs.freqs <- (1/yrs.period) * (1/omega)  # Convert annual period --> annual freq --> monthly freq
  # # ggplot(data = s.df %>% filter(freq>1/12), aes(x=freq, y=spec)) + geom_line() +
  # #   geom_vline(xintercept = s.df %>% filter(f <= 0.05) %>% pull(freq), lty=2, col='red') +
  # #   scale_x_log10("Period (years)", breaks = yrs.freqs, labels = yrs.labels) +
  # #   scale_y_log10()
  
  # 4. Upsample original signal using FFT + 10th harmonic for best fit --------------------
  vals_up <- nff(vals_filled, n=10, upsample=TRUE) 
  
  # Create new signal from upsampled time series
  vals_up <- ts(vals_up$y,
                start = c(year(timedf$DATE[1]), month(timedf$DATE[1])),
                frequency = omega/0.125)
  
  # 5. Singular Spectrum Analysis -------------------
  
  vals_ssa <- ssa(vals_up)
  
  # plot(vals_ssa, type='vectors') # Plot eigenvectors
  # plot(vals_ssa, type='series') # Plot reconstructed series
  
  # Dominant seasonality (where present) is captured in first eigenvectors
  EigenDF_Wide <- data.frame(vals_ssa$U) %>%
    rename_all(~gsub('X', 'E', .x)) %>% 
    mutate(time = 1:nrow(.)) 
  
  # ggplot(EigenDF_Wide, aes(x=E4, y=E5)) + geom_point(lwd=2)
  
  EigenDF_Long <- EigenDF_Wide %>% reshape2::melt(c("time"))
  
  # ggplot(EigenDF_Long %>% filter(variable %in% c("E1", "E2")), 
  #        aes(x=time, y=value, group=variable, col=variable)) + 
  #   geom_line(lwd=2)
  
  # Loop over first 10 eigenvector pairs to detect potential seasonality
  Esp <- lapply(1:10, function(e){
    
    print(paste0("Processing eigenvector pairs ", e, " and ", e+1))
    
    ex <- EigenDF_Wide %>% pull(!!paste0("E", e))
    ey <- EigenDF_Wide %>% pull(!!paste0("E", e+1))
    
    # ggplot(data = data.frame(ex = ex, ey = ey), aes(x=ex, y=ey)) + geom_point(lwd=2)
    
    # Convert Cartesian spiral to polar coordinates
    e_p <- cart2polar(ex, ey) %>% mutate(index = 1:nrow(.))
    
    # Count how many spirals are present based on # of jumps in theta
    n_spirals <- e_p %>%
      mutate(theta = abs(theta)) %>%
      # Calculate jumps
      mutate(JUMPS = c(NA, (diff(sign(diff(theta))<0)), NA )) %>% 
      filter(JUMPS!=0) %>% 
      # Make sure jump is sustained for at least a quarter
      mutate(TDIFF = c(diff(index), NA) ) %>%
      filter(TDIFF > (omega/2/0.125) ) %>% 
      nrow()
    
    return( data.frame(EigenPairs = e, N_Spirals = ceiling(n_spirals/2)) )
    
  }) %>%
    bind_rows() %>%
    filter(N_Spirals > 0)
  
  print(paste0(paste(sort(unique(Esp$N_Spirals)), collapse = ", "), " spirals found!"))
  
  # 6. Generate harmonic regression model ------------------
  
  if(max(Esp$N_Spirals)>1){
    m_pref <- glm(detrend ~ SIN2PI + COS2PI + SIN4PI + COS4PI + 
                    INDEX + INDEX2 + INDEX3, 
                  data = timedf)
  } else if(max(Esp$N_Spirals)==1){
    m_pref <- glm(detrend ~ SIN2PI + COS2PI + INDEX + INDEX2 + INDEX3, 
                  data = timedf)
  }  
  
  # 7. Calculate peak timings --------------------

  if( exists("m_pref") ){
    
    PT <- peaktimecalc(m_pref, omega)
    
    if(isTRUE(transform_ln)){
      PT[[1]] <- PT[[1]] %>% 
        mutate(PEAKVALUE = exp(PEAKVALUE), NADIRVALUE = exp(NADIRVALUE),
               PEAK2VALUE = exp(PEAK2VALUE), NADIR2VALUE = exp(NADIR2VALUE))
    }
    
    # ggplot(PT[[2]], aes(x=INDEX, y=exp(PRED))) + geom_line()
    
    return(PT)
    
  }

  # End function
}