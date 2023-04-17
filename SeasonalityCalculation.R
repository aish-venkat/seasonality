# Function to calculate harmonic peaks and nadirs from weekly, monthly, or daily time series
# Last Updated: April 12, 2023

# Use pacman package to load/install needed packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, reshape2, multitaper, broom, Rssa, sf, imputeTS)

######################################################
# REFERENCES #########################################
######################################################

# Naumova, E. N. & MacNeill, I. B. (2007). 
# Seasonality assessment for biosurveillance systems. 
# In Advances in Statistical Methods for the Health Sciences (eds Mesbah, M. et al.) 
# 437–450 (Birkhäuser, Basel, 2007).

# Naumova, E. N., Jagai, J. S., Matyas, B., DeMaria, A., 
# MacNeill, I. B., & Griffiths, J. K. (2007). 
# Seasonality in six enterically transmitted diseases and 
# ambient temperature. Epidemiology & Infection, 135(2), 281-292.

# Alsova, O. K., Loktev, V. B., & Naumova, E. N. (2019). 
# Rotavirus seasonality: An application of singular spectrum analysis and 
# polyharmonic modeling. IJERPH 16(22), 4309.

# Falconi, T. A., Cruz, M. S. & Naumova, E. N. (2018).
# The shift in seasonality of legionellosis in the USA. 
# Epidemiol. Infect. 146, 1824–1833.

# Simpson, R. B., Zhou, B., & Naumova, E. N. (2020). Seasonal synchronization of 
# foodborne outbreaks in the United States, 1996–2017. Scientific reports, 10(1), 17500.
# https://www.nature.com/articles/s41598-020-74435-9
# Equations for negative binomial regression in Supplementary Table 3

# Ramanathan, K., Thenmozhi, M., George, S., Anandan, S., 
# Veeraraghavan, B., Naumova, E. N., & Jeyaseelan, L. (2020). 
# Assessing seasonality variation with harmonic regression: 
# accommodations for sharp peaks. International journal of environmental research and public health, 17(4), 1318.

######################################################
# SUPPORTING FUNCTIONS ###############################
######################################################

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
peaktimecalc <- function(mod, f, omega, vals_con, timedf, linkpar){
  
  # print(str(linkpar$family))
  
  # New data for which to predict
  preddf <- data.frame(INDEX=seq(1, f+1, by=f/100))
  preddf <- preddf %>%
    mutate(SIN2PI = sin(2*pi*omega*INDEX),
           COS2PI = cos(2*pi*omega*INDEX),
           SIN4PI = sin(4*pi*omega*INDEX),
           COS4PI = cos(4*pi*omega*INDEX)) %>%
    mutate(INDEX2 = INDEX^2, INDEX3 = INDEX^3) 

  # Remove linear, quadratic, cubic trend from model
  mod_notrend <- update(mod, as.formula(paste(c(formula(mod), " - INDEX2 - INDEX3"),
                                              collapse = " ") ) )
  # Add prediction from no trend model
  preddf <- preddf %>%
    mutate(PRED = as.vector(predict(mod_notrend, newdata=preddf, type="response"))) %>%
    # Adjust index of prediction to match omega for polar plots
    mutate(INDEX = INDEX-1)
  
  # Add original data back into prediction dataframe
  originaldf <- mod$model %>% data.frame() %>% 
    mutate(INDEX = (INDEX %% f),
           INDEX = ifelse(INDEX==0, f, INDEX)) %>% 
    group_by(INDEX) %>% summarize(ORIGINAL = mean(value, na.rm=T)) %>%
    ungroup()
  
  preddf <- preddf %>% left_join(originaldf, by="INDEX")

  #ggplot(preddf, aes(x=INDEX, y=PRED)) + geom_line() + scale_x_continuous(breaks=1:12)
  
  if(any(grepl("SIN4PI", mod$call))){
    
    preddf <- preddf %>%
      # Calculate derivatives
      mutate(dydx_1 = c(NA, diff(PRED)), dydx_2 = c(NA, NA, diff(diff(PRED)))) %>% 
      # Calculate peak values
      mutate(PEAKS = c(NA, (diff(sign(diff(PRED))<0)), NA )) %>% 
      mutate(MAXIMA = ifelse(PEAKS==1, 1, 0), MINIMA = ifelse(PEAKS==-1, 1, 0)) 
    
    peaktiming <- preddf %>% filter(MAXIMA==1) %>% arrange(-PRED) %>% pull(INDEX)
    peakvalue <- preddf %>% filter(MAXIMA==1) %>% arrange(-PRED) %>% pull(PRED)
    
    nadirtiming <- preddf %>% filter(MINIMA==1) %>% arrange(-PRED) %>% pull(INDEX)
    nadirvalue <- preddf %>% filter(MINIMA==1) %>% arrange(-PRED) %>% pull(PRED)
    
    # Calculate amplitudes from simulated time series ----------------
    fitted_model <- auto.arima(vals_con)
    
    # Simulate the time series and extract amplitudes several times
    SIMS <- lapply(1:99, function(z){
      
      #print(paste0("Running simulation: ", z))
      
      # Simulate one full cycle
      vals_sim <- simulate(fitted_model, n.sim = 3*length(vals_con))
      N <- length(vals_sim)
      
      fft_sim <- data.frame(coef = fft(vals_sim / N, inverse=T), 
                            freqindex = 1:length(vals_sim),
                            freq = (1:N)*(f)/(2*N)) %>% 
        mutate(value_original = vals_sim,
               amp = Mod(coef), freqindex = 1:N,
               value = Re(coef), time = freqindex %% f) %>% 
        slice(-1) %>% 
        mutate(SIM = z) %>% 
        # Calculate derivatives
        mutate(dydx_1 = c(NA, diff(amp)), dydx_2 = c(NA, NA, diff(diff(amp)))) %>% 
        # Calculate peak values
        mutate(PEAKS = c(NA, (diff(sign(diff(amp))<0)), NA )) %>% 
        mutate(MAXIMA = ifelse(PEAKS==1, 1, 0), MINIMA = ifelse(PEAKS==-1, 1, 0)) 
      
      return(fft_sim)
      
    }) %>% bind_rows()
    
    AMP <- SIMS %>% group_by(SIM) %>% 
      summarize(MEAN = mean(amp, na.rm=T),
                RANGE = abs(max(amp) - min(amp)),
                SD = sd(amp, na.rm=T))
    
    pt <- data.frame(MODEL = "4PI", 
                     PEAKTIMING = peaktiming[1], PEAKVALUE = peakvalue[1],
                     NADIRTIMING = nadirtiming[1], NADIRVALUE = nadirvalue[1],
                     PEAK2TIMING = peaktiming[2], PEAK2VALUE = peakvalue[2],
                     NADIR2TIMING = nadirtiming[2], NADIR2VALUE = nadirvalue[2],
                     AMPLITUDE = mean(AMP$MEAN), AMP_CI = mean(AMP$SD))
    
  } else {
    
    # Extract regression estimate
    coefs <- tidy(mod)
    
    intercept <- coefs %>% filter(term=="(Intercept)") %>% pull(estimate)
    
    coef_sin <- coefs %>% filter(term=="SIN2PI") %>% pull(estimate)
    coef_cos <- coefs %>% filter(term=="COS2PI") %>% pull(estimate)
    
    # Extract terms from variance-covariance matrix
    # Own-term variance is on diagonal, other terms are covariances
    X <- unname(vcov(mod))
    sigma_sin <- X[2,2]; sigma_cos <- X[3,3]; sigma_sincos <- X[2,3];

    # Peak timing (same across specifications) -----------
    
    # Define cycle
    period <- f/(2*pi);
    
    # Angle
    ang <- coef_sin / coef_cos
    
    # Phi or Shift 
    #shift <- atan2(coef_cos, coef_sin) # range in -pi to pi values
    shift <- atan(ang)
    
    # Variance of phi
    var_phi <- ( (coef_sin^2 * sigma_cos) + (coef_cos^2 * sigma_sin) - 
                   (2* sigma_sincos * coef_sin * coef_cos) )    / ( (coef_sin^2 + coef_cos^2)^2 )
    
    # Peak Timing
    peaktiming <- ifelse(coef_sin>0 & coef_cos>0, (shift*period),
                         ifelse(coef_sin>0 & coef_cos<0, (shift+pi)*period,
                                ifelse(coef_sin<0 & coef_cos<0, (shift+pi)*period,
                                       ifelse(coef_sin<0 & coef_cos>0, (shift + 2*pi)*period))))

    # CI of peak timing
    ci_peaktiming <- 1.96 * sqrt(var_phi) * period
    
    # Nadir Timing
    nadirtiming <- ifelse(coef_cos < 0, shift + (pi/2),
                          ifelse(coef_cos > 0 & coef_sin > 0, (shift/2), (shift + (2 * pi))/2 ))
    
    # Values different by model specification ------------
    if( any(linkpar$family %in% c("poisson", "quasipoisson")) ) {
      
      # Amplitude
      amp <- exp(sqrt(coef_sin^2 + coef_cos^2))
      
      # CI of amplitude
      var_amp <- ((coef_cos^2 * sigma_cos) + (coef_sin^2 * sigma_sin) + 
                    (2* sigma_sincos * coef_sin * coef_cos)) / (coef_sin^2 + coef_cos^2)
      var_amp <- (amp^2) * var_amp
      
      ci_amp <- 1.96 * sqrt(var_amp)
      
      # Peak Value
      peakvalue = exp(intercept) + amp
      
      # Nadir Value
      nadirvalue = exp(intercept) - amp
      
      # Absolute Intensity 
      intensity = peakvalue - nadirvalue
      var_intensity = 4*(var_amp)
      ci_intensity = 1.96 * sqrt((var_intensity))
      
      # Variance of relative intensity
      var_intensity_rel <- ( (exp(2*intercept) / (exp(intercept) - amp)^4) ) * var_intensity
      
    } else {
      
      # Amplitude
      amp <- sqrt(coef_sin^2 + coef_cos^2)
      
      # CI of amplitude
      var_amp <- ((coef_cos^2 * sigma_cos) + (coef_sin^2 * sigma_sin) + 
                    (2* sigma_sincos * coef_sin * coef_cos)) / (coef_sin^2 + coef_cos^2)
      ci_amp <- 1.96 * sqrt(var_amp)
      
      # Peak Value
      peakvalue = intercept + amp
      
      # Nadir Value
      nadirvalue = intercept - amp
      
      # Absolute Intensity 
      intensity = peakvalue - nadirvalue
      var_intensity = 4*(var_amp)
      ci_intensity = 1.96 * sqrt((var_intensity))
      
      # Variance of relative intensity
      var_intensity_rel <- ( (intercept^2) / (intercept - amp)^4 )  * var_intensity
      
    }

    # Relative Intensity 
    intensity_rel = peakvalue / nadirvalue
    var_intensity_rel = 4*(var_amp)
    ci_intensity_rel = 1.96 * sqrt((var_intensity_rel))
    
    pt <- data.frame(MODEL = "2PI", 
                      PEAKTIMING = peaktiming, PEAKVALUE = peakvalue,
                      NADIRTIMING = nadirtiming, NADIRVALUE = nadirvalue,
                      PEAK2TIMING = NA, PEAK2VALUE = NA,
                      NADIR2TIMING = NA, NADIR2VALUE = NA,
                      AMPLITUDE = amp, AMP_CI = ci_amp, 
                      INTENSITY = intensity, INTENSITY_CI = ci_intensity,
                     INTENSITY_REL = intensity_rel, INTENSITY_REL_CI = ci_intensity_rel)
    
  }
  
  # Return estimated peak timings and predicted values
  return(list(pt, preddf))
  
}

# Function to calculate seasonality 
seasonalitycalc <- function(df, tfield, f, outcome, 
                            linkpar = gaussian(link="identity"),
                            fspec = FALSE, fsing = FALSE){
  
  ##########################
  # EXPLANATION OF ARGUMENTS
  
  # df: data frame
  # tfield: string title of column containing time field; can be date or numeric
  # outcome: string title of column containing outcome of interest
  # linkpar : FAMILY AND LINK FUNCTION parameters for expected values and linear predictors
  #           specify as you would for GLM function (e.g. binomial(link = 'logit'))
  # fspec : whether spectral analysis should be implemented; results will be returned in list element [[3]]
  # fsing : whether singular spectral analysis (eigenvector analysis) should be implemented;
  #            results will be returned in list element [[4]]
  
  ##########################

  # Process time fields
  if(f==12){
    fstr = "months"
  } else if(f==365 | f==365.25){
    fstr = "days"
  } else if(f==52 | f==53){
    fstr = "weeks"
  }
  
  timefield <- df %>% pull(tfield)
  
  if(class(timefield)=="Date"){
    
    # Round to beginning of year
    mindate <- min(timefield, na.rm=T); mindate <- floor_date(mindate, unit='year');
    maxdate <- max(timefield, na.rm=T); maxdate <- ceiling_date(maxdate, unit='year');
    
    # Create sequential series
    omega <- 1/f
    indexdf <- data.frame(DATE = seq.Date(mindate, maxdate, by=fstr)) %>%
      mutate(INDEX = 1:nrow(.), INDEX2 = INDEX^2, INDEX3 = INDEX^3,
             SIN2PI = sin(2*pi*omega*INDEX), 
             COS2PI = cos(2*pi*omega*INDEX),
             SIN4PI = sin(4*pi*omega*INDEX), 
             COS4PI = cos(4*pi*omega*INDEX)) %>%
      mutate(YEAR = year(DATE))
    
    # Merge sequence into time series
    timedf <- data.frame(DATE = timefield, value = df %>% pull(outcome)) %>% 
      full_join(indexdf, by="DATE") %>% 
      # Simple summary of outcome by date for time series
      group_by(DATE, YEAR, INDEX, INDEX2, INDEX3, SIN2PI, COS2PI, SIN4PI, COS4PI) %>% 
      summarize(value = mean(value, na.rm=T)) %>%
      ungroup() %>%
      arrange(DATE) %>%
      distinct()
    
  } else if(class(timefield)=="Numeric"){
    
    # Create sequential series
    
    mindate <- min(timefield, na.rm=T); maxdate <- max(timefield, na.rm=T)
    indexdf <- data.frame(DATE = seq(mindate, maxdate, by=f),
                          value = df %>% pull(outcome)) %>%
      mutate(INDEX = 1:nrow(.)) %>%
      mutate(INDEX = 1:nrow(.), INDEX2 = INDEX^2, INDEX3 = INDEX^3,
             SIN2PI = sin(2*pi*omega*INDEX), 
             COS2PI = cos(2*pi*omega*INDEX),
             SIN4PI = sin(4*pi*omega*INDEX), 
             COS4PI = cos(4*pi*omega*INDEX)) %>%
      mutate(YEAR = year(DATE))
    
    # Merge sequence into time series
    timedf <- data.frame(DATE = timefield, value = df %>% pull(outcome)) %>% 
      full_join(indexdf, by="DATE") %>% 
      arrange(DATE) %>% distinct()
    
  }
  
  ######################################################
  # PERIODICITY ASSESSMENT #############################
  ######################################################
  
  # 1. Create time series for complete cycles -----------
  vals <- ts(timedf$value,
             start = c(year(timedf$DATE[1]), month(timedf$DATE[1])),
             deltat = 1/f)
  
  N <- sum(!is.na(timedf$value))
  
  # 2. If time series is sufficiently long, 
  # 2. Subset time series for complete cycles ----------
  # 2. If not, fill in missing values using kalman filter / interpolation
  year_totals <- timedf %>%
    drop_na(value) %>%
    group_by(YEAR) %>% tally() %>% ungroup() %>%
    # Subset only complete observations
    filter(n==f) %>% pull(YEAR)

  if(length(year_totals)>2){
    vals_con <- ts(timedf %>% filter(YEAR %in% year_totals) %>% pull(value),
                   start = c(year(timedf$DATE[1]), month(timedf$DATE[1])),
                   deltat = 1/f)
  } else if(isTRUE(fspec) | isTRUE(fsing)) {
    # Only use Kalman smoothing on time series when spectral analyses are needed
    vals_con <- na_seadec(vals, algorithm="kalman")
  } else if(N < f/2){
    # If less than half of a cycle is available, ignore 
    vals_con <- NA
  } else{
    vals_con <- vals
  }

  # 3. Spectral Analysis ----------------------
  
  if(fspec){
    
    # Use complete cycles to look at spectra
    s <- spec.mtm(vals_con, Ftest = TRUE, jackknife = T)
    
    # calculate spectrum
    spec1 <- s$spec[1]
    s.df <- data.frame(freq = s$freq,  spec = s$spec, spec_scaled = s$spec/spec1,
                       f = s$mtm$Ftest, var_jk = s$mtm$jk$varjk,
                       ci_lower = s$mtm$jk$lowerCI, ci_upper = s$mtm$jk$upperCI) %>%
      # Multiply spectra by 2
      #mutate(spec = 2*spec) %>%
      # Remove zero value as it is meaningless
      dplyr::filter(freq!=0) %>%
      # Convert from frequency to periods
      mutate(period = freq / (1/f),
             time = round(period - (f * (period %/% f)), 0) )
    
    # Which frequencies are statistically significant?
    sig_months <- s.df %>% dplyr::filter(f <= 0.05)
    
    # Are same frequencies statistically significant across cycles?
    # ggplot(data = s.df, aes(x=period, y=spec)) + geom_line() +
    #   geom_vline(xintercept = s.df %>% filter(f <= 0.05) %>% pull(period), lty=2, col='red') +
    #   scale_x_continuous("Period (years)", breaks = seq(0, length(vals)*2, by=f)) +
    #   scale_y_log10()
    # yrs.period <- rev(c(1/12, 1/10, 1/8, 1/6, 1/5, 1/4, 1/3, 0.5, 1, 2))
    # yrs.labels <- rev(c("1/12", "1/10", "1/8", "1/6", "1/5", "1/4", "1/3", "1/2", "1", "2"))
    # yrs.freqs <- (1/yrs.period) * (1/f)  # Convert annual period --> annual freq --> monthly freq
    # ggplot(data = s.df %>% filter(freq>1/12), aes(x=freq, y=spec)) + geom_line() +
    #   geom_vline(xintercept = s.df %>% filter(f <= 0.05) %>% pull(freq), lty=2, col='red') +
    #   scale_x_log10("Period (years)", breaks = yrs.freqs, labels = yrs.labels) +
    #   scale_y_log10()
    
    
  } else{
    s.df <- NA
  }

  # 4. Singular Spectrum Analysis -------------------
  
  if(fsing){
    
    vals_ssa <- ssa(vals_con)
    
    # plot(vals_ssa, type='vectors') # Plot eigenvectors
    # plot(vals_ssa, type='series') # Plot reconstructed series
    
    # Dominant seasonality (where present) is captured in first eigenvector pairs
    EigenDF <- data.frame(vals_ssa$U) %>%
      rename_all(~gsub('X', 'E', .x)) %>%
      mutate(time = 1:nrow(.))
    
    # ggplot(EigenDF, aes(x=E4, y=E5)) + geom_point(lwd=2)
    
    # Loop over first 10 eigenvector pairs to detect potential seasonality
    num_eigens <- grep("E", names(EigenDF))
    Esp <- lapply(2:(length(num_eigens)-1), function(e){
      
      print(paste0("Processing eigenvector pairs ", e, " and ", e+1))
      
      ex <- EigenDF %>% pull(!!paste0("E", e))
      ey <- EigenDF %>% pull(!!paste0("E", e+1))
      
      p <- ggplot(data = data.frame(ex = ex, ey = ey), aes(x=ex, y=ey)) + geom_point(lwd=2)
      print(p)
      
      # Convert Cartesian spiral to polar coordinates
      e_p <- cart2polar(ex, ey) %>% mutate(index = 1:nrow(.)) %>%
        mutate(theta_abs = abs(theta)) %>%
        # Calculate jumps from derivatives
        mutate(JUMPS = c(NA, (diff(sign(diff(theta_abs))<0)), NA ))
      
      # Count how many spirals are present based on # of jumps in theta
      n_spirals <- e_p %>%
        filter(JUMPS!=0) %>%
        # Make sure jump is sustained for at least half cycle
        mutate(TDIFF = c(diff(index), NA) ) %>%
        filter(TDIFF > f/2 ) %>%
        nrow()
      
      # Assign shape index to distinguish shapes into groups
      coords <- data.frame(ex, ey, e_p) %>%
        mutate(JUMPS = ifelse(is.na(JUMPS)|JUMPS<0, 0, JUMPS)) %>%
        mutate(GID = cumsum(JUMPS) + JUMPS*0)
      
      # Calculate approximate centroid and area of each spiral
      if(n_spirals > 0){
        
        # Create shape
        sf::sf_use_s2(FALSE)
        spiralshape <- st_as_sf(coords, coords = c("ex", "ey"), crs = 4326) %>%
          group_by(GID) %>%
          tally() %>%
          summarise(geometry = st_combine(geometry)) %>%
          st_cast("POLYGON") %>%
          ungroup() %>%
          st_make_valid() %>%
          mutate(A = st_area(.), C = st_centroid(.)) %>%
          # Look only at shapes defined by # of spirals
          arrange(-A) %>% slice(1:n_spirals) %>%
          # Convert from polar to cartesian coordinates
          mutate(C_rad = st_coordinates(C)[,1],
                 C_theta = st_coordinates(C)[,2]) %>%
          mutate(C_x = polar2cart(C_rad, C_theta)[,1],
                 C_y = polar2cart(C_rad, C_theta)[,2])
        
        # Phase shift in time units
        if(transformation == "log"){
          phaseshift <- diff(exp(spiralshape %>% st_drop_geometry() %>% pull(C_x)))[1]
        } else{
          phaseshift <- diff(spiralshape %>% st_drop_geometry() %>% pull(C_x))[1]
        }
        
      } else{
        phaseshift <- NA
      }
      
      return( data.frame(EigenPairs = e,
                         N_Spirals = ceiling(n_spirals/2),
                         Phase = phaseshift) )
      
    })
    E_n <- Esp %>% bind_rows() %>% filter(N_Spirals > 0)
    
    print(paste0(paste(sort(unique(E_n$N_Spirals)), collapse = ", "), " spirals found!"))
    
  } else{
    E_n <- NA
  }

  ######################################################
  # MODEL FIT ASSESSMENT #############################
  ######################################################
  
  # Calculate models if non-NA values total at least half a cycle
  if( N > (f/2) ){
    
    # Remove any infinite values from timedf
    timedf <- timedf %>% filter(!is.na(value) & !is.infinite(value))
    
    m1 <- glm(value ~ SIN2PI + COS2PI + INDEX + INDEX2 + INDEX3, 
              data = timedf, family = linkpar )
    m2 <- glm(value ~ SIN2PI + COS2PI + SIN4PI + COS4PI + 
                INDEX + INDEX2 + INDEX3, data = timedf, family = linkpar )
    
    # Are 4pi terms significant?
    psig_4 <- tidy(m2) %>% filter(grepl("4PI", term)) %>% filter(p.value<0.05) %>% nrow()
    psig_2 <- tidy(m1) %>% filter(grepl("2PI", term)) %>% filter(p.value<0.05) %>% nrow()
    
    # Compare models
    m_pref <- glance(m1) %>% mutate(MODEL = "2PI") %>% 
      bind_rows(glance(m2) %>% mutate(MODEL = "4PI")) %>% 
      melt(c("MODEL")) %>%
      filter(variable=="BIC") %>% 
      slice(which.min(value)) %>%
      pull(MODEL)
    
    if(psig_4>0){
      m_pref <- m2
    } else if(psig_2>0){
      m_pref <- m1
    } else{
      m_pref <- NA
    } 
    
    if( !all(is.na(m_pref)) ){
      
      PT <- peaktimecalc(m_pref, f, omega, vals_con, timedf, linkpar)
      
    } else{
      PT <- NA
    } # End m_pref check
    
  }  else{
    PT <- NA
  } # End vals_con length check

  ######################################################
  # HARMONIC CHARACTERISTICS CALCULATION ###############
  ######################################################
  
  # Add s.df and E_n as list elements
  PT <- append(PT, s.df); PT <- append(PT, E_n)
  
  return( PT )
}
