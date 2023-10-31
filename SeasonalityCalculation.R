# Function to calculate harmonic peaks and nadirs from weekly, monthly, or daily time series
# Last Updated: April 23, 2023

# Use pacman package to load/install needed packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, reshape2, multitaper, broom, Rssa, sf, imputeTS, boot)

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

# Function to extract measures for multiple harmonic regression from predicted curves
calc_multiple_harmonics <- function(prediction_frame){
  
  prediction_frame <- prediction_frame %>%
    # Calculate derivatives
    mutate(dydx_1 = c(NA, diff(PRED)), 
           dydx_2 = c(NA, NA, diff(diff(PRED))),
           dydx_3 = c(NA, NA, NA, diff(diff(diff(PRED)))),
           dydx_4 = c(NA, NA, NA, NA, diff(diff(diff(diff(PRED))))) ) %>% 
    # Calculate signs of differences
    mutate(SIGDIFF_1 = c(NA, (diff(sign(PRED)<0))),
           SIGDIFF_2 = c(NA, (diff(sign(dydx_1)<0))),
           SIGDIFF_3 = c(NA, (diff(sign(dydx_2)<0))),
           SIGDIFF_4 = c(NA, (diff(sign(dydx_3)<0))))  %>%
    # Calculate peak values 
    mutate(PEAKS = SIGDIFF_2 ) %>% 
    mutate(MAXIMA = ifelse(PEAKS==1, 1, 0), 
           MINIMA = ifelse(PEAKS==-1, 1, 0))
  
  # Visual check: confirm alignment of 2nd/3rd derivatives
  # vizdf <- prediction_frame %>% mutate(PRED = rescale(PRED, to=c(0,1)),
  #                            dydx_1 = rescale(dydx_1, to=c(0,1)),
  #                            dydx_2 = rescale(dydx_2, to=c(0,1)),
  #                            dydx_3 = rescale(dydx_3, to=c(0,1)),
  #                            dydx_4 = rescale(dydx_4, to=c(0,1)))
  # ggplot(vizdf, aes(x=INDEX)) +
  #   geom_line(aes(y=PRED), col='black', lwd=2) +
  #   # First deriv (no sign diff in PRED)
  #   geom_line(aes(y=dydx_1), col='red', lty=3, lwd=2) +
  #   geom_point(data = vizdf %>% filter(SIGDIFF_1!=0), aes(x=INDEX, y=dydx_1),
  #              inherit.aes = F, col='red', size=8, pch=3) +
  #   # Second deriv
  #   geom_line(aes(y=dydx_2), col='green', lty=3, lwd=2) +
  #   geom_point(data = vizdf %>% filter(SIGDIFF_2!=0), aes(x=INDEX, y=dydx_2),
  #              inherit.aes = F, col='green', size=8, pch=3) +
  #   # Third deriv
  #   geom_line(aes(y=dydx_3), col='blue', lty=3, lwd=2) +
  #   geom_point(data = vizdf %>% filter(SIGDIFF_3!=0), aes(x=INDEX, y=dydx_3),
  #              inherit.aes = F, col='blue', size=8, pch=3) +
  #   # Fourth deriv
  #   geom_line(aes(y=dydx_4), col='purple', lty=3, lwd=2) +
  #   geom_point(data = vizdf %>% filter(SIGDIFF_3!=0), aes(x=INDEX, y=dydx_3),
  #              inherit.aes = F, col='purple', size=8, pch=3)

  # Amplitude (empirical) --------
  
  # Get list of breaks from 2nd and 3rd derivatives
  DERIV_BREAKS <- prediction_frame %>% 
    dplyr::select(INDEX, PRED, contains("SIGDIFF"), MAXIMA, MINIMA) %>%
    mutate(across(where(is.numeric), ~na_if(., 0))) %>% 
    mutate_at(vars(starts_with("SIG")), as.character) %>%
    mutate(SIGDIFF_2 = gsub("1|-1", "2nd", SIGDIFF_2),
           SIGDIFF_3 = gsub("1|-1", "3rd", SIGDIFF_3)) %>%
    mutate(SIG_BRK = coalesce(SIGDIFF_2, SIGDIFF_3)) %>% 
    filter(!is.na(SIG_BRK)) %>%
    dplyr::select(INDEX, PRED, SIG_BRK, MAXIMA, MINIMA) %>% 
    arrange(INDEX)

  # If no inflection points are found, break out of loop
  if( nrow(DERIV_BREAKS)==0){
    return(NULL) 
  } else{
    # Repeat first index at the end to support looping over end of calendar year
    DERIV_BREAKS <- DERIV_BREAKS %>% bind_rows(DERIV_BREAKS %>% slice(1))
    
    DERIV_BREAKS <- lapply(1:(nrow(DERIV_BREAKS)-1), function(x){
      data.frame(from = DERIV_BREAKS$INDEX[x], 
                 to = DERIV_BREAKS$INDEX[x+1],
                 fromtype = DERIV_BREAKS$SIG_BRK[x],
                 totype = DERIV_BREAKS$SIG_BRK[x+1], 
                 # Retain fields from previous
                 PRED = DERIV_BREAKS$PRED[x+1],
                 MAXIMA = DERIV_BREAKS$MAXIMA[x+1],
                 MINIMA = DERIV_BREAKS$MINIMA[x+1])
    }) %>% bind_rows()
    
    # Calculate amplitudes for each interval
    AMPS <- lapply(1:nrow(DERIV_BREAKS), function(x){
      
      if( DERIV_BREAKS$from[x] < DERIV_BREAKS$to[x] ){
        derivchunk <- prediction_frame %>%
          filter(INDEX >= DERIV_BREAKS$from[x] & INDEX <=  DERIV_BREAKS$to[x]) 
      } else{
        derivchunk <- prediction_frame %>% 
          filter(INDEX >= DERIV_BREAKS$from[x] | INDEX <=  DERIV_BREAKS$to[x])
      }

      # Get slope of dydx_3 and dydx_4
      if(nrow(derivchunk)>3){
        SL3 <- sign(coef(lm(dydx_3 ~ INDEX, data = derivchunk))[[2]])
        SL4 <- sign(coef(lm(dydx_4 ~ INDEX, data = derivchunk))[[2]])
      } else{
        SL3 <- NA; SL4 <- NA;
      }
      
      derivchunk %>% summarize(AMP = max(PRED) - min(PRED)) %>% 
        bind_cols(DERIV_BREAKS %>% slice(x))  %>%
        mutate(SLOPE_3 = SL3, SLOPE_4 = SL4) %>%
        arrange(from, to, AMP)
      
    }) %>% bind_rows() 
  
    # PULL breakpoints where 2nd and 3rd derivatives have the same sign / 
    # are moving in the same direction
    AMPS <- AMPS %>%  filter(SLOPE_3 == SLOPE_4) 
    
    # IF only two rows are present with the same amplitudes, we have a square curve
    # so the 'PEAK' is technically in between the two zero values
    if(nrow(AMPS)>0){
      
      if( all( as.character(AMPS$AMP) %in% c("1", "0") ) ){
        print(AMPS) %>% data.frame()
        AMPS <- AMPS %>% 
          mutate(MINIMA = ifelse(AMP==0, 1, NA), MAXIMA = ifelse(AMP==1, 1, NA)) %>% 
          rowwise() %>% mutate(NEWMEAN = mean(from, to, na.rm=T)) %>% ungroup() %>% 
          mutate(from = NEWMEAN, to = NEWMEAN) %>% 
          dplyr::select(-NEWMEAN)
      }
      
      # Extract values, but
      peaktiming <- AMPS %>% filter(MAXIMA==1) %>% arrange(-AMP) %>% pull(to)
      peakvalue <- AMPS %>% filter(MAXIMA==1) %>% arrange(-AMP) %>% pull(PRED)
    
      nadirtiming <- AMPS %>% filter(MINIMA==1) %>% arrange(-AMP) %>% pull(to)
      nadirvalue <- AMPS %>% filter(MINIMA==1) %>% arrange(-AMP) %>% pull(PRED)
      
    
      return ( data.frame(MODEL = "4PI", 
                          PEAKTIMING = peaktiming[1], PEAKVALUE = peakvalue[1],
                          NADIRTIMING = nadirtiming[1], NADIRVALUE = nadirvalue[1],
                          PEAK2TIMING = peaktiming[2], PEAK2VALUE = peakvalue[2],
                          NADIR2TIMING = nadirtiming[2], NADIR2VALUE = nadirvalue[2],
                        AMP1 = AMPS$AMP[1], AMP2 = AMPS$AMP[2] ) )
  
    } # End nrow check for AMPS
  } # End check for zero inflection points
} # End calc_multiple_harmonics function 

# Function to calculate peak/nadir timings and values 
peaktimecalc <- function(mod, f, omega, vals_con, timedf, linkpar){
  
  # print(str(linkpar$family))
  
  # Calculate predicted values from preferred model -----------------
  
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
    mutate(PRED = as.vector(predict(mod_notrend, newdata=preddf, type='response'))) %>%
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
    
    # Calculate actual peak/nadir values from model prediction ----------
    PRED_M <- calc_multiple_harmonics(preddf)

    # Run subsequent lines if PRED_M is not NULL (i.e. enough observations to calculate peaks)
    
    if( !is.null(PRED_M)  ){

        # Calculate amplitudes and CIs from simulated time series ----------------
        fitted_model <- auto.arima(vals_con)
        
        # Simulate the time series 999 times and retain those highly correlated with original series
        SIMS <- lapply(1:999, function(z){
          
          #print(paste0("Running simulation: ", z))
          vals_sim <- simulate(fitted_model, nsim=length(vals_con), bootstrap=T)
          N <- length(vals_sim)
          
          # Only run downstream blocks if cross-correlation between
          # FIRST-DIFFERENCED original and simulated time series at lag 0 is statsig
          
          ccftest <- ccf(diff(vals_con %>% as.numeric()), # 1st diff original TS
                         diff(vals_sim %>% as.numeric()), # 1st diff simulated TS
                         lag.max = f, type='correlation', plot=F)
          ss_ccf <- data.frame(LAG = ccftest$lag, ACF = ccftest$acf) %>%
            filter(LAG == 0) %>%
            mutate(THRESHOLD_UPPER = qnorm((1 + 0.95)/2)/sqrt(length(vals_con)),
                   THRESHOLD_LOWER = -THRESHOLD_UPPER) %>%
            filter(ACF>=THRESHOLD_UPPER | ACF<=THRESHOLD_LOWER)
          
          if( nrow(ss_ccf) > 0 ){
            
            # Create new simulation data frame
            simdf <- data.frame(INDEX=1:N, value=vals_sim %>% as.vector() )
            simdf <- simdf %>%
              mutate(SIN2PI = sin(2*pi*omega*INDEX),
                     COS2PI = cos(2*pi*omega*INDEX),
                     SIN4PI = sin(4*pi*omega*INDEX),
                     COS4PI = cos(4*pi*omega*INDEX)) %>%
              mutate(INDEX2 = INDEX^2, INDEX3 = INDEX^3) 
            
            # Wrap simulation regression in tryCatch in case
            # there are errors with the simulated model
            MODELSIM_GLM <- function(simdf) {
              out <- tryCatch(
                {
                  glm(value ~ SIN2PI + COS2PI + SIN4PI + COS4PI,
                      data = simdf, family = linkpar )
                },
                error=function(cond) { return(NA) },
                warning=function(cond) { return(NA) },
                finally=function(cond) { }
              )    
              return(out)
            }
            
            m_sim <- MODELSIM_GLM(simdf)
            
            if(! all(is.na(m_sim)) ){
              
              # Predict values for time series
              simpreddf <- data.frame(INDEX=seq(1, f+1, by=f/100))
              simpreddf <- simpreddf %>%
                mutate(SIN2PI = sin(2*pi*omega*INDEX),
                       COS2PI = cos(2*pi*omega*INDEX),
                       SIN4PI = sin(4*pi*omega*INDEX),
                       COS4PI = cos(4*pi*omega*INDEX)) %>%
                mutate(INDEX2 = INDEX^2, INDEX3 = INDEX^3) 
              
              # Add prediction from no trend model
              simpreddf <- simpreddf %>%
                mutate(PRED = as.vector(predict(m_sim, newdata=simpreddf, type='response'))) %>%
                # Adjust index of prediction to match omega for polar plots
                mutate(INDEX = INDEX-1)
              
              mh <- calc_multiple_harmonics(simpreddf)
              print(mh); return(mh); 
              
            } else{
              return(NULL)
            }      
          } # End correlation check
          
        }) %>% bind_rows()
    
        # Run following if simulations yielded at least one simulated time series
        if(nrow(SIMS)>0){
      
            # Function to bootstrap CIs of variables
            btvar <- function(var){
              varvec <- SIMS %>% pull(get(!!var))
              if( (length(varvec)>0) & (!all(is.na(varvec))) ){
                varbt <- boot( na.omit(varvec) , function(u,i) mean(u[i]), R=999)
                varci <- boot.ci(varbt, conf=0.99, type='bca')
                return( mean(abs(varci$bca[4:5] - varci$t0)) )
              } else{
                return(NA)
              }
            }
        
            pt <- data.frame(MODEL = "4PI", 
                             PEAKTIMING = PRED_M$PEAKTIMING, PEAKTIMING_CI = btvar("PEAKTIMING"),
                             PEAKVALUE = PRED_M$PEAKVALUE, PEAKVALUE_CI = btvar("PEAKVALUE"),
                             NADIRTIMING = PRED_M$NADIRTIMING, NADIRTIMING_CI = btvar("NADIRTIMING"),
                             NADIRVALUE = PRED_M$NADIRVALUE, NADIRVALUE_CI = btvar("NADIRVALUE"),
                             PEAK2TIMING = PRED_M$PEAK2TIMING, PEAK2TIMING_CI = btvar("PEAK2TIMING"),
                             PEAK2VALUE = PRED_M$PEAK2VALUE, PEAK2VALUE_CI = btvar("PEAK2VALUE"),
                             NADIR2TIMING = PRED_M$NADIR2TIMING, NADIR2TIMING_CI = btvar("NADIR2TIMING"),
                             NADIR2VALUE = PRED_M$NADIR2VALUE, NADIR2VALUE_CI = btvar("NADIR2VALUE"),
                             AMP1 = PRED_M$AMP1, AMP1_CI = btvar("AMP1"),
                             AMP2 = PRED_M$AMP2, AMP2_CI = btvar("AMP2")) %>%
              rowwise() %>%
              mutate(INTENSITY =  (PEAKVALUE - NADIRVALUE) , 
                     INTENSITY2 =  (PEAK2VALUE - NADIR2VALUE),
                     INTENSITY_REL = (PEAKVALUE / NADIRVALUE),
                     INTENSITY2_REL = (PEAK2VALUE / NADIR2VALUE) ) %>%
              ungroup()
            # Not calculating CIs of intensity for 4pi specifications
    
        } else{
            pt <- data.frame(MODEL = "4PI", 
                             PEAKTIMING = PRED_M$PEAKTIMING, PEAKTIMING_CI = NA,
                             PEAKVALUE = PRED_M$PEAKVALUE, PEAKVALUE_CI = NA,
                             NADIRTIMING = PRED_M$NADIRTIMING, NADIRTIMING_CI = NA,
                             NADIRVALUE = PRED_M$NADIRVALUE, NADIRVALUE_CI = NA,
                             PEAK2TIMING = PRED_M$PEAK2TIMING, PEAK2TIMING_CI = NA,
                             PEAK2VALUE = PRED_M$PEAK2VALUE, PEAK2VALUE_CI = NA,
                             NADIR2TIMING = PRED_M$NADIR2TIMING, NADIR2TIMING_CI = NA,
                             NADIR2VALUE = PRED_M$NADIR2VALUE, NADIR2VALUE_CI = NA,
                             AMP1 = PRED_M$AMP1, AMP1_CI = NA,
                             AMP2 = PRED_M$AMP2, AMP2_CI = NA) %>%
              rowwise() %>%
              mutate(INTENSITY =  (PEAKVALUE - NADIRVALUE) , 
                     INTENSITY2 =  (PEAK2VALUE - NADIR2VALUE),
                     INTENSITY_REL = (PEAKVALUE / NADIRVALUE),
                     INTENSITY2_REL = (PEAK2VALUE / NADIR2VALUE) ) %>%
              ungroup()
        }
    } else{
      pt <- NA
    } # End PRED_M check
    
  } else{
    
    # Extract regression estimate
    coefs <- tidy(mod_notrend)
    
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
    ci_peaktiming <- 1.96 * sqrt(var_phi)
    
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
      
      # CI of peak value
      ci_peakvalue <- mean(amp + ci_amp, amp - ci_amp)
      
      # Nadir Value
      nadirvalue = exp(intercept) - amp
      
      # CI of nadir value
      ci_nadirvalue <- mean(amp + ci_amp, amp - ci_amp)
      
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
      
      # CI of peak value
      ci_peakvalue <- mean(amp + ci_amp, amp - ci_amp)
      
      # Nadir Value
      nadirvalue = intercept - amp
      
      # CI of nadir value
      ci_nadirvalue <- mean(amp + ci_amp, amp - ci_amp)
      
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
                     PEAKTIMING = peaktiming, PEAKTIMING_CI = ci_peaktiming,
                     PEAKVALUE = peakvalue, PEAKVALUE_CI = ci_peakvalue,
                     NADIRTIMING = nadirtiming, NADIRTIMING_CI = ci_peaktiming,
                     NADIRVALUE = nadirvalue,  NADIRVALUE_CI = ci_nadirvalue,
                     PEAK2TIMING = NA, PEAK2VALUE = NA,
                     NADIR2TIMING = NA, NADIR2VALUE = NA,
                     AMP1 = amp, AMP1_CI = ci_amp, AMP2 = NA, AMP2_CI = NA,
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
  
  print(df %>% slice(1))
  
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
  
  # If the model family is measuring continuous variables, they can be averaged
  if(! (linkpar$family %in% c("binomial", "poisson", "quasi")) ){
    timedf <- timedf %>%
      # Simple summary of outcome by date for time series
      group_by(DATE, YEAR, INDEX, INDEX2, INDEX3, SIN2PI, COS2PI, SIN4PI, COS4PI) %>% 
      summarize(value = mean(value, na.rm=T)) %>%
      ungroup()
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
  
  # Count number of unique months for which data is available
  month_totals <- length(which(!is.na(vals)) %% f %>% unique())
  
  # Decide how to create continuous subsets
  if(length(year_totals)>0){
    vals_con <- ts(timedf %>% filter(YEAR %in% year_totals) %>% pull(value),
                   start = c(year(timedf$DATE[1]), month(timedf$DATE[1])),
                   deltat = 1/f)
    
  } else if(length(year_totals)==0 & (isTRUE(fspec) | isTRUE(fsing)) ) {
    # Only use Kalman smoothing on time series when spectral analyses are needed
    vals_con <- na_seadec(vals, algorithm="kalman")
    
  } else if( month_totals >= (f/2) ){
    # If less than one cycle is available, generate filled cycle of averages based on moving average
    vals_avg <- tapply(vals, cycle(vals), mean, na.rm=T)
    vals_con <- ts(vals_avg %>% as.vector(),
                   start = c(year(timedf$DATE[1]), month(timedf$DATE[1])),
                   deltat = 1/f)
    vals_con <- na_ma(vals_con, k=2)
    
  } else{
    vals_con <- NA
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
  if( !all(is.na(vals_con)) ){
    
    # Remove any infinite values from timedf
    timedf <- timedf %>% filter(!is.na(value) & !is.infinite(value))
    
    m1 <- glm(value ~ SIN2PI + COS2PI + INDEX + INDEX2 + INDEX3, 
              data = timedf, family = linkpar )
    m2 <- glm(value ~ SIN2PI + COS2PI + SIN4PI + COS4PI + 
                INDEX + INDEX2 + INDEX3, data = timedf, family = linkpar )
    
    # Are 4pi terms significant?
    psig_4 <- tidy(m2) %>% filter(grepl("4PI", term)) %>% filter(p.value<0.05) %>% nrow()
    psig_2 <- tidy(m1) %>% filter(grepl("2PI", term)) %>% filter(p.value<0.05) %>% nrow()
    
    # BIC penalizes distance away from TRUE model; therefore, not appropriate for this task
    # We use AIC here instead
    # Compare models
    # m_pref <- glance(m1) %>% mutate(MODEL = "2PI") %>%
    #   bind_rows(glance(m2) %>% mutate(MODEL = "4PI")) %>%
    #   melt(c("MODEL")) %>%
    #   filter(variable=="AIC") %>%
    #   slice(which.min(value)) %>%
    #   pull(MODEL)
    
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
