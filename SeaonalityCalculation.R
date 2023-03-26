# Function to calculate harmonic peaks and nadirs from weekly, monthly, or daily time series
# Last Updated: March 25, 2023

library(tidyverse)

# Function to calculate peak timings from regression models 
peaktimecalc <- function(mod){
  
  if(any(grepl("SIN4PI", mod$call))){
    
    # Ref: https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
    maxpred<-preddf[(which(diff(sign(diff(preddf$PRED)))<0)+1),] %>% 
      arrange(PRED) %>% mutate_if(is.numeric, ~round(., 5)) %>%
      distinct(INDEX, SIN2PI, COS2PI, SIN4PI, COS4PI, .keep_all = TRUE)
    
    minpred<-preddf[(which(diff(sign(diff(preddf$PRED)))>0)+1),] %>% 
      arrange(PRED) %>% mutate_if(is.numeric, ~round(., 5)) %>%
      distinct(INDEX, SIN2PI, COS2PI, SIN4PI, COS4PI, .keep_all = TRUE)
    
    peaktiming<-maxpred %>% arrange(-PRED) %>% slice(1) %>% pull(INDEX)
    peakvalue <- maxpred %>% arrange(-PRED) %>% slice(1) %>% pull(PRED)
    
    nadirtiming<-minpred %>% arrange(PRED) %>% slice(1) %>% pull(INDEX)
    nadirvalue <- minpred %>% arrange(PRED) %>% slice(1) %>% pull(PRED)
    
    if(length(nadirtiming)==0){
      secpeaktiming <- NA; secpeakvalue <- NA;
      nadirtiming <- NA; nadirvalue <- NA;
      secnadirtiming <- NA; secnadirvalue <- NA;
    } else{
      if(nrow(maxpred)>1 & nrow(minpred)>1){
        secpeaktiming <- maxpred %>% arrange(-PRED) %>% slice(2) %>% pull(INDEX)
        secpeakvalue <- maxpred %>% arrange(-PRED) %>% slice(2) %>% pull(PRED)
        
        secnadirtiming <- minpred %>% arrange(PRED) %>% slice(2) %>% pull(INDEX) 
        secnadirvalue <- minpred %>% arrange(PRED) %>% slice(2) %>% pull(PRED)      
        
      } else if(nrow(maxpred)>1){
        secpeaktiming <- maxpred %>% arrange(-PRED) %>% slice(2) %>% pull(INDEX)
        secpeakvalue <- maxpred %>% arrange(-PRED) %>% slice(2) %>% pull(PRED)
        
        secnadirtiming <- NA; secnadirvalue <- NA;
      } else if(nrow(minpred)>1){
        secpeaktiming <- NA; secpeakvalue <- NA;
        
        secnadirtiming <- minpred %>% arrange(PRED) %>% slice(2) %>% pull(INDEX)
        secnadirvalue <-  minpred %>% arrange(PRED) %>% slice(2) %>% pull(PRED)
        
      } else{
        secpeaktiming <- NA; secpeakvalue <- NA;
        secnadirtiming <-NA; secnadirvalue <- NA;
      }
    }
    
    return( data.frame(MODEL = "4PI", 
                       PEAKTIMING = peaktiming, PEAKVALUE = peakvalue,
                       NADIRTIMING = nadirtiming, NADIRVALUE = nadirvalue,
                       PEAK2TIMING = secpeaktiming, PEAK2VALUE = secpeakvalue,
                       NADIR2TIMING = secnadirtiming, NADIR2VALUE = secnadirvalue) )
    
  } else {
    
    M <- (12/(2*pi))
    coefs<-as.vector(coef(mod))
    intercept<-coefs[grep("(Intercept)", names(coef(mod)))[1]] 
    coef_sin<-coefs[grep("sin", names(coef(mod)), ignore.case=TRUE)[1]] 
    coef_cos<-coefs[grep("cos", names(coef(mod)), ignore.case=TRUE)[1]] 
    
    # Amplitude
    amp <- sqrt(coef_sin^2 + coef_cos^2)
    
    # Angle
    ang <- coef_sin / coef_cos
    
    # Phi or Shift 
    shift <- atan(ang)
    
    # Peak Timing
    peaktiming <- ifelse(coef_cos < 0, ((shift + pi)*M),
                         ifelse(coef_cos > 0 & coef_sin > 0, (shift*M),
                                ((shift + (2 * pi))*M) ))
    
    # Peak Value
    peakvalue = intercept + amp
    
    # Nadir Timing
    nadirtiming <- ifelse(coef_cos < 0, shift + (pi/2),
                          ifelse(coef_cos > 0 & coef_sin > 0, (shift/2), (shift + (2 * pi))/2 ))
    
    # Nadir Value
    nadirvalue = intercept - amp
    
    return(data.frame(MODEL = "2PI", 
                      PEAKTIMING = peaktiming, PEAKVALUE = peakvalue,
                      NADIRTIMING = nadirtiming, NADIRVALUE = nadirvalue,
                      PEAK2TIMING = NA, PEAK2VALUE = NA,
                      NADIR2TIMING = NA, NADIR2VALUE = NA))
    
  } 
}

# Function to calculate seasonality 
seasonalitycalc <- function(df, tfield, omega, outcome){

  ############ PROCESS TIME FIELDS ################
  
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
      left_join(indexdf, by="DATE")
    
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
      left_join(indexdf, by="DATE")
  }
  
  ############ CREATE PREDICTION DATA FRAME ################
  
  # Generate a prediction data frame
  preddf <- data.frame(INDEX=seq(1, omega+1, by=1/omega)) %>% 
    mutate(SIN2PI = sin(2*pi*(1/omega)*INDEX),
           COS2PI = cos(2*pi*(1/omega)*INDEX),
           SIN4PI = sin(4*pi*(1/omega)*INDEX),
           COS4PI = cos(4*pi*(1/omega)*INDEX)) %>%
    mutate(INDEX2 = INDEX^2, INDEX3 = INDEX^3)
  
  
  ########### DECIDE WHETHER 2PI OR 4PI MODEL IS APPROPRIATE ###########

  # Start with 4pi harmonics
  m2 <- glm(value ~ SIN2PI + COS2PI + SIN4PI + COS4PI + INDEX + 
              INDEX2 + INDEX3, data = timedf)
  
  # If any 4pi terms are significant, use 4pi specification
  mcheck <- broom::tidy(m2) %>%
    filter(grepl(".*4PI$", term)) %>%
    filter(p.value<0.05)
  
  if(nrow(mcheck)>0){
    
    # Calculate predicted seasonal curve
    preddf <- preddf %>% mutate(PRED = predict(m2, newdata=preddf, type="response")) 
    
    # Calculate peak timing
    pt <- peaktimecalc(m2) 
    
  } else{
    
    m1 <- glm(value ~ SIN2PI + COS2PI + TSINDEX, data = FG_SUB)
    
    # Calculate peak timing
    pt <- peaktimecalc(m1) 
    
    # Calculate predicted seasonal curve
    preddf_new <- preddf %>% 
      mutate(PREDICTED = predict(m1, newdata=preddf, type="response")) 
    
  }
  
  # Graphical check
  #ggplot(preddf, aes(x=INDEX, y=PRED)) + geom_line()
  
  # Get average of values at end of seasonal cycle / overlapping time periods
  preddf <- preddf %>%
    mutate(INDEX = ifelse(INDEX>omega, INDEX-omega, INDEX)) %>%
    group_by(INDEX) %>%
    summarize(PRED = mean(PRED, na.rm=T)) %>%
    ungroup()

  # Return estimated peak timings and predicted values
  return(list(pt, preddf))
  
}