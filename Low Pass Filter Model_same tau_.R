#  Run model for different parameters
#make dataframe of time, voltage, integral for simple RC circuit. 
#plot these.
library(gridExtra)
install.packages("flux")
library(flux)
library(deSolve)

setwd("/Users/Lukasz/Documents/Lim Lab/Manuscripts/Cancer profiling/Model")

#function to make models and their integrals of first order RC circuits
#enter in values for period, pulse lenghts, resistor, and capacitor
#outputs a list of the RC model for each condition, , its value above a threshold, its digital threshold response, and integral of both
firstOrderRC <-function(inputPeriodVals, inputPulseVals, resistor, capacitor, threshold){
  times <- seq(0, 500, by = 1)
  signal <- as.data.frame(list(times = times, voltage = rep(0, length(times))))
  
  traceList = list() # create list to store all traces
  
  #inputPeriod = 60
  #inputPulse = 10
  
  
  for (i in 1:length(inputPeriodVals)){ # for each period
    inputPeriod = inputPeriodVals[i]
    for (j in 1:length(inputPulseVals)){ # for each pulse
      inputPulse = inputPulseVals[j]
      
      #define input signal
      signal$voltage <- ifelse((trunc(signal$times) %% inputPeriod> inputPulse ), 0, 1)
      volts = signal$voltage
      
      #make input function to put into solver
      input <- approxfun(signal, rule = 2)
      
      ##make matrix holder for the model output
      out = matrix(0,length(times),2)
  
      #params = c(R = 2, C = 2)
      
      #placeholder for parameters
      params = c(R = resistor, C = capacitor)
      
      #set "initial condition"
        y = c(Q = 0)
      
      dynamicModel <- function(time, y, params){
        with(as.list(c(y,params)),{
          import = input(time)
          dQ = -(Q/(R*C))  + (import/R)
          list(dQ)
        })
      }
      
      out <- ode(y, times, dynamicModel, params)
     # out = results[,c(1:2)]
      colnames(out) = c("time", "Q")
    
    
 
    
      #set threshold to interpret ppErk above
      thresh = threshold
   
      

      #make a dataframe to hold times, input signal, model, and integral

      traceFrame <-  data.frame(time = times,
                               inputs = signal$voltage,
                               models = out[,"Q"], 
                               anAbThresh = ifelse(out[,"Q"] > thresh, out[,"Q"], 0),
                               anAbThreshIntegral = cumsum(ifelse(out[,"Q"] > thresh, out[,"Q"], 0)),
                               digAbThresh = ifelse(out[,"Q"] > thresh, 1, 0),
                               digAbThreshIntegral = cumsum(ifelse(out[,"Q"] > thresh, 1, 0))
                               )
      
      #traceList[[paste("int", inputPeriod, "pulse", inputPulse, "R1", R1, "C1", C1, "R0", R0, "C0", C0, sep = "")]] = traceFrame
      traceList[[paste("int", inputPeriod, "pulse", inputPulse, "Resistor", resistor, "Capacitor", capacitor, sep = "")]] = traceFrame
    }
  }
  return(traceList)
}

pulseList = firstOrderRC(inputPeriodVals = 60,
                           inputPulseVals = c(10,20),
                           resistor = 2,
                           capacitor = 2,
                           threshold = 0.4)




#automate a scanning R, C values, and pulse widths. make 3x3 matrix
totModel = list()
periodVals = 60
pulseVals = c(10,20,30,40)
rVals = c(2, 20, 40)
cVal = 1
thresh = 0.4

for (i in 1:length(rVals)){
    
  model = firstOrderRC(inputPeriodVals = periodVals,
                                       inputPulseVals = pulseVals,
                                       resistor = rVals[i],
                                       capacitor = cVal,
                                       threshold = thresh)
    
    #model[[paste("int 60 pulse 10 R1_", rRise[i], " R0_", rFall[j] sep = "")]] = model
    namedModel = list()
    namedModel[[paste("R_", rVals[i], "_C_", cVal, sep = "")]] = model
    totModel = c(totModel, model)
}

# now  graph that shit
plotModels = function(list, columns, rows){
  par(mfrow=c(columns,rows))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"models"], type = "l", col = "red", ylim = c(0,1), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotModels(totModel, columns = length(rVals), rows = length(pulseVals))

plotThreshModels = function(list, columns, rows){
  par(mfrow=c(columns, rows))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"anAbThresh"], type = "l", col = "red", ylim = c(0,1), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotThreshModels(totModel, columns = length(rVals), rows = length(pulseVals))

plotDigModels = function(list, columns, rows){
  par(mfrow=c(columns, rows))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"digAbThresh"], type = "l", col = "red", ylim = c(0,1), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotDigModels(totModel, columns = length(rVals), rows = length(pulseVals))

plotDigIntegrals = function(list, columns, rows){
  par(mfrow=c(rows,columns))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"digAbThreshIntegral"], type = "l", col = "red", ylim = c(0,600), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotDigIntegrals(totModel, columns = length(pulseVals), rows = length(rVals))


#now make interpolated surfaces
#first change model to have input of ON pulse and OFF pulse, not ON pulse and period
firstOrderRC_ON_OFF <-function(ONpulseVal, OFFpulseVal, resistor, capacitor, threshold){
  times <- seq(0, 500, by = 1)
  signal <- as.data.frame(list(times = times, voltage = rep(0, length(times))))
  
  traceList = list() # create list to store all traces

  for (i in 1:length(ONpulseVal)){ # for each ON pulse
    ONpulse = ONpulseVal[i]
    for (j in 1:length(inputPulseVals)){ # for each OFF pulse
      OFFpulse = OFFpulseVal[j]
      
      #define input signal
      signal$voltage <-  rep_len(c(rep(1,ONpulse), rep(0,OFFpulse)), length.out = length(times))
      volts = signal$voltage
      
      #make input function to put into solver
      input <- approxfun(signal, rule = 2)
      
      ##make matrix holder for the model output
      out = matrix(0,length(times),2)
      
      #params = c(R = 2, C = 2)
      
      #placeholder for parameters
      params = c(R = resistor, C = capacitor)
      
      #set "initial condition"
      y = c(Q = 0)
      
      dynamicModel <- function(time, y, params){
        with(as.list(c(y,params)),{
          import = input(time)
          dQ = -(Q/(R*C))  + (import/R)
          list(dQ)
        })
      }
      
      out <- ode(y, times, dynamicModel, params)
      # out = results[,c(1:2)]
      colnames(out) = c("time", "Q")
      
      
      
      
      #set threshold to interpret ppErk above
      thresh = threshold
      
      
      
      #make a dataframe to hold times, input signal, model, and integral
      
      traceFrame <-  data.frame(time = times,
                                inputs = signal$voltage,
                                models = out[,"Q"], 
                                modelIntegral = cumsum(out[,"Q"]),
                                anAbThresh = ifelse(out[,"Q"] > thresh, out[,"Q"], 0),
                                anAbThreshIntegral = cumsum(ifelse(out[,"Q"] > thresh, out[,"Q"], 0)),
                                digAbThresh = ifelse(out[,"Q"] > thresh, 1, 0),
                                digAbThreshIntegral = cumsum(ifelse(out[,"Q"] > thresh, 1, 0))
      )
      
      #traceList[[paste("int", inputPeriod, "pulse", inputPulse, "R1", R1, "C1", C1, "R0", R0, "C0", C0, sep = "")]] = traceFrame
      traceList[[paste("on_", ONpulse, "_off_", OFFpulse, "Resistor", resistor, "Capacitor", capacitor, sep = "")]] = traceFrame
    }
  }
  return(traceList)
}

#test
firstOrderRC_ON_OFF(c(10,30),20,2,2,.5)

#simulate across ON pulses with 10 min off. Make a vector of integrated thresh after 500 mins
totModel = list()
ONpulses = c(5,10,20,30,40,50,60,90,120)
OFFpulses = 10
rVals = c(2, 20, 40)
cVal = 1
thresh = 0.4

#simulate different ON pusles, OF Fpulses
for (i in 1:length(rVals)){
  
  model = firstOrderRC_ON_OFF(ONpulseVal = ONpulses,
                            OFFpulseVal = OFFpulses,
                            resistor = rVals[i],
                            capacitor = cVal,
                            threshold = thresh)
 
  totModel = c(totModel, model)
}
digitalIntegrated = matrix(nrow = length(ONpulses), ncol = length(rVals))
analogIntegrated = matrix(nrow = length(ONpulses), ncol = length(rVals))

#extract the digital integrated signal above a threshold
for(i in 1:length(rVals)){
  for (j in 1:length(ONpulses)){
    digitalIntegrated[j,i] = totModel[[(i-1)*length(ONpulseVals)+j]]$digAbThreshIntegral[501]
  }
}

#extract the analog integrated signal
for(i in 1:length(rVals)){
  for (j in 1:length(ONpulses)){
    analogIntegrated[j,i] = totModel[[(i-1)*length(ONpulseVals)+j]]$modelIntegral[501]
  }
}

###NEXT, expand to make 10x10 matrix for several values of R. 

pdf("3_3 params model_CORRECT_60OFF_10ON.pdf",width=8,height=8)
plotModels(totModel)
dev.off()


