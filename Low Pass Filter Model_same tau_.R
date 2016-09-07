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
                           threshold = 0.65)




#automate a scanning R, C values, and pulse widths. make 3x3 matrix
totModel = list()
periodVals = 60
pulseVals = c(10, 30, 60)
rVals = c(0.2, 2, 20)
cVal = 1

for (i in 1:length(rVals)){
    model = firstOrderRC(inputPeriodVals = periodVals,
                                       inputPulseVals = pulseVals,
                                       resistor = rVals[i],
                                       capacitor = cVal,
                                       threshold = 0.65)
    
    #model[[paste("int 60 pulse 10 R1_", rRise[i], " R0_", rFall[j] sep = "")]] = model
    namedModel = list()
    namedModel[[paste("R_", rVals[i], "_C_", cVal, sep = "")]] = model
    totModel = c(totModel, model)
}

# now  graph that shit

plotModels = function(list1){
  par(mfrow=c(2,1))
  for (i in 1:length(list1)){
    plot(list1[[i]][[1]][,"time"],list1[[i]][[1]][,"model"], type = "l", col = "red", ylim = c(0,2), main = names(list1[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotIntegrals = function(list1){
  par(mfrow=c(2,1))
  for (i in 1:length(list1)){
    plot(list1[[i]][[1]][,"time"],list1[[i]][[1]][,"integral"], type = "l", col = "red", ylim = c(1,1000), main = names(list1[i]))
    #lines(list2[[i]][[1]][,"time"],list2[[i]][[1]][,"integral"], col = "purple")
  }
}

pdf("3_3 params model_CORRECT_60OFF_10ON.pdf",width=8,height=8)
plotModels(totModel)
dev.off()

pdf("3_3 params integral_60OFF_30ON.pdf",width=8,height=8)
plotIntegrals(totModel)
dev.off()

pdf("2_1 params model_60OFF_10ON.pdf",width=8,height=8)
plotModels(totModel)
dev.off()

pdf("2_1 params integral_60OFF_10ON.pdf",width=8,height=8)
plotIntegrals(totModel)
dev.off()

save.image("3x3 model matrix_correct model.RData")




pdf("diffParam_r2c2r2c2.pdf",width=6,height=8)
plotGrid(r2c2r2c2)
dev.off()




head(r2c2r2c2)
plotGrid = function(list){
  par(mfrow=c(3,2))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"model"], type = "l", ylim = c(0,1.4), main = names(list[i]))
    plot(list[[i]][,"time"],list[[i]][,"integral"], type = "l", ylim = c(0,800), main = names(list[i]))
  }
}

plotModels = function(list1, list2){
  par(mfrow=c(3,1))
  for (i in 1:length(list1)){
    plot(list1[[i]][,"time"],list1[[i]][,"model"], type = "l", col = "red", main = names(list1[i]))
    lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}

plotIntegrals = function(list1, list2){
  par(mfrow=c(3,1))
  for (i in 1:length(list1)){
    plot(list1[[i]][,"time"],list1[[i]][,"integral"], type = "l", col = "red", ylim = c(1,800), main = names(list1[i]))
    lines(list2[[i]][,"time"],list2[[i]][,"integral"], col = "purple")
  }
}

plotThresh = function(list1, list2){
  par(mfrow=c(3,1))
  for (i in 1:length(list1)){
    plot(list1[[i]][,"time"],list1[[i]][,"abThresh"], type = "l", col = "red", ylim = c(0,1.2), main = names(list1[i]))
    lines(list2[[i]][,"time"],list2[[i]][,"abThresh"], col = "purple")
  }
}

