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
#outputs a list of the RC model for each condition, and also its integral
firstOrderRC_d_Rise_d_fall <-function(inputPeriodVals, inputPulseVals, resistor0, capacitor0, resistor1, capacitor1, threshold){
  times <- seq(0, 500, by = 1)
  signal <- as.data.frame(list(times = times, voltage = rep(0, length(times))))
  
  traceList = list() # create list to store all traces
  for (i in 1:length(inputPeriodVals)){ # for each period
    inputPeriod = inputPeriodVals[i]
    for (j in 1:length(inputPulseVals)){ # for each pulse
      inputPulse = inputPulseVals[j]
      
      #inputPeriod = 40
      #inputPulse = 10
      
      ##TEST
      
      #signal <- as.data.frame(list(times = times, voltage = rep(0, length(times))))
      #inputPeriod = 20
      #inputPulse = 10
      #resistor = 1
      #capacitor = 2
      
      #define input signal
      signal$voltage <- ifelse((trunc(signal$times) %% inputPeriod> inputPulse ), 0, 1)
      volts = signal$voltage
      
      #make input function to put into solver
      input <- approxfun(signal, rule = 2)
      
      ##make a diff model if signal is 1 or 0.
      out = matrix(0,length(times),2)
      
      #placeholder for parameters
      params = c(R = 0, C = 0)
      newTimes = 0
      #start indexing time and volt vector at 1
      i = 1
      
      #keep track of how may times have been in voltage of 1 or 0
      
      while (i < length(times)){
        if(volts[i] == 1){
          newTimes = seq(times[i], times[length(times)])
          newVolts = volts[c(i: length(volts))]
          
          timeLim = newTimes[which(newVolts == 0)[1]]
          
          #to prevent an NA at end of timecourse
          if (is.na(timeLim)){
            timeLim = length(times)
          }
          
          newTimes = seq(times[i], times[timeLim])
          
          R1 = resistor1
          C1 = capacitor1
          params <- c(R = R1, C = C1)
        }
        else if(volts[i] == 0){
          newTimes = seq(times[i], times[length(times)])
          newVolts = volts[c(i: length(volts))]
          
          timeLim = newTimes[which(newVolts == 1)[1]]
          
          #to prevent an NA at end of timecourse
          if (is.na(timeLim)){
            timeLim = length(times)
          }
          
          newTimes = seq(times[i], times[timeLim])
          
          R0 =resistor0
          C0 = capacitor0
          params <- c(R = R0, C = C0)
        }
        
        #set "initial condition"
        if (i == 1){
          y = c(Q = 0)
        }
        else{
          y = c(Q = out[i-1,2])
        }
        
        
        dynamicModel <- function(time, y, params){
          with(as.list(c(y,params)),{
            import = input(time)
            dQ = -(Q/(R*C))  + (import/R)
            list(dQ)
          })
        }
        
        outSegment <- ode(y, newTimes, dynamicModel, params)
        out[c(i:(i+length(newTimes)-1)),] = outSegment
        
        i = timeLim+1
      }
      
      colnames(out) = c("time", "Q")
      #trace = approxfun(out[,"time"], out[,"Q"])
      
      #function to calculate integral
      #intFun = Vectorize(function(x) integrate(Vectorize(trace), 0, x, stop.on.error = FALSE)$value)
      
      
      times_minus = times[c(1:length(times-1))] 
      
      #make a dataframe to hold times, input signal, model, and integral
      
      #first set threshold to interpret ppErk above
      thresh = threshold
      
      #now make the dataframe
      traceFrame =  data.frame(time= times_minus,
                               input = signal$voltage[c(1:length(times_minus))],
                               model = out[c(1:length(times_minus)),"Q"], 
                               integral = cumsum(out[c(1:length(times_minus)),"Q"]),
                               abThresh = ifelse(out[c(1:length(times_minus)),"Q"] > thresh, 1, 0))
      
      #traceList[[paste("int", inputPeriod, "pulse", inputPulse, "R1", R1, "C1", C1, "R0", R0, "C0", C0, sep = "")]] = traceFrame
      traceList[[paste("int", inputPeriod, "pulse", inputPulse, sep = "")]] = traceFrame
      
    }
  }
  return(traceList)
}


rRise = c(1)
rFall = c(2,20)

#automate a scan of R_rise and R_fall. make 3x3 matrix
totModel = list()
for (i in 1:length(rRise)){
  for (j in 1:length(rFall)){
    model = firstOrderRC_d_Rise_d_fall(inputPeriodVals = 60,
                                       inputPulseVals = c(10),
                                       resistor0 = rFall[j],
                                       capacitor0 = 2,
                                       resistor1 = rRise[i],
                                       capacitor1 = 2,
                                       threshold = 0.65)
    #model[[paste("int 60 pulse 10 R1_", rRise[i], " R0_", rFall[j] sep = "")]] = model
    totModel[[paste("R1_", rRise[i], " R0_", rFall[j], sep = "")]] = model
    
  }
}
totModel = c(totModel)

###Next to automate: do this over 3 different pulse timings
# but first graph that shit

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

