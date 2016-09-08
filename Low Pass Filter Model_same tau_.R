#  Run model for different parameters
#make dataframe of time, voltage, integral for simple RC circuit. 
#plot these.
library(gridExtra)
install.packages("flux")
library(flux)
install.packages("deSolve")
library(deSolve)
install.packages("ggplot2")
library(ggplot2)
library(reshape2)
install.packages("RColorBrewer")
library(RColorBrewer)
library(mgcv)
install.packages("akima")
library(akima)


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

#model that takes ON pulse and OFF pulse as inputs
firstOrderRC_ON_OFF <-function(ONpulseVal, OFFpulseVal, resistor, capacitor, threshold){
  times <- seq(0, 500, by = 1)
  signal <- as.data.frame(list(times = times, voltage = rep(0, length(times))))
  
  traceList = list() # create list to store all traces
  
  for (i in 1:length(ONpulseVal)){ # for each ON pulse
    ONpulse = ONpulseVal[i]
    for (j in 1:length(OFFpulseVal)){ # for each OFF pulse
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


#automate a scanning R, C values, and pulse widths. make 3x3 matrix
totModel = list()
periodVals = 40
pulseVals = c(0,20,40)
rVals = c(2, 32)
cVal = 1
thresh = 0.3

#simulate
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

# now  graph!
#output traces
plotModels = function(list, columns, rows){
  par(mfrow=c(columns,rows))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"models"], type = "l", col = "red", ylim = c(0,1), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotModels(totModel, columns = length(rVals), rows = length(pulseVals))

pdf("2Rvals0_20_40_onModels.pdf",width=10,height=5)
plotModels(totModel, columns = length(rVals), rows = length(pulseVals))
dev.off()

#output traces above threshod
plotDigModels = function(list, columns, rows){
  par(mfrow=c(columns, rows))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"digAbThresh"], type = "l", col = "red", ylim = c(0,1), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotDigModels(totModel, columns = length(rVals), rows = length(pulseVals))

pdf("2Rvals0_20_40_onDigModels.pdf",width=10,height=5)
plotDigModels(totModel, columns = length(rVals), rows = length(pulseVals))
dev.off()


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
    plot(list[[i]][,"time"],list[[i]][,"digAbThreshIntegral"], type = "l", col = "red", ylim = c(0,500), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}
plotDigIntegrals(totModel, columns = length(pulseVals), rows = length(rVals))

pdf("2Rvals0_20_40_onDigModelsIntegral.pdf",width=10,height=5)
plotDigIntegrals(totModel, columns = length(pulseVals), rows = length(rVals))
dev.off()


plotModelThreshIntegral = function(list, column,rows){
  par(mfrow=c(rows,columns))
  for (i in 1:length(list)){
    plot(list[[i]][,"time"],list[[i]][,"model"], type = "l", col = "red", ylim = c(0,600), main = names(list[i])
    plot(list[[i]][,"time"],list[[i]][,"digAbThresh"], type = "l", col = "red", ylim = c(0,600), main = names(list[i])
    plot(list[[i]][,"time"],list[[i]][,"digAbThreshIntegral"], type = "l", col = "red", ylim = c(0,600), main = names(list[i]))
    #lines(list2[[i]][,"time"],list2[[i]][,"model"], col = "purple")
  }
}

#now make interpolated surfaces
#first change model to have input of ON pulse and OFF pulse, not ON pulse and period
firstOrderRC_ON_OFF <-function(ONpulseVal, OFFpulseVal, resistor, capacitor, threshold){
  times <- seq(0, 500, by = 1)
  signal <- as.data.frame(list(times = times, voltage = rep(0, length(times))))
  
  traceList = list() # create list to store all traces
  
  for (i in 1:length(ONpulseVal)){ # for each ON pulse
    ONpulse = ONpulseVal[i]
    for (j in 1:length(OFFpulseVal)){ # for each OFF pulse
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
  
  
    digitalIntegrated = matrix(nrow = length(ONpulses), ncol = length(rVals))
    analogIntegrated = matrix(nrow = length(ONpulses), ncol = length(rVals))
    
    #extract the digital integrated signal above a threshold
    for(i in 1:length(rVals)){
      for (j in 1:length(ONpulses)){
        digitalIntegrated[j,i] = totModel[[(i-1)*length(ONpulses)+j]]$digAbThreshIntegral[501]
      }
    }
    
    #extract the analog integrated signal
    for(i in 1:length(rVals)){
      for (j in 1:length(ONpulses)){
        analogIntegrated[j,i] = totModel[[(i-1)*length(ONpulses)+j]]$modelIntegral[501]
      }
    }


###NEXT, expand to make 10x10 matrix for several values of R. 
##First do for just multiple thresholds, multiple values of R

threshholds = c(0.1,0.2,0.3, 0.4, 0.6, 0.8)    
for (a in 1:length(threshholds)){   
  ONpulses = c(1,5,10,20,30,40,50,60,90,120)
  OFFpulses = c(1,5,10,20,30,40,50,60,90,120)
  rVals = c(2,4,8,16,32, 64)
  cVal = 1
  thresh = threshholds[a]
  
  #simulate different ON pusles, OF Fpulses, for different R
  on_off_Digital = list()
  on_off_Analog = list()
  
  #simulate
  totModel = list()
  for (z in 1:length(rVals)){
    
    model = firstOrderRC_ON_OFF(ONpulseVal = ONpulses,
                                OFFpulseVal = OFFpulses,
                                resistor = rVals[z],
                                capacitor = cVal,
                                threshold = thresh)
    
    digitalIntegrated = matrix(nrow = length(ONpulses), ncol = length(OFFpulses))
    analogIntegrated = matrix(nrow = length(ONpulses), ncol = length(OFFpulses))
    
    
    #extract the digital integrated signal above a threshold
    for(i in 1:length(ONpulses)){
      for (j in 1:length(OFFpulses)){
        digitalIntegrated[j,i] = model[[(i-1)*length(ONpulses)+j]]$digAbThreshIntegral[501]
      }
    }
    on_off_Digital[[paste("R_",rVals[z],sep="")]] = digitalIntegrated
    
    #extract the analog integrated signal
    for(i in 1:length(ONpulses)){
      for (j in 1:length(OFFpulses)){
        analogIntegrated[j,i] = model[[(i-1)*length(ONpulses)+j]]$modelIntegral[501]
      }
    }
    on_off_Analog[[paste("R_",rVals[z],sep="")]] = analogIntegrated
  }
  
  ## Make list of interpolated arrays 
  onOffDigInterp = list()
  for (i in 1:length(on_off_Digital)){
    mat = on_off_Digital[[i]]
    colnames(mat) = ONpulses #label columns with ON times
    tempdf = data.frame(mat, check.names = FALSE)
    tempdf["OFFtimes"] = OFFpulses #put OFF times for melting purposes
    
    #melt this to interpolate
    meltdf = melt(tempdf, id = "OFFtimes")
    colnames(meltdf) = c("OFFtimes", "ONtimes", "value")
    meltdf$ONtimes = as.numeric(as.character(meltdf$ONtimes)) # need to make ONtimes be numeric, not factors
    
    akInterp = interp(meltdf[,"ONtimes"],
                          meltdf[,"OFFtimes"],
                          meltdf[,"value"], 
                          xo = seq(0, 
                                   max(meltdf[,"ONtimes"]), length = 121), 
                          yo = seq(0, max(meltdf[,"OFFtimes"]), length = 121))
    
    interpDf = data.frame(akInterp[[3]]) # extract the interpolated matrix
    colnames(interpDf) = seq (0,120,by =1)
    rownames(interpDf) = seq (0,120,by =1)
    interpDf["OnTime"] = seq(0,120,1) #label ON times for melthing purposes
    
    #melt for ggplot purposes
    meltInterp = melt(interpDf, id = "OnTime")
    colnames(meltInterp) = c("OnTime", "OffTime", "value")
    meltInterp["Rval"] = names(on_off_Digital[i])
    
    onOffDigInterp[[names(on_off_Digital)[i]]] = meltInterp
  }
  
  
  #combine all maps so that can ggplot them together
  totInterpMaps = data.frame() 
  for (i in 1:length(onOffDigInterp)){
    totInterpMaps = rbind(totInterpMaps, onOffDigInterp[[i]])
  }
  
  totInterpMaps$Rval = factor(totInterpMaps$Rval, levels = c("R_2", "R_4", "R_8","R_16","R_32", "R_64"))
  
 
  
  #plot!
  plotInterp(subset(totInterpMaps))
  ggsave(paste("LPF model digital_Thresh ", thresh,".pdf",sep = ""), height = 5, width = 8)

}

#plotting function
plotInterp = function(df){
  ggplot(data = df, aes(x = OnTime, y = OffTime))+
    geom_raster(aes(fill = value), interpolate = TRUE)+
    facet_wrap(~Rval)+ 
    scale_y_discrete(breaks = c(0,20,40,60,80,100,120))+
    scale_x_continuous(breaks = c(20,40,60,80,100,120))+
    scale_fill_distiller(type = "div", palette = "PuOr")+
    theme_bw()+ 
    labs(title = paste("Threshold is", thresh, sep = " " ))+
    theme(
      axis.text.x= element_text(angle = 0, face = "bold", size = 15),
      axis.text.y= element_text(angle = 0, face = "bold", size = 15),
      axis.title.x= element_text(angle = 0, face = "bold", size = 20),
      axis.title.y= element_text(angle = 90, face = "bold", size = 20),
      #  axis.line= element_line(size = 10),
      axis.line = element_line(size = 1.5),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())
}
