#1) Clear memory and set working directory
{
  #Clear memory
  rm(list=ls())
  
  #Set working directory
  setwd("/media/felix/Elements/Felix_Data_Backup/Science/PostDoc Project/Experiments/1_Experiment_1_Frequency_of_sex/")
}

#2) Load packages and scripts
{
  library(tidyverse)
  library(nlsMicrobio)
  # library(devtools)
  # install_github("femoerman/PBPGM")
  #library(PBPGM)
  library(nlme)
  library(MuMIn)
  library(car)
  library(ggfortify)
  source("/media/felix/Elements/Felix_Data_Backup/Science/PostDoc Project/Experiments/0_r_functions/process_xpt.R")
}

#3) Read in the data, and combine with the treatment data from the experiment
{
  xpt.files <- list.files("2_data/", full.names = T, pattern = ".txt")
  
  OD.750 <- process_xpt(xpt.files[1])
  for (i in 2:length(xpt.files)){
    OD.750 <- rbind(OD.750, process_xpt(xpt.files[i]))
  }
  
  #Load data with treatment information
  treatments <- read.csv(file="1_setup/SampleList_FinalAssays_redone.csv", sep=",")
  treatments$Well <- paste0(treatments$Row, treatments$Column)
  treatments$Barcode <- treatments$Plate
  treatments <- dplyr::select(treatments, -Row, -Column, -Plate)
  
  #Combine all the data
  OD.data <- left_join(OD.750, treatments, by=c("Barcode", "Well"))
}

#4) Standardize OD measurements for the blank values (median per plate), and add a time since start value
{
  #Create standardized OD values and add a time variable
  OD.data$OD650_Stand <- 0
  OD.data$OD750_Stand <- 0
  OD.data$time <- 0
  OD.data$timepoint <- 0
  
  for (i in unique(OD.data$Plate.ID)){
    temp <- filter(OD.data, Plate.ID==i)
    for (j in unique(temp$hours)) {
      temp2 <- filter(temp, hours==j)
      blank <- filter(temp, Population_No=="Blank", hours==j)
      OD.data[which(OD.data$Plate.ID==i & OD.data$hours==j), "OD650_Stand"] <- temp2$OD650 - min(blank$OD650)
      OD.data[which(OD.data$Plate.ID==i & OD.data$hours==j), "OD750_Stand"] <- temp2$OD750 - min(blank$OD750)
      OD.data[which(OD.data$Plate.ID==i & OD.data$hours==j), "time"] <- temp2$hours - min(temp$hours)
      
    }
    OD.data[which(OD.data$Plate.ID==i), "timepoint"] <- match(OD.data[which(OD.data$Plate.ID==i), "time"], sort(unique(OD.data[which(OD.data$Plate.ID==i), "time"])))
    
  }
  
  OD.data$OD650_Stand <- ifelse(OD.data$OD650_Stand<0, 0, OD.data$OD650_Stand)
  OD.data$OD750_Stand <- ifelse(OD.data$OD750_Stand<0, 0, OD.data$OD750_Stand)
  
  #Add an ID for each replicate measurement
  OD.data$ID <- paste0(OD.data$Population_No, OD.data$Assay_stress, OD.data$Replicate, OD.data$Batch)
  
  #Filter the data to exclude the blanks
  OD.data.filtered <- filter(OD.data, Population_No != "Blank")
  
  #Order the sexual reproduction treatment by frequency, and create a numeric value for it
  OD.data.filtered$Reproduction <- ifelse(OD.data.filtered$Reproduction == "biweekly", "Biweekly", OD.data.filtered$Reproduction)
  OD.data.filtered$Reproduction <- factor(OD.data.filtered$Reproduction, levels = c("None", "Once", "Monthly", "Biweekly", "Ancestor"))
  OD.data.filtered$Reproduction_Numeric <- ifelse(OD.data.filtered$Reproduction=="None", 0, ifelse(OD.data.filtered$Reproduction=="Once", 1, ifelse(OD.data.filtered$Reproduction=="Monthly", 4, ifelse(OD.data.filtered$Reproduction=="Biweekly", 8, -1))))
}

#5) Plot the OD values for the populations
{
  ggplot(OD.data.filtered, aes(x = time, y = OD650, ident = ID, colour=Reproduction)) + geom_point() + geom_line() + facet_grid(Assay_stress~Historical_Stress)
  
  ggplot(OD.data.filtered, aes(x = time, y = OD750, ident = ID, colour=Reproduction)) + geom_point() + geom_line() + facet_grid(Assay_stress~Historical_Stress)
  
  ggplot(OD.data.filtered, aes(x = time, y = OD650_Stand, ident = ID, colour=Reproduction)) + geom_point() + geom_line() + facet_grid(Assay_stress~Historical_Stress)
  
  ggplot(OD.data.filtered, aes(x = time, y = OD750_Stand, ident = ID, colour=Reproduction)) + geom_point() + geom_line() + facet_grid(Assay_stress~Historical_Stress)
  
}

#6) Calculate metrics (max r; population growth model fits) from the raw data, and save the output
{
  
  {
    #Order by timepoint
    OD.data.filtered <- arrange(OD.data.filtered, timepoint)
    
    #Create an output dataframe
    data.output <- filter(OD.data.filtered, time==0) %>% dplyr::select(-OD650, -OD750, -OD650_Stand, -OD750_Stand)
    data.output$max_r_OD750 <- 0
    data.output$max_r_OD650 <- 0
    data.output$max_K_OD750 <- 0
    data.output$max_K_OD650 <- 0
    
    #Method 1: Loop over the timepoints and calculate the max r0 and K
    for (i in data.output$ID){
      temp <- filter(OD.data.filtered, ID==i)
      OD750 <- (log(temp[2:nrow(temp), "OD750_Stand"]) - log(temp[1:nrow(temp)-1, "OD750_Stand"])) / (temp[2:nrow(temp), "time"]-temp[1:nrow(temp)-1, "time"])
      OD650 <- (log(temp[2:nrow(temp), "OD650_Stand"]) - log(temp[1:nrow(temp)-1, "OD650_Stand"])) / (temp[2:nrow(temp), "time"]-temp[1:nrow(temp)-1, "time"])
      OD750 <- OD750[which(! is.nan(OD750) & ! is.infinite(OD750) & OD750 != -Inf)]
      OD650 <- OD650[which(! is.nan(OD650) & ! is.infinite(OD650) & OD650 != -Inf)]
      
      
      data.output[which(data.output$ID==i), "max_r_OD750"] <- max(OD750)
      data.output[which(data.output$ID==i), "max_r_OD650"] <- max(OD650)
      data.output[which(data.output$ID==i), "max_K_OD750"] <- max(temp$OD750_Stand)
      data.output[which(data.output$ID==i), "max_K_OD650"] <- max(temp$OD650_Stand)
    }
    
    #Replace -INF with NA
    data.output[data.output == -Inf] <- NA
  }
  
}

#7) Plot the output data based on the max r0 and K values
{
  #Plot all data
  ggplot(data.output, aes(y=max_r_OD750, x = Historical_Stress, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  ggplot(data.output, aes(y=max_r_OD650, x = Historical_Stress, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  
  ggplot(data.output, aes(y=max_K_OD750, x = Historical_Stress, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  ggplot(data.output, aes(y=max_K_OD650, x = Historical_Stress, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  
  
  #Plot only for the environment where populations evolved
  data.output.filt <- data.output[which(data.output$Population_No=="Ancestor" | data.output$Assay_stress==data.output$Historical_Stress),]
  
  ggplot(data.output.filt, aes(y=max_r_OD750, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  ggplot(data.output.filt, aes(y=max_r_OD650, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  ggplot(data.output.filt, aes(y=max_K_OD750, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  ggplot(data.output.filt, aes(y=max_K_OD650, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  
  
  #Standardize by the ancestor value (median), but now for the complete dataset, not just the local adaptation
  data.output.stand <- mutate(data.output, r_OD750_stand = 0, r_OD650_stand = 0, K_OD750_stand = 0, K_OD650_stand = 0)
  for (i in unique(data.output.stand$ID)){
    line <- filter(data.output.stand, ID==i)
    temp <- filter(data.output.stand, Population_No == "Ancestor", Assay_stress == line$Assay_stress)
    data.output.stand[which(data.output.stand$ID==i), "r_OD750_stand"] <- line$max_r_OD750 / median(temp$max_r_OD750)
    data.output.stand[which(data.output.stand$ID==i), "r_OD650_stand"] <- line$max_r_OD650 / median(temp$max_r_OD650)
    data.output.stand[which(data.output.stand$ID==i), "K_OD750_stand"] <- line$max_K_OD750 / median(temp$max_K_OD750)
    data.output.stand[which(data.output.stand$ID==i), "K_OD650_stand"] <- line$max_K_OD650 / median(temp$max_K_OD650)
  }
  
  #Add clearer labels for the facets
  data.output.stand$Evolution <- ifelse(data.output.stand$Historical_Stress=="NaCl", "Salt lines", ifelse(data.output.stand$Historical_Stress=="None", "No salt lines", "Ancestor"))
  data.output.stand$Environment2 <- ifelse(data.output.stand$Assay_stress=="NaCl", "Salt lines", ifelse(data.output.stand$Assay_stress=="None", "No salt lines", "Ancestor"))
  data.output.stand$Environment <- ifelse(data.output.stand$Assay_stress=="NaCl", "Salt environment", ifelse(data.output.stand$Assay_stress=="None", "No salt environment", "Ancestor"))
  
  #Get the local adaptation data
  data.output.filt.stand <- filter(data.output.stand, Historical_Stress==Assay_stress)
  

  
  #Plot the standardized values for the local adaptation
  {
    ggplot(data.output.filt.stand, aes(y=r_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(data.output.filt.stand, aes(y=r_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(data.output.filt.stand, aes(y=K_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(data.output.filt.stand, aes(y=K_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=r_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=r_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=K_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=K_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=r_OD750_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=r_OD650_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=K_OD750_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
    ggplot(filter(data.output.filt.stand, Population_No != "Ancestor"), aes(y=K_OD650_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_wrap(~Assay_stress)
  }
  
  #Plot the standardozed values for the complete dataset
  {
    ggplot(data.output.stand, aes(y=r_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(data.output.stand, aes(y=r_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(data.output.stand, aes(y=K_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(data.output.stand, aes(y=K_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=r_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=r_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=K_OD750_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=K_OD650_stand, x = Reproduction, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=r_OD750_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=r_OD650_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=K_OD750_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
    ggplot(filter(data.output.stand, Population_No != "Ancestor"), aes(y=K_OD650_stand, x = Reproduction_Numeric, colour=Reproduction, group=Historical_Stress)) + geom_point() + geom_boxplot(aes(group=Population_No)) + facet_grid(Evolution~Environment)
  }
 
}

#8) Statistically analyze the data
{
  #Prepare the data for model fitting
  data.stats <- filter(data.output.filt.stand, Reproduction != "Ancestor")
  
  #8.1) First look at r0
  {
    #8.1.1) With reproduction as numeric variable
    {
      #Fit a full interaction model (with assay stress)
      full.r0.1 <-  lme(data=data.stats, r_OD750_stand~poly(Reproduction_Numeric,2)*Evolution, random =  ~ 1|Population_No, method = "ML")
      
      #Look at all possible models
      dredge(full.r0.1)
      
      #Make the best fitting model
      best.r0.1 <- lme(data=data.stats, r_OD750_stand~poly(Reproduction_Numeric,2)*Evolution, random =  ~ 1|Population_No, method = "REML")
      
      #Look at model output
      summary(best.r0.1)
      Anova(best.r0.1, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
      plot(best.r0.1, type=c("p","smooth"), col.line=1)
      plot(best.r0.1,
           sqrt(abs(resid(.)))~fitted(.),
           type=c("p","smooth"), col.line=1)
      qqnorm(resid(best.r0.1))
      
      qqline(resid(best.r0.1))
      
      
      #Create model predictions
      {
        predict.data <- expand.grid(Reproduction_Numeric=seq(from = min(data.stats$Reproduction_Numeric), to = max(data.stats$Reproduction_Numeric), length = 100),
                                    Evolution=unique(data.stats$Evolution))
        predictions <- predict(newdata=predict.data, best.r0.1, se.fit=T, level=0)
        predict.data$r0.mean <-  (predictions$fit)
        predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
        predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
      }
      
      #Plot model predictions with the raw data on top
      {
        r0.numeric.la <- ggplot(data.stats, aes(y=r_OD750_stand, x=Reproduction_Numeric)) + geom_point() + facet_wrap(~Evolution) +
          geom_line(data=predict.data, inherit.aes = F, mapping = aes(x=Reproduction_Numeric, y = r0.mean))+
          geom_ribbon(data=predict.data, inherit.aes = F, 
                       mapping = aes(x=Reproduction_Numeric, ymin=r0.lower, ymax=r0.upper, y = r0.mean),
                    stat = "identity", alpha=0.3)+ theme_bw() +
          theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), axis.text=element_text(size=14), legend.text=element_text(size=14), 
                strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
                axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
          ylab(expression(paste("Relative change in intrinsic rate of increase ", italic("r₀")))) + xlab("Number of sexual reproduction events")

        r0.numeric.la
        ggsave(r0.numeric.la, filename = "4_results/Figures/r0_numeric_local_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2, width=67*3)
        
      }
    }

    #8.1.2) With reproduction as categorical variable
    {
      #Fit a full interaction model (with assay stress)
      full.r0.2 <-  lme(data=data.stats, r_OD750_stand~Reproduction*Evolution, random =  ~ 1|Population_No, method = "ML")
      
      #Look at all possible models
      dredge(full.r0.2)
      
      #Make the best fitting model
      best.r0.2 <- lme(data=data.stats, r_OD750_stand~Reproduction*Evolution, random =  ~ 1|Population_No, method = "REML")
      
      #Look at model output
      summary(best.r0.2)
      Anova(best.r0.2, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
      plot(best.r0.2, type=c("p","smooth"), col.line=1)
      plot(best.r0.2,
           sqrt(abs(resid(.)))~fitted(.),
           type=c("p","smooth"), col.line=1)
      qqnorm(resid(best.r0.2))
      qqline(resid(best.r0.2))
      
      #Create model predictions
      {
        predict.data <- expand.grid(Reproduction=unique(data.stats$Reproduction),
                                    Evolution=unique(data.stats$Evolution))
        predictions <- predict(newdata=predict.data, best.r0.2, se.fit=T, level=0)
        predict.data$r0.mean <-  (predictions$fit)
        predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
        predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
      }
      
      #Plot model predictions with the raw data on top
      {
        r0.categorical.la <- ggplot(data.stats, aes(y=r_OD750_stand, x=Reproduction)) + geom_point() + facet_wrap(~Evolution) +
          geom_boxplot(data=predict.data, inherit.aes = F, 
                      mapping = aes(x=Reproduction, ymin=r0.lower, ymax=r0.upper, middle=r0.mean, lower=r0.lower, upper=r0.upper, fill = Reproduction),
                      stat = "identity", alpha=0.3)+ theme_bw() +
          theme(axis.text=element_text(size=14), legend.text=element_text(size=14), 
                strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
                axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
          ylab(expression(paste("Relative change in intrinsic rate of increase ", italic("r₀")))) + xlab("Frequency of sexual reproduction") + 
        scale_color_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) + 
          scale_fill_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) 
        
        r0.categorical.la
        ggsave(r0.categorical.la, filename = "4_results/Figures/r0_categorical_local_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2, width=67*3.3)
        
        
        r0.categorical.la.2 <- ggplot(predict.data, aes(y=(r0.mean-1)*100, x=Reproduction, colour=Reproduction, fill=Reproduction)) + geom_bar(stat="identity") + facet_wrap(~Evolution) +
          theme_bw() +
          theme(axis.text=element_text(size=14), legend.text=element_text(size=14), 
                strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
                axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
          ylab("Evolutionary change (%)") + xlab("Frequency of sexual reproduction") + 
          scale_color_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) + 
          scale_fill_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) 
        
        r0.categorical.la.2
        ggsave(r0.categorical.la.2, filename = "4_results/Figures/r0_categorical_local_adaptation_PintOS.png", device = "png", dpi = 500, units = "mm", height=74*2, width=67*3.3)
        
      }
      
    }
    
  }
  
  #8.2) Next look at K
  {
    #8.1.1) With reproduction as numeric variable
    {
      #Fit a full interaction model (with assay stress)
      full.K.1 <-  lme(data=data.stats, K_OD750_stand~poly(Reproduction_Numeric,2)*Evolution, random =  ~ 1|Population_No, method = "ML")
      
      #Look at all possible models
      dredge(full.K.1)
      
      #Make the best fitting model
      best.K.1 <- lme(data=data.stats, K_OD750_stand~Evolution, random =  ~ 1|Population_No, method = "REML")
      
      #Look at model output
      summary(best.K.1)
      Anova(best.K.1, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
      plot(best.K.1, type=c("p","smooth"), col.line=1)
      plot(best.K.1,
           sqrt(abs(resid(.)))~fitted(.),
           type=c("p","smooth"), col.line=1)
      qqnorm(resid(best.K.1))
      qqline(resid(best.K.1))
      
      
      #Create model predictions
      {
        predict.data <- expand.grid(Reproduction_Numeric=seq(from = min(data.stats$Reproduction_Numeric), to = max(data.stats$Reproduction_Numeric), length = 100),
                                    Evolution=unique(data.stats$Evolution))
        predictions <- predict(newdata=predict.data, best.K.1, se.fit=T, level=0)
        predict.data$r0.mean <-  (predictions$fit)
        predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
        predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
      }
      
      #Plot model predictions with the raw data on top
      {
        K.numeric.la <- ggplot(data.stats, aes(y=K_OD750_stand, x=Reproduction_Numeric)) + geom_point() + facet_wrap(~Evolution) +
          geom_line(data=predict.data, inherit.aes = F, mapping = aes(x=Reproduction_Numeric, y = r0.mean))+
          geom_ribbon(data=predict.data, inherit.aes = F, 
                      mapping = aes(x=Reproduction_Numeric, ymin=r0.lower, ymax=r0.upper, y = r0.mean),
                      stat = "identity", alpha=0.3)+ theme_bw() +
          theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), axis.text=element_text(size=14), legend.text=element_text(size=14), 
                strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
                axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
          ylab(expression(paste("Relative change in Equilibrium density ", italic("K")))) + xlab("Number of sexual reproduction events")
        
        K.numeric.la
        ggsave(K.numeric.la, filename = "4_results/Figures/K_numeric_local_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2, width=67*3)
        
      }
    }
    
    #8.1.2) With reproduction as categorical variable
    {
      #Fit a full interaction model (with assay stress)
      full.K.2 <-  lme(data=data.stats, K_OD750_stand~Reproduction*Evolution, random =  ~ 1|Population_No, method = "ML")
      
      #Look at all possible models
      dredge(full.K.2)
      
      #Make the best fitting model
      best.K.2 <- lme(data=data.stats, K_OD750_stand~Evolution, random =  ~ 1|Population_No, method = "REML")
      
      #Look at model output
      summary(best.K.2)
      Anova(best.K.2, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
      plot(best.K.2, type=c("p","smooth"), col.line=1)
      plot(best.K.2,
           sqrt(abs(resid(.)))~fitted(.),
           type=c("p","smooth"), col.line=1)
      qqnorm(resid(best.K.2))
      qqline(resid(best.K.2))
    }
    
    
    #Create model predictions
    {
      predict.data <- expand.grid(Reproduction=unique(data.stats$Reproduction),
                                  Evolution=unique(data.stats$Evolution))
      predictions <- predict(newdata=predict.data, best.K.2, se.fit=T, level=0)
      predict.data$r0.mean <-  (predictions$fit)
      predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
      predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
    }
    
    #Plot model predictions with the raw data on top
    {
      K.categorical.la <- ggplot(data.stats, aes(y=K_OD750_stand, x=Reproduction)) + geom_point() + facet_wrap(~Evolution) +
        geom_boxplot(data=predict.data, inherit.aes = F, 
                     mapping = aes(x=Reproduction, ymin=r0.lower, ymax=r0.upper, middle=r0.mean, lower=r0.lower, upper=r0.upper, fill = Reproduction),
                     stat = "identity", alpha=0.3)+ theme_bw() +
        theme(axis.text=element_text(size=14), legend.text=element_text(size=14), 
              strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
              axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
        ylab(expression(paste("Relative change in intrinsic rate of increase ", italic("K")))) + xlab("Frequency of sexual reproduction") + 
        scale_color_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) + 
        scale_fill_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) 
      
      K.categorical.la
      ggsave(K.categorical.la, filename = "4_results/Figures/K_categorical_local_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2, width=67*3.3)
      
    }
    
  }
}

#9) Perform the statistical analyses to detect whether there are trade-offs/how adaptation of populations is affected in the environment where
#they did not evolve (that is, look at the complete data, all at once)

#Here, I use quadratic numeric and categorical reproduction variables (because clear non linear trends in some groups)
{
  #Prepare the data for model fitting
  data.stats.2 <- filter(data.output.stand, Reproduction != "Ancestor")
  
  #9.1) First look at r0
  {
    
    #9.1.1) With reproduction as numerical variable
    {
      #Fit a full interaction model (with assay stress)
      full.r0.3 <- lme(data=data.stats.2, r_OD750_stand~poly(Reproduction_Numeric,2)*Evolution*Environment, random =  ~ 1|Population_No, method = "ML")
      
      #Look at all possible models
      dredge(full.r0.3)
      
      #Make the best fitting model
      best.r0.3 <- lme(data=data.stats.2, r_OD750_stand~poly(Reproduction_Numeric,2)*Environment+Evolution, random =  ~ 1|Population_No, method = "REML")
      
      #Look at model output
      summary(best.r0.3)
      Anova(best.r0.3, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
      plot(best.r0.3, type=c("p","smooth"), col.line=1)
      plot(best.r0.3,
           sqrt(abs(resid(.)))~fitted(.),
           type=c("p","smooth"), col.line=1)
      qqnorm(resid(best.r0.3))
      qqline(resid(best.r0.3))
      
    }
    
    #Create model predictions
    {
      predict.data <- expand.grid(Reproduction_Numeric=seq(from = min(data.stats.2$Reproduction_Numeric), to = max(data.stats.2$Reproduction_Numeric), length = 100),
                                  Evolution=unique(data.stats.2$Evolution),
                                  Environment = unique(data.stats.2$Environment))
      predictions <- predict(newdata=predict.data, best.r0.3, se.fit=T, level=0)
      predict.data$r0.mean <-  (predictions$fit)
      predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
      predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
    }
    
    #Plot model predictions with the raw data on top
    {
      r0.numeric.all <- ggplot(data.stats.2, aes(y=r_OD750_stand, x=Reproduction_Numeric, colour = Environment)) + geom_point() + facet_wrap(~Evolution) +
        geom_line(data=predict.data, inherit.aes = F, mapping = aes(x=Reproduction_Numeric, y = r0.mean, colour = Environment))+
        geom_ribbon(data=predict.data, inherit.aes = F, 
                    mapping = aes(x=Reproduction_Numeric, ymin=r0.lower, ymax=r0.upper, y = r0.mean, fill = Environment),
                    stat = "identity", alpha=0.3)+ theme_bw() +
        theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), axis.text=element_text(size=14), legend.text=element_text(size=14), 
              strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
              axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
        ylab(expression(paste("Relative change in intrinsic rate of increase ", italic("r₀")))) + xlab("Number of sexual reproduction events") + 
        scale_color_manual("Assay environment", values = c("#d8b365", "#5ab4ac")) + 
        scale_fill_manual("Assay environment", values = c("#d8b365", "#5ab4ac")) 
      
      r0.numeric.all
      ggsave(r0.numeric.all, filename = "4_results/Figures/r0_numeric_all_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2, width=67*3)
      
    }
    
    #9.1.2) With reproduction as categorical variable
    {
      #Fit a full interaction model (with assay stress)
      full.r0.4 <- lme(data=data.stats.2, r_OD750_stand~Reproduction*Evolution*Environment, random =  ~ 1|Population_No, method = "ML")
      
      #Look at all possible models
      dredge(full.r0.4)
      
      #Make the best fitting model
      best.r0.4 <- lme(data=data.stats.2, r_OD750_stand~Evolution+Environment*Reproduction, random =  ~ 1|Population_No, method = "REML")
      
      #Look at model output
      summary(best.r0.4)
      Anova(best.r0.4, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
      plot(best.r0.4, type=c("p","smooth"), col.line=1)
      plot(best.r0.4,
           sqrt(abs(resid(.)))~fitted(.),
           type=c("p","smooth"), col.line=1)
      qqnorm(resid(best.r0.4))
      qqline(resid(best.r0.4))
      
    }
    
    #Create model predictions
    {
      predict.data <- expand.grid(Reproduction = unique(data.stats.2$Reproduction),
                                  Evolution=unique(data.stats.2$Evolution),
                                  Environment = unique(data.stats.2$Environment))
      predictions <- predict(newdata=predict.data, best.r0.4, se.fit=T, level=0)
      predict.data$r0.mean <-  (predictions$fit)
      predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
      predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
    }
    
    #Plot model predictions with the raw data on top
    {
      r0.categorical.all <- ggplot(data.stats.2, aes(y=r_OD750_stand, x=Reproduction)) + geom_point() + facet_grid(Environment~Evolution) +
        geom_boxplot(data=predict.data, inherit.aes = F, 
                     mapping = aes(x=Reproduction, ymin=r0.lower, ymax=r0.upper, middle=r0.mean, lower=r0.lower, upper=r0.upper, fill = Reproduction),
                     stat = "identity", alpha=0.3)+ theme_bw() +
        theme(axis.text=element_text(size=14), legend.text=element_text(size=14), 
              strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
              axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
        ylab(expression(paste("Relative change in intrinsic rate of increase ", italic("r₀")))) + xlab("Frequency of sexual reproduction") + 
        scale_color_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) + 
        scale_fill_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) 
      
      r0.categorical.all
      ggsave(r0.categorical.all, filename = "4_results/Figures/r0_categorical_all_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2.5, width=67*3.3)
      
    }
    
  }
  
  #9.2) Next look at K
  {
    
    #9.2.1) With reproduction as numerical variable
    {
      #Fit a full interaction model (with assay stress)
      full.K.3 <- lme(data=data.stats.2, K_OD750_stand~poly(Reproduction_Numeric,2)*Evolution*Environment, random =  ~ 1|Population_No, method = "ML")
      
      #Look at all possible models
      dredge(full.K.3)
      
      #Make the best fitting model
      best.K.3 <- lme(data=data.stats.2, K_OD750_stand~poly(Reproduction_Numeric,2)*Evolution*Environment, random =  ~ 1|Population_No, method = "REML")
      
      #Look at model output
      summary(best.K.3)
      Anova(best.K.3, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
      plot(best.K.3, type=c("p","smooth"), col.line=1)
      plot(best.K.3,
           sqrt(abs(resid(.)))~fitted(.),
           type=c("p","smooth"), col.line=1)
      qqnorm(resid(best.K.3))
      qqline(resid(best.K.3))
      
    }
    
    #Create model predictions
    {
      predict.data <- expand.grid(Reproduction_Numeric=seq(from = min(data.stats.2$Reproduction_Numeric), to = max(data.stats.2$Reproduction_Numeric), length = 100),
                                  Evolution=unique(data.stats.2$Evolution),
                                  Environment = unique(data.stats.2$Environment))
      predictions <- predict(newdata=predict.data, best.K.3, se.fit=T, level=0)
      predict.data$r0.mean <-  (predictions$fit)
      predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
      predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
    }
    
    #Plot model predictions with the raw data on top
    {
      K.numeric.all <- ggplot(data.stats.2, aes(y=K_OD750_stand, x=Reproduction_Numeric, colour = Environment)) + geom_point() + facet_wrap(~Evolution) +
        geom_line(data=predict.data, inherit.aes = F, mapping = aes(x=Reproduction_Numeric, y = r0.mean, colour = Environment))+
        geom_ribbon(data=predict.data, inherit.aes = F, 
                    mapping = aes(x=Reproduction_Numeric, ymin=r0.lower, ymax=r0.upper, y = r0.mean, fill = Environment),
                    stat = "identity", alpha=0.3)+ theme_bw() +
        theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), axis.text=element_text(size=14), legend.text=element_text(size=14), 
              strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
              axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
        ylab(expression(paste("Relative change in equilibrium density ", italic("K")))) + xlab("Number of sexual reproduction events") + 
        scale_color_manual("Assay environment", values = c("#d8b365", "#5ab4ac")) + 
        scale_fill_manual("Assay environment", values = c("#d8b365", "#5ab4ac")) 
      
      K.numeric.all
      ggsave(K.numeric.all, filename = "4_results/Figures/K_numeric_all_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2, width=67*3)
      
    }
    
  }
  
  #9.2.2) With reproduction as categorical variable
  {
    #Fit a full interaction model (with assay stress)
    full.K.4 <- lme(data=data.stats.2, K_OD750_stand~Reproduction*Evolution*Environment, random =  ~ 1|Population_No, method = "ML")
    
    #Look at all possible models
    dredge(full.K.4)
    
    #Make the best fitting model
    best.K.4 <- lme(data=data.stats.2, K_OD750_stand~Evolution*Environment, random =  ~ 1|Population_No, method = "REML")
    
    #Look at model output
    summary(best.K.4)
    Anova(best.K.4, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
    plot(best.K.4, type=c("p","smooth"), col.line=1)
    plot(best.K.4,
         sqrt(abs(resid(.)))~fitted(.),
         type=c("p","smooth"), col.line=1)
    qqnorm(resid(best.K.4))
    qqline(resid(best.K.4))
    
  }
  
  #Create model predictions
  {
    predict.data <- expand.grid(Reproduction = unique(data.stats.2$Reproduction),
                                Evolution=unique(data.stats.2$Evolution),
                                Environment = unique(data.stats.2$Environment))
    predictions <- predict(newdata=predict.data, best.K.4, se.fit=T, level=0)
    predict.data$r0.mean <-  (predictions$fit)
    predict.data$r0.upper <- (predictions$fit + 1.96*predictions$se.fit)
    predict.data$r0.lower <- (predictions$fit - 1.96*predictions$se.fit)
  }
  
  #Plot model predictions with the raw data on top
  {
    K.categorical.all <- ggplot(data.stats.2, aes(y=K_OD750_stand, x=Reproduction)) + geom_point() + facet_grid(Environment~Evolution) +
      geom_boxplot(data=predict.data, inherit.aes = F, 
                   mapping = aes(x=Reproduction, ymin=r0.lower, ymax=r0.upper, middle=r0.mean, lower=r0.lower, upper=r0.upper, fill = Reproduction),
                   stat = "identity", alpha=0.3)+ theme_bw() +
      theme(axis.text=element_text(size=14), legend.text=element_text(size=14), 
            strip.text.x=element_text(size=14), strip.text.y=element_text(size=14), legend.title=element_text(size=16), legend.position = "bottom",
            axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5), legend.box="vertical", legend.key = element_rect(fill = "white")) + 
      ylab(expression(paste("Relative change in equilibrium density ", italic("K")))) + xlab("Frequency of sexual reproduction") + 
      scale_color_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) + 
      scale_fill_manual("Frequency of sexual reproduction", values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) 
    
    K.categorical.all
    ggsave(K.categorical.all, filename = "4_results/Figures/K_categorical_all_adaptation.pdf", device = cairo_pdf, dpi = 500, units = "mm", height=74*2.5, width=67*3.3)
    
  }
  
}

#10) Check correlations to see how trade-offs are affected for r-K and for one environment versus the other
{
  ggplot(data.stats.2, aes(x = r_OD750_stand, y = K_OD750_stand, colour = Reproduction)) + geom_point() + geom_smooth(method="lm") + facet_wrap(Evolution~Environment)
  plot(data.stats.2$r_OD750_stand, data.stats.2$K_OD750_stand)
  
  data.stats.3 <- filter(data.stats.2, Environment == "No salt environment") 
  data.stats.3 <- mutate(data.stats.3, r_OD750_stand_salt = 0, data.stats.3, K_OD750_stand_salt = 0)
  for (i in unique(data.stats.3$ID)){
    temp <- dplyr::filter(data.stats.3, ID==i)
    temp2 <- filter(data.stats.2, Population_No == temp$Population_No, Replicate == temp$Replicate, Environment=="Salt environment")
    data.stats.3[which(data.stats.3$ID==i), "r_OD750_stand_salt"] <- temp2$r_OD750_stand
    data.stats.3[which(data.stats.3$ID==i), "K_OD750_stand_salt"] <-temp2$K_OD750_stand
  }
  
  ggplot(data.stats.3, aes(x = r_OD750_stand, y = r_OD750_stand_salt, colour = Reproduction)) + geom_point() + geom_smooth(method="lm") + facet_grid(~Evolution)
  
}

#11) Compare coefficient of variation between the populations
{
  data.cv <- data.stats.2 %>% group_by(Reproduction, Evolution, Environment) %>% 
    summarize(CV_r = sd(r_OD750_stand)/mean(r_OD750_stand), CV_K = sd(K_OD750_stand)/mean(K_OD750_stand))
  
  ggplot(data.cv, aes(x = Reproduction, y = CV_r)) + geom_point() + facet_grid(Environment~Evolution)
  cor.test(data.cv$CV_r, data.cv$CV_K)
  model <- lm(data=data.cv, CV_K~Environment , na.action = "na.fail")
  dredge(model)
  summary(model)
}