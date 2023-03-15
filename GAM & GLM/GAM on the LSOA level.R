###################################################################################

library("ggplot2")
library("dplyr")
library("mgcv")

rm(list = ls())
gc()
dev.off()

# load data
asb.dataset <- read.csv("ASB.csv")
# sort the arrangement so that its ordered by lsoa and month
asb.dataset <- asb.dataset[order(asb.dataset$lsoa, asb.dataset$month),]

# create unique group id for the lsoa codes in asb.dataset
asb.dataset$lsoa1 <- asb.dataset$lsoa
asb.dataset$number <- as.numeric(as.factor(asb.dataset$lsoa1))

# create a loop
for (i in 1:182) {
	# perform it on first LSOA when i == 1
	if (i == 1) {
	data.filtered <- asb.dataset[asb.dataset$number==i,]
	# run the gam model
	gam.result <- gam(rate~s(mon), data=data.filtered, method="REML")
	# compute summary into a table object
	table.object = summary(gam.result)
	# extract edf value and p-value - this needs to be plotted on all graphs as evidence you made the plots
	table.object$s.table[1]
	table.object$s.table[4]
	
	# create custom graphs for each LSOA
	png(paste0("LSOA_",data.filtered$lsoa[1],".png"), width=956, height=845)
	plot(gam.result, shade = TRUE, shade.col = "azure2", xlab = "Months", ylab = "Effect of Lockdown on ASB", main = bquote("LSOA:"~.(data.filtered$lsoa[1]) ~ " EDF Value =" ~.(table.object$s.table[1]) ~ " p-value =" ~.(table.object$s.table[4])))
	abline(h=0, col="black", lty=2)
	abline(v=11, col="red", lty=2)
	abline(v=15, col="red", lty=2)
	abline(v=19, col="red", lty=2)
	abline(v=24, col="red", lty=2)
	dev.off()
	
	# extract values to determine significance at each time point
	estimates.object <- plot(gam.result, residuals = TRUE)
	estimates.object <- estimates.object[[1]]
	smoothed_estimates.object <- as.data.frame(estimates.object[c("x", "se", "fit")])
	smoothed_estimates.object$lower.limit <- smoothed_estimates.object$fit - smoothed_estimates.object$se
	smoothed_estimates.object$upper.limit <- smoothed_estimates.object$fit + smoothed_estimates.object$se
	
	smoothed_estimates.object$significance <- NA
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit>0] <- 0    # NOT SIGNIFICANT
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit==0 | smoothed_estimates.object$upper.limit==0] <- 0  # NOT SIGNIFICANT
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit>0 & smoothed_estimates.object$upper.limit>0] <- 1    # SIGNIFICANT INCREASE
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit<0] <- -1   # SIGNIFICANT DECREASE
	
	# keep pre-lockdown phase
	PRE = smoothed_estimates.object[(smoothed_estimates.object$x >= 1) & (smoothed_estimates.object$x <= 11),]
	# keep lockdown Phase 1
	LKD1 = smoothed_estimates.object[(smoothed_estimates.object$x >= 11) & (smoothed_estimates.object$x <= 15),]
	# keep post lockdown Phase 1
	LKD1_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 15) & (smoothed_estimates.object$x <= 19),]
	# keep lockdown Phase 2 & 3
	LKD23 = smoothed_estimates.object[(smoothed_estimates.object$x >= 20) & (smoothed_estimates.object$x <= 24),]
	# keep post lockdown
	LKD23_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 24) & (max(smoothed_estimates.object$x)),]
	
	# housekeeping
	#::: create empty data frame to store final significance class outputs
	spatial_data_base = data.frame(lsoa=character(), pre_lkd = integer(), lkd_phase1 = integer(), lkd_post_phase1 = integer(), lkd_phase23=integer(), lkd_postphase=integer())
	spatial_data_base <- spatial_data_base[1:1,]
	row.names(spatial_data_base) <- 1:1
	#::: create empty data frame to store significance class outputs and pass to spatial_data_base object
	#::: the magic will take place in this object
	sig_value = data.frame(Sig.R=character(), Not.Sig = character(), Sig.I = character())
	sig_value<- sig_value[1:1,]
	row.names(sig_value) <- 1:1
	
	# attach lsoa label into spatial_data_base and start populating results
	spatial_data_base$lsoa[1] = data.filtered$lsoa[1]
	
	# determine consist estimates that were sigificant reduction, not significant or significant increase.
	#::: only one can be "TRUE"
	sig_value$Sig.R[1] = all(PRE$significance == -1)
	sig_value$Not.Sig[1] = all(PRE$significance == 0)
	sig_value$Sig.I[1] = all(PRE$significance == 1)
	
	# print the name of the column that matches the "TRUE" statement
	# :::assign significance value to LKD1 column
	spatial_data_base$pre_lkd[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_phase1
	sig_value$Sig.R[1] = all(LKD1$significance == -1)
	sig_value$Not.Sig[1] = all(LKD1$significance == 0)
	sig_value$Sig.I[1] = all(LKD1$significance == 1)
	spatial_data_base$lkd_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_post_phase1
	sig_value$Sig.R[1] = all(LKD1_post$significance == -1)
	sig_value$Not.Sig[1] = all(LKD1_post$significance == 0)
	sig_value$Sig.I[1] = all(LKD1_post$significance == 1)
	spatial_data_base$lkd_post_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_phase23
	sig_value$Sig.R[1] = print(all(LKD23$significance == -1))
	sig_value$Not.Sig[1] = print(all(LKD23$significance == 0))
	sig_value$Sig.I[1] = print(all(LKD23$significance == 1))
	spatial_data_base$lkd_phase23[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_postphase
	sig_value$Sig.R[1] = print(all(LKD23_post$significance == -1))
	sig_value$Not.Sig[1] = print(all(LKD23_post$significance == 0))
	sig_value$Sig.I[1] = print(all(LKD23_post$significance == 1))
	spatial_data_base$lkd_postphase[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	initial_base <- spatial_data_base
	}
	
	# perform it on rest of all LSOAs when i = 2, 3, 4 and so...
	if (i != 1) {
	data.filtered <- asb.dataset[asb.dataset$number==i,]
	# run the gam model
	gam.result <- gam(rate~s(mon), data=data.filtered, method="REML")
	# compute summary into a table object
	table.object = summary(gam.result)
	# extract edf value and p-value - this needs to be plotted on all graphs as evidence you made the plots
	table.object$s.table[1]
	table.object$s.table[4]
			
	# create custom graphs for each LSOA
	png(paste0("LSOA_",data.filtered$lsoa[1],".png"), width=956, height=845)
	plot(gam.result, shade = TRUE, shade.col = "azure2", xlab = "Months", ylab = "Effect of Lockdown on ASB", main = bquote("LSOA:"~.(data.filtered$lsoa[1]) ~ " EDF Value =" ~.(table.object$s.table[1]) ~ " p-value =" ~.(table.object$s.table[4])))
	abline(h=0, col="black", lty=2)
	abline(v=11, col="red", lty=2)
	abline(v=15, col="red", lty=2)
	abline(v=19, col="red", lty=2)
	abline(v=24, col="red", lty=2)
	dev.off()
			
	# extract values to determine significance at each time point
	estimates.object <- plot(gam.result, residuals = TRUE)
	estimates.object <- estimates.object[[1]]
	smoothed_estimates.object <- as.data.frame(estimates.object[c("x", "se", "fit")])
	smoothed_estimates.object$lower.limit <- smoothed_estimates.object$fit - smoothed_estimates.object$se
	smoothed_estimates.object$upper.limit <- smoothed_estimates.object$fit + smoothed_estimates.object$se
			
	smoothed_estimates.object$significance <- NA
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit>0] <- 0    # NOT SIGNIFICANT
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit==0 | smoothed_estimates.object$upper.limit==0] <- 0  # NOT SIGNIFICANT
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit>0 & smoothed_estimates.object$upper.limit>0] <- 1    # SIGNIFICANT INCREASE
	smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit<0] <- -1   # SIGNIFICANT DECREASE
			
	# keep pre-lockdown phase
	PRE = smoothed_estimates.object[(smoothed_estimates.object$x >= 1) & (smoothed_estimates.object$x <= 11),]
	# keep lockdown Phase 1
	LKD1 = smoothed_estimates.object[(smoothed_estimates.object$x >= 11) & (smoothed_estimates.object$x <= 15),]
	# keep post lockdown Phase 1
	LKD1_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 15) & (smoothed_estimates.object$x <= 19),]
	# keep lockdown Phase 2 & 3
	LKD23 = smoothed_estimates.object[(smoothed_estimates.object$x >= 19) & (smoothed_estimates.object$x <= 24),]
	# keep post lockdown
	LKD23_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 24) & (max(smoothed_estimates.object$x)),]
	
	# housekeeping
	#::: create empty data frame to store final significance class outputs
	spatial_data_base = data.frame(lsoa=character(), pre_lkd = integer(), lkd_phase1 = integer(), lkd_post_phase1 = integer(), lkd_phase23=integer(), lkd_postphase=integer())
	spatial_data_base <- spatial_data_base[1:1,]
	row.names(spatial_data_base) <- 1:1
	#::: create empty data frame to store significance class outputs and pass to spatial_data_base object
	#::: the magic will take place in this object
	sig_value = data.frame(Sig.R=character(), Not.Sig = character(), Sig.I = character())
	sig_value<- sig_value[1:1,]
	row.names(sig_value) <- 1:1
	
	# attach lsoa label into spatial_data_base and start populating results
	spatial_data_base$lsoa[1] = data.filtered$lsoa[1]
	
	# determine consist estimates that were sigificant reduction, not significant or significant increase.
	#::: only one can be "TRUE"
	sig_value$Sig.R[1] = all(PRE$significance == -1)
	sig_value$Not.Sig[1] = all(PRE$significance == 0)
	sig_value$Sig.I[1] = all(PRE$significance == 1)
	
	# print the name of the column that matches the "TRUE" statement
	# :::assign significance value to LKD1 column
	spatial_data_base$pre_lkd[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_phase1
	sig_value$Sig.R[1] = all(LKD1$significance == -1)
	sig_value$Not.Sig[1] = all(LKD1$significance == 0)
	sig_value$Sig.I[1] = all(LKD1$significance == 1)
	spatial_data_base$lkd_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_post_phase1
	sig_value$Sig.R[1] = all(LKD1_post$significance == -1)
	sig_value$Not.Sig[1] = all(LKD1_post$significance == 0)
	sig_value$Sig.I[1] = all(LKD1_post$significance == 1)
	spatial_data_base$lkd_post_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_phase23
	sig_value$Sig.R[1] = print(all(LKD23$significance == -1))
	sig_value$Not.Sig[1] = print(all(LKD23$significance == 0))
	sig_value$Sig.I[1] = print(all(LKD23$significance == 1))
	spatial_data_base$lkd_phase23[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
	
	# repeat for lkd_postphase
	sig_value$Sig.R[1] = print(all(LKD23_post$significance == -1))
	sig_value$Not.Sig[1] = print(all(LKD23_post$significance == 0))
	sig_value$Sig.I[1] = print(all(LKD23_post$significance == 1))
	spatial_data_base$lkd_postphase[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
		"Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])

			
	append_base <- spatial_data_base
	initial_base <- rbind(initial_base, append_base)

	}
}



initial_base$pre_lkd[initial_base$pre_lkd == 'Not.Sig'] <- "Not Significant"
initial_base$pre_lkd[initial_base$pre_lkd == 'Sig.I'] <- "Significant Increase"
initial_base$pre_lkd[initial_base$pre_lkd == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_phase1[initial_base$lkd_phase1 == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_phase1[initial_base$lkd_phase1 == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_phase1[initial_base$lkd_phase1 == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_phase23[initial_base$lkd_phase23 == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_phase23[initial_base$lkd_phase23 == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_phase23[initial_base$lkd_phase23 == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_postphase[initial_base$lkd_postphase == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_postphase[initial_base$lkd_postphase == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_postphase[initial_base$lkd_postphase == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_post_phase1[initial_base$lkd_post_phase1 == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_post_phase1[initial_base$lkd_post_phase1 == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_post_phase1[initial_base$lkd_post_phase1 == 'Sig.R'] <- "Significant Decrease"

LSOA_ASB<- merge(data_not,initial_base, by.x="code", by.y="lsoa")

color2 <- c("floralwhite","cadetblue3","coral1")

tm_shape(LSOA_ASB) + tm_borders(lwd = .2)+ tm_polygons(col="pre_lkd", palette=color2, midpoint = 0)+tmap_options(check.and.fix = TRUE) +tm_layout(frame = FALSE)+ 
  tm_compass(type = "8star", position = c("right", "top")) +tm_scale_bar(breaks = c(0, 1, 2), text.size = 1, position = c("right", "bottom"))+ tm_layout( legend.outside=TRUE,legend.text.size = 1, legend.title.size=1.5)




View(asb.dataset)

View(VSA_effect)

VSA_effect$asb <- VSA_effect$rate

VSA_effect$rate <- VSA_effect$vsa

##

VSA_effect_1 <- data.frame(VSA_effect$lsoa, VSA_effect$month.y,VSA_effect$rate,VSA_effect$mon)

VSA_effect$month <- VSA_effect$month.y

# load data
asb.dataset <- read.csv("ASB.csv")

asb.dataset <- Notting_VAS



# sort the arrangement so that its ordered by lsoa and month
asb.dataset <- asb.dataset[order(asb.dataset$lsoa, asb.dataset$month),]

# create unique group id for the lsoa codes in asb.dataset
asb.dataset$lsoa1 <- asb.dataset$lsoa
asb.dataset$number <- as.numeric(as.factor(asb.dataset$lsoa1))

# create a loop
for (i in 1:182) {
  # perform it on first LSOA when i == 1
  if (i == 1) {
    data.filtered <- asb.dataset[asb.dataset$number==i,]
    # run the gam model
    gam.result <- gam(rate~s(mon), data=data.filtered, method="REML")
    # compute summary into a table object
    table.object = summary(gam.result)
    # extract edf value and p-value - this needs to be plotted on all graphs as evidence you made the plots
    table.object$s.table[1]
    table.object$s.table[4]
    
    # create custom graphs for each LSOA
    plot(gam.result, shade = TRUE, shade.col = "azure2", xlab = "Months", ylab = "Effect of Lockdown on ASB", main = bquote("LSOA:"~.(data.filtered$lsoa[1]) ~ " EDF Value =" ~.(table.object$s.table[1]) ~ " p-value =" ~.(table.object$s.table[4])))
    abline(h=0, col="black", lty=2)
    abline(v=11, col="red", lty=2)
    abline(v=15, col="red", lty=2)
    abline(v=19, col="red", lty=2)
    abline(v=24, col="red", lty=2)
    dev.off()
    
    # extract values to determine significance at each time point
    estimates.object <- plot(gam.result, residuals = TRUE)
    estimates.object <- estimates.object[[1]]
    smoothed_estimates.object <- as.data.frame(estimates.object[c("x", "se", "fit")])
    smoothed_estimates.object$lower.limit <- smoothed_estimates.object$fit - smoothed_estimates.object$se
    smoothed_estimates.object$upper.limit <- smoothed_estimates.object$fit + smoothed_estimates.object$se
    
    smoothed_estimates.object$significance <- NA
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit>0] <- 0    # NOT SIGNIFICANT
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit==0 | smoothed_estimates.object$upper.limit==0] <- 0  # NOT SIGNIFICANT
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit>0 & smoothed_estimates.object$upper.limit>0] <- 1    # SIGNIFICANT INCREASE
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit<0] <- -1   # SIGNIFICANT DECREASE
    
    # keep pre-lockdown phase
    PRE = smoothed_estimates.object[(smoothed_estimates.object$x >= 1) & (smoothed_estimates.object$x <= 11),]
    # keep lockdown Phase 1
    LKD1 = smoothed_estimates.object[(smoothed_estimates.object$x >= 11) & (smoothed_estimates.object$x <= 14),]
    # keep post lockdown Phase 1
    LKD1_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 14) & (smoothed_estimates.object$x <= 19),]
    # keep lockdown Phase 2 & 3
    LKD23 = smoothed_estimates.object[(smoothed_estimates.object$x >= 20) & (smoothed_estimates.object$x <= 24),]
    # keep post lockdown
    LKD23_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 24) & (max(smoothed_estimates.object$x)),]
    
    # housekeeping
    #::: create empty data frame to store final significance class outputs
    spatial_data_base = data.frame(lsoa=character(), pre_lkd = integer(), lkd_phase1 = integer(), lkd_post_phase1 = integer(), lkd_phase23=integer(), lkd_postphase=integer())
    spatial_data_base <- spatial_data_base[1:1,]
    row.names(spatial_data_base) <- 1:1
    #::: create empty data frame to store significance class outputs and pass to spatial_data_base object
    #::: the magic will take place in this object
    sig_value = data.frame(Sig.R=character(), Not.Sig = character(), Sig.I = character())
    sig_value<- sig_value[1:1,]
    row.names(sig_value) <- 1:1
    
    # attach lsoa label into spatial_data_base and start populating results
    spatial_data_base$lsoa[1] = data.filtered$lsoa[1]
    
    # determine consist estimates that were sigificant reduction, not significant or significant increase.
    #::: only one can be "TRUE"
    sig_value$Sig.R[1] = all(PRE$significance == -1)
    sig_value$Not.Sig[1] = all(PRE$significance == 0)
    sig_value$Sig.I[1] = all(PRE$significance == 1)
    
    # print the name of the column that matches the "TRUE" statement
    # :::assign significance value to LKD1 column
    spatial_data_base$pre_lkd[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                          "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_phase1
    sig_value$Sig.R[1] = all(LKD1$significance == -1)
    sig_value$Not.Sig[1] = all(LKD1$significance == 0)
    sig_value$Sig.I[1] = all(LKD1$significance == 1)
    spatial_data_base$lkd_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                             "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_post_phase1
    sig_value$Sig.R[1] = all(LKD1_post$significance == -1)
    sig_value$Not.Sig[1] = all(LKD1_post$significance == 0)
    sig_value$Sig.I[1] = all(LKD1_post$significance == 1)
    spatial_data_base$lkd_post_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                                  "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_phase23
    sig_value$Sig.R[1] = print(all(LKD23$significance == -1))
    sig_value$Not.Sig[1] = print(all(LKD23$significance == 0))
    sig_value$Sig.I[1] = print(all(LKD23$significance == 1))
    spatial_data_base$lkd_phase23[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                              "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_postphase
    sig_value$Sig.R[1] = print(all(LKD23_post$significance == -1))
    sig_value$Not.Sig[1] = print(all(LKD23_post$significance == 0))
    sig_value$Sig.I[1] = print(all(LKD23_post$significance == 1))
    spatial_data_base$lkd_postphase[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                                "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    initial_base <- spatial_data_base
  }
  
  # perform it on rest of all LSOAs when i = 2, 3, 4 and so...
  if (i != 1) {
    data.filtered <- asb.dataset[asb.dataset$number==i,]
    # run the gam model
    gam.result <- gam(rate~s(mon), data=data.filtered, method="REML")
    # compute summary into a table object
    table.object = summary(gam.result)
    # extract edf value and p-value - this needs to be plotted on all graphs as evidence you made the plots
    table.object$s.table[1]
    table.object$s.table[4]
    
    # create custom graphs for each LSOA
    png(paste0("LSOA_",data.filtered$lsoa[1],".png"), width=956, height=845)
    plot(gam.result, shade = TRUE, shade.col = "azure2", xlab = "Months", ylab = "Effect of Lockdown on ASB", main = bquote("LSOA:"~.(data.filtered$lsoa[1]) ~ " EDF Value =" ~.(table.object$s.table[1]) ~ " p-value =" ~.(table.object$s.table[4])))
    abline(h=0, col="black", lty=2)
    abline(v=11, col="red", lty=2)
    abline(v=15, col="red", lty=2)
    abline(v=19, col="red", lty=2)
    abline(v=24, col="red", lty=2)
    dev.off()
    
    # extract values to determine significance at each time point
    estimates.object <- plot(gam.result, residuals = TRUE)
    estimates.object <- estimates.object[[1]]
    smoothed_estimates.object <- as.data.frame(estimates.object[c("x", "se", "fit")])
    smoothed_estimates.object$lower.limit <- smoothed_estimates.object$fit - smoothed_estimates.object$se
    smoothed_estimates.object$upper.limit <- smoothed_estimates.object$fit + smoothed_estimates.object$se
    
    smoothed_estimates.object$significance <- NA
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit>0] <- 0    # NOT SIGNIFICANT
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit==0 | smoothed_estimates.object$upper.limit==0] <- 0  # NOT SIGNIFICANT
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit>0 & smoothed_estimates.object$upper.limit>0] <- 1    # SIGNIFICANT INCREASE
    smoothed_estimates.object$significance[smoothed_estimates.object$lower.limit<0 & smoothed_estimates.object$upper.limit<0] <- -1   # SIGNIFICANT DECREASE
    
    # keep pre-lockdown phase
    PRE = smoothed_estimates.object[(smoothed_estimates.object$x >= 1) & (smoothed_estimates.object$x <= 11),]
    # keep lockdown Phase 1
    LKD1 = smoothed_estimates.object[(smoothed_estimates.object$x >= 11) & (smoothed_estimates.object$x <= 15),]
    # keep post lockdown Phase 1
    LKD1_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 15) & (smoothed_estimates.object$x <= 19),]
    # keep lockdown Phase 2 & 3
    LKD23 = smoothed_estimates.object[(smoothed_estimates.object$x >= 19) & (smoothed_estimates.object$x <= 24),]
    # keep post lockdown
    LKD23_post = smoothed_estimates.object[(smoothed_estimates.object$x >= 24) & (max(smoothed_estimates.object$x)),]
    
    # housekeeping
    #::: create empty data frame to store final significance class outputs
    spatial_data_base = data.frame(lsoa=character(), pre_lkd = integer(), lkd_phase1 = integer(), lkd_post_phase1 = integer(), lkd_phase23=integer(), lkd_postphase=integer())
    spatial_data_base <- spatial_data_base[1:1,]
    row.names(spatial_data_base) <- 1:1
    #::: create empty data frame to store significance class outputs and pass to spatial_data_base object
    #::: the magic will take place in this object
    sig_value = data.frame(Sig.R=character(), Not.Sig = character(), Sig.I = character())
    sig_value<- sig_value[1:1,]
    row.names(sig_value) <- 1:1
    
    # attach lsoa label into spatial_data_base and start populating results
    spatial_data_base$lsoa[1] = data.filtered$lsoa[1]
    
    # determine consist estimates that were sigificant reduction, not significant or significant increase.
    #::: only one can be "TRUE"
    sig_value$Sig.R[1] = all(PRE$significance == -1)
    sig_value$Not.Sig[1] = all(PRE$significance == 0)
    sig_value$Sig.I[1] = all(PRE$significance == 1)
    
    # print the name of the column that matches the "TRUE" statement
    # :::assign significance value to LKD1 column
    spatial_data_base$pre_lkd[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                          "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_phase1
    sig_value$Sig.R[1] = all(LKD1$significance == -1)
    sig_value$Not.Sig[1] = all(LKD1$significance == 0)
    sig_value$Sig.I[1] = all(LKD1$significance == 1)
    spatial_data_base$lkd_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                             "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_post_phase1
    sig_value$Sig.R[1] = all(LKD1_post$significance == -1)
    sig_value$Not.Sig[1] = all(LKD1_post$significance == 0)
    sig_value$Sig.I[1] = all(LKD1_post$significance == 1)
    spatial_data_base$lkd_post_phase1[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                                  "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_phase23
    sig_value$Sig.R[1] = print(all(LKD23$significance == -1))
    sig_value$Not.Sig[1] = print(all(LKD23$significance == 0))
    sig_value$Sig.I[1] = print(all(LKD23$significance == 1))
    spatial_data_base$lkd_phase23[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                              "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    # repeat for lkd_postphase
    sig_value$Sig.R[1] = print(all(LKD23_post$significance == -1))
    sig_value$Not.Sig[1] = print(all(LKD23_post$significance == 0))
    sig_value$Sig.I[1] = print(all(LKD23_post$significance == 1))
    spatial_data_base$lkd_postphase[1] = ifelse(sig_value[,1]=="FALSE" & sig_value[,2]=="FALSE" & sig_value[,3]=="FALSE", 
                                                "Not.Sig", names(sig_value)[max.col(sig_value == "TRUE")])
    
    
    append_base <- spatial_data_base
    initial_base <- rbind(initial_base, append_base)
    
  }
}


initial_base$pre_lkd[initial_base$pre_lkd == 'Not.Sig'] <- "Not Significant"
initial_base$pre_lkd[initial_base$pre_lkd == 'Sig.I'] <- "Significant Increase"
initial_base$pre_lkd[initial_base$pre_lkd == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_phase1[initial_base$lkd_phase1 == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_phase1[initial_base$lkd_phase1 == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_phase1[initial_base$lkd_phase1 == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_phase23[initial_base$lkd_phase23 == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_phase23[initial_base$lkd_phase23 == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_phase23[initial_base$lkd_phase23 == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_postphase[initial_base$lkd_postphase == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_postphase[initial_base$lkd_postphase == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_postphase[initial_base$lkd_postphase == 'Sig.R'] <- "Significant Decrease"

initial_base$lkd_post_phase1[initial_base$lkd_post_phase1 == 'Not.Sig'] <- "Not Significant"
initial_base$lkd_post_phase1[initial_base$lkd_post_phase1 == 'Sig.I'] <- "Significant Increase"
initial_base$lkd_post_phase1[initial_base$lkd_post_phase1 == 'Sig.R'] <- "Significant Decrease"

LSOA_ASB<- merge(data_not,initial_base, by.x="code", by.y="lsoa")

View(initial_base)






