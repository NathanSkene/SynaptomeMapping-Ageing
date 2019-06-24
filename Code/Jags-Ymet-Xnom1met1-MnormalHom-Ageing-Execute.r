# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
a.n <- function(x){return(as.numeric(as.character(x)))}
#------------------------------------------------------------------------------- 
library(reshape)


##########################################################################################
# SET PARAMETERS:
##########################################################################################
# << THE BELOW LINE WILL NEED EDITING FOR ANYONE ELSE RUNNING THIS CODE ON WINDOWS >>
baseDir = "~/ownCloud/Synaptome Map Statistics/Ageing/"
setwd(baseDir)
graphFileType = "pdf" 

##########################################################################################
# FOR EACH DATA FILE:
##########################################################################################

# << EDIT BELOW LINE TO CONTAIN NAMES OF THE DATA FILES >>
for(file_tag in c("PSD95_Intensity","PSD95_Density","PSD95_Size","SAP102_Intensity","SAP102_Density","SAP102_Size")){
	
	##########################################################################################
	# SETUP FOLDERS TO STORE THE RESULTING DATA
	##########################################################################################
	variable_type = gsub(".*_","",file_tag)[1]
	synapse_data=data.frame()
	colony = gsub("_.*","",file_tag)
	folder_ref = gsub(":","_",sprintf("%s_%s",paste(file_tag,collapse="-"),Sys.time()))
	analysis_dir = sprintf("%sOutput/%s",baseDir,folder_ref)
	dir.create(file.path(sprintf("%sOutput",baseDir), folder_ref))
	dir.create(file.path(analysis_dir, "Figs_DiagMCMC"))
	dir.create(file.path(analysis_dir, "Figs_PlotMCMC"))
	dir.create(file.path(analysis_dir, "Figs_RegionPlots"))
	dir.create(file.path(analysis_dir, "mcmcChains"))
	dir.create(file.path(analysis_dir, "ResultsTables"))
	
	##########################################################################################
	# STEP 1: LOAD AND PREPARE THE DATA
	##########################################################################################
	
	FULL_RES = data.frame()
	raw_synapse_data = read.csv(sprintf("Data/%s.csv",file_tag))
	raw_synapse_data = t(raw_synapse_data)
	colnames(raw_synapse_data)=raw_synapse_data[1,]
	raw_synapse_data = raw_synapse_data[-1,]
	raw_synapse_data = cbind(mID=rownames(raw_synapse_data),raw_synapse_data)
	syn_dat = melt(data.frame(raw_synapse_data),id.vars=c("mID","Age"))
	syn_dat$value = as.numeric(as.character(syn_dat$value))
	colnames(syn_dat)[colnames(syn_dat)=="variable"]="subregion"
	syn_dat$subregion = gsub("main_","",syn_dat$subregion)
	syn_dat$region = gsub("_.*","",syn_dat$subregion)
	syn_dat = syn_dat[!is.na(syn_dat$value),]
	syn_dat$value = as.numeric(syn_dat$value)
	syn_dat$Age = a.n(syn_dat$Age)

	########################################################################################################
	# STEP 2: PREPARE ARGUMENTS AND SCRIPTS FOR THE MODELLING
	########################################################################################################
	yName="value" 
	xNomName="subregion" 
	xMetName="Age"             
	graphFileType = "pdf"
	source("Code/Jags-Ymet-Xnom1met1-MnormalHom-Ageing.r")	

	########################################################################################################
	# STEP 3: FOR EACH FULL REGION (e.g. hippocampus/cortex/striatum), SEQUENTIALLY PERFORM THE ANALYSIS
	########################################################################################################
	for(region in unique(syn_dat$region)){
		######################################################################################
		### SEPERATE THE DATA RELEVANT TO THAT REGION ########################################
		######################################################################################	
	    myDataFrame = syn_dat[syn_dat$region==region,]
	    myDataFrame$subregion = as.factor(as.character(myDataFrame$subregion))
	   
	    
	    fileNameRoot = file_tag      
	        
	    ##########################################################################################
	    # LOAD AND RUN THE BAYESIAN MODELLING FOR THIS REGION
	    ##########################################################################################        
	    
	    # Generate the MCMC chain:
	    mcmc_fileNameRoot = sprintf("%s/mcmcChains/mcmc_%s_%s_",analysis_dir,region,variable_type) 
	    mcmcCoda = genMCMC( datFrm=myDataFrame , 
	                        yName=yName , xNomName=xNomName , xMetName=xMetName ,
	                        numSavedSteps=11000 , thinSteps=10 , saveName=mcmc_fileNameRoot )
	                        # <-- set numSavedSteps=11000
	    save(mcmcCoda,file=sprintf("%s/mcmcChains/mcmcCoda_%s_%s.RData",analysis_dir,region,variable_type))
	
	    
	    ##########################################################################################
	    # COMMENTED SECTIONS ARE FOR CHECKING THE MODELLING OF PARAMETERS WORKED WELL
	    ##########################################################################################    
	    DiagMCMC_fileNameRoot = sprintf("%s/Figs_DiagMCMC/%s_%s",analysis_dir,region,variable_type)
	    # Display diagnostics of chain, for specified parameters:
	    parameterNames = varnames(mcmcCoda) 
	    show( parameterNames ) # show all parameter names, for reference
	    for ( parName in parameterNames ) {
	        diagMCMC( codaObject=mcmcCoda , parName=parName , 
	                saveName=DiagMCMC_fileNameRoot , saveType=graphFileType )
	    }
	    graphics.off()
	    #------------------------------------------------------------------------------- 
	    # Get summary statistics of chain:
	    summaryInfo = smryMCMC( mcmcCoda , datFrm=myDataFrame , xNomName=xNomName , 
	                            xMetName=xMetName , #contrasts=contrasts , 
	                            saveName=fileNameRoot )
	    show(summaryInfo)
	    # Display posterior information:
		PlotMCMC_fileNameRoot = sprintf("%s/Figs_PlotMCMC/%s_%s",analysis_dir,region,variable_type)
	    plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xNomName=xNomName , 
	              xMetName=xMetName , #contrasts=contrasts , 
	              saveName=PlotMCMC_fileNameRoot , saveType=graphFileType )
	    graphics.off()
	    #------------------------------------------------------------------------------- 
	    
	    ##########################################################################################
	    # HDI ANALYSIS FOR THIS REGION
	    ##########################################################################################    
	    aMet_pLessZero = 1-sum(as.matrix(mcmcCoda)[,"aMet"]<0)/length(as.matrix(mcmcCoda)[,"aMet"])
	    aMet_pGtZero = 1-sum(as.matrix(mcmcCoda)[,"aMet"]>0)/length(as.matrix(mcmcCoda)[,"aMet"])
	    aMet_mean = mean(as.matrix(mcmcCoda)[,"aMet"])
	    region_dat = data.frame(region=region,ageing_mean_monthly_change=aMet_mean,aMet_pLessZero=aMet_pLessZero,aMet_pGtZero=aMet_pGtZero)
	    if(dim(FULL_RES)[1]==0){
	        FULL_RES = region_dat
	    }else{
	        FULL_RES = rbind(FULL_RES,region_dat)
	    }
	    graphics.off() # This closes all of R's graphics windows.
	}
	write.csv(FULL_RES,file=sprintf("%s/ResultsTables/FULL_RES_%s.csv",analysis_dir,variable_type))
}