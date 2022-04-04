#' # Generate in silico IMC measurement
#' 	Author: Attila Gabor
#'	Email: attila.gabor@iuni-hidelberg.de
#'  	Date: 17.01.2022
#'
#' main functions in this script:
#' * import a cell layout file (segmented image)
#' * build a model using a configurations file
#' * simulate the reaction diffusion system
#' * plot/show the data
#' * export measurements

# Geometry from Jovan: Tissues.zip
# 

segmentationFile = "./IMC/segmentation/Tissues/Tissue1_grid_preproc.csv"

# config_v7: smaller (20 nodes) case studies
reacDiff_configFile = "./IMC/insilico_data_generation/general_config_v8_correctIC.R"
outputFolder = "./data/IMC/insilico_data_generation/test_folder"
exp_id = "manual_test"


insilico_IMC_generator = function(segmentationFile, reacDiff_configFile, exp_id, outputFolder){
	
    library(RColorBrewer)
    library(dplyr)
    library(ggplot2)
    library(CellNOptR)
    library(igraph)
    source("./IMC/insilico_data_generation/IMC_model.R")
    
	# check inputs
	if(!file.exists(segmentationFile)) stop(paste("segmentation file: ", segmentationFile," does not exist."))
	if(!file.exists(reacDiff_configFile)) stop(paste("segmentation file: ", reacDiff_configFile," does not exist."))
	if(!dir.exists(outputFolder)) dir.create(outputFolder,recursive = T)
	
	
	# read the geometry of the image (segmented image). in the file we store the cell-types on each grid point. 
	cells = read.delim(segmentationFile,sep = ",",colClasses=c("character"),header = F)
	if(any(cells=="")) cells[cells==""] = "0"
	
	# An object from the IMC class is initialised using the configuration file.
	configFile = reacDiff_configFile
	M = IMC$new(cells,configFile=configFile,randomizeIC=TRUE)
	
	# Export the connectivity
	celltypes <- M$celltypes

	
	get_interaction_for_cell_type <- function(cell_type="CT1",inited_IMC_model ){
		
		cell_types <- lapply(inited_IMC_model$nodeTable,function(x)x$type) %>% unlist()
		
		node_indx = which(cell_types==cell_type)[[1]]
		node <- inited_IMC_model$nodeTable[[node_indx]]
		
		inter_mat <- node$prodCoeff %>% as_tibble() %>% add_column(target = rownames(node$prodCoeff),.before = 1)
		
		interactions <- inter_mat %>% gather(from,Vmax,-target) %>% 
			mutate(cell_type = cell_type )
	}
	
	intra_interactions <- list("CT1","CT2","CT3","CT4") %>% 
		lapply(.,get_interaction_for_cell_type,inited_IMC_model=M)  %>%
		bind_rows()
	write_csv(intra_interactions,file = file.path(outputFolder,paste0(exp_id,"_cell_type_interactions.csv")))
	
	
	# start the simulation and measure how much time it takes:
	simulation_max_steps = 50
	cat("Simulating steady state:\n")
	t0 = Sys.time()
	M$simulate(Nstep = simulation_max_steps,dt = 0.2)
	elapsed = Sys.time() - t0
	print(elapsed)
	
	# plot the geometry
	pdf(file.path(outputFolder,paste0(exp_id,"_geometry.pdf")))
	M$plotCells()
	dev.off()
	
	#M$recordGIF(folderName = "./data/IMC/insilico_data_generation/random_IMC1_healthy_v4/gif")
	# visual check that steady state achieved: 
	if(FALSE){
		M$plotTotal(ispecies = 1:5)
		M$plotTotal(ispecies = 6:10)
		M$plotTotal(ispecies = c(4,9,16))
		M$plotTotal(ispecies = c(1:29)) + scale_y_log10()
		gg = M$plotTotal(ispecies = 1:29)
		
		sim_table %>% filter(row<20, col < 20) %>%
			ggplot() + geom_point(aes(time,L4),col="red") + 
			 geom_point(aes(time,R4),,col="blue")+
			 geom_point(aes(time,X23),col="black") + facet_wrap(~type)
		
		sim_table %>% filter(time == 50) %>%
			ggplot() + geom_point(aes(L4,R4),col="red") + 
			 facet_wrap(~type,scales = "free")
		sim_table %>% filter(time == 50) %>%
			ggplot() + geom_point(aes(R4,X23),col="red") + 
			facet_wrap(~type,scales = "free")
		
		
		gg + facet_wrap(~variable,scales = "free")
		M$plotsimStates(itime = 50, ispecies = 1)
		M$plotsimStates(itime = 50, ispecies = 2)
		M$plotsimStates(itime = 50, ispecies = 3)
		M$plotsimStates(itime = 50, ispecies = 4)
		M$plotsimStates(itime = 50, ispecies = 5)
		M$plotsimStates(itime = 10, ispecies = 4)
		M$plotsimStates(itime = 10, ispecies = 9)
		M$plotsimStates(itime = 10, ispecies = 16)
		
		
		M$plotsimStates(itime = 50, ispecies = 8)
		
		M$plotsimStates(itime = 10, ispecies = 19)
		M$plotsimStates(itime = 1, ispecies = 18)
		M$states
		M$plotsimStates(itime = 10, ispecies = 11)
		M$plotsimStates(itime = 2, ispecies = 5)
		M$plotsimStates(itime = 2, ispecies = 5)
		M$plotsimStates(itime = 3, ispecies = 5)
		M$plotsimStates(itime = 1, ispecies = 18)
		
		M$plotsimStates(itime = 3, ispecies = 10)
		M$plotsimStates(itime = 5, ispecies = 29)
		M$plotsimStates(itime = 10, ispecies = 2)
		
	}
	
	# Extract and then plot the final distribution of the species:
	M_simTable = M$toTable() %>% gather(species, value,-row,-col,-time)
	
	M_simTable <- M_simTable %>% as_tibble() %>% 
		dplyr::group_by(species) %>% 
		dplyr::mutate(rel_value = value/max(value))
	
	
	myPalette <- colorRampPalette(RColorBrewer:::brewer.pal(9, "BuPu"))
	sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0, 1))
	gg = dplyr::filter(M_simTable,time==simulation_max_steps) %>%
		select(-value) %>%
		tidyr::spread(species,rel_value) %>%
		tidyr::gather(species,rel_value,-1:-3) %>%
		ggplot() + 
		geom_tile(aes(x=col,y=row,fill=rel_value)) +
		facet_wrap(~species) + sc + theme_bw( ) +
		theme(axis.text.x = element_blank(), 
			  axis.text.y = element_blank(), 
			  axis.ticks = element_blank())
	
	dplyr::filter(M_simTable,time==50) %>%
		ggplot() + 
		geom_tile(aes(x=col,y=row,fill=value)) +
		facet_wrap(~species) + theme_bw( ) +
		theme(axis.text.x = element_blank(), 
			  axis.text.y = element_blank(), 
			  axis.ticks = element_blank())
	
	print(gg)
	
	ggsave(file.path(outputFolder,paste0(exp_id,"_simulatedDsitributions_rel.png")),plot = gg, width = 8,height = 8)
	
	#' save the model to a file
	saveRDS(M,file=file.path(outputFolder,paste0(exp_id,"_IMCmodelObject.RDS")))
	
	#' extract the species, that we want to use later in the analysis. 
	M_measTable = dplyr::filter(M_simTable,time==simulation_max_steps)
	
	
	#' export data to file
	write.csv(M_measTable,file = file.path(outputFolder,paste0(exp_id,"_simulatedData.csv")))
	
	#' export data to file, where we also record the neighbour cells
	#M_measTable_neighbr = M$exportWithNeighbours()
	
	#write.csv(M_measTable_neighbr,file =  file.path(outputFolder,paste0(exp_id,"_simulatedData_with_neighbours.csv")))
	
	
	#' export for Jovan
	# export each pixel
	M_measTable_wide = M_measTable %>% select(row,col,species, value) %>%
		tidyr::spread(key = species ,value="value")
	M_measTable_wide$row = as.numeric(M_measTable_wide$row)
	M_measTable_wide$col = as.numeric(M_measTable_wide$col)
	write_csv(M_measTable_wide,file = file.path(outputFolder,paste0(exp_id,"_position_expression.csv")))
	
	# export only pixel/nodes occupied by a real cell
	M_measTable_wide = M_measTable %>% select(row,col,species, value) %>% 
		tidyr::spread(key = species ,value="value")
	M_measTable_wide$row = as.numeric(M_measTable_wide$row)
	M_measTable_wide$col = as.numeric(M_measTable_wide$col)
	
	type_table = do.call("rbind.data.frame",
						 lapply(M$nodeTable, function(x){
						 	c(x$coord, x$type)
						 }))
	colnames(type_table) <- c("row","col","type")
	type_table <- type_table %>% mutate(row=as.numeric(row),col=as.numeric(col))
	M_measTable_wide = full_join(M_measTable_wide,type_table, by = c("row", "col"))
	
	M_measTable_wide = M_measTable_wide[M_measTable_wide$type != "0",]
	write.csv(M_measTable_wide,file = file.path(outputFolder,paste0(exp_id,"_position_expression_real_cells.csv")),row.names = F)
	
	
	heatPlot <- M_measTable_wide %>%
		select(-row,-col,-type) %>% pheatmap::pheatmap(show_rownames = FALSE)
	
	heatPlot$gtable %>% ggsave(plot = .,filename = file.path(outputFolder,paste0(exp_id,"_heatmap_expression_differences.pdf")),width = 10, height = 10)
	
}



## process all image in the segmentation folder

# segmentation files were generated by Jovan/Denis
segmentationFile_list = list.files("./IMC/segmentation/Tissues",pattern = ".*preproc.csv",full.names = T)

reacDiff_configFile = "./IMC/insilico_data_generation/general_config_v8_correctIC.R"
outputFolder = "./data/IMC/insilico_data_generation/TIGER_general_network_v8_IC_small"



### Process paralllel

library(foreach)

parallel::detectCores()
# n.cores <- parallel::detectCores() - 1
n.cores = 2
#create the cluster
my.cluster <- parallel::makeCluster(
	n.cores, 
	type = "PSOCK",
	setup_strategy = "sequential"
)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

foreach::getDoParWorkers()
ifile  = 2

x <- foreach(ifile = 1:length(segmentationFile_list)) %dopar% {
	segmentationFile = segmentationFile_list[[ifile]]
	
	outputFolder_i = paste0(outputFolder,"_IMC",ifile)
	
	exp_id = "small_network"
	insilico_IMC_generator(segmentationFile, reacDiff_configFile, exp_id, outputFolder_i)
	TRUE	
}


parallel::stopCluster(cl = my.cluster)


for(ifile in  1:length(segmentationFile_list)){
	segmentationFile = segmentationFile_list[[ifile]]
	
	outputFolder_i = paste0(outputFolder,"_IMC",ifile)
	
	exp_id = "small_network"
	insilico_IMC_generator(segmentationFile, reacDiff_configFile, exp_id, outputFolder_i)
	TRUE	
}
