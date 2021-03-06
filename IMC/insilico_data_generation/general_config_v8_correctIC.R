# CONFIG file
# 	Author: Attila Gabor
#	Email: gabor.attila87@gmail.com
#  	Date: 25.01.2022

# 1) Basic model structure
# we generate a general random network, with given number of input nodes (receptors) and 
# intermediate signaling nodes. The network structure is randomly generated using 
# a tree-gnerating algorithm from igraph (Rpackage) and some heuristics, which adds
# connections between the branches of the tree. 
#
# 2) TFs: 
# terminal nodes (nodes without outgoing edges) are considered as transcription factors. 
# Some of these transcription factors are randomly selected and connected to ligand
# productions. We add these TF-ligand interactions to the network. 
#
# 3) ligand - receptro interactions: 
# we assume receptor specificity and ligands (LX) are connected to therr respective
# receptors (RX).
#
# 4) cell-type specific networks
# we generate four cell types. 
# we want cell-type specific interactions, but, we also want consistency between different
# cell types: this means that if protein X1 and X3 interacted in cell type A then 
# they can also interact within cell type B, if both X1 and X3 are expressed. 
#
# Therefore, each cell type are generated by eliminating some receptors, 
# ligands and signaling nodes from the general network. The selection of the removed nodes 
# are random. Interactions of removed nodes are also removed. 

# generate a random network --------------------------------------------------
# using CellNOpt network generator:
# - 100 nodes, 5 stimuli, 95 measured
# - contains only activation for simplicity
# this serves as a general model. Subsets of this will be considered in each cell type. 


# Kinetics are linear -> avoide saturation

source("./SIF/util_functions_insilico_generation.R")

model_output_folder = "./data/IMC/insilico_data_generation/TIGER_general_network_v8_IC_small_models"
dir.create(model_output_folder,recursive = TRUE,showWarnings = FALSE)
visualize = TRUE
### ------  in-silico network generation parameters ---- 

set.seed(1237)
N_nodes = 30
p_neg = 0
N_measured = 25
N_stimulated = 5
N_inhibited = 0
cross_talk_ratio = 0.7
barabasi_power= 1
barabasi_m = 1
n_cell_types = 4

# reactions coefficients:
# change the parameters of the random generator
generate_prodCoeff <- function(N){
	runif(N, min = 0.5, max = 1.5)
}

# Initial condition parameters
# applied only to signaling nodes;
# receptors and ligands are initally set to 0. 
generate_initial_condition <- function(N){
	runif(N, min = 1, max = 2)
}
init_cond_receptors = 0
init_cond_ligands = 0

# diffusion coefficients:
generate_diffusion_coeffs <- function(N){
	runif(N,0.8,1)	
}

## ---- Basic model structure ----
syn1 <- create_model(N_nodes=N_nodes,
					 N_stimulated=N_stimulated,
					 N_inhibited=N_inhibited,
					 N_measured=N_measured,
					 p_neg=p_neg,
					 barabasi_power = barabasi_power,
					 barabasi_m = barabasi_m,
					 cross_talk_ratio=cross_talk_ratio)

syn1$model$namesSpecies %>% length()
if(visualize) CellNOptR::plotModel(syn1$model)

# receptors: nodes that do not have incoming edges. They are stimuli in CellNopt.
receptors = colnames(syn1$cnolist@stimuli)
names(receptors) = paste("R",1:length(receptors),sep = "")


# select TFs from the terminal nodes
intMat = syn1$model$interMat
intMat[intMat==-1] = 0
n_inputs = intMat %>% t() %>% colSums()

intMat = syn1$model$interMat
intMat[intMat==1] = 0
n_outgoing_edges = intMat %>% t() %>% colSums()
TFs = n_outgoing_edges[n_outgoing_edges==0] %>% names()


# introduce ligands. They are produced based on TF activity. 
ligands = paste("L",1:length(receptors),sep = "")

lig_production = list()

# ---- Ligand production ----
# select randomly the TFs regulating the production of each ligand
for(ilig in 1:length(ligands)){
	n_TF_partner = round(runif(1,1,3))
	partners = sample(TFs,n_TF_partner)
	
	lig_production[[ilig]] = tibble(from = partners,sign = 1, to = ligands[[ilig]])
}

lig_network = do.call("rbind",lig_production)
TFs <- lig_network$from %>% unique()


# ---- Ligand-receptor signaling ----
# each ligand is specific to it's receptor
lig_rec_interactions <- tibble(from = ligands,sign = 1, to = names(receptors))


# working with the network representation
tmp_file <- tempfile(fileext = ".sif")
CellNOptR::writeSIF(syn1$model, tmp_file)

network_tbl <- read_tsv(tmp_file,col_names = FALSE,col_types = "cdc")
colnames(network_tbl) <- c("from","sign","to")


# rename receptors in the network
network <- network_tbl %>% 
	dplyr::rowwise() %>%
	# rename receptors
	dplyr::mutate(from = ifelse(from %in% receptors, names(receptors)[which(from == receptors)],from)) %>%
	# add TF-ligand interactions
	rbind(lig_network) 
# %>%
# 	ungroup() %>%
# 	mutate(keep = !to %in% TFs) %>%
# 	filter(keep) %>% select(-keep)

if(visualize) {
	network %>% write_tsv(tmp_file,col_names = FALSE)
	global_network <- CellNOptR::readSIF(tmp_file)
	pdf(file = file.path(model_output_folder,"global_network.pdf"),width = 10,height = 10)
	CellNOptR::plotModel(global_network)
	dev.off()
}


# adding ligand- receptor interactions would mess up the cellnopt plotting, because 
# then the receptors are no longer located on the top of the figure and ligands
# at the bottom. 
network <- network %>%
	# add lig-rec signaling
	rbind(lig_rec_interactions)


# summarize the nodes: 
ligands
receptors = names(receptors)
signaling_nodes <- c(network$from,network$to) %>% unique() 
signaling_nodes <- signaling_nodes %>% grep("X[0-9]+",x = ., value = TRUE)


# ---- generate cell-type specific models -----
# - remove some receptors in each cell type
# - remove ligand production in each cell type
# - remove some signaling nodes and their interactions in each cell type



ct_parameters <- vector("list",n_cell_types)
ct_network <- vector("list",n_cell_types)

for (iCT  in 1:(n_cell_types)){
	
	# random number of remaining nodes per category:
	ct_parameters[[iCT]]$n_receptors = round(runif(1,2,3))
	ct_parameters[[iCT]]$n_signaling_nodes = round(0.9*length(signaling_nodes))
	ct_parameters[[iCT]]$n_ligands = round(runif(1,2,3))
	
	# randomly selected nodes to be removed
	ct_parameters[[iCT]]$KO_receptors = sample(receptors,length(receptors)-ct_parameters[[iCT]]$n_receptors)
	ct_parameters[[iCT]]$KO_ligands = sample(ligands,length(ligands)-ct_parameters[[iCT]]$n_ligands)
	ct_parameters[[iCT]]$KO_signaling = sample(signaling_nodes,length(signaling_nodes)-ct_parameters[[iCT]]$n_signaling_nodes)
	
	KO_nodes = c(ct_parameters[[iCT]]$KO_receptors,
				 ct_parameters[[iCT]]$KO_ligands,
				 ct_parameters[[iCT]]$KO_signaling)
	
	# remove nodes
	ct_network[[iCT]] <- network %>%
		ungroup() %>% 
		dplyr::filter(!from %in% KO_nodes ) %>%
		dplyr::filter(!to   %in% KO_nodes )
	
	out_file <- file.path(model_output_folder,paste0("CT",iCT,".sif"))
	ct_network[[iCT]] %>% filter(!grepl("^R",to)) %>% 
		write_tsv(.,
									file = out_file,
									col_names = FALSE)
	
	ct_network_cno <- CellNOptR::readSIF(out_file)
	ct_network_cno %>% write_rds(.,file = file.path(model_output_folder,paste0("cellnoptObj_CT",iCT,".rds")))
	if(visualize){
		pdf(file = file.path(model_output_folder,paste0("CT",iCT,"_network.pdf")),width = 10,height = 10)
		CellNOptR::plotModel(ct_network_cno)
		dev.off()
	}
}

# we should check that each receptor and ligand appears in at least one network
n_KO_receptors = lapply(ct_parameters, function(ct)ct$KO_receptors) %>% unlist() %>% table()
n_KO_ligands= lapply(ct_parameters, function(ct)ct$KO_ligands) %>% unlist()%>% table()
n_KO_signaling= lapply(ct_parameters, function(ct)ct$KO_signaling) %>% unlist()%>% table()

if(any(n_KO_receptors==n_cell_types)){
	stop("Error, at least one receptor was KO-d in all cell types.")
}

if(any(n_KO_ligands==n_cell_types)){
	stop("Error, at least one ligand was KO-d in all cell types.")
}

if(any(n_KO_signaling==n_cell_types)){
	stop("Error, at least one signaling node was KO-d in all cell types.")
}

# End of network _structure_ generation for cell-types.
ct_network


# ---- Quantitative model generation -------------------------------------------
state = c(ligands, receptors,signaling_nodes)

celltype = c(paste0("CT",1:n_cell_types))
nodetypes = c("0",celltype)
n_node_types = length(nodetypes)


prodKinetics = "MA"  # new feature. 


prodCoeff = vector("list",n_node_types)
names(prodCoeff) = nodetypes

# Empty node type: there is no source term
prodCoeff[[1]] = matrix(0,
						nrow=length(state),
						ncol=length(state),
						dimnames = list(state,state))


# setup interaction coefficients for each cell type
for (iNT in 2:n_node_types){
	iCT = iNT-1
	prodCoeff[[iNT]] = matrix(0,
							  nrow=length(state),
							  ncol=length(state),
							  dimnames = list(state,state))
	
	
	ct_network[[iCT]]$prodCoeff <- generate_prodCoeff(nrow(ct_network[[iCT]])) 
	
	# no acrtivation of TFs
	ct_network[[iCT]] <- ct_network[[iCT]] %>% mutate(prodCoeff = ifelse(to %in% all_of(TFs), prodCoeff/100,prodCoeff))
	
	for (iInter in 1:nrow(ct_network[[iCT]])){
		
		from = ct_network[[iCT]]$from[iInter]
		to = ct_network[[iCT]]$to[iInter]
		prodCoeff[[iNT]][to,from] = ct_network[[iCT]]$prodCoeff[iInter]
		
	}
	
	# add ligand binding reactions:
	# check which receptors are in the cell type and add the corresponding 
	# ligand-receptor interaction
	existing_receptors = c(ct_network[[iCT]]$from, ct_network[[iCT]]$to) %>%
		unique() %>% 
		grep("^R",x = .,value = TRUE)
	binding <- tibble(receptor = existing_receptors ) %>% 
		mutate(receptor_type = substr(receptor,2,2)) %>%
		mutate(ligand = paste0("L",receptor_type)) %>%
		mutate(prodCoeff = generate_prodCoeff(n()))
	
	for (iInter in 1:nrow(binding)){
		
		from = binding$ligand[iInter]
		to = binding$receptor[iInter]
		prodCoeff[[iNT]][to,from]  = binding$prodCoeff[iInter]
		
	}
	
}

# Degradation coefficients -----------------------
degradCoeff = vector("list",n_node_types)
names(degradCoeff) = nodetypes

for (iNT in 1:n_node_types){
	degradCoeff[[iNT]]  = runif(length(state),0.05,0.1)
	names(degradCoeff[[iNT]])  = state
	
}

# Initial state activation levels ---------------

initState = vector("list",n_node_types)
names(initState) = nodetypes
# CT0 -> no IC
initState[[1]] = rep(0,length(state))
names(initState[[1]]) <- state

for (iNT in 2:n_node_types){
	iCT = iNT -1 
	
	initState[[iNT]]  = generate_initial_condition(length(state))
	initState[[iNT]][grepl("^L",state)] = init_cond_ligands
	initState[[iNT]][grepl("^R",state)] = init_cond_receptors
	# cell states which are not expressed should have no IC
	initState[[iNT]][state %in% ct_parameters[[iCT]]$KO_signaling] = 0
	
	names(initState[[iNT]])  = state
}


# diffusion coefficients ----------------------------
diffCoeff = rep(0,length(state))
names(diffCoeff) = state

diffCoeff[names(diffCoeff) %in% ligands] = generate_diffusion_coeffs(length(ligands))

