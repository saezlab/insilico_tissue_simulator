# coded: A.Gabor
# last modification: 24.04.2020
#
# node class: stores the cell-type specific information
# IMC class: simulates the image.
#
# last modifications: 
# - 24.04.2020: juxtacrine signaling disabled due to a bug 
node =	R6::R6Class("node",
				   public = list(
				   	type = list(),
				   	coord = c(x=0,y=0),
				   	nStates = 0,
				   	state=c(),
				   	
				   	prodCoeff = c(),
				   	prodKinetics = "MA", # alternatively, MM
				   	prodMM_Kconst = c(), # matrix like to the prodCoeff
				   	prodMM_Vmax = c(), # matrix like to the prodCoeff
				   	degradCoeff = c(),
				   	juxtaCoeff = c(),
				   	
				   	
				   	initialize=function(type="0",coord=c(x=0,y=0),state,prodCoeff,
				   						degradCoeff,juxtaCoeff,
				   						prodKinetics = "MA",prodMM_Kconst=NULL,prodMM_Vmax=NULL){
				   		
				   		self$type = type
				   		self$coord = coord
				   		self$prodCoeff = prodCoeff
				   		self$degradCoeff = degradCoeff
				   		self$juxtaCoeff = juxtaCoeff
				   		self$prodKinetics = prodKinetics
				   		if(prodKinetics == "MM"){
				   			self$prodMM_Kconst = prodMM_Kconst
				   			self$prodMM_Vmax = prodMM_Vmax
				   		}
				   		
				   		self$state = state
				   		self$nStates = length(state)
				   		
				   	},
				   	
				   	production = function(state=self$state){
				   		
				   		if( self$type == "0") {
				   			prod = rep(0,self$nStates)
				   			return(prod)
				   		}
				   		
				   		if(self$prodKinetics == "MA"){
				   			
				   			prod = self$prodCoeff %*% state	
				   			
				   			
				   		}else if(self$prodKinetics == "MM"){
				   			# browser()
				   			# prod = rep(1,self$nStates)
				   			
				   			stateMat = matrix(rep(state,length(state)),nrow = length(state),byrow = TRUE)
				   			prod = rowSums(stateMat * self$prodMM_Vmax / (self$prodMM_Kconst + stateMat + 0.0001))
				   			
				   			# for(iState in 1:self$nStates){
				   			# 	prod[iState] = sum(self$prodMM_Vmax[iState,]*state / (self$prodMM_Kconst[iState,] + state + 0.0001))
				   			# }
				   			
				   		}else(stop("invalid production kinetics"))
				   		
				   		return(prod)
				   		
				   	},
				   	
				   	degradation = function(state=self$state){
				   		deg = self$degradCoeff * state
				   		return(deg)
				   		
				   	},
				   	# This is wrong: it returns a scalar instead of a vector
				   	# Not an essentinal feature, there it is removed. 
				   	# juxtacrine = function(NN){
				   	# 	coeff=self$juxtaCoeff
				   	# 	
				   	# 	juxtacrine = 0
				   	# 	for(N in NN){
				   	# 		juxtacrine = juxtacrine + sum(coeff[,N$type] * N$state)
				   	# 	}
				   	# 	
				   	# 	return(juxtacrine)
				   	# },
				   	
				   	plotProduction = function(){
				   		#self = I_healthy$nodeTable[[2]]
				   		library(igraph)
				   		
				   		G = graph_from_adjacency_matrix(t(self$prodCoeff),
				   										mode="directed",weighted=TRUE)
				   		plot(G,main=paste0("Production in celltype: ",self$type),edge.label=E(G)$weight)
				   		
				   	}
				   )
)


IMC =
	R6::R6Class("IMC",
				public = list(
					nodeTable = matrix(),
					nStates= c(),
					states=c(),
					celltypes=c(),
					dims=c(),
					diffCoeff=c(),
					boundaryCondition = "reflective",
					simulatedStates = list(), # contains the states for each time. 
					# Each list item is a matrix: columns are state variables, each 
					# row is a gridpoint: columns first convention
					# (1,1), (2,1),...(i,1), then (1,2), (2,2),... (j,2) so on. 
					prodKinetics = "MA", # alternatve MM for Michaelis Menten 
					
					initialize=function(image,configFile,randomizeIC=FALSE){
						
						self$dims = dim(image)
						
						# load the configFile
						config = self$checkConfig(configFile)
						
						self$nStates = length(config$state)
						self$states = config$state
						self$celltypes = config$celltype
						self$diffCoeff = config$diffCoeff
						
						if(!is.null(config$prodKinetics)){
							self$prodKinetics = config$prodKinetics
						}
						
						
						# build the grid
						tmp=list()
						n = 0
						pb <- progress::progress_bar$new(total = nrow(image)*ncol(image),format = "  Initializing the model [:bar] :percent eta: :eta",)
						
						for(j in 1:nrow(image)){
							for(i in 1:ncol(image)){
								pb$tick()
								n=n+1
								type=image[i,j]
								initState = config$initState[[type]]
								
								if(randomizeIC) initState = initState + 0.3 * initState * runif(length(initState),min = -1,max = 1)
								
								tmp[[n]] = node$new(type=type,
													coord=c("x"=i,"y"=j),
													state=initState,
													prodCoeff=config$prodCoeff[[type]],
													degradCoeff=config$degradCoeff[[type]],
													juxtaCoeff=config$juxtaCoeff[[type]],
													prodKinetics=self$prodKinetics,
													prodMM_Kconst=config$prodMM_Kconst[[type]],
													prodMM_Vmax=config$prodMM_Vmax[[type]])
								
							}
						}
						dim(tmp) = dim(image)
						self$nodeTable = tmp
					},
					
					# checks and loads the configuration file
					checkConfig = function(configFile){
						conf = new.env()
						source(configFile,local = conf)
						cat("processing config file....\n")
						
						requiredConfs = c( "diffCoeff",   "degradCoeff", "celltype",   "state" ,  "initState")
						
						if(any(! requiredConfs %in% names(conf) )){
							print("missing variables in the configuration file: ")
							print (requiredConfs[which(!requiredConfs %in% names(conf))])
							stop('add missing variables to Conf file')
						}
						cat("\tall required fields found. \n")
						
						nStates = length(conf$state)
						cat("\tnumber of states: ", nStates,"\n")
						
						nCellType = length(conf$celltype)
						cat("\tnumber of cell types: ", nCellType,"\n")
						
						# checking diffusion coefficients
						if(!all( conf$state == names(conf$diffCoeff))){
							cat("!!! no diffusion coefficient found for the following state(s):\n")
							print(conf$state[which(!conf$state %in% names(conf$diffCoeff))])
							print("diffCoeff must be a named vector corresponding to each state variable, in the same order as defined in the state variable!")
							stop()
						}
						
						# checking degradation coefficients
						
						if(!all( conf$celltype %in% names(conf$degradCoeff))){
							cat("!!! no degradCoeff found for the following celltype(s):\n")
							print(conf$celltype[which(!conf$celltype %in% names(conf$degradCoeff))])
							print("degradCoeff must be a list of vectors corresponding to each celltype. Each vector is named, must follow the order of the state variable !")
							stop("see above")
						}else{
							nameCheck = unlist(lapply(conf$degradCoeff, function(m) all(conf$state ==  names(m) ) ))
							if(!all(nameCheck)){
								cat("\nameing check failed for degradCoeff[['")
								cat(names(nameCheck)[which(!nameCheck)])
								cat("']]\n")
								stop()
							}
							
						}
						
						
						# checking production coefficients
						# By default: Mass action (MA)
						# alternative: Michelis Menten (MM)
						
						if(is.null(conf$prodKinetics) | conf$prodKinetics=="MA"){
							if(is.null(conf$prodCoeff)) stop("prodCoeff is required in the config file for MA kinetics")
							if(!all( conf$celltype %in% names(conf$prodCoeff))){
								cat("!!! no production coefficient found for the following celltype(s):\n")
								print(conf$celltype[which(!conf$celltype %in% names(conf$prodCoeff))])
								print("prodCoeff must be a list of matrices corresponding to each celltype. The Matrix' dimension must be nStates x nStates !")
								stop("see above")
							}else{
								dimCheck = unlist(lapply(conf$prodCoeff, function(m) all(dim(m)==c(nStates,nStates)) ))
								if(!all(dimCheck)){
									cat("\tdimensionality check failed for prodCoeff[['")
									cat(names(dimCheck)[which(!dimCheck)])
									cat("']]\n")
									stop()
								}
								
							}
							
						}else if(conf$prodKinetics=="MM"){
							if(is.null(conf$prodMM_Kconst) | is.null(conf$prodMM_Vmax)) {
								stop("prodMM_Vmax and prodMM_Kconst are required in the config file for MM kinetics")
							}
							# Checking K constant
							if(!all( conf$celltype %in% names(conf$prodMM_Kconst))){
								cat("!!! no prodMM_Kconst coefficient found for the following celltype(s):\n")
								print(conf$celltype[which(!conf$celltype %in% names(conf$prodMM_Kconst))])
								print("prodMM_Kconst must be a list of matrices corresponding to each celltype. The Matrix' dimension must be nStates x nStates !")
								stop("see above")
							}else{
								dimCheck = unlist(lapply(conf$prodMM_Kconst, function(m) all(dim(m)==c(nStates,nStates)) ))
								if(!all(dimCheck)){
									cat("\tdimensionality check failed for prodMM_Kconst[['")
									cat(names(dimCheck)[which(!dimCheck)])
									cat("']]\n")
									stop()
								}
								
							}
							# checking the Vmax parameter
							if(!all( conf$celltype %in% names(conf$prodMM_Vmax))){
								cat("!!! no prodMM_Vmax coefficient found for the following celltype(s):\n")
								print(conf$celltype[which(!conf$celltype %in% names(conf$prodMM_Vmax))])
								print("prodMM_Vmax must be a list of matrices corresponding to each celltype. The Matrix' dimension must be nStates x nStates !")
								stop("see above")
							}else{
								dimCheck = unlist(lapply(conf$prodMM_Vmax, function(m) all(dim(m)==c(nStates,nStates)) ))
								if(!all(dimCheck)){
									cat("\tdimensionality check failed for prodMM_Vmax[['")
									cat(names(dimCheck)[which(!dimCheck)])
									cat("']]\n")
									stop()
								}
								
							}
							
						}
						
						# checking initial state
						
						if(!all( conf$celltype %in% names(conf$initState))){
							cat("!!! no initial state found for the following celltype(s):\n")
							print(conf$celltype[which(!conf$celltype %in% names(conf$initState))])
							print("initState must be a list of vectors corresponding to each celltype. Each vector is named, must follow the order of the state variable !")
							stop("see above")
						}else{
							nameCheck = unlist(lapply(conf$initState, function(m) all(conf$state ==  names(m) ) ))
							if(!all(nameCheck)){
								cat("\nameing check failed for initState[['")
								cat(names(nameCheck)[which(!nameCheck)])
								cat("']]\n")
								stop()
							}
							
						}
						
						
						
						return(conf)
					},
					
					plotCells = function(){
						# TODO bring the colors to the users configFile
						r3 <- structure(vapply(self$nodeTable, function(x)factor(as.character(x$type),levels=self$celltypes), factor(1,levels=self$celltypes)), dim=dim(self$nodeTable))
						image(1:nrow(r3),1:ncol(r3),t(r3),col=c("#F99349","#C2BEA1", "#3A5E8F","#FAFF40" ),ylim=c(nrow(r3)+0.5,0.5),xlab="Y",ylab="X")
					},
					get_Cells = function(){
						r3 <- structure(vapply(self$nodeTable, function(x)factor(as.character(x$type),levels=self$celltypes), factor(1,levels=self$celltypes)), dim=dim(self$nodeTable))
					},
					plotState = function(nodeTable = self$nodeTable, species=self$states,showText = F){
						library(RColorBrewer)
						
						redcols <- brewer.pal(9, 'Reds')
						
						species = match.arg(species)
						r3 <- structure(vapply(nodeTable, function(x)x$state[species], numeric(1)), dim=dim(self$nodeTable))
						image(1:nrow(r3),1:ncol(r3),t(r3),col=redcols,ylim=c(nrow(r3)+0.5,0.5),xlab="Y",ylab="X")
						if(showText){
							for(i in 1:nrow(r3)){
								text(i,1:ncol(r3),labels = format( t(r3[i,]),digits=1),ylim=c(nrow(r3)+0.5,0.5))
							}
						}
					},
					get_state = function(species=self$states){
						species = match.arg(species)
						r3 <- structure(vapply(self$nodeTable, function(x)x$state[species], numeric(1)), dim=dim(self$nodeTable))
					},
					get_neighbours = function(node){
						coord = node$coord
						# print("get_neighbourhood, dim of nodeTable:")
						# print(dim(self$nodeTable))
						
						if(self$boundaryCondition=="reflective"){
							
							if(coord["x"]+1 < nrow(self$nodeTable))	southNB = self$nodeTable[[coord["x"]+1,coord["y"]]]
							else southNB = self$nodeTable[[coord["x"]-1,coord["y"]]]
							
							if(coord["x"]-1 > 0)	northNB = self$nodeTable[[coord["x"]-1,coord["y"]]]
							else northNB = self$nodeTable[[coord["x"]+1,coord["y"]]]
							
							if(coord["y"]+1 < ncol(self$nodeTable))	eastNB = self$nodeTable[[coord["x"],coord["y"]+1]]
							else eastNB = self$nodeTable[[coord["x"],coord["y"]-1]]
							
							
							if(coord["y"]-1 > 0)	westNB = self$nodeTable[[coord["x"],coord["y"]-1]]
							else westNB = self$nodeTable[[coord["x"],coord["y"]+1]]
							
							
							
							
						}else stop('boundary condition unknown.')
						list("north"=northNB,"east"=eastNB,"south"=southNB,"west"=westNB)
						
						
					},
					diffusion = function(node){
						#browser()
						NN = self$get_neighbours(node)
						
						
						# coeff = c(ECM =0,
						# 		  TGFb=0.2,
						# 		  PDGF=0.2,
						# 		  WNT =0.2,
						# 		  CTGF=0.2,
						# 		  bcatenin=0,
						# 		  SMAD=0,
						# 		  TGFb_TF = 0,
						# 		  PDGF_TF = 0,
						# 		  WNT_TF = 0,
						# 		  CTGF_TF = 0)
						coeff=self$diffCoeff
						
						diffusion = coeff *(NN[["north"]]$state + NN[["east"]]$state + NN[["south"]]$state +NN[["west"]]$state - 4*node$state)
						
						
					},
					# dt: step size
					# rel_noise: parameter to control stochasticity. in each step 
					# noise with normal distribution, mean = 0 and sd = rel_noise*State is added to the states. 
					# keep it small (0.01)
					simulate = function(Nstep, dt=1, rel_noise = 0, run_parallel = FALSE ){
						# simulate the diffusion-reaction system
						# M is the node table with the cells
						# V is a state array: 20*20*7 (dimX*dimY*nStates)
						#
						nodeTable = self$nodeTable
						
						V =  t(vapply(nodeTable, function(x)x$state, rep(1,self$nStates)))
						
						# dim(V) = c(self$dims[[1]],self$dims[[2]],self$nStates)
						# dimnames(V) = list(paste0("R",1:self$dims[[1]]),paste0("C",1:self$dims[[1]]),self$states)
						self$simulatedStates[[1]] = V
						
						if(run_parallel){
							# multi thread
							library(foreach)
							
							pb <- progress::progress_bar$new(total = Nstep,format = "  simulating [:bar] :percent eta: :eta",)
							for(itime in 1:Nstep){
								pb$tick()
								
								grid_points = expand.grid(1:nrow(nodeTable),1:ncol(nodeTable))
								t0 = Sys.time()
								
								fV = foreach(i = grid_points$Var1,
											 j = grid_points$Var2) %do% {
											 	
											 	increment = dt*(
											 		nodeTable[[i,j]]$production() - 
											 			nodeTable[[i,j]]$degradation() +
											 			self$diffusion(nodeTable[[i,j]]))
											 	
											 	state_update = nodeTable[[i,j]]$state + increment + rnorm(self$nStates,mean = 0,sd = nodeTable[[i,j]]$state*rel_noise)
											 	# make sure we dont go to negative. 
											 	state_update[state_update<0] = 0 
									return(state_update)
								}
								
								# Now we want to trainsform the list of vectors to a 
								# 3D array, where the first 2 dimensions are the spatial
								# coordinates and the species are along the 3rd dimension.
								# 
								# Here is a small study, why we have to do unlist and then 
								# permuting the dimensions: 
								# 
								# grid_points = expand.grid(0:9,0:9)
								# fv = foreach(i = grid_points$Var1,
								# 			 j = grid_points$Var2) %do% {
								# 			 	
								# 			 	state_update = (1:22) +0.1*i + 0.01*j 
								# 			 	return(state_update)
								# 			 }
								# unlist(fv) %>% array(data = .,dim = c(22,10,10)) %>% aperm(.,perm = c(2,3,1)) 
								# 
								# After unlist, we will have a vector, where species follow each other,
								# of different coordinates. Forming them into an array
								# will put the species along the columns, so dim=1 is species and 
								# dim 2 and 3 are spacial coordinates. So we need a permutation that 
								# puts the 1st dim into 3, 2nd into 1 and 3rd into 2. 
								
								V <- unlist(fV) %>%
									array(data = .,dim = c(self$nStates,nrow(nodeTable),ncol(nodeTable))) %>%
									aperm(.,perm = c(2,3,1)) 
								dt <- Sys.time() - t0
								print("parallel step done")
								print(dt)
								
								
								# store the updated state:
								self$simulatedStates[[itime+1]] = V
								
								# update nodeTable's state
								for(i in 1:nrow(nodeTable)){
									for(j in 1:ncol(nodeTable)){
										nodeTable[[i,j]]$state = V[i,j,]
									}
								}
							}
							
						}else{
							# single thread
							pb <- progress::progress_bar$new(total = Nstep,format = "  simulating [:bar] :percent eta: :eta",)
							for(itime in 1:Nstep){
								pb$tick()
								
								nodeTable_dims = dim(nodeTable)
								
								dim(nodeTable) <- c(prod(dim(nodeTable)),1)
								nStates = self$nStates
								
								calc_increment <- function(N){
									#	browser()
									#if(all(N[[1]]$coord == c(7,1))) browser()
									N = N[[1]]
									increment <- dt*(
										N$production() - 
											N$degradation() +
											self$diffusion(N) +  
											rnorm(nStates,mean = 0,sd = N$state*rel_noise) )
									
									state_update = N$state + increment 
									
									# make sure we dont go to negative. 
									state_update[state_update<0] = 0 
									
									return(state_update)
									
								}
								
								state_update = t(apply(nodeTable,1,calc_increment))
								colnames(state_update) <- self$states
								
								# dim(state_update) <- c(nStates,nodeTable_dims[[1]],nodeTable_dims[[2]])
								# V <- aperm(state_update,perm = c(2,3,1)) 
								
								
								# store the updated state:
								self$simulatedStates[[itime+1]] = state_update
								
								# update nodeTable's state
								nodeTable <- lapply(1:length(nodeTable), function(i,NT){
									NT[[i]]$state = state_update[i,]
									return(NT[[i]])
									
								},nodeTable ) %>% as.matrix()
								dim(nodeTable) <- nodeTable_dims
								
								# for(i in 1:nrow(nodeTable)){
								# 	for(j in 1:ncol(nodeTable)){
								# 		nodeTable[[i,j]]$state = V[i,j,]
								# 	}
								# }
							}
						}
						
					},
					plotsimStates =  function(V = self$simulatedStates,itime=1, ispecies=1){
						# plots the state at the defined timepoint
						species = self$states[[ispecies]]
						slice <- V[[itime]][,ispecies] %>% 
							matrix(data =  .,nrow = self$dims[[1]],ncol = self$dims[[2]])
						
						maxValue = max(1,max(slice))
						
						#Vplot = apply(t(V[[itime]][,,ispecies]), 2, rev)
						
						lattice::levelplot(slice,
										   col.regions = rev(heat.colors(101)),  at=seq(0,maxValue,length.out = 100),
										   main= paste0(species,"@",itime))
						
						
					},
					
					plotTotal =  function(V = self$simulatedStates, ispecies=1){
						# plots the state at the defined timepoint
						
						
						# apply sum on each species: 
						timecourse <- lapply(V,colSums) %>%
							bind_rows() %>%
							dplyr::mutate(n_step = 1:n(),.before=1)
						
						timecourse %>%
							select(n_step,ispecies+1) %>%
							gather(variable, value,-n_step) %>%
							ggplot() + geom_line(aes(n_step,value,col=variable))
						
						
					},
					
					plotNode = function(V = self$simulatedStates,coord=c(1,1), ispecies=1){
						# plots the state at the defined timepoint
						
						species = self$states[ispecies]
						dims = self$dims
						
						totalValue = unlist(lapply(V,function(x)x[(coord[[2]]-1)*dims[[1]] + coord[[1]], ispecies]))
						
						plot(1:length(V),totalValue,main=paste0(species," @x:",coord[[1]],",y:",coord[[2]]),xlab="time",ylab=species,type = "b")
						
					},
					
					
					toTable = function(){
						V = self$simulatedStates
						
						add_meta <- function(itime,ldf){
							ldf[[itime]] %>% as_tibble() %>%
								mutate(row = (1:n()-1) %% self$dims[[1]]+1,
										  col = floor((1:n()-1)/self$dims[[1]])+1,
										  time = itime, .before=1)
						}
						
						table_out <- lapply(1:length(V),add_meta,V) %>% bind_rows()
						
						return(table_out)
					},
					# export the node information with info on the neighbours.
					exportWithNeighbours = function(measuredSpecies = self$states){
						out = list()
						
						M = self$nodeTable
						counter = 1
						for(i in 1:nrow(M)){
							for(j in 1:ncol(M)){
								# i = 2
								# j = 2
								v = data.frame(xCoord = i,yCoord = j)
								
								node = M[[i,j]]
								v$type = node$type
								
								
								v = cbind(v,as.data.frame(t(node$state[measuredSpecies])))
								
								
								
								NN = self$get_neighbours(node)
								NN_state = plyr::ldply(NN,function(x)x$state,.id = "direction")
								rownames(NN_state) = NN_state$direction
								NN_state$direction = NULL
								
								NN_state = NN_state[,colnames(NN_state) %in% measuredSpecies]
								
								NN_state$type = plyr::ldply(NN,function(x)x$type,.id = "direction")$V1
								
								m=c()
								for(irow in 1:nrow(NN_state)){
									tmp = NN_state[irow,]
									colnames(tmp) = paste(rownames(tmp),colnames(tmp),sep=".")
									m= c(m,tmp)
								}
								v = cbind(v,as.data.frame(m))
								out[[counter]] = v
								counter = counter + 1
							}
						}
						
						out = do.call("rbind",out)
						
						
						
					},
					
					
					
					recordGIF=function(frames = 1:length(self$simulatedStates),folderName){
						
						V = self$simulatedStates
						nSpecies = self$nStates
						system(paste0("mkdir ", folderName))
						pb = progress::progress_bar$new(total=nSpecies)
						cat('>>> generating pictures...\n')
						
						# calculates the maximum value for each species
						maxValue_time = do.call("rbind",lapply(M,function(Mi){
							# calculates the maximum for each species across all the nodes on a plate
							apply(Mi,MARGIN = 3,max)}
						))
						maxValue_time = lapply(V,function(x)apply(x,2,max)) %>% bind_rows()
							
											# take max over time too:
						maxValue = apply(maxValue_time,2,max)
						
						for(ispecies in 1:nSpecies){
							pb$tick()
							species = self$states[ispecies]
							png(file=paste0(folderName,"/",species,"_%02d.png"), width=600, height=600)
							for(itime in frames){
								#plot.new()
								slice <- V[[itime]][,ispecies] %>% 
									matrix(data =  .,nrow = self$dims[[1]],ncol = self$dims[[2]])
								
								print(lattice::levelplot(slice,col.regions = rev(heat.colors(101)),  at=seq(0,max(maxValue[[ispecies]],1),length.out = 100),main=paste0(species," @",itime)))
								
							}
							dev.off()
							
						}
						cat('>>> generating GIFs...\n')
						for(ispecies in 1:nSpecies){
							species = self$states[ispecies]
							system(paste0("convert -delay 10 ",folderName,"/",species,".png ",folderName,"/",species,".gif"))
							cat(folderName,"/",species,".gif is ready.\n",sep = "")
						}
					}
					
				)
	)
