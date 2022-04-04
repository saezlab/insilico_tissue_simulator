

library(tidyverse)
library(igraph)
library(CellNOptR)


# plot a network using cellnopt layout
show_by_cellnopt <- function(graph){
    
    sif = igraph::as_data_frame(graph,what = "edges") %>% 
        as_tibble() %>%
        select(from, sign,to)
    
    # write to temp file and read back with CellNopt
    tmp_file <- tempfile(fileext = ".sif")
    sif %>%
        write_tsv(tmp_file,col_names = FALSE)
    
    pkn <- CellNOptR::readSIF(tmp_file)
    
    CellNOptR::plotModel(pkn)
}

# convert the igraph graph object to a CellNOptR prior knowledge network object
get_pkn_from_graph <- function(graph){
    
    sif = igraph::as_data_frame(graph,what = "edges") %>% 
        as_tibble() %>%
        select(from, sign,to)
    
    # write to temp file and read back with CellNopt
    tmp_file <- tempfile(fileext = ".sif")
    sif %>%
        write_tsv(tmp_file,col_names = FALSE)
    
    pkn <- CellNOptR::readSIF(tmp_file)
}

# reverse the directions of interactions: 
# source: https://lists.nongnu.org/archive/html/igraph-help/2013-07/msg00085.html
graph_reverse <- function (graph) {
    if (!is.directed(graph))
        return(graph)
    e <- get.data.frame(graph, what="edges")
    ## swap "from" & "to"
    neworder <- 1:length(e)
    neworder[1:2] <- c(2,1)
    e <- e[neworder]
    names(e) <- names(e)[neworder]
    graph.data.frame(e, vertices = get.data.frame(graph, what="vertices"))
}



# generate a suitable graph:
# uses barabasi scale free model then does some heuristic modifications: 
# - removes the largest hub
# - reverse the edge directions (many to few nodes --> few nodes to many)
# - this would still result trees. To make the structure more realistic, we add
# random connections between nodes" "crosstalk"
graph_template <- function(N_nodes,barabasi_power = 1, barabasi_m = 1, p_neg = 0.1,cross_talk_ratio=0.5){
    # starting from scale free network
    
    g = igraph::barabasi.game(n = N_nodes,power = barabasi_power,m = barabasi_m)
    
    E(g)$sign = ifelse(runif(length(E(g)))>p_neg,1,-1)
    V(g)$name = make.names(1:length(V(g)))
    
    # we dont need to find largest component, becuase the barabasi algorithm always 
    # give connected graphs
    
    # reverse the graph edges of the graph 
    g = graph_reverse(g)
    
    # then find and remove the biggest hub that have only outgoing edges
    node_deg <- igraph::degree(graph = g,mode = "out")
    g = g - names(which.max(node_deg))
    
    # remove unconnected nodes 
    node_deg <- igraph::degree(graph = g,mode = "all")
    g = g - names(node_deg)[node_deg==0]
    
    
    
    # At this point we have trees (we detect by components). We add some random edges to complicated the 
    # structure "simulate" cross-talk effects
    
    add_random_edges <- function(g,n_crosstalks){
        node_deg_out <- igraph::degree(graph = g,mode = "out")
        node_deg_in <- igraph::degree(graph = g,mode = "in")
        
        
        # connect from node, that is not at the bottom and to node that is not at the top
        # from_node = sample(names(node_deg_out)[node_deg_out>0],size = n_crosstalks,replace = TRUE)
        from_node = sample(names(node_deg_out),size = n_crosstalks,replace = TRUE)
        to_node = sample(names(node_deg_in)[node_deg_in>0],size = n_crosstalks,replace = TRUE)
        
        # remove potential self-loops
        is_self_loop = from_node == to_node
        from_node = from_node[!is_self_loop]
        to_node = to_node[!is_self_loop]
        
        g2 = add_edges(g,c(rbind(from_node,to_node)),sign = 1)
    }
    
    
    g = add_random_edges(g,n_crosstalks = floor(length(E(g))*cross_talk_ratio) )
    
}

# generate a suitable graph:
# uses erdos renyi model then finds the largest component. 
# - this approach is outdated. The uniform edge-degree distribution resulted in 
# unrealistic scenarios
graph_template_erdos <- function(N_nodes,N_edges,p_neg = 0.1){
    # erdos renyi graph: edges are uniformly distributed
    # add sign with p_neg = .1 probability
    
    g = igraph::erdos.renyi.game(n = N_nodes,p.or.m = N_edges,type = "gnm",directed = TRUE)
    # g = igraph::barabasi.game(n = N_nodes,power = 1,m = 1)
    
    
    E(g)$sign = ifelse(runif(length(E(g)))>p_neg,1,-1)
    V(g)$name = make.names(1:length(V(g)))
    
    
    # find the largest component
    clu <- igraph::components(g, mode = c("weak"))
    largest_component <- clu$membership %>% table() %>% which.max()
    g <- igraph::induced_subgraph(g,clu$membership==largest_component)
}


# use some heuristics to select nodes for stimulated, inhibited and measured:
# stimulated: a node that has no incoming edge; or randomly selected nodes with 
#	outgoing edges and then incoming edges removed
# we remove all the nodes that are not reachable from any of the stimulated nodes
# measured nodes: randomly selected nodes that are not stimulated 
# inhibited nodes: randomly selected nodes that are not stimulated and have outgoing edge
# returns: an igraph object, with a type attribute. 
assign_node_types <- function(g,N_stimulated, N_inhibited, N_measured){
    
    V(g)$in_deg <- igraph::degree(graph = g,mode = "in")
    V(g)$out_deg <- igraph::degree(graph = g,mode = "out")
    
    
    # designing the stimulated nodes: 
    # stimulated nodes should have indegree 0 (no incoming edges)
    # we need a reasonable number of stimulated nodes : N_stimulated
    
    stimulated_nodes = V(g)$name[V(g)$in_deg == 0]
    
    if(length(stimulated_nodes) >= N_stimulated){
        # if there are too many, sample from them
        stimulated_nodes = sample(stimulated_nodes, size = N_stimulated)
        
    }else if(length(stimulated_nodes) < N_stimulated) {
        # add nodes, but be careful, we need to remove incoming edges
        # stimulated node should have downstream neighbors, so positive out_deg
        possible_nodes <- V(g)$name[V(g)$out_deg > 0]
        
        new_stimulated_nodes = sample(possible_nodes,size = N_stimulated-length(stimulated_nodes) )
        
        stimulated_nodes = c(stimulated_nodes,new_stimulated_nodes)
        
        # remove the incoming edges
        edgeSeq <- unlist(igraph::incident_edges(g, v = V(g)[V(g)$name %in% new_stimulated_nodes], mode = c("in")))
        g <- delete_edges(g,edgeSeq)
    }
    
    # make sure that edges from stimulated nodes activates their target:
    stim_out_edges <- unlist(igraph::incident_edges(g, v = V(g)[V(g)$name %in% stimulated_nodes], mode = c("out")))
    E(g)$sign[stim_out_edges] = 1
    
    
    # keep only nodes of the network that are reachable from the inputs
    
    # take the neighbourhood that can be reached from the input nodes
    stimuli_downstreams = igraph::make_ego_graph(g,order=30,nodes = stimulated_nodes, mode = "out")
    
    # stimuli_downstreams is a list, we construct the set of all nodes that are reachable
    neighbourhood_names <- lapply(stimuli_downstreams,function(x)V(x)) %>% unlist() %>% names() %>% unique() %>% sort()
    neighbourhood_ids = V(g)[V(g)$name %in% neighbourhood_names]
    
    # take the subgraph that are reachable from the inputs
    g_new <- induced_subgraph(g,neighbourhood_ids)
    
    #show_by_cellnopt(g_new)
    
    # design measured nodes:
    # select one from each set of nodes that are downstream of stimuli, which
    # makes sure that each stimuli is obeservable at least by one measurement,
    # then select the rest randomly
    
    if(N_measured<N_stimulated) stop("at least as many measured nodes are needed as the number of stimuli.")
    
    measured_nodes <- lapply(stimuli_downstreams,function(x){
        nodes <- names(V(x))
        possible_nodes <- setdiff(nodes,stimulated_nodes)
        measured <- sample(possible_nodes,1)
        }) %>% unlist()
    N_already_measured <- length(measured_nodes)
    # select the rest randomly: 
    possible_nodes <- V(g_new)$name
    possible_nodes <- possible_nodes[!possible_nodes %in% c(stimulated_nodes,measured_nodes)]
    
    if(length(possible_nodes)<N_measured) warning("number of measured nodes decreased")
    
    measured_nodes = c(measured_nodes, 
                       sample(possible_nodes,
                              size = min(length(possible_nodes),N_measured-N_already_measured)))
    
    # design inhibited nodes: 
    # nodes, that are not stimulated and has outgoing edge 
    possible_nodes <- V(g_new)$name
    possible_nodes <- possible_nodes[!possible_nodes %in% stimulated_nodes]
    possible_nodes <- possible_nodes[!possible_nodes %in% V(g_new)$name[V(g_new)$out_deg == 0]]
    inhibited_nodes = sample(possible_nodes,size = N_inhibited)
    
    
    V(g_new)$peturbation = case_when(V(g_new)$name %in% stimulated_nodes ~ "stimulated",
                                     V(g_new)$name %in% inhibited_nodes ~ "inhibited",
                                     TRUE ~ "node")
    V(g_new)$is_measured = V(g_new)$name %in% measured_nodes
    
    return(g_new)
}


# cnodata_template_from_graph
#
# takes a igraph object g with attributes name and type
# creates a CNOdata object:
# exp conditions are setup according to the stimulated, inhibited and measured nodes
# exp data is setup with data dependeing which nodes are reachable from stimulated and inhibited nodes
# returns a CellNOptR::CNOdata object

cnodata_template_from_graph <- function(g){
    
    if(is.null(V(g)$peturbation) || is.null(V(g)$name) || is.null(V(g)$is_measured)) stop("the input graph must have name, peturbation and is_measured attributes")
    
    files<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
    cnolist  = CellNOptR::CNOlist(files[[1]])
    
    node_tbl <- tibble(node = V(g)$name, perturbation=V(g)$peturbation, measured =V(g)$is_measured)
    cue_tbl <- node_tbl %>% filter(perturbation %in% c("stimulated","inhibited"))
    stim_tbl <- node_tbl %>% filter(perturbation %in% c("stimulated"))
    inhib_tbl <- node_tbl %>% filter(perturbation %in% c("inhibited"))
    measured_tbl<- node_tbl %>% filter(measured)
    
    n_exp = 1 + nrow(stim_tbl) +  nrow(stim_tbl)*nrow(inhib_tbl)
    
    cues = matrix(0, nrow=n_exp,ncol = nrow(cue_tbl))
    colnames(cues) = c(stim_tbl$node, inhib_tbl$node)
    
    # first experiment: all inputs are zero
    cues[1,] = 0
    filled_rows <- 1
    
    # independent stimulation: 
    cues[filled_rows+seq_along(stim_tbl$node),seq_along(stim_tbl$node)] = diag(1,length(stim_tbl$node))
    filled_rows <- filled_rows + length(stim_tbl$node)
    
    for(i in seq_along(inhib_tbl$node)){
        cues[filled_rows+seq_along(stim_tbl$node),seq_along(stim_tbl$node)] = diag(1,length(stim_tbl$node))
        cues[filled_rows+seq_along(stim_tbl$node),length(stim_tbl$node)+i] = 1
        filled_rows = filled_rows +length(stim_tbl$node)
        
    }
    
    
    inhibitors = matrix(nrow = nrow(cues),ncol = nrow(inhib_tbl))
    if(nrow(inhib_tbl)>0){
        inhibitors = cues[,(length(stim_tbl$node)+1):ncol(cues),drop=FALSE]    
    }
    stimuli = cues[,1:length(stim_tbl$node),drop=FALSE]
    # cues done	
    
    
    # create meaningful measurements
    measured = measured_tbl$node
    
    infer_measurements <- function(g,stimulated, inhibited, measured, time){
        
        if(time == 0 || length(stimulated)+length(inhibited)==0){
            
            measurement = matrix(0,nrow = 1,ncol = length(measured))
            colnames(measurement) = measured
        }else{
            
            pos_influenced <- igraph::ego(graph = g,order = 10,nodes = stimulated,mode = "out") %>%
                unlist() %>% names()
            
            neg_influenced <- igraph::ego(graph = g,order = 3,nodes = inhibited,mode = "out") %>%
                unlist() %>% names()
            
            measurement = matrix(0,nrow = 1,ncol = length(measured))
            colnames(measurement) = measured
            
            measurement[colnames(measurement) %in% pos_influenced & !colnames(measurement) %in% neg_influenced ] = 1
            
        }
        return(measurement)
    }	
    
    
    signals = list()
    signals$`0` = matrix(0,nrow = n_exp,ncol = length( measured_tbl$node))
    colnames(signals$`0`)  = measured_tbl$node
    signals$`10` = matrix(0,nrow = n_exp,ncol = length( measured_tbl$node))
    colnames(signals$`10`)  = measured_tbl$node
    variances = signals
    
    # perturbed state: 
    time = 1
    for(i_cue in 1:nrow(cues)){
        stimulated = colnames(stimuli)[as.logical(stimuli[i_cue,])]
        inhibited = colnames(inhibitors)[as.logical(inhibitors[i_cue,])]
        measurements <- infer_measurements(g,stimulated, inhibited, measured, time)
        signals$`10`[i_cue,] = measurements
    }
    
    
    
    timepoints = c(0,10)
    
    cnolist@cues = cues
    cnolist@inhibitors = inhibitors
    cnolist@stimuli = stimuli
    cnolist@signals = signals
    cnolist@variances = variances
    cnolist@timepoints = timepoints
    
    return(cnolist)
    
}

#' generate CNOdata object by model simulation
#'
#' simulates model and stores the data in cnodata object
#' @param cnodata CellNOptR data object
#' @param model CellNOptR model object
measurements_by_simulation <- function(cnodata, model){
    
    simList = CellNOptR::prep4sim(model)
    indexList = CellNOptR::indexFinder(cnodata, model)
    bString = rep(1,length(model$reacID))
    
    modelCut <- CellNOptR::cutModel(model, bString)
    simListCut <- CellNOptR::cutSimList(simList, bString)
    Sim0 <- CellNOptR::simulatorT0(CNOlist = cnodata, model = modelCut, 
                        simList = simListCut, indexList = indexList)
    simRes0 <- as.matrix(Sim0[, indexList$signals, drop = FALSE])
    Sim <- CellNOptR::simulatorT1(CNOlist = cnodata, model = modelCut, 
                       simList = simListCut, indexList = indexList)
    simRes <- as.matrix(Sim[, indexList$signals, drop = FALSE])
    
    simResults <- list(t0 = simRes0, t1 = simRes)
    
    colnames(simRes0) = colnames(cnodata@signals[[1]])
    colnames(simRes) = colnames(cnodata@signals[[2]])
    
    new_cnodata <- cnodata
    new_cnodata@signals = list(`0` = simRes0, `10`=simRes)
    return(new_cnodata)
}

#' insilico model generator
#' 
#' this is a wrapper around the different tasks to generate the insilico model:
#' 
#' @param N_nodes number of nodes in the initial network (resulting network will be smaller)
#' @param N_stimulated number of stimulated nodesZ
#' @param N_inhibited number of inhibited nodes
#' @param N_measured number of measured nodes
#' @param barabasi_power parameter for barabasi random network, see igraph::barabasi.game
#' @param barabasi_m parameter for barabasi random network igraph::barabasi.game
#' @param p_neg fraction of negative edges (0-1)
#' @param cross_talk_ratio by default 0.5. How many random edges we add to the original
#' network that is generated by the `barabasi.game` function 
#' 
#' @return Returns a list of elmenets: 
#' 
#' - `cnolist` CellNOptR data object, contains the simulated data
#' - `model` CellNOptR model object, contains the model that produced the data
#' - `ext_model` CellNOptR model object, contains the model with additional random edges
#' - `ext_model_BString` BitString with which ext_model fits perfecetly to the data
#' 
create_model <- function(N_nodes, N_stimulated, N_inhibited, N_measured,
                         barabasi_power = 1, barabasi_m = 1,
                         p_neg=0.05, cross_talk_ratio=0.5){
    
    # create graph
    g <- graph_template(N_nodes = N_nodes,
                        barabasi_power = barabasi_power,
                        barabasi_m = barabasi_m, 
                        p_neg = p_neg,
                        cross_talk_ratio = cross_talk_ratio)
    
    # select stimulated, inhibited and measured nodes
    g <- assign_node_types(g = g,
                           N_stimulated = N_stimulated,
                           N_inhibited = N_inhibited, 
                           N_measured = N_measured)
    
    # create experiment template
    cnodata = cnodata_template_from_graph(g)
    
    pkn = get_pkn_from_graph(g)
    model = CellNOptR::preprocessing(cnodata,pkn,expansion = FALSE)
    
    # simulate data
    cnodata_exp <- measurements_by_simulation(cnodata,model)
    
    
    # add edges to the model to challenge the optimization
    model_extended <- add_random_edges_to_pkn(model = model,
                                              cnodata = cnodata,
                                              n_random_edges = floor(length(model$reacID)*0.75),
                                              p_neg = p_neg)
    
    # the model generates perfect fir with a 1-string (all elements are 1). We 
    # determine the bitstring that gives perfect fit to the extended model. 
    
    optim_bitrstring <- as.numeric(model_extended$reacID %in% model$reacID)
    
    
    
    return(list(cnolist = cnodata_exp,
                model = model,
                ext_model = model_extended,
                ext_model_BString = optim_bitrstring))
    
}

#' Add random edges to the PKN object
#'  
#' it is difficult to add edges directly to the model object. 
#' we export to a SIF file, then import as a table where each 
#' row is an interaction. It is easier to add interactions there. 
#' Finally we export to file and load as a model object. 
add_random_edges_to_pkn <- function(model,cnodata, n_random_edges,p_neg){
    
    tmp_file <- tempfile(fileext = ".sif")
    CellNOptR::writeSIF(model, tmp_file)
    
    network_tbl <- read_tsv(tmp_file,col_names = FALSE,col_types = "cdc")
    colnames(network_tbl) <- c("from","sign","to")
    
    # add random edges:
    # we can draw from any node, but only to nodes that are not stimulated. 
    possible_ending_nodes <- model$namesSpecies[!model$namesSpecies %in% colnames(cnodata@stimuli)]
    possible_starting_nodes <- model$namesSpecies
    
    
    from_node = sample(possible_starting_nodes,size = n_random_edges, replace = TRUE)
    to_node = sample(possible_ending_nodes,size = n_random_edges,replace = TRUE)
        
    # remove potential self-loops
    is_self_loop = from_node == to_node
    from_node = from_node[!is_self_loop]
    to_node = to_node[!is_self_loop]
        
    random_edge_table <- tibble(from = from_node,sign =1, to = to_node) %>%
        dplyr::mutate(sign = ifelse(runif(n())>p_neg,1,-1))

    # bind and shuffle the rows, then write to tmp file to read it back
    bind_rows(network_tbl,random_edge_table) %>%
        sample_frac() %>%
        unique() %>%
        write_tsv(tmp_file,append = FALSE,col_names = FALSE)
    
    model <- CellNOptR::readSIF(tmp_file)
    return(model)
}

# exports the model and data and creates h5 description to the target directory
export_case_study <- function(syn_data,target_dir){
    dir.create(target_dir)
    saveRDS(syn_data,file = file.path(target_dir,"insilico_model_data.RDS"))
    export_model_to_hdf5(syn_data$cnolist,syn_data$ext_model, file.path(target_dir,"insilico_model_data.h5"))
    
}

