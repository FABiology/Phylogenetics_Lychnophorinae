#The script below is freely available and does not offer any guarantee. The use of it is entirely the responsibility of the
#user, as well as any possible eventuality resulting from its misuse.
#This code was tested using R version 3.6.3 and castor 1.6.6, on Windows 10.

#The script was created for phylogenetic analyzes using the "castor" package in order to maintain only the essentials for
#its execution and visualization of the results. Several notes were made in order to help the user. Curly braces "{}" were
#used when necessary to speed up the process. A loop in parallel computing was applied in some cases to automate certain
#processes, however, individual execution is also available. The script was built so that the user can make as few manual
#adjustments as possible, however, feel free to improve it and share it with the community.

#Script author: Fábio Alves.

{memory.limit(size = 2*memory.limit()) #Increased limit of available memory
  Packages <- c("ape", "dplyr", "mgsub", "castor", "strap", "phyloch", "phytools", "geoscale", "doParallel",
                "ggplot2", "gganimate", "ggthemes", "hrbrthemes", "dygraphs", "tidyverse", "htmlwidgets",
                "phylobase", "phylosignal", "caper")
  lapply(Packages, library, character.only = T)
  numCores <- detectCores() #Detect the number of CPU cores on the current host. It is important for parallel computing
  registerDoParallel(numCores)} #Register the parallel backend with the foreach package

#Read the tree for analysis ("tree" object)
#Read the caracter states document ("tip.states" and "tip_states" objects. Both have different purposes along the script)
{tree = read_tree( file = "./tree.tre",
                   edge_order = "cladewise",
                   include_edge_lengths = T,
                   look_for_edge_labels = T,
                   look_for_edge_numbers = T,
                   include_node_labels = T,
                   underscores_as_blanks = F,
                   check_label_uniqueness = F,
                   interpret_quotes = F,
                   trim_white = T)
  
  tip.states = read.table("./tip_states.csv", h = T, sep = ";", dec = ",")
  tip_states = lapply(tip.states, as.numeric)
  Ntips = length(tree$tip.label) #Number of tips
  Nnodes = tree$Nnode #Number of nodes
  root_age = get_tree_span(tree)$max_distance} #Root age
#node.depth(tree) #Return the depths of nodes and tips. Count number of species of nodes and tips
#branching.times(tree) #Computes the branching times of a tree, that is the distance from each node to the tips


#Plot tree with posterior probabilities, ages 95% HPD and geologic time scale
tree.plotter <- function(basepath, xmin, xmax) {
  base_tree <- read.nexus(paste(basepath, "tree_beast.tre", sep = ""))
  base_tree$root.time <- max(nodeHeights(base_tree))
  annot_tree <- read.beast(paste(basepath, "tree_beast.tre", sep = ""))
  annot_tree$root.time <- max(annot_tree$height)
  #age_table <- read.table(paste(basepath, "age_ranges.txt", sep = ""), stringsAsFactors = F)
  params <- read.table(paste(basepath, "loganalyser_params.txt", sep = ""), header = T, stringsAsFactors = F)
  if (is.null(annot_tree$`CAheight_95%_HPD_MIN`)) {
    annot_tree$min_ages <- annot_tree$`height_95%_HPD_MIN`
    annot_tree$max_ages <- annot_tree$`height_95%_HPD_MAX`
  } else {
    annot_tree$min_ages <- annot_tree$`CAheight_95%_HPD_MIN`
    annot_tree$max_ages <- annot_tree$`CAheight_95%_HPD_MAX`
  }
  if (length(params$mean[params$statistic == "offset"] != 0)) {
    offset <- params$mean[params$statistic == "offset"]
    base_tree$root.time <- base_tree$root.time + offset
    annot_tree$min_ages <- annot_tree$min_ages + offset
    annot_tree$max_ages <- annot_tree$max_ages + offset
  } else {
    base_tree$root.time <- base_tree$root.time
    annot_tree$min_ages <- annot_tree$min_ages
    annot_tree$max_ages <- annot_tree$max_ages
  }
  #age_mat <- cbind(as.numeric(age_table[,2]), as.numeric(age_table[,3]))
  #rownames(age_mat) <- age_table[,1]
  #colnames(age_mat) <- c("FAD", "LAD")  #ages = age_mat,
  geoscalePhylo(ladderize(base_tree, right = T), show.tip.label = F, x.lim = c(xmin, xmax), y.lim = c(2, 77), #quat.rm = T,
                units = c("Period", "Epoch", "Age"), cex.tip = 0.5, cex.age = 1, cex.ts = 0.85, width = 1.5, tick.scale = 5)
  nodelabels(round(annot_tree$posterior, 2), adj = c(1.1, 1.1), frame = "n", cex = 0.7, col = "firebrick2", font = 2)
  tiplabels(sub("_", " ", base_tree$tip.label), adj = 0, frame = "n", cex = 0.7, col = "black", font = 3, offset = 0.1)
  T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  for(i in (Ntip(base_tree) + 1):(base_tree$Nnode + Ntip(base_tree))) {
    lines(x = c(T1$root.time - annot_tree$min_ages[i - Ntip(base_tree)],
                T1$root.time - annot_tree$max_ages[i - Ntip(base_tree)]),
          y = rep(T1$yy[i], 2), lwd = 6, lend = 0, col = make.transparent("blue", 0.4))}}

#tiff("TREE.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
tree.plotter("./Scripts/", -0.8, 11)
#dev.off() #Activate to export

#Phylogenetic Signal (Fritz and Purvis' D)
d.phy = phylo.d(tip.states, tree, bladder_like, names.col = sp, permut = 10000)
tiff("01.Caducous_D.tiff", units = "px", width = 1920*2, height = 1080*2, res = 300) #Activate to export
plot(d.phy, main = paste0("Pappus Caducous D Value = ",round(d.phy$DEstimate,4)))
legend(legend = c(paste0("P-bownian: ", d.phy$Pval0), paste0("P-random: ", d.phy$Pval1)),"topleft", inset = 0.05)
dev.off() #Activate to export

#Phylogenetic Autocorrelation Function Of A Numeric (continuous) Trait ####
ACF = get_trait_acf(tree, tip_states$floret, Npairs=1e7, Nbins=6)
plot(ACF$distances, ACF$autocorrelations, type="l", xlab="distance", ylab="ACF")


#Trait Depth (All) - Phylogenetic signal ####
#Parallel loop process for calculating the trait depth of each character state
trait_depth_all = foreach(i = seq_along(tip_states), .packages = "castor", .combine = cbind) %dopar% {
  trait_depth = get_trait_depth(tree,
                                tip_states[[i]],
                                min_fraction = 0.9,
                                count_singletons = T,
                                singleton_resolution= 0,
                                weighted = F,
                                Npermutations = 1000000)}

#Write .csv document with all trait depths
write.table(trait_depth_all, "./trait_depth.csv", sep = ";", dec = ".")

#Trait Depth (Individual) - Phylogenetic signal ####
#Process for calculating the trait depth of one character state at a time
trait_depth = get_trait_depth(tree,
                              tip_states$setose, #Select here the column corresponding to the desired character state
                              min_fraction = 0.9,
                              count_singletons = T,
                              singleton_resolution= 0,
                              weighted = F,
                              Npermutations = 1000000)

#Print the mean depth and standard deviation of the individual trait depth
#cat(sprintf("Mean depth = %g, std = %g\n",trait_depth$mean_depth,sqrt(trait_depth$var_depth)))

#Plot Trait Depth ####
#$positive_clades #Indices of tips and nodes (from 1 to Ntips+Nnodes) that were found to be positive in the trait
#You need to change the title and some limits according to your tree tip/node value. The rest is optional.
#The green number shows the positive species for the corresponding node. The number in red is just a tip/node ID
{plot.phylo(ladderize(tree), y.lim = 77, x.lim = 12, show.tip.label = F)
  title(main = "Setose")
  legend("left", legend = c("Mean Depth", trait_depth$mean_depth, "P Value", trait_depth$P), box.col = "white")
  nodelabels(node = trait_depth$positive_clades, adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
  tiplabels(text = tree$node.label[trait_depth$positive_clades], trait_depth$positive_clades,
            pch = 19, frame = "n", cex = 0.7, col = "black")
  tiplabels(text = trait_depth$positives_per_clade[trait_depth$positive_clades], trait_depth$positive_clades, 
            adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "darkgreen", font = 2)
  b = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(trait_depth$positives_per_clade)))) %>% 
                             filter(as.data.frame(trait_depth$positives_per_clade) == 1 %in% 
                                      row.names(as.data.frame(trait_depth$positives_per_clade))) %>% filter(. < 80)))
  tiplabels(text = sub("_", " ", tree$tip.label[b]), b, 
            adj = 0, frame = "n", cex = 0.7, col = "royalblue1", font = 3, offset = 0.2)
  c = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(trait_depth$positives_per_clade)))) %>% 
                             filter(as.data.frame(trait_depth$positives_per_clade) == 0 %in% 
                                      row.names(as.data.frame(trait_depth$positives_per_clade))) %>% filter(. < 80)))
  tiplabels(text = sub("_", " ", tree$tip.label[c]), c, 
            adj = 0, frame = "n", cex = 0.7, col = "black", font = 3, offset = 0.2)
  axisPhylo(1)}

#$positives_per_clade #Number of descending tips per clade (tip or node) that were positive in the trait
#Use "!" before the filter, inverts the selection
{plot.phylo(ladderize(tree), tip.color = "royalblue1", cex = 0.7, label.offset = 0.2, y.lim = 77, x.lim = 12)
  t = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(trait_depth$positives_per_clade)))) %>% 
                             filter(as.data.frame(trait_depth$positives_per_clade) > 0 %in% 
                                      row.names(as.data.frame(trait_depth$positives_per_clade))) %>% filter(. > 79)))
  
  nodelabels(node = t, adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
  tiplabels(text = tree$node.label[t], t, pch = 1, frame = "n", cex = 0.7, col = "black")
  tiplabels(text = trait_depth$positives_per_clade[t], t, 
            adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "darkgreen", font = 2)
  axisPhylo(1)}


#Ancestral State Reconstruction (ASR) step ####
#Map states of a discrete trait to integers (All)
#Parallel loop process for map all states
state_map_all = foreach(i = seq_along(tip.states), .packages = "castor", .combine = rbind) %dopar% {
  state_map = map_to_state_space(tip.states[[i]])}

#ASR with Markov (Mk) models (All)
#Parallel loop process for calculating the ASR of each character with all known tip states
amk_all = foreach(i = c(17, 22, 26, 29), .packages = "castor", .combine = rbind) %dopar% {
  amk = asr_mk_model(tree, state_map_all[,3][[i]], state_map_all[,1][[i]], rate_model = "ARD", store_exponentials = T, 
                     include_ancestral_likelihoods = T, Ntrials = numCores, Nthreads = numCores)}

#ASR with hidden state prediction via Mk model max-likelihood (All)
#Parallel loop process for calculating the ASR of each character with tip states absence, i.e. hidden state
hmk_all = foreach(i = c(2, 5, 8, 12), .packages = "castor", .combine = rbind) %dopar% {
  hmk = hsp_mk_model(tree, state_map_all[,3][[i]] -1, state_map_all[,1][[i]] -1, rate_model = "ARD", store_exponentials = T, 
                     include_likelihoods = T, Ntrials = numCores, Nthreads = numCores)}


#Map states of a discrete trait to integers (Individual)
#Process for map one state at a time. You need to specify the character column from tip.states object changing the number
state_map = map_to_state_space(tip.states[,2])

#ASR with Mk models and rerooting (Individual)
#Process for calculating the ASR of one character at a time with all known tip states
{amk = asr_mk_model(tree, state_map$mapped_states, state_map$Nstates, rate_model = "ARD",include_ancestral_likelihoods = T,
                    store_exponentials = T, Ntrials = numCores, Nthreads = numCores)$ancestral_likelihoods
  amk_estimated = max.col(amk[1:(Nnodes),])
  amk_estimated} #Print estimated node states

#ASR with hidden state prediction via Mk model max-likelihood (Individual)
#Process for calculating the ASR of one character at a time with tip states absence, i.e. hidden state
{hmk = hsp_mk_model(tree, state_map$mapped_states -1, state_map$Nstates -1, rate_model = "ARD", store_exponentials = T,
                    include_likelihoods = T, Ntrials = numCores, Nthreads = numCores)$likelihoods
  hmk_estimated = max.col(hmk[1:(Ntips + Nnodes),])
  hmk_estimated} #Print estimated tip and node states


#Plot ASR for AMK approach with pie chart likelihoods and geologic time scale
amk.plotter <- function(basepath, xmin, xmax) {
  tree$root.time <- root_age
  geoscalePhylo(tree, show.tip.label = F, x.lim = c(xmin, xmax), y.lim = c(2, 77),
                units = c("Period", "Epoch", "Age"), cex.tip = 0.5, cex.age = 1, cex.ts = 0.85, width = 1.5, tick.scale = 5)
  colors = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c", "#ffffbf")
  nodelabels(round(amk, 2), adj = c(0.5, 0.5), frame = "n", cex = 0.25, col = "#00000000", font = 2, piecol = colors, pie = amk)
  tiplabels(sub("_", " ", tree$tip.label), adj = 0, frame = "n", cex = 0.7, col = "black", font = 3, offset = 0.1)
  rr = data.frame(V1 = map_to_state_space(tip.states[,18])$mapped_states -1,
                  V2 = map_to_state_space(tip.states[,19])$mapped_states -1,
                  V3 = map_to_state_space(tip.states[,20])$mapped_states -1,
                  V4 = map_to_state_space(tip.states[,21])$mapped_states -1)
  tiplabels(rr, adj = c(0.5, 0.5), frame = "n", cex = 0.12, col = "#00000000", font = 2, piecol = colors, pie = rr)
  legend(x = "bottomleft", title = "Habit", legend = state_map$state_names, pch = 21, pt.cex = 2.5, pt.bg = colors,
         cex = 1.5, bty = "n")}

#tiff("AMK.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
amk.plotter("./Scripts/", -0.8, 11)
#dev.off() #Activate to export

#Plot ASR for HMK approach with pie chart likelihoods and geologic time scale
hmk.plotter <- function(basepath, xmin, xmax) {
  tree$root.time <- root_age
  geoscalePhylo(ladderize(tree, right = T), show.tip.label = F, x.lim = c(xmin, xmax), y.lim = c(2, 77),
                units = c("Period", "Epoch", "Age"), cex.tip = 0.5, cex.age = 1, cex.ts = 0.85, width = 1.5, tick.scale = 5)
  colors = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c", "#ffffbf")
  nodelabels(round(hmk[((max(tree$Nnode)+2):max(tree$edge)),1], 2), adj = c(0.5, 0.5), frame = "n", cex = 0.25,
             col = "#00000000", font = 2, piecol = colors, pie = hmk[((max(tree$Nnode)+2):max(tree$edge)),])
  tiplabels(sub("_", " ", tree$tip.label), adj = 0, frame = "n", cex = 0.7, col = "black", font = 3, offset = 0.1)
  tiplabels(round(hmk[(1:(max(tree$Nnode)+1)),1], 2), adj = c(0.5, 0.5), frame = "n", col = "#00000000", font = 3,
            cex = 0.12, piecol = colors, pie = hmk[(min(tree$edge):(max(tree$Nnode)+1)),])
  legend(x = "bottomleft", title = "Pappus type", legend = state_map$state_names[-1], pch = 21, pt.cex = 2.5, pt.bg = colors,
         cex = 1.5, bty = "n")}

#tiff("HMK.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
hmk.plotter("./Scripts/", -0.8, 11)
#dev.off() #Activate to export

#Plot ASR for AMK approach
#You need to change the title, legend and some limits according to your tree tip/node value. The rest is optional.
{plot.phylo(ladderize(tree), y.lim = 77, x.lim = 12, show.tip.label = F)
  title(main = "A.S.R.: Habit")
  legend("left", legend = c("C = Caulirosula", "S = Shrub", "T = Trees"),
         text.col = c("firebrick2", "royalblue1", "darkgreen"), box.col = "white")
  #nodelabels(((Ntips + 1):(Ntips + Nnodes)), adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
  tiplabels(tip.states[,12], adj = 0, frame = "n", cex = 0.6, col = "darkgreen", font = 2)
  nodelabels(mgsub(amk_estimated[1:Nnodes], state_map$name2index, state_map$state_names),
             adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "purple3", font = 2)
  b = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(state_map$mapped_states)))) %>% 
                             filter(as.data.frame(state_map$mapped_states) == 1 %in% 
                                      row.names(as.data.frame(state_map$mapped_states)))))
  tiplabels(text = sub("_", " ", tree$tip.label[b]), b, 
            adj = 0, frame = "n", cex = 0.7, col = "firebrick2", font = 3, offset = 0.2)
  c = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(state_map$mapped_states)))) %>% 
                             filter(as.data.frame(state_map$mapped_states) != 1 %in% 
                                      row.names(as.data.frame(state_map$mapped_states)))))
  tiplabels(text = sub("_", " ", tree$tip.label[c]), c, 
            adj = 0, frame = "n", cex = 0.7, col = "royalblue1", font = 3, offset = 0.2)
  axisPhylo(1)}


#Plot ASR for HMK approach
#You need to change the title, legend and some limits according to your tree tip/node value. The rest is optional.
{plot.phylo(ladderize(tree), y.lim = 77, x.lim = 12, show.tip.label = F)
  title(main = "A.S.R.: Pappus type")
  legend("left", legend = c("S = Setose", "P = Paleaceous"), text.col = c("firebrick2", "royalblue1"), box.col = "white")
  #nodelabels(((Ntips + 1):(Ntips + Nnodes)), adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
  tiplabels(mgsub(hmk_estimated[1:Ntips], state_map$name2index -1, state_map$state_names),
            adj = 0, frame = "n", cex = 0.6, col = "darkgreen", font = 2)
  nodelabels(mgsub(hmk_estimated[80:(Ntips + Nnodes)], state_map$name2index -1, state_map$state_names),
             adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "purple3", font = 2)
  b = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(hmk_estimated)))) %>% 
                             filter(as.data.frame(hmk_estimated) == 1 %in% 
                                      row.names(as.data.frame(hmk_estimated))) %>% filter(. < 80)))
  tiplabels(text = sub("_", " ", tree$tip.label[b]), b, 
            adj = 0, frame = "n", cex = 0.7, col = "royalblue1", font = 3, offset = 0.2)
  c = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(hmk_estimated)))) %>% 
                             filter(as.data.frame(hmk_estimated) != 1 %in% 
                                      row.names(as.data.frame(hmk_estimated))) %>% filter(. < 80)))
  tiplabels(text = sub("_", " ", tree$tip.label[c]), c, 
            adj = 0, frame = "n", cex = 0.7, col = "firebrick2", font = 3, offset = 0.2)
  axisPhylo(1)}

#Mirror Cladogram
{syn_txt = read.csv2("./syn.txt")
  syn = as.matrix(syn_txt)
  syn = as.numeric(syn)}
{flor_txt = read.csv2("./flor.txt")
  flor = as.matrix(flor_txt)
  flor = as.numeric(flor)}

mirror.plotter <- function(basepath) {
  layout(matrix(1:3,1,3),widths=c(0.45,0.10,0.45))
  colors = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c", "#ffffbf")
  plotBranchbyTrait(ladderize(tree), syn, "nodes", palette=colorRampPalette(colors[1:2]), show.tip.label = F,
                    legend = F, prompt=F, cex=1)
  rr = data.frame(V1 = map_to_state_space(tip.states[,27])$mapped_states -1,
                  V2 = map_to_state_space(tip.states[,28])$mapped_states -1)
  tiplabels(rr, adj = c(0.58, 0.5), frame = "n", cex = 0.35, col = "#00000000", font = 2, piecol = colors, pie = rr)
  nodelabels(round(amk, 2), adj = c(0.5, 0.5), frame = "n", cex = 0.45, col = "#00000000", font = 2, piecol = colors, pie = amk)
  state_map_syn = map_to_state_space(tip.states[,26])
  legend(x = "bottomleft", title = "Syncephalum", legend = state_map_syn$state_names, pch = 21, pt.cex = 3.7, pt.bg = colors,
         cex = 2.2, bty = "n", inset = 0.05)
  plot.new()
  plot.window(xlim=c(-0.1,0.1),ylim=c(1, length(tree$tip.label)))
  par(cex=1)
  text(rep(0,length(tree$tip.label)), 1:length(tree$tip.label), sub("_", " ", tree$tip.label), cex=0.8, font = 3)
  plotBranchbyTrait(ladderize(tree), flor, "nodes", palette=colorRampPalette(colors[1:4]), show.tip.label = F,
                    legend = F, prompt=F, cex=1, direction="leftwards")
  tiplabels(round(hmk[(1:(max(tree$Nnode)+1)),1], 2), adj = c(0.42, 0.5), frame = "n", col = "#00000000", font = 3,
            cex = 0.35, piecol = colors, pie = hmk[(min(tree$edge):(max(tree$Nnode)+1)),])
  nodelabels(round(hmk[((max(tree$Nnode)+2):max(tree$edge)),1], 2), adj = c(0.5, 0.5), frame = "n", cex = 0.45,
             col = "#00000000", font = 2, piecol = colors, pie = hmk[((max(tree$Nnode)+2):max(tree$edge)),])
  state_map_flor = map_to_state_space(tip.states[,12])
  legend(x = "bottomright", title = "Florets", legend = state_map_flor$state_names[-1], pch = 21, pt.cex = 2.5, pt.bg = colors,
         cex = 1.5, bty = "n", inset = 0.05)}
#tiff("Mirror.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
mirror.plotter("./Desktop/")
tiff("Mirror.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
mirror.plotter("./Desktop/")
dev.off() #Activate to export

#Lineages Through Time (LTT) step ####
#Calculate the number of lineages represented in the tree at various time points
{ltt = count_lineages_through_time(tree, Ntimes = 1000, include_slopes = T)
ltt_values = read.table("./ltt_values.csv", h = T, sep = ";", dec = ".")}

#LTT Plot
#tiff("LTT.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
{plot(ltt$times, ltt$lineages, type = "l", xlab = "Time", ylab = "Species")
  title("Lineages Through Time (LTT)")}
#dev.off() #Activate to export

#LTT Static
#tiff("LTT_static.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
ggplot(aes(-(times[nrow(ltt_values):0]), lineages), data = ltt_values) +
  theme_ipsum(axis_title_size = 12, grid_col = "gray30") +
  geom_line(size = 1) +
  geom_area(fill = "blue", alpha = 0.25) +
  labs(x = "Time (Ma)", y = "Species", title = "Lineages Through Time (LTT)")
#dev.off() #Activate to export

#LTT Animate
ltt_animate = ltt_values %>%
  ggplot(aes(-(times[nrow(ltt_values):0]), lineages)) +
  geom_line() +
  geom_point() +
  labs(x = "Time (Ma)", y = "Species", title = "Lineages Through Time (LTT)") +
  theme_ipsum() +
  transition_reveal(times)

animate(ltt_animate, end_pause = 15, width = 960, height = 540)

anim_save("./LTT_animate.gif") #Export animated LTT

#LTT Dygraph
ltt_dygraph = dygraph(ltt_values, main = "Lineages Through Time (LTT)", xlab = "Time (Ma)", ylab = "Species") %>%
  dyOptions(labelsUTC = T, fillGraph = T, fillAlpha = 1, drawGrid = F, colors = "royalblue") %>%
  dyRangeSelector() %>%
  dyCrosshair(direction = "vertical") %>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = F)  %>%
  dyRoller(rollPeriod = 1)

saveWidget(ltt_dygraph, file = "./LTT_dygraph.html") #Export interactive LTT

#Guesstimate start lambdas & mus based on the LTT
#Lambda can be inferred from the present-day slope of the LTT, correcting for the sampling fraction
{guess = list()
  guess$birth_rates = tail(ltt$relative_slopes,1)
  guess$death_rates = guess$birth_rates} #mu is typically close to lambda

#Specify upper bounds for the model parameters
#This helps constrain the fitting within reasonable parameter space, preventing excessive computation time
upper = list(transition_matrix=1e4/root_age, birth_rates=1e4/root_age, death_rates=1e4/root_age)

#Limit the time (seconds) spent on each likelihood evaluation
#This further helps avoid extreme parameter regions where likelihood calculations become very slow
max_model_runtime = max(1,Ntips/1e3)

#Fit model (BiSSE/MuSSE/HiSSE) ####
system.time({fit = fit_musse(tree,
  #Nstates = state_map$Nstates,
  Nstates = (state_map$Nstates-1)*(state_map$Nstates-1), #Activate for HiSSE
  NPstates = state_map$Nstates-1, #Activate for HiSSE
  proxy_map = c(as.integer(state_map$name2index[-1]-1), as.integer(state_map$name2index[-1]-1)), #Activate for HiSSE
  #state_names	= state_map$state_names,
  state_names = 1:((state_map$Nstates-1)*(state_map$Nstates-1)), #Activate for HiSSE
  tip_pstates = state_map$mapped_states-1, #For HiSSE use -1
  oldest_age = root_age,
  root_prior = "likelihoods",
  root_conditioning	= "madfitz",
  transition_rate_model = "ARD",
  birth_rate_model = "ARD",
  death_rate_model = "ARD",
  Ntrials = numCores,
  Nthreads = numCores,
  max_start_attempts = 10,
  optim_algorithm = "subplex",
  optim_max_iterations = 10000,
  #optim_max_evaluations = 10000,
  #include_ancestral_likelihoods = T,
  Nbootstraps = numCores,
  Ntrials_per_bootstrap = numCores,
  D_temporal_resolution = root_age,
  #first_guess = guess, #Deactivate for HiSSE
  upper = upper,
  max_model_runtime = max_model_runtime,
  verbose = T,
  verbose_prefix = "")
if(!fit$success){
  cat(sprintf("ERROR: Fit failed, with the following error: %s\n",fit$error))
}else{
  #Fitting succeeded, so print results
  cat(sprintf("Fit finished successfully\nLog-likelihood: %.10g\nFitted parameters:\n",fit$loglikelihood))
  print(fit$parameters)
  #Compare fitted birth rates to true values
  errors = (fit$parameters$birth_rates - guess$birth_rates)
  relative_errors = errors/guess$birth_rates
  cat(sprintf("Fit relative birth-rate errors:\n"))
  print(relative_errors)
  #Print 95%-confidence interval
  cat(sprintf("CI95 for lambda: %g-%g",
              fit$CI95lower$birth_rates,
              fit$CI95upper$birth_rates))
  cat(sprintf("CI95 for mu: %g-%g",
              fit$CI95lower$death_rates,
              fit$CI95upper$death_rates))}})


#Fit pulled speciation rates (PSR) of birth-death models on a time grid ####
{Ngrid = root_age
age_grid = seq(from=0,to=root_age,length.out=Ngrid)
psr = fit_hbd_psr_on_grid(tree,
                          oldest_age = root_age,
                          age_grid = age_grid,
                          min_PSR = 0,
                          max_PSR = +100,
                          condition = "auto",
                          Ntrials = numCores,
                          Nthreads = numCores,
                          max_model_runtime = 1)
if(!psr$success){
  cat(sprintf("ERROR: Fitting failed: %s\n",psr$error))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",psr$loglikelihood))
  #Plot fitted PSR
  plot( x = psr$age_grid,
        y = psr$fitted_PSR,
        main = 'Fitted PSR',
        xlab = 'Age',
        ylab = 'PSR',
        type = 'l',
        xlim = c(root_age,0))
  #Plot deterministic LTT of fitted model
  plot( x = psr$age_grid,
        y = psr$fitted_LTT,
        main = 'Fitted dLTT',
        xlab = 'Age',
        ylab = 'Lineages',
        type = 'l',
        #log = 'y',
        xlim = c(root_age,0))
  #Get fitted PSR as a function of age
  PSR_fun = approxfun(x=psr$age_grid, y=psr$fitted_PSR)}}


#Fit pulled diversification rates (PDR) of birth-death models on a time grid ####
{Ngrid = root_age
age_grid = seq(from=0,to=root_age,length.out=Ngrid)
pdr = fit_hbd_pdr_on_grid(tree,
                          oldest_age = root_age,
                          age_grid = age_grid,
                          min_PDR = -100,
                          max_PDR = +100,
                          guess_rholambda0 = 0.50,
                          condition = "auto",
                          Ntrials = numCores,
                          Nthreads = numCores,
                          max_model_runtime = 1)
if(!pdr$success){
  cat(sprintf("ERROR: Fitting failed: %s\n",pdr$error))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",pdr$loglikelihood))
  #Plot fitted PDR
  plot( x = pdr$age_grid,
        y = pdr$fitted_PDR,
        main = 'Fitted PDR',
        xlab = 'Age',
        ylab = 'PDR',
        type = 'l',
        xlim = c(root_age,0))
  #Get fitted PDR as a function of age
  PDR_fun = approxfun(x=pdr$age_grid, y=pdr$fitted_PDR)}}










#BiSSE model ####
bisse = fit_musse(tree,
                  Nstates = 2,
                  tip_pstates = tip_states,
                  transition_rate_model = "ARD",
                  birth_rate_model = "ARD",
                  death_rate_model = "ARD",
                  Ntrials = 12,
                  Nthreads = 12,
                  max_start_attempts = 10,
                  optim_algorithm = "subplex",
                  #optim_max_iterations = 10000,
                  #optim_max_evaluations = 10000,
                  include_ancestral_likelihoods = T,
                  Nbootstraps = 10,
                  Ntrials_per_bootstrap = 10,
                  D_temporal_resolution = 100,
                  verbose = T)
if(!bisse$success){
  cat(sprintf("ERROR: Fitting failed"))
}else{
  # compare fitted birth rates to true values
  errors = (bisse$parameters$birth_rates - parameters$birth_rates)
  relative_errors = errors/parameters$birth_rates
  cat(sprintf("BiSSE relative birth-rate errors:\n"))
  print(relative_errors)}


#HiSSE model ####
hisse = fit_musse(tree,
                  Nstates = 4,
                  NPstates = 2,
                  proxy_map = c(1,2,1,2),
                  tip_pstates = state_map$mapped_states-1,
                  transition_rate_model = "ARD",
                  birth_rate_model = "ARD",
                  death_rate_model = "ARD",
                  Ntrials = 12,
                  Nthreads = 12,
                  max_start_attempts = 10,
                  optim_algorithm = "subplex",
                  #optim_max_iterations = 10000,
                  #optim_max_evaluations = 10000,
                  include_ancestral_likelihoods = T,
                  Nbootstraps = 10,
                  Ntrials_per_bootstrap = 10,
                  D_temporal_resolution = 100,
                  verbose = T)
if(!hisse$success){
  cat(sprintf("ERROR: Fitting failed"))
}else{
  # compare fitted birth rates to true values
  errors = (hisse$parameters$birth_rates - parameters$birth_rates)
  relative_errors = errors/parameters$birth_rates
  cat(sprintf("HiSSE relative birth-rate errors:\n"))
  print(relative_errors)
  # print 95%-confidence interval for first birth rate
  cat(sprintf("CI95 for lambda1: %g-%g",
              hisse$CI95lower$birth_rates[1],
              hisse$CI95upper$birth_rates[1]))}
