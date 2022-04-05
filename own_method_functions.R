#Only works if objects have same number of reduced dims
centroid_match <- function(embeddings, clustering){
  match_df <- data.frame("clus1.1", "clus2.1")
  centroids <- matrix(ncol=length(embeddings[[1]][1,]), nrow=1)
  mtx_row_id <- c()
  #calculate centroid embedding
  for (i in 1:length(embeddings)){
    current_embed <- embeddings[[i]]
    current_cluster <- clustering[[i]]
    
    for (j in unique(current_cluster[,1])){
      cell_list <- current_cluster[,2][which(current_cluster[,1] == j)]
      centroids <- rbind(centroids,colMeans( current_embed[ which( row.names(current_embed) %in% cell_list ), ] ))
      mtx_row_id <- append(mtx_row_id, paste0(j,"_",i))
      
    }
    
  }
  centroids <- centroids[-1,]
  row.names(centroids) <- mtx_row_id
  
  #create distance matrix
  dist_mtx <-as.matrix( dist(centroids, diag = T, upper = T) )
  
  
  #remove distances between centroids on the same dataset
  dist_mtx <- dist_mtx[-which( row.names(dist_mtx) %in% paste0( unique(clustering[[1]][,1]), "_1" ) ), ]
  
  dist_mtx <- dist_mtx[ , -which( colnames(dist_mtx) %in% paste0( unique(clustering[[2]][,1]), "_2" ) )]
  
  return(dist_mtx)
  
  #minimum combination of row values
}

#Assigns cells to a node based on euclidean distance for some metric
#First column is names of individuals
#Second column is value (pseudotime)
cell_capture <- function(cell_embed, node_embed){
  output_list <- list()
  
  node_mtx <- t( replicate(length(cell_embed[,1]), node_embed[,1]) )
  cell_mtx <- replicate( length(node_embed[,1]) , cell_embed[,1])
  
  dist_mtx <- ( node_mtx - cell_mtx )^2
  
  colnames(dist_mtx) <- row.names(node_embed)
  row.names(dist_mtx) <- row.names(cell_embed)
  
  row_min <- rowMins(dist_mtx)
  
  closest_node <- as.matrix(apply( dist_mtx, 1, which.min))
  
  for (i in 1:length( colnames(dist_mtx) )){
    output_list[[colnames(dist_mtx)[i]]] <- row.names(closest_node)[which(closest_node == i)]
    
  }
  
  return(output_list)
  
}

monocle_pseudobulk <- function(mon_obj){
  output <- list()
  cell_embed <- t(mon_obj@reducedDimS)
  node_embed <- t(mon_obj@reducedDimK)
  exp_mtx <- mon_obj@assayData[["exprs"]]
  
  capture_output <- cell_capture(cell_embed, node_embed)
  
  #return(capture_output)
  
  pseudobulk_output<- matrix(ncol =1, nrow = length(exp_mtx[,1]))
  for (i in 1:length(capture_output)){
    filtered_expr_mtx <- exp_mtx[, which( colnames(exp_mtx) %in% capture_output[[i]] ) ]
    
    if (length(capture_output[[i]]) == 1){
      pseudobulk_exp <- filtered_expr_mtx
    }
    
    else{
      pseudobulk_exp <- rowSums(filtered_expr_mtx)
    }
    
    
    pseudobulk_output <- cbind(pseudobulk_output,as.matrix(pseudobulk_exp))
  }
  
  pseudobulk_output <- pseudobulk_output[,-1]
  colnames(pseudobulk_output) <- names(capture_output)
  
  output[[ paste0("pseudobulk") ]] <- pseudobulk_output
  output[[ paste0("capture") ]] <- capture_output
  
  return(output)
  
}

slingshot_pseudo_bulk <- function(exp_mtx, sce_obj, slingshot_obj,node_info){
  output <- list()
  
  for (lineage in 1:length(slingshot_obj@curves)){
    node_pseudotime <- node_info[[ paste0("lineage_", lineage, "_pseudotime") ]]
    cell_pseudotime <- data.frame(sce_obj@colData@listData[[paste0("slingPseudotime_", lineage)]],sce_obj@colData@rownames)
    cell_pseudotime <- cell_pseudotime[ order(cell_pseudotime[,1]), ]
    
    #Removes all items from variable if there is only 1 lineage
    if(length(slingshot_obj@curves) > 1){
      cell_pseudotime <- cell_pseudotime[ -which( is.na(cell_pseudotime) == T ), ]
    }

    cell_pseudotime <- data.frame(cell_pseudotime[,1], row.names = cell_pseudotime[,2])
    node_pseudotime <- data.frame(node_pseudotime[,1], row.names = row.names(node_pseudotime))
    capture_output <- cell_capture(cell_pseudotime, node_pseudotime)
    
    #if(lineage == 2){return(capture_output)}
    
    pseudobulk_output<- matrix(ncol =1, nrow = length(exp_mtx[,1]))
    
    for (i in 1:length(capture_output)){
      filtered_expr_mtx <- as.matrix(exp_mtx[, which( colnames(exp_mtx) %in% capture_output[[i]] ) ])
      
      pseudobulk_exp <- rowSums(filtered_expr_mtx)
      
      pseudobulk_output <- cbind(pseudobulk_output,as.matrix(pseudobulk_exp))
    }
    
    pseudobulk_output <- pseudobulk_output[,-1]
    colnames(pseudobulk_output) <- names(capture_output)
    
    output[[ paste0("lineage_", lineage, "_pseudobulk") ]] <- pseudobulk_output
    output[[ paste0("lineage_", lineage, "_capture") ]] <- capture_output
  }
  
  return(output)
}

dis_mtx_calculator <- function(exp_mtx1, exp_mtx2, dis_method){
  if (dis_method != "euclidean"){
    dis_matrix <- cor(cbind(exp_mtx1, exp_mtx2), method=dis_method)
    dis_matrix <- 1- dis_matrix
    
  }
  
  if (dis_method == "euclidean"){
    dis_matrix <- as.matrix(dist(rbind(t(exp_mtx1), t(exp_mtx2)), diag = T, upper = T)) 
    
  }
  
  #process the matrix to remove unnecessary rows and columns
  dis_matrix <- dis_matrix[ -which( row.names(dis_matrix) %in% colnames(exp_mtx1) ) , ]
  dis_matrix <- dis_matrix[ , -which( colnames(dis_matrix) %in% colnames(exp_mtx2) ) ]
  dis_matrix <- as.matrix(dis_matrix)

  print( pheatmap::pheatmap(as.matrix(t(dis_matrix)), cluster_rows = F, cluster_cols = F, scale="none") )
  
  
  return(dis_matrix)
}

#Estimate gene expression values for the nodes weighted by their euclidean distance in pseudotime to the cells
cell_align_node_GEV <- function(exp_mtx, node_pseudo, cell_pseudo, window_size){
  
  output <- list()
  
  for ( lineage in 2:length(cell_pseudo[1,]) ){
    #node_exp_mtx <- matrix(nrow = length(exp_mtx[,1]), ncol = length(node_pseudo[,1]))
    node_exp_mtx <- matrix(nrow = 50, ncol = length(node_pseudo[,1]))
    
    cell_exp_mtx <- exp_mtx[, which( is.na(cell_pseudo[,lineage]) != T )]
    cell_pseudo_cut <- cell_pseudo[ which( is.na(cell_pseudo[,lineage]) != T ), ]
    
    
    weight_vector <- c()
    weight_mtx <- matrix(nrow = length(cell_pseudo_cut[,1]), ncol = length(node_pseudo[,1]))
    
    for (cell in 1:length(cell_pseudo_cut[,1])){
      cell_weight_vector <- c()
      
      for (node in 1:length(node_pseudo[,1])){
        weight <- exp( -1 * ( (cell_pseudo_cut[cell,lineage] - node_pseudo[node,2])^2 / window_size^2 ) )
        cell_weight_vector <- append(cell_weight_vector, weight)
        
      }
      weight_mtx[cell,] <- cell_weight_vector
      
    }
    
    row.names(weight_mtx) <- cell_pseudo_cut[,1]
    colnames(weight_mtx) <- node_pseudo[,1]
    
    
    #mtx multiplication way
    node_exp_mtx <- cell_exp_mtx %*% weight_mtx 
    
    normalizer <- 1/(colSums(weight_mtx))
    
    for (node in 1:length(node_exp_mtx[1,])){
      node_exp_mtx[,node] <- node_exp_mtx[,node] * normalizer[node]
    }
    
    output[[paste("lineage_", lineage, "_exp_mtx")]] <- node_exp_mtx
    
  }

  return(output)
  
}

#If win size is too low, it creates 0's (not really 0's but very very small numbers) which fucks up everything
#Pseudotimes are 
cell_align_node_branch <- function(tree, end_nodes, root_node, exp_mtx, node_embed, cell_embed, node_pseudotime, cell_pseudotime, win_size){
  capture_output <- cell_capture(cell_embed, node_embed)
  
  #return(capture_output)
  
  node_exp_mtx <- matrix(nrow = length(exp_mtx[,1]), ncol = 1)
  row.names(node_exp_mtx) <- row.names(exp_mtx)
  
  one_hot_mtx <- matrix(0, nrow =length(node_embed[,1]), ncol = length(end_nodes))
  row.names(one_hot_mtx) <- row.names(node_embed)
  
  weight_mtx <- matrix(nrow = length(cell_pseudotime[,1]), ncol = length(node_pseudotime[,1]))
  
  #Calculate weight matrix
  for (cell in 1:length(cell_pseudotime[,1])){
    cell_weight_vector <- c()
    
    for (node in 1:length(node_pseudotime[,1])){
      weight <- exp( -1 * ( (cell_pseudotime[cell,2] - node_pseudotime[node,2])^2 / win_size^2 ) )
      
      if(weight == 0){ 
        print(paste("Node: ", node_pseudotime[node,1], " and cell: ", cell_pseudotime[node,1]))
        print(paste("Node: ", node_pseudotime[node,2], " and cell: ", cell_pseudotime[node,2]))
        return(c(cell_pseudotime[cell,2], node_pseudotime,2))
        }
      #print(paste( "Cell pseudo: ", cell_pseudotime[cell,2] ))
      #print(paste( "Node pseudo: ", node_pseudotime[node,2] ))
      #print("")
      
      cell_weight_vector <- append(cell_weight_vector, weight)
      
    }
    weight_mtx[cell,] <- cell_weight_vector
    
  }
  
  row.names(weight_mtx) <- cell_pseudotime$node
  colnames(weight_mtx) <- node_pseudotime$node
  
  #Find what leaf to root paths each node is on
  for (end in 1:length(end_nodes)){
    
    selection <- c(end_nodes[end])
    current_node <- end_nodes[end]
    
    while (current_node != root_node){
      current_node <- tree[which( tree[,1] == current_node), 2]
      selection <- append(selection, current_node)
      
    }
    
    one_hot_mtx[ which(row.names(one_hot_mtx) %in% selection), end] <- one_hot_mtx[ which(row.names(one_hot_mtx) %in% selection), end] + 1
    
  }
  
  for (node in 1:length(one_hot_mtx[,1])){
    one_hot_mtx_cp <- one_hot_mtx
    weight_mtx_cp <- weight_mtx
    exp_mtx_cp <- exp_mtx
    
    r_t_l_index <- which(one_hot_mtx[node,] == 1)
    
    one_hot_mtx_cp <- as.matrix(one_hot_mtx_cp[, r_t_l_index])
    one_hot_mtx_cp <- as.matrix(one_hot_mtx_cp[ which(rowSums(one_hot_mtx_cp) != 0), ])
    
    cell_list_index <- which(row.names(one_hot_mtx_cp) %in% names(capture_output))
    
    cell_list <- unname(unlist(capture_output[cell_list_index]))
    
    weight_mtx_cp <- weight_mtx_cp[ which(row.names(weight_mtx_cp) %in% cell_list) , which( colnames(weight_mtx_cp) == row.names(one_hot_mtx)[node ]) ]
    
    exp_mtx_cp <- exp_mtx_cp[ , which( colnames(exp_mtx_cp) %in% cell_list ) ]
    
    node_exp_vector <- exp_mtx_cp %*% weight_mtx_cp 
    normalizer <- 1/(sum(weight_mtx_cp))
    
    normalized_node_exp_vector <- node_exp_vector * normalizer
    
    #print(row.names(one_hot_mtx)[node])
    #print(sum(normalized_node_exp_vector))
    
    
    node_exp_mtx <- cbind(node_exp_mtx, normalized_node_exp_vector)
    #print(sum(node_exp_mtx))
    #print("")
  }
  
  node_exp_mtx <- node_exp_mtx[,-1]
  colnames(node_exp_mtx) <- row.names(one_hot_mtx)
  return(node_exp_mtx)
}

#pseudobulk or wilcoxon plus permutation or t-test
slingshot_sliding_window_de <- function(expr_mtx_1, expr_mtx_2, pseudo_1, pseudo_2, pseudo_window, method, logfc_threshold, pval_threshold){
  output <- list()
  
  #get the windows of pseudotime we need for comparing the cells 
  window_intervals <- c(0)
  while( max(window_intervals) < min( c( max(pseudo_1[,2]), max(pseudo_2[,2]) ) ) ){
    window_intervals <- append(window_intervals, ( max(window_intervals) + pseudo_window) )
  }
  
  all_grouping <- list()
  #group cells based on their pseudotime window collation 
  for (i in 1:2){
    if(i == 1){pseudotime <- pseudo_1}
    
    else{pseudotime <- pseudo_2}
    
    current_grouping <- list()
    for (window in 2:length(window_intervals)){
      
      group <- pseudotime[ which(pseudotime[,2] > window_intervals[(window-1)] & pseudotime[,2] < window_intervals[window]) ,1]
      #ensures we don't try to make a comparison with only a few cells 
      if (length(group) < 2){
        current_grouping[[ paste0(window_intervals[(window-1)] , " to ", window_intervals[window]) ]] <- c()
      }
      
      else{
        current_grouping[[ paste0(window_intervals[(window-1)] , " to ", window_intervals[window]) ]] <- group
        
      }
      
    }
    
    
    all_grouping[[paste0("condition_", i)]] <- current_grouping
    
  }
  
  
  
  if (method == "t-test"){ t_test_de(expr_mtx_1, expr_mtx_2, all_grouping, logfc_threshold, pval_threshold) }
  
  else { permutation_de(expr_mtx_1, expr_mtx_2, all_grouping) }
}

#Only includes aligned cells in analysis 
slingshot_align_only_window_de <- function(expr_mtx_1, expr_mtx_2, pseudo_1, pseudo_2, pseudo_window ,method, logfc_threshold, pval_threshold, all.genes, logfc.metric){
  minimum_limit <- max( c( min(pseudo_1[,1]), min(pseudo_2[,1]) ) )
  
  maximum_limit <- min( c( max(pseudo_1[,1]), max(pseudo_2[,1]) ) )
  
  print(minimum_limit)
  print(maximum_limit)
  
  cell_list_1 <- row.names(pseudo_1)[ which(pseudo_1[,1] < maximum_limit & pseudo_1[,1] > minimum_limit)]
  cell_list_2 <- row.names(pseudo_2)[ which(pseudo_2[,1] < maximum_limit & pseudo_2[,1] > minimum_limit)]
  
  pseudo_1 <- pseudo_1[ which(pseudo_1[,1] < maximum_limit & pseudo_1[,1] > minimum_limit) ,]
  pseudo_2 <- pseudo_2[ which(pseudo_2[,1] < maximum_limit & pseudo_2[,1] > minimum_limit) ,]
  
  
  expr_mtx_1 <- expr_mtx_1[ , which(colnames(expr_mtx_1) %in% cell_list_1) ]
  expr_mtx_2 <- expr_mtx_2[ , which(colnames(expr_mtx_2) %in% cell_list_2) ]
  
  expr_mtx_1_sum <- rowSums(expr_mtx_1)
  expr_mtx_2_sum <- rowSums(expr_mtx_2)
  
  expr_mtx_1 <- expr_mtx_1[which(expr_mtx_1_sum != 0),]
  expr_mtx_2 <- expr_mtx_2[which(expr_mtx_2_sum != 0),]
  
  gene_list <- unique(c(row.names(expr_mtx_1), row.names(expr_mtx_2)))
  gene_list <- gene_list[which(gene_list %in% row.names(expr_mtx_1))]
  gene_list <- gene_list[which(gene_list %in% row.names(expr_mtx_2))]
  
  expr_mtx_1 <- expr_mtx_1[which(row.names(expr_mtx_1) %in% gene_list),]
  expr_mtx_2 <- expr_mtx_2[which(row.names(expr_mtx_2) %in% gene_list),]
  
  output <- data.frame(gene_name = "name", p_val = 0, adj_p_val = 0 ,logfc = 0)
  
  
  #get the windows of pseudotime we need for comparing the cells 
  window_intervals <- c(minimum_limit, (minimum_limit + pseudo_window))
  window_intervals <- data.frame(minimum_limit, (minimum_limit+pseudo_window))
  print(window_intervals[1,])
  print(maximum_limit)
  print(pseudo_window)
  
  #while( max(window_intervals) < maximum_limit ){
  #  window_intervals <- append(window_intervals, ( ( max(window_intervals) - (pseudo_window/2) ) + pseudo_window ) )
  #}
  
  while( max(window_intervals) < maximum_limit){
    interval <- c( (max(window_intervals) - (pseudo_window/2)) , ( (max(window_intervals) - (pseudo_window/2)) + pseudo_window) )
    window_intervals <- rbind(window_intervals, interval)
    
  }
  
  print(window_intervals)
  
  all_grouping <- list()
  
  #group cells based on their pseudotime window collation 
  for (i in 1:2){
    if(i == 1){pseudotime <- pseudo_1}
    
    else{pseudotime <- pseudo_2}
    
    current_grouping <- list()
    for (window in 1:length(window_intervals[,1])){
      
      group <- row.names(pseudotime)[ which(pseudotime[,1] > window_intervals[window,1] & pseudotime[,1] < window_intervals[window,2])]
      
      #ensures we don't try to make a comparison with only a few cells 
      if (length(group) < 2){
        current_grouping[[ paste0(window_intervals[window,1] , " to ", window_intervals[window,2]) ]] <- c()
      }
      
      else{
        current_grouping[[ paste0(window_intervals[window,1] , " to ", window_intervals[window,2]) ]] <- group
        
      }
      
    }
    
    all_grouping[[paste0("condition_", i)]] <- current_grouping
    
  }
  
  #some cells can be present in a pseudotime interval for one trajectory but not the other
  #We remove intervals where no comparison is possible here
  
  cond1_intervals <- names(all_grouping[[1]])
  cond2_intervals <- names(all_grouping[[2]])
  
  common_intervals <- intersect(cond1_intervals, cond2_intervals)
  
  all_grouping[[1]] <- all_grouping[[1]][which(names(all_grouping[[1]]) %in% common_intervals)]
  all_grouping[[2]] <- all_grouping[[2]][which(names(all_grouping[[2]]) %in% common_intervals)]
  
  #return(all_grouping)
  
  if (method == "t-test"){ output <- window_t_test_de(expr_mtx_1, expr_mtx_2, all_grouping, logfc_threshold, pval_threshold, all.genes, logfc.metric) }
  
  else { output <- permutation_de(expr_mtx_1, expr_mtx_2, all_grouping) }
  
  #Bonferoni for the windows and genes
  
  # for (i in 1:length(output)){
  #   output[[i]][,"adj_pval"] <- output[[i]][,"adj_pval"] * length(output)
  # 
  #   output[[i]] <- output[[i]][output[[i]][,"adj_pval"] < 0.05,]
  # 
  # }
  
  return(output)
}

window_t_test_de <- function(expr_mtx_1, expr_mtx_2, cell_groupings, logfc_threshold, pval_threshold, all.genes, logfc.metric){
  output <- list()
  for (i in 1:length(cell_groupings[[1]])){
    print(i)
    #print(cell_groupings[[1]][i])
    #print(cell_groupings[[2]][i])
    expr_mtx_1_window <- as.matrix(expr_mtx_1[ , which(colnames(expr_mtx_1) %in% cell_groupings[[1]][[i]]) ])
    expr_mtx_2_window <- as.matrix(expr_mtx_2[ , which(colnames(expr_mtx_2) %in% cell_groupings[[2]][[i]]) ])
    
    gene_list <- intersect(row.names(expr_mtx_1_window), row.names(expr_mtx_2_window))
    
    expr_mtx_1_window <- expr_mtx_1_window[ which( row.names(expr_mtx_1_window) %in% gene_list) ,]
   
    expr_mtx_2_window <- expr_mtx_2_window[ which( row.names(expr_mtx_2_window) %in% gene_list) ,]

    expr_mtx_1_window <- expr_mtx_1_window[match(row.names(expr_mtx_1_window), gene_list),]
    expr_mtx_2_window <- expr_mtx_2_window[match(row.names(expr_mtx_2_window), gene_list),]

    #remove genes whose mean is zero in both conditions
    mean_mtx1 <- rowMeans(expr_mtx_1_window)
    names(mean_mtx1) <- row.names(expr_mtx_1_window)
    
    mean_mtx2 <- rowMeans(expr_mtx_2_window)
    names(mean_mtx2) <- row.names(expr_mtx_2_window)
    
    gene_df <- data.frame("gene", 0, 0)
    colnames(gene_df) <- c("gene_name", "p_val", "logfc.change")
    
    if( logfc.metric == "Seurat" ){
      condition_1 <- mean.fxn(expr_mtx_1_window)
      condition_2 <- mean.fxn(expr_mtx_2_window)
      
      logfc_change <- condition_1 - condition_2
      print("Seurat method")
    }
    
    else{
      condition_1 <- rowMeans(expr_mtx_1_window)
      condition_2 <- rowMeans(expr_mtx_2_window)
      
      logfc_change <- log2(condition_1 / condition_2)
      print("Own method")
    }
    
   
    test_result_vector <- c()
    
    for (gene in gene_list){
      test_result <- t.test(expr_mtx_1_window[gene,], expr_mtx_2_window[gene,])
      
      test_result_vector <- append(test_result_vector, test_result$p.value)
    }
    
    gene_df <- data.frame(gene_name = gene_list, logfc = logfc_change, adj_pval = test_result_vector)
    
    library(stats)
    
    gene_df[,3] <- p.adjust(gene_df[,3], method = "bonferroni", length(gene_list))
    if (all.genes == F){
      #Remove genes with non-significant p-values and logfc not in the threshold
      gene_df <- gene_df[which(gene_df[,3] < pval_threshold),]
      
      gene_df <- gene_df[which(gene_df[,2] > logfc_threshold | gene_df[,2] < (-1*logfc_threshold) ) , ]
    }
    
    output[[ names(cell_groupings[[1]])[i] ]] <- gene_df
  }
  
  return(output)
  
}

align_de <- function(expr_mtx_1, expr_mtx_2, pseudo_1, pseudo_2 ,method, logfc_threshold, pval_threshold, all.genes, logfc.metric){
  
  
  cells_1 <- pseudo_1$ID[which(pseudo_1$status == "align")]
  cells_2 <- pseudo_2$ID[which(pseudo_2$status == "align")]
  
  expr_1_subset <- expr_mtx_1[, cells_1]
  expr_2_subset <- expr_mtx_2[, cells_2]
  
  output <- t_test_de_no_grouping(expr_1_subset, expr_2_subset, logfc_threshold, pval_threshold, all.genes,logfc.metric)
  
  return(output)
}

#Did have fold change based on sum of the genes expression, but realised those could be thrown off if one of the conditions had loads more cells than the other
t_test_de <- function(expr_mtx_1, expr_mtx_2, cell_groupings, logfc_threshold, pval_threshold, all.genes){
  
  output <- list()
  
  for (i in 1:length(cell_groupings[[1]])){
    #print(cell_groupings[[1]][i])
    #print(cell_groupings[[2]][i])
    expr_mtx_1_window <- as.matrix(expr_mtx_1[ , which(colnames(expr_mtx_1) %in% cell_groupings[[1]][[i]]) ])
    expr_mtx_2_window <- as.matrix(expr_mtx_2[ , which(colnames(expr_mtx_2) %in% cell_groupings[[2]][[i]]) ])
    
    #remove genes whose sum is zero
    sum_mtx1 <- rowSums(expr_mtx_1_window)
    names(sum_mtx1) <- row.names(expr_mtx_1_window)
    
    sum_mtx2 <- rowSums(expr_mtx_2_window)
    names(sum_mtx2) <- row.names(expr_mtx_2_window)
    
    mean_mtx1 <- rowMeans(expr_mtx_1_window)
    names(sum_mtx1) <- row.names(expr_mtx_1_window)
    
    mean_mtx2 <- rowMeans(expr_mtx_2_window)
    names(sum_mtx2) <- row.names(expr_mtx_2_window)
    
    mean_mtx1 <- mean_mtx1[which(mean_mtx1 >0)]
    mean_mtx2 <- mean_mtx2[which(mean_mtx2 >0)]
    
    sum_mtx1 <- sum_mtx1[which(sum_mtx1 >0)]
    sum_mtx2 <- sum_mtx2[which(sum_mtx2 >0)]
    
    sum_mtx1 <- mean_mtx1
    sum_mtx2 <- mean_mtx2
    
    gene_list <- intersect(names(sum_mtx1), names(sum_mtx2))
    
    gene_df <- data.frame("gene", 0, 0)
    colnames(gene_df) <- c("gene_name", "p_val", "logfc.change")
    
    mean.fxn <- function(x) {
      return(log(x = mean(x = expm1(x)) + pseudocount.use, base = base))
    }
    
    for (gene in gene_list){
      test1 <- mean.fxn(expr_mtx_1[gene,])
      test2 <- mean.fxn(expr_mtx_2[gene,])
      
      logfc_change <- test1 - test2
      
      #old way
      #change <- sum_mtx1[gene] / sum_mtx2[gene]
      #logfc_change <- log2(change)
      
      if (logfc_change > logfc_threshold | logfc_change < (-1*logfc_threshold)){
        test_result <- t.test(expr_mtx_1_window[gene,], expr_mtx_2_window[gene,])
        
        row_output <- data.frame(gene, test_result$p.value, logfc_change)
        colnames(row_output) <- colnames(gene_df)
        
        gene_df <- rbind(gene_df, row_output)
      }
    }
    gene_df <- gene_df[-1,]
    
    library(stats)
    
    gene_df[,2] <- p.adjust(gene_df[,2], method = "bonferroni", length(row.names(expr_mtx_1)))
    
    if (all.genes == F){
      #Remove genes with non-significant p-values
      gene_df <- gene_df[which(gene_df[,2] < pval_threshold),]
    }
    
    output[[ names(cell_groupings[[1]])[i] ]] <- gene_df
  }
  
 
  
  return(output)
}

#Seurat way of getting LogFC 
mean.fxn <- function(x) {
  return(log(x = rowMeans(x = expm1(x)) + 1, base = 2))
  
}

t_test_de_no_grouping <- function(expr_mtx_1, expr_mtx_2, logfc_threshold, pval_threshold, all.genes,logfc.metric){
  
  bon_correction <- max( c( length(row.names(expr_mtx_1)), length(row.names(expr_mtx_2)) ) )
  
  gene_list <- intersect(row.names(expr_mtx_1), row.names(expr_mtx_2))
  

  expr_mtx_1 <- expr_mtx_1[ which( row.names(expr_mtx_1) %in% gene_list) ,]
    
  expr_mtx_2 <- expr_mtx_2[ which( row.names(expr_mtx_2) %in% gene_list) ,]
  
  expr_mtx_1 <- expr_mtx_1[match(row.names(expr_mtx_1), gene_list),]
  expr_mtx_2 <- expr_mtx_2[match(row.names(expr_mtx_2), gene_list),]
  
  #remove genes whose mean is zero in both conditions
  mean_mtx1 <- rowMeans(expr_mtx_1)
  names(mean_mtx1) <- row.names(expr_mtx_1)
    
  mean_mtx2 <- rowMeans(expr_mtx_2)
  names(mean_mtx2) <- row.names(expr_mtx_2)
    
  if( logfc.metric == "Seurat" ){
    print("in")
    condition_1 <- mean.fxn(expr_mtx_1)
    condition_2 <- mean.fxn(expr_mtx_2)
    
    logfc_change <- condition_1 - condition_2
  }
    
  else{
    print("in_")
    condition_1 <- rowMeans(expr_mtx_1)
    condition_2 <- rowMeans(expr_mtx_2)
      
    logfc_change <- log2(condition_1 / condition_2)
  }
    
    
  test_result_vector <- c()
    
  for (gene in gene_list){
    test_result <- t.test(expr_mtx_1[gene,], expr_mtx_2[gene,])
      
    test_result_vector <- append(test_result_vector, test_result$p.value)
  }
    
  gene_df <- data.frame(gene_name = gene_list, logfc = logfc_change, adj_pval = test_result_vector)

  library(stats)
    
  gene_df[,3] <- p.adjust(gene_df[,3], method = "bonferroni", bon_correction)
  
  if (all.genes == F){
    #Remove genes with non-significant p-values and logfc not in the threshold
    gene_df <- gene_df[which(gene_df[,3] < pval_threshold),]
      
    gene_df <- gene_df[which(gene_df[,2] > logfc_threshold | gene_df[,2] < (-1*logfc_threshold) ) , ]
  }
    
  output <- gene_df

  return(output)
}
  
gene_over_pseudo_display <- function(gene_list, exp_mtx1, exp_mtx2, pseudo1, pseudo2, interval){
  #Order cells by their pseudotime so we can plot easily
  
  pseudo1 <- pseudo1[pseudo1$status == "align",]
  pseudo2 <- pseudo2[pseudo2$status == "align",]
  
  pseudo1 <- pseudo1[ order(pseudo1[,1]), ]
  pseudo2 <- pseudo2[ order(pseudo2[,1]), ]
  
  for (gene in gene_list){
    
    pseudo_gene1 <- row.names(pseudo1)
    pseudo_gene2 <- row.names(pseudo2)
    
    print(length(pseudo_gene1))
    
    gene_exp1 <- exp_mtx1[gene,pseudo_gene1]
    gene_exp2 <- exp_mtx2[gene,pseudo_gene2]
    
    print(length(gene_exp1))
    
    plot1 <- cbind(pseudo1[,1] ,gene_exp1)
    plot2 <- cbind(pseudo2[,1], gene_exp2)
    
    max_pseudo <- max(c(pseudo1[,1], pseudo2[,1]))
    print(max_pseudo)
    
    print(ggplot() + geom_point(aes(plot1[,1], plot1[,2], col="Condition 1")) + geom_point(aes(plot2[,1], plot2[,2], col = "Condition 2")) + scale_x_continuous(breaks = seq(0, max_pseudo , by = signif(interval, 2))))
    #print(ggplot() + geom_point(aes(plot1[,1], plot1[,2], col="Condition 1")) + scale_x_continuous(breaks = seq(0, max_pseudo , by = interval)))
    #print(ggplot() + geom_point(aes(plot2[,1], plot2[,2], col = "Condition 2")) + scale_x_continuous(breaks = seq(0, max_pseudo , by = interval)))
    
  }
  
  
}

scale_nodes <- function(exp_mtx_1, exp_mtx_2){
  
  max_val <- max( c( max( exp_mtx_1), max(exp_mtx_2) ) )
  min_val <- min( c( min( exp_mtx_1), min(exp_mtx_2) ) )
  
  exp_mtx_1 <- ( exp_mtx_1 - min_val ) / (max_val - min_val)
  exp_mtx_2 <- ( exp_mtx_2 - min_val ) / (max_val - min_val)
  
  #exp_mtx_1 <- ( exp_mtx_1 - min(exp_mtx_1) ) / (max(exp_mtx_1) - min(exp_mtx_1))
  #exp_mtx_2 <- ( exp_mtx_2 - min(exp_mtx_2) ) / (max(exp_mtx_2) - min(exp_mtx_2))
  
  return(list(exp_mtx_1, exp_mtx_2))
  
}

#Aligning pseudotimes of matched trajectories

#Align pseudotimes
#First need to project the initial points backwards then after that squish them between the aligned nodes
#Version of this is present in the Traverse_monocle_trajectory, but works with output from slingshot_lineage_to_lineage rather than slingshot_node_maker
pseudo_align_own <- function(condition_1, condition_2,  alignment, dims){
  condition_1_order <- row.names(condition_1)
  condition_2_order <- row.names(condition_2)
  
  condition_1_pseudotime <- condition_1[,1]
  condition_2_pseudotime <- condition_2[,1]
  
  
  cond_1_counter <- 1
  cond_2_counter <- 1
  for (i in 1:length(alignment[,1])){
    x <- F
    cond_1_interval <- 1
    cond_2_interval <- 1
    
    #Initial alignment
    if(i == 1){
      while (x != T){
        #print(paste0("1: ",cond_1_counter))
        #print(paste0("2: ",cond_2_counter))
        if (condition_1_order[cond_1_counter] == alignment[i,1]){
          #print("hmm1")
          cond_1_interval <- 0
        }
        
        if (condition_2_order[cond_2_counter] == alignment[i,2]){
          #print("hmm2")
          cond_2_interval <- 0
        }
        
        cond_1_counter <- cond_1_counter + cond_1_interval
        cond_2_counter <- cond_2_counter + cond_2_interval
        #print(paste0("1 interval: ", cond_1_interval))
        #print(paste0("2 interval: ", cond_2_interval))
        if(cond_1_interval == 0 & cond_2_interval == 0){
          x <- T
        }
      }
      #Align node pseudotimes
      #Wondering what happens if condition 2 has more nodes between but a smaller pseudotime than the aligned node on condition 1
      #Assume that both trajectories span a similar pseudotime and node axis (and share same feature space)
      
      if(cond_1_counter > cond_2_counter){
        #if (condition_2[cond_2_counter, (dims+1)] < condition_1[cond_1_counter, (dims+1)]) 
        #Condense pseudotime of cond_1 down
        diff <- condition_1_pseudotime[cond_1_counter] - condition_2_pseudotime[cond_2_counter]
        #condition_1[seq(1,(cond_1_counter-1), by = 1), (dims+1)] <- condition_1[seq(1,(cond_1_counter-1), by = 1), (dims+1)] - diff
        condition_1_pseudotime <- condition_1_pseudotime- diff
        
        
        condition_1_pseudotime[cond_1_counter] <- condition_2_pseudotime[cond_2_counter]
      }
      
      
      else{
        #condense pseudotime of cond_2 down
        diff <- condition_2_pseudotime[cond_2_counter] - condition_1_pseudotime[cond_1_counter]
        #condition_2[seq(1,(cond_2_counter-1), by = 1), (dims+1)] <- condition_2[seq(1,(cond_2_counter-1), by = 1), (dims+1)] - diff
        
        condition_2_pseudotime <- condition_2_pseudotime- diff
        
        condition_2_pseudotime[cond_2_counter] <- condition_1_pseudotime[cond_1_counter]
      }
    }
    #Problem will arise if the trajectory with lots of nodes spans a larger pseudotime than the other trajectory with fewer nodes between the aligned ones
    
    #Count how many nodes it takes to reach an aligned node in each of the conditions
    #The condition with the highest number of nodes between aligned nodes needs to be condensed in terms of pseudotime
    
    #else start
    else{
      previous_cond1 <- cond_1_counter
      previous_cond2 <- cond_2_counter
      x <- F
      cond_1_interval <- 1
      cond_2_interval <- 1
      while (x != T){
        if (condition_1_order[cond_1_counter] == alignment[i,1]){
          cond_1_interval <- 0
        }
        
        if (condition_2_order[cond_2_counter] == alignment[i,2]){
          cond_2_interval <- 0
        }
        
        
        cond_1_counter <- cond_1_counter + cond_1_interval
        cond_2_counter <- cond_2_counter + cond_2_interval
        
        if(cond_1_interval == 0 & cond_2_interval == 0){
          x <- T
        }
      }
      print("problem")
      if(cond_1_counter > cond_2_counter){
        #Condense pseudotime of cond_1 down
        #Make pseudotimes the same for the aligned nodes 
        condition_1_pseudotime[cond_1_counter] <- condition_2_pseudotime[cond_2_counter]
        inbetween <- cond_1_counter - previous_cond1
        
        #If there is intervening nodes between aligned nodes
        if(inbetween != 1){
          diff <- condition_1_pseudotime[cond_1_counter] - condition_1_pseudotime[which(alignment[i-1,1] == row.names(condition_1))]
          
          interval <- diff / (inbetween)
          pseudo <- condition_1_pseudotime[which(alignment[i-1,1] == row.names(condition_1))]
          node_pseudo <- c()
          
          
          #For every intervening node that is aligned, calculate it's new pseudotime
          for (node in 1:(inbetween-1)){
            pseudo <- pseudo + interval
            node_pseudo <- append(node_pseudo, pseudo)
          }
          print("problemstart")
          print(node_pseudo)
          print(previous_cond1)
          print(cond_1_counter)
          condition_1_pseudotime[(previous_cond1+1):(cond_1_counter-1)] <- node_pseudo
          print("problemend")
        }
        
      }
      
      else{
        #condense pseudotime of cond_2 down
        condition_2_pseudotime[cond_2_counter] <- condition_1_pseudotime[cond_1_counter]
        inbetween <- cond_2_counter - previous_cond2
        
        if (inbetween != 1){
          diff <- condition_2_pseudotime[cond_2_counter] - condition_2_pseudotime[which(alignment[i-1,2] == row.names(condition_2))]
          
          interval <- diff / (inbetween)
          pseudo <- condition_2_pseudotime[which(alignment[i-1,2] == row.names(condition_2))]
          node_pseudo <- c()
          
          #For every intervening node that is aligned, calculate it's new pseudotime
          for (node in 1:(inbetween-1)){
            print(pseudo)
            pseudo <- pseudo + interval
            node_pseudo <- append(node_pseudo, pseudo)
          }
          
          
          condition_2_pseudotime[(previous_cond2+1):(cond_2_counter-1)] <- node_pseudo
          
        }
        
      }
    }
    
  } 
  
  #As we might have negative pseudotime we want to scale/normalise the pseudotime by making the lowest pseduotime value of the conditions equal to 0
  #The other pseudotimes are then boosted up by the absolute value of the lowest pseudotime
  if (min(condition_1_pseudotime) < 0 | min(condition_2_pseudotime) < 0){
    if (min(condition_1_pseudotime) < min(condition_2_pseudotime)){
      condition_1_pseudotime <- condition_1_pseudotime + abs(min(condition_1_pseudotime))
      condition_2_pseudotime <- condition_2_pseudotime + abs(min(condition_1_pseudotime))
    }
    
    else{
      condition_1_pseudotime <- condition_1_pseudotime + abs(min(condition_2_pseudotime))
      condition_2_pseudotime <- condition_2_pseudotime + abs(min(condition_2_pseudotime))
    }
  }
  
  condition_1_pseudotime <- as.matrix(condition_1_pseudotime)
  condition_2_pseudotime <- as.matrix(condition_2_pseudotime)
  
  row.names(condition_1_pseudotime) <- row.names(condition_1)
  row.names(condition_2_pseudotime) <- row.names(condition_2)
  
  x <- list(condition_1_pseudotime, condition_2_pseudotime)
  names(x) <- c("condition_1", "condition_2")
  return(x)
}

#Only works with linear trajectories at the moment 
#Version of this is present in the Traverse_monocle_trajectory, but works with output from slingshot_lineage_to_lineage rather than slingshot_node_maker
#This version doesn't use euclidean distance for the kmeans part 
#Depracated
pseudo_change_own <- function(output, cond1_sce,cond1_sling,cond2_sce, cond2_sling, dims, lin1, lin2, cond1_old_pseudo, cond2_old_pseudo){
  kmean_output <- list()
  
  for (condition in 1:2){
    if(condition == 1){
      #cond_cells <- cond1_sce@int_colData@listData[["reducedDims"]]@listData[["PHATE"]]
      cond_cells <- cond1_sce@colData@rownames
      cond_orig_pseudo <- cond1_sce@colData@listData[[paste0("slingPseudotime_", lin1)]]
      
      #Makes sure we only deal with cells which are on the lineage
      erase <- which(is.na(cond_orig_pseudo)==T)
      keep <- which(is.na(cond_orig_pseudo)==F)
      #print(cond_orig_pseudo)
      
      #print(cond_orig_pseudo)
      if(length(erase) != 0){
        cond_cells <- cond_cells[-erase]
        cond_orig_pseudo <- cond_orig_pseudo[-erase]
        
      }
      
      #print(cond_cells)
      cond_node_old <- cond1_old_pseudo
    }
    
    else{
      #cond_cells <- cond2_sce@int_colData@listData[["reducedDims"]]@listData[["PHATE"]] 
      cond_cells <- cond2_sce@colData@rownames
      cond_orig_pseudo <- cond2_sce@colData@listData[[paste0("slingPseudotime_", lin2)]]
      erase <- which(is.na(cond_orig_pseudo)==T)
      keep <- which(is.na(cond_orig_pseudo)==F)
      
      
      if(length(erase) != 0){
        cond_cells <- cond_cells[-erase]
        cond_orig_pseudo <- cond_orig_pseudo[-erase]
      }
      cond_node_old <- cond2_old_pseudo
    }
    
    cond_node_new <- output[[paste0("condition_", condition)]]
    
    
    
    kmean_cell_input <- as.matrix(as.matrix(cond_orig_pseudo))
    row.names(kmean_cell_input) <- cond_cells
    
    #return(list(cond_orig_pseudo, cond_node_old))
    
    kmean_node <- cell_capture(kmean_cell_input, cond_node_old)
    #for (cell in 1:length(cond_cells[,1])){
      #print(cell)
    #  previous_dist <- 1000
    #  for (node in 1:length(cond_node_new[, (dims+2)])){
    #    current_dist <- dist(rbind(cond_cells[cell,], cond_node_new[node,1:dims]), "euclidean")
        
    #    if (current_dist < previous_dist){
    #      previous_dist <- current_dist
    #      closest_node <- cond_node_new[node, (dims+2)]
    #    }
    #  }
    #  kmean_node[[closest_node]] <- append(kmean_node[[closest_node]], row.names(cond_cells)[cell])
      
    #}
    
    #future error might occur if one of the nodes doesn't have any cells assigned to it. Current dataset doesn't have this scenario
    #for (node in 1:length(cond_node_new[,(dims+2)])){
    #diff <- cond_node_new[node, (dims+1)] - cond_node_old[node, (dims+1)]
    
    #if(condition ==2 & node %in% c(1,2,3)){
    #print(paste0("node: ", node))
    #print(paste0("orig+diff: ",cond_orig_pseudo[which(row.names(cond_cells) %in% kmean_node[[node]])] + diff))
    #print(paste0("orig: ", cond_orig_pseudo[which(row.names(cond_cells) %in% kmean_node[[node]])]))
    #print(kmean_node[[node]])
    #print("")
    #}
    
    #cond_orig_pseudo[which(row.names(cond_cells) %in% kmean_node[[node]])] <- cond_orig_pseudo[which(row.names(cond_cells) %in% kmean_node[[node]])] + diff
    #}
    
    for (node in 1:length(kmean_node)){
      #diff <- cond_node_new[which(cond_node_new[,(dims+2)] == names(kmean_node)[node]), (dims+1)] - cond_node_old[which(row.names(cond_node_old == names(kmean_node)[node]),(dims+1)]
      diff <- cond_node_new[,1][which( row.names(cond_node_new) == names(kmean_node)[node] )] - cond_node_old[,1][which( row.names(cond_node_old) == names(kmean_node)[node] )]
      cond_orig_pseudo[which(cond_cells %in% kmean_node[[node]])] <- cond_orig_pseudo[which(cond_cells %in% kmean_node[[node]])] + diff
      
      if(diff != 0){
        print(which(cond_cells %in% kmean_node[[node]]))
        print(condition)
        print(diff)
      }
    }
    
    #Put new pseudotimes into object 
    if(condition ==1){
      int <- cond1_sce@colData@listData[[paste0("slingPseudotime_", lin1)]]
      int[keep]  <- cond_orig_pseudo
      cond1_sce@colData@listData[[paste0("newslingPseudotime_", lin1)]] <- int
    }
    
    else{
      int <- cond2_sce@colData@listData[[paste0("slingPseudotime_", lin2)]]
      int[keep]  <- cond_orig_pseudo
      cond2_sce@colData@listData[[paste0("newslingPseudotime_", lin2)]] <- int
    }
    kmean_output[[paste0("kmean_condition_",condition)]] <- kmean_node
  }
  output_final <- list(cond1_sce, cond2_sce, kmean_output[[paste0("kmean_condition_",(condition-1))]], kmean_output[[paste0("kmean_condition_",condition)]])
  names(output_final) <- c("condition_1", "condition_2", "kmean_condition_1", "kmean_condition_2")
  return(output_final)
}

#Takes the pseudotimes of the node and gives it a space in the embedding
#depcracated
pseudo_node_to_embed_node <- function(node_pseudotime, cell_pseudotime, cell_embed){
  node_embed <- matrix(nrow = 1, ncol = length(cell_embed[1,]))
  
  for (node in node_pseudotime[,1]){
    previous_dist <- 1000
    
    for (cell in 1:length(cell_pseudotime)){
      dist <- abs(node - cell_pseudotime[cell])
      
      if (dist < previous_dist){
        cell_choice <- row.names(cell_embed)[cell]
        previous_dist <- dist
      }
      
    }
    node_embed <- rbind(node_embed, cell_embed[ which( row.names(cell_embed) == cell_choice ), ])
 
  }
  node_embed <- node_embed[-1,]
  row.names(node_embed) <- row.names(node_pseudotime)
  return(node_embed)
  
}

#assumes matrix has colnames(mtx)[1] == root of T1 and row.names(mtx)[1] == root of T2
#Assumes that the starting point of one of the processes correlates well with somewhere on the other process. I.e. can't align converging processes
#Assumes that the end of one of the processes will correlate well with somewhere on the other trajectory. I.e. can't align diverging processes
#Need to change it so it cuts nodes where the correlation isn't high. Don't remove unmatched nodes from the alignment output, just make them so they have no match in the other column
find_surface_path <- function(error_surface){
  
  choice_string <- c("diag", "vertical", "horizontal")
  
  choice <- c(min(error_surface[,1]), min(error_surface[1,]))
  
  index_choice <- which(choice == min(choice))

  if (index_choice == 1){
    loc <- c(which(error_surface[,1] == choice[index_choice]), 1)
  }
  
  else{
    loc <- c(1, which(error_surface[1,] == choice[index_choice]))
  }
  
  x <- T
  
  match_frame <- data.frame(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]])
  
  diag <- error_surface[(loc[1]+1), (loc[2]+1)]
  vertical <- error_surface[loc[1], (loc[2]+1)]
  horizontal <- error_surface[(loc[1]+1), loc[2]]
  
  options <- c(diag, vertical, horizontal)
  
  move_choice <- which( options == min(options) )
  
  previous_choice = move_choice
  while (x == T){
    #Three movement options
    #1. Diagonally - a 1 to 1 match
    #2 vertically - a many to 1 match 
    #3 horizontally a many to 1 match (but to a different tree from above)
    
    diag <- error_surface[(loc[1]+1), (loc[2]+1)]
    vertical <- error_surface[loc[1], (loc[2]+1)]
    horizontal <- error_surface[(loc[1]+1), loc[2]]
    
    options <- c(diag, vertical, horizontal)
    
    move_choice <- which( options == min(options) )
    
    print(loc)
    print(move_choice)
    print(previous_choice)
    
    #Stops us from creating a massive chain of aligned cells between the two processes
    if ((move_choice + previous_choice) == 5 ){
      print("enter")
      print(options)
      move_choice <- diag
    }
    
    if (move_choice == 1){
      loc <- loc + 1
    }
    
    else if (move_choice == 2){
      loc[2] <- loc[2] + 1
    }
    
    else{
      loc[1] <- loc[1] + 1
    }
    
    print("seperate")
    print(loc)
    print(move_choice)
    print(previous_choice)
    
    print("")
    
    
    
    match_frame <- rbind(match_frame, c(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
   
    previous_choice <- move_choice
    
    if (loc[1] == length(error_surface[,1]) | loc[2] == length(error_surface[1,])){
      x <- F
    }
    
    print("")
    
  }
  return(match_frame)
}

#Trying to sort the problem of big gaps in pseudotime
pseudo_gap_fix <- function(condition_1, condition_2){
  
  pseudo1 <- condition_1[,1]
  pseudo2 <- condition_2[,1]
  
  pseudotime_points <- length(unique(c(pseudo1, pseudo2)))
  
  pseudotime_order_old <- sort(unique(c(pseudo1, pseudo2)))
  
  interval <- max(c(pseudo1, pseudo2)) / (pseudotime_points+1)
  
  pseudotime_order_new <- c(0)
  
  for (point in 1:length(pseudotime_order_old)){
    
    pseudotime_order_new <- append(pseudotime_order_new, (interval + pseudotime_order_new[point]))
    
  }
  
  pseudotime_order_new <- append(pseudotime_order_new, max(c(pseudo1, pseudo2)))
  
  pseudo1_new <- c()
  pseudo2_new <- c()
  
  for (i in pseudo1){
    
    pseudo1_new <- append(pseudo1_new, pseudotime_order_new[which(pseudotime_order_old == i)])
    
  }
  
  for (i in pseudo2){
    pseudo2_new <- append(pseudo2_new, pseudotime_order_new[which(pseudotime_order_old == i)])
    
  }
  
  pseudo1_new <- data.frame(pseudo1_new)
  pseudo2_new <- data.frame(pseudo2_new)
  
  row.names(pseudo1_new) <- row.names(condition_1)
  row.names(pseudo2_new) <- row.names(condition_2)
  
  output <- list(pseudo1_new, pseudo2_new)
  
  names(output) <- c("condition_1", "condition_2")
  
  return(output)
  
}

#Get the most common cell type for the nodes
commonID <- function(captures, clusters){
  
  output <- data.frame(rep("id", length(names(captures))), row.names= names(captures))
  
  for (i in 1:length(captures)){
    cluster_c <- clusters[captures[[i]],]
    
    occurence <- table(cluster_c)
    print(table(cluster_c))
    print(names(occurence)[which(occurence == max(occurence))])
    
    choice <- names(occurence)[which(occurence == max(occurence))]
    
    if (length(choice) > 1){
      choice <- paste0(choice, collapse = "/")
      
    }
    
    output[names(captures)[i],1] <- choice
      
  
    
    print("")
    
    
  }
  
  print(occurence)
  
  return(output)
}

#Creates nodes across the trajectory, 
slingshot_node_maker <- function(sce_obj, slingshot_obj, total_nodes, prefix, dims, id){
  #We don't need to give the nodes a place on the phate map, they just need pseudotimes to do pseudobulk or cellalign
  output <- list()
  
  for (lineage in 1:length(slingshot_obj@curves)){
    print(length(slingshot_obj@curves))
    #pseudotime <- data.frame(sce_obj@colData@rownames,sce_obj@colData@listData[["pseudotime"]])
    pseudotime <- data.frame(sce_obj@colData@listData[[paste0("slingPseudotime_", lineage)]], row.names = sce_obj@colData@rownames)
    pseudotime <- data.frame(pseudotime[ order(pseudotime[,1]), ], row.names = sce_obj@colData@rownames[ order(pseudotime[,1]) ] )
    max_pseudo <- max(pseudotime[,1], na.rm=T)
    
    print(max_pseudo)
    
    #Get pseudotime of nodes
    interval <- max_pseudo / total_nodes[lineage]
    print(interval)
    node_pseudotime_vector <- c(0)
    #Include 0 as a node but not the maximum
    for (node in 2:total_nodes[lineage]){ 
      node_pseudotime_vector <- append(node_pseudotime_vector, (node_pseudotime_vector[(node-1)] + interval))
    }
    
    node_names <- paste0(prefix ,seq(1, total_nodes[lineage], 1) )
    
    #make the tree
    node_tree <- data.frame("Child", "Parent")
    
    #Don't want to make the last node a parent so we'll stop the loop before then
    for (node in 1:(length(node_names)-1)){
      node_tree <- rbind( node_tree , c( node_names[node+1], node_names[node]) )
    }
    
    node_pseudotime <- data.frame(node_pseudotime_vector,row.names=node_names)
    colnames(node_pseudotime) <- "pseudotime"
    
    #Get most common cell type for node
    captures <- cell_capture(pseudotime, node_pseudotime)
    cluster_info <- data.frame(sce_obj@colData@listData[[id]], row.names = sce_obj@colData@rownames)
    node_cluster <- commonID(captures, cluster_info)
    
    colnames(node_tree) <- node_tree[1, ]
    node_tree <- node_tree[-1, ]
    
    output[[paste0("lineage_", lineage, "_ID" )]] <- node_cluster
    output[[paste0("lineage_", lineage, "_tree" )]] <- node_tree
    output[[paste0("lineage_", lineage, "_pseudotime")]] <- node_pseudotime
    
  }
  
  return(output)
}

fill_mtx <- function(penalty_mtx){
  
  choice <- c(min(penalty_mtx[,1]), min(penalty_mtx[1,]))
  
  index_choice <- which(choice == min(choice))
  
  if (index_choice == 1){
    loc <- c(which(penalty_mtx[,1] == choice[index_choice]), 1)
  }
  
  else{
    loc <- c(1, which(penalty_mtx[1,] == choice[index_choice]))
  }
  
  loc <- c(1,1)
  penalty_mtx <- penalty_mtx[ loc[1]:length(penalty_mtx[,1]) , loc[2]:length(penalty_mtx[1,]) ]
  
  traceback_mtx <- matrix(0, nrow = length(row.names(penalty_mtx)), ncol = length(colnames(penalty_mtx)))
  
  traceback_mtx[1,1] <- penalty_mtx[1,1]
  
  traceback_mtx[1,] <- own_fib(penalty_mtx[1,])
  
  traceback_mtx[,1] <- own_fib(penalty_mtx[,1])
  
  for (i in 2:length(traceback_mtx[,1])){
    
    for (j in 2:length(traceback_mtx[1,])){
      choice1 <- traceback_mtx[(i-1), j] + penalty_mtx[i,j]
      choice2 <- traceback_mtx[i,(j-1)] + penalty_mtx[i,j]
      choice3 <- traceback_mtx[(i-1),(j-1)] + penalty_mtx[i,j]
      
      final_choice <- min(c(choice1, choice2, choice3))
      
      traceback_mtx[i,j] <- final_choice
      
    }
    
  }
  row.names(traceback_mtx) <- row.names(penalty_mtx)
  colnames(traceback_mtx) <- colnames(penalty_mtx)
  
  return(traceback_mtx)
  
}

ROC_fill_mtx <- function(penalty_mtx, cut_type){
  
  penalty_mtx_cp <- penalty_mtx
  
  if(cut_type == "minimum"){
    choice <- c(min(penalty_mtx[,1]), min(penalty_mtx[1,]))
    
    index_choice <- which(choice == min(choice))
    
    if (index_choice == 1){
      loc <- c(which(penalty_mtx[,1] == choice[index_choice]), 1)
    }
    
    else{
      loc <- c(1, which(penalty_mtx[1,] == choice[index_choice]))
    }
    
    start <- loc
    
    last_row <- length(penalty_mtx[,1])
    last_column <- length(penalty_mtx[1,])
    
    #Find minimum at end to start from
    choice <- c( min( penalty_mtx[ last_row , ] ), min( penalty_mtx[ , last_column ] ) )

    index_choice <- which(choice == min(choice))
    

    if (index_choice == 1){
      loc <- c( last_row , which( penalty_mtx[ last_row , ] == choice[index_choice] ) )
    }
    
    else{
      loc <- c( which( penalty_mtx[ , last_column ] == choice[index_choice] ), last_column )
    }
    
    end <- loc
    
    output <- find_start_end(penalty_mtx, start, end)
    
    penalty_mtx <- penalty_mtx[ output[[1]][1]:output[[2]][1] , output[[1]][2]:output[[2]][2] ]
    
    
  }
  
  
  else{
    loc <- c(1,1)
    
    penalty_mtx <- penalty_mtx[ loc[1]:length(penalty_mtx[,1]) , loc[2]:length(penalty_mtx[1,]) ]
    
  }
  #Above section is if we want to start the scoring matrix from the lowest point of initial alignment, at the moment this doesn't happen
  
  
  traceback_mtx <- matrix(0, nrow = length(row.names(penalty_mtx)), ncol = length(colnames(penalty_mtx)))
  
  traceback_mtx[1,1] <- penalty_mtx[1,1]
  
  traceback_mtx[1,] <- own_fib_ROC(penalty_mtx[1,])
  
  traceback_mtx[,1] <- own_fib_ROC(penalty_mtx[,1])
  
  for (i in 2:length(traceback_mtx[,1])){
    
    for (j in 2:length(traceback_mtx[1,])){
      choice1 <- traceback_mtx[(i-1), j] + ( (penalty_mtx[i,j] - penalty_mtx[(i-1), j] ))
      choice2 <- traceback_mtx[i,(j-1)] + ((penalty_mtx[i,j] - penalty_mtx[i, (j-1)])  )
      choice3 <- traceback_mtx[(i-1),(j-1)] + (penalty_mtx[i,j] - penalty_mtx[(i-1), (j-1)] )
      
      #choice1 <- traceback_mtx[(i-1), j] + ( (penalty_mtx[i,j]))
      #choice2 <- traceback_mtx[i,(j-1)] + ((penalty_mtx[i,j] )  )
      #choice3 <- traceback_mtx[(i-1),(j-1)] + (penalty_mtx[i,j]  )
      
      final_choice <- min(c(choice1, choice2, choice3))
      
      traceback_mtx[i,j] <- final_choice
      
    }
    
  }
  row.names(traceback_mtx) <- row.names(penalty_mtx)
  colnames(traceback_mtx) <- colnames(penalty_mtx)
  
  return(traceback_mtx)
}

#Find optimal starting point, taking into consideration the correlation landscape of the error surface
find_start_end <- function(penalty_mtx, start, end){
  
  #Need to find where to take the slice
  #determined by finding the maximum index based off the context of the metaCell total of each condition
  index_1 <- which(start == max(start))
  index_1 <- which( c(start[1]/length(penalty_mtx[,1]), start[2]/length(penalty_mtx[1,])) == max(  c(start[1]/length(penalty_mtx[,1]), start[2]/length(penalty_mtx[1,]))  )     )
  print("section")
  print(start)
  print(paste0("index:",index_1))
  
  #get rate of change for the start
  if(index_1 == 1){
    slice_1 <- penalty_mtx[,1]
    
  }
  
  else{
    slice_1 <- penalty_mtx[1,]
  }
  
  #print(slice_1)
  
  distances_1 <- c()
  
  for (i in 1:length(slice_1)){
    distances_1 <- append(distances_1, abs(slice_1[start[index_1]] - slice_1[i]) )
  }
  
  change_v1 <- c()
  
  for(i in 2:length(slice_1)){
    change_v1 <- append(change_v1, abs(slice_1[i]-abs(slice_1[i-1])))
    
  }

  print(distances_1)
  
  #cutoff_1 <- median(change_v1)
  cutoff_1 <- mean(change_v1)
  
  print(cutoff_1)
  
  selection_1 <- slice_1[1:start[index_1]]
  
  #print(which(slice_1 < cutoff_1))
  
  #print(cutoff_1)
  
  start[index_1] <- min(which(distances_1 < cutoff_1))
  
  index_2 <- which(end == min(end))
  index_2 <- which( c(end[1]/length(penalty_mtx[,1]), end[2]/length(penalty_mtx[1,])) == min(  c(end[1]/length(penalty_mtx[,1]), end[2]/length(penalty_mtx[1,]))  )     )
  
  #print(end[index_2])
  
  if(index_2 == 1){
    slice_2 <- penalty_mtx[,length(penalty_mtx[1,])]
    
  }
  
  else{
    slice_2 <- penalty_mtx[length(penalty_mtx[,1]),]
  }
  
  distances_2 <- c()
  
  for (i in 1:length(slice_1)){
    distances_2 <- append(distances_2, abs(slice_2[end[index_2]] - slice_2[i]) )
  }
  
  change_v2 <- c()
  
  for(i in 2:length(slice_2)){
    change_v2 <- append(change_v2, abs(slice_2[i]-abs(slice_2[i-1])))
    
  }
  print("")
  print(distances_2)
  print(median(change_v2))
  print("")
  
  #cutoff_2 <- median(change_v2)
  cutoff_2 <- mean(change_v2)
  
  end[index_2] <- max(which(distances_2 < cutoff_2))
  
  #print(max(which(distances_2 < cutoff_2)))
  
  return(list(start, end))
  
}

#Successively sums together the scores from our penalty matrix to make the initial row and column for our traceback matrix
own_fib <- function(number_vector){
  output <- c(number_vector[1])
  
  for (number in 2:length(number_vector)){
    result <- output[(number-1)] + number_vector[number]
    
    output <- append(output, result)
  }
  
  return(output)
  
}

own_fib_ROC <- function(number_vector){
  output <- c(number_vector[1])
  
  for (number in 2:length(number_vector)){
    result <- output[(number-1)] + (number_vector[number] - number_vector[(number-1)]) 
    
    #result <- output[(number-1)] + (number_vector[number]) 
    
    output <- append(output, result)
  }
  
  return(output)
}

traceforward_pathfind <- function(penalty_mtx, cut_type){
  penalty_mtx_cp <- penalty_mtx
  error_surface <- ROC_fill_mtx(penalty_mtx, cut_type)
  
  penalty_mtx_cut <- penalty_mtx[row.names(error_surface), colnames(error_surface)]
  print("start")
  print(row.names(error_surface))
  print(colnames(error_surface))
  
  #add high score ends to error_surface
  error_surface <- rbind(error_surface, rep(1000,length(colnames(error_surface))))
  error_surface <- cbind(error_surface, rep(1000,length(row.names(error_surface))))
  
  #Find minimum at end to start from
  choice <- c( min( error_surface[ 1 , ] ), min( error_surface[ , 1 ] ) )
  #print(choice)
  
  index_choice <- which(choice == min(choice))
  
  if (index_choice == 1){
    loc <- c( 1 , which( error_surface[ 1 , ] == choice[index_choice] ) )
  }
  
  else{
    loc <- c( which( error_surface[ , 1 ] == choice[index_choice] ), 1 )
  }
  
  loc <- c(1,1)
  
  #diag, vertical, horizontal
  movement_frame <- data.frame(cond_1 = c(1, 1, 0), cond_2 = c(1, 0, 1))

  match_frame <- data.frame(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]])
  
  score_vector <- c(error_surface[loc[1], loc[2]])
  score_vector <- c(penalty_mtx_cut[loc[1], loc[2]])
  
  
  previous_choice <- 1
  
  x <- T
  
  while (x == T){
    choice_vector <- c(error_surface[loc[1]+1, loc[2]+1], error_surface[loc[1]+1 , loc[2]], error_surface[loc[1], loc[2]+1]  )
    
    choice <- which(choice_vector == min(choice_vector))
    # print(loc)
    print(paste0("Choice: ", choice))
    # print(paste0("Previous choice: ", previous_choice))
    print("problem")
    
    if (choice %in% c(2,3)){
      print(loc)
      
      temp_loc <- loc + as.numeric(movement_frame[choice,])
      
      temp_choice_vector <- c(error_surface[temp_loc[1]+1, temp_loc[2]+1], error_surface[temp_loc[1]+1 , temp_loc[2]], error_surface[temp_loc[1], temp_loc[2]+1]  )
      
      temp_choice <- which(temp_choice_vector == min(temp_choice_vector))
      #print(paste0("temp_choice: ", temp_choice))
      if ( !(temp_choice %in% c(choice,1)) ){
        print("in")
        loc <- loc + c(1,1)
        
        choice <- 1
        
        match_frame <- rbind(match_frame,data.frame(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
        
        score_vector <- append(score_vector , c(penalty_mtx_cut[loc[1], loc[2]]) )
        
      }
      
      else{
        loc <- loc + as.numeric(movement_frame[choice,])
        
        match_frame <- rbind(match_frame,data.frame(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
        
        score_vector <- append(score_vector , c(penalty_mtx_cut[loc[1], loc[2]]) )
        
      }
      
    }
    
    else{
      loc <- loc + as.numeric(movement_frame[choice,])
      
      match_frame <- rbind(match_frame,data.frame(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
      
      score_vector <- append(score_vector , c(penalty_mtx_cut[loc[1], loc[2]]) )
      
      
    }
    
    #Exit when we reach the end of one of the processes
    if (loc[1] == (length(error_surface[,1])-1) | loc[2] == (length(error_surface[1,])-1) ){
      x <- F
    }
    
    #print(loc[2] == length(error_surface[1,]))
    
    previous_choice <- choice
    #print("")
  }
  

  #add the steps that take us to the final matrix position
  finish_choice <- c( abs(loc[1] - length(row.names(error_surface))  )-1  , abs(loc[2] - length(colnames(error_surface)) )  -1 ) 
  
  finish_index <- which(finish_choice == max(finish_choice))
  
  loc[finish_index] <- loc[finish_index] + 1
  
  finish_choice[finish_index] <- finish_choice[finish_index] - 1
  
  last_index_1 <- seq(loc[1], loc[1] + finish_choice[1])
  last_index_2 <- seq(loc[2], loc[2] + finish_choice[2])
  
  print(last_index_1)
  print(last_index_2)
  
  print(finish_choice)
  
  finish_choice <- finish_choice + 1
  
  print(finish_choice)
  
  final_addition <- data.frame(rep(row.names(error_surface)[last_index_1], finish_choice[2]), rep(colnames(error_surface)[last_index_2], finish_choice[1]) )
  colnames(final_addition) <- colnames(match_frame)
  print(final_addition)
  
  print(colnames(match_frame))
  
  match_frame <- rbind(match_frame,final_addition)
  
  print(penalty_mtx_cut[last_index_1, last_index_2])
  
  score_vector <- append(score_vector , c(penalty_mtx_cut[last_index_1, last_index_2]) )
  
  
  output <- data.frame('X' = match_frame[,2],"Y" = match_frame[,1] , score_vector)
  
  return(output)
  
  
}


#Traceback from the lowest point on the last row or column
#If you do that it follows the path of just moving through the minimum like in my first way of finding the path
traceback_pathfind <- function(penalty_mtx, cut_type){
  
  penalty_mtx_cp <- penalty_mtx
  error_surface <- ROC_fill_mtx(penalty_mtx, cut_type)
  
  
  last_row <- length(error_surface[,1])
  last_column <- length(error_surface[1,])
    
  factor <- -1 
  
  penalty_mtx_cut <- penalty_mtx[row.names(error_surface), colnames(error_surface)]
  
  
  #Find minimum at end to start from
  choice <- c( min( error_surface[ last_row , ] ), min( error_surface[ , last_column ] ) )
  print(choice)
  
  index_choice <- which(choice == min(choice))
  
  if (index_choice == 1){
    loc <- c( last_row , which( error_surface[ last_row , ] == choice[index_choice] ) )
  }
  
  else{
    loc <- c( which( error_surface[ , last_column ] == choice[index_choice] ), last_column )
  }
  
  loc <- c(last_row, last_column)
  
  match_frame <- data.frame(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]])
  
  score_vector <- c(penalty_mtx_cut[loc[1], loc[2]])
  
  diag <- error_surface[(loc[1]+factor), (loc[2]+factor)]
  vertical <- error_surface[(loc[1]+factor), loc[2]]
  horizontal <- error_surface[loc[1], (loc[2]+factor)]
  
  options <- c(diag, vertical, horizontal)
  
  move_choice <- which( options == min(options) )
  
  print(move_choice)
  
  previous_choice = move_choice
  
  if (1 == move_choice){
    loc <- loc + factor
    print("diag")
  }
  
  else if (2 == move_choice){
    loc[1] <- loc[1] + factor
    print("vertical")
  }
  
  else{
    loc[2] <- loc[2] + factor
    #print("horizontal")
  }
  
  #print(paste("match: ", c(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]])))
  
  match_frame <- rbind(match_frame, c(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
  score_vector <- append(score_vector,penalty_mtx_cut[loc[1], loc[2]])
  
  x <- T
  counter <- 1
  while (x == T){
    #print(loc)
    #Three movement options
    #1. Diagonally - a 1 to 1 match
    #2 vertically - a many to 1 match 
    #3 horizontally a many to 1 match (but to a different tree from above)
    
    diag <- error_surface[(loc[1]+factor), (loc[2]+factor)]
    vertical <- error_surface[(loc[1]+factor), loc[2]]
    horizontal <- error_surface[loc[1], (loc[2]+factor)]
    
    options <- c(diag, vertical, horizontal)
    
    #print(paste("diag: ", diag))
    #print(paste("horizontal: ", horizontal))
    #print(paste("vertical: ", vertical))
    #print(options)
    
    move_choice <- which( options == min(options, na.rm=T) )
    
    #print(paste("move_choice:" , c("diag", "vertical", "horizontal")[move_choice]))
    #print("")
    
    #Stops us from creating a massive chain of aligned cells between the two processes
    if ((move_choice + previous_choice) == 5 ){
      
      #We want to move diagonally, but not for the next move, but for the move we just made, so we retroactively change the previous result.
      match_frame <- match_frame[-length(match_frame[,1]),]
      score_vector <- score_vector[-length(score_vector)]
      
      if(move_choice == 2){
        loc[2] <- loc[2] - factor
      }
      
      else{
        loc[1] <- loc[1] - factor
      }
      
      move_choice <- 1
      
    }
    
    if (1 == move_choice){
      loc <- loc + factor
      #print("diag")
    }
    
    else if (2 == move_choice){
      loc[1] <- loc[1] + factor
      #print("vertical")
    }
    
    else{
      loc[2] <- loc[2] + factor
      #print("horizontal")
    }
    
    #print(paste("match: ", c(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]])))
    
    match_frame <- rbind(match_frame, c(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
    score_vector <- append(score_vector,penalty_mtx_cut[loc[1], loc[2]])
    
    previous_choice <- move_choice
    
    #print("")
    
    if (loc[1] == 2 | loc[2] == 2){
      print(paste("final: ", c(loc[1], loc[2])))
      
      #If we've reached the end of the alignment already then we don't want to add it again
      if(sum(loc) != 2){
        
        loc[1] <- loc[1] + factor
        loc[2] <- loc[2] + factor
        
        final_addition <- cbind(rep(row.names(error_surface)[loc[1]:1], loc[2]) , rep(colnames(error_surface)[loc[2]:1], loc[1]))
        score_vector <- append(score_vector,penalty_mtx_cut[loc[1]:1, loc[2]:1])
        
        print(final_addition)
        
        colnames(final_addition) <- colnames(match_frame)
        
        match_frame <- rbind(match_frame, final_addition)
      }
      
      
      x <- F
    }
    #print(options)
    counter <- counter + 1
    
  }
  
  output <- data.frame(rev(match_frame[,2]), rev(match_frame[,1]), rev(score_vector))
  colnames(output) <- c("X", "Y")
  
  return(output)
  
}


pathfind <- function(penalty_mtx, cut_type){
  choice_list <- list(traceforward_pathfind(penalty_mtx, cut_type), traceback_pathfind(penalty_mtx, cut_type) )
  
  scores <- c(sum(choice_list[[1]][,3]), sum(choice_list[[2]][,3]))
  
  #return whichever path has the lowest overall score
  
  output <- choice_list[[which(scores == min(scores))]]
  
  return(output)
}

#Finds the median of the matched and unmatched nodes. Then calculates the difference between each of the matched pairs and these two thresholds
#If the match is closer to median of the unmatched nodes than the matched nodes, then we designate it to be 'cut' from the alignment
cut_deviate <- function(alignment, dis_mtx){
  
  threshold_output <- get_thresholds(alignment, dis_mtx)
  
  median_match <- threshold_output[[1]]
  median_unmatch <- threshold_output[[2]]
  match_cost <- threshold_output[[3]]
  
  
  score_frame <- data.frame( match_median = abs(match_cost[["score"]] - median_match), unmatch_median = abs(match_cost[["score"]] - median_unmatch) )
  
  final_match <- cbind(alignment, rep("match", length(alignment[,1])))
  
  final_match[which(score_frame[,1] > score_frame[,2]),4] <- "cut"
  
  colnames(final_match) <- c("X", "Y", "Score","Status")
  
  #return(score_frame)
  return(final_match)
  
}

#Gets the thresholds for the cut_deviate() function
get_thresholds <- function(alignment, penalty_mtx){
  match_cost <- c()
  
  for (match in 1:length(alignment[,1])){
    loc <- c(which(row.names(penalty_mtx) == alignment[match,2]), which(colnames(penalty_mtx) == alignment[match,1]))
    match_cost <- append(match_cost, penalty_mtx[loc[1], loc[2]])
    
  }
  match_cost <- as.matrix(match_cost)
  
  match_cost <- data.frame(score = match_cost, group = rep("match", length(match_cost)))
  
  match_cost_unmatched <- penalty_mtx
  
  rows <- as.numeric(str_replace(alignment[,2], "Y", ""))
  cols <- as.numeric(str_replace(alignment[,1], "X", ""))
  
  coord <- as.matrix(data.frame(rows, cols))
  
  match_cost_unmatched[coord] <- NA
  
  match_cost_unmatched <- as.vector(match_cost_unmatched)
  match_cost_unmatched <- data.frame(score = match_cost_unmatched[which(is.na(match_cost_unmatched) == F )], group = rep("unmatch", length(match_cost_unmatched[which(is.na(match_cost_unmatched) == F )])))
  
  median_match <- median(match_cost[["score"]])
  median_unmatch <- median(match_cost_unmatched[["score"]])
  
  return(list(median_match, median_unmatch, match_cost))
}

#Permutation t-test function
permutation_test <- function(group_1, group_2){
  combination <- c(group_1, group_2)
  
  print(length(combination))
  
  permutations <- length(combination)/2
  
  initial <- mean(group_1) - mean(group_2)
  
  perm_test <- c()
  
  for (i in 1:permutations){
    choice <- sample(1:2, 1)
    
    if(choice == 1){
      split_n <- floor(permutations)
    }
    
    else{
      split_n <- ceiling(permutations)
    }
    
    sample_1 <- sample(combination, split_n)
    
    print(sample_1)
    
    perm_1 <- combination[which(combination %in% sample_1)]
    perm_2 <- combination[which(!(combination %in% sample_1))]
    
    perm_test <- append(perm_test,( mean(perm_1) - mean(perm_2) ) )
    
  }
  
  result <- length(which(perm_test >= initial)) / permutations
  
  print(result)
  
  return(list(initial, perm_test))
  
}

find_node_position <- function(cell_embedding, cell_pseudo, node_pseudo){
  node_embedding <- matrix(nrow = 1, ncol=length(cell_embedding[1,]))
  colnames(node_embedding) <- colnames(cell_embedding)
  
  for (node in node_pseudo[,1]){
    choice <- row.names(cell_pseudo)[which( abs(cell_pseudo-node) == min( abs(cell_pseudo-node) ) )]
    node_embedding<- rbind(node_embedding, cell_embedding[choice,])
  }
  node_embedding <- as.data.frame(node_embedding[-1,])
  node_embedding <- cbind(node_embedding, rep("node", length(row.names(node_embedding))))
  row.names(node_embedding) <- row.names(node_pseudo)
  colnames(node_embedding) <- c(colnames(cell_embedding), "type")
  
  return(node_embedding)
  
}

full_alignment_function <- function(alignment, penalty_mtx){
  alignment_tracker <- list()
  
  penalty_mtx_orig <- penalty_mtx
  
  path <- traceback_pathfind(penalty_mtx)
  
  threshold_output <- get_thresholds(alignment, penalty_mtx)
  
  median_match <- threshold_output[[1]]
  median_unmatch <- threshold_output[[2]]
  
  cut_path <- cut_deviate(path, penalty_mtx)
  
  index <- which(cut_path[,3] == "cut")
  
  rows <- as.numeric(str_replace(alignment[index,2], "Y", ""))
  cols <- as.numeric(str_replace(alignment[index,1], "X", ""))
  
  coord <- as.matrix(data.frame(rows, cols))
  
  penalty_mtx[coord] <- penalty_mtx[coord] + 0.5
  
  path <- traceback_pathfind(penalty_mtx)
  
  previous_cuts <- length(index)
  
  alignment_tracker[[previous_cuts]] <- path
  
  #return(penalty_mtx)
  
  counter <- 1
  chance <- 5
  while(chance != 0){
    
    cut_path <- cut_deviate(path, penalty_mtx)
    
    index <- which(cut_path[,3] == "cut")
    
    rows <- as.numeric(str_replace(path[index,2], "Y", ""))
    cols <- as.numeric(str_replace(path[index,1], "X", ""))
    
    coord <- as.matrix(data.frame(rows, cols))
    
    penalty_mtx[coord] <- penalty_mtx[coord] + 0.5
    
    print(penalty_mtx[coord])
    
    path <- traceback_pathfind(penalty_mtx)
    
    if (length(index) > previous_cuts){
      chance <- chance - 1
    }
    
    previous_cuts <- length(index)
    
    alignment_tracker[[paste(counter, " - ", length(index))]] <- path
    
    print(previous_cuts)
    
    #chance <- chance -1
    counter <- counter + 1
    
  }
  
  return(list(alignment_tracker, penalty_mtx))
  
}

align_nodes_linear_1_to_1 <- function(alignment, condition_1, condition_2){
  new_pseudo1 <- condition_1
  new_pseudo2 <- condition_2
  
  #As we're assuming that the starting node of either cond_1 or cond2 must be the starting point for the alignment, we need to adjust based on this
  
  #Initial aligning point is at start of condition 1
  if (alignment[1,1] == "X1"){
    new_root_pseudotime <- condition_2[alignment[1,2],]
    
    #as pseudotime of the first node is 0 we can just add the new_root_pseudotime to get our shunted pseudotime
    new_pseudo1 <- new_pseudo1 + new_root_pseudotime
  }
  
  
  else{
    new_root_pseudotime <- condition_1[alignment[1,1],]
    
    new_pseudo2 <- new_pseudo2 + new_root_pseudotime
  }
  
  for (match in 1:length(alignment[,1])){
    
    choices <- c(new_pseudo1[alignment[match,1],1], new_pseudo2[alignment[match,2],1])
    
    min_choice <- which(choices == min(choices))
    
    if (min_choice == 1){
      diff <- new_pseudo2[alignment[match,2],1] - new_pseudo1[alignment[match,1],1]
      
      pulldown_index <- which(row.names(new_pseudo2) == alignment[match,2])
      
      print(pulldown_index)
      print(length(new_pseudo2[,1]))
      
      new_pseudo2[pulldown_index:length(new_pseudo2[,1]),1] <- new_pseudo2[pulldown_index:length(new_pseudo2[,1]),1] - diff
    
    }
    
    else{
      diff <- new_pseudo1[alignment[match,1],1] - new_pseudo2[alignment[match,2],1]
      
      pulldown_index <- which(row.names(new_pseudo1) == alignment[match,1])

      new_pseudo1[pulldown_index:length(new_pseudo1[,1]),1] <- new_pseudo1[pulldown_index:length(new_pseudo1[,1]),1] - diff
      
      
    }
    
    
  }
  
  row.names(new_pseudo1) <- row.names(condition_1)
  row.names(new_pseudo2) <- row.names(condition_2)
  
  output <- gap_fix(new_pseudo1, new_pseudo2)
  return(output)
  
}

align_nodes_linear_seperate <- function(alignment, condition_1, condition_2){
  cut_sequences <- paste(alignment[,3], collapse = "")
  
  cut_sequences <- str_split(cut_sequences, "match")
  
  cut_sequences <- cut_sequences[[1]][which(cut_sequences[[1]] != "")]
  length_cuts <- c()
  
  for (i in cut_sequences){
    length_cuts <- append(length_cuts, (length(str_split(i, "t")[[1]])-1))
  }
  
  new_pseudo1 <- condition_1
  new_pseudo2 <- condition_2
  
  #As we're assuming that the starting node of either cond_1 or cond2 must be the starting point for the alignment, we need to adjust based on this
  
  #Initial aligning point is at start of condition 1
  if (alignment[1,1] == "X1"){
    new_root_pseudotime <- condition_1[alignment[1,1],]
    
    #as pseudotime of the first node is 0 we can just add the new_root_pseudotime to get our shunted pseudotime
    new_pseudo1 <- new_pseudo1 + new_root_pseudotime
  }
  
  
  else{
    new_root_pseudotime <- condition_2[alignment[1,2],]
    
    new_pseudo2 <- new_pseudo2 + new_root_pseudotime
  }
  
  cut_counter <- 1
  
  for (match in 1:length(alignment[,1])){
    
    choices <- c(new_pseudo1[alignment[match,1],1], new_pseudo2[alignment[match,2],1])
    
    min_choice <- which(choices == min(choices))
    
    if(alignment[match,3] == "cut"){
      start1 <- new_pseudo1[alignment[match,1],1]
      end1 <- new_pseudo1[alignment[(match+length_cuts[cut_counter]),1],1]
      
      start2 <- new_pseudo2[alignment[match,2],1]
      end2 <- new_pseudo2[alignment[(match+length_cuts[cut_counter]),2],1]
      
      match <- match + length_cuts[cut_counter]
      cut_counter <- cut_counter + 1
    }
    
    if (min_choice == 1){
      diff <- new_pseudo2[alignment[match,2],1] - new_pseudo1[alignment[match,1],1]
      
      pulldown_index <- which(names(new_pseudo2) == alignment[match,2])
      
      new_pseudo2[pulldown_index:length(new_pseudo2[,1]),1] <- new_pseudo2[pulldown_index:length(new_pseudo2[,1]),1] - diff
      
    }
    
    else{
      diff <- new_pseudo1[alignment[match,1],1] - new_pseudo2[alignment[match,2],1]
      
      pulldown_index <- which(names(new_pseudo1) == alignment[match,1])
      
      new_pseudo1[pulldown_index:length(new_pseudo1[,1]),1] <- new_pseudo1[pulldown_index:length(new_pseudo1[,1]),1] - diff
      
      
    }
    
    
  }
  
  
}

align_nodes_linear_remove <- function(alignment, condition_1, condition_2){}

error_graph_node <- function(alignment, aligned_pseudo1, aligned_pseudo2){
  
  previous_pseudopoint <- aligned_pseudo1[alignment[1,1],1]
  
  output_df <- data.frame(score = 0, pseudotime = 0)
  score_vector <- c(alignment[1,3])
  
  for (align in 2:length(alignment[,3])){
    
    pseudo_point <- aligned_pseudo1[alignment[align,1],1]
    
    print(pseudo_point)
    
    if(pseudo_point != previous_pseudopoint){
      print("enter")
      score <- mean(score_vector)
      score_vector <- c()
      output_df <- rbind(output_df, data.frame(score=score, pseudotime = previous_pseudopoint))
      previous_pseudopoint <- pseudo_point
      
    }
    
    score_vector <- append(score_vector, alignment[align,3])
    
  }
  
  score <- mean(score_vector)
  output_df <- rbind(output_df, data.frame(score=score, pseudotime = pseudo_point))
  
  return(output_df)
}

#Finds the 'error' (1-correlation) for each of the nodes in our alignment. i.e. how correlated the gene expression of each of the nodes in the match is
error_node <- function(alignment){
  error_list <- list()
  
  for (col in 1:2){
    error_df <- data.frame(node = "X", score = 0)
    
    for ( node in unique( alignment[,col] ) ){
      #print(which( alignment[,col] == node))
      
      error <- mean(alignment[ which( alignment[,col] == node) , 3 ])
      
      #print(error)
      
      error_df <- rbind(error_df, data.frame(node = node, score = error))
      
    }
    
    error_df <- error_df[-1,]
    
    error_list[[paste("Condition_", col)]] <- error_df
    
  }
  
  return(error_list)
  
}

#Using the error of the nodes, we can estimate the error of the individual cells in the dataset, 
#based on their distance (in terms of pseudotime) to the nodes
error_cells <- function(node_score, node_pseudo, cell_pseudo, window_size){
  output_list <- list()
  
  for (i in 1:2){
    cell_pseudo_current <- cell_pseudo[[i]]
    node_pseudo_current <- node_pseudo[[i]]
    node_score_current <- node_score[[i]]
    
    #remove nodes that have never been matched
    node_pseudo_current <- node_pseudo_current[ which(node_pseudo_current[,1] %in% node_score_current[,1]) , ]
    
    node_mtx <- t( replicate(length(cell_pseudo_current[,2]), node_pseudo_current[,2]) )
    cell_mtx <- replicate( length(node_pseudo_current[,2]) , cell_pseudo_current[,2])
    
    #remove cells that are outside the window (in terms of pseudotime) from our aligned nodes
    dist_filter <- (node_mtx - cell_mtx)^2
    
    avg_weight <- rowMins(dist_filter)
    
    #return(avg_weight)
    
    #print(avg_weight)
  
    #weight. N x C matrix 
    weight <-  exp( -1 * ( (node_mtx - cell_mtx)^2 / (window_size)^2 )) 
    
    row.names(weight) <- cell_pseudo_current[,1]
    
    cells_na <- row.names(weight[ which(avg_weight > ((window_size)^2) ) , ])

    #NxC dot 1xN
    
    #print(t(node_score_current[,2]))
    
    #return(weight)
    
    score <- weight %*% node_score_current[,2]
    normalizer <- 1/rowSums(weight)
    score <- score * normalizer
    row.names(score) <- row.names(weight)
    
    score <- cellAlign_function(node_pseudo_current, node_score_current, cell_pseudo_current, window_size)
    
    score[which(row.names(score) %in% cells_na)] <- NA
    
    output_list[[paste("condition_", i)]] <- score
  }
  
  return(output_list)
  
}

#Function which estimates variable 2 for set2 based off of variable 2 for set 1 and weighted by the distance between variable 1 for set 1 and set 2.
cellAlign_function <- function(set1_variable1, set1_variable2, set2_variable1 ,window_size){
  set1_variable1_mtx <- t( replicate(length(set2_variable1[,1]), set1_variable1[,1]) )
  set2_variable1_mtx <- replicate( length(set1_variable1[,1]) , set2_variable1[,1])
  
  #weight. N x C matrix 
  weight <-  exp( -1 * ( (set1_variable1_mtx - set2_variable1_mtx)^2 / ((window_size)^2) )) 
  
  row.names(weight) <- row.names(set2_variable1)
  
  #NxC dot 1xN
  
  set2_variable2 <- weight %*% set1_variable2[,2]
  
  normalizer <- 1/rowSums(weight)
  
  set2_variable2 <- set2_variable2 * normalizer
  row.names(set2_variable2) <- row.names(weight)
  
  return(set2_variable2)
}


#Finds equivalent output for nodes based on distance in values between the nodes and cells
node_estimator <- function(node_values, cell_values, cell_output){
  
  node_mtx <- matrix(rep(node_values[,2], each = length(cell_values[,2])), ncol=length(cell_values[,2]), byrow=T)
  
  cell_mtx <- matrix(rep(cell_values[,2], each = length(node_values[,2])), nrow = length(node_values[,2]))
  diff_mtx <- abs(node_mtx - cell_mtx)
  
  row.names(diff_mtx) <- row.names(node_values)
  colnames(diff_mtx) <- cell_values[,1]
  
  #return(diff_mtx)
  
  minimun_row <- rowMin(diff_mtx) 
  
  node_output <- matrix(nrow=1, ncol=length(colnames(cell_output)))
  
  for (row in 1:length(diff_mtx[,1])){
    
    cell_choice <- colnames(diff_mtx)[which(diff_mtx[row,] == minimun_row[row])]
    
    if(length(cell_choice) > 1){
      cell_choice <- cell_choice[1]
    }
    
    node_output <- rbind(node_output, cell_output[cell_choice,])
    
  }
  
  node_output <- node_output[-1,]
  return(node_output)
  
}

#Finds equivalent output for cells based on distance in values between the nodes and cells
cell_estimator <- function(node_values, cell_values, node_output){
  
  node_mtx <- matrix(rep(node_values[["pseudotime"]], each = length(cell_values[["pseudotime"]])), ncol=length(cell_values[["pseudotime"]]), byrow=T)
  
  cell_mtx <- matrix(rep(cell_values[["pseudotime"]], each = length(node_values[["pseudotime"]])), nrow = length(node_values[["pseudotime"]]))
  diff_mtx <- abs(node_mtx - cell_mtx)
  
  row.names(diff_mtx) <- row.names(node_values)
  colnames(diff_mtx) <- cell_values[["cell_ID"]]
  
  #return(diff_mtx)
  
  minimun_row <- rowMin(diff_mtx) 
  
  node_output <- matrix(nrow=1, ncol=length(colnames(cell_output)))
  
  for (row in 1:length(diff_mtx[,1])){
    
    cell_choice <- colnames(diff_mtx)[which(diff_mtx[row,] == minimun_row[row])]
    
    if(length(cell_choice) > 1){
      cell_choice <- cell_choice[1]
    }
    
    node_output <- rbind(node_output, cell_output[cell_choice,])
    
  }
  
  node_output <- node_output[-1,]
  return(node_output)
}

#Not been tested on alignments with gaps 
#Nodes that come after the final alignment won't have their pseudotimes pulled down
pseudo_align_multi_match <- function(condition_1, condition_2,  alignment, dims){
  new_pseudo1 <- condition_1
  new_pseudo2 <- condition_2
  
  #As we're assuming that the starting node of either cond_1 or cond2 must be the starting point for the alignment, we need to adjust based on this
  
  #Initial aligning point is at start of condition 1
  if (alignment[1,1] == "X1"){
    print("in")
    new_root_pseudotime <- condition_2[alignment[1,2],]
    
    #as pseudotime of the first node is 0 we can just add the new_root_pseudotime to get our shunted pseudotime
    new_pseudo1 <- new_pseudo1 + new_root_pseudotime
    
    print(new_pseudo1[1,1])
    print(new_pseudo2[1,1])
  }
  
  
  else{
    new_root_pseudotime <- condition_1[alignment[1,1],]
    
    new_pseudo2 <- new_pseudo2 + new_root_pseudotime
  }
  
  #Work from end of the alignment backwards, changing the pseudotime of nodes 
  #Need to account for unaligned nodes
  between_nodes <- c(0,0)
  for (align in 1:length(alignment[,1])){
    
    if ( alignment[align,2] == "cut"){
      
      loc <- which(alignment[align,3] != "cut")
      
      if (loc == 1){
        between_nodes[1] <- between_nodes[1] + 1
      }
      
      else{
        between_nodes[2] <- between_nodes[2] + 1
      }
      
    }
    
    else{
      pseudopoint_1 <- new_pseudo1[alignment[align,1],]
      pseudopoint_2 <- new_pseudo2[alignment[align,2],]
      
      #Whichever node has the highest pseudotime should be condensed down to the other
      choice <- which( c(pseudopoint_1, pseudopoint_2) == max( c(pseudopoint_1, pseudopoint_2) ) )
      #print(pseudopoint_1)
      #print(pseudopoint_2)
      #print(choice)
      #print("")
      
      if (choice == 1){
        #print("problem")
        new_pseudo1[alignment[align,1],] <- pseudopoint_2
        
        #difference <- pseudopoint_1 - pseudopoint_2
        
        #print(new_pseudo1[ which(row.names(new_pseudo1) == alignment[align,1]):length(new_pseudo1[,1]) , ])
        
        #new_pseudo1[ which(row.names(new_pseudo1) == alignment[align,1]):length(new_pseudo1[,1]) , ] <- new_pseudo1[ which(row.names(new_pseudo1) == alignment[align,1]):length(new_pseudo1[,1]) , ] - difference
        
        
        #Squish intervening nodes
        if(sum(between_nodes) != 0){
          squish_min <- new_pseudo1[alignment[ ( align - between_nodes[1] ) , 1]]
          squish_max <- pseudopoint_2
          
          #If we're going backwards use squish max, if it's forwards from root use squish min
          cond1_squish <- c(squish_max)
          cond2_squish <- c(squish_max)
          
          cond1_interval <- (squish_max-squish_min) / (between_nodes[1] +1)
          cond1_interval <- (squish_max-squish_min) / (between_nodes[2] +1)
          
          
          for (i in 1:between_nodes[1]){cond1_squish <- append(cond1_squish, (cond1_squish[i] + cond1_interval))}
          
          for (i in 1:between_nodes[2]){cond2_squish <- append(cond2_squish, (cond1_squish[i] + cond2_interval))}
          
          if (between_nodes[1] > 0){
            new_pseudo1[(align - between_nodes[1]):align] <- cond1_squish[2:length(cond1_squish)]
          }
          
          
          if (between_nodes[2] > 0){
            new_pseudo2[(align - between_nodes[2]):align] <- cond2_squish[2:length(cond2_squish)]
            
          }
          
          
          between_nodes <- c(0,0)
        }
        
      }
      
      else{
        new_pseudo2[alignment[align,2],] <- pseudopoint_1
        
        #difference <- pseudopoint_2 - pseudopoint_1
        
        #new_pseudo2[ which(row.names(new_pseudo2) == alignment[align,2]):length(new_pseudo2[,1]) , ] <- new_pseudo2[ which(row.names(new_pseudo2) == alignment[align,2]):length(new_pseudo2[,1]) , ] - difference
        ƒnode
        #Squish intervening nodes
        if(sum(between_nodes) != 0){
          squish_min <- new_pseudo1[alignment[ ( align - between_nodes[1] ) , 1]]
          squish_max <- pseudopoint_2
          
          cond1_squish <- c(squish_min)
          cond2_squish <- c(squish_min)
          
          cond1_interval <- (squish_max-squish_min) / (between_nodes[1] +1)
          cond1_interval <- (squish_max-squish_min) / (between_nodes[2] +1)
          
          
          for (i in 1:between_nodes[1]){cond1_squish <- append(cond1_squish, (cond1_squish[i] + cond1_interval))}
          
          for (i in 1:between_nodes[2]){cond2_squish <- append(cond2_squish, (cond1_squish[i] + cond2_interval))}
          
          if (between_nodes[1] > 0){
            new_pseudo1[(align - between_nodes[1]):align] <- cond1_squish[2:length(cond1_squish)]
          }
          
          
          if (between_nodes[2] > 0){
            new_pseudo2[(align - between_nodes[2]):align] <- cond2_squish[2:length(cond2_squish)]
            
          }
          
          
          between_nodes <- c(0,0)
        }
        
      }
      
    }
    
  }
  
  row.names(new_pseudo1) <- row.names(condition_1)
  row.names(new_pseudo2) <- row.names(condition_2)
  
  output <- pseudo_gap_fix(new_pseudo1, new_pseudo2)
  return(output)
}

#Aligns nodes, with initial alignment being the point of first match
pseudo_align_nodes <- function(condition_1, condition_2, alignment){
  
  alignment_cp <- alignment
  
  new_pseudo1 <- condition_1
  new_pseudo2 <- condition_2
  
  alignment_initial <- subset(alignment, Status == "match")

  alignment <- alignment[min(which(alignment[,4] == "match")):length(alignment[,1]),]
  
  unmatched <- subset(alignment_cp, Status == "cut")
  unmatched <- unique( c(unmatched[,1], unmatched[,2])  )
  
  unmatched <- unmatched[which( !(unmatched %in%  alignment_initial[,1]))]
  unmatched <- unmatched[which( !(unmatched %in%  alignment_initial[,2]))]
  
  unmatched <- append(unmatched, row.names(new_pseudo1)[ which( !( row.names(new_pseudo1) %in% alignment_initial[,1] ) ) ])
  unmatched <- append(unmatched, row.names(new_pseudo2)[ which( !( row.names(new_pseudo2) %in% alignment_initial[,2] ) ) ])
  
  cond1_align <- data.frame(node = "node", pseudotime = 0)
  cond1_noalign <- data.frame(node = "node", pseudotime = 0)
  cond2_align <- data.frame(node = "node", pseudotime = 0)
  cond2_noalign <- data.frame(node = "node", pseudotime = 0)
  
  choices <- c(as.numeric(str_replace(alignment_initial[1,1], "X", "")) , as.numeric(str_replace(alignment_initial[1,2], "Y", "")) ) 
  max_choice <- which(choices == max(choices))
  
  print(new_pseudo1[1,1])
  print(new_pseudo2[1,1])
 
  print(max_choice)
  
  #Initial aligning point is at start of condition 1
  if (max_choice == 1){
    new_root_pseudotime <- condition_1[alignment_initial[1,1],]
    diff <- abs(new_root_pseudotime - condition_2[alignment_initial[1,2],])
    
    #as pseudotime of the first node is 0 we can just add the new_root_pseudotime to get our shunted pseudotime
    new_pseudo2 <- new_pseudo2 + diff
    
    print(new_pseudo1[1,1])
    print(new_pseudo2[1,1])
  }
  
  else{
    new_root_pseudotime <- condition_2[alignment_initial[1,2],]
    diff <- abs(new_root_pseudotime - condition_1[alignment_initial[1,1],])
    
    #as pseudotime of the first node is 0 we can just add the new_root_pseudotime to get our shunted pseudotime
    new_pseudo1 <- new_pseudo1 + diff
    
    print(new_pseudo1[1,1])
    print(new_pseudo2[1,1])
    
  }
  
  print("help")
  
  
  for (match in 2:length(alignment[,1])){
    choices <- c( new_pseudo1[alignment[match, 1],], new_pseudo2[alignment[match,2], ] )
    
    min_choice <- which(choices == min(choices))
    
    if (min_choice == 1){
      diff <- max(choices) - min(choices)
      
      node_change <- which(row.names(new_pseudo2) == alignment[match,2])
  
      new_pseudo2[node_change:length(new_pseudo2[,1]),1] <- new_pseudo2[node_change:length(new_pseudo2[,1]),1] - diff
      
    }
    
    else{
      
      diff <- max(choices) - min(choices)
      
      node_change <- which(row.names(new_pseudo1) == alignment[match,1])
      
      new_pseudo1[node_change:length(new_pseudo1[,1]),1] <- new_pseudo1[node_change:length(new_pseudo1[,1]),1] - diff
      
    }
    
    
  }
  
  recursion_list <- list(new_pseudo1, new_pseudo2)
  output <- list()
  for (i in 1:2){
    
    obj <- recursion_list[[i]]
    vector <- c()
    
    for (node in row.names(obj)){
      
      if (node %in% unmatched){
        #vector <- append(vector, paste0("condition_",i,"_noalign"))
        vector <- append(vector, "noalign")
      }
      
      else{
        #vector <- append(vector, paste0("condition_",i,"_align"))
        vector <- append(vector, "align")
      }
      
    }
    #levels_vector <- c( "condition_1_noalign", "condition_1_align", "condition_2_align", "condition_2_noalign" )
    #vector <- factor(vector, levels = levels_vector)
    output[[paste0("condition_", i)]] <- cbind(recursion_list[[i]], vector)
    
  }
  return(output)
  
}

pseudo_align_nodes_squish <- function(condition_1, condition_2, alignment){
  alignment_cp <- alignment
  
  new_pseudo1 <- condition_1
  new_pseudo2 <- condition_2
  
  alignment_initial <- subset(alignment, Status == "match")
  alignment <- alignment[min(which(alignment[,4] == "match")):length(alignment[,1]),]
  
  unmatched <- subset(alignment_cp, Status == "cut")
  unmatched <- unique( c(unmatched[,1], unmatched[,2])  )
  
  unmatched <- unmatched[which( !(unmatched %in%  alignment_initial[,1]))]
  unmatched <- unmatched[which( !(unmatched %in%  alignment_initial[,2]))]
  
  unmatched <- append(unmatched, row.names(new_pseudo1)[ which( !( row.names(new_pseudo1) %in% alignment_initial[,1] ) ) ])
  unmatched <- append(unmatched, row.names(new_pseudo2)[ which( !( row.names(new_pseudo2) %in% alignment_initial[,2] ) ) ])
  
  cond1_align <- data.frame(node = "node", pseudotime = 0)
  cond1_noalign <- data.frame(node = "node", pseudotime = 0)
  cond2_align <- data.frame(node = "node", pseudotime = 0)
  cond2_noalign <- data.frame(node = "node", pseudotime = 0)
  
  choices <- c(as.numeric(str_replace(alignment_initial[1,1], "X", "")) , as.numeric(str_replace(alignment_initial[1,2], "Y", "")) ) 
  max_choice <- which(choices == max(choices))
  
  #print(new_pseudo1[1,1])
  #print(new_pseudo2[1,1])
  
  #print(max_choice)
  
  #Initial aligning point is at start of condition 1
  if (max_choice == 1){
    new_root_pseudotime <- condition_1[alignment_initial[1,1],]
    diff <- abs(new_root_pseudotime - condition_2[alignment_initial[1,2],])
    
    #as pseudotime of the first node is 0 we can just add the new_root_pseudotime to get our shunted pseudotime
    new_pseudo2 <- new_pseudo2 + diff
    
    #print(new_pseudo1[1,1])
    #print(new_pseudo2[1,1])
  }
  
  else{
    new_root_pseudotime <- condition_2[alignment_initial[1,2],]
    diff <- abs(new_root_pseudotime - condition_1[alignment_initial[1,1],])
    
    #as pseudotime of the first node is 0 we can just add the new_root_pseudotime to get our shunted pseudotime
    new_pseudo1 <- new_pseudo1 + diff
    
    #print(new_pseudo1[1,1])
    #print(new_pseudo2[1,1])
    
  }
  
  col1 <- unique(alignment[duplicated(alignment[,1]),1])
  col2 <- unique(alignment[duplicated(alignment[,2]),2])
  
  squish_list <- list()
  for (i in col1){
    squish_list[[i]] <- alignment[which(alignment[,1] == i),2]
  }
  for (i in col2){
    squish_list[[i]] <- alignment[which(alignment[,2] == i),1]
  }
  
  
  newPseudo_list <- list(new_pseudo1, new_pseudo2)
  #print(length(newPseudo_list[[1]]))
  
  previous_match <- c(alignment[1,1], alignment[1,2])
  
  cond_counter <- c(0,0)
  
  condense <- F
  match <- 1
  
  while(match < length(alignment[,1])){
    current_match <- c(alignment[match,1], alignment[match,2])
    
    #print("")
    #print(previous_match)
    #print(current_match)
    #print(paste0("match: ", match))
    squish_index <- which(current_match %in% names(squish_list))
    #print(paste0("squish index: ", squish_index))
    
    
    #If there is multiple matching
    if(sum(squish_index) > 0){
      
      if(current_match[1] == alignment[length(alignment[,1]),1] |  current_match[2] == alignment[length(alignment[,1]),2] ){
        #return(newPseudo_list)
        #print("end")
        #if it's already at the end of the nodes, then create a new node (pseudotime based on the pseudotime interval of the nodes from the old pseudotime)
        #and squish values between the current node and this new node
        if (current_match[squish_index] == row.names(newPseudo_list[[squish_index]])[length(newPseudo_list[[squish_index]][,1])]){
          #print("if")
          pseudo_interval <- abs(new_pseudo1[1,1] - new_pseudo1[2,1])
          
          squish1 <- newPseudo_list[[squish_index]][current_match[squish_index],1]
          squish2 <- squish1 + pseudo_interval
          
          values_id <- squish_list[[current_match[squish_index]]]
          values_index <- which(!(current_match %in% names(squish_list)))
          newPseudo_list[[values_index]][values_id,1] <- squish_function(newPseudo_list[[values_index]][values_id,1], squish2, squish1)
          
          
          
        }
        
        #squish cells between the node that comes after the one at the end of the alignment 
        else{
          #print("else")
          #print(current_match[squish_index])
          index <- which(row.names(newPseudo_list[[squish_index]]) == current_match[squish_index])
          squish2_id <- row.names(newPseudo_list[[squish_index]])[(index+1)]
          #print(squish2_id)
          
          squish1 <- newPseudo_list[[squish_index]][current_match[squish_index],1]
          squish2 <- newPseudo_list[[squish_index]][squish2_id,1]
          
          values_id <- squish_list[[ current_match[squish_index] ]]
          values_index <- which(!(current_match %in% names(squish_list)))
          newPseudo_list[[values_index]][values_id,1] <- squish_function(newPseudo_list[[values_index]][values_id,1], squish2, squish1)
          
          
        }
        match <- length(alignment[,1])
        
      }
      
      else{
        current_id <- which(row.names(newPseudo_list[[squish_index]]) == current_match[squish_index])
        #Could be that the next node along from the current one has multiple matches, thus we pick the first position
        match <- which(alignment[,squish_index] == row.names(newPseudo_list[[squish_index]])[current_id+1])[1]
        print(match)
        choice <- c(newPseudo_list[[1]][alignment[match,1],] , newPseudo_list[[2]][alignment[match,2],])
        diff <- abs(choice[1] - choice[2])
        #print(paste0("difference:", diff))
        index <- which(choice == max(choice))[1]
        
        newPseudo_list[[index]][which(row.names(newPseudo_list[[index]]) == alignment[match, index]):length(newPseudo_list[[index]][,1]),1] <- newPseudo_list[[index]][which(row.names(newPseudo_list[[index]]) == alignment[match, index]):length(newPseudo_list[[index]][,1]),1] - diff
        
        
        squish1 <- newPseudo_list[[squish_index]][current_match[squish_index],1]
        squish2 <- min(choice)
        values_id <- squish_list[[current_match[squish_index]]]
        #print(values_id)
        values_index <- which(!(current_match %in% names(squish_list)))
        #print(paste0("squish output: ", squish_function(newPseudo_list[[values_index]][values_id,1], squish2, squish1)))
        newPseudo_list[[values_index]][values_id,1] <- squish_function(newPseudo_list[[values_index]][values_id,1], squish2, squish1)
        
      }
      
      
    }
      
    else{
      choice <- c(newPseudo_list[[1]][alignment[match,1],] , newPseudo_list[[2]][alignment[match,2],])
      diff <- abs(choice[1] - choice[2])
      #print(paste0("difference:", diff))
      index <- which(choice == max(choice))[1]
        
      newPseudo_list[[index]][which(row.names(newPseudo_list[[index]]) == alignment[match, index]):length(newPseudo_list[[index]][,1]),1] <- newPseudo_list[[index]][which(row.names(newPseudo_list[[index]]) == alignment[match, index]):length(newPseudo_list[[index]][,1]),1] - diff
        
      #print(paste0("test1: ", newPseudo_list[[1]][which(row.names(newPseudo_list[[1]]) == alignment[match, 1]),1]))
      #print(paste0("test2: ", newPseudo_list[[2]][which(row.names(newPseudo_list[[2]]) == alignment[match, 2]),1]))
      print(index)
      match <- match + 1
        
    }
    previous_match <- current_match
    cond_counter[1] <- which(row.names(newPseudo_list[[1]]) == alignment[match,1])
    cond_counter[2] <- which(row.names(newPseudo_list[[2]]) == alignment[match,2])
    
  }
  
  output <- node_status(newPseudo_list, unmatched)
  
  return(output)
}

#Identify unaligned nodes
node_status <- function(newPseudo_list, unmatched){
  recursion_list <- newPseudo_list
  output <- list()
  for (i in 1:2){
    
    obj <- recursion_list[[i]]
    vector <- c()
    
    for (node in row.names(obj)){
      
      if (node %in% unmatched){
        #vector <- append(vector, paste0("condition_",i,"_noalign"))
        vector <- append(vector, "noalign")
      }
      
      else{
        #vector <- append(vector, paste0("condition_",i,"_align"))
        vector <- append(vector, "align")
      }
      
    }
    #levels_vector <- c( "condition_1_noalign", "condition_1_align", "condition_2_align", "condition_2_noalign" )
    #vector <- factor(vector, levels = levels_vector)
    output[[paste0("condition_", i)]] <- cbind(recursion_list[[i]], vector)
    
  }
  
  return(output)
}

new_squish_node <- function(cond1_pseudo, cond2_pseudo, tree){
  
  pseudo_list <- list(cond1_pseudo, cond2_pseudo)
  
  #tree <- tree[which(tree[,4] != "cut"),]
  
  #Initial alignment
  start_pseudo_choice <- c(pseudo_list[[1]][tree[1,1],1], pseudo_list[[2]][tree[1,2],1])
  
  #The trajectory which has the lowest initial pseudotime needs to have their pseudotimes moved up to come in line with the point of 
  #initial alignment with the other trajectory 
  
  choice_index <- which(start_pseudo_choice == min(start_pseudo_choice))
  
  diff <- abs(start_pseudo_choice[1] - start_pseudo_choice[2])
  
  pseudo_list[[choice_index]][,choice_index] <- pseudo_list[[choice_index]][,choice_index] + diff
  
  unique_nodes <- c(names(table(output_solution_cut$Y))[which(table(output_solution_cut$Y)==1)], names(table(output_solution_cut$X))[which(table(output_solution_cut$X)==1)])
  
  #tree <- rbind(tree, data.frame(X="end", Y="end", Score = 0, Status = "end"))
  item <- 1
  #for (item in 1:length(tree[,1])){
  while( item != (length(tree[,1])+1) ){
   
    pseudo_choice <- c(pseudo_list[[1]][tree[item,1],1], pseudo_list[[2]][tree[item,2],1])
    
    if(pseudo_choice[1] != pseudo_choice[2]){
      
      choice_index <- which(pseudo_choice != min(pseudo_choice))
      
      diff <- abs(pseudo_choice[1] - pseudo_choice[2])
      
      node_position <- which(row.names(pseudo_list[[choice_index]]) == tree[item,choice_index])[1]
      
      pseudo_list[[choice_index]][node_position:(length(pseudo_list[[choice_index]][,1])),1] <- pseudo_list[[choice_index]][node_position:(length(pseudo_list[[choice_index]][,1])),1] - diff
      
    }
  
    #if the item graph is a subset of another larger graph
    
    if (!(tree[item,1] %in% unique_nodes) | !(tree[item,2] %in% unique_nodes)){
     
      #get index of the common node for the multi-matching
      common_index <- which(!(c(tree[item,1], tree[item,2]) %in% unique_nodes))
      
      lower_squish <- pseudo_list[[1]][tree[item,1],1]
      
      
      
      #Find the upper squish value
      pseudo_position <- which(row.names(pseudo_list[[common_index]]) == tree[item,common_index]) 
      
     
      tree_nodes <- tree[ which(tree[,common_index] == row.names(pseudo_list[[common_index]])[pseudo_position+1])[1], 1:2]
     
      #upper_squish <- min( c( pseudo_list[[1]][ tree_nodes[[1]] , 1 ], pseudo_list[[2]][ tree_nodes[[2]] , 1 ] ) )
      squish_nodes <- tree[ which(tree[,common_index]==tree[item,common_index]) ,(3-common_index)]
      
      #Finds the nodes that occur after the final matching in this particular multimatch
      after_nodes <-c( row.names(pseudo_list[[3-common_index]])[ ( which( row.names(pseudo_list[[3-common_index]]) == squish_nodes[length(squish_nodes)] ) +1 ) ],
                       row.names(pseudo_list[[common_index]])[pseudo_position+1])
     
      upper_squish <- min( c(  pseudo_list[[3-common_index]][after_nodes[1],], pseudo_list[[common_index]][after_nodes[2],] ), na.rm = T )
      
    
      
      output <- squish_function(values = pseudo_list[[3-common_index]][squish_nodes,1], upper_squish, lower_squish)
      
     
      
      pseudo_list[[3-common_index]][squish_nodes,1] <- output
      
      item <- item + length(squish_nodes) - 1
      
    }
    item <- item + 1
 
  }
  
  return(pseudo_list)
  
}

chunk_generator <- function(tree){
  
  output <- list()
  
  previous <- tree[1,4]
  
  vector_c <- c(1)
  
  list_counter <- 1
  
  for (i in 2:length(tree[,1])){
    
    current <- tree[i,4]
    
    if (current != previous){
      
      output[[list_counter]] <- vector_c
      
      vector_c <- c()
      
      list_counter <- list_counter + 1
      
    }
    
    vector_c <- append(vector_c, i)
    
    previous <- current
    
  }
  
  output[[list_counter]] <- vector_c
  
  return(output)
  
}

chunk_node <- function(cond1_pseudo, cond2_pseudo, tree){
  
  #get list of nodes which are only present on cut parts
  
  cut_nodes <- c(unique(tree[tree[,4] == "cut",1]), unique(tree[tree[,4] == "cut",2]))
  match_nodes <- c(unique(tree[tree[,4] == "match",1]), unique(tree[tree[,4] == "match",2]))
  
  cut_nodes <- cut_nodes[ which( !(cut_nodes %in% match_nodes)) ]
  
  pseudo_list <- list(cond1_pseudo, cond2_pseudo)
  
  cut_tree <- tree[which(tree[,4] != "cut"),]
  
  #Initial alignment
  start_pseudo_choice <- c(pseudo_list[[1]][cut_tree[1,1],1], pseudo_list[[2]][cut_tree[1,2],1])
  
  #The trajectory which has the lowest initial pseudotime needs to have their pseudotimes moved up to come in line with the point of 
  #initial alignment with the other trajectory 
  
  choice_index <- which(start_pseudo_choice == min(start_pseudo_choice))[1]
  
  diff <- abs(start_pseudo_choice[1] - start_pseudo_choice[2])
  
  print(diff)
  
  pseudo_list[[choice_index]][,1] <- pseudo_list[[choice_index]][,1] + diff
  
  #remove the points before the initial alignment 
  cut_index <- which(tree[,4] == "match")
  tree <- tree[cut_index[1]:length(tree[,1]),]
  
  chunks <- chunk_generator(tree)

  seperate_ <- function(input){min(input)}
  
  for (i in 1:length(chunks)){
    tree_chunk <- tree[chunks[[i]],]
    
    item <- 1
    
    unique_nodes <- c(names(table(tree_chunk[,1]))[which(table(tree_chunk[,1])==1)], names(table(tree_chunk[,2]))[which(table(tree_chunk[,2])==1)])
    
    MC_choice <- tree_chunk[1,4] == "match"

    if (MC_choice){
      
       while( item != (length(tree_chunk[,1])+1) ){
         
         pseudo_choice <- c( pseudo_list[[1]][tree_chunk[item,1],1], pseudo_list[[2]][tree_chunk[item,2],1] )
         
         if(pseudo_choice[1] != pseudo_choice[2]){
           
           choice_index <- which(pseudo_choice != seperate_(pseudo_choice))
           
           diff <- pseudo_choice[3-choice_index] - pseudo_choice[(choice_index)]
           
           node_position <- which(row.names(pseudo_list[[choice_index]]) == tree_chunk[item,choice_index])[1]
           
           pseudo_list[[choice_index]][node_position:(length(pseudo_list[[choice_index]][,1])),1] <- pseudo_list[[choice_index]][node_position:(length(pseudo_list[[choice_index]][,1])),1] + diff
           
           seperate_ <- function(input){min(input)}
           
         }
         
         #if the item graph is a subset of another larger graph
         if (!(tree_chunk[item,1] %in% unique_nodes) | !(tree_chunk[item,2] %in% unique_nodes)){
           
           #get index of the common node for the multi-matching
           common_index <- which(!(c(tree_chunk[item,1], tree_chunk[item,2]) %in% unique_nodes))
           
           print(tree_chunk[item,common_index])
           
           lower_squish <- pseudo_list[[common_index]][tree_chunk[item,common_index],1]
           
           #Find the upper squish value
           pseudo_position <- which(row.names(pseudo_list[[common_index]]) == tree_chunk[item,common_index]) 
           
           tree_chunk_nodes <- tree_chunk[ which(tree_chunk[,common_index] == row.names(pseudo_list[[common_index]])[pseudo_position+1])[1], 1:2]
           
           #upper_squish <- min( c( pseudo_list[[1]][ tree_chunk_nodes[[1]] , 1 ], pseudo_list[[2]][ tree_chunk_nodes[[2]] , 1 ] ) )
           squish_nodes <- tree_chunk[ which(tree_chunk[,common_index]==tree_chunk[item,common_index]) ,(3-common_index)]
           print(paste0("squish: ", squish_nodes))
           
           #Finds the nodes that occur after the final matching in this particular multimatch
           after_nodes <-c( row.names(pseudo_list[[3-common_index]])[ ( which( row.names(pseudo_list[[3-common_index]]) == squish_nodes[length(squish_nodes)] ) +1 ) ],
                            row.names(pseudo_list[[common_index]])[pseudo_position+1])
           
           #If we're doing a multi-align at the end of the process, there will be no nodes that occur after it to squish between
           if (is.na(after_nodes) == T){
             print("argh")
             upper_squish <- lower_squish + mean( abs( pseudo_list[[common_index]][,1] - mean(pseudo_list[[common_index]][,1]) ) )

           }
           
           else{
             upper_squish <- min( c(  pseudo_list[[3-common_index]][after_nodes[1],], pseudo_list[[common_index]][after_nodes[2],] ), na.rm = T )
             
           }
           
           upper_squish <- upper_squish - ((upper_squish-lower_squish)/ length(squish_nodes))
           
           print(paste0("upper squish", upper_squish))
           print(paste0("lower squish", lower_squish))
           
           
           
           output <- squish_function(values = pseudo_list[[3-common_index]][squish_nodes,1], upper_squish, lower_squish)
           output <- scale_function(pseudo_list[[3-common_index]][squish_nodes,1], upper_squish, lower_squish)
           print(paste0("output: ", output))
           diff <- abs(pseudo_list[[3-common_index]][ squish_nodes[length(squish_nodes)] , 1 ] - output[length(output)])
           
           pseudo_list[[3-common_index]][squish_nodes,1] <- output
           
           #Don't want to pulldown if we've already reached the end of the process we're moving
           if (!( !( tree[length(tree[,1]),1] %in% squish_nodes ) | !( tree[length(tree[,1]),2]  %in% squish_nodes ) )){
             pulldown_index <- which(row.names(pseudo_list[[3-common_index]]) == squish_nodes[length(squish_nodes)]) + 1
             
             pseudo_list[[3-common_index]][ pulldown_index : length(pseudo_list[[3-common_index]][,1]) , 1] <- pseudo_list[[3-common_index]][ pulldown_index : length(pseudo_list[[3-common_index]][,1]) , 1]  - diff
             
           }
           
           item <- item + length(squish_nodes) - 1
           
         }
         item <- item + 1
         
       }
       
    }
    
    #If the current chunk is unaligned then we want the next aligned bit to be separated in terms of pseudotime
    else{
      seperate_ <- function(input){max(input)}
    }
  
  }
  
  matrix_1 <- pseudo_list[[1]]
  matrix_1 <- data.frame("pseudotime" = matrix_1, alignment = rep("align", length(matrix_1[,1])))
  
  matrix_2 <- pseudo_list[[2]]
  matrix_2 <- data.frame("pseudotime" = matrix_2, alignment = rep("align", length(matrix_2[,1])))
  
  matrix_1[which(!(row.names(matrix_1) %in% match_nodes)),2] <- "noalign"
  matrix_2[which(!(row.names(matrix_2) %in% match_nodes)),2] <- "noalign"
  
  return(pseudo_list <- list(condition_1 = matrix_1, condition_2 = matrix_2))
  
}

#find the best starting point and ending point for the trajectories
#Starts from highest correlated point and works backwards 
cut_point <- function(cor_mtx){
  
  #find start point
  
  choice <- c(min(cor_mtx[,1]), min(cor_mtx[1,]))
  
  index_choice <- which(choice == min(choice))
  
  if (index_choice == 1){
    loc <- c(which(cor_mtx[,1] == choice[index_choice]), 1)
  }
  
  else{
    loc <- c(1, which(cor_mtx[1,] == choice[index_choice]))
  }
  
  cor_mtx <- cor_mtx[ loc[1]:length(cor_mtx[,1]) , loc[2]:length(cor_mtx[1,]) ]
  
  
  
  #Find end point
  
  last_row <- length(cor_mtx[,1])
  last_column <- length(cor_mtx[1,])
  
  choice <- c( min( cor_mtx[ last_row , ] ), min( cor_mtx[ , last_column ] ) )
  
  index_choice <- which(choice == min(choice))
  
  if (index_choice == 1){
    loc <- c( last_row , which( cor_mtx[ last_row , ] == choice[index_choice] ) )
  }
  
  else{
    loc <- c( which( cor_mtx[ , last_column ] == choice[index_choice] ), last_column )
  }
  
  
  cor_mtx <- cor_mtx[ 1:loc[1] , 1:loc[2] ]
  
  
}

#scales values
scale_function <- function(vals, scale_max, scale_min){
  print(paste0("values scale: ",vals))
  print(paste0("upper squish", scale_max))
  print(paste0("lower squish", scale_min))
  output <- ( (vals - min(vals)) / (max(vals) - min(vals)) ) *(  (scale_max - scale_min)) + scale_min
    
  return(output)
  
}

#Output
x <- function(){
  previous_match <- c(alignment[1,1], alignment[1,2])
  
  cond_counter <- c(0,0)
  
  condense <- F
  for (match in 1:length(alignment[,1])){
    current_match <- c(alignment[match,1], alignment[match,2])
    
    
    print("")
    #If there is multiple matching
    pass_condition <- which(current_match %in% previous_match)
    print(previous_match)
    print(current_match)
    print(paste0("match: ", match))
    print(paste0("problemstart: ", newPseudo_list[[1]][which(row.names(newPseudo_list[[1]]) == alignment[3, 1]),1]))
    
    if (sum(pass_condition) > 0){
      print("hi")
      condense <- T
    }
    
    else{
      
      
      choice <- c(newPseudo_list[[1]][alignment[match,1],] , newPseudo_list[[2]][alignment[match,2],])
      diff <- abs(choice[1] - choice[2])
      print(paste0("difference:", diff))
      index <- which(choice == max(choice))[1]
      
      newPseudo_list[[index]][which(row.names(newPseudo_list[[index]]) == alignment[match, index]):length(newPseudo_list[[index]][,1]),1] <- newPseudo_list[[index]][which(row.names(newPseudo_list[[index]]) == alignment[match, index]):length(newPseudo_list[[index]][,1]),1] - diff
      
      print(paste0("test1: ", newPseudo_list[[1]][which(row.names(newPseudo_list[[1]]) == alignment[match, 1]),1]))
      print(paste0("test2: ", newPseudo_list[[2]][which(row.names(newPseudo_list[[2]]) == alignment[match, 2]),1]))
      print(index)
      
      if (condense == T){
        print(paste0("cond_counter: ", cond_counter))
        squish_index <- which(cond_counter == max(cond_counter))
        print(paste0("index: ", index))
        
        squish1 <- newPseudo_list[[squish_index]][(min(cond_counter)),1]
        squish2 <- newPseudo_list[[squish_index]][(min(cond_counter)+1),1]
        
        print(paste0("squish 1: ", squish1))
        print(paste0("squish 2: ", squish2))
        
        print(paste0("match - 1: ", (match-1)))
        print(paste0("match - max(con d_counter): ", (match - max(cond_counter))))
        newPseudo_list[[which(cond_counter == min(cond_counter))]][min(cond_counter):max(cond_counter),1] <- squish_function(newPseudo_list[[which(cond_counter == min(cond_counter))]][min(cond_counter):max(cond_counter),1], squish2, squish1)
      }
      print(paste0("problem: ", newPseudo_list[[1]][which(row.names(newPseudo_list[[1]]) == alignment[3, 1]),1]))
      
      condense <- F
    }
    
    previous_match <- current_match
    cond_counter[1] <- which(row.names(newPseudo_list[[1]]) == alignment[match,1])
    cond_counter[2] <- which(row.names(newPseudo_list[[2]]) == alignment[match,2])
  }
}

squish_function <- function(values, max_lim, min_lim){
  diff <- abs(min_lim - max_lim)
  val <- diff / (length(values))
  output <- c(min_lim)
  value_progress <- min_lim
  for (i in 2:length(values)){
    value_progress <- value_progress + val
    output <- append(output, value_progress)
    
  }
    
    
  return(output)
  
}

#Cell new pseudoptime calculated based on the singular closest node
pseudo_align_cells <- function(cell_pseudo, node_pseudo_new, node_pseudo_old){
  #As we're dealing with slingshot we need to remove cells which have a pseudotime of NA
  cell_pseudo[,2] <- cell_pseudo[ which(is.na(cell_pseudo[,2]) == F),2]
  
  #Do cell capture based on pseudotime between node_pseudo_old and node_pseudo_new
  kmean_output <- cell_capture(cell_pseudo, node_pseudo_old)
  
  for ( node in names(kmean_output)){
    print(node_pseudo_new[node,2])
    print(node_pseudo_old[node,2])
    node_change <- node_pseudo_new[node,2] - node_pseudo_old[node,2]
    print(node_change)
    
    cell_pseudo[kmean_output[[node]],2] <- cell_pseudo[kmean_output[[node]],2] + node_change
    
  }
  
  return(cell_pseudo)
}

#Was doing the correlation of the cells closest to the unaligned nodes which were next to aligned ones so make sure we weren't chucking out anything useful
#but, basically, if it's closer in pseudotime to an unaligned node, it will be more correlated to that node than to a correlated one
pseudo_cell_align_deprecated <- function(cell_pseudo, node_pseudo_new, node_pseudo_old, window_size, node_expr_mtx, cell_expr_mtx){
  
  #Make sure genes in cell_expr_mtx are the genes in node_expr_mtx
  cell_expr_mtx <- cell_expr_mtx[which(row.names(cell_expr_mtx) %in% row.names(node_expr_mtx)),]
  
  cell_expr_mtx <- cell_expr_mtx[order(match(row.names(cell_expr_mtx), row.names(node_expr_mtx))),]
  
  cell_pseudo_new <- cellAlign_function(node_pseudo_old, node_pseudo_new, cell_pseudo, window_size)
  
  #get intervals of non-alignment
  
  kmean_output <- cell_capture(cell_pseudo, node_pseudo_old)
  
  
  noalign_pseudo <- data.frame(ID = "ID", pseudotime = 0)
  align_pseudo <- data.frame(ID = "ID", pseudotime = 0)
  previous_node <- node_pseudo_new[1,1]
  previous_item <- node_pseudo_new[1,3]
  
  for (i in 1:length(node_pseudo_new[,3])){
    item <- node_pseudo_new[i,3]
    node <- node_pseudo_new[i,1]
    
    #If we need to decide which cells should not be on the aligned pseudotime
    if(item != previous_item){
      choice <- which(c(item, previous_item) == "noalign")
      
      if(choice == 1){
        cells_decide <- kmean_output[[node]]
        
        node_expr_vector <- as.matrix(node_expr_mtx[,node])
        
        print(node)
        print(previous_node)
        return(list(node_expr_mtx[,c(node, previous_node)], cell_expr_mtx[,cells_decide]))
        
        cor_mtx <- dis_mtx_calculator(data.frame(node = node_expr_mtx[,node]), cell_expr_mtx[,cells_decide], "spearman", "median", "linear")
        

      }
      
      else{
        cells_decide <- kmean_output[[previous_node]]
        cor_mtx <- dis_mtx_calculator(node_expr_mtx[,c(node, previous_node)], cell_expr_mtx[,cells_decide], "spearman", "median", "linear")
        
      }
      
    }
    
    else{
      
      if(item == "align"){
        df <- data.frame(ID = kmean_output[[item]], pseudotime = cell_pseudo_new[ kmean_output[[item]] , 1])
        align_pseudo <- rbind(align_pseudo, df)
        
      }
      
      #We know with a degree of certainty that the cells are either aligned or unaligned because the adjacent node isn't different to the current alignment status
      else{
        df <- data.frame(ID = kmean_output[[item]], pseudotime = cell_pseudo_new[ kmean_output[[item]] , 1])
        noalign_pseudo <- rbind(noalign_pseudo, df)
      }
      
    }
    
    
    previous_node <- node
    previous_item <- item
    
  }
  
  
  return(cell_pseudo_new)
}

#Cell pseudotime calculated using the CellAlign function
#Slightly increase correlation of the new cell pseudotime with the old cell pseudotime (preserves the original order better)
pseudo_cell_align <- function(cell_pseudo, node_pseudo_new, node_pseudo_old, window_size){
  
  cell_pseudo_new <- cellAlign_function(node_pseudo_old, node_pseudo_new, cell_pseudo, window_size)
  
  
  kmean_output <- cell_capture(cell_pseudo, node_pseudo_old)
  #Some nodes have no nodes close to it and have to be removed
  index <- which(lengths(kmean_output) != 0)
  kmean_output <- kmean_output[names(index)]
  node_pseudo_new <- node_pseudo_new[which(node_pseudo_new[,1] %in% names(index)),]
  
  #return(cell_pseudo_new)
  
  #If the closest node to a cell is unaligned, that cell becomes unaligned
  noalign_pseudo <- data.frame(ID = "ID", pseudotime = 0, status = "hi")
  align_pseudo <- data.frame(ID = "ID", pseudotime = 0, status = "hi")
  for (i in 1:length(node_pseudo_new[,3])){
    item <- node_pseudo_new[i,3]
    node <- node_pseudo_new[i,1]
    print(node)
    if (item == "noalign"){
      df <- data.frame(ID = kmean_output[[node]], pseudotime = cell_pseudo_new[ kmean_output[[node]] , 1], status = "noalign")
      noalign_pseudo <- rbind(noalign_pseudo, df)
      
    }
    
    else{
      #return(kmean_output)
      df <- data.frame(ID = kmean_output[[node]], pseudotime = cell_pseudo_new[ kmean_output[[node]] , 1], status = "align")
      align_pseudo <- rbind(align_pseudo, df)
    }
    
  }
  
  noalign_pseudo <- noalign_pseudo[-1,]
  align_pseudo <- align_pseudo[-1,]
  
  output <- rbind(noalign_pseudo, align_pseudo)
  
  #Keep order of cells the same as it is in the sce and seurat objects
  output <- output[order(match(output[,1], cell_pseudo[,1])),]
  output <- output[,-1]
  
  return(output)
}



#DE methods
#Get the rate of change in gene expression
gROC <- function(exp_mtx_1, exp_mtx_2, pseudo_1, pseudo_2, window, min.cell){
  
  #remove cells that aren't aligned
  pseudo_1 <- pseudo_1[which(pseudo_1$status == "align"), ]
  pseudo_2 <- pseudo_2[which(pseudo_2$status == "align"), ]
  
  exp_mtx_1 <- exp_mtx_1[, pseudo_1$ID]
  exp_mtx_2 <- exp_mtx_2[, pseudo_2$ID]
  
  output <- list("condition_1" = list(), "condition_2" = list())
  
  #get window size that fits the aligned pseudotime axis
  lower <- min( c(min(pseudo_1[pseudo_1$status == "align", 2]), min(pseudo_2[pseudo_2$status == "align", 2] )) )
  upper <- max( c(max(pseudo_1[pseudo_1$status == "align", 2]), max(pseudo_2[pseudo_2$status == "align", 2]) ) )
  
  diff <- abs(lower-upper)
  window <-  diff / floor(diff/window)

  common_genes <- intersect(row.names(exp_mtx_1), row.names(exp_mtx_2))
  print(length(common_genes))
  
  exp_mtx_1 <- exp_mtx_1[common_genes,]
  exp_mtx_2 <- exp_mtx_2[common_genes,]
  
  exp_row_sum_1 <- rowSums(exp_mtx_1 != 0)
  exp_row_sum_2 <- rowSums(exp_mtx_2 != 0)

  genes_keep_1 <- names(exp_row_sum_1)[exp_row_sum_1 > min.cell]
  genes_keep_2 <- names(exp_row_sum_2)[exp_row_sum_2 > min.cell]
  
  genes_keep <- intersect(genes_keep_1, genes_keep_2)
  
  exp_mtx_1 <- exp_mtx_1[genes_keep,]
  exp_mtx_2 <- exp_mtx_2[genes_keep,]
  
  exp_list <- list(exp_mtx_1, exp_mtx_2)
  pseudo_list <- list(pseudo_1, pseudo_2)
  print(length(genes_keep))
  
  #for each condition
  
  test_list <- list( "condition_1" = list(), "condition_2" = list() )
  
  for (i in 1:length(exp_list)){
    
    #exp_list[[i]] <- exp_list[[i]][,aligned_cells]
    #pseudo_list[[i]] <- pseudo_list[[i]][aligned_cells,]
    
    current_pseudo <- lower
    
    for (j in 1:( diff/window ) ){
      
      lower_capture <- current_pseudo
      
      upper_capture <- lower_capture + window
      
      cell_window_1 <- pseudo_list[[i]][pseudo_list[[i]][,2] < upper_capture & pseudo_list[[i]][,2] > lower_capture, 1 ]
      
      window_exp_list_1 <- exp_list[[i]][,cell_window_1]
      
      #find the cells in the next part of the sliding window
      
      lower_capture <- upper_capture - (window/2)
      upper_capture <- lower_capture + window
      
      cell_window_2 <- pseudo_list[[i]][pseudo_list[[i]][,2] < upper_capture & pseudo_list[[i]][,2] > lower_capture, 1 ]
      
      print(paste0("cells one: ", length(cell_window_1)))
      
      print(paste0("cells two: ", length(cell_window_2)))
      
      window_exp_list_2 <- exp_list[[i]][,cell_window_2]
      
      if ( length(window_exp_list_1[1,]) > 10 & length(window_exp_list_2[1,]) > 10 ){
        #return(list(window_exp_list_1, window_exp_list_2))
        window_1_pseudobulk <- pseudobulk(window_exp_list_1, ( length(window_exp_list_1[1,]) /10 ) )
        
        window_2_pseudobulk <- pseudobulk(window_exp_list_2, ( length(window_exp_list_2[1,]) /10 ) )
      }
      
      else{
        window_1_pseudobulk <- window_exp_list_1
        
        window_2_pseudobulk <- window_exp_list_2
      }
      
      #window_1_pseudobulk <- window_exp_list_1
      
      #window_2_pseudobulk <- window_exp_list_2
      
      #store_mtx structure <- rows = genes, columns = cells, values = average change 
      
      store_mtx <- matrix(nrow =1 , ncol = length(window_2_pseudobulk[1,]))
      
      for (gene in row.names(window_1_pseudobulk)){
        window_mtx_1 <- matrix( rep( window_1_pseudobulk[gene,], each= length(window_2_pseudobulk[1,]) ),nrow= length(window_2_pseudobulk[1,]) )

        window_mtx_2 <- matrix( rep( window_2_pseudobulk[gene,], each= length(window_1_pseudobulk[1,]) ), ncol= length(window_1_pseudobulk[1,]), byrow=TRUE )
          
        difference <- window_mtx_2 - window_mtx_1
        
        #difference <- colSums(difference) / length(difference[,1])
        difference <- rowSums(difference) / length(difference[1,])
        
        #return(list(window_mtx_1, window_mtx_2, difference))
        
        store_mtx <- rbind(store_mtx, difference)
      }
      #print(dim(store_mtx))
      store_mtx <- store_mtx[-1,]
      
     
      
      row.names(store_mtx) <- row.names(window_exp_list_1)
      #print(dim(store_mtx))
      
     
      
      #return(store_mtx)
      
      output[[paste0("condition_", i)]][[paste0("slide_", j)]] <- store_mtx
      
      if(j == 1){
        #return(list(window_1_pseudobulk, window_2_pseudobulk, output))
      }
      

      current_pseudo <- lower_capture
        
    }
    
    
  }
  #print("problem")
  
  #collect cells that are in pseudo_window_i and pseudo_window_i+1
  
  #remove genes whose sum is zero across all the cells in both windows
  
  #for every gene
  
  #distance matrix of cells in window_i and window_i+1
  
  #average for each cell in window for the cells in the other window
  
  #compare cells change average between the conditions with t-test
  
  #return(output)
  
  
  output <- test_function(output, "t.test")
  
  return(output)
  
  
}

#random pseudobulk of cells
pseudobulk <- function(exp_mtx, n.cells){
  print(n.cells)
  n.cells <- round(n.cells)
  
  #print(n.cells)
  
  partitions <- round(runif(length(exp_mtx[1,]),1, n.cells))
  
  output <- matrix(0, nrow = length(exp_mtx[,1]), ncol = 1)
  
  #print(table(partitions))
  #print(paste0("max: ", max(partitions)))
  
  for (i in 1:n.cells){
    #print(paste0("partition i: ", i))
    bulk_mtx <- exp_mtx[, which(partitions == i)]
    
    #print(dim(bulk_mtx))
    
    if(length(which(partitions==i)) > 1){
      bulk_mtx <- rowSums(bulk_mtx) / (length(which(partitions == i)))
      
    }
    
    output <- cbind(output, bulk_mtx)
    
  }
  output <- output[,-1]
  
  print(dim(output))

  row.names(output) <- row.names(exp_mtx)
  colnames(output) <- paste0("c", seq(1, length(output[1,]) , 1) )
  
  
  return(output)
  
}

test_function <- function(list_input, test.use){
  
  #gene name, p_val, adj.p_val, change
  output <- data.frame("gene_name" = "g", "p_val" = 0, "cond1_change" = 0, "cond2_change" = 0,"slide" = 0, "bon_pval" = 0)
  
  
  stats_output <- data.frame("gene_name" = "g", "p_val" = 0, "cond1_change" = 0, "cond2_change" = 0,"slide" = 0)
  
  for ( i in 1:length(list_input[[1]]) ){
    
    for (gene in row.names(list_input[[1]][[1]])){
      
      if (test.use == "t.test"){
        #print(gene)
        
        print(list_input[[1]][[i]][ gene, ])
        
        print(list_input[[2]][[i]][ gene, ])
        
        #if ( length(unique(list_input[[1]][[i]][ gene, ])) != 1 & length(unique(list_input[[2]][[i]][ gene, ])) != 1 ){
          test_out <- t.test(list_input[[1]][[i]][ gene, ], list_input[[2]][[i]][ gene, ])
          
          p_val <- test_out$p.value
          
        #}
        
        #else{
          #p_val <- -1
          
        #}
        
        
        stats_output <- rbind(stats_output, data.frame("gene_name" = gene, "p_val" = p_val, "cond1_change" =  mean(list_input[[1]][[i]][ gene, ]), "cond2_change" = mean(list_input[[2]][[i]][ gene, ])  , "slide" = i) )
        
        
      }
      
      
    }
    
    stats_output <- stats_output[-1,]
    
    output <- rbind( output, data.frame(stats_output, "bon_pval" = stats::p.adjust(stats_output[,"p_val"], method = "bonferroni", length(stats_output[,"gene_name"])) ) )
    
    stats_output <- data.frame("gene_name" = "g", "p_val" = 0, "cond1_change" = 0, "cond2_change" = 0,"slide" = 0)
    
    
  }
  
  output <- output[-1,]
  
  return(output)
  
}

#Plot's the output of alignment
PlotOutput <- function(pseudo_1, pseudo_2, ID_1, ID_2, alignment, cond_names){
  
  cut_alignment <- alignment[alignment$Status == "match",]
  
  unaligned_nodes <- append(row.names(pseudo_1)[which( !(row.names(pseudo_1) %in% cut_alignment[,1]) )] , row.names(pseudo_2)[which( !(row.names(pseudo_2) %in% cut_alignment[,2]) )] )
  
  plot_df <- data.frame("node" = "x", "group" = 0, "pseudotime" = 0, "condition" = "ah", "ID" = "ah")
  
  plot_counter <- 1
  
  for (i in 1:length(cut_alignment[,1])){
    print(cut_alignment[i,1])
    plot_df <- rbind(plot_df, data.frame("node" = cut_alignment[i,1], "group" = plot_counter, "pseudotime" = pseudo_1[cut_alignment[i,1],], "condition" = cond_names[1] , 
                     "ID" = ID_1[cut_alignment[i,1],]) )
    
    plot_df <- rbind(plot_df, data.frame("node" = cut_alignment[i,2], "group" = plot_counter, "pseudotime" = pseudo_2[cut_alignment[i,2],], "condition" = cond_names[2] , 
                     "ID" = ID_2[cut_alignment[i,2],]) )
    
    plot_counter <- plot_counter + 1
  }
  
  plot_df <- plot_df[-1,]
  #add unaligned nodes
  unaligned_frame <- data.frame("node" = unaligned_nodes,
                                "group" = seq(plot_counter, (plot_counter + length(unaligned_nodes)-1), 1),
                                "pseudotime" = c(pseudo_1[unaligned_nodes,1], pseudo_2[unaligned_nodes,1])[ is.na(c(pseudo_1[unaligned_nodes,1], pseudo_2[unaligned_nodes,1])) == F ],
                                "condition" = c( rep(cond_names[1], length( which( grepl("X",unaligned_nodes) == T) ) ), rep(cond_names[2], length( which( grepl("Y",unaligned_nodes) == T ) ) )),
                                                 "ID" = c( ID_1[unaligned_nodes[which( grepl("X",unaligned_nodes) == T)],], ID_2[unaligned_nodes[which( grepl("Y",unaligned_nodes) == T)],] ) ) 
  
  
  plot_df <- rbind(plot_df, unaligned_frame)
  
  print(max(plot_df$pseudotime))
  
  print( ggplot(data = plot_df, aes(pseudotime, condition)) + geom_point(aes(col = ID)) + geom_line(aes(group = group)) + theme_classic() )#+scale_x_continuous(breaks = seq(0, max(plot_df$pseudotime), max(plot_df$pseudotime)/5)) )
  
}

#plot heatmap of gene logfc over windows
gene_window_heatmap <- function(de_output, conditions = 0, window, top_n, genes = 0){
  
  genes_vector <- c()
  heatmap_cols <-c()
  for (i in 1:length(de_output)){
    genes_vector <- append(genes_vector, de_output[[i]][,"gene_name"])
    heatmap_cols <- append( heatmap_cols, paste("logfc", i, sep="_") )
    
  }
  
  genes_vector <- unique(genes_vector)
  len_genes <- length(genes_vector)
  
  heatmap_mtx <- matrix(0, nrow = len_genes, ncol = length(heatmap_cols), dimnames = list(genes_vector, heatmap_cols))
  
  for (i in genes_vector){
    
    for (j in 1:length(de_output)){
      
      if (i %in% de_output[[j]][,"gene_name"]){
        heatmap_mtx[i,paste("logfc", j, sep="_")] <- de_output[[j]][i,"logfc"]
        
      }
      
      else{
        heatmap_mtx[i,paste("logfc", j, sep="_")] <- NA
        
      }
      
    }
    
  }
  
  heatmap_melt <- melt(heatmap_mtx)
  colnames(heatmap_melt) <- c("Gene", "Window", "Logfc")
  
  if(genes != 0){
    heatmap_melt <- heatmap_melt[which(heatmap_melt$Gene %in% genes),]
    
  }
  
  #return(heatmap_melt)
  
  heatmap_melt_subset <- subset(heatmap_melt, Window %in% paste("logfc",window, sep="_") )
  
  #split logfc + and -
  heatmap_melt_pos <- subset(heatmap_melt_subset, Logfc > 0)
  heatmap_melt_neg <- subset(heatmap_melt_subset, Logfc < 0)
  
  obj_list <- list(heatmap_melt_subset, heatmap_melt_pos, heatmap_melt_neg)
  
  #remove NA values
  #heatmap_melt_subset <- heatmap_melt_subset[which(!(is.na(heatmap_melt_subset[,"Logfc"]))),]
  
  for (i in 1:length(obj_list)){
    #heatmap_melt <- heatmap_melt[heatmap_melt[,"Gene"] %in% heatmap_melt_subset[,"Gene"],]
    #heatmap_melt <- heatmap_melt_subset
    obj_current <- obj_list[[i]]
    
    obj_current <- obj_current %>%
      group_by(Window) %>%
      slice_max(n = top_n, order_by = abs(Logfc))
    
    gene_list <- unique(obj_current$Gene)
    #heatmap_melt_old <- heatmap_melt
    
    obj_current <- heatmap_melt[heatmap_melt$Gene %in%gene_list,]
    
    order_frame <- c()

    for ( j in unique(obj_current$Gene) ){

      order_frame[j] <- mean(obj_current[obj_current$Gene == j, "Logfc"], na.rm =T)
      print(j)
      print(obj_current[obj_current$Gene == j, "Logfc"])
      print("")
    }

    order_frame <- sort(order_frame, decreasing = T)

    obj_current$Gene <- as.factor(as.character(obj_current$Gene) )
    obj_current$Gene <- factor(obj_current$Gene, levels = names(order_frame))
    print(length(obj_current$Gene))
    

    obj_list[[i]] <- obj_current
    
    print("")
    
    
  }
  

  #heatmap_melt <- heatmap_melt[heatmap_melt[,"Gene"] %in% heatmap_melt_subset[,"Gene"],]
  #heatmap_melt <- heatmap_melt_subset
  # heatmap_melt_subset <- heatmap_melt_subset%>%
  #   group_by(Window) %>%
  #   slice_max(n = top_n, order_by = abs(Logfc))
  # 
  # gene_list <- unique(heatmap_melt_subset$Gene)
  # print(gene_list)
  # heatmap_melt_old <- heatmap_melt
  # 
  # heatmap_melt <- heatmap_melt[heatmap_melt$Gene %in%gene_list,]
  
  #print(ggplot(data =heatmap_melt, aes(x= Window, y =Gene, fill = Logfc)) + geom_tile()+scale_fill_gradient2())
  return(obj_list)
  
}



PlotAlignment <- function(alignment, score_mtx){
  cut_alignment <- alignment[alignment$Status == "match",]
  
  score_mtx <- t(score_mtx)
  score_mtx[as.matrix(cut_alignment[,1:2])] <- NA
  
  pheatmap::pheatmap(as.matrix(score_mtx), cluster_rows = F, cluster_cols = F, scale="none")
}


