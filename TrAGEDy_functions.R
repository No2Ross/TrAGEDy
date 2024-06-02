

#Assigns cells to a node based on euclidean distance for some metric
#Functions as one round k-means clustering where k is the number of nodes
#First column is names of cells/nodes
#Second column is value (pseudotime)
cell_capture <- function(cell_embed, node_embed){
  output_list <- list()
  
  node_mtx <- t( replicate(length(cell_embed), node_embed) )
  cell_mtx <- replicate( length(node_embed) , cell_embed)
  
  dist_mtx <- ( node_mtx - cell_mtx )^2
  
  colnames(dist_mtx) <- names(node_embed)
  row.names(dist_mtx) <- names(cell_embed)
  
  row_min <- rowMins(dist_mtx)
  
  closest_node <- as.matrix(apply( dist_mtx, 1, which.min))
  
  for (i in 1:length( colnames(dist_mtx) )){
    output_list[[colnames(dist_mtx)[i]]] <- row.names(closest_node)[which(closest_node == i)]
    
  }

  return(output_list)
  
}


nodeExpressionEstimate <- function(cell_exp_mtx, tree, window_size, adjust.window){
  
  #initialize variables
  cell_pseudo <- tree$cell_pseudotime
  node_pseudo <- tree$node_pseudotime
  
  output <- list()
  
  if(adjust.window){
    #Get window sizes for each interpolated point
    hist_output <- hist(cell_pseudo, breaks = seq(min(cell_pseudo) , max(cell_pseudo), length.out = (length(node_pseudo) + 1 ) ) )
    
    #If it's 0, then it won't calculate it
    #hist_output$counts <- hist_output$counts + 1
    
    hist_mean <- mean(hist_output$breaks)
    
    alt_hist_mean <- hist_output$counts[min(which(hist_output$breaks > mean(cell_pseudo)))]
    
    #If the interpolated point is in an area of pseudotime with a cell density less than the mean, we increase the window size. 
    #If it's an area of high cell density, we decrease the window
    window_vector <- window * (1 - ( ( hist_output$counts- alt_hist_mean ) / max(hist_output$counts)) )
    
  }
  
  else{
    window_vector <- rep( window, length(node_pseudo) )
  }
  
  
  node_exp_mtx <- matrix(nrow = 50, ncol = length(node_pseudo))
  
  weight_vector <- c()
  weight_mtx <- matrix(nrow = length(cell_pseudo), ncol = length(node_pseudo))
  
  for (cell in 1:length(cell_pseudo)){
    cell_weight_vector <- c()
    
    for (node in 1:length(node_pseudo)){
      weight <- exp( -1 * ( (cell_pseudo[cell] - node_pseudo[node])^2 / window_vector[node]^2 ) )
      cell_weight_vector <- append(cell_weight_vector, weight)
      
    }
    weight_mtx[cell,] <- cell_weight_vector
    
  }
  
  row.names(weight_mtx) <- names(cell_pseudo)
  colnames(weight_mtx) <- names(node_pseudo)
  
  
  #mtx multiplication way
  node_exp_mtx <- cell_exp_mtx %*% weight_mtx 
  
  normalizer <- 1/(colSums(weight_mtx))
  
  for (node in 1:length(node_exp_mtx[1,])){
    node_exp_mtx[,node] <- node_exp_mtx[,node] * normalizer[node]
  }
  
  output[["exp_mtx"]] <- node_exp_mtx
  
  
  
  return(node_exp_mtx)
  
}

pct_calculator <- function(expr_mtx_1, expr_mtx_2, min.pct, all.genes, own.genes = "no"){
  genes <- intersect(row.names(expr_mtx_1), row.names(expr_mtx_2))
  
  final_genes <- genes
  
  #print(final_genes)
  
  count_mtx_1 <- 1 - (rowCounts(expr_mtx_1, value = 0) /length(colnames(expr_mtx_1)))
  names(count_mtx_1) <- row.names(expr_mtx_1)
  
  count_mtx_2 <- 1 - (rowCounts(expr_mtx_2, value = 0) /length(colnames(expr_mtx_2)))
  names(count_mtx_2) <- row.names(expr_mtx_2)
  
  for (i in 1:length(genes)){
    
    
    if (all.genes == T){
      "all.genes"
    }
    
    else if(count_mtx_1[genes[i]] < min.pct & count_mtx_2[genes[i]] < min.pct ){
      final_genes <- final_genes[-(which(final_genes == genes[i]))]
      
    }
    
  }
  
  if(own.genes == "no"){
    count_mtx_1 <- count_mtx_1[final_genes]
    count_mtx_2 <- count_mtx_2[final_genes]
  }
  
  else{
    count_mtx_1 <- count_mtx_1[own.genes]
    count_mtx_2 <- count_mtx_2[own.genes]
  }
  
  
  return(list(count_mtx_1, count_mtx_2))
}

#Seurat way of getting LogFC 
mean.fxn <- function(x) {
  return(log(x = (rowSums(x = expm1(x = x)) + 1)/NCOL(x), base = 2))
  
}

commonID <- function(captures, clusters){
  
  output <- data.frame(rep("id", length(names(captures))), row.names= names(captures))
  
  #If we have more than 1 meta data variable
  if( length(clusters[1,]) > 1){
    for(i in 2:length(clusters[1,])){ output <- cbind(output, data.frame(rep("id", length(names(captures))), row.names= names(captures)) ) }
  }
  
  name_vector <- c()
  
  for (i in 1:length(captures)){
    cluster_c <- as.matrix(clusters[captures[[i]],])

    for (j in 1:length( clusters[1,] )){
  
      occurence <- table(cluster_c[,j])
      choice <- names(occurence)[which(occurence == max(occurence))]
      
      if (length(choice) > 1){
        # choice <- paste0(choice, collapse = "/")
        choice <- choice[1]
        
      }
      
      #If no cells are closer to it than another interpolated point
      if(sum(table(cluster_c)) ==0){
        choice <- "None"
      }
      
      name_vector <- append(name_vector, choice)
      # output[names(captures)[i],1] <- choice
      output[i,j] <- choice
      
    }
    
  }
  
  
  for (j in 1:length( clusters[1,] )){
    
    none_which <- which(output[,j] == "None")
    
    for (i in none_which){
      output[i,j] <- unique( c( output[i-1,j], output[i+1,j] ) )[1]
    }
    
  }
  

  colnames(output) <- colnames(clusters)
  
  return(output)
}


nodePseudotime <- function(sc_obj,pseudotime_ID,meta_IDs, total_nodes, prefix){
  
  cell_pseudotime <- matrix(sc_obj[[pseudotime_ID]], dimnames =list(colnames(sc_obj), "pseudotime"))
  cell_pseudotime <- sc_obj[[pseudotime_ID]]
  names(cell_pseudotime) <- colnames(sc_obj)
  
  cell_ID <- data.frame(sc_obj[[ meta_IDs[1] ]], row.names = colnames(sc_obj) )
  
  #If the user has multiple meta data columns 
  if(length(meta_IDs) > 1){
    
    for(i in 2:length(meta_IDs)){
      
      cell_ID <- cbind(cell_ID, data.frame( sc_obj[[ meta_IDs[i] ]], row.names = colnames(sc_obj) ) )
      
    }
    
  }
  
  colnames(cell_ID) <- meta_IDs

  output <- list()
  
  # cell_pseudotime <- cell_pseudotime[ order(cell_pseudotime), ]
  cell_pseudotime_ordered <- sort(cell_pseudotime)

  max_pseudo <- max(cell_pseudotime, na.rm=T)
  
  
  #Get pseudotime of interpolated points
  interval <- max_pseudo / (total_nodes + 1)
  
  node_pseudotime_vector <- seq(interval, max_pseudo - interval, interval)
  
  node_names <- paste0(prefix, "_" ,seq(1, total_nodes, 1) )
  
  #make the tree
  node_tree <- data.frame("Child", "Parent")
  
  #Don't want to make the last interpolated point a parent so we'll stop the loop before then
  for (node in 1:(length(node_names)-1)){
    node_tree <- rbind( node_tree , c( node_names[node+1], node_names[node]) )
  }
  
  # node_pseudotime <- data.frame(node_pseudotime_vector,row.names=node_names)
  # colnames(node_pseudotime) <- "pseudotime"
  
  node_pseudotime <- node_pseudotime_vector
  names(node_pseudotime) <- node_names

  #Get most common cell type for each inteprolated point
  captures <- cell_capture(cell_pseudotime_ordered, node_pseudotime)
  node_cluster <- commonID(captures, cell_ID)
  # node_cluster <- node_cluster[,1]
  # names(node_cluster) <- node_names
  
  colnames(node_tree) <- node_tree[1, ]
  
  #Remove the initial values we used to create the data frame
  node_tree <- node_tree[-1, ]
  
  output[["ID"]] <- node_cluster
  output[["tree"]] <- node_tree
  output[["cell_pseudotime"]] <- cell_pseudotime
  output[["node_pseudotime"]] <- node_pseudotime
  
  
  return(output)
}

ROC_fill_mtx <- function(penalty_mtx, cut_type, method){
  
  penalty_mtx_cp <- penalty_mtx 
  
  #Make the score matrix
  traceback_mtx <- matrix(0, nrow = length(row.names(penalty_mtx)), ncol = length(colnames(penalty_mtx)))
  
  traceback_mtx[1,1] <- penalty_mtx[1,1]
  
  traceback_mtx[1,] <- own_fib_ROC(penalty_mtx[1,])
  
  traceback_mtx[,1] <- own_fib_ROC(penalty_mtx[,1])
  
  #Fill in the score matrix
  for (i in 2:length(traceback_mtx[,1])){
    
    for (j in 2:length(traceback_mtx[1,])){
      choice1 <- traceback_mtx[(i-1), j] + ( (penalty_mtx[i,j] - (penalty_mtx[(i-1), j] ) ))
      choice2 <- traceback_mtx[i,(j-1)] + ((penalty_mtx[i,j] - (penalty_mtx[i, (j-1)]) )  )
      choice3 <- traceback_mtx[(i-1),(j-1)] + (penalty_mtx[i,j] - (penalty_mtx[(i-1), (j-1)] ) )
    
      choice1 <- traceback_mtx[(i-1), j] + ( (penalty_mtx[i,j]))
      choice2 <- traceback_mtx[i,(j-1)] + ((penalty_mtx[i,j] )  )
      choice3 <- traceback_mtx[(i-1),(j-1)] + (penalty_mtx[i,j]  )
  
      final_choice <- min(c(choice1, choice2, choice3))
      
      traceback_mtx[i,j] <- final_choice
      
    }
    
  }
  row.names(traceback_mtx) <- row.names(penalty_mtx)
  colnames(traceback_mtx) <- colnames(penalty_mtx)
  
  return(traceback_mtx)
}


bootstrap_path <- function(input, bootstrap_iterate, bootstrap_len){
  
  result <- c()
  
  for(i in bootstrap_iterate){
    
    output <- input[sample(seq(1, length(input), 1), size = bootstrap_len, replace = T)]
    
    result <- append(result, mean(output))
  }
  
  return(mean(result))
  
}

#As some paths will have more connections than others (i.e. longer paths), we will bootstrap each path to be the same length

bootstrap_pathfind <- function(sequence_1,sequence_2, similarity_method = "spearman",threshold_method, bootstrap_iterate = 10){
  
  if(threshold_method == "mean"){
    cutoff_metric <- function(x){return(mean(x))}
    
  }
  
  else if (threshold_method == "median"){
    cutoff_metric <- function(x){return(median(x))}
    
  }
  
  penalty_mtx <- dis_mtx_calculator(sequence_1, sequence_2, similarity_method)

  #Get the most minimum start and end points
  #Get the last possible start and end points
  #Bootstrap paths (n parameter times, default 100) for all the possible start and end points and get their average score
  #Use the average scores to get the threshold for choosing where the start and end point are (like i've done already)
  
  #Find the point of best (i.e. lowest dissimilarity) initial alignment of the two processes
  choice <- c(min(penalty_mtx[,1]), min(penalty_mtx[1,]))
  
  ###########       NOTE        #################
  #Won't work if there are points with the same minimum value
  choice_index <- c( which( penalty_mtx[,1] == min(penalty_mtx[,1]) ), 
                     which( penalty_mtx[1,] == min(penalty_mtx[1,]) ) )
  
  start_index_choice <- which(choice == min(choice))
  
  if(start_index_choice == 1){
    min_slice <-penalty_mtx[,1]
  }
  
  else if(start_index_choice == 2){
    min_slice <- penalty_mtx[1,]
  }
  
  start_min <- c(1,1)
  start_min[2] <- choice_index[start_index_choice]
  
  start_min[1] <- find_optimal(min_slice, start_min[1], position = "start")
  
  #Find the point of best (i.e. lowest dissimilarity) last alignment of the two processes
  last_row <- length(penalty_mtx[,1])
  last_column <- length(penalty_mtx[1,])
  
  # Find minimum at end to start from
  choice <- c(min(penalty_mtx[,last_column]), min(penalty_mtx[last_row,]))
  
  choice_index <- c( which( penalty_mtx[,last_column] == min(penalty_mtx[,last_column]) ),
                     which( penalty_mtx[last_row,] == min(penalty_mtx[last_row,]) ) )
  
  end_index_choice <- which(choice == min(choice))
  
  if(end_index_choice == 1){
    min_slice <-penalty_mtx[,last_row]
  }
  
  else if(end_index_choice == 2){
    min_slice <- penalty_mtx[last_column,]
  }
  
  
  end_min <- c(last_row, last_column)
  
  end_min[1] <- choice_index[end_index_choice]
  
  end_min[2] <- find_optimal(min_slice, end_min[1], position = "end")

  #last_start_end <- find_start_end(penalty_mtx, start_min, end_min, threshold_method)
  
  print("Start and ends")
  print(start_min)
  print(end_min)
  print("")
  # print(last_start_end)
  
  #We add the [1] to the end of the index choices, just in case the alignment runs from the very beginning of both sequences,
  # runs to the very end of both sequences, or both.

  start_locs <- seq( start_min[1], start_min[ 2 ], 1 )
  
  end_locs <- seq( end_min[ 1 ], end_min[ 2 ],1 )
  
  print(start_locs)
  print(end_locs)
  
  #Iterate through the combinations
  
  paths <- list()
  
  for(i in 1:length(start_locs)){
    start_current <- start_locs[i]
    
    for(j in 1:length(end_locs)){
      end_current <- end_locs[j]
      
      start_loc <- c(start_current,start_current)
      end_loc <- c(end_current,end_current)
      
      start_loc[3-start_index_choice] <- 1
      
      #If the first sequence is the one we are cutting
      if (end_index_choice == 1){
        end_loc[3-end_index_choice] <- length(penalty_mtx_cut[1,])
      }
      
      #If the second sequence is the one we are cutting
      else{
        end_loc[3-end_index_choice] <- length(penalty_mtx_cut[,1])
      }
      
      
      sequence_1_cut <- sequence_1[,start_loc[1]:end_loc[1]]
      sequence_2_cut <- sequence_2[,start_loc[2]:end_loc[2]]
      
      penalty_mtx <- dis_mtx_calculator(sequence_1_cut, sequence_2_cut, similarity_method)
      
      path <- pathfind(penalty_mtx, cut_type = "both", threshold_method)
    
      # if(start_current == 3 & end_current == 53){
      #   return(path)
      # }
      
      paths[[paste0(start_current, "_", end_current)]] <- path[,3]
      

    }
    
  }
  
  print("Possible start and ends")
  print(start_locs)
  print(end_locs)
  print("")
  
  score_mtx <- matrix(0, nrow = length(start_locs), ncol = length(end_locs))
  
  colnames(score_mtx) <- as.character(end_locs)
  row.names(score_mtx) <- as.character(start_locs)
  
  bootstrap_n <- max(lengths(paths))
  
  for(i in names(paths)){
    
    current_path <- paths[[i]]
    
    current_index <- str_split(i, "_")[[1]]
    
    score_mtx[current_index[1],current_index[2]] <- bootstrap_path(current_path, bootstrap_iterate,  bootstrap_n)
    
  }
  
  #Run the new score_mtx through the find_optimal function 
  

  #Return the path which had the lowest average score across the bootstraps
  
  ###########       NOTE        #################
  #Won't work if there are points with the same minimum value
  
  ###New section start###
  
  #We create a threshold around the current chosen start/end combination based on the mean or median distance between the bootstrapped path scores.
  #We then select the path that falls within that threshold and maximises being the furthest back at the start and the furthest forward at the end

  #In the case that the score mtx is 1D
  if(1 %in% dim(score_mtx)){
    score_vector <- c()
    
    for (value in score_mtx){score_vector <- append(score_vector, value)}
    
    threshold_value <- min(score_mtx) + cutoff_metric(dist(score_vector))
    
    
  }
  
  else{
    threshold_value <- min(score_mtx) + cutoff_metric(dist(score_mtx))
  }

  valid_paths_index <- which(score_mtx <= threshold_value, arr.ind = TRUE)
  
  valid_path_starts <- as.integer(row.names(score_mtx)[valid_paths_index[,1]])
  valid_path_ends <-  as.integer(colnames(score_mtx)[valid_paths_index[,2]])
  
  valid_path <- which(valid_path_ends - valid_path_starts == max(valid_path_ends - valid_path_starts))
  
  row_index <- valid_path_starts[valid_path]
  col_index <- valid_path_ends[valid_path]
  ###New section end###
  
  # row_index <- row.names(score_mtx)[which(rowMins(score_mtx) == min(rowMins(score_mtx)))]
  # col_index <- colnames(score_mtx)[which(colMins(score_mtx) == min(colMins(score_mtx)))]
  
  start_loc <- c(row_index,row_index)
  end_loc <- c(col_index,col_index)
  
  print("why")
  print(valid_path)
  print(start_loc)
  print(end_loc)
  print("")
  
  start_loc[3-start_index_choice] <- 1
  
  #If the first sequence is the one we are cutting
  if (end_index_choice == 1){
    end_loc[3-end_index_choice] <- length(penalty_mtx_cut[,1])
  }
  
  #If the second sequence is the one we are cutting
  else{
    end_loc[3-end_index_choice] <- length(penalty_mtx_cut[1,])
  }
  
  print(start_loc)
  print(end_loc)
  
  sequence_1_cut <- sequence_1[,start_loc[1]:end_loc[1]]
  sequence_2_cut <- sequence_2[,start_loc[2]:end_loc[2]]
  
  penalty_mtx <- dis_mtx_calculator(sequence_1_cut, sequence_2_cut, similarity_method)
  
  path <- pathfind(penalty_mtx, cut_type = "both", threshold_method)
  
  return(list(path, score_mtx, paths))
  
}

#Redo this function so it can take the start and end seperately. Will make the function simpler. 
#Also add a for statement before entering to make sure the very start and very end are not the initial start and ends

#Find optimal starting point, taking into consideration the correlation landscape of the error surface
find_start_end <- function(penalty_mtx, start, end, method){
  
  if(method == "mean"){
    cutoff_metric <- function(x){return(mean(x))}
    
  }
  
  else if (method == "median"){
    cutoff_metric <- function(x){return(median(x))}
    
  }
  
  #If the minimum score start point is at the very first position of the matrix, there is no need to do this step
  
  if(start[1] != 1 | start[2] != 1 ){
    #Need to find where to take the slice
    #determined by finding the maximum index based off the context of the number of metaCells in each condition
    index_1 <- which(start != 1)

    
    #get rate of change for the start
    if(index_1 == 1){
      slice_1 <- penalty_mtx[,1]
      
    }
    
    else{
      slice_1 <- penalty_mtx[1,]
    }
    
    distances_1 <- c()
    
    for (i in 1:length(slice_1)){
      distances_1 <- append(distances_1, abs(slice_1[start[index_1]] - slice_1[i]) )
    }
    
    change_v1 <- c()
    
    for(i in 2:length(slice_1)){
      change_v1 <- append(change_v1, abs(slice_1[i]-abs(slice_1[i-1])))
      
    }
    
    #get histogram of distances from the current start point
    hist_distance <- hist(distances_1)
    
    choices_1 <- which(distances_1 < hist_distance$breaks[2])
    
    cutoff_1 <- cutoff_metric(change_v1)

    selection_1 <- slice_1[1:start[index_1]]
    
    start[index_1] <- min(which(distances_1 < cutoff_1))
    start[index_1] <- min(choices_1)
  }
  
  #If the minimum score end point is at the very last position of the matrix, there is no need to do this step
  
  if(end[1] != 1 | end[1] != 1){
    index_2 <- which(end == min(end))
    index_2 <- which( c(end[1]/length(penalty_mtx[,1]), end[2]/length(penalty_mtx[1,])) == min(  c(end[1]/length(penalty_mtx[,1]), end[2]/length(penalty_mtx[1,]))  )     )[1]
    
    if(index_2 == 1){
      slice_2 <- penalty_mtx[,length(penalty_mtx[1,])]
      
    }
    
    else{
      slice_2 <- penalty_mtx[length(penalty_mtx[,1]),]
    }
    
    distances_2 <- c()
    
    for (i in 1:length(slice_2)){
      distances_2 <- append(distances_2, abs(slice_2[end[index_2]] - slice_2[i]) )
    }
    
    change_v2 <- c()
    
    for(i in 2:length(slice_2)){
      change_v2 <- append(change_v2, abs(slice_2[i]-abs(slice_2[i-1])))
      
    }
    
    cutoff_2 <- cutoff_metric(change_v2)
    
    hist_distance <- hist(distances_2)
    
    choices_2 <- which(distances_2 < hist_distance$breaks[2])
    
    end[index_2] <- max(which(distances_2 < cutoff_2))
    end[index_2] <- max(choices_2)
  }
  
  return(list(start, end))

}

find_optimal <- function(score_slice, point,position, method = "mean"){
  
  #point - numeric vector giving the row (first sequence) and column (second sequence) position of the current optimal point for alignment beginning/end
  #Position - string which tells the function whether we are trying to find the optimal start or end point
  
  if(method == "mean"){
    cutoff_metric <- function(x){return(mean(x))}
    
  }
  
  else if (method == "median"){
    cutoff_metric <- function(x){return(median(x))}
    
  }
  
  if(position == "start"){
    direction_metric <- function(x){return(min(x))}
    
  }
  
  else if(position == "end"){
    direction_metric <- function(x){return(max(x))}
    
  }
  
  #Find the distance between the score of the current optimal point and all the other possible points
  distances_2 <- c()
  for (i in 1:length(score_slice)){
    distances_2 <- append(distances_2, abs(score_slice[point] - score_slice[i]) )
  }
  
  #Find the rate of change in score as we go across the possible optimal points
  change_v2 <- c()
  for(i in 2:length(score_slice)){
    change_v2 <- append(change_v2, abs(score_slice[i]-abs(score_slice[i-1])))
    
  }
  
  cutoff_2 <- cutoff_metric(change_v2)
  
  hist_distance <- hist(distances_2)
  
  choices_2 <- which(distances_2 < hist_distance$breaks[2])
  
  point <- direction_metric(choices_2)
  
  return(point)
  
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
    
    result <- output[(number-1)] + (number_vector[number]) 
    
    output <- append(output, result)
  }
  
  return(output)
}

#Traceback from the lowest point on the last row or column
#If you do that it follows the path of just moving through the minimum like in my first way of finding the path
traceback_pathfind <- function(penalty_mtx, cut_type, method){
  
  penalty_mtx_cp <- penalty_mtx
  error_surface <- ROC_fill_mtx(penalty_mtx, cut_type, method)
  
  last_row <- length(error_surface[,1])
  last_column <- length(error_surface[1,])
    
  factor <- -1 
  
  penalty_mtx_cut <- penalty_mtx[row.names(error_surface), colnames(error_surface)]
  
  #Find position to start tracing back from
  loc <- c(last_row, last_column)
  
  match_frame <- data.frame(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]])
  
  score_vector <- c(penalty_mtx_cut[loc[1], loc[2]])
  
  #Perform first traceback move to intialise the output variables
  diag <- error_surface[(loc[1]+factor), (loc[2]+factor)]
  vertical <- error_surface[(loc[1]+factor), loc[2]]
  horizontal <- error_surface[loc[1], (loc[2]+factor)]
  
  options <- c(diag, vertical, horizontal)
  
  move_choice <- which( options == min(options) )
  
  previous_choice = move_choice
  
  if (1 == move_choice){
    loc <- loc + factor
  }
  
  else if (2 == move_choice){
    loc[1] <- loc[1] + factor
  }
  
  else{
    loc[2] <- loc[2] + factor
  }
  
  match_frame <- rbind(match_frame, c(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
  score_vector <- append(score_vector,penalty_mtx_cut[loc[1], loc[2]])
  
  x <- T
  counter <- 1
  while (x == T){
    #Three movement options
    #1. Diagonally - a 1 to 1 match
    #2 vertically - a many to 1 match 
    #3 horizontally a many to 1 match (but to a different tree from above)
    
    diag <- error_surface[(loc[1]+factor), (loc[2]+factor)]
    vertical <- error_surface[(loc[1]+factor), loc[2]]
    horizontal <- error_surface[loc[1], (loc[2]+factor)]
    
    options <- c(diag, vertical, horizontal)
    
    move_choice <- which( options == min(options, na.rm=T) )
    
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
    }
    
    else if (2 == move_choice){
      loc[1] <- loc[1] + factor
    }
    
    else{
      loc[2] <- loc[2] + factor
    }
    
    match_frame <- rbind(match_frame, c(row.names(error_surface)[loc[1]], colnames(error_surface)[loc[2]]))
    score_vector <- append(score_vector,penalty_mtx_cut[loc[1], loc[2]])
    
    previous_choice <- move_choice
    
    #Because we've decided the start point already, we leave when we are one interpolated point away from the end on either process
    if (loc[1] == 2 | loc[2] == 2){

      #If we've reached the end of the alignment already then we don't want to add it again
      if(sum(loc) != 2){
        
        loc[1] <- loc[1] + factor
        loc[2] <- loc[2] + factor
        
        final_addition <- cbind(rep(row.names(error_surface)[loc[1]:1], loc[2]) , rep(colnames(error_surface)[loc[2]:1], loc[1]))
        score_vector <- append(score_vector,penalty_mtx_cut[loc[1]:1, loc[2]:1])
        
        colnames(final_addition) <- colnames(match_frame)
        
        match_frame <- rbind(match_frame, final_addition)
      }
      
      
      x <- F
    }

    counter <- counter + 1
    
  }
  
  output <- data.frame(rev(match_frame[,1]), rev(match_frame[,2]), rev(score_vector))
  colnames(output) <- c("row", "col")
  

  return(output)
  
}

pathfind <- function(penalty_mtx, cut_type = "both", method = "mean", dis_metric = "spearman"){
  
  output <- traceback_pathfind(penalty_mtx, cut_type, method)
  
  if(mean(output[,3]) > 0.5 & dis_metric %in% c("spearman", "pearson")){
    warning(paste0("
             ################################################################
             IMPORTANT MESSAGE TO FOLLOW. PLEASE READ THE MESSAGE THOROUGHLY:
             ################################################################
            
             
             The mean dissimilarity value of the alignment path is ",mean(output[,3]), "
             This, likely, means that there is no shared process between the two datasets. Continue at your own risk    "
            ))
    
  }
  
  return(output)
}

#Finds the median of the matched and unmatched nodes. Then calculates the difference between each of the matched pairs and these two thresholds
#If the match is closer to median of the unmatched nodes than the matched nodes, then we designate it to be 'cut' from the alignment
cut_deviate <- function(alignment, dis_mtx, method){
  
  threshold_output <- get_thresholds(alignment, dis_mtx, method)
  
  median_match <- threshold_output[[1]]
  median_unmatch <- threshold_output[[2]]
  match_cost <- threshold_output[[3]]
  
  
  score_frame <- data.frame( match_median = abs(match_cost[["score"]] - median_match), unmatch_median = abs(match_cost[["score"]] - median_unmatch) )
  
  final_match <- cbind(alignment, rep("match", length(alignment[,1])))
  
  final_match[which(score_frame[,1] > score_frame[,2]),4] <- "cut"
  
  colnames(final_match) <- c("row", "col", "Score","Status")
  
  return(final_match )
  
}

#Gets the thresholds for the cut_deviate() function
get_thresholds <- function(alignment, penalty_mtx, method){
  match_cost <- c()
  
  if(method == "mean"){
    score_method <- function(x){ return( mean(x) ) }
    
  }
  
  else if(method == "median"){
    score_method <- function(x){ return( median(x) ) }
  }
  
  for (match in 1:length(alignment[,1])){
    loc <- c(which(row.names(penalty_mtx) == alignment[match,1]), which(colnames(penalty_mtx) == alignment[match,2]))
    match_cost <- append(match_cost, penalty_mtx[loc[1], loc[2]])
    
  }
  match_cost <- as.matrix(match_cost)
  
  match_cost <- data.frame(score = match_cost, group = rep( "match", length(match_cost) ))
  
  match_cost_unmatched <- penalty_mtx
  
  rows <- as.numeric(str_replace(alignment[,1], "Y", ""))
  cols <- as.numeric(str_replace(alignment[,2], "X", ""))
  
  coord <- as.matrix(data.frame(rows, cols))
  
  match_cost_unmatched[coord] <- NA
  
  match_cost_unmatched <- as.vector(match_cost_unmatched)
  match_cost_unmatched <- data.frame(score = match_cost_unmatched[which(is.na(match_cost_unmatched) == F )], group = rep("unmatch", length(match_cost_unmatched[which(is.na(match_cost_unmatched) == F )])))
  
  median_match <- score_method(match_cost[["score"]])
  median_unmatch <- score_method(match_cost_unmatched[["score"]])
  
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

#Function which estimates variable 2 for set2 based off of variable 2 for set 1 and weighted by the distance between variable 1 for set 1 and set 2.
cellAlign_function <- function(set1_variable1, set1_variable2, set2_variable1 ,window_size){
  set1_variable1_mtx <- t( replicate(length(set2_variable1), set1_variable1) )
  set2_variable1_mtx <- replicate( length(set1_variable1) , set2_variable1)
  
  #weight. N x C matrix 
  weight <-  exp( -1 * ( (set1_variable1_mtx - set2_variable1_mtx)^2 / ((window_size)^2) )) 
  
  row.names(weight) <- names(set2_variable1)
  
  #NxC dot 1xN
  
  set2_variable2 <- weight %*% set1_variable2
  
  normalizer <- 1/rowSums(weight)
  
  set2_variable2 <- set2_variable2 * normalizer
  names(set2_variable2) <- row.names(weight)
  
  return(set2_variable2)
}

###
#Splits the alignment up into segments split by where there is a cut in the alignment space.
chunk_generator <- function(tree){
  
  output <- list()
  
  previous <- tree[1,4]
  
  vector_c <- c(1)
  
  list_counter <- 1
  
  for (i in 2:length(tree[,1])){
    
    current <- tree[i,4]
    
    #If we've reached a match that's been cut:
    #End the current segment and start a new one
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

chunk_node <- function(cond1_pseudo, cond2_pseudo, tree){
  
  #get list of nodes which are only present on cut parts
  cut_nodes <- c(unique(tree[tree[,4] == "cut",1]), unique(tree[tree[,4] == "cut",2]))
  
  #Get list of nodes which are only present on matched parts
  match_nodes <- c(unique(tree[tree[,4] == "match",1]), unique(tree[tree[,4] == "match",2]))
  
  #Some nodes may be multi matched to a matched and cut node, we want to remove the ones that are matched to a cut node 
  cut_nodes <- cut_nodes[ which( !(cut_nodes %in% match_nodes)) ]
  
  pseudo_list <- list(cond1_pseudo, cond2_pseudo)
  
  cut_tree <- tree[which(tree[,4] != "cut"),]
  
  #Initial alignment
  start_pseudo_choice <- c(pseudo_list[[1]][cut_tree[1,1]], pseudo_list[[2]][cut_tree[1,2]])
  
  #The trajectory which has the lowest initial pseudotime needs to have their pseudotimes moved up to come in line with the point of 
  #initial alignment with the other trajectory 
  
  choice_index <- which(start_pseudo_choice == min(start_pseudo_choice))[1]
  
  print(choice_index)
  print(start_pseudo_choice)
  
  diff <- abs(start_pseudo_choice[1] - start_pseudo_choice[2])
  
  #If we have unaligned sections at the beginning, it means the initial alignment will accidently overlap with it
  pseudo_list[[choice_index]] <- pseudo_list[[choice_index]] + diff
  
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
        
        pseudo_choice <- c( pseudo_list[[1]][tree_chunk[item,1]], pseudo_list[[2]][tree_chunk[item,2]] )
        
        if(pseudo_choice[1] != pseudo_choice[2]){
          
          choice_index <- which(pseudo_choice != seperate_(pseudo_choice))
          
          diff <- pseudo_choice[3-choice_index] - pseudo_choice[(choice_index)]
          
          node_position <- which(names(pseudo_list[[choice_index]]) == tree_chunk[item,choice_index])[1]
          print(node_position)
          
          pseudo_list[[choice_index]][node_position:(length(pseudo_list[[choice_index]]))] <- pseudo_list[[choice_index]][node_position:(length(pseudo_list[[choice_index]]))] + diff
          
          seperate_ <- function(input){min(input)}
          
        }
        
        #if the item graph is a subset of another larger graph
        if (!(tree_chunk[item,1] %in% unique_nodes) | !(tree_chunk[item,2] %in% unique_nodes)){
          
          #get index of the common node for the multi-matching
          common_index <- which(!(c(tree_chunk[item,1], tree_chunk[item,2]) %in% unique_nodes))
          
          #print(tree_chunk[item,common_index])
          
          lower_squish <- pseudo_list[[common_index]][tree_chunk[item,common_index]]
          
          #Find the upper squish value
          pseudo_position <- which(names(pseudo_list[[common_index]]) == tree_chunk[item,common_index]) 
          
          tree_chunk_nodes <- tree_chunk[ which(tree_chunk[,common_index] == names(pseudo_list[[common_index]])[pseudo_position+1])[1], 1:2]
          
          #upper_squish <- min( c( pseudo_list[[1]][ tree_chunk_nodes[[1]] , 1 ], pseudo_list[[2]][ tree_chunk_nodes[[2]] , 1 ] ) )
          squish_nodes <- tree_chunk[ which(tree_chunk[,common_index]==tree_chunk[item,common_index]) ,(3-common_index)]
          print(paste0("squish: ", squish_nodes))
          
          #Finds the nodes that occur after the final matching in this particular multimatch
          after_nodes <-c( names(pseudo_list[[3-common_index]])[ ( which( names(pseudo_list[[3-common_index]]) == squish_nodes[length(squish_nodes)] ) +1 ) ],
                           names(pseudo_list[[common_index]])[pseudo_position+1])
          
          #If we're doing a multi-align at the end of the process, there will be no nodes that occur after it to squish between
          #If one of the values is NA, the sum of after nodes will be more than 0, thus we enter this if statement
          
          if (sum(is.na(after_nodes)) > 0){
            print("aaaah")
            #The new upper squish becomes the lower squish the average distance between the previous aligned matches
            upper_squish <- lower_squish + median( abs( c(pseudo_list[[common_index]][1:length(pseudo_list[[common_index]])-1]) - c(pseudo_list[[common_index]][2:length(pseudo_list[[common_index]])], pseudo_list[[common_index]][1]) ) )
            
            print(lower_squish)
            print(upper_squish)
            # upper_squish <- lower_squish + mean( abs( pseudo_list[[common_index]] - mean(pseudo_list[[common_index]]) ) )
          }
          
          else{
            upper_squish <- min( c(  pseudo_list[[3-common_index]][after_nodes[1]], pseudo_list[[common_index]][after_nodes[2]] ), na.rm = T )
            
          }
          
          #We don't want the last point in the multi match to have the pseudotime of the next pair, so we take it's pseudotime and move slightly back from it
          upper_squish <- upper_squish - ((upper_squish-lower_squish)/ length(squish_nodes))
          
          #print(paste0("upper squish", upper_squish))
          #print(paste0("lower squish", lower_squish))
          
          output <- squish_function(values = pseudo_list[[3-common_index]][squish_nodes], upper_squish, lower_squish)
          output <- scale_function(pseudo_list[[3-common_index]][squish_nodes], upper_squish, lower_squish)
          #print(paste0("output: ", output))
          diff <- abs(pseudo_list[[3-common_index]][ squish_nodes[length(squish_nodes)]  ] - output[length(output)])
          
          pseudo_list[[3-common_index]][squish_nodes] <- output
          
          #Don't want to pulldown if we've already reached the end of the process we're moving
          if (!( !( tree[length(tree[,1]),1] %in% squish_nodes ) | !( tree[length(tree[,1]),2]  %in% squish_nodes ) )){
            pulldown_index <- which(names(pseudo_list[[3-common_index]]) == squish_nodes[length(squish_nodes)]) + 1
            
            pseudo_list[[3-common_index]][ pulldown_index : length(pseudo_list[[3-common_index]]) ] <- pseudo_list[[3-common_index]][ pulldown_index : length(pseudo_list[[3-common_index]]) ]  - diff
            
          }
          
          item <- item + length(squish_nodes) - 1
          
        }
        item <- item + 1
        
      }
      
    }
    
    #If the current chunk is unaligned/cut, we want to make sure it does not overlap with cells from the other condition (in terms of pseudotime). 
    #Therefore, the next aligned nodes should be moved up to the node in the match with the highest pseudotime, so the cut node remains seperate.
    else{
      seperate_ <- function(input){max(input)}
    }
    
  }
  
  matrix_1 <- pseudo_list[[1]]
  matrix_1 <- data.frame("pseudotime" = matrix_1, alignment = rep("align", length(matrix_1)))
  
  matrix_2 <- pseudo_list[[2]]
  matrix_2 <- data.frame("pseudotime" = matrix_2, alignment = rep("align", length(matrix_2)))
  
  matrix_1[which(!(row.names(matrix_1) %in% match_nodes)),2] <- "noalign"
  matrix_2[which(!(row.names(matrix_2) %in% match_nodes)),2] <- "noalign"
  
  return(pseudo_list <- list(condition_1 = matrix_1, condition_2 = matrix_2))
  
}

#scales values
scale_function <- function(vals, scale_max, scale_min){
  #print(paste0("values scale: ",vals))
  #print(paste0("upper squish", scale_max))
  #print(paste0("lower squish", scale_min))
  
  
  
  output <- ( (vals - min(vals)) / (max(vals) - min(vals)) ) *(  (scale_max - scale_min)) + scale_min
  
  
  
    
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
  
  dis_matrix <- t(dis_matrix)
  
  print( pheatmap::pheatmap(as.matrix(dis_matrix), cluster_rows = F, cluster_cols = F, scale="none", border_color = NA) )
  
  return(dis_matrix)
}

pseudo_cell_align <- function(cell_pseudo, node_information, node_pseudo , window_size){
  
  node_alignment <- node_information$alignment
  
  node_pseudo_new <- node_information$pseudotime
  names(node_pseudo_new) <- row.names(node_information)
  
  cell_pseudo_new <- cellAlign_function(node_pseudo, node_pseudo_new, cell_pseudo, window_size)
  
  print(length(cell_pseudo))
  
  kmean_output <- cell_capture(cell_pseudo, node_pseudo)
  
  #Some nodes have no cells close to it and have to be removed
  index <- which(lengths(kmean_output) != 0)
  kmean_output <- kmean_output[names(index)]

  node_pseudo_new <- node_pseudo_new[which(names(node_pseudo_new) %in% names(index))]
  
  #If the closest node to a cell is unaligned, that cell becomes unaligned
  noalign_pseudo <- data.frame(ID = "ID", pseudotime = 0, status = "hi")
  align_pseudo <- data.frame(ID = "ID", pseudotime = 0, status = "hi")
  
  for (i in 1:length(node_pseudo_new)){
    
    item <- node_alignment[i]
    node <- names(node_pseudo_new)[i]
    
    if (item == "noalign"){
      df <- data.frame(ID = kmean_output[[node]], pseudotime = cell_pseudo_new[ kmean_output[[node]] ], status = "noalign")
      noalign_pseudo <- rbind(noalign_pseudo, df)
      
    }
    
    else{
      df <- data.frame(ID = kmean_output[[node]], pseudotime = cell_pseudo_new[ kmean_output[[node]] ], status = "align")
      align_pseudo <- rbind(align_pseudo, df)
    }
    
  }
  
  # return(list(noalign_pseudo, align_pseudo))
  
  noalign_pseudo <- noalign_pseudo[-1,]
  align_pseudo <- align_pseudo[-1,]
  
  output <- rbind(noalign_pseudo, align_pseudo)
  
  print(output)
  print(cell_pseudo)
  
  #Keep order of cells the same as it is in the sce and seurat objects
  output <- output[ order( match( output[,1], names(cell_pseudo) ) ),]
  output <- output[,-1]
  
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

PlotOutput<- function(tree_1, tree_2, alignment){
  
  cut_alignment <- alignment[alignment$Status == "match",]
  
  unaligned_nodes <- append(names(tree_1$node_pseudotime)[which( !(names(tree_1$node_pseudotime) %in% cut_alignment[,1]) )] ,names(tree_2$node_pseudotime)[which( !(names(tree_2$node_pseudotime) %in% cut_alignment[,2]) )] )
  
  plot_df <- data.frame("node" = "x", "group" = 0, "pseudotime" = 0, "condition" = "ah", "ID" = "ah")
  
  plot_counter <- 1
  
  cond_names <- c( str_split( names(tree_1$node_pseudotime)[1], pattern = "_" )[[1]][1]  ,
                   str_split( names(tree_2$node_pseudotime)[1], pattern = "_" )[[1]][1] )
  
  
  for (i in 1:length(cut_alignment[,1])){
    plot_df <- rbind(plot_df, data.frame("node" = cut_alignment[i,1], "group" = plot_counter, "pseudotime" = tree_1$node_pseudotime[cut_alignment[i,1]], "condition" = cond_names[1] , 
                                         "ID" = tree_1$ID[cut_alignment[i,1],]) )
    
    plot_df <- rbind(plot_df, data.frame("node" = cut_alignment[i,2], "group" = plot_counter, "pseudotime" = tree_2$node_pseudotime[cut_alignment[i,2]], "condition" = cond_names[2] , 
                                         "ID" = tree_2$ID[cut_alignment[i,2],]) )
    
    plot_counter <- plot_counter + 1
  }
  
  #return(plot_df)
  
  plot_df <- plot_df[-1,]
  
  print(length(unaligned_nodes))
  print(unaligned_nodes)
  print(length(seq(plot_counter, (plot_counter + length(unaligned_nodes)-1), 1)))
  print(length(c(tree_1$node_pseudotime[unaligned_nodes], tree_2$node_pseudotime[unaligned_nodes])[ is.na(c(tree_1$node_pseudotime[unaligned_nodes], tree_2$node_pseudotime[unaligned_nodes])) == F ]))
  print(c( rep(cond_names[1], length( which( grepl(cond_names[1],unaligned_nodes) == T) ) ), rep(cond_names[2], length( which( grepl(cond_names[2],unaligned_nodes) == T ) ) )))
  
  #add unaligned nodes
  unaligned_frame <- data.frame("node" = unaligned_nodes,
                                "group" = seq(plot_counter, (plot_counter + length(unaligned_nodes)-1), 1),
                                "pseudotime" = c(tree_1$node_pseudotime[unaligned_nodes], tree_2$node_pseudotime[unaligned_nodes])[ is.na(c(tree_1$node_pseudotime[unaligned_nodes], tree_2$node_pseudotime[unaligned_nodes])) == F ],
                                "condition" = c( rep(cond_names[1], length( which( grepl(cond_names[1],unaligned_nodes) == T) ) ), rep(cond_names[2], length( which( grepl(cond_names[2],unaligned_nodes) == T ) ) )),
                                "ID" = c( tree_1$ID[unaligned_nodes[which( grepl(cond_names[1],unaligned_nodes) == T)],], tree_2$ID[unaligned_nodes[which( grepl(cond_names[2],unaligned_nodes) == T)],] ) ) 
  
  plot_df <- rbind(plot_df, unaligned_frame)

  #print(max(plot_df$pseudotime))
  
  print( ggplot(data = plot_df, aes(pseudotime, condition)) + geom_point(aes(col = ID), size=3) + geom_line(aes(group = group)) + theme_classic() )#+scale_x_continuous(breaks = seq(0, max(plot_df$pseudotime), max(plot_df$pseudotime)/5)) )
  
}

PlotWindowHeatmap <- function(DE_all, DE_output, DE_genes, rowClus = F, colClus = F ){
  
  heatmap_mat <- matrix(0, nrow = length(DE_genes), ncol =  length(DE_output), dimnames = list(c(DE_genes) , paste0("window_", seq(1, length(DE_output), 1)) ) )
  
  for (i in 1:length(DE_all)){
    input_vector <- DE_all[[i]][DE_genes,"logfc"]

    heatmap_mat[,i] <-input_vector 
  }
  
  input_mat <- matrix("", nrow = length( DE_genes), ncol =  length(DE_output), dimnames = list(c( DE_genes) , paste0("window_", seq(1,  length(DE_output), 1)) ) )
  
  #Go through the original matrix and see if the gene is present in the DE table. If so, give it a *
  
  for (i in 1:length(DE_all)){
    
    for (j in 1:length(DE_genes)){
      
      if(DE_genes[j] %in% DE_output[[i]][,"gene_name"]){
        
        input_mat[j,i] <- "*"
        
      }
      
    }
    
  }
  
  
  library(gplots)
  myColor <- colorRampPalette(c("red", "white", "blue"))(200)
  return(heatmap.2(heatmap_mat, cellnote = input_mat,Rowv = rowClus, Colv = colClus, notecex = 2, tracecol = NA, col = myColor, notecol = "black", cexCol = 1.5 ,margins = c(8,8) ,dendrogram = "none", key.title = "", scale = "none",
                   cexRow = 1))
  
}

PlotAlignment <- function(alignment, score_mtx){
  cut_alignment <- alignment[alignment$Status == "match",]
  
  score_mtx[as.matrix( cbind(cut_alignment[,1], cut_alignment[,2]))] <- NA
  
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(200)
  
  breaks_plot <- seq(0, max(score_mtx, na.rm=T), (max(score_mtx, na.rm=T)/200))
  
  return(pheatmap::pheatmap(as.matrix(score_mtx), cluster_rows = F, cluster_cols = F, scale="none", border_color = NA, margins = c(8,8),  breaks = breaks_plot, color = myColor, cexCol = 0.2))
}

#Assign metaCells based on cell density across pseudotime
#DE sliding window on cells closet to matched metaCells
#If the number of cells don't reach a minimum cut off then we include the next match 
#ALternatively user defines how many unique matches to include in window 
#Window size is how many windows you want in total

TrajDE <- function(sce_obj, traj_info, matches ,n_windows, overlap, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes, own.genes = "no",test_use, correct = T){
  
  if(test_use == "t"){
    test_function <- function(x,y){t.test(x,y)}
  }
  
  else if (test_use == "wilcox"){
    test_function <- function(x,y){wilcox.test(x,y)}
    
  }
  
  else if (test_use == "perm"){
    test_function <- function(x,y){permTS(x,y)}
    
  }
  
  de_genes_output <- list()
  
  capture_output <- list()
  
  for (i in 1:2){
    #Assign cells to the closet metaCell in terms of pseudotime
    cell_embed <- traj_info[[i]]$cell_pseudotime
    names(cell_embed) <- names(traj_info[[i]]$cell_pseudotime)
    
    node_embed <- traj_info[[i]]$node_pseudotime
    names(node_embed) <- names(traj_info[[i]]$node_pseudotime)
    
    capture_output[[i]] <- cell_capture(cell_embed, node_embed)
    
  }
  
  matches <- matches[which(matches$Status == "match"),]
  
  #get window size 
  
  match_counter <- 0
  
  multi_match <- c( names( table(matches[,1]) )[ table(matches[,1])>1 ], names( table(matches[,2]) )[ table(matches[,2])>1 ] )
  
  comparisons <- list()
  
  #FInd groups for comparisons
  i <- 1
  while ( i < (length(matches[,1]) +1) ){
    
    current <- c(matches[i,1], matches[i,2])
    
    if(current[1] %in% multi_match | current[2] %in% multi_match){
      multi_index <- which(current %in% multi_match)
      
      match_counter <- match_counter + 1
      
      #Depending on which condition has the multi match we assign to the comparisons variable
      
      if (multi_index == 1){
        comparisons <- append(comparisons, list( list(current[multi_index] , matches[ which(matches[,multi_index] == current[multi_index]) , 3-multi_index])))
        
      }
      
      else{
        comparisons <- append(comparisons, list( list( matches[ which(matches[,multi_index] == current[multi_index]) , 3-multi_index] , current[multi_index] )))
        
        
      }
      
      i <- i + length(matches[ which(matches[,multi_index] == current[multi_index]) , 3-multi_index])
      
      
    }
    
    else{
      
      comparisons <- append(comparisons, list( list(matches[i,1], matches[i,2])))
      
      match_counter <- match_counter + 1
      i <- i + 1
    }
    
  }
  
  #overlap + n_windows <= matches
  matches <- length(comparisons)
  
  total_matches <- ( matches + (overlap * n_windows)  ) - overlap
  
  #number_matches_per_window <- (length(comparisons) * 2) / number_windows
  
  final_comparisons <- c()
  
  n_windows_cp <- n_windows
  
  for (i in 1:n_windows){
    length_window <- round(total_matches/n_windows_cp)
    
    n_windows_cp <- n_windows_cp - 1
    total_matches <- total_matches - length_window
    
    final_comparisons <- append(final_comparisons, length_window)
    
  }
  
  
  exp_mtx_list <- list()
  
  match_counter <- 1
  
  comparisons_cp <- comparisons
  
  for (i in 1:length(final_comparisons)){
    # print("")
    # print(i)
    move <- final_comparisons[i] - overlap
    
    cells_1 <- c()
    cells_2 <- c()
    
    for (j in 1:final_comparisons[i]){
      # print(paste0("j = ", j))
      # print(comparisons[[j]][[1]])
      cells_1 <- append(cells_1, unlist( capture_output[[1]][ comparisons[[j]][[1]] ] ) )  
      cells_2 <- append(cells_2, unlist( capture_output[[2]][ comparisons[[j]][[2]] ] ) )
      
    }
    
    expr_mtx_1 <- sce_obj[[1]]@assays@data$logcounts[ , cells_1 ]
    expr_mtx_2 <- sce_obj[[2]]@assays@data$logcounts[ , cells_2 ]
    
    exp_mtx_list[[i]] <- list(expr_mtx_1, expr_mtx_2)
    
    output <- de_metrics_return(expr_mtx_1, expr_mtx_2, logfc, p_val, min.pct, all.genes, test_function, own.genes, correct)
    
    #return(output)
    
    de_genes_output <- append(de_genes_output, list(output))
    
    match_counter <- match_counter + move
    
    comparisons <- comparisons[-(1:move)]
    
  }
  
  #return(exp_mtx_list)
  
  print(final_comparisons)
  return(list(de_genes_output, exp_mtx_list, capture_output, comparisons_cp))
}


de_metrics_return <- function(expr1, expr2, logfc_threshold, pval_threshold, min.pct, all.genes, test_function, own.genes = "no", correct = T){
  
  gene_list_og <- intersect(row.names(expr1), row.names(expr2))
  
  gene_list <- intersect(row.names(expr1), row.names(expr2))
  gene_all <- gene_list
  
  print(length(gene_list))
  
  expr1 <- expr1[gene_list,]
  expr2 <- expr2[gene_list,]
  
  temp <- cbind(expr1, expr2)
  
  temp <- rowSums(temp)
  
  if(all.genes == F & own.genes == "no"){
    gene_list <- names(temp)[which(temp > 0)]
    expr1 <- expr1[gene_list,]
    expr2 <- expr2[gene_list,]
  }
  

  pct_output <- pct_calculator(expr1, expr2, min.pct = min.pct, all.genes, own.genes )
  
  gene_list <-names(pct_output[[1]])
  
  expr1 <- expr1[gene_list,]
  expr2 <- expr2[gene_list,]
  
  gene_df <- data.frame("gene", 0, 0, 0, 0)
  colnames(gene_df) <- c("gene_name", "p_val", "logfc.change", "min.pct1", "min.pct2")
  
  condition_1 <- mean.fxn(expr1)
  condition_2 <- mean.fxn(expr2)
    
  logfc_change <- condition_1 - condition_2
  
  test_result_vector <- c()
  
  # print(length(logfc_change))
  # 
  # print(length(which(logfc_change < logfc_threshold & logfc_change > (-1*logfc_threshold) )))
  # 
  if(all.genes == F & own.genes == "no" & logfc_threshold != 0){
    gene_list <- gene_list[ -( which(logfc_change < logfc_threshold & logfc_change > (-1*logfc_threshold) ) ) ]
    
  }
  
  
  logfc_change <- logfc_change[ gene_list ]
    
  pct_output[[1]] <- pct_output[[1]][gene_list]
  pct_output[[2]] <- pct_output[[2]][gene_list]

  # print(length(gene_list))
  # 
  # print("")

  if(all.genes == T & own.genes == "no"){
    print("help")
    gene_list <- gene_all
  }
  
  else if(all.genes == F & own.genes != "no"){
    #print("somebody")
    gene_list <- own.genes
  }
  
  print(length(gene_list))
  
  for (gene in gene_list){
    test_result <- test_function(expr1[gene,], expr2[gene,])
    
    test_result_vector <- append(test_result_vector, test_result$p.value)
  }
  
  gene_df <- data.frame(gene_name = gene_list, logfc = logfc_change, p_val = test_result_vector ,adj_pval = test_result_vector, min.pct1 = pct_output[[1]], min.pct2 = pct_output[[2]])
  
 
    library(stats)
    print("corrected P-value")
    print("length before")
    print(length(gene_df[,1]))
    print("")
    print(length(gene_list_og))
    gene_df[,4] <- p.adjust(gene_df[,4], method = "bonferroni", length(gene_list_og))
    if (all.genes == F & own.genes == "no"){
      #Remove genes with non-significant p-values and logfc not in the threshold
      
      if(correct){
        gene_df <- gene_df[which(gene_df[,4] < pval_threshold),]
      }
      
      else{
        gene_df <- gene_df[which(gene_df[,3] < pval_threshold),]
      }
      
      gene_df <- gene_df[which(gene_df[,2] > logfc_threshold | gene_df[,2] < (-1*logfc_threshold) ) , ]
    }
    print("Length after")
    print(length(gene_df[,1]))
  
  
  
   
  return(gene_df) 
}


window_t_test_de <- function(expr_mtx_1, expr_mtx_2, cell_groupings, logfc_threshold, pval_threshold, all.genes, logfc.metric, min.pct){
  
  output <- list()
  for (i in 1:length(cell_groupings[[1]])){
    print(i)
    #print(cell_groupings[[1]][i])
    #print(cell_groupings[[2]][i])
    expr_mtx_1_window <- as.matrix(expr_mtx_1[ , which(colnames(expr_mtx_1) %in% cell_groupings[[1]][[i]]) ])
    expr_mtx_2_window <- as.matrix(expr_mtx_2[ , which(colnames(expr_mtx_2) %in% cell_groupings[[2]][[i]]) ])
    
    pct_output <- pct_calculator(expr_mtx_1_window, expr_mtx_2_window, min.pct = min.pct, all.genes )
    
    final_genes <-names(pct_output[[1]])
    
    expr_mtx_1_window <- expr_mtx_1_window[final_genes,]
    expr_mtx_2_window <- expr_mtx_2_window[final_genes,]
    
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
    
    gene_df <- data.frame("gene", 0, 0, 0, 0)
    colnames(gene_df) <- c("gene_name", "p_val", "logfc.change", "min.pct1", "min.pct2")
    
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
    
    gene_df <- data.frame(gene_name = gene_list, logfc = logfc_change, adj_pval = test_result_vector, min.pct1 = pct_output[[1]], min.pct2 = pct_output[[2]])
    
    
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



