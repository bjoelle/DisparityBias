#' Make joined trees from several cbd subtrees
#'
#' @param n_sim number of replicates (final trees)
#' @param n_tips number of total tips
#' @param n_join number of transition events (i.e. n_join+1 subtrees are simulated)
#' @param lambda birth rate
#' @param mu death rate
#'
#' @return list of phylo trees
joined_trees = function(n_sim, n_tips, n_join, lambda, mu) {
  trees = lapply(1:n_sim, function(i) .one_joined_tree(n_tips, n_join, lambda, mu))
  trees
}

.one_joined_tree = function(n_tips, n_join, lambda, mu) {
  library(data.tree)
  
  tips = rep(round(n_tips / (n_join + 1)), n_join + 1)
  tips[1] = tips[1] + n_tips - sum(tips)
  subtrees = lapply(1:(n_join + 1), function(t) .cbd_tree(tips[t], lambda, mu))
  
  for(i in 1:length(subtrees)) {
    subtrees[[i]]$Set(tree_idx = i)
    subtrees[[i]]$Do(function(n) n$name = paste0(n$name, "_", i))
  }
  for(ii in 1:n_join) {
    ages = sapply(subtrees, function(t) t$phylo_height)
    merge_idxs = order(ages)[1:2]
    new_tree = .merge_trees(subtrees[[merge_idxs[1]]], subtrees[[merge_idxs[2]]])
    subtrees = c(subtrees[- merge_idxs], new_tree)
  }
  root_area = sample(2, 1)
  tree = .set_areas(subtrees[[1]], root_area)
  
  tree = .convert_phylo(tree)
  return(tree)
}

.cbd_tree = function(n_tips, lambda, mu)  {
  leaves = sapply(1:n_tips, function(i) Node$new(i, phylo_height = 0))
  time = 0
  node_counter = n_tips + 1
  
  while (length(leaves) > 0) {
    timestep <- rexp(1, (length(leaves) * (lambda + mu)))
    time = time + timestep
    specevent <- runif(1, 0, 1)
    new_node = Node$new(node_counter, phylo_height = time)
    node_counter = node_counter + 1
    if ((lambda/(lambda + mu)) >= specevent) {
      if (length(leaves) > 1) {
        species <- sample(length(leaves), 2)
        new_node$AddChildNode(leaves[[species[1]]])
        new_node$AddChildNode(leaves[[species[2]]])
        leaves <- c(leaves[-species], new_node)
      }
      else {
        new_node$AddChildNode(leaves[[1]])
        return(new_node)
      }
    }
    else {
      leaves <- c(leaves, new_node)
    }
  }
}

.merge_trees = function(t1, t2) {
  attacht = t1$phylo_height
  find_nodes = function(node) {
    if(node$phylo_height < attacht) return(node)
    res = c()
    for(n in node$children) res = c(res, find_nodes(n))
    return(res)
  }
  
  attach_nodes = find_nodes(t2)
  attach_node = sample(attach_nodes, 1)[[1]]
  parent = attach_node$parent
  
  parent$RemoveChild(attach_node$name)
  t1$AddChildNode(attach_node)
  t1$tree_idx = attach_node$tree_idx
  parent$AddChildNode(t1)
  
  return(t2)
}

.set_areas = function(t, root_area) {
  set_area = function(node, area) {
    node$area = area
    for(ch in node$children) {
      if(ch$tree_idx == node$tree_idx) set_area(ch, area)
      else {
        new_area = if(area == 1) 2 else 1
        set_area(ch, new_area)
      }
    }
  }
  
  set_area(t, root_area)
  return(t)
}

.convert_phylo = function(tree) {
  if(length(tree$children) < 2) {
    child = tree$children[[1]]
    phylo = data.tree::as.phylo.Node(child, heightAttribute = "phylo_height")
    phylo$root.edge = tree$phylo_height - child$phylo_height
  }
  else {
    phylo = data.tree::as.phylo.Node(tree, heightAttribute = "phylo_height")
  }
  
  phylo$area = rep(0, length(phylo$edge.length))
  labels = c(phylo$tip.label, phylo$node.label)
  set_area = function(node, phylo) {
    lab = which(labels == node$name)
    edge = which(phylo$edge[, 2] == lab)
    if(length(edge) > 0) phylo$area[edge] = node$area
    for(ch in node$children) phylo = set_area(ch, phylo)
    return(phylo)
  }
  
  phylo = set_area(tree, phylo)
  return(phylo)
}