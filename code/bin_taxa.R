bin.taxa = function(taxa, nbins, max.age) {
  if(nbins%%1 != 0 || nbins == 0 || nbins < 0) {
    stop("Number of bins must be a positive integer, check nbins")
  }
  
  bin.ages = seq(0, max.age, max.age/nbins)
  bin.min = bin.ages[-length(bin.ages)]
  bin.max = bin.ages[-1]
  
  fs = data.frame()
  for(tax in 1:nrow(taxa)) {
    overlap = which(bin.min <= taxa$start[tax] & bin.max >= taxa$end[tax])
    if(length(overlap) == 0) {
      warning("Taxa overlaps with no bins, check max.age")
      next
    }
    
    for(bin.index in overlap) {
      min.age = max(bin.min[bin.index], taxa$end[tax])
      max.age = min(bin.max[bin.index], taxa$start[tax])
      mid.age = (min.age + max.age) / 2
      fs = rbind(fs, data.frame(sp = taxa$sp[tax], edge = taxa$edge[tax], hmin = mid.age, hmax = mid.age))
    }
  }
  
  FossilSim::fossils(fs)
}