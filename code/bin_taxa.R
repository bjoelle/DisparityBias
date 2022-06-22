bin.taxa = function(taxa, nbins, max.age) {
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
    
    for(ov in overlap) {
      min.age = max(bin.min[ov], taxa$end[tax])
      max.age = min(bin.max[ov], taxa$start[tax])
      fs = rbind(fs, data.frame(sp = taxa$sp[tax], edge = taxa$edge[tax], hmin = min.age, hmax = max.age))
    }
  }
  
  FossilSim::fossils(fs)
}