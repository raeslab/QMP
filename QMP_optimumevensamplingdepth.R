# Rarefaction to optimum even sampling depth #
##############################################
# Contributers: Doris Vandeputte, Gunther Kathagen, Mireia Valles-Colomer
# function rarefy_even_sampling_depth
# with cnv_corrected_abundance_table: a copy number variation corrected abundance table with sample-identifiers as rows, copy number corrected taxa-abundances as columns
# with cell_counts_table: a table with sample-identifiers as rows, cell counts as columns 
# makes use of phyloseq function rarefy_even_depth
library(phyloseq)
rarefy_even_sampling_depth_opt <- function(cnv_corrected_abundance_table, cell_counts_table, minimum_nr_reads) 
{  
  try(if(identical(sort(row.names(cnv_corrected_abundance_table)), sort(row.names(cell_counts_table)))==FALSE) stop("cnv_corrected_abundance_table and cell_counts_table do not have the same sample names. Please check!")) 
  
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations 
  
  cell_counts_table = t(cell_counts_table[order(match(row.names(cnv_corrected_abundance_table),row.names(cell_counts_table))),,drop=F]) #order cell_counts_table as cnv_corrected_abundance_table
  
  sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  
  minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
  
  rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  if (all(rarefy_to > minimum_nr_reads)){ #if all samples above min nr reads --> no problem
    samples_to_exclude = c()
    minimum_sampling_depth_opt = minimum_sampling_depth
  } else { #try next sampling depth (and exclude x)
    # find optimal sampling depth (losing the minimum amount of samples in order to reach the n reads threshold)
    lost_max = length(which(rarefy_to < minimum_nr_reads)) #n samples that are below nr reads
    lost_due_to_exclusion = 1:(length(sampling_depths)) #vector with all n
    lost_after_rarefaction = c()
    ## check how many samples are lost after rarefaction (because they did not reach the minimum nr of reads specified by the treshold, for each sample with a low sampling depth that is removed)
    for (i in 1:(length(sampling_depths))){
      #print(i)
      minimum_sampling_depth = sort(sampling_depths, decreasing=TRUE)[length(sampling_depths)-i] # take sampling depth of next sample
      rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
      lost_after_rarefaction = c(lost_after_rarefaction, length(which(rarefy_to < minimum_nr_reads))) #for each sampling depth --> n samples that are lost
    }
    ## define the optimal sampling depth (depth at which the minimum of samples are lost in total)
    lost_in_total = lost_after_rarefaction + lost_due_to_exclusion 
    nr_samples_to_exclude = min(which(lost_in_total == min(lost_in_total))) #CAUSES PROBLEMS IF SEVERAL WITH MIN 
    
    plot(lost_in_total, xlim= c(0,lost_max), ylim=c(0,lost_max), pch=16,
         ylab = "Total nr of samples lost",
         xlab= "Nr. of excluded samples",
         main = paste0("Lost using optimal sampling depth: ", min(lost_in_total), " (instead of ", lost_max, ")"))
    
    mtext(paste(nr_samples_to_exclude, "due to exclusion and ", lost_after_rarefaction[nr_samples_to_exclude], " because they didn't make the rarefaction threshold."))
    abline(h=min(lost_in_total), v=nr_samples_to_exclude)
    minimum_sampling_depth_opt = sort(sampling_depths, decreasing = TRUE)[length(sampling_depths)-nr_samples_to_exclude] # minimum of all sampling depths at which the least samples are lost to reach the minimum nr of reads. 
    # define the samples that need to be excluded for this optimal sampling depth
    samples_to_exclude = which(sampling_depths<minimum_sampling_depth_opt)
    
  }
  
  # make QMP matrix with optimal sampling depth
  rarefy_to_opt = round(cell_counts_table * minimum_sampling_depth_opt) # nr of reads to rarefy in each sample in order to get to the optimum even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  out=NULL
  samples_included = c()
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq)){
    if (!(i %in% samples_to_exclude)){# skip the samples that need to be excluded
      if (rarefy_to_opt[i]>minimum_nr_reads){ # only include the samples that pass the treshold
        print(paste("sample",row.names(cnv_corrected_abundance_table_phyloseq)[i], "rarefied to", rarefy_to_opt[i], "reads."))
        x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to_opt[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
        out=rbind(out,x)
        samples_included = c(samples_included, i) # make a final list of samples that get included
      }
    }
  }
  print("_____________________________________________________________")
  print("Optimal sampling depth:")
  print(minimum_sampling_depth_opt)
  print("_____________________________________________________________")
  print(paste("The following",length(samples_to_exclude),"samples were excluded due to exclusion:"))
  print(row.names(cnv_corrected_abundance_table)[samples_to_exclude])
  print(paste("The following",lost_after_rarefaction[nr_samples_to_exclude],"samples were excluded because of not making the rarefaction threshold:"))
  print(row.names(cnv_corrected_abundance_table)[!(row.names(cnv_corrected_abundance_table)) %in% c(row.names(cnv_corrected_abundance_table)[samples_to_exclude],row.names(cnv_corrected_abundance_table)[samples_included])])
  
  rarefied_matrix = as.matrix(out)
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,samples_included] 
  return(QMP)
}