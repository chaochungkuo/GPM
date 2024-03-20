add_QC <- function(tag, description, file_dir, save_dir) {
  scripts  <- readLines("scRNAseq_QC_template.Rmd")
  scripts <- gsub(pattern = "SAMPLE_DIR", replace = file_dir, x = scripts )
  scripts <- gsub(pattern = "DESCRIPTION", replace = description, x = scripts )
  scripts  <- gsub(pattern = "FILE_TAG", replace = tag, x = scripts)
  filename <- paste0("scRNAseq_QC_",tag)
  writeLines(scripts, con=paste0(save_dir,"/",filename,".Rmd"))
}




###################################################################################################################################################
################################################# Helper functions ################################################################################
###################################################################################################################################################

#Function to list directories recursively up to n-depth 
#https://stackoverflow.com/questions/48297440/list-files-recursive-up-to-a-certain-level-in-r
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}


sample_name_ndepth = function(directory, n = 2){
  
  in_dir = directory
  name = ""
  
  for(i in seq(n)){
    if(i == 1){
      name = paste0(basename(in_dir), name)
    }else{
      name = paste0(basename(in_dir), "_", name)
    }
    in_dir = dirname(in_dir)
    
  }
  return(name)
}
