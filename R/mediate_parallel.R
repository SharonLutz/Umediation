mediate_parallel.unix <- function(list_of_job_args, num_jobs=getOption("mediate.jobs", parallel::detectCores() - 1)) {
  options(mc.cores = num_jobs)
  
  result <- pbapply::pblapply(list_of_job_args, perform_mediation, cl = num_jobs)
  attr(result, "dim") <- dim(list_of_job_args)
  
  return(result)
}

mediate_parallel.non_unix <-function(list_of_job_args, num_jobs=getOption("mediate.jobs", parallel::detectCores() - 1)) {
  options(cl.cores = num_jobs)
  snow::setDefaultClusterOptions(type="SOCK")
  this.cluster <- snow::makeCluster(num_jobs)
  on.exit(snow::stopCluster(this.cluster))
  #make sure reverseC is loaded on the nodes
  if(pkgload::is_dev_package("Umediation")){
    #we're in a dev environment, need to load with load_all
    snow::clusterCall(cl=this.cluster,function(){suppressMessages(library(devtools));suppressMessages(load_all())})
  } else {
    #we're being used from an installed copy of reverseC load package explicitly
    snow::clusterCall(cl=this.cluster,function(){suppressMessages(library(Umediation))})
  }
  
  result <- pbapply::pblapply(list_of_job_args, perform_mediation, cl=this.cluster)
  dim(result) = dim(list_of_job_args)
  
  return (result)
}

mediate_parallel <- function(list_of_job_args, num_jobs=getOption("mediate.jobs", parallel::detectCores() - 1)){
  pbapply::pboptions(type="timer", style=1)
  if(.Platform$OS.type == "unix") {
    result <- mediate_parallel.unix(list_of_job_args=list_of_job_args, num_jobs=num_jobs)
  } else {
    result <- mediate_parallel.non_unix(list_of_job_args=list_of_job_args, num_jobs=num_jobs)
  }
  return(result)
}