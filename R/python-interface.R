# script for matlab connection evaluation:
# evaluate(matlab, "A=1+2;")
# data <- getVariable(matlab, c("A"))
# 
#   evaluate(matlab, "run('D:/Dropbox/Git/DPHigh-dim/JunctionTree/matlab/r_matlab.m')")
# #evaluate(matlab, "exit")
# close(matlab)
# showConnections()
# closeAllConnections()
# print("finish running matlab script !")

write_clique_python <- function(out.dir, filename, domain) {
    
  filepath <- paste(out.dir, filename, ".clique", sep = '')
  clique.lines <- readLines(filepath)
  num.of.clique <- length(clique.lines)
  clique.out <- lapply(seq_along(clique.lines), function(ii) {
    attrs <- unlist(strsplit(clique.lines[[ii]], " "))
    indices <- which(domain$name %in% attrs)
    line <- paste(paste(indices, collapse = " "), "\n", sep = "")
    return(line)
  })
  curr.path <- getwd()
  
  if (num.of.clique > 2){
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    command.file <- paste(out.dir, "/r_matlab", "_", timestamp, ".m", sep = "")
    outfile <- paste(out.dir, filename, ".in", sep = '')
    domain.out <- lapply(seq_along(domain$name), function(ii) {paste(ii, " ", domain$dsize[ii], "\n", sep = "")})
    cat(unlist(domain.out), file = outfile, sep = "")
    cat(paste(paste(rep("-", 20), collapse=""), "\n", sep = "")
        , file = outfile, sep = "", append = TRUE)
    
    cat(unlist(clique.out), file=outfile, sep = "", append = TRUE)
    cat(paste("cd(\'", curr.path, "\')\n", sep = "")
        , file = command.file
        , sep = "")
    cat(paste("addpath(\'", curr.path, "/matlab\')\n", sep = "")
        , file = command.file
        , sep = "", append=TRUE)    
    fin <- paste("\'",out.dir, filename, ".in\'", sep= '')

    merge.function.path <- paste(python.dir, "marginal-optimization.py", sep="/")
    python.command <- "python"
    python.lines <- lapply(seq.int(from = 2, to = num.of.clique - 1)
                                   , function(ii){
                                      fout <- paste("\'", out.dir, filename,  ".merge", ii, "\'", sep="")
                                      allArgs <- c(merge.function.path, fin, ii, 100, fout)
				      system2(python.command, args=allArgs, stdou=TRUE)
                                      if (file.exists(fout)) file.remove(fout)
                                    })  
  }
  merge.1 <- paste(out.dir, filename,  ".merge", 1, sep="")
  cat(paste(paste(seq.int(length(domain$name)), collapse=" "), "\n", sep="")
      , file = merge.1, append = FALSE)
  merge.last <- paste(out.dir, filename,  ".merge", num.of.clique, sep="")
  cat(unlist(clique.out), file=merge.last, sep = "", append = FALSE) 
  return()
}

compute_best_clusters <- function(out.dir, filename, domain) {
  filepath <- paste(out.dir, filename,  ".clique", sep = '')
  outpath <- paste(out.dir, filename,  ".cluster", sep = '')
  clique.lines <- readLines(filepath)
  num.of.clique <- length(clique.lines)
  total.variances <- rep(0, num.of.clique)
  for (ii in seq.int(num.of.clique)) {
    merge.file <- paste(out.dir, filename,  ".merge", ii, sep = '')
    cliques <- readLines(merge.file, warn=FALSE)
    clique.sizes <- lapply(cliques, function(x) {
      indices <- do.call(as.integer,(strsplit(str_trim(x), " ")))
      sizes <- domain$dsize[indices]
      return(sizes)
      
    })
    products <- lapply(clique.sizes, function(x) prod(x))
    total.variances[[ii]] <- 2 * (ii ^ 2) * sum(unlist(products))
  }
  best.merge <- which.min(total.variances)
  best.merge.file <- paste(out.dir, filename,  ".merge", best.merge, sep= '')
  best.cliques <- readLines(best.merge.file, warn=FALSE)
  cluster.out <- lapply(best.cliques, function(x){
    indices <- do.call(as.integer, (strsplit(str_trim(x), " ")))
    attrs <- domain$name[indices]
    line <- paste(paste(attrs, collapse=" "),"\n", sep = "")
    return(line)
  })
  cluster.file <- paste(out.dir, filename,  ".cluster", sep = '')
  cat(unlist(cluster.out), file = cluster.file, append=FALSE, sep="")
}


load_margin <- function(filepath){
  lines <- readLines(filepath)
  num.of.lines <- length(lines)
  margins <- lapply(seq_along(lines), function(ii) {
    attrs <- unlist(strsplit(lines[[ii]], " "))
    return(attrs)
  })
  return(margins)
  
}

load_clusters <- function(out.dir, filename){
  filepath <- paste(out.dir, filename,  ".cluster", sep = '')
  clusters <- load_margin(filepath) 
  return(clusters)
}

load_cliques <- function(out.dir, filename){
  filepath <- paste(out.dir, filename,  ".clique", sep = '')
  cliques <- load_margin(filepath)
  return(cliques)
}


