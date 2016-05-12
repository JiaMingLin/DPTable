# get current directory
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                       error=function(e) # works when using R CMD
                         normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
this.dir <- dirname(full.fpath)
code.dir <- paste(dirname(this.dir), "R", sep = "/")
data.dir <- paste(dirname(this.dir), "data", sep = "/")
python.dir <- paste(dirname(this.dir), "python", sep = "/")

# set the working directory to the main folder containing all the directories
wd.dir <- dirname(this.dir)
setwd(wd.dir)                                       #設定工作路徑
#options(error = recover)
options(error = NULL)                               #如果產生錯誤則不做為
source(paste(code.dir, "init.R", sep = "/"))        #啟動INIT

main <- function(args) {                            #args(參數列表)
  exp.specs <- parse_args(args)                     #把args存入exp.specs?  # FLAG類似一個開關 也許有預設 
  flag.sample <- exp.specs$flag.sample
  data.name <- exp.specs$data.name
  print(exp.specs)
  epsilon.1 <- exp.specs$epsilon.1                   #構建noisy junction tree的ep
  epsilon.2 <- exp.specs$epsilon.2                  #向marginal中加入noise的ep
  CV.thresh <- exp.specs$CV                         #numeric value in (0, 1), threshold Cramer's V value for picking correlated attributes pairs, 0.2 as default,typically choose 0.2 for weakly correlated datasets; 0.3 for highly correlated datasets
  nrun<-exp.specs$nrun                              #run次數，預設10
  flag.sim <- exp.specs$flag.sim                    #布林型 只能取TRUE或是FALSE，表示是否要生成或模擬synthetic data 預設是F
  flag.process.query <- exp.specs$flag.process.query
#   random_kways = c(4, 6)                          
#   all_kways = c(2, 3)
  random_kways = c(2)                                 #k for k way marginal????????????????
  all_kways = c(2)                                    #?????????????

#adding benchmark
op <- options(digits.secs=6)
print( paste("loading data. start time:", Sys.time(), sep=" ") )         #paste連接字符串，每一個sys.time中間用 隔開
  cat("load data: ", data.name, "\n")                                   #轉化為字符，再列印
  curr.data <- Data$new(data.name)                                              #new??????????TAG?????
  tag.sample <- as.character(epsilon.1)                        #把e1從數值行轉換成字符型(字符串不能演算 會直接輸出)
  out.dir <- paste('./output/', data.name                               #out.dir是一個顯示文檔名，門檻值，時間的字符串
                   , "_CV_", as.character(CV.thresh), "_"
                   ,format(Sys.time(), "%Y%m%d_%H%M%S"), "/"
                   , sep = ""
  ) 

  dir.create(out.dir)
  errors <- ErrorStats(data.name, epsilon.1, epsilon.2, out.dir)   #創資料夾
  
  for (i in 1:nrun) {                                                #第i圈跑到NRUN圈
    tag.run <- paste("-run-",i, sep = "")                         #顯示跑第幾圈
    tag.out <- paste(tag.sample, tag.run, sep = "")               #顯示資訊
    sample.filename <- paste(out.dir, data.name                        #貼出資訊
                             , '-eps1-', tag.out, '.dat', sep = "")
    if (!file.exists(sample.filename)|| flag.sample) {                    #確認資料存在
      beta <- compute_best_sampling_rate_with_Gtest(data.name, curr.data$DB.size              #β：取樣率(在util.R裡面有函式計算)
                                                    , epsilon.1
                                                    , curr.data$domain)
      data.file <- curr.data$sample_data(out.dir, rate = beta                          #輸出資訊     
                                         , out.tag = paste('-eps1-', tag.out, sep = "")
      )  
    } else {                                                          #如果資料不存在?????
      data.file <- paste(data.name, '-eps1-', tag.out, sep = "")
    }
    
    
    sample.data <- curr.data$load_sample_data(out.dir, filename = data.file)
    sample.info <- curr.data$load_sample_info(out.dir, filename = data.file)

#adding benchmark
op <- options(digits.secs=6)
print( paste("start dependency graph time:", Sys.time(), sep=" ") )


    beta <- as.numeric(sample.info$sample.rate)                                         #取樣率
    sample.depgraph <- DependenceGraph$new(sample.data                                  #相依圖(DGRAPH.R裡面)
                                           , flag.sample = TRUE
                                           , flag.noise = TRUE                             #執行dgraph.r並更新參數
                                           , beta = beta
                                           , epsilon = epsilon.1
                                           , thresh.CV = CV.thresh
                                           , flag.debug = FALSE
    )
  

#adding benchmark
op <- options(digits.secs=6)
print( paste("start junctiontree time:", Sys.time(), sep=" ") )                     
  
    types <- c('CV', 'chi2', 'CV2.noisy', 'Gtest.noisy')                               #將一堆參數轉化為向量
    jtrees <- lapply(types, function(x){                                          #lapply是把東西都套入之後回傳LIST
      jtree <- JunctionTree$new(out.dir, edges=sample.depgraph$edges[[x]]
                                , nodes = sample.depgraph$nodes
                                , data.filename = data.file
                                , type=x
                                , flag.debug = TRUE)   
      return(jtree)
    })                              
    names(jtrees) <- types                                       #?????
    plot(jtrees[['CV2.noisy']]$jtree)                              
    curr.jtree <- jtrees[['CV2.noisy']]
    type <- "CV2.noisy" 
    # if using exiting junction tree structure
    #     curr.jtree <- JunctionTree$new(flag.build = FALSE
    #                                    , jtree.file = paste("output/",data.file, "-", type, "-jtree.Rdata",sep=""))

#adding benchmark
op <- options(digits.secs=6)
print( paste("start inference with merge time:", Sys.time(), sep=" ") )               #form attribute clusters

    curr.jtree$do_inference_with_merge(                                      #呼叫MAT計算MERGE
      out.dir
      , curr.data$origin
      , curr.data$domain
      , data.filename = paste(data.file, "-", type, sep="")
      , flag.noise = TRUE
      , do.consistent = TRUE
      , epsilon.2=epsilon.2
      , flag.debug = FALSE
      , flag.matlab = TRUE
    )

    if (flag.sim){

#adding benchmark
op <- options(digits.secs=6)
print( paste("start generate data time:", Sys.time(), sep=" ") )                           #模擬出和DB一樣大的synthetic data

      data.sim <- curr.jtree$simulate(curr.data$DB.size)                              #生成SIMULATE.DAT(在JTREE.R)
      data.sim <- data.sim[, colnames(curr.data$origin)]
      data.sim.file <- paste(out.dir, data.file, "-eps2-", epsilon.2, "-sim.dat",sep="")

#adding benchmark
op <- options(digits.secs=6)
print( paste("start write data time:", Sys.time(), sep=" ") )

      write.table(data.sim                         #輸出多變數資料檔
                  , data.sim.file
                  , row.names = FALSE
                  , col.names = FALSE
                  , quote=FALSE
                  , sep = ",")
      
    }
    
#adding benchmark
op <- options(digits.secs=6)
print( paste("finish time:", Sys.time(), sep=" ") )   
 
    #random query
    print("finish")
  }
}
args <- commandArgs(trailingOnly = TRUE)                                  #?????????????
main(args)
