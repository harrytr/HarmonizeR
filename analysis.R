MNR <-function(regression_data,perm,disease_filename, reg_type, resDir, screening_method, key, LM,folds, GLM_user_word)  {
  # install.packages(c("glmnet","Matrix","parallel","doParallel","foreach","stats","utils","matrixTests","graphics"))
  #  install.packages(c("testthat","knitr","rmarkdown","plotly","htmltools","kableExtra","DT"))


  library(plotly)
  library(htmltools)
  folds <-as.integer(folds)
  set.seed(
    #A seed
    seed = 5381L,                   #a randomly chosen integer value
    #The kind of RNG to use
    kind = "Mersenne-Twister",      #we make explicit the current R default value
    #The kind of Normal generation
    normal.kind = "Inversion"       #we make explicit the current R default value
  )
  #setwd(resDir)
  print(paste0(reg_type," classification starting learning with ", LM, " and classifying on: ", key))
  #get hyperparameters
  get_hp = function(id, y){

    #Generalised Linear Models with Penalisation
    lambda = 10^seq(3, -2, length=100)
    # alpha = seq(0.1, 0.9, length = 9)
    alpha = seq(0.1, 0.9, length = 5)
    gamma = c(0, 0.25, 0.5, 0.75, 1)

    #Random Forest
    ntree = c(10, 50, 100, 250, 500)

    #Generalised Boosted Regression Modelling
    eta = c(0.3, 0.1, 0.01, 0.001)

    #Support Vector Machines
    cost      = 2^seq(from = -5, to = 15, length.out = 5)
    svm.gamma = 2^seq(from = -15, to = 3, length.out = 4)
    degree    = seq(from = 1, to = 3, length.out = 3)
    #Note that for classification nu must be
    #nu * length(y)/2 <= min(table(y)). So instead of
    #fixing it as
    # nu        = seq(from = 0.1, to = 0.6, length.out = 5)
    #we do
    nu.to = floor((min(table(y)) * 2/length(y)) * 10) / 10
    nu = seq(from = 0.1, to = nu.to, length.out = 5)

    #kNN
    k = seq(from = 1, to = 9, length.out = 5)

    #Nearest Shrunken Centroid
    threshold = seq(0, 2, length.out = 30)

    #hyperparameters
    out = switch(
      id,
      'lasso'              = list(lambda = lambda),
      'ridge'              = list(lambda = lambda),
      'elasticnet'         = list(lambda = lambda, alpha = alpha),
      'relaxed_lasso'      = list(lambda = lambda, gamma = gamma),
      'relaxed_ridge'      = list(lambda = lambda, gamma = gamma),
      'relaxed_elasticnet' = list(lambda = lambda, gamma = gamma, alpha = alpha),
      'randomForest'       = list(ntree = ntree),
      'gbm'                = list(eta = eta, ntree = ntree),
      'linear_SVM'         = list(cost = cost),
      'polynomial_SVM'     = list(cost = cost, gamma = svm.gamma, degree = degree),
      'radial_SVM'         = list(cost = cost, gamma = svm.gamma),
      'sigmoid_SVM'        = list(cost = cost, gamma = svm.gamma),
      'linear_NuSVM'       = list(nu = nu),
      'polynomial_NuSVM'   = list(nu = nu, gamma = svm.gamma, degree = degree),
      'radial_NuSVM'       = list(nu = nu, gamma = svm.gamma),
      'sigmoid_NuSVM'      = list(nu = nu, gamma = svm.gamma),
      'gknn'               = list(k = k),
      'nsc'                = list(threshold = threshold)
    )

    return(out)
  }
  #Set some directories
  #Input file must be in the directory defined here
  #for the following code to be run without changes
  print("Setting up the data")
  data = readRDS(regression_data)

  colnames(data) <- sapply(strsplit(colnames(data), split='..', fixed=TRUE),function(x) (x[1]))

  analysis.type = reg_type

  x = as.matrix(x = data[,-(1),drop=F])

  rownames(x) = data$Variant_Classification


  if(identical(analysis.type, "binomial")){
    print("Binomial selected")
    y = data[,1]

    resp.type = analysis.type
    if (!is.null(GLM_user_word)) {y = ifelse(stringr::str_detect(y,GLM_user_word),"1","0")}
    else {y = ifelse(stringr::str_detect(y,key),"1","0")}


    names(y) = data$Variant_Classification

    x = x[names(y),,drop=F]

  } else if(identical(analysis.type, "multinomial"))

    {
    print("multinomial selected")
    resp.type = "multinomial"

    dataXX <- data
    write.csv(dataXX[,1], "dataYY.csv")

    dataXX$Variant_Classification <- str_sub(dataXX$Variant_Classification, start = 1, end=2)
    dataXX <- dataXX %>% group_by(Variant_Classification) %>% mutate(Variant_Classification = cur_group_id())

    dataXX <- dataXX %>% group_by(Variant_Classification) %>% filter(n() > folds)

    write.csv(dataXX,"dataXX.csv")
    dataY <- dataXX$Variant_Classification

    names(dataY) = dataXX$Variant_Classification

    x = x[names(y),,drop=F]
    #invisible(readline(prompt="Press [enter] to continue"))




  }

  #as factor
  y = as.factor(y)
  write.csv(x,"X_GLM.csv")
  write.csv(y,"Y_GLM.csv")
  print("Setting parameters...")
  #-------------------------------------------------------------------------------#
  # NEW SET UP
  #-------------------------------------------------------------------------------#
  #learning method
 # learning.methods = list_supported_learning_methods(x = "classification")
  learning.methods.ids = LM

  #metric for tuning
  performance.metric.id.tuning = "acc"

  #metrics for evaluation
  performance.metric.ids.evaluation = c("acc", "precision","class")

  #sampling for tuning
  sampling.method.id.tuning = "cv"

  #sampling for evaluation
  sampling.method.id.evaluation = "random"

  resp.type = analysis.type

  #container
  learners = list()

  tuner = Tuner(
    id = "grid.search",
    sampler = Sampler(
      method = sampling.method.id.tuning,
      k = folds,
      n = integer(),
      strata = y
    ),
    looper   = Looper(cores = 1L),
    logger   = Logger(verbose = T, level = "INFO")
  )

  #loop
  for(learning.method.id in learning.methods.ids){
    #manual setup
    learners[[learning.method.id]] = Learner(
      tuner      = tuner,
      trainer    = Trainer(id = learning.method.id),
      forecaster = Forecaster(id = learning.method.id),
      scorer     = ScorerList(Scorer(id = performance.metric.id.tuning)),
      selector   = Selector(id = learning.method.id),
      recorder   = Recorder(id = learning.method.id),
      marker     = Marker(id = learning.method.id),
      logger     = Logger(level = "ALL")
    )
  }

  #Evaluator
  evaluator = Evaluator(
    #Sampling strategy: stratified random sampling without replacement
    sampler = Sampler(
      method = "random",
      k = folds,
      strata = y,
      N = as.integer(length(y))
    ),

    #Performance metric
    scorer  = ScorerList(
      Scorer(id = performance.metric.ids.evaluation[1]),
      Scorer(id = performance.metric.ids.evaluation[2]),
      Scorer(id = performance.metric.ids.evaluation[3])
    )
  )


  #define path
  outdir = file.path(resDir, "data-raw", "benchmark", "classification", disease_filename, "analysis")

  #create if not existing
  if(!dir.exists(outdir)){dir.create(path = outdir, showWarnings = F, recursive = T)}

  #container list
  resl = list()

  #loop
  for(learning.method.id in learning.methods.ids){

    #Each analysis can take hours, so we save data
    #for future faster load

    #path to file
    fp.obj = file.path(outdir, paste0(learning.method.id,".rds"))
    fp.sum = file.path(outdir, paste0("st_",learning.method.id,".rds"))

    #check if exists
    if(file.exists(fp.sum)){
      #load
      cat(paste0("Reading ", learning.method.id, "..."), sep = "")
      resl[[learning.method.id]] = readRDS(file = fp.sum)
      cat("DONE", sep = "\n")
    } else {

      cat(paste("Learning method:", learning.method.id), sep = "\n")

      #Set a seed for RNG
      set.seed(
        #A seed
        seed = 5381L,                   #a randomly chosen integer value
        #The kind of RNG to use
        kind = "Mersenne-Twister",      #we make explicit the current R default value
        #The kind of Normal generation
        normal.kind = "Inversion"       #we make explicit the current R default value
      )

      resl[[learning.method.id]] = renoir(
        # filter,

        #Training set size
        npoints = 3,
        # ngrid,
        nmin = round(nrow(x)/2),

        #Loop
        looper = Looper(),

        #Store
        filename = "renoir",
        outdir   = NULL,
        restore  = TRUE,

        #Learn
        learner   = learners[[learning.method.id]],

        #Evaluate
        evaluator = evaluator,

        #Log
        logger    = Logger(level = "ALL", verbose = T),

        #Data for training
        hyperparameters = get_hp(id = learning.method.id, y = y),
        x         = x,
        y         = y,
        weights   = NULL,
        offset    = NULL,
        resp.type = resp.type,

        #Free space
        rm.call = FALSE,
        rm.fit  = FALSE,

        #Group results
        grouping = TRUE,

        #No screening
        screening = NULL,

        #Remove call from trainer to reduce space
        keep.call = F
      )

      #save obj
      saveRDS(object = resl[[learning.method.id]], file = fp.obj)

      #create summary table
      resl[[learning.method.id]] = renoir:::summary_table.RenoirList(resl[[learning.method.id]], key = c("id", "config"))

      #save summary table
      saveRDS(object = resl[[learning.method.id]], file = fp.sum)

      cat("\n\n", sep = "\n")
    }
  }

  #create summary table
  resl = do.call(what = rbind, args = c(resl, make.row.names = F, stringsAsFactors = F))
  #plot
  g1<- renoir:::plot.RenoirSummaryTable(
    x = resl[resl$config == "opt",,drop=F],
    measure     = "precision",
    set         = "train",
    interactive = T,
    add.boxplot = F,
    add.scores  = F,
    add.best    = F,
    key         = c("id", "config")
  )
  print(g1)


  #plot
  g3<-renoir:::plot.RenoirSummaryTable(
    x = resl[resl$config == "opt",,drop=F], #select opt config
    measure     = "precision",
    set         = "test",
    interactive = T,
    add.boxplot = F,
    add.scores  = F,
    add.best    = F,
    key         = c("id", "config")
  )
  print(g3)


  #plot
  g5<-renoir:::plot.RenoirSummaryTable(
    x = resl[resl$config == "opt",,drop=F], #select opt config
    measure     = "precision",
    set         = "full",
    interactive = T,
    add.boxplot = F,
    add.scores  = F,
    add.best    = F,
    key         = c("id", "config")
  )
  print(g5)


  #plot
  g7 <- renoir:::plot.RenoirSummaryTable(
    x = resl[resl$config == "opt",,drop=F], #select opt config
    measure     = "acc",
    set         = "train",
    interactive = T,
    add.boxplot = F,
    add.scores  = F,
    add.best    = F,
    key         = c("id", "config")
  )
  print(g7)


  #plot
  g9<-renoir:::plot.RenoirSummaryTable(
    x = resl[resl$config == "opt",,drop=F], #select opt config
    measure     = "acc",
    set         = "test",
    interactive = T,
    add.boxplot = F,
    add.scores  = F,
    add.best    = F,
    key         = c("id", "config")
  )
  #plot
  print(g9)

  g10<- renoir:::plot.RenoirSummaryTable(
    x = resl[resl$config == "opt",,drop=F], #select opt config
    measure     = "acc",
    set         = "full",
    interactive = T,
    add.boxplot = F,
    add.scores  = F,
    add.best    = F,
    key         = c("id", "config")
  )
  print(g10)



} #end
