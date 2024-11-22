options(warn=-1)
graphics.off()
print("Loading libraries required...")
list.of.packages <- c("shiny","shinyWidgets","shinyjs","igraph","stringr", "dplyr")
invisible(lapply(list.of.packages, library, character.only = TRUE))

working_directory <- getwd()
library(renoir)
library(devtools)
library(ggbiplot)

learning.methods = list_supported_learning_methods()
LM <- learning.methods$id


source(paste0(working_directory,'//CCLE2.R'))
inputs_dir = paste0(working_directory,"//inputs")
setwd(inputs_dir)

mut_matrix_csv_small <- mut_matrix_csv[,sapply(mut_matrix_csv, function(x) length(unique(x))) <= 100]
molecular_features <- names(mut_matrix_csv_small)


mutations <- unique(mut_matrix_csv$VariantInfo)

genes_HGNC <- sort(unique(colnames(expr_matrix_csv)))
genes_HGNC <- unique(genes_HGNC)
genes_HGNC <- str_sort(genes_HGNC)
##########NEW 2024##############

setwd(working_directory)

models <- as.data.frame(read.csv(paste0(inputs_dir,"/Model.csv"), header=TRUE))
models_bkp <-models


CCLE_models<- unique(models$DepmapModelType)
CCLE_models<- c(CCLE_models,"ALL_MODELS")


clinical_small <- models[,sapply(models, function(x) length(unique(x))) <= 20]
clinical_features <- names(clinical_small)

PMA <- unique(models$PrimaryOrMetastasis)
PMA <- c(PMA, "ALL")

AK <- unique(models$AgeCategory)
AK <- c(AK, "ALL")
Sex <- unique(models$Sex)
Sex<- c(Sex, "ALL")
result <- ""
Lineage <- unique(models$OncotreeLineage)

clinical_features_col <-""

molecular_features_col <- ""

library(bslib)

ui <- fluidPage(

  tags$head(
    tags$style(HTML('#run{background-color:white}'))
  ),
  useShinyjs(),
  setBackgroundColor(
    color = c("#F7FBFF", "#2171B5"),
    gradient = "linear",
    direction = "top"
  ),
  titlePanel("HARMONIZER"),
  sidebarLayout(
    sidebarPanel(
      actionButton("run_ccle", "RUN EXPERIMENT", icon("play"),
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),

      ########################################################################################################
      shinyjs::hidden(p(id = "running_text", "Processing...")),
      selectInput("violin_gene", "Select gene of interest (usually a TF):", choice = NULL, selected = "TP53..7157.", selectize = TRUE),
      selectInput('selectfile','Select custom celllines from a file ("inputs>>CELL_LINES"):',choice = tools::file_path_sans_ext(list.files('./inputs/CELL_LINES/')),selected = "TNBC"),
      checkboxInput("MODELS_CCLE", label = "Choose from predefined models in CCLE", FALSE),
      selectInput("OncotreeLineage", "Select lineage:", choice = Lineage, selected = "Breast", selectize = TRUE),
      selectInput('selectfile2','Select DepMap Model Type:',choice = sort(CCLE_models),selected = "BRCA"),
      selectInput("primary_metastasis", "Select either Primary, Metastatic, or all cell lines:", choice = PMA, selected = "ALL", selectize = TRUE),
      selectInput("age", "Select age group:", choice = AK, selected = "ALL", selectize = TRUE),
      selectInput("sex", "Select Sex:", choice = Sex, selected = "ALL", selectize = TRUE),

      radioButtons(inputId="choice", label="Cell line dataset selection:",
                   choices=list("Run only for the selected disease" = 1, "Run for all cancer cell line samples in single run" = 2,
                                "Run per cancer types separately" = 3),selected = 1),
      tags$hr(style="border-color: blue;"),
      #######################################################################################################


    ),
    mainPanel(

      tags$hr(style="border-color: blue;"),

      ######################################################################################################
      checkboxInput("CARNIVAL_flag", label = "Run CARNIVAL optimization", FALSE),


      checkboxInput("molecular_flag", label = "Select from molecular features", FALSE),
      selectInput("key_opt_molecular", "Enter the network comparison molecular feature:", choice = NULL, selected = NULL, selectize = TRUE),

      checkboxInput("clinical_flag", label = "Select from clinical features", FALSE),
      selectInput("key_opt_clinical", "Enter the network comparison clinical feature:", choice = NULL, selected = NULL, selectize = TRUE),

      selectInput("feature", "Select resulting feature :", choice = NULL, selected = result[1], selectize = TRUE),
      selectInput("comparison_type", "Select comparison type :", choice = c("Equal","VS all other"), selected = "Equal", selectize = TRUE),

      textInput("top", "Number of top measurements from Dorothea to use (0 = all):", 0),
      checkboxInput("FEM_flag", label = "Use whole expression matrix?", TRUE),
      tags$hr(style="border-color: blue;"),
      ######################################################################################################

      checkboxInput("GLM_flag", label = "Perform Classification", FALSE),
      checkboxGroupInput("LM", "Select learning method:",LM, selected = "elasticnet"),
      radioButtons(inputId="GLM", label="GLM options:",
                   choices=list("Binomial" = 1,"Multinomial" = 2),selected = 1),
      textInput("folds", "Enter the number of Folds :", 5),
      textInput("hotspots_user", "User defined GLM classification delimiter:", "p.R175H|p.R248Q|p.R273H|p.R248W|p.R273C|p.R282W|p.G245S"),
      radioButtons(inputId="binomial_key", label="Classify across:",
                   choices=list("WT/MT" = 1, "Deleterious" = 2,
                                "User defined (change above)" = 3),selected = 4),
      checkboxInput("mutation_type_flag", label = "or Select mutation type", TRUE),
      selectInput("mutation_type", "Select mutation to classify on:", choice = mutations , selected = "missense_variant", selectize = TRUE),
      checkboxInput("GLM_all_u", label = "Use all genes for GLM", FALSE),
      tags$hr(style="border-color: blue;"),
      textInput("dataset_version_u", "CCLE_version:", "24Q2"),
      checkboxInput("new_version_u", label = "Read new CCLE version?", FALSE),

      img(src='logo.png', align = "center")


    )

  )
)

######################################################################################

server <- function(input, output,session) {

  updateSelectizeInput(session, 'violin_gene', choices = genes_HGNC,selected = "TP53..7157.", server = TRUE)
  updateSelectizeInput(session, 'key_opt_molecular', choices = molecular_features, selected = molecular_features[1], server = TRUE)
  updateSelectizeInput(session, 'key_opt_clinical', choices = clinical_features, selected = clinical_features[1], server = TRUE)
  updateSelectizeInput(session, 'feature', choices = result, selected = "PrimaryOrMetastasis", server = TRUE)


  observeEvent(input$key_opt_molecular, {
    col <- which( colnames(mut_matrix_csv)== input$key_opt_molecular )
    result <- unique(mut_matrix_csv[,col])
    result <- result[result != ""]
    if (input$molecular_flag == TRUE) {
      updateSelectInput(session, "feature",
                        choices = result[!is.na(result)])
    }
  })

  observeEvent(input$molecular_flag, {
    col <- which( colnames(mut_matrix_csv)== input$key_opt_molecular )
    result <- unique(mut_matrix_csv[,col])
    result <- result[result != ""]
    if (input$molecular_flag == TRUE) {
      updateSelectInput(session, "feature",
                        choices = result[!is.na(result)])
    }
  })

  observeEvent(input$key_opt_clinical, {
    col <- which( colnames(models_bkp)== input$key_opt_clinical )
    result <- unique(models_bkp[,col])
    result <- result[result != ""]
    if (input$clinical_flag == TRUE) {
      updateSelectInput(session, "feature",
                        choices = result[!is.na(result)])
    }
  })

  observeEvent(input$clinical_flag, {
    col <- which( colnames(models_bkp)== input$key_opt_clinical )
    result <- unique(models_bkp[,col])
    result <- result[result != ""]
    if (input$clinical_flag == TRUE) {
      updateSelectInput(session, "feature",
                        choices = result[!is.na(result)])
    }
  })


  observeEvent(input$OncotreeLineage, {
    models<-models_bkp
    models <- models %>% dplyr::filter(models$OncotreeLineage == input$OncotreeLineage)
    CCLE_models <- unique(models$DepmapModelType)
    CCLE_models<- c("ALL_MODELS",CCLE_models)
    CCLE_models<-sort(CCLE_models)
    updateSelectInput(session, "selectfile2",
                      choices = CCLE_models
    )
  })

  observe({

    shinyjs::toggleState("MODELS_CCLE", input$choice != 2)
    shinyjs::toggleState("selectfile", input$choice != 2)


    shinyjs::toggleState("GLM", input$GLM_flag == "TRUE")
    shinyjs::toggleState("selectfile", input$MODELS_CCLE == "FALSE" &  input$choice == 1)
    shinyjs::toggleState("OncotreeLineage", input$MODELS_CCLE == "TRUE")
    shinyjs::toggleState("selectfile2", input$MODELS_CCLE == "TRUE")
    shinyjs::toggleState("primary_metastasis", input$MODELS_CCLE == "TRUE")
    shinyjs::toggleState("age", input$MODELS_CCLE == "TRUE")
    shinyjs::toggleState("sex", input$MODELS_CCLE == "TRUE")
    shinyjs::toggleState("OncotreeLineage", (input$choice == 1 | input$choice == 3) && input$MODELS_CCLE == "TRUE")
    shinyjs::toggleState("selectfile2", (input$choice == 1 | input$choice == 3) && input$MODELS_CCLE == "TRUE")

    shinyjs::toggleState("GLM_all_u", input$GLM_flag == "TRUE" & input$GLM == 1)
    shinyjs::toggleState("LM", input$GLM_flag == "TRUE" & input$GLM == 1)
    shinyjs::toggleState("folds", input$GLM_flag == "TRUE" & input$GLM == 1)
    shinyjs::toggleState("binomial_key", input$GLM_flag == "TRUE" & input$GLM == 1)
    shinyjs::toggleState("mutation_type", input$GLM_flag == "TRUE" & input$GLM == 1)
    shinyjs::toggleState("mutation_type_flag", input$GLM_flag == "TRUE" & input$GLM == 1)
    shinyjs::toggleState("hotspots_user", input$GLM_flag == "TRUE" & input$GLM == 1 & input$mutation_type_flag == "FALSE")

    shinyjs::toggleState("FEM_flag", input$CARNIVAL_flag == "TRUE")


    shinyjs::toggleState("clinical_flag", input$CARNIVAL_flag == "FALSE")
    shinyjs::toggleState("molecular_flag", input$CARNIVAL_flag == "FALSE")

    shinyjs::toggleState("key_opt_clinical", input$CARNIVAL_flag == "TRUE")
    shinyjs::toggleState("key_opt_molecular", input$CARNIVAL_flag == "TRUE")

    shinyjs::toggleState("key_opt_molecular", input$molecular_flag == "TRUE")
    shinyjs::toggleState("key_opt_clinical", input$clinical_flag == "TRUE")


    shinyjs::toggleState("molecular_flag", input$clinical_flag == "FALSE" & input$CARNIVAL_flag == "TRUE")
    shinyjs::toggleState("clinical_flag", input$molecular_flag == "FALSE" & input$CARNIVAL_flag == "TRUE")



    shinyjs::toggleState("top", input$CARNIVAL_flag == "TRUE")
    shinyjs::toggleState("comparison_type", input$CARNIVAL_flag == "TRUE")

    shinyjs::toggleState("feature", input$CARNIVAL_flag == "TRUE")

    shinyjs::toggleState("binomial_key", input$mutation_type_flag == "FALSE")
    shinyjs::toggleState("mutation_type", input$mutation_type_flag == "TRUE" & input$GLM_flag == "TRUE" & input$GLM == 1)


  })
  v <- reactiveValues()
  plotReady <- reactiveValues(ok = FALSE)


  observeEvent(input$run_ccle, {

    shinyjs::disable("run_ccle")
    shinyjs::show("running_text")

    top_score <- input$score_user
    GLM_user_word <- NULL

    if ((input$GLM_flag == TRUE) & (input$mutation_type_flag == TRUE)){
      if (input$binomial_key == 3)
      {
        GLM_user_word <- input$hotspots_user
      }
    }

    hotspots_default <- c("p.R175H","p.R248Q","p.R273H","p.R248W", "p.R273C", "p.R282W", "p.G245S")

    condition <- input$conditions

    dataset_version <- input$dataset_version_u

    plotReady$ok <- FALSE
    c_flag = FALSE
    GLM_all = FALSE
    new_version = FALSE
    FEM_user = FALSE
    load_other_GLM = FALSE

    if (input$GLM_all_u == TRUE){GLM_all = TRUE}
    if (input$CARNIVAL_flag == TRUE) {c_flag = TRUE}
    if (input$new_version_u == TRUE) {new_version=TRUE}
    if (input$FEM_flag == TRUE) {FEM_user <- TRUE}
    GLM_flag <- input$GLM_flag

    if (input$GLM_flag == TRUE) {
      if (input$GLM == 1) {reg_type = "binomial"}
      else if (input$GLM == 2) {reg_type = "multinomial"}
    }

    GLM_user_word <- stringr::str_replace_all(GLM_user_word, fixed(" "), "")


    PMA_user <- input$primary_metastasis
    Age_user <- input$age
    Sex_user <- input$sex
    OncotreeLineage_user <- input$OncotreeLineage

    folds <- input$folds


    key_opt <- input$feature
    key_opt_type <- input$comparison_type

    if (input$mutation_type_flag == FALSE) {
      if (input$binomial_key == 1) {key = "WT"}
      else if (input$binomial_key == 2) {key = "NOTisDeleterious"}
      else if (input$binomial_key == 3) {key = input$hotspots_user}
    }
    else {
      key <- input$mutation_type

    }

    if(input$selectfile2 == "ALL_MODELS" & input$choice != 2 ) {cell_lines_model <- OncotreeLineage_user}


    if (input$choice == 1 )
    {
      if (input$MODELS_CCLE == "FALSE")
      {
        cell_lines_model <- tools::file_path_sans_ext(input$selectfile)
      }
      else{
        cell_lines_model <- input$selectfile2
      }
    }


    if (input$choice == 2)
    {
      if (input$MODELS_CCLE == "FALSE")
      {
        cell_lines_model <- "all_cell_lines"
      }
      else{
        cell_lines_model <- "all_cell_lines"
      }
    }


    if (input$choice == 3)
    {
      if (input$MODELS_CCLE == "FALSE")
      {
        cell_lines_model <- "all_custom"
      }
      else{
        cell_lines_model <- "all_CCLE"
      }
    }

    LM <- input$LM

    if (input$molecular_flag == TRUE) {
      molecular_features_col <- input$key_opt_molecular
    }

    else if(input$clinical_flag == TRUE){
      clinical_features_col <- input$key_opt_clinical
    }

    v$res <- CCLE2(cell_lines_model,
                   PMA_user,
                   Age_user,
                   Sex_user,
                   OncotreeLineage_user,
                   input$violin_gene,
                   GLM_flag,
                   reg_type,
                   c_flag,
                   GLM_all,
                   new_version,
                   dataset_version,
                   GLM_user_word,
                   key,
                   load_other_GLM,
                   FEM_user,
                   LM,
                   folds,
                   key_opt,
                   key_opt_type,molecular_features_col,clinical_features_col
    )

    setwd(working_directory)
    stopApp(returnValue = invisible())


  })

  Sys.sleep(2)
  plotReady$ok <- TRUE

}

######################################################################################

shinyApp(ui = ui, server = server)
