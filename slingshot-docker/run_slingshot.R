#!/usr/bin/env Rscript

# Load libraries
options(rgl.useNULL = TRUE)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("logging"))
suppressPackageStartupMessages(library("slingshot"))
suppressPackageStartupMessages(library("SingleCellExperiment"))

# Set up logging
# Logging level choices
C_LEVEL_CHOICES <- names(loglevels)
#initialize to info setting.
logging::basicConfig(level = "INFO")

# Expected test suite path for validation in slingshot docker container
TEST_SUITE <-
  "/usr/local/lib/R/site-library/slingshot/tests/testthat/test_slingshot.R"

# Command line arguments
pargs <- optparse::OptionParser(usage = paste("%prog [options]",
                                            "--input file"))

# input (string) path to SingleCellExperiment object
pargs <- optparse::add_option(pargs, c("--input"),
            type = "character",
            default = NULL,
            action = "store",
            dest = "input",
            help = paste("Input file for analysis. Matrix (tsv)",
                         "or R data object (rds) ",
                         "\n\t\tRequires --input-type declaration.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--input-type"),
            type = "character",
            default = NULL,
            action = "store",
            dest = "input_type",
            metavar = "input_type",
            help = paste("Specify format [ rds | matrix ] of input file",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--reduced-dim"),
            type = "character",
            default = NULL,
            action = "store",
            dest = "reduced_dim",
            metavar = "reduced_dim",
            help = paste("Identifier to use if reducedDim(data) has multiple",
                        "elements. \n\t\tIf NULL, will use first",
                        "element in reducedDim(data)",
                        "\n\t\t[Default %default]"))

pargs <- optparse::add_option(pargs, c("--cluster-labels"),
            type = "character",
            default = NULL,
            action = "store",
            dest = "cluster_labels",
            metavar = "cluster_labels",
            help = paste("For matrix-type input files, required cluster-label",
                        "file \n\t\t\tcontains a vector of length n denoting",
                        "cluster labels.\n\t\tFor SingleCellExperiment objects",
                        ", provides cluster label identifier. \n\t\tIf",
                        "SlingshotDataSet, cluster labels will be taken",
                        "from it. \n\t\tIf cluster-labels=NULL, one",
                        "cluster is assumed [Default %default]"))


pargs <- optparse::add_option(pargs, c("--start-clus"),
            type = "character",
            default = NULL,
            action = "store",
            dest = "start_clus",
            metavar = "start_clus",
            help = paste("Indicates the cluster(s) *from*",
                         "which lineages will be drawn",
                         "\n\t\t(optional), [Default %default]"))
                         
pargs <- optparse::add_option(pargs, c("--end-clus"),
            type = "character",
            default = NULL,
            action = "store",
            dest = "end_clus",
            metavar = "end_clus",
            help = paste("Indicates the cluster(s) forced to",
                         "be leaf nodes in their trees",
                         "\n\t\t(optional), [Default %default]"))

pargs <- optparse::add_option(pargs, c("--validate"),
                              type = "logical",
                              default = FALSE,
                              action = "store_true",
                              dest = "validate",
                              help = paste("Run slingshot validation,",
                                          "(optional), [Default %default]"))

# Check arguments to ensure user input meets certain requirements    
check_arguments <- function(arguments){
  # Require an input-type if input file provided
  if ("input" %in% names(arguments)) {
    if ( ! ("input_type" %in% names(args))) {
      logging::logerror(paste(":: --input-type: must declare input type."))
      stop(paste("error, --input-type: must declare input type for input file"))
    } else {
      # require cluster-labels if input is matrix file
      if (arguments$input_type == "matrix") {
        if (( ! ("cluster_labels" %in% names(arguments) )) ||
          (arguments$cluster_labels == "") ||
          (is.na(arguments$cluster_labels))) {
          error_message <- paste("error, no --cluster-labels provided.",
                                "--cluster-labels file required if input of",
                                "matrix type")
          logging::logerror(paste("::", error_message))
          stop(error_message)
        } else {
          # Make sure the clusterLabel file exists
          if ( ! file.exists(arguments$cluster_labels)){
            error_message <- paste0("cluster-labels file, ",
              arguments$cluster_labels, ", does not exist. ",
              "Please check for input files and try again")
            logging::logerror(error_message)
            stop(error_message)
          }
        }
      }
    }
    # Make sure the input data file exists if provided
    if ( ! file.exists(arguments$input)){
      error_message <- paste0("Provided input file, ",
      arguments$input, ", does not exist. ",
      "Please check your input and try again")
      logging::logerror(error_message)
      stop(error_message)
    }
    # Parse input filename
    input_file <- unlist(strsplit(args$input, "[.]"))
    # Define output name based on input name.
    # If path to file is provided, path is maintained
    # Set the default name of an output rds file
    arguments$output_file <- paste0(input_file[1:length(input_file) - 1], ".rds")
    input_folder <-  unlist(strsplit(args$input, '/'))
    arguments$pt_output_file <- paste(c(input_folder[1:length(input_folder)-1],"/SlingshotPT.csv"), collapse='/')
    arguments$pt_curve_file <- paste(c(input_folder[1:length(input_folder)-1],"/curves.csv"), collapse='/')

  }
  return(arguments)
}


#process user-submitted inputs and set defaults
args_parsed <- optparse::parse_args(pargs, positional_arguments = TRUE)
args <- args_parsed$options
logging::loginfo(paste("Checking input arguments", sep = "" ))
args <- check_arguments(args)

#read in data and run slingshot
if ( args$validate ){
    # check that test suite is in expected location
    if  ( ! file.exists(TEST_SUITE)){
    error_message <- paste0("Test suite file, ", TEST_SUITE, ", is missing. ",
    "Please confirm the correct docker container is in use")
    stop(error_message)
    } else {
      test_result <- tools::testInstalledPackage("slingshot",
                                                 types = c("tests"))
    }
    if (test_result == 0L) {
      logging::loginfo(paste("Validation successful. For details consult",
                              "output file ./slingshot-tests/testthat.Rout\n"))
    } else {
      error_message <- paste0("Validation failed, please consult file",
                              "./slingshot-tests/testthat.Rout.fail")
      logging::logerror(error_message)
      stop(error_message)
    }
} else if ( args$input_type == "rds" ) {
  sce <- readRDS(args$input)
  if ( "cluster_labels" %in% names(args) ) {
    result <- slingshot(sce, clusterLabels = args$cluster_labels,
                        reducedDim = args$reduced_dim, 
                        start.clus = args$start_clus, end.clus = args$end_clus)
  } else {
    result <- slingshot(sce, reducedDim = args$reduced_dim, 
                        start.clus = args$start_clus, end.clus = args$end_clus)
  }
} else if ( args$input_type == "plain" ) {
    coordinates <- as.matrix(read.delim(args$input, header = FALSE,
                            row.names = NULL, sep = "\t", skip = 0))
    labels <- as.vector(as.numeric(read.delim(args$cluster_labels,
                        header = FALSE, row.names = NULL, sep = "\t",
                        stringsAsFactors = TRUE, skip = 0)[, 1]))
    result <- slingshot(coordinates, labels,
                        start.clus = args$start_clus, end.clus = args$end_clus)
} else if ( args$input_type == "matrix" ) {
  coordinates <- as.matrix(read.delim(args$input, header = TRUE,
                                      row.names = 1, sep = "\t"))
  labels <- as.vector(as.numeric(read.delim(args$cluster_labels,
                                            header = TRUE, row.names = 1, sep = "\t",
                                            stringsAsFactors = TRUE)[, 1]))

  # print(as.character(strsplit(args$end_clus, ",")[[1]]))
  if (!is.null(args$end_clus)){
      result <- slingshot(coordinates, labels,
                          start.clus = args$start_clus, end.clus = as.character(strsplit(args$end_clus, ",")[[1]]))
      }else{
      result <- slingshot(coordinates, labels,
                          start.clus = args$start_clus)
      }
} else {
    error_message <- paste0("Unknown input type, ", args$input_type, ". ",
    "Please check your input-type and try again")
    logging::logerror(error_message)
    stop(error_message)
}
if ("input" %in% names(args)) {
  #save data to current working directory
  saveRDS(result, args$output_file)
  ptValues <- slingPseudotime(result)
  colnames(ptValues) <- gsub("curve", "PseudoTime",colnames(ptValues))
  write.table(ptValues, args$pt_output_file, quote = F, sep = ',')
ds <- list()
for (col in colnames(ptValues)){
idx <- !is.na(ptValues[,col])
x <- as.numeric(coordinates[idx,1])
y <- as.numeric(coordinates[idx,2])
lo <- loess(y~x)

ypred <- predict(lo,x)
ds[paste0('y',col)] <- list(as.character(ypred))
ds[paste0('x',col)] <- list(as.character(x))
}

cat(sapply(ds, toString), file = args$pt_curve_file, sep="\n")

} else {
  cat("Validation log in slingshot-tests/testthat.Rout\n")
}


