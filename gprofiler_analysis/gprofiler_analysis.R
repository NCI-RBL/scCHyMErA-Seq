
## Developed by Carl McIntosh
## Date: 5 March 2025
## Version 0.5.0

############################## Extract Arguements ##############################
## https://cran.r-project.org/web/packages/argparse/vignettes/argparse.html
library(argparse)

print("############################# gprofiler_analysis.R #############################")
parser <- ArgumentParser()

parser$add_argument("--output", default="./", help = "Specify output directory")
parser$add_argument("--use_bkg_known", nargs='?', type="logical", help="Run against all known.", default=FALSE)
args <- parser$parse_args()

print(paste("Ouput results directory:", args$output, sep=" "))

################################################################################
# https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
## Required Library Packages
# install.packages("gprofiler2")
# install.packages("argparse")
# install.packages("TAF")

## Install if necessary required packages.
required_packages <- c("gprofiler2", "argparse", "plotly", "TAF")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

## Load required packages.
library(gprofiler2)
library(plotly)
library(TAF)

################################### Functions ##################################
## Extract Bkg List Name - the first line
gene_list_name <- function(line) {
  x <- unlist(strsplit(toString(line), ";"))
  my_name <- unlist(strsplit(toString(x[1]), ";"))
  return(my_name)
}

## Extract Bkg List - the first line
bkg_gene_list <- function(line) {
  x <- unlist(strsplit(toString(line), ";"))
  my_list <- unlist(strsplit(toString(x[2]), ","))
  return(my_list)
}

## Extact lines 2-N as named list
named_gene_list <- function(line) {
  x <- unlist(strsplit(toString(line), ";"))
  my_list <- strsplit(toString(x[2]), ",")
  my_name <- unlist(strsplit(toString(x[1]), ";"))
  names(my_list) <- "name"
  list_length <- length(my_list$name)

  names(my_list) <- my_name
  return(my_list)
}

## Get length of gene list
named_gene_list_length <- function(line) {
  x <- unlist(strsplit(toString(line), ";"))
  my_list <- strsplit(toString(x[2]), ",")
  names(my_list) <- "name"
  list_length <- length(my_list$name)
  return(list_length)
}

## Extact lines 2-N as named list
gene_list <- function(line) {
  x <- unlist(strsplit(toString(line), ";"))
  my_list <- strsplit(toString(x[2]), ",")
  return(my_list)
}

## Open CSV File
setwd(args$output)
gene_lists <- as.list(readLines("GeneLists/Group_All_GeneList.csv"))

# Pull in background list
bkg_list_name <- gene_list_name(gene_lists[1])
bkg_list <- bkg_gene_list(gene_lists[1])

custom_bg <- NULL

if (bkg_list_name == "background") { custom_bg <- bkg_list }

mkdir("results_custom_bkg")

custom_urls <- c()
custom_url_names <- c()


for (list_idx in 2:length(gene_lists)) {
  gList_length <- named_gene_list_length(gene_lists[list_idx])
  gList <- named_gene_list(gene_lists[list_idx])
  name_list <- gene_list_name(gene_lists[list_idx])

  print("--------------------------------------------------------------------------------")
  print(paste("Group Name:      ", name_list, sep = " "))
  print(paste("Gene List Length:", gList_length, sep = " "))

  tmp_list_name <- names(gList)

  skip_to_next <- TRUE

  if(args$use_bkg_known)
  {
    mkdir("results_know_bkg")
    tryCatch(
      {

        ################################## Known background ##################################
        print("++++++ Known Option Background ++++++")
        names(gList) <- paste(tmp_list_name, "- Known Background", sep=" ")

        gostres = gost(
          query = gList,
          organism = "hsapiens",
          significant = TRUE,
          ordered_query = FALSE,
          multi_query = FALSE,
          correction_method = "g_SCS",
          user_threshold = 0.05,
          domain_scope = "known",
          numeric_ns = ""
          )

        # head(gostres$result, 3)
        graph <- gostplot(gostres, capped = TRUE, interactive = TRUE)

        # Save html
        file_name <- gsub(":", "_", name_list)
        file_name <- gsub("-", "_", file_name)

        html_file_name <- paste(file_name, "index.html", sep="_")
        htmlwidgets::saveWidget(as_widget(graph), paste("results_know_bkg", html_file_name, sep="/"))

        # Save table
        results_file_name <- paste(file_name, "results_know_bkg.csv", sep="_")
        gost_results_file_name <- paste("results_know_bkg", results_file_name, sep="/")
        gost_results_df <- apply(gostres$result,2,as.character)
        write.csv(gost_results_df, gost_results_file_name)
      },
      error = function(msg){
        print(paste("ERROR with ", name_list, sep=" "))
        print(paste("    ", msg, sep=""))
      },
      finally = {
      }
    )
  }

  tryCatch(
    {
      ################################## Custom Background ##################################
      print("++++++ Custom Background ++++++")
      names(gList) <- paste(tmp_list_name, "- Custom Background", sep=" ")
        gostres_bkg = gost(
          query = gList,
          organism = "hsapiens",
          ordered_query = FALSE,
          multi_query = FALSE,
          correction_method = "g_SCS",
          user_threshold = 0.05,
          domain_scope = "custom",
          custom_bg = custom_bg,
          numeric_ns = "",
          highlight = TRUE,
        )

      # head(gostres_bkg$result, 3)
      graph_bkg <- gostplot(gostres_bkg, capped = TRUE, interactive = TRUE)

      # Save html
      file_name <- gsub(":", "_", name_list)
      file_name <- gsub("-", "_", file_name)

      html_bkg_file_name <- paste(file_name, "bkg_custom.html", sep="_")
      htmlwidgets::saveWidget(as_widget(graph_bkg), paste("results_custom_bkg", html_bkg_file_name, sep="/"))

      # Save table
      results_bkg_file_name <- paste(file_name, "bkg_custom.csv", sep="_")
      gost_results_bkg_file_name <- paste("results_custom_bkg", results_bkg_file_name, sep="/")
      gost_results_bkg_df <- apply(gostres_bkg$result,2,as.character)
      write.csv(gost_results_bkg_df, gost_results_bkg_file_name)

      ## Get URL
      gostres_bkg_url = gost(
        query = gList,
        organism = "hsapiens",
        ordered_query = FALSE,
        multi_query = FALSE,
        correction_method = "g_SCS",
        user_threshold = 0.05,
        domain_scope = "custom",
        custom_bg = custom_bg,
        source= c("GO:BP","GO:MF","KEGG"),
        numeric_ns = "",
        highlight = TRUE,
        as_short_link = TRUE
      )
      custom_urls <- c(custom_urls, gostres_bkg_url)
      custom_url_names <- c(custom_url_names, name_list)
    },
    error = function(msg){
      print(paste("ERROR with ", name_list, sep=" "))
      print(paste("    ", msg, sep=""))
    },
    finally = {
    }
  )

  url_df <- data.frame(custom_url_names, custom_urls)
  write.csv(url_df, "custom_urls.txt")

}

