# import functions
library(shiny)
library(sendmailR)
library(shinyjs)
library(tools)
library(ggplot2)
# set path variable
if (file.exists("/home/eden/eden.sh")) {
  # we are inside the docker
  packrat::on()
  # folders to store input faa and ffn files
  faa.path <<- "/home/eden/data/faa"
  ffn.path <<- "/home/eden/data/ffn"
  sample.path <<- "/home/eden/data/samples.txt"
  hmm.path <<- "/home/eden/data/model.hmm"
  tar.path <<- "/home/eden/data/tar"
  fasta.path <<- "/home/eden/data/fasta"
  folder.path <<- "data"
  log.path <<- "/home/eden/shinylog.txt"
  lock.file <<- "/home/eden/lock.txt"
} else {
  # we are online/local hosted
  # folders to store input faa and ffn files
  log.path <<- "log.txt"
  faa.path <<- "faa"
  ffn.path <<- "ffn"
  fasta.path <<- "fasta"
  tar.path <<- "tar"
  sample.path <<- "samples.txt"
  hmm.path <<- "model.hmm"
  log.path <<- "log.txt"
  lock.file <<- "lock.txt"
}
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

if (!file.exists(fasta.path)) {
  dir.create(fasta.path)
}
if (!file.exists(faa.path)) {
  dir.create(faa.path)
}
if (!file.exists(ffn.path)) {
  dir.create(ffn.path)
}
Sys.chmod(fasta.path, mode = "0777", use_umask = TRUE)
Sys.chmod(ffn.path, mode = "0777", use_umask = TRUE)
Sys.chmod(faa.path, mode = "0777", use_umask = TRUE)


system2("echo", paste('";;server ready" >> ', log.path, sep = ""))

# define header
headerPanel_2 <- function(title, h, windowTitle = title) {
  tagList(tags$head(tags$title(windowTitle)),
          h(title))
}

shinyUI(fluidPage(
  theme = "edentheme.css",
  headerPanel_2(HTML(' '), h1, ""),
  fluidRow(
    useShinyjs(),
    column(
      8,
      
      #     wellPanel(
      tabsetPanel(
        tabPanel("File upload",
                 wellPanel(uiOutput("upload_ui_head"))
                 ,
                 value = "start"),
        
        tabPanel("File grouping",
                 wellPanel(uiOutput("upload_ui_mid")),
                 value = "start"),
        
        tabPanel("Gene family definition",
                 wellPanel(uiOutput("upload_ui_hmm")),
                 value = "start"),
        
        tabPanel("Parameters",
                 wellPanel(uiOutput("upload_ui_bottom")),
                 value = "start"),
        
        #
        # tabPanel(
        #   "start",
        #   ,
        #   value = "start"
        # ),
        #
        
        id = "tsp"
      )
      
    ),
    
    column(
      4,
      htmlOutput("statusrunning"),
      htmlOutput("statusfinished"),
      htmlOutput("statuscheck"),
      tableOutput("log"),
      htmlOutput("text1"),
      # textOutput("text2"),
      uiOutput("upload_ui_buttons")
      
    )
    
  )
))
