# ==============================================================================================
# global options
# ==============================================================================================
options(shiny.maxRequestSize = 2024 * 1024 ^ 2) # set max file size for file upload
options(shiny.deprecation.messages = FALSE)
options(shiny.sanitize.errors = FALSE)

shinyServer(function(input, output, session) {
  # ==============================================================================================
  # initialize status
  # ==============================================================================================
  status <-
    reactiveValues() # save the run status as a reactive object
  status$orf_finished <- FALSE # TRUE after start_orf.sh is executed
  status$eden_finished <-
    FALSE # TRUE if a finished signal detected in the log file
  status$eden_failed <-
    FALSE # is TRUE if a error signal detected in the log file
  status$samples_provided <-
    FALSE # TRUE if a samples.txt is provided
  status$faa <- FALSE
  status$ffn <- FALSE
  status$fasta <- FALSE
  status$files <- NULL
  # ==============================================================================================
  # check procedure for the detection which files are provided by the user
  # ==============================================================================================
  
  ### running status
  # function to evaluate if a process is running in the background
  update_running_status <- function() {
    my_running_status <<-
      get_running_status() # get_new_status returns the html-alert box visible on the first page
  }
  output$statusrunning = renderText({
    invalidateLater(millis = 1000, session) # update every 1 seconds
    update_running_status() # get a new status line
  })
  
  ### finished status
  # function to evaluate if a process is running in the background
  update_finished_status <- function() {
    my_finished_status <<-
      get_finished_status() # get_new_status returns the html-alert box visible on the first page
  }
  output$statusfinished = renderText({
    invalidateLater(millis = 1000, session) # update every 1 seconds
    update_finished_status() # get a new status line
  })
  
  
  # ==============================================================================================
  # procedure to update log file / progress message
  # ==============================================================================================
  # loads the log file and extract the informations for the check process bar
  get_new_log <- function() {
    print(" get_new_log")
    msg <- ""
    data <- read.table(log.path, header = F, sep = ";")
    colnames(data) <- c("step", "steps", "message")
    
    last_event <- data[nrow(data), ]$message # current message
    steps <- data[nrow(data), ]$steps # current step
    step <-
      data[nrow(data), ]$step # total number of steps till finished
    # print(last_event) # for debugging, print the last event
    if (last_event == "error") {
      # error signal
      status$check_failed <<- TRUE
      msg <-
        paste("Last log message: <span class='label label-danger'>check success</span>")
    } else if (last_event == "finished") {
      # finished signal for eden
      msg <-
        paste("Last log message: <span class='label label-success'>eden finished</span>")
      status$eden_finished <- TRUE
    } else {
      msg <-
        paste("Last log message: <span class='label label-success'>",
              last_event,
              "</span>")
    }
    return(msg)
  }
  
  # ==============================================================================================
  # user interface
  # ==============================================================================================
  
  
  ### ui head part
  output$upload_ui_head <- renderUI({
    conditionalPanel(
      condition = "input.tsp=='start'",
      helpText(
        "Step 1: specify input file format. You can either upload .faa and .ffn files of open reading frames (ORF) or the nucleotide .fasta file (in this case the ORFs will be predicted inside the pipeline)."
      ),
      
      fluidRow(column(
        5,
        fileInput(
          'files_faa',
          'upload .faa files',
          accept = c('.faa'),
          multiple = TRUE
        )
      ), column(
        5,
        fileInput(
          'files_ffn',
          'upload .ffn files',
          accept = c('.ffn'),
          multiple = TRUE
        )
      )),
      actionButton('upload_orf', label = "upload ORF"),
      hr(),
      helpText("OR you can provide .fasta files"),
      
      fileInput(
        'files_fasta',
        'select one or multiple fasta files',
        accept = c('.fasta'),
        multiple = TRUE
      ),
      actionButton('upload_fasta', label = "upload fasta")
    )
  })
  
  ### ui mid part
  output$upload_ui_mid <- renderUI({
    conditionalPanel(
      condition = "input.tsp=='start'",
      helpText(
        "Step 2: specify file handling: If you want to perform a comparative analysis you have to specify which samples are get pooled together. On default all samples will be pooled together."
      ),
      uiOutput('uploadsamplestxt')
    )
  })
  
  output$upload_ui_bottom <- renderUI({
    conditionalPanel(
      condition = "input.tsp=='start'",
      helpText("Step 3: Specify name and thresholds"),
      textInput("eden_run_name", label = "give your analysis a unique name", value = "eden_run_1"),
      #  textInput("eden_run_cpus", label = "number of CPUs used for analysis", value = "4"),
      numericInput("eden_run_cpus", label = "number of CPUs", value = 4),
      helpText(
        "Low-confidence postions can be automatically filtered for the dN/dS analysis. On deault, positions with more than 80% gaps in the alignment will not be used for dN/dS caluclation."
      ),
      sliderInput(
        "eden_run_gap",
        label = "gap proportion to filter out low confidence positions",
        min = 0,
        max = 100,
        value = 80
      )
    )
  })
  
  output$upload_ui_hmm <- renderUI({
    conditionalPanel(
      condition = "input.tsp=='start'",
      helpText(
        "Step 4: select group definition. You can select precalculated hidden markov models (HMM) or upload a .HMM file which may contain multiple hmm models for the gene families of interest."
      ),
      
      radioButtons(
        "radio",
        label = "Group definition",
        choices = list(
          "uploaded HMM model" = 1,
          "TIGRFAM models" = 2
        ),
        selected = 2
      ),
      
      fileInput(
        'hmmfile',
        'upload a hidden markov model file',
        accept = c('.HMM'),
        multiple = FALSE
      ),
      actionButton('upload_hmm', label = "upload HMM file")
    )
    # uiOutput("hmmui"))
  })
  
  output$upload_ui_buttons <- renderUI({
    conditionalPanel(
      condition = "input.tsp=='start",
      actionButton('checkButton', label = "start analysis"),
      actionButton('deletefiles', label = "reset files")
    )
  })
  
  
  # show upload section for samples.txt
  output$uploadsamplestxt <- renderUI({
    conditionalPanel(
      condition = "input.tsp=='start'",
      helpText("Please define which samples are pooled together:"),
      numericInput("groupnum", label = "number of groups", value = 2),
      uiOutput("groupboxes"),
      uiOutput("updategroups")
    )
  })
  
  output$groupboxes <- renderUI({
    members <<- as.integer(input$groupnum) # default 2
    max_pred <- as.integer(20)
    # faa mode
    #   if (input$uploadtype == "orf"){
    
    lapply(1:members, function(i) {
      fluidRow(column(
        5,
        textInput(
          inputId = paste0("name", i),
          label = paste("name", i),
          width = "100%",
          value = paste("group_", i, sep = "")
        )
      ), column(
        5,
        selectInput(
          inputId = paste0("ind", i),
          label = paste("group", i),
          choices = status$files,
          multiple = TRUE,
          width = "100%"
        )
      ))
    })
    
  })
  
  output$updategroups <- renderUI({
    actionButton('updategrouping', label = "generate grouping file")
  })
  
  # ==============================================================================================
  # button press events
  # ==============================================================================================
  
  # delete files on button press
  observeEvent(input$deletefiles, {
    isolate({
      unlink(faa.path, recursive = T, force = T)
      unlink(ffn.path, recursive = T, force = T)
      unlink(fasta.path, recursive = T, force = T)
      unlink(hmm.path, recursive = T, force = T)
      unlink(sample.path, recursive = T, force = T)
      unlink("/home/eden/data/groups.txt",
             recursive = T,
             force = T)
      unlink(log.path, recursive = T, force = T)
      system2("echo", paste('";;server ready" >> ', log.path, sep = ""))
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
      
      status <-
        reactiveValues() # save the run status as a reactive object
      input <-
        reactiveValues() # save the run status as a reactive object
      
      #
      # input <- NULL
      # js$reset()
    })
  })
  
  
  # faa upload
  observeEvent(input$upload_orf, {
    print("upload_faa using observedEvent")
    infiles_faa <- as.data.frame(input$files_faa)
    infiles_faa$dest <-
      paste(faa.path, infiles_faa$name, sep = "/")
    for (i in 1:nrow(infiles_faa)) {
      cmd <-
        paste("mv ", infiles_faa$datapath[i], " ", infiles_faa$dest[i], sep = "")
      err <- system(cmd,  intern = TRUE)
    }
    system2("echo", paste('";;faa files added" >> ', log.path, sep = ""))
    status$faa <- TRUE
    
    # process ffn
    infiles_ffn <- as.data.frame(input$files_ffn)
    infiles_ffn$dest <-
      paste(ffn.path, infiles_ffn$name, sep = "/")
    for (i in 1:nrow(infiles_ffn)) {
      cmd <-
        paste("mv ", infiles_ffn$datapath[i], " ", infiles_ffn$dest[i], sep = "")
      err <- system(cmd,  intern = TRUE)
    }
    system2("echo", paste('";;ffn files added" >> ', log.path, sep = ""))
    status$ffn <- TRUE
    
    # update files
    status$files <- list.files(path =  faa.path,
                               full.names = FALSE,
                               recursive = FALSE)
    
    
    
  })
  
  # fasta upload
  observeEvent(input$upload_fasta, {
    print("upload_fasta using observedEvent")
    
    infiles_fasta <- as.data.frame(input$files_fasta)
    infiles_fasta$dest <-
      paste(fasta.path, infiles_fasta$name, sep = "/")
    for (i in 1:nrow(infiles_fasta)) {
      cmd <-
        paste("mv ",
              infiles_fasta$datapath[i],
              " ",
              infiles_fasta$dest[i],
              sep = "")
      err <- system(cmd,  intern = TRUE)
    }
    out <- paste(err)
    system2("echo", paste('";;fasta files added" >> ', log.path, sep = ""))
    status$fasta <- TRUE
    
    status$files <- list.files(path =  fasta.path,
                               full.names = FALSE,
                               recursive = FALSE)
    
    
  })
  
  
  
  # hmm upload
  observeEvent(input$upload_hmm, {
    print("upload_hmm using observedEvent")
    # process ffn
    hmmfile <- as.data.frame(input$hmmfile)
    hmmfile$dest <- hmm.path
    cmd <-
      paste("mv ", hmmfile$datapath, " ", hmmfile$dest, sep = "")
    err <- system(cmd,  intern = TRUE)
    out <- paste(err)
    system2("echo", paste('";;hmmfile added" >> ', log.path, sep = ""))
    hmmfile_success <<- TRUE
    
  })
  
  
  # iterate over input$name_ and input$ind_ and create a samples table
  observeEvent(input$updategrouping, {
    sam <-  rep("unknown", members)
    nam <- rep("unknown", members)
    for (i in 1:members) {
      temp <- paste0("name", i)
      temp2 <<- paste0("ind", i)
      nam[i] <- input[[temp]]
      print(paste(basename(file_path_sans_ext(input[[temp2]])),  collapse =
                    "+"))
      #  sam[i] <- paste(substr(input[[temp2]], 1, nchar(input[[temp2]])-4), collapse="+")
      
      #  sam[i] <- paste(substr(input[[temp2]], 1, nchar(input[[temp2]])-4), collapse="+")
      sam[i] <-
        paste(basename(file_path_sans_ext(input[[temp2]])),  collapse = "+")
    }
    dat <- data.frame(name = nam, sample = sam)
    write.table(
      dat,
      file = sample.path,
      row.names = F,
      quote = F,
      col.names = F,
      sep = ";"
    )
    status$samples_provided <- TRUE
    system2("echo",
            paste('";;grouping file updated" >> ', log.path, sep = ""))
    
  })
  
  # ==============================================================================================
  # function that updates message based on files prsovided
  # ==============================================================================================
  
  get_running_status <- function() {
    if (file.exists(lock.file)) {
      msg <-
        "</br><div class='alert alert-dismissible alert-info'>
      <button type='button' class='close' data-dismiss='alert'>&times;</button>
      <p><strong>A process is running in the background! Please wait</strong> until the process is finished.</p>"
      # print ok when criteria are met
      return(msg)
    } else {
      return(NULL)
    }
  }
  
  
  get_finished_status <- function() {
    if (isolate(status$eden_finished)) {
      msg <-
        "</br><div class='alert alert-dismissible alert-info'>
      <button type='button' class='close' data-dismiss='alert'>&times;</button>
      <p><strong>Eden finished!</stong>  <a href='../eden-visualizer' class='alert-link'> download and visualize the results</p>"
      # clean up
      unlink(faa.path, recursive = T, force = T)
      unlink(ffn.path, recursive = T, force = T)
      unlink(fasta.path, recursive = T, force = T)
      unlink(hmm.path, recursive = T, force = T)
      unlink(sample.path, recursive = T, force = T)
      unlink("/home/eden/data/groups.txt",
             recursive = T,
             force = T)
      unlink(log.path, recursive = T, force = T)
      system2("echo", paste('";;server ready" >> ', log.path, sep = ""))
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
      status <-
        reactiveValues() # save the run status as a reactive object
      input <-
        reactiveValues() # save the run status as a reactive object
      return(msg)
    } else {
      return(NULL)
    }
  }
  
  # ==============================================================================================
  # execution of shell scripts and update of statusmsg print
  # ==============================================================================================
  
  # use observe event to check if button pressed
  observeEvent(input$checkButton, {
    #  file.create(lock.file)
    
    system(
      paste(
        "/home/eden/start_check.sh",
        faa.path,
        ffn.path,
        input$eden_run_cpus,
        input$eden_run_name,
        input$eden_run_gap / 100,
        hmm.path,
        sample.path
      ),
      wait = FALSE
    )
    
    # unlink(faa.path, recursive = T, force = T)
    # unlink(ffn.path, recursive = T, force = T)
    # unlink(fasta.path, recursive = T, force = T)
    # unlink(hmm.path, recursive = T, force = T)
    # unlink(sample.path, recursive = T, force = T)
    # unlink("/home/eden/data/groups.txt", recursive = T, force = T)
    # unlink(log.path, recursive = T, force = T)
    # system2("echo", paste('";;server ready" >> ', log.path, sep = ""))
    # if (!file.exists(fasta.path)) {
    #   dir.create(fasta.path)
    # }
    # if (!file.exists(faa.path)) {
    #   dir.create(faa.path)
    # }
    # if (!file.exists(ffn.path)) {
    #   dir.create(ffn.path)
    # }
    # Sys.chmod(fasta.path, mode = "0777", use_umask = TRUE)
    # Sys.chmod(ffn.path, mode = "0777", use_umask = TRUE)
    # Sys.chmod(faa.path, mode = "0777", use_umask = TRUE)
    #
  })
  
  # Initialize log
  my_log <<- get_new_log()
  
  # Function to update my_data
  update_log <- function() {
    my_log <<- get_new_log()
  }
  
  output$log = renderText({
    invalidateLater(millis = 1000, session)
    update_log()
  })
  
  output$text1 <-  renderText({
    "</br>"
  })
  
  output$text2 <- renderText({
    paste(status$files)
  })
  
  #
  # ==============================================================================================
  # clean up on close
  # ==============================================================================================
  
  ### clean up routine
  cancel.onSessionEnded <- session$onSessionEnded(function() {
    print("SessionEnded")
    unlink(faa.path, recursive = T, force = T)
    unlink(ffn.path, recursive = T, force = T)
    unlink(fasta.path, recursive = T, force = T)
    unlink(hmm.path, recursive = T, force = T)
    unlink(sample.path, recursive = T, force = T)
    unlink("/home/eden/data/groups.txt",
           recursive = T,
           force = T)
    unlink(log.path, recursive = T, force = T)
    system2("echo", paste('";;server ready" >> ', log.path, sep = ""))
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
    status <-
      reactiveValues() # save the run status as a reactive object
    input <-
      reactiveValues() # save the run status as a reactive object
    #    js$reset()
  })
})