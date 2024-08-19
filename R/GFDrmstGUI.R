#source("Hilfsfunktionen.R")
#source("test_functions.R")
#source("plot_function.R")
#source("summary.GFDrmst.R")

#library('shinyMatrix')
#library('shinyWidgets')

GFDrmstGUI <- function(){
  requireNamespace("shiny", quietly = TRUE)
  if (!("package:shiny" %in% search())) {
    attachNamespace("shiny")
  }
  requireNamespace("shinyMatrix", quietly = TRUE)
  if (!("package:shinyMatrix" %in% search())) {
    attachNamespace("shinyMatrix")
  }
  requireNamespace("tippy", quietly = TRUE)
  if (!("package:tippy" %in% search())) {
    attachNamespace("tippy")
  }
  ui <- fluidPage(theme = shinythemes::shinytheme("cerulean"),
                  shinyjs::useShinyjs(), titlePanel("Tests for GFDrmst"),
                  sidebarLayout(sidebarPanel(splitLayout(cellWidths = c("40%","15%","20%","20%","0%"), fileInput("infile",
                                                                   "Choose CSV File", accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
                                                         actionButton("infoButton2", "", icon = icon("info-circle")),
                                                         checkboxInput("header", "Header", TRUE),
                                                         selectInput("sep", "Separator in csv", c(",", ";", ".", "|")),
                                                         tippy::tippy_this("infoButton2", "The csv file should be structured in the following way: different rows must contain the data for the independent individuals, different columns contain the different values (at least the event time, censoring status, and one factor variable).   ")
                  ),
                                             tags$head(tags$style(HTML("\n .shiny-split-layout > div {\n overflow: visible;\n }\n"))),
                                             tags$style(HTML("\n input[type=number] {\n -moz-appearance:textfield;\n }\n input[type=number]::{\n -moz-appearance:textfield;\n }\n
                  input[type=number]::-webkit-outer-spin-button,\n input[type=number]::-webkit-inner-spin-button {\n -webkit-appearance: none;\n margin: 0;\n }\n")),
                                             h3(id = "titleLoadData", "Load dataset first!", style = "color:red"),
                                             splitLayout(cellWidths = c("80%"), shinyjs::hidden(selectInput("Method", "Select Testing Method:", c(`Asymptotic` = "asymptotic", `Groupwise Bootstrap` = "groupwise", `Permutation` = "permutation"), selected = "groupwise"))),
                                             splitLayout(cellWidths = c("50%", "50%"), uiOutput(outputId = "dynamicInput"), uiOutput(outputId = "dynamicInput2")),
                                             splitLayout(cellWidths = c("50%"), shinyjs::hidden(checkboxInput("stepwise", "Stepwise extension", TRUE))),
                                             splitLayout(cellWidths = c("50%"), shinyjs::hidden(checkboxInput("plots", "Plot the confidence intervals", FALSE))),
                                             splitLayout(cellWidths = c("50%", "5%", "45%"), uiOutput(outputId = "dynamicInput3"),
                                                         #shinyjs::hidden(textInput("formula", "Formula ", "timeFactor ~ FactorA * FactorB")),
                                                         shinyjs::hidden(actionButton("infoButton", "", icon = icon("info-circle"))),
                                                         tippy::tippy_this("infoButton", "Example:<br><br>\n - 1 Factor: <br>\n  time ~ factorA <br><br>\n
                                    - 2 Factors:\n <br> time ~ factorA * factorB \n <br> or\n <br> time ~ factorA + factorB  ",
                                                                           placement = "right")),
                                             shinyjs::hidden(h5(id = "titleWeights", strong("Further Arguments"), style = "color:grey")),
                                             splitLayout(cellWidths = c("30%"), shinyjs::hidden(selectInput("hyp_mat", "Select contrast matrix:", c(`Dunnett` = "Dunnett", `Tukey` = "Tukey", `Center` = "center", `crossed factorial` = "crossed factorial", `Other` = "Other"), selected = "Tukey"))),
                                             # Other
                                             splitLayout(cellWidths = c("70%","30%"), uiOutput(outputId = "out_mat"), uiOutput(outputId = "out_vec")),
                                             conditionalPanel(condition = "input.hyp_mat == 'Other'",
                                                              splitLayout(cellWidths = c("50%","50%"),
                                                                          actionButton(inputId = "add", label = "Add matrix"),
                                                                          conditionalPanel(condition = "input.add > input.del",actionButton(inputId = "del", label = "Delete matrix"))
                                                              )),

                                             splitLayout(cellWidths = c("15%"), shinyjs::hidden(numericInput("tau", "Endpoint tau of the relevant time window [0,tau]", 1))),
                                             splitLayout(cellWidths = c("40%", "40%"), shinyjs::hidden(autonumericInput("alpha", "Level of significance alpha", value = 0.05)),
                                                         shinyjs::hidden(numericInput("nres", "Number of resampling repetitions", value = 4999))),
                                             shinyjs::hidden(actionButton("process", "Calculate", class = "btn-primary")), width = 6),
                                mainPanel(verbatimTextOutput("group"), verbatimTextOutput("result"), plotOutput("result_plot"), width = 6)))

  server <- function(input, output, session){
    datasetInput <- reactive({
      req(input$infile)
      if (is.null(input$infile))
        return(NULL)
      read.csv(input$infile$datapath, header = input$header,
               sep = as.character(input$sep))
    })
    observeEvent(input$infile, {
      if (is.null(input$infile)) {
        shinyjs::hide(id = "Method")
        #shinyjs::hide(id = "formula0")
        shinyjs::hide(id = "formula")
        shinyjs::hide(id = "infoButton")
        shinyjs::hide(id = "hyp_mat")
        shinyjs::hide(id = "add")
        shinyjs::hide(id = "del")
        shinyjs::hide(id = "alpha")
        shinyjs::hide(id = "nres")
        shinyjs::hide(id = "tau")
        shinyjs::hide(id = "process")
        shinyjs::hide(id = "titleWeights")
        shinyjs::hide(id = "stepwise")
        shinyjs::hide(id = "plots")
      }
      else {
        shinyjs::show(id = "Method")
        shinyjs::show(id = "formula")
        shinyjs::show(id = "infoButton")
        shinyjs::show(id = "stepwise")
        shinyjs::show(id = "process")
        shinyjs::show(id = "titleWeights")
        shinyjs::show(id = "hyp_mat")
        shinyjs::show(id = "tau")
        shinyjs::show(id = "alpha")
        shinyjs::show(id = "nres")
        shinyjs::hide(id = "titleLoadData")
        observeEvent(input$stepwise, {
          if (input$stepwise == TRUE){
            shinyjs::hide(id = "plots")
            updateCheckboxInput(inputId = "plots", value = FALSE)
          }
          else{
            shinyjs::show(id = "plots")
          }
        })
      }
    })
    values <- reactiveValues()
    output$dynamicInput <- renderUI({
      selectInput(inputId = "dynamic", label = "Name of censoring status variable",
                  choices = colnames(datasetInput()))
    })
    observe({ if (input$Method == "asymptotic" || input$Method == "groupwise" || input$Method == "permutation") {
      values$dyn <- input$dynamic
    }
      else {
        values$dyn <- NULL
      }
    })
    values2 <- reactiveValues()
    output$dynamicInput2 <- renderUI({
      selectInput(inputId = "dynamic2", label = "Label of censored variable",
                  choices = unique(datasetInput()[, values$dyn]),
                  selected = suppressWarnings(min(unique(datasetInput()[, values$dyn]))))
    })
    output$dynamicInput3 <- renderUI({
      if(length(grep("time", colnames(datasetInput())[colnames(datasetInput()) != input$dynamic])) > 0){
        textInput(inputId = "dynamic3", label = "Formula",
                  value = suppressWarnings(paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][(grep("time", colnames(datasetInput())[colnames(datasetInput()) != input$dynamic])[1])],
                                                 "~", paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][-(grep("time", colnames(datasetInput())[colnames(datasetInput()) != input$dynamic])[1])], collapse = " * ", sep = " ")) ))
      }
      else{
        textInput(inputId = "dynamic3", label = "Formula",
                  value = suppressWarnings(paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][1],
                                                 "~", paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][-1], collapse = " * ", sep = " ")) ))
      }
    })
    # Define a reactive value to track the formula
    formula_reactive <- reactive({
      req(input$dynamic3)  # Ensure the formula input is available

      # Safely return the formula if it's valid
      if (nchar(input$dynamic3) > 0) {
        tryCatch({
          as.formula(input$dynamic3)  # Convert to formula
        }, error = function(e) {
          NULL  # Return NULL if invalid formula
        })
      } else {
        NULL
      }
    })

    observe({
      values2$dyn <- input$dynamic2
      values2$dyn3 <- input$dynamic3
    })
    values3 <- reactiveValues()

    output$out_mat <- renderUI({
      if (input$hyp_mat == "Other") {
        data <- as.data.frame(datasetInput())
        event <- data[, input$dynamic]
        data[, input$dynamic] <- ifelse(event == input$dynamic2, 0, 1)
        k    <- length(unique(formula2input(input$dynamic3, input$dynamic, data)$group))
        lapply(1:100, function(i){
          shinyjs::hidden( matrixInput(inputId = paste0("new_mat", i), label = paste0("H_",i," ="),
                                       value = matrix("", 1, k), rows = list(extend = TRUE, names = FALSE),
                                       cols = list(n = k, names = FALSE, extend = TRUE, delta = 0), class = "numeric"))
        })
      }else {
        return(NULL)
      }
    })


    output$out_vec <- renderUI({
      vecinputs <- lapply(1:100, function(i){
        if (input$hyp_mat == "Other") {
          return(matrixInput(inputId = paste0("new_vec", i), label = paste0("c_",i," ="),
                             value = matrix(rep("",1)),
                             rows = list(extend = TRUE, names = FALSE),
                             cols = list(n = 1, names = FALSE, extend = TRUE, delta = 0), class = "numeric") )     }
        else {
          return(NULL)
        }
      })
      if(input$hyp_mat == "Other") shinyjs::hidden(vecinputs[!sapply(vecinputs, is.null)])
    })


    observeEvent(eventExpr = input$add, handlerExpr = {
      shinyjs::show(id = paste0("new_mat", input$add - input$del))
      shinyjs::show(id = paste0("new_vec", input$add  - input$del))
    })
    observeEvent(eventExpr = input$del, handlerExpr = {
      shinyjs::hide(id = paste0("new_mat", input$add - input$del + 1))
      shinyjs::hide(id = paste0("new_vec", input$add  - input$del + 1))
    })

    observe({
      if(input$hyp_mat == "Other"){
        for(i in 1:(input$add - input$del)){
          values3[[paste0("new_mat", i)]] <- input[[paste0("new_mat", i)]]
          values3[[paste0("new_vec", i)]] <- input[[paste0("new_vec", i)]]
        }}
      else {
        for(i in 1:(input$add - input$del)){
          values3[[paste0("new_mat", i)]] <- NULL
          values3[[paste0("new_vec", i)]] <- NULL
        }
      }
    })

    observeEvent(input$process, {
      if (length(unique(as.data.frame(datasetInput())[, input$dynamic])) != 2) {
        output$result <- renderPrint({
          "ERROR: More or less than two censoring types"
        })
      }else{
        if(input$hyp_mat == "Other" && input$add - input$del == 0){
          output$result <- renderPrint({
            "ERROR: Specify at least one contrast matrix."
          })}
        else {
          if(input$hyp_mat == "Other" && any(unlist(sapply(1:(input$add - input$del),
                                                    function(i){
                                                      (apply(input[[paste0("new_mat",i)]], 1, function(row){
                                                      if(length(row) > 1){
                                                        return(anyNA(row) && !all(is.na(row)))
                                                      }else{return(FALSE)}
                                                    }) )
                                                    })))){
            output$result <- renderPrint({
              "ERROR: Matrices are not completely filled in."
            })}else{
              if(input$hyp_mat == "Other" && any((sapply(1:(input$add - input$del), function(i) (nrow(na.omit(input[[paste0("new_mat",i)]])) == 0))))){
                output$result <- renderPrint({
                  "ERROR: Delete empty matrices."
                })}
              else {
                if(input$hyp_mat == "Other" && any(sapply(1:(input$add - input$del),
                                                          function(i){
                                                            (nrow(na.omit(input[[paste0("new_mat",i)]])) !=
                                                                       length(na.omit(input[[paste0("new_vec",i)]])))}))  ){
                  output$result <- renderPrint({
                    "ERROR: The vectors c must have the same length as number of rows in H."
                  })}else{
                    if(input$hyp_mat == "Other"){
                      my_hyp_mat <- lapply(1:(input$add - input$del), function(i) na.omit(input[[paste0("new_mat",i)]]))
                      my_hyp_vec <- lapply(1:(input$add - input$del), function(i) na.omit(input[[paste0("new_vec",i)]]))
                      solution <- logical(length(my_hyp_mat))
                      for(i in 1:length(my_hyp_mat)){
                        solution[i] <- existence(my_hyp_mat[[i]], my_hyp_vec[[i]], input$tau)
                      }
                    }
                    else{
                      my_hyp_mat <- input$hyp_mat
                      my_hyp_vec <- NULL
                      solution <- TRUE
                    }
                    if(any(!solution)){
                      output$result <- renderText({
                        paste0("ERROR: Hypothesis", which(!solution),
                               " does not have a possible solution in [0,tau]^k.")
                      })
                    }
                    else{
                      data <- as.data.frame(datasetInput())
                      event <- data[, input$dynamic]
                      data[, input$dynamic] <- ifelse(event == input$dynamic2, 0,1)

                      output$group <- renderPrint({
                        cat("Group assignment: \n")
                        print(formula2groups(formula = isolate(input$dynamic3),
                                             event = input$dynamic, data = isolate(data)),
                              row.names = FALSE)
                      })

                      showModal(modalDialog("Calculating!"))

                      if (input$Method == "asymptotic") {
                        output_asy <- RMST.asymptotic.test(formula = isolate(input$dynamic3), event = input$dynamic, Nres = input$nres,
                                                           data = isolate(data), hyp_mat = my_hyp_mat, hyp_vec = my_hyp_vec,
                                                           tau = input$tau, stepwise = input$stepwise, alpha = input$alpha)
                        removeModal()
                        output$result <- renderPrint({
                          summary.GFDrmst(output_asy)
                        })
                        if (input$plots) {
                          output$result_plot <- renderPlot({
                            plot.GFDrmst(output_asy)
                          })
                        }
                      }
                      if (input$Method == "groupwise") {
                        output_grp <- RMST.groupwise.test(formula = isolate(input$dynamic3),
                                                          event = input$dynamic, Nres = input$nres,
                                                          data = isolate(data), hyp_mat = my_hyp_mat, hyp_vec = my_hyp_vec,
                                                          tau = input$tau, stepwise = input$stepwise, alpha = input$alpha)

                        removeModal()
                        output$result <- renderPrint({
                          summary.GFDrmst(output_grp)
                        })
                        if (input$plots) {
                          output$result_plot <- renderPlot({
                            plot.GFDrmst(output_grp)
                          })
                        }
                      }
                      if (input$Method == "permutation") {
                        output_perm <- RMST.permutation.test(formula = isolate(input$dynamic3),
                                                             event = input$dynamic, Nres = input$nres,
                                                             data = isolate(data), hyp_mat = my_hyp_mat, hyp_vec = my_hyp_vec,
                                                             tau = input$tau, stepwise = input$stepwise, alpha = input$alpha)

                        removeModal()
                        output$result <- renderPrint({
                          summary.GFDrmst(output_perm)
                        })
                        if (input$plots) {
                          output$result_plot <- renderPlot({
                            plot.GFDrmst(output_perm)
                          })
                        }
                      }
                    }
                  }}}}}
      #}
    })
  }
  shinyApp(ui = ui, server = server)
}
GFDrmstGUI()

