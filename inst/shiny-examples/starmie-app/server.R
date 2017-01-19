library(shiny)
library(starmie)

function(input, output, session) {

  app_data <- reactive({
    if(input$program == "structure") {
      if (is.null(input$out_files$datapath) &
          is.null(input$log_files$datapath)) {

        app_data <- NULL

      } else if (!is.null(input$out_files$datapath) &
                 is.null(input$log_files$datapath)) {

        app_data <- structList(mapply(loadStructure, input$out_files$datapath,
                                                SIMPLIFY = FALSE, USE.NAMES = FALSE))

      } else if (!is.null(input$out_files) & !is.null(input$log_files)) {
        if (length(input$out_files$datapath) == length(input$log_files$datapath)) {

          app_data <- structList(mapply(loadStructure,
                                        sort(input$out_files$datapath), sort(input$log_files$datapath),
                                        SIMPLIFY = FALSE, USE.NAMES = FALSE))
        }
        else {

          stop(safeError("Invalid inputs. Please try again."))

        }
      }
    }

    else if (input$program == "admixture") {
      if (is.null(input$Q_files) & is.null(input$P_files) &
          is.null(input$alog_files) & is.null(input$Q_bias_files) &
          is.null(input$Q_se_files)) {

        app_data <- NULL

      } else if (!is.null(input$Q_files) & !is.null(input$P_files) &
                 is.null(input$alog_files) & is.null(input$Q_bias_files) &
                 is.null(input$Q_se_files)) {

        P_files <- sort(input$P_files$datapath)
        Q_files <- sort(input$Q_files$datapath)
        if (length(P_files) == length(Q_files)) {

          app_data <- admixList(mapply(loadAdmixture, Q_files, P_files,
                                       SIMPLIFY = FALSE, USE.NAMES = FALSE))
        } else {
          stop(safeError("Number of Q files must equal number of P files."))
        }

      } else if (!is.null(input$Q_files) & !is.null(input$P_files) &
                 !is.null(input$alog_files) & is.null(input$Q_bias_files) &
                 is.null(input$Q_se_files)) {

        P_files <- sort(input$P_files$datapath)
        Q_files <- sort(input$Q_files$datapath)
        alog_files <- sort(input$alog_files$datapath)
        if (length(alog_files) == length(Q_files) &
            length(alog_files) == length(P_files)) {
          app_data <- admixList(mapply(loadAdmixture, Q_files, P_files,
                                       alog_files, SIMPLIFY = FALSE,
                                       USE.NAMES = FALSE))
        } else {
          stop(safeError("Number of log files must match Q and P files."))
        }


      } else if (!is.null(input$Q_files) & !is.null(input$P_files) &
                 !is.null(input$alog_files) & !is.null(input$Q_bias_files) &
                 !is.null(input$Q_se_files)) {

        P_files <- sort(input$P_files$datapath)
        Q_files <- sort(input$Q_files$datapath)
        alog_files <- sort(input$alog_files$datapath)
        Qb_files <- sort(input$Q_bias_files$datapath)
        Qs_files <- sort(input$Q_se_files$datapath)

        app_data <- admixList(mapply(loadAdmixture, Q_files, P_files,
                                     alog_files, Qb_files, Qs_files,
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE))

      } else if (!is.null(input$Q_files) & !is.null(input$P_files) &
                 is.null(input$alog_files) & !is.null(input$Q_bias_files) &
                 !is.null(input$Q_se_files)) {

        P_files <- sort(input$P_files$datapath)
        Q_files <- sort(input$Q_files$datapath)
        Qb_files <- sort(input$Q_bias_files$datapath)
        Qs_files <- sort(input$Q_se_files$datapath)
        app_data <- admixList(mapply(loadAdmixture, qfile = Q_files,
                                     pfile = P_files, bootstrap_bias = Qb_files,
                                     bootstrap_se =  Qs_files,
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE))
      } else {
        stop(safeError("Invalid inputs. The following combinations of admixture files are valid uploads:
                       1. P and Q
                       2. P and Q and log
                       3. P and Q and log and all bootstrap files.
                       4. P and Q and all bootstrap files."))
      }
    }

    else if (input$program == "fineStructure") {
      app_data <- NULL
    }

    else if (input$program == "fastStructure") {
      app_data <- NULL
    }
  })

  summary_values <- reactiveValues()
  summary_values$list_size <- length(isolate(app_data()))
  summary_values$min_k <- min(unlist(lapply(isolate(app_data()), getK)))
  summary_values$max_k <- max(unlist(lapply(isolate(app_data()), getK)))

  output$summary <- renderPrint({ print(app_data())})

  output$bestKcontrols <- renderUI({
    if (input$program == "structure") {
      method <- c("evanno", "structure")
      radioButtons("method", "Choose Method:", method)
    }

  })

  output$inference <- renderPlot({
    method <- input$method
    if (length(method) > 0 & input$program == "structure") {
      bestK(app_data(), method)
    } else if (input$program == "admixture") {
      bestK(app_data())
    } else {
      stop(safeError("Invalid input for bestK function."))
    }

  })

  output$visualiseControls <- renderUI({
    common_plots <- c("bar", "tree", "multiple", "mds")
    if (input$program == "structure") {
      radioButtons("plot_type", "Choose a plot type:", c(common_plots, "mcmc_diagnostics"))
    }
    else if (input$program == "admixture") {
      radioButtons("plot_type", "Choose a plot type:", common_plots)
    }
  })

  output$plotOptions <- renderUI({
    if (input$plot_type %in% c("bar", "tree", "mds")) {
      numericInput("K", "Select K value:", summary_values$min_k,
                   min = summary_values$min_k, max = summary_values$max_k)
    }

  })

  output$plot_out <- renderUI({
    switch(input$plot_type,
           "bar" = plotOutput("bar"),
           "mcmc_diagnostic" = plotOutput("mcmc_diagnostic"),
           "tree" =  plotOutput("tree"),
           "multiple" = plotOutput("multiple"),
           "mds" = plotOutput("mds"))
  })

  output$bar <- renderPlot({
    index_out <- unlist(lapply(app_data(), getK)) == input$K
    if (sum(index_out) > 1) {
      stop(safeError("more than 1 struct object with K-value, use multiple plot instead"))
    }

    input_obj <- app_data()[index_out]

    plotBar(input_obj)

  })

  output$mcmc_diagnostic <- renderPlot({
    can_plot <- any(unlist(lapply(app_data(), function(y) is.null(y[["nonburn_df"]]))))

    if (can_plot) {
      stop(safeError("Some struct objects have no MCMC data. Please upload to plot MCMC chains."))
    }

    plotMCMC(app_data())

  })

  output$tree <- renderPlot({
    index_out <- unlist(lapply(app_data(), getK)) == input$K
    if (sum(index_out) > 1) {
      stop(safeError("more than one object with same K-value, use multiple plot instead"))
    }

    print(class(app_data()))
    input_obj <- app_data()[index_out]
    print(class(input_obj))
    print(index_out)
    plotTreeBar(input_obj)

  })

  output$mds <- renderPlot({
    index_out <- unlist(lapply(app_data(), getK)) == input$K
    if (sum(index_out) > 1) {
      stop(safeError("more than one object with K-value, use multiple plot instead"))
    }

    input_obj <- app_data()[index_out]

    plotMDS(input_obj)


  })

  output$multiple <- renderPlot({
    plotMultiK(app_data())
  })




}
