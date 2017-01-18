library(shiny)
library(starmie)

function(input, output) {

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
          is.null(input$log_files) & is.null(input$Q_bias_files) &
          is.null(input$Q_se_files)) {

        app_data <- NULL

      } else if (!is.null(input$Q_files) & !is.null(input$P_files) &
                 is.null(input$log_files) & is.null(input$Q_bias_files) &
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
                 !is.null(input$log_files) & is.null(input$Q_bias_files) &
                 is.null(input$Q_se_files)) {

        P_files <- sort(input$P_files$datapath)
        Q_files <- sort(input$Q_files$datapath)
        log_files <- sort(input$log_files$datapath)
        if (length(log_files) == length(Q_files) &
            length(log_files) == length(P_files)) {
          app_data <- admixList(mapply(loadAdmixture, Q_files, P_files,
                                       log_files, SIMPLIFY = FALSE,
                                       USE.NAMES = FALSE))
        } else {
          stop(safeError("Number of log files must match Q and P files."))
        }


      } else if (!is.null(input$Q_files) & !is.null(input$P_files) &
                 !is.null(input$log_files) & !is.null(input$Q_bias_files) &
                 !is.null(input$Q_se_files)) {

        P_files <- sort(input$P_files$datapath)
        Q_files <- sort(input$Q_files$datapath)
        log_files <- sort(input$log_files$datapath)
        Qb_files <- sort(input$Q_bias_files$datapath)
        Qs_files <- sort(input$Q_se_files$datapath)

        app_data <- admixList(mapply(loadAdmixture, Q_files, P_files,
                                     log_files, Qb_files, Qs_files,
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE))

      } else if (!is.null(input$Q_files) & !is.null(input$P_files) &
                 is.null(input$log_files) & !is.null(input$Q_bias_files) &
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
      app_data <- "coming soon"
    }

    else if (input$program == "fastStructure") {
      app_data <- "comming soon"
    }
  })
  output$contents <- renderText({ print(app_data())})

}
