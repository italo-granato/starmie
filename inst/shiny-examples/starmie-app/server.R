library(shiny)
library(starmie)

function(input, output) {


  app_data <- reactive({
    if(input$program == "structure") {
      if (is.null(input$out_files$datapath) & is.null(input$log_files$datapath)) {
        app_data <- NULL
      } else if (!is.null(input$out_files$datapath) & is.null(input$log_files$datapath)) {
        app_data <- structList(mapply(loadStructure, input$out_files$datapath,
                                                SIMPLIFY = FALSE, USE.NAMES = FALSE))
      } else if (!is.null(input$out_files) & !is.null(input$log_files)) {
        if (length(input$out_files$datapath) == length(input$log_files$datapath)) {
          app_data <- structList(mapply(loadStructure, sort(input$out_files$datapath), sort(input$log_files$datapath),
                                                  SIMPLIFY = FALSE, USE.NAMES = FALSE))
        }
        else {
          stop("Invalid inputs. Please try again.")
        }
      }
    }

    else if (input$program == "admixture") {
      app_data <- "coming soon"
    }

    else if (input$program == "fineStructure") {
      app_data <- "coming soon"
    }
  })
  output$contents <- renderText({ print(app_data())})

}
