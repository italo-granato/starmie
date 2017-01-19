library(shiny)

fluidPage(
  # Application title
  titlePanel("starmie - interactive analysis of population structure models"),
  sidebarLayout(
    sidebarPanel(
      selectInput('program',
                  'Which program did you run?',
                  c("structure", "admixture", "fastStructure", "fineStructure")),

      conditionalPanel(condition = "input.program == 'structure'",
                       fileInput("out_files", "out_f files", multiple = TRUE, accept = c(".out_f")),
                       fileInput("log_files", "STRUCTURE log files (optional)", multiple = TRUE, accept = c(".out", ".log"))
      ),
      conditionalPanel(condition = "input.program == 'admixture'",
                       fileInput("Q_files", "Q files", multiple = TRUE, accept = ".Q"),
                       fileInput("P_files", "P files", multiple = TRUE, accept = ".P"),
                       fileInput("alog_files", "ADMIXTURE log files (optional)", multiple = TRUE, accept = c("text/plain", ".out", ".log")),
                       fileInput("Q_bias_files", "bootstrap bias files (optional)", multiple = TRUE, accept = ".Q_bias"),
                       fileInput("Q_se_files", "bootstrap standard error files (optional)", multiple = TRUE, accept = ".Q_se")
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("starmie-summary", verbatimTextOutput("summary")),
        tabPanel("starmie-inference", plotOutput("inference"), uiOutput("bestKcontrols")),
        tabPanel("starmie-visualise",  uiOutput("visualiseControls"), uiOutput("plotOptions"), uiOutput("plot_out"))
      )
    )
  )
)
