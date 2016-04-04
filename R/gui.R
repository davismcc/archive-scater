#' scater GUI function
#'
#' scater shiny app GUI for workflow for less programmatically inclined users
#'
#' @param sce_set SCESet object after running \code{\link{calculateQCMetrics}} 
#' on it
#'
#' @return Opens a browser window with an interactive shiny app and visualize
#' all possible plots included in the scater
#'
#' @import shiny shinydashboard
#'
#' @export
#' 
#' @examples 
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
#' example_sceset <- example_sceset[!drop_genes, ]
#' example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:40)
#' \dontrun{
#' scater_gui(example_sceset)
#' }
scater_gui <- function(sce_set) {
    
    pd <- names(pData(sce_set))
    pd.plot <- pd[!grepl("filter_", pd) & !grepl("is_", pd)]
    
    exprs_values <- c("counts", "exprs", "tpm", "fpkm")
    
    shinyApp(
        ui <- dashboardPage(
            dashboardHeader(title = "scater"),
            dashboardSidebar(
                sidebarMenu(
                    menuItem("plot", tabName = "plot"),
                    menuItem("plotQC", tabName = "plotQC")
                )
            ),
            dashboardBody(
                tabItems(
                    tabItem(tabName = "plot",
                            fluidRow(
                                box(HTML("<h4>Overview of expression for each cell</h4>
                                         Plots produced by this function are intended 
                                         to provide an overview of large-scale 
                                         differences between cells. For each cell, 
                                         the features are ordered from most-expressed 
                                         to least-expressed and the cumulative 
                                         proportion of the total expression for 
                                         the cell is computed across the top 
                                         nfeatures features. These plots can flag 
                                         cells with a very high proportion of the 
                                         library coming from a small number of features; 
                                         such cells are likely to be problematic 
                                         for analyses. Using the colour and blocking 
                                         arguments can flag overall differences in 
                                         cells under different experimental conditions 
                                         or affected by different batch and other variables."),
                                    width = 12,
                                    status = "success")
                            ),
                            fluidRow(
                                column(width = 8,
                                       box(plotOutput("plot", height = 700),
                                           width = NULL
                                       )
                                ),
                                column(width = 4,
                                       selectInput("block1", "block1:",
                                                   pd.plot,
                                                   selected = pd.plot[2]),
                                       selectInput("block2", "block1:",
                                                   pd.plot,
                                                   selected = pd.plot[3]),
                                       selectInput("colour_by", "colour_by:",
                                                   pd.plot,
                                                   selected = pd.plot[4]),
                                       selectInput("exprs_values", "exprs_values:",
                                                   exprs_values)
                                )
                            )
                    ),
                    tabItem(tabName = "plotQC",
                            fluidRow(
                                box(HTML("<h4>General plots</h4>
                                         <b>highest-expression</b> shows features with 
                                         highest expression<br>
                                         <b>explanatory-variables</b> shows a set of 
                                         explanatory variables plotted against each other, 
                                         ordered by marginal variance explained<br>
                                         <b>exprs-mean-vs-freq</b> plots the mean expression 
                                         levels against the frequency of expression for a 
                                         set of features"),
                                    width = 12,
                                    status = "success")
                            ),
                            fluidRow(
                                box(plotOutput("plotQC", height = 600), width = 8),
                                box(
                                    radioButtons("QCtype",
                                                 label = "Choose a type of QC plot",
                                                 choices = c("highest-expression",
                                                             "explanatory-variables",
                                                             "exprs-freq-vs-mean"),
                                                 selected = "highest-expression"),
                                    width = 4
                                )
                            ),
                            fluidRow(
                                box(HTML("<h4>Find PCs</h4>
                                         This plot shows the most important principal 
                                         components for a given variable"),
                                    width = 12,
                                    status = "success"),
                                box(plotOutput("plotQCfindpc", height = 600), width = 8),
                                box(
                                    radioButtons("QCvar",
                                                 label = "Choose a variable of interest",
                                                 choices = pd.plot,
                                                 selected = "total_features"),
                                    width = 4
                                )
                            )
                    )
                )
            )
        ),
        server <- function(input, output, session) {
            output$plot <- renderPlot({
                plot(sce_set, exprs_values = input$exprs_values,
                     block1 = input$block1,
                     block2 = input$block2,
                     colour_by = input$colour_by)
            })
            output$plotQC <- renderPlot({
                plotQC(sce_set, type = input$QCtype)
            })
            output$plotQCfindpc <- renderPlot({
                plotQC(sce_set, type = "find-pcs", variable = input$QCvar)
            })
            session$onSessionEnded(function() {
                stopApp()
            })
        },
        options = list(launch.browser = TRUE)
    )
}
