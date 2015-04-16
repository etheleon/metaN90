library(shiny)
library(ggvis)
load("graph_prototype.rda")
pathwayChoice = do.call(c,sapply(data, function(x) paste(x$name, x$fullName, sep=": ")))

shinyUI(pageWithSidebar(
    headerPanel("Metabolic"),
    sidebarPanel(
    ############
        selectInput('pathwayID', label=h4('Select the Pathway of your choice'), pathwayChoice),
        checkboxInput('koID', 'Display KOID', value=FALSE),
        checkboxInput('ko', 'Display KO label', value=TRUE),
        checkboxInput('cpd', 'Display cpd label', value=TRUE),
        sliderInput('kotextsize', 'KO label size', value=1, min=0.1, max=10, step=0.5),
        sliderInput('cpdtextsize', 'Cpd label size', value = 0.1, min=0.1, max=10, step=0.1),
        sliderInput('kosize', 'KO size', value=4, min=0.1, max=10, step=0.5),
        sliderInput('cpdsize', 'Cpd size', value = 0.1, min=0.1, max=10, step=0.1),
        selectInput('Nvar', label=h4('Select NXX for table'), seq(0.1,0.95,0.05), selected=0.9),
        uiOutput("plot_ui")
      ),
    mainPanel(
    ###########
              h1(textOutput("SelectedPathway")),
              tabsetPanel(
          tabPanel("NXX", ggvisOutput("ggvis")),
          tabPanel("Metabolic Graph", 
                    p("Below shows a subset of selected pathway"),
                    plotOutput("graphObj", height=1000, width=1000)),
          tabPanel("Summary", 
                    h3("Summary"),
                    p("NXX against total expression"),
                   plotOutput("summary", height = 1000,width=1000)), 
          tabPanel("Table", dataTableOutput("table"))
                )
    )
))
