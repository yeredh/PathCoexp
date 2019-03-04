reactiveNetwork <- function (outputId) 
{
  HTML(paste("<div id=\"", outputId, "\" class=\"shiny-network-output\"><svg /></div>", sep=""))
}

# changed account to track [http://spark.rstudio.com/yeredh/PathCoexp01/]
googleAnalytics <- function(account="UA-43248788-1"){
  HTML(paste("<script type=\"text/javascript\">

    var _gaq = _gaq || [];
  _gaq.push(['_setAccount', '",account,"']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

  </script>", sep=""))
}


shinyUI(
  pageWithSidebar(
    headerPanel(title="Pathway Coexpression Network",windowTitle="PathCoexp"),
    
    sidebarPanel(  
      # Select pathway
      includeHTML("www/js/tools.js"),
      selectInput("tool1", label = "Select pathway:", choices = path.names, selected = NULL, multiple = FALSE),
      
      br(),
      
      # Select the number top n connected pathways
      sliderInput("top.n", "Pathways to display:", 
                  min=0, max=50, value=15),

      br(),
      
      # select the BIC cut-off
      sliderInput("BIC.cut","BIC cut-off:",
                  min=0, max=1, value=0.50,step=0.05),

      br(),
      
      # select the correlation cut-off
      sliderInput("cor.cut","PathCor cut-off:",
                  min=0, max=1, value=0.25,step=0.05),
      
      
      br(),
      
      selectInput("method", "Top edges by:", 
                  choices = c("Absolute Value", "Decreasing", "Increasing")),
      br(),
      # Select edges to display in graph
      selectInput("edge.weight", "Edge weight:", 
                  list("Adjusted Correlation"="PathCor",
                       "Raw Correlation"="Cor")),
      
      
      br(),
      
      submitButton("Submit"),
      
      br(),
      
      h4("Legend"),
      img(src = "legend.png")


    ),
    
    mainPanel(
      h3(textOutput("caption")),
      includeHTML("graph.js"),
      googleAnalytics(),
      tabsetPanel(
        tabPanel("Network", plotOutput("graphplot")),
        tabPanel("Table", tableOutput("view")),
        #tabPanel("Interactive Network",reactiveNetwork("net.list")),
        tabPanel("Coefficient Selection",plotOutput("BICplot"),p(coeffSelString))
      )
      
    )
  )
)