library(ggvis)
library(gridExtra)
library(igraph)
library(ggplot2)
library(tidyr)
library(dplyr)
library(shiny)
theme_set(theme_bw())
load("pathway.rda")
load("graph_prototype.rda")

data = data[which(!sapply(data, function(x) is.null(x)))]
pathways = sapply(data, function(x) x$name)

everythingDF = do.call(rbind,lapply(data, function(x) x$df)) %>% unique

shinyServer(
  function(input, output){
    eve = reactive({
      everythingDF %>% filter(percentage==input$Nvar)
    })
    output$SelectedPathway = renderText({
         pathwayID = unlist(strsplit(input$pathwayID, ": "))[2]
    })

    df = reactive({
         pathwayID = unlist(strsplit(input$pathwayID, ": "))[1]
          df = data[[which(pathways == pathwayID)]]
          df$df = df$df %>% filter(percentage==input$Nvar)

        label = V(df$graph)$label

        if(input$ko==FALSE){
            V(df$graph)$label[grep("ko:", V(df$graph)$name)] = ""
        }
        if(input$ko==TRUE){
            V(df$graph)$label[grep("ko:", V(df$graph)$name)] = label[grep("ko:", V(df$graph)$name)]
        }
        if(input$koID==TRUE){
            V(df$graph)$label[grep("ko:", V(df$graph)$name)] = V(df$graph)$name[grep("ko:", V(df$graph)$name)]
        }
        if(input$cpd==FALSE){
            V(df$graph)$label[grep("cpd:", V(df$graph)$name)] = ""
        }
        V(df$graph)$label.cex=rep(input$cpdtextsize, vcount(df$graph))
        V(df$graph)$label.cex[grep("ko", V(df$graph)$name)] = input$kotextsize

        V(df$graph)$size=rep(input$cpdsize, vcount(df$graph))
        V(df$graph)$size[grep("ko", V(df$graph)$name)] = input$kosize
          df
    })


    #bar = reactive({
        #pathwayID = unlist(strsplit(input$pathwayID, ": "))[1]
        #df = data[[which(pathways == pathwayID)]]
        #bar = df$df                                 %>%
        #filter(percentage==input$Nvar)              %>%
        #select(-nPerc, -percentage, -X, -totalMRNA) %>%
        #gather(type, count, -ko)                    %>%
        #mutate(color=paste(ko, type, sep="_"))
        #bar
    #})
    ##################################################
    #Graph
    ##################################################
    output$graphObj  = renderPlot({
        plot(df()$graph)
      },
      height = 1000,
      width = 1000
    )
    ##################################################
    #Plot
    ##################################################
    output$summary = renderPlot({
#        p0 = ggplot(df()$df,aes(totalMRNA, contigsRequired))+
#        geom_point(aes(color=ko, size=sqrt(contigsRequired)),alpha=1)+
#        geom_point(aes(color=ko, size=sqrt(totalNum)),alpha=0.5)+
#        #geom_text(aes(label=totalNum), size=3)+
#        xlab("mRNA count")+
#        ylab(
#             sprintf("N%s",
#        gsub("\\.", "",format(round(as.numeric(input$Nvar), 2), nsmall=2))
#             ))+
#        scale_color_discrete("KO")+
#        scale_size_continuous("Sq root number of contigs")

        #p1 =    ggplot(bar(),aes(ko,log10(count), group=type, alpha=type))+
                #geom_bar(stat="identity")+
                #theme(axis.text.x=element_text(angle=90))+
                #scale_alpha_discrete(labels=c("Required for NXX", "Total number of contigs"))
                #ggtitle("")

        raw =   ggplot(eve(), aes(x=contigsRequired, y=totalMRNA))+
                geom_point()+
                ggtitle("NXX against expression (All KOs)")
        pathwayType = mapply(function(type, colorDot){
                KOI = unique(subset(newPath, cat.type == type)$ko)
                doi = filter(eve(), ko %in% KOI)
                ggplot(doi,aes(y=contigsRequired, x=totalMRNA)) +
                    geom_point(color=colorDot)+
                    ggtitle(sprintf("%s", type))
        }, SIMPLIFY=FALSE,
        type     =  unique(newPath$cat.type),
        colorDot =  c("#e41a1c","#377eb8","#4daf4a","#984ea3")
        )
        pathwayType$raw = raw
        do.call(grid.arrange,c(pathwayType, list(nrow=3)))
      })
    ##################################################
    #Table
    ##################################################
      output$table = renderDataTable({
        koDF = df()$df
        koDF %>% filter(percentage==input$Nvar)
          })

    ##################################################
    #ggvis
    ##################################################

    vis = reactive({
            path  = df()$df                             %>%
                    select(-percentage, -nPerc, -X)     %>%
                    gather(type, size, -ko, -totalMRNA)
            path2 = df()$df %>%
                    select(ko, contigsRequired)         %>%
                    unique
            pathDone = merge(path, path2, all=T)
            pathDone$id = 1:nrow(pathDone)
            #print(pathDone)
            lb = linked_brush(keys = 1:nrow(pathDone), "red")
            pathDone %>%
            ggvis(x=~totalMRNA, y=~contigsRequired, fill=~ ko, size=~size, size.brush:= 1000)                                                                        %>%
            layer_points(opacity= ~type)                           %>%
            scale_numeric("size", range=c(0,1000))                                                                                                                                %>%
            scale_ordinal("opacity", range=c(1,0.2))                                                                                                                                %>%
            add_tooltip(function(data){paste0("Expression: ", data$totalMRNA, "<br>", "ContigsRequired: ",filter(pathDone[which(pathDone$ko == data$ko),], type == 'contigsRequired')$size,  "<br>","TotalNum: ",filter(pathDone[which(pathDone$ko == data$ko),], type == 'totalNum')$size, "<br>", "KO: ", data$ko)}, "hover") %>%
            hide_legend("fill")                                                                                                                                                   %>%
            lb$input()
        })
    vis %>% bind_shiny("ggvis", "plot_ui")
  })
