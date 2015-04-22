library(MetamapsDB)
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

data          =  data[which(!sapply(data, function(x) is.null(x)))]
pathways      =  sapply(data, function(x) x$name)
everythingDF  =  do.call(rbind,lapply(data, function(x) x$df)) %>% unique

all_values <- function(x) {
if(is.null(x)) return(NULL)
paste0(names(x), ": ", format(x), collapse = "<br />")
}

shinyServer(
  function(input, output){
    output$SelectedPathway = renderText({
         pathwayID = unlist(strsplit(input$pathwayID, ": "))[2]
    })

    eve = reactive({
      everythingDF %>% filter(percentage==input$Nvar)
    })

    df = reactive({
         pathwayID = unlist(strsplit(input$pathwayID, ": "))[1]
          df = data[[which(pathways == pathwayID)]]
          df$df = df$df %>% filter(percentage==input$Nvar)
        #label = V(df$graph)$label

        #if(input$ko==FALSE){
            #V(df$graph)$label[grep("ko:", V(df$graph)$name)] = ""
        #}
        #if(input$ko==TRUE){
            #V(df$graph)$label[grep("ko:", V(df$graph)$name)] = label[grep("ko:", V(df$graph)$name)]
        #}
        #if(input$koID==TRUE){
            #V(df$graph)$label[grep("ko:", V(df$graph)$name)] = V(df$graph)$name[grep("ko:", V(df$graph)$name)]
        #}
        #if(input$cpd==FALSE){
            #V(df$graph)$label[grep("cpd:", V(df$graph)$name)] = ""
        #}
        #V(df$graph)$label.cex=rep(input$cpdtextsize, vcount(df$graph))
        #V(df$graph)$label.cex[grep("ko", V(df$graph)$name)] = input$kotextsize

        #V(df$graph)$size=rep(input$cpdsize, vcount(df$graph))
        #V(df$graph)$size[grep("ko", V(df$graph)$name)] = input$kosize
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
    #Graph (Igraph Plot)
    ##################################################
    #output$graphObj  = renderPlot({
        #plot(df()$graph)
      #},
      #height = 1000,
      #width = 1000
    #)
    ##################################################
    #Summary Plot:: All KOs
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
        raw =   ggplot(eve(), aes(x=contigsRequired, y=totalMRNA)) +
                geom_point(aes(color=readType), alpha = 0.7)                     +
                scale_color_manual(values=c("grey", "black"))    +
                ggtitle("NXX against expression (All KOs)")

        pathwayType = mapply(function(type, colorDot){

                KOI = unique(subset(newPath, cat.type == type)$ko)
                doi = filter(eve(), ko %in% KOI)

                ggplot(doi,aes(y=contigsRequired, x=totalMRNA))     +
                    geom_point(aes(color=readType), alpha=0.7)      +
                    scale_color_manual(values=c(colorDot, "black")) +
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

    pathDone = reactive({
        #path0 = data[[1]]$df %>% filter(percentage == 0.8)
        path  = df()$df                               %>%
        #path = path0                                          %>%
        filter(readType == input$scatter)        %>%
        select(-percentage, -nPerc, -X, -readType)               %>%
        gather(type, size, -ko, -totalMRNA)

        path2 = df()$df                       %>%
        #path2 = path0                                  %>%
        filter(readType == input$scatter)        %>%
        select(ko, contigsRequired) %>%
        unique

        pathDone = merge(path, path2, all=T)
        pathDone$id = 1:nrow(pathDone)
        pathDone 
    })

    lb = linked_brush(keys = 1:nrow(pathDone()), "red")
    #TODO: need to separate gDNA and cDNA (and i need to know how much mRNA was actually attached to those contigs; second thought, will probably have to dig upstream again ... SIGH
    pathDone                                                                                                                                                                                                                                                                                                            %>% 
#    #ggvis automatically takes the reactive function; so no need to pass it double brackets ()
    ggvis(x=~totalMRNA, y=~contigsRequired, stroke.brush:="red", fill=~ ko, size=~size, size.brush:= 1000)                                                                                                                                                                                                              %>%
    layer_points(opacity= ~type)                                                                                                                                                                                                                                                                                        %>%
    scale_numeric("size", range=c(0,1000))                                                                                                                                                                                                                                                                              %>%
    scale_ordinal("opacity", range=c(1,0.2))                                                                                                                                                                                                                                                                            %>%
    add_tooltip(function(data){paste0("Expression: ", data$totalMRNA, "<br>", "ContigsRequired: ",filter(pathDone()[which(pathDone()$ko == data$ko),], type == 'contigsRequired')$size,  "<br>","TotalNum: ",filter(pathDone()[which(pathDone()$ko == data$ko),], type == 'totalNum')$size, "<br>", "KO: ", data$ko)}, "hover") %>%
    add_axis("x", title = "Total Mappable Reads") %>%
    hide_legend("fill")                                                                                                                                                                                                                                                                                                 %>%
    lb$input() %>%
    bind_shiny("ggvis1")

    selected = lb$selected

    graphObj = reactive({
         ig2ggvis(df()$graph)
    })

    graphHighlight = reactive({
        #cat(
            #as.character(graphObj()[which(graphObj()$name %in% gsub("^","ko:", unique(pathDone()[selected(), ]$ko))),]$name)
            #)
        graphObj()[graphObj()$name %in% gsub("^","ko:", unique(pathDone()[selected(), ]$ko)),]
    })

base = graphObj %>%
            ggvis(~x, ~y)                                                                                                                                                %>%
            group_by(row)                                                                                                                                                           %>%
            layer_paths()                                                                                                                                                           %>%
            layer_points(size= ~type, fill=~type)                                                                                                                                   %>%
            add_axis("x", title = "", properties = axis_props(axis = list(strokeWidth = 0), grid = list(strokeWidth = 0),ticks = list(strokeWidth = 0), labels = list(fontSize=0))) %>%
            add_axis("y", title = "", properties = axis_props(axis = list(strokeWidth = 0), grid = list(strokeWidth = 0),ticks = list(strokeWidth = 0), labels = list(fontSize=0))) %>%
            scale_ordinal("size", range=c(20,100)) %>%
            scale_ordinal("fill", range=c("grey","red")) %>%
            add_data(graphHighlight)                                      %>%   #same here ggvis takes the reactive directly
            layer_points(size:=1000, fill:="yellow")                          %>%
            add_tooltip(all_values, "hover") %>%
            set_options(width = 1000, height = 1000) #padding = padding(20, 20, 20, 20))

       reactive({
            if(input$koShow){
           if(input$koID){
               base_layer1 = base %>% layer_text(data=subset(graphObj(),type == 'ko'),text:=~name, fontSize:= input$kotextsize)
                if(input$cpdShow){
                    if(input$cpd){
                        base_layer1 %>% layer_text(data=subset(graphObj(),type == 'cpd'),text:=~name, fontSize:= input$cpdtextsize)
                    }else{
                        base_layer1 %>% layer_text(data=subset(graphObj(),type == 'cpd'),text:=~label, fontSize:= input$cpdtextsize)
                    }
                }else{base_layer1}
            }else{
                base_layer1 = base %>% layer_text(data=subset(graphObj(),type == 'ko'),text:=~label, fontSize:=input$kotextsize)
                if(input$cpdShow){
                    if(input$cpd){
                        base_layer1 %>% layer_text(data=subset(graphObj(),type == 'cpd'),text:=~name, fontSize:= input$cpdtextsize)
                    }else{
                        base_layer1 %>% layer_text(data=subset(graphObj(),type == 'cpd'),text:=~label, fontSize:= input$cpdtextsize)
                    }
                }else{base_layer1}
           }
            }else{base}
       }) %>% bind_shiny("ggvis2")
#TODO: the graph brushing stops working at pathK01100, not sure why, maybe too many datapoints


#else{
    #reactive({
        #baseG %>%
        #layer_text(data=subset(graphObj(),type == 'ko'),text:=~label, fontSize:=input$kotextsize)
    #}) %>%
    #bind_shiny("ggvis2")
#}
        #layer_text(data=subset(graphObj(),type == 'cpd'),text:=~name, fontSize:= input$cpdtextsize) %>%
                #add_tooltip(function(data){
                    #paste0(
                    #"ID: ", graphObj()[graphObj()$id == data$id,]$name, "<br>",
                    #"label: ", graphObj()[graphObj()$id == data$id,]$label)
                #}) %>%
  })
