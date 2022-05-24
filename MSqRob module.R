Genelist_UImodule <- function(id, label = "Counter"){
  
  ns <- NS(id)
  tagList(
    # navbarPage("Genelist tool",
    fluidRow(
      column(width = 12,
             tabBox(title = "Genelist tool", width = 12,
                    tabPanel("Generate gene list by cutoff",
                             column(width = 6,
                                    h3("Significant gene lists by threshold"),
                                    # fluidRow(
                                    #   conditionalPanel(condition = paste("input['", id, "-genelistFrom']", " ==  'DEG' || input.threshold_method == 'intersect'" , sep = ""), #ns = ns, #input.threshold_method == 'intersect' ||  paste("input['", id, "-genelistFrom']", " ==  'DEG'" , sep = "")
                                    #                    box(numericInput(ns("p_genelist"),
                                    #                                     "adjusted P cutoff",
                                    #                                     min = 0.0001, max = 0.1, value = 0.05),
                                    #                        width = 6),
                                    #                    box(numericInput(ns("lfc_genelist"),
                                    #                                     "log2 fold change cutoff",
                                    #                                     min = 0, max = 10, value = 1),
                                    #                        width = 6)  
                                    #   )
                                    #   
                                    #   # box(numericInput(ns("p_genelist"),
                                    #   #                       "adjusted P cutoff",
                                    #   #                       min = 0.0001, max = 0.1, value = 0.05),
                                    #   #          width = 6),
                                    #   # box(numericInput(ns("lfc_genelist"),
                                    #   #                       "log2 fold change cutoff",
                                    #   #                       min = 0, max = 10, value = 1),
                                    #   #          width = 6)                            
                                    # ),
                                    fluidRow(
                                      box(selectizeInput(ns("genelistFrom"),
                                                         "Generate gene list from",
                                                         choices=c("MSqRob"," "), selected = " ",
                                                         multiple = F), width = 6),
                                      # box(selectizeInput(ns("selected_compare"),
                                      #                    "Select groups",
                                      #                    choices=NULL,
                                      #                    multiple = F),width = 4),
                                      box(uiOutput(ns("selected_compare_ui")),width = 6)
                                    ),
                                    fluidRow(
                                      
                                      box(selectizeInput(ns("type_genelist"),
                                                         "Significant type",
                                                         choices=c("Significant(Up&Down)","Upregulated(Up)","Downregulated(Down)"),
                                                         multiple = F), width = 6)                            
                                    ),
                                    fluidRow(
                                      column(width = 12,
                                             actionBttn(inputId = ns("Generate_genelist"),
                                                        label = "Generate",
                                                        style = "jelly",
                                                        color = "primary")                              
                                      )                             
                                    )
                             ),
                             column(width = 5,
                                    br(),
                                    br(),
                                    br(),
                                    h4("existing gene lists"),
                                    htmlOutput(ns("existed_lists"))
                             )
                    ),
                    tabPanel("Import gene list",
                             h3("Import gene list"),
                             fluidRow(
                               column(h4("Import by pasting gene list"),
                                      textAreaInput(inputId = ns("text_input_genelist"),
                                                    label = "Please paste your gene list",
                                                    placeholder = "TP53\nMAGED2\n",
                                                    rows = 8, width = "350px"),
                                      textInput(ns("name_genelist"),"the name of imported gene list","Imported list 1", width = "350px"),
                                      actionBttn(
                                        inputId = ns("Import_genelist"),
                                        label = "Import",
                                        style = "jelly",
                                        color = "primary"
                                      ), 
                                      br(),
                                      br(),
                                      textOutput(ns("input_genelist")),
                                      tableOutput(ns("dataframe"))                                                               
                                      ,width = 5),
                               column(h4("existing imported gene lists"),
                                      htmlOutput(ns("existed_imported_lists")),
                                      width = 7)
                             ),
                             # tags$br(),
                             # textInput(ns("name_genelist"),"the name of import Genelist","Imported list 1", width = "350px"),
                             
                             # br(),
                             # br(),
                             # textOutput(ns("input_genelist")),
                             # tableOutput(ns("dataframe"))
                    )         
                    
             )
      )
    )
    
    # )
    
  )
  
}


Genelist_server_module <- function(id,alpha,object,lfc,condition_matched_frame) {
  moduleServer(
    id,
    function(input, output, session,id=id) {
      
      # rm all objects of the envir that we build, in order to make each run app rm objects last run
      rm(list = ls(envir = get_msqrob_genelist_Env()),envir = get_msqrob_genelist_Env())
      rm(list = ls(envir = get_imported_genelist_Env()),envir = get_imported_genelist_Env())
      rm(list = ls(envir = get_all_genelist_Env()),envir = get_all_genelist_Env())
      rm(list = ls(envir = get_gsea_genelist_Env()),envir = get_gsea_genelist_Env())
      
      msqrob_frame <- reactive({
        msqrob_frame <- result_from_object(object=object(),condition_matched_frame())
        # test_frame <<- msqrob_frame
      })
      
      
      lists_msqrob <- reactive({
        lists_msqrob <- all_msqrob_siglist(frame = msqrob_frame(),alpha = alpha(),lfc = lfc())
        # testlist <<- lists_msqrob
      })
      
      
      output$existed_lists = renderUI({
        validate(need(!is.null(msqrob_frame()), "Please go to the MSqRob options, and do the differential analysis first"))
        if(!is.null(msqrob_frame())) {
          lis <- unique(c(lists_msqrob()))
          # lis <- unique(ls(envir = get_msqrob_genelist_Env()))
        } 
        HTML(paste(lis, collapse = "</br>"))
      })
      
      # 
      # output$existed_lists <- renderUI({
      #   HTML(paste(ls(envir = get_msqrob_genelist_Env()), collapse = "</br>"))
      # })
      
      
      output$input_genelist = renderText({
        thename <- input$name_genelist
      })
      
      #
      #       output$Vennplot = renderPlot({NULL})
      
      
      cat("bb")
      
      output$selected_compare_ui <- renderUI({
        froms <- input$genelistFrom
        if(froms == "MSqRob"){
          validate(need(!is.null(msqrob_frame()), "Please go to the MSqRob options, and do the differential analysis first"))
          frame <- msqrob_frame()
          cols <- grep(".significant", colnames(frame))
          names <- colnames(frame)[cols]
          names <- gsub(".significant", "", names)
          selectizeInput(session$ns("selected_compare"),"Select groups",choices=names, multiple = F)#getDefaultReactiveDomain()
        }
      })
      
      
      ## Event generate  genelist
      observeEvent(input$Generate_genelist,{
        # if(length(input$selected_compare)==0){
        #   sendSweetAlert(
        #     session = session,
        #     title = "Warning !",
        #     text = "please choose compare",
        #     type = "warning"
        #   )
        # }else{
        # aaa <- all_msqrob_siglist(object = msqrob_object,alpha = alpha,lfc = lfc)
        froms = input$genelistFrom
        list_type = input$type_genelist
        # test_type <<- list_type
        type <- switch(list_type,
                       `Significant(Up&Down)` = "Sig",
                       `Upregulated(Up)`= "Up",
                       `Downregulated(Down)` ="Down")
        
        
        # if(froms == "Msqrob"){
        thelist <- get(paste0(input$selected_compare,".",alpha(),"_",lfc(),
                              "_",type),envir = get_msqrob_genelist_Env())
        # testthelist <<- paste0(input$selected_compare,".",alpha,"_",lfc,"_",type) 
        the_name <- paste0(input$selected_compare,".",alpha(),"_",lfc(),
                           "_",type)
        
        if(nrow(thelist)==0){
          sendSweetAlert(
            session = session,
            title = "Warning !",
            text = "no expected protein under this threshold\nPlease back to Msqrob option to do hypothesis test",
            type = "warning"
          )
        }
        # thelist <- get_msqrob_siglist(msqrob_object = msqrob(),
        #                                 alpha = input$p_genelist , lfc = input$lfc_genelist,
        #                                 compare = input$selected_compare,type = list_type) 
        else{
          sendSweetAlert(
            session = session,
            title = "Success",
            text = paste("the list saved as: ", the_name,
                         "  (","total ",nrow(thelist)," proteins",")",
                         sep = ""),
            type = "success"
          )
        }
        
      })
      
      ## Event import  genelist
      observeEvent(input$Import_genelist, {
        texts <- input$text_input_genelist
        ## warning if
        if(nchar(texts) == 0){
          sendSweetAlert(
            session = session,
            title = "Warning !",
            text = "your input text is empty",
            type = "warning"
          )
          import_table = NULL
        }else{
          import_table = try(read.table(text = texts,header = F,sep = "\t",stringsAsFactors = F),silent = T)
          
          ## read in talbe
          if(nrow(import_table) ==1 ){
            ## warning if only have one gene
            sendSweetAlert(
              session = session,
              title = "Warning !",
              text = "please input list with multiple genes",
              type = "warning"
            )
            import_table = NULL
          }else{
            if(ncol(import_table)==1){
              colnames(import_table) = "symbol"
            }else if(ncol(import_table)==2){
              colnames(import_table) = c("symbol","FC")
            }else if(ncol(import_table)>2){
              import_table = import_table[,1:2]
              colnames(import_table) = c("symbol","FC")
            }
            ## assign table
            assign(x = input$name_genelist,value = import_table,envir = get_imported_genelist_Env())
            assign(x = input$name_genelist,value = import_table,envir = get_all_genelist_Env())
            
            sendSweetAlert(
              session = session,
              title = "Success",
              text = paste("list import success: (total ",nrow(import_table),")",sep = ""),
              type = "success"
            )
          }
        }
        output$dataframe <- renderTable({
          import_table
        })
        
        output$existed_imported_lists <- renderUI({
          HTML(paste(ls(envir = get_imported_genelist_Env()), collapse = "</br>"))
        })
        
      })
    }
  )
}