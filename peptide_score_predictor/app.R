#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/

#Shiny app -- predictor

#libraries

library(shiny)
library(shinythemes)
library(tidyverse)
library(glmnet)

#model1 <- readRDS(file = "../data/NN_1.rds")

#Encoding function

encode_peptide = function(x, encoder){
    X_enc = x %>%
        str_split('') %>%
        lapply(function(x_i){
            encoder[x_i,] %>%
                as.matrix %>%
                t() %>%
                matrix(nrow = 1, byrow = TRUE) %>%
                return
        })
    X_enc = do.call(rbind, X_enc)
    rownames(X_enc) = x
    return(X_enc)
}


# Define UI
#------------------------------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(theme = shinytheme("readable"),
                
                #Neural network predictor
                navbarPage(
                    "Shiny peptide fitness",
                    tabPanel("Elastic net predictor",
                             
                             sidebarPanel(
                                 tags$h3("Input:"),
                                 textInput("txt1", "Peptide sequence:", ""),
                                 selectInput("encoding1", h4("Encoding"), 
                                             choices = list("blosum62", "blosum45", "blosum50", "blosum62", "blosum80", "blosum90", "pam30", "pam70",  "pam250","Z-scales", selected = 1)),
                                 selectInput("protein1", h4("Select protein"), 
                                             choices = list("BRCA1", "ERK2", "DLRAP1", "Pab1" ,selected = 1)),
                                 actionButton("submitbutton", "Submit", class = "btn btn-primary")
                             ), # sidebarPanel
                             
                             mainPanel(
                                 h1("Fitness prediction"),
                                 h4("Results"),
                                 verbatimTextOutput("txtout1"),verbatimTextOutput("txtout2"),verbatimTextOutput("txtout3")
                                 ,) # mainPanel
                             
                             
                    ), # Elastic net predictor
                    tabPanel("Artificial neural network predictor",
                             sidebarPanel(
                                 tags$h3("Input:"),
                                 textInput("txt2", "Peptide sequence:", ""),
                                 selectInput("encoding1", h4("Encoding"), 
                                             choices = list("blosum62", "blosum45", "blosum50", "blosum62", "blosum80", "blosum90", "pam30", "pam70",  "pam250","Z-scales", selected = 1)),
                                 selectInput("protein1", h4("Select protein"), 
                                             choices = list("BRCA1", "ERK2", "DLRAP1", "Pab1" ,selected = 1)),
                             ), # sidebarPanel
                             
                             
                             mainPanel(
                                 h1("Fitness prediction"),
                                 h4("Results"),
                                 h5("Service under development -- not available"),) # mainPanel
                             
                    )
                ) # navbarPage
) # fluidPage


# Define server function
#------------------------------------------------------------------------------------------------------------------------------------------------

server <- function(input, output) {
    
    output$txtout1 <- renderText({
        
        #selected sequence
        paste("Sequence :",input$txt1)
        
    })
    
    output$txtout2 <- renderText({
        
        #selected parameters 
        paste("Selected parameters:",input$encoding1,",", input$protein1)
        
        
    })
    
    encoding_computation <- reactive({  
        
        #Format sequence to character
        sequence <-  as.character(input$txt1)
        
        #EXAMPLES of sequences
        
        #sequence_BRCA1 <- "DLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEF"
        #sequence_ERK2 <- "DALKSAGRALIRSPSLAKQSWGGGGRHRKLPENWTDTRETLLEGMLFSLKYLGMTLVEQPKGEELSAAAIKRIVATAKASGKKLQKVTLKVSPRGIILTDNLTNQLIENVSIYRISYCTADKMHDKVFAYIAQSQHNQSLECHAFLCTKRKMAQAVTLTVAQAFKVAFEFWQVSKEEKEKRDKASQEGGDVLGARQDCTPSLKSLVATGNLLDLEETAKAPLSTVSANTTNMDEVPRPQALSGSSVVWELDDGLDEAFSRLAQSRTNPQVLDTGLTAQDMHYAQCLSPVDWDKPDSSGTEQDDLFQ"
        #sequence_DLRAP1 <- "DALKSAGRALIRSPSLAKQSWGGGGRHRKLPENWTDTRETLLEGMLFSLKYLGMTLVEQPKGEELSAAAIKRIVATAKASGKKLQKVTLKVSPRGIILTDNLTNQLIENVSIYRISYCTADKMHDKVFAYIAQSQHNQSLECHAFLCTKRKMAQAVTLTVAQAFKVAFEFWQVSKEEKEKRDKASQEGGDVLGARQDCTPSLKSLVATGNLLDLEETAKAPLSTVSANTTNMDEVPRPQALSGSSVVWELDDGLDEAFSRLAQSRTNPQVLDTGLTAQDMHYAQCLSPVDWDKPDSSGTEQDDLFQ"
        #sequence_Pab1 <- "ANLHPDIDNKALYDTFSVFGDILSSKIATDENGKSKGFGFVHFEEEGAAKEAIDALNGMLLNGQEIY"
        
        #Encode sequence
        m <- input$encoding1
        encoder_matrix <- read.table(file = "z_scales.txt", row.names = 1, header =TRUE)
        encoded_sequence <- encode_peptide(sequence,encoder_matrix)
        
        #Load model
        model_protein <- input$protein1
        model_protein <- paste("elasticnet_",model_protein,".RData", sep="")
        model <- load(model_protein)
        
        #Compute prediction on the sequence input
        Y_pred_ElasticNet <- predict(fit.elasticnet, newx = encoded_sequence, s = s_, alpha=alpha_)
        
        #print output
        return(Y_pred_ElasticNet)
    })
    
    output$txtout3 <- renderText({
        
        encoding_computation()
        
    })
    
} # server




# Create Shiny object
shinyApp(ui = ui, server = server)
