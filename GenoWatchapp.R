# This is APP for monitor the genotyping quality and intergrating with pedigree information. 
library(shiny)
library(ggplot2) 
library(scales)
library(DT)
library(ggpubr)
library(purrr) 
library(dplyr) 
library(vtable)
library(gtsummary)
library(gt)
library(data.table)
library(viridis)
library("gridExtra")
library(stringr)
library(plotly)
library(tidyr) 
library(GGally)
ui <- fluidPage(
  title="GenoWatch: Real Time Genotype Quality Monitor",
  titlePanel(
    h1("GenoWatch: Real Time Genotype Quality Monitor", align = "center")
    ),

  sidebarLayout(

    sidebarPanel(
  selectInput("sampleName", "FarmLine",choices = c("BH-12",
"L-1",
"L-2",
"L-3",
"L-4",
"L-5",
"L-6",
"L-7",
"L-8",
"L-9",
"L-10",
"L-11"
)),
   tags$hr(),
  radioButtons("rb", "Choose one (QC):",
               choiceNames = list(
                 "All Samples",
				 "Good Call Samples"
               ),
               choiceValues = list(
                 "all", "good"
               ))
  ,textInput("regIDupdate", "Registration ID:", ""),
   textAreaInput("commentBox", "Enter your QC comment:", ""),
   actionButton("submitBtn1", "Submit QC comment"),actionButton("resetButton", "Reset QC"),
   tags$hr(),
  # actionButton("updateBtnPT", "Update PT"),
    radioButtons("donePT", "Choose one (PT) :",
               choiceNames = list(
                 "Not Yet",
				 "Done"
               ),
               choiceValues = list(
                 "FALSE", "TRUE"
               ))
  ,
    textInput("ptIDupdate", "PT run ID:", ""),
   textAreaInput("PTcommentBox", "Enter your PT comment:", ""),
   actionButton("PTsubmitBtn", "Submit PT comment"),actionButton("PTresetButton", "Reset PT")			   
	  ,width = 3),
    mainPanel(
            tabsetPanel(
        tabPanel("summaryByRegID", align="center",gt_output(outputId = "mytable4"),column(12,align="left",verbatimTextOutput("savedComments")),column(12,align="left",verbatimTextOutput("newComments")),plotOutput("hetgraph1", "auto", "auto"),plotOutput("hetgraph2", "auto", "auto"),plotOutput("hetgraph3", "auto", "auto"),plotOutput("hetgraph4", "auto", "auto"),plotOutput("bargraph1", "auto", "auto"),tableOutput("mytableTop")
		),

        tabPanel("summaryQC", DT::dataTableOutput("mytable"),plotOutput("graph1", "auto", "auto"),plotOutput("QCgraph2", "auto", "auto")),
        tabPanel("summaryByPT",align="center",gt_output(outputId = "PTmytable1"),column(12,align="left",verbatimTextOutput("PTsavedComments")),column(12,align="left",verbatimTextOutput("PTnewComments")),plotOutput("PThetgraph1", "auto", "auto"),plotOutput("PTbargraph1", "auto", "auto")),
		tabPanel("summaryByMG", plotOutput("barMGgraph1"),plotlyOutput('plotMGnoG'), 
		    HTML("<h4 style='text-align: center;'>Table 1. The Sample Quality Distribution by MG</h4>"),
		column(12,align="center",tableOutput("mytableMG")),tableOutput("mytableMGtest"),gt_output("summaryByMGtable1")),
        tabPanel("summaryStatistics", gt_output("mytable2"), gt_output(outputId = "mytable3"),plotOutput("summaryGraph1", "auto", "auto")),
		tabPanel("about", imageOutput("demo_image",inline = TRUE),tableOutput("mytableQCrule"))#,
      )
    )
  )
)
server <- function(input, output, session) {
  sampleset<- fread(file="Line_Summary.txt",head=TRUE,sep="\t")
    fileNameSplit <- do.call(rbind, strsplit(sampleset$regID, split = "_"))

	sampleset$Farm=fileNameSplit[,1]
	sampleset$Line=fileNameSplit[,2]
	sampleset$MG=fileNameSplit[,3]
	sampleset$Hatch=fileNameSplit[,4]
	sampleset$DateT=fileNameSplit[,9]   # add the date. 
	
	sampleset$posX=substring(sampleset$Well.on.DNA.plate,2,3)   # add the date. 
	sampleset$posY=factor(substring(sampleset$Well.on.DNA.plate,1,1))   # add the date. 
	sampleset$farmLinePlate=paste(sampleset$farmline,sampleset$Sample.Plate,sep=":")
	sampleset$Sample.Plate=as.factor(sampleset$Sample.Plate)
	sampleset$MGHatch=paste(sampleset$MG,sampleset$Hatch,sep="-")
	sampleset$lowCall=ifelse(sampleset$Call.Rate<0.85,"true","false")

    # Specify the breaks for grouping
    breaks <- c(0.0, 5.0, 10.0,20.0,40.0, Inf)  # Inf represents infinity
    # Use cut() to create groups
    sampleset$DNAConcentrationGroups <- cut(sampleset$dnaConcentration, breaks, labels = c("[0,5]", "[5,10]", "[10,20]", "[20,40]", ">40"))
    sampleset$DNAConcentrationGroups <- factor(sampleset$DNAConcentrationGroups) 
	sampleset$indexID <- paste(sampleset$Sample.ID,sampleset$regID,sep=":")
# debug

  vipSNPabout<- read.table(file="GQ_about.txt",head=TRUE,sep="\t")
  qcRule<- read.table(file="QC_rule.txt",head=TRUE,sep="\t")

    pedigreeSet<- fread(file="Pedigree.txt",head=TRUE,sep="\t")
	

# PT data.
Data_PT= fread("PT_Summary.txt", header = TRUE ,sep="\t")
Data_PT$indexID <- paste(Data_PT$Sample.ID,Data_PT$regID,sep=":")
Data_PT=Data_PT[,c(-1,-2)]
# comments part
  observeEvent(input$resetButton, {
    updateTextAreaInput(session, "commentBox", value = "")
    updateTextInput(session, "regIDupdate", value = "")
	comments$data <- ""
  })
  
  comments <- reactiveValues(data = character(0))
  allcomments <-  fread(file="comments.txt",head=FALSE,sep="\t")
  observeEvent(input$submitBtn1, {
    newComment <- isolate(input$commentBox)
    if (nchar(newComment) > 0 && input$regIDupdate !="") {
      date1 <- format(Sys.time(), "%Y-%m-%d %H:%M:%S") 
      newComment <- paste(date1,newComment,sep=" comment is ")
      comments$data <- c(comments$data, newComment)
	  newComment2 <- paste(input$regIDupdate,newComment,sep="\t")
      write(newComment2, "comments.txt", append = TRUE)
    }
  })

  output$newComments <- renderPrint({
    cat(comments$data, sep = "\n")
  })

goodcomments <- reactive({
  if (input$regIDupdate !="") {
  oldcomments=subset(allcomments,V1 == str_trim(input$regIDupdate))
  }
 return(oldcomments)		
      })

  output$savedComments <- renderPrint({  
  if (input$regIDupdate !="") {
  oldcomments=goodcomments()
  if (dim(oldcomments)[1] != 0){
    cat(oldcomments$V2, sep = "\n")
  }
  }
  })
  
  observeEvent(input$PTresetButton, {
    updateTextAreaInput(session, "PTcommentBox", value = "")
    updateTextInput(session, "ptIDupdate", value = "")
	PTcomments$data <- ""
  })
  
  PTcomments <- reactiveValues(data = character(0))
  PTallcomments <-  fread(file="PTcomments.txt",head=FALSE,sep="\t")
  observeEvent(input$PTsubmitBtn, {
    PTnewComment <- isolate(input$PTcommentBox)
    if (nchar(PTnewComment) > 0 && input$ptIDupdate !="" && input$donePT=="TRUE") {
      date1 <- format(Sys.time(), "%Y-%m-%d %H:%M:%S") 
      PTnewComment <- paste(date1,PTnewComment,sep=" comment is ")
      PTcomments$data <- c(PTcomments$data, PTnewComment)
	  PTnewComment2 <- paste(input$ptIDupdate,PTnewComment,sep="\t")
      write(PTnewComment2, "PTcomments.txt", append = TRUE)
    }
  })
  # Display saved PTcomments
  output$PTnewComments <- renderPrint({
    cat(PTcomments$data, sep = "\n")
  })

PTgoodcomments <- reactive({
  if (input$ptIDupdate !="" && input$donePT=="TRUE") {
  PToldcomments=subset(PTallcomments,V1 == str_trim(input$ptIDupdate))
  }
 return(PToldcomments)		
      })

  output$PTsavedComments <- renderPrint({  
  if (input$ptIDupdate !="" && input$donePT=="TRUE" ) {
  PToldcomments=PTgoodcomments()
  if (dim(PToldcomments)[1] != 0){
    cat(PToldcomments$V2, sep = "\n")
  }
  }
  })
  
# end PT comments part 

  goodSample <- reactive({
       namePlot=input$sampleName
	   goodSampleplot1=subset(sampleset, farmline %in% namePlot)
	   goodSampleplot1<- goodSampleplot1 %>% group_by(Sample.ID) %>% mutate(numCount=n())
	   if(input$rb=="good") {
	    goodSampleplot=subset(goodSampleplot1,Call.Rate >=0.85)
		} else {
		goodSampleplot=goodSampleplot1
		}
		if(input$regIDupdate !="") {
	    goodSampleplot=subset(goodSampleplot,regID == str_trim(input$regIDupdate))
		} else {
		goodSampleplot=goodSampleplot
		}
      return(goodSampleplot)		
      })

  PTgoodSample <- reactive({
       namePlot=input$sampleName
	   goodSampleplot=subset(Data_PT, farmline %in% namePlot)
	   goodSampleplot$Difference[which(is.na(goodSampleplot$Difference))]<-"noError"
	   goodSampleplot$Reason[which(is.na(goodSampleplot$Reason))]<-"noError"
	   goodSampleplot$AssignGender[which(is.na(goodSampleplot$AssignGender))]<-"notAssign"

		if(input$ptIDupdate !="") {
	    goodSampleplot=subset(goodSampleplot,ptID == str_trim(input$ptIDupdate))
		goodSampleplot=left_join(goodSampleplot,sampleset)
		} else {
		goodSampleplot=goodSampleplot
		goodSampleplot=left_join(goodSampleplot,sampleset)
		}
      return(goodSampleplot)		
      })
	  
	  
  calcHeight <- reactive({
  namePlot=input$sampleName
  numSample <- dim(goodSample())[1]
  newheight=500

  if(numSample>=1600) newheight=800
  if(numSample>=2600) newheight=1000
    # c(100 * length(input[["scoresTable_rows_selected"]]))
	return(newheight)
  })

  PTcalcHeight <- reactive({
  namePlot=input$sampleName
  numSample <- dim(PTgoodSample())[1]
  newheight=500

  if(numSample>=1600) newheight=800
  if(numSample>=2600) newheight=1000
    # c(100 * length(input[["scoresTable_rows_selected"]]))
	return(newheight)
  })
  
  pedSample <- reactive({
		namePlot=input$sampleName

		subsampleset=select(goodSample(),-c("MG","Hatch","MGHatch","farmline"))
 	    # goodPedplot=subset(pedigreeSet, farmline %in% namePlot) 
 	    subpedigreeSet=subset(pedigreeSet, farmline %in% namePlot) 
	    subData_PT=select(PTgoodSample(),c("Geno.Sire","Geno.Dam","Difference","Reason","sire","dam","ptID","AssignGender","indexID"))		
      if(input$regIDupdate =="" && input$ptIDupdate =="") {

        uniqsubsampleset <- subsampleset %>%
                              group_by(Sample.ID) %>%
                               slice(which.max(Call.Rate)) %>%
                               ungroup()

        submerger1 <- left_join(uniqsubsampleset,subData_PT,by="indexID")
		submerger2 <- left_join(subpedigreeSet,submerger1,by="Sample.ID")
		submerger2$Geno="G"
        submerger2$Geno[is.na(submerger2$indexID)]="No_G"		
	    submerger2$MGHatch <- paste(submerger2$MG,submerger2$Hatch,sep="-")
		
		return(submerger2)
      }
	  
	  })

  
  output$plotMGnoG <- renderPlotly({
  
  if(input$regIDupdate =="" && input$ptIDupdate =="" ) {
    output2=pedSample() %>% group_by(MGHatch,Geno) %>% tally()
    output2$Geno=as.factor(output2$Geno)
    output2$MGHatch=as.factor(output2$MGHatch)		
	outputwide=spread(output2, Geno, n)

	fig <- plot_ly(outputwide,x = ~MGHatch, y=~G, type = 'bar', name = 'G') %>% add_trace(y=~No_G,name = 'No_G') %>% 
       layout(yaxis = list(title = 'Count'),
              barmode = 'stack',title="Number of Genotyping vs. Non-Genotyping Animals")  %>% config(displayModeBar = FALSE)							   
# fig
}
})


  output$mytableMG <-  renderTable({
        goodSampleplot=pedSample()
		if(input$regIDupdate =="" && input$ptIDupdate =="" ) {
        summary_result <- goodSampleplot %>%
        group_by(MGHatch) %>%
         summarise(

           NumTotal = n(),
           NumNoGeno = paste(length(which(Geno=="No_G"))," (",percent(length(which(Geno=="No_G")) / n(), accuracy=0.01),")",sep=""),
           NumLowCall = paste(length(which(lowCall=="true"))," (",percent(length(which(lowCall=="true")) / n(), accuracy=0.01),")",sep=""),
           NumHetAb = paste(length(which(HetAbnormal=="true"))," (",percent(length(which(HetAbnormal=="true")) / n(), accuracy=0.01),")",sep=""),
           NumPTerror = paste(length(which(Reason=="PT"))," (",percent(length(which(Reason=="PT")) / n(), accuracy=0.01),")",sep=""),
           NumAssignGender = paste(length(which(AssignGender=="M"))+length(which(AssignGender=="F"))," (",percent((length(which(AssignGender=="M"))+length(which(AssignGender=="F"))) / n(), accuracy=0.01),")",sep=""),
           NumRerun=paste(sum(numCount>=2,na.rm = TRUE)," (",percent(sum(numCount>=2,na.rm = TRUE) / n(), accuracy=0.01),")",sep="")
               )
  summary_result
		} 
	}, type = "html", bordered = TRUE, striped = TRUE, align = "c")

  output$mytableQCrule <- renderTable({
    qcRule
  }, caption = "Table 1. The Line, Farm and Rule", caption.placement = "top",type = "html", bordered = TRUE, striped = TRUE, align = "c")
  
  output$mytable2 <- render_gt({
    sampleset[,c("farmline","Sex.Estimate","SexPredict")]  %>% tbl_summary(by = farmline) %>%
        # CONVERT TO A {gt} TABLE! VERY IMPORTANT STEP!
        as_gt() %>%
        tab_header(md("**Table 1. Gender Summary **"))	%>%
        gt::tab_source_note(gt::md("*The summary is based on the 4 recent MG *"))
  })
  
  output$mytable3 <-  render_gt({

    sampleset[,c("farmline","Call.Rate","HetRate","AveReadDepth","dnaConcentration")]  %>% tbl_summary(statistic =  all_continuous() ~ "{mean} ({sd})",
	by = farmline,missing = "no") %>%
        # CONVERT TO A {gt} TABLE! VERY IMPORTANT STEP!
        as_gt() %>%
        tab_header(md("**Table 2. Quality Summary **"))	%>%
        gt::tab_source_note(gt::md("*The summary is based on the 4 recent MG *"))	
  })

 output$summaryGraph1 <- renderPlot({
		# if(input$regIDupdate =="") {

plot1 <- ggplot(sampleset, aes(x = farmline, y = Call.Rate, fill = farmline)) +
  geom_boxplot() +
  theme_minimal() +ggtitle("The boxplot for Call Rate by Farm Line")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20),axis.text.x = element_text(size = 15),,axis.text.y = element_text(size = 15)) # + geom_text(aes(label=lowCall),position = position_fill(vjust = -0.5))

plot2 <- ggplot(sampleset, aes(x = farmline, y = HetRate, fill = farmline)) +
  geom_boxplot() +
  theme_minimal() +ggtitle("The boxplot for Het Rate by Farm Line")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20),axis.text.x = element_text(size = 15),,axis.text.y = element_text(size = 15)) # + geom_text(aes(label=lowCall),position = position_fill(vjust = -0.5))
	
grid.arrange(plot1, plot2, nrow= 2)   # ncol
		# }
} , height = 1000) 

 
  output$mytable4 <-  render_gt({
        goodSampleplot=goodSample()
		if(input$regIDupdate !="") {
		regID = str_trim(input$regIDupdate)
		ouputMytable4Text=paste("*This data is only for Registration ID ( ",regID," ) from the text box*",sep="")

        goodSampleplot[,c("SexPredict","lowCall","HetAbnormal","GenoDuplicate","Call.Rate","HetRate","AveReadDepth","dnaConcentration")]  %>% tbl_summary(statistic =  all_continuous() ~ "{mean} ({sd})") %>%
        # CONVERT TO A {gt} TABLE! VERY IMPORTANT STEP!
        as_gt() %>%
        tab_header(md("**Table 1. One Registration Quality Summary**")) %>%
        gt::tab_source_note(gt::md(ouputMytable4Text))
        # gt::tab_source_note(gt::md("*This data is only for one Registration ID from the text box*"))
		} else {
		}
		
  })


  output$mytableTop <-  renderTable({
        goodSampleplot=goodSample()
		if(input$regIDupdate !="") {
        result_dplyr1 <- goodSampleplot %>%
        group_by(Sample.Plate) %>%
        summarise(LowCallNumber = length(which(lowCall=="true")),LowCallPercentage = percent(length(which(lowCall=="true")) / n(), accuracy=0.01))
        result_dplyr2 <- goodSampleplot %>%
        group_by(Sample.Plate) %>%
        summarise(HetAbnormalNumber = length(which(HetAbnormal=="true")),HetAbnormalPercentage = percent(length(which(HetAbnormal=="true")) / n(), accuracy=0.01))	
		result_dplyr3 <- goodSampleplot %>%
		group_by(Sample.Plate) %>% summarise(NumberOfRerun=sum(numCount>=2))		
        # summarise(Percentage = n() / nrow(sampleset) * 100)
        result_dplyr12=merge(result_dplyr1,result_dplyr2)
        result_dplyr=left_join(result_dplyr12,result_dplyr3)		
		result_dplyr
		} else {
		}
  }, caption = "Table 2. One Registration Quality Percentage Summary", caption.placement = "top",type = "html", bordered = TRUE, striped = TRUE, align = "c")
  # html_caption_str <- as.character(shiny::tags$b(style = "color: red", "a styled caption"))

  output$mytable <- DT::renderDataTable(server = FALSE, {
      goodSampleplot=goodSample()
	   # goodSampleplot=sampleset
	   # DT::datatable(vipSNPplot[,])
	  DT::datatable(
	    goodSampleplot[,c(1:12,14)],
	    extensions = c('Buttons'), options = list(
	      dom = 'Bfrtip',
          buttons = list(
          list(extend = "csv", text = "Download Current Page", filename = "page",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv", text = "Download Full Results", filename = "data",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          ))
	    ),filter = 'top',rownames = F	
	  )	   
  })

  output$graph1 <- renderPlot({
  # output$graph1 <- renderPlotly({
	  goodSampleplot=goodSample()
	namePlot=input$sampleName
if(namePlot=="BH-12" || namePlot=="DC-12" || namePlot=="HV-12" ){	
	  
plot1 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=1-Call.Rate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("Call.Rate distribution for samples") + xlab("") + ylab("1-Call.Rate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.15, linetype="dashed", color = "red") 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=HetRate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("HetRate distribution for samples") + xlab("") + ylab("HetRate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.40, linetype="dashed", color = "red") 

plot3 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=AveReadDepth,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("AveReadDepth distribution for samples") + xlab("Animal") + ylab("AveReadDepth") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

plot4 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=dnaConcentration,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("dnaConcentration distribution for samples") + xlab("Animal") + ylab("dnaConcentration") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

grid.arrange(plot1, plot2,plot3,plot4, nrow= 4)
}

if(namePlot=="BH-58" || namePlot=="DC-58" || namePlot=="HV-58" ){	
	  
plot1 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=1-Call.Rate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("Call.Rate distribution for samples") + xlab("") + ylab("1-Call.Rate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.15, linetype="dashed", color = "red") 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=HetRate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("HetRate distribution for samples") + xlab("") + ylab("HetRate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.45, linetype="dashed", color = "red") 

plot3 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=AveReadDepth,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("AveReadDepth distribution for samples") + xlab("Animal") + ylab("AveReadDepth") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

plot4 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=dnaConcentration,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("dnaConcentration distribution for samples") + xlab("Animal") + ylab("dnaConcentration") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

grid.arrange(plot1, plot2,plot3,plot4, nrow= 4)
}

if(namePlot=="BH-61" || namePlot=="DC-61" || namePlot=="HV-61" ){	
	  
plot1 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=1-Call.Rate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("Call.Rate distribution for samples") + xlab("") + ylab("1-Call.Rate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.15, linetype="dashed", color = "red") 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=HetRate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("HetRate distribution for samples") + xlab("") + ylab("HetRate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.50, linetype="dashed", color = "red") 

plot3 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=AveReadDepth,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("AveReadDepth distribution for samples") + xlab("Animal") + ylab("AveReadDepth") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

plot4 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=dnaConcentration,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("dnaConcentration distribution for samples") + xlab("Animal") + ylab("dnaConcentration") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

grid.arrange(plot1, plot2,plot3,plot4, nrow= 4)
}

if(namePlot=="HV-18" ){	
	  
plot1 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=1-Call.Rate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("Call.Rate distribution for samples") + xlab("") + ylab("1-Call.Rate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.15, linetype="dashed", color = "red") 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=HetRate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("HetRate distribution for samples") + xlab("") + ylab("HetRate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.43, linetype="dashed", color = "red") 

plot3 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=AveReadDepth,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("AveReadDepth distribution for samples") + xlab("Animal") + ylab("AveReadDepth") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

plot4 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=dnaConcentration,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("dnaConcentration distribution for samples") + xlab("Animal") + ylab("dnaConcentration") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

grid.arrange(plot1, plot2,plot3,plot4, nrow= 4)
}

if(namePlot=="HV-15" ){	
	  
plot1 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=1-Call.Rate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("Call.Rate distribution for samples") + xlab("") + ylab("1-Call.Rate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.15, linetype="dashed", color = "red") 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=HetRate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("HetRate distribution for samples") + xlab("") + ylab("HetRate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.45, linetype="dashed", color = "red") 

plot3 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=AveReadDepth,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("AveReadDepth distribution for samples") + xlab("Animal") + ylab("AveReadDepth") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

plot4 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=dnaConcentration,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("dnaConcentration distribution for samples") + xlab("Animal") + ylab("dnaConcentration") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

grid.arrange(plot1, plot2,plot3,plot4, nrow= 4)
}

if(namePlot=="BH-70" || namePlot=="DC-70" || namePlot=="HV-70" ){	
	  
plot1 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=1-Call.Rate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("Call.Rate distribution for samples") + xlab("") + ylab("1-Call.Rate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.15, linetype="dashed", color = "red") 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=HetRate,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("HetRate distribution for samples") + xlab("") + ylab("HetRate") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan")) + geom_hline(yintercept=0.41, linetype="dashed", color = "red") 

plot3 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=AveReadDepth,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("AveReadDepth distribution for samples") + xlab("Animal") + ylab("AveReadDepth") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

plot4 <- ggplot(goodSampleplot, aes(x=Sample.ID ,y=dnaConcentration,color=MG,shape =Sex.Estimate)) + geom_point(size = 2) + ggtitle("dnaConcentration distribution for samples") + xlab("Animal") + ylab("dnaConcentration") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))  + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan"))

grid.arrange(plot1, plot2,plot3,plot4, nrow= 4)
}


}, height = 1000 )

  output$QCgraph2 <- renderPlot({
	  goodSampleplot=goodSample()
 
plot1 <- ggpairs(goodSampleplot, columns = c("Call.Rate","HetRate","AveReadDepth","dnaConcentration"), aes(colour = lowCall)) + theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))

plot1
}, height = 600 )

  output$hetgraph1 <- renderPlot({
		if(input$regIDupdate !="") {
			  goodSampleplot=goodSample()
			  goodSampleplot$posY <- factor(goodSampleplot$posY, levels = c("H","G","F","E","D","C","B","A"))
ggplot(goodSampleplot, aes(x = reorder(posX,as.numeric(posX)) , y = posY, fill = Call.Rate)) + 
  geom_tile() + scale_fill_viridis() + 
  labs(title = paste("Call Rate By Plate for this RegID",sep=""), 
       x = "posX",
       y = "posY") + facet_wrap(~farmLinePlate, ncol = 4,scales = "free") + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))   # strip.text.x = element_text(size = 30)  https://statisticsglobe.com/change-font-size-of-ggplot2-facet-grid-labels-in-r
		}
		
 }, height = function(){calcHeight()})
  
# https://groups.google.com/g/shiny-discuss/c/dkZxTvfHOvo
  output$hetgraph2 <- renderPlot({
		if(input$regIDupdate !="") {
			  goodSampleplot=goodSample()
			  goodSampleplot$posY <- factor(goodSampleplot$posY, levels = c("H","G","F","E","D","C","B","A"))
ggplot(goodSampleplot, aes(x = reorder(posX,as.numeric(posX)) , y = posY, fill = HetRate)) +
  geom_tile() + scale_fill_viridis() + 
  labs(title = paste("HetRate By Plate for this RegID",sep=""), 
       x = "posX",
       y = "posY") + facet_wrap(~farmLinePlate, ncol = 4,scales = "free") + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))
}
}, height = function(){calcHeight()})
  

  output$hetgraph3 <- renderPlot({
		if(input$regIDupdate !="") {
			  goodSampleplot=goodSample()
			  goodSampleplot$posY <- factor(goodSampleplot$posY, levels = c("H","G","F","E","D","C","B","A"))
ggplot(goodSampleplot, aes(x = reorder(posX,as.numeric(posX)) , y = posY, fill = AveReadDepth)) +
  geom_tile() + scale_fill_viridis() + 
  labs(title = paste("AveReadDepth  By Plate for this RegID",sep=""), 
       x = "posX",
       y = "posY") + facet_wrap(~farmLinePlate, ncol = 4,scales = "free") + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))
	}
}, height = function(){calcHeight()})

  output$hetgraph4 <- renderPlot({
		if(input$regIDupdate !="") {
			  goodSampleplot=goodSample()
			  goodSampleplot$posY <- factor(goodSampleplot$posY, levels = c("H","G","F","E","D","C","B","A"))
ggplot(goodSampleplot, aes(x = reorder(posX,as.numeric(posX)) , y = posY, fill = DNAConcentrationGroups)) +
  # geom_tile() + scale_fill_viridis() + 
  geom_tile() + 
  labs(title = paste("dnaConcentration  By Plate for this RegID",sep=""), 
       x = "posX",
       y = "posY") + facet_wrap(~farmLinePlate, ncol = 4,scales = "free") + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))
	}
}, height = function(){calcHeight()})
 
  output$bargraph1 <- renderPlot({
		if(input$regIDupdate !="") {
			  goodSampleplot=goodSample()
plot1 <- ggplot(goodSampleplot, aes(x=Sample.Plate, fill=lowCall))+
  geom_bar(stat="count", width=0.7)+
  theme_minimal() +ggtitle("Number of Animals by Plate")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20)) # + geom_text(aes(label=lowCall),position = position_fill(vjust = -0.5))
 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.Plate, fill=HetAbnormal))+
  geom_bar(stat="count", width=0.7)+
  theme_minimal() +ggtitle("Number of Animals by Plate")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))	
grid.arrange(plot1, plot2, nrow= 2)   # ncol
		}
} , height = 700)   #


   output$barMGgraph1 <- renderPlot({
		if(input$regIDupdate =="") {
			  goodSampleplot=goodSample()
ggplot(goodSampleplot, aes(x=MG, fill=Hatch))+
  geom_bar(stat="count", width=0.7)+
  theme_minimal() +ggtitle("Number of Genotyping Animals by MG")+
 theme(plot.title = element_text(hjust = 0.5,size=20,face = "plain"),axis.title=element_text(size=15))
 }
 }) 
  
    output$demo_image <- renderImage({
    list(src = "copyright.png",
         width = 800,
         height = 200,
         alt = "A flower")
  }, deleteFile = FALSE)
  
# PT part

  output$PTmytable1 <-  render_gt({
        goodSampleplot=PTgoodSample()
		if(input$ptIDupdate !="" && input$donePT=="TRUE") {
		regID = str_trim(input$ptIDupdate)
		ouputMytable4Text=paste("*This data is only for PT ID ( ",regID," ) from the text box*",sep="")
#  goodSampleplot
        goodSampleplot[,c("Difference","Reason","AssignGender")]  %>% tbl_summary(statistic =  all_continuous() ~ "{mean} ({sd})") %>%
        as_gt() %>%
        tab_header(md("**Table 1. One PT run ID Quality Summary**")) %>%
        gt::tab_source_note(gt::md(ouputMytable4Text))
		} else {
		}
  })

  output$PThetgraph1 <- renderPlot({
		if(input$ptIDupdate !="" && input$donePT=="TRUE" ) {
			  goodSampleplot=PTgoodSample()
			  goodSampleplot$posY <- factor(goodSampleplot$posY, levels = c("H","G","F","E","D","C","B","A"))
ggplot(goodSampleplot, aes(x = reorder(posX,as.numeric(posX)) , y = posY, fill = Difference)) + 
  geom_tile() + 
  labs(title = paste("PT By Plate for this PT ID",sep=""), 
       x = "posX",
       y = "posY") + facet_wrap(~farmLinePlate, ncol = 4,scales = "free") + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))   # strip.text.x = element_text(size = 30)  https://statisticsglobe.com/change-font-size-of-ggplot2-facet-grid-labels-in-r
		}
		
 }, height = function(){PTcalcHeight()})
 
 output$PTbargraph1 <- renderPlot({
		if(input$ptIDupdate !="" && input$donePT=="TRUE") {
			  goodSampleplot=PTgoodSample()
plot1 <- ggplot(goodSampleplot, aes(x=Sample.Plate, fill=Difference))+
  geom_bar(stat="count", width=0.7)+
  theme_minimal() +ggtitle("Number of Animals by Plate")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20)) # + geom_text(aes(label=lowCall),position = position_fill(vjust = -0.5))
 

plot2 <- ggplot(goodSampleplot, aes(x=Sample.Plate, fill=Reason))+
  geom_bar(stat="count", width=0.7)+
  theme_minimal() +ggtitle("Number of Animals by Plate")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(plot.title = element_text(hjust = 0.5,size=20),axis.title=element_text(size=20))	
grid.arrange(plot1, plot2, nrow= 2)   # ncol
		}
} , height = 700) 


  
}

shinyApp(ui, server)
