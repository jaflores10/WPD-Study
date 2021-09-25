# fillmap code
fillmap<-function(map, figtitle, y , n.col, bk="e", cuts,legendtxt="",
                  leg.loc="bottomright",leg.cex=1.5,main.cex=1.5,main.line=-2,leg.horiz=F,map.lty=1){
  
  if(bk=="q"){if (min(y)<min(y) | max(y)>max(y)){
    print("The minimum or maximum	values of y fall outside of those for y")
  } else {p <- seq(0,1, length=n.col+1)
  br <- round(quantile(y, probs=p),2)}}
  if(bk=="e"){if (min(y)<min(y) | max(y)>max(y)){
    print("The minimum or maximum values of y fall outside of those for y")
  } else {br <- round(seq(min(y), max(y), length=n.col+1),6)}}
  if(bk=="c"){if (length(cuts)!= (n.col+1)) {cat("Cut off and color categories 
	do not match. ", "\n")
    break}  else {br <- cuts}  }
  
  # 0: dark 1: light light Current shading ranges from darkest to light gray white (to distinguish with lakes).
  shading<-gray(rev(0:(n.col-1)/(n.col-1)))
  #shading<-hsv(.6,alpha=0:(n.col-1)/(n.col-1))
  y.grp<-findInterval(y, vec=br, rightmost.closed = TRUE, all.inside = TRUE)
  y.shad<-shading[y.grp]
  
  plot(map,col=y.shad,axes=F, lty=map.lty)
  title(main=figtitle,cex.main=main.cex,line=main.line) 
  
  br<-round(br, 2)
  
  if (is.na(legendtxt[1])){print("No legend specifed")
  } else if (legendtxt[1]==""){
    leg.txt<-paste("[",br[n.col],",",br[n.col+1],"]",sep="")
    for(j in (n.col-1):1){ 
      leg.txt<-append(leg.txt,paste("[",br[j],",",br[j+1],")",sep="")) }
    leg.txt<-rev(leg.txt)
    legend(leg.loc,legend=leg.txt,fill=shading,cex=leg.cex,ncol=1,bty="n",
           horiz=leg.horiz)
  } else if (length(legendtxt) != n.col){cat("Length of lengendtxt must equal 
		n.col", "\n")
    break
  } else {leg.txt<-legendtxt
  legend(leg.loc,legend=leg.txt,fill=shading,cex=leg.cex,ncol=1,bty="n",
         horiz=leg.horiz)}
}

#fillmap2
fillmap2<-function(map, figtitle, y , leg.loc="beside", y.scl=NULL,
                   main.cex=1.5,main.line=0,map.lty=1,leg.rnd=0,
                   leg.cex=1){
  
  # 0: dark 1: light light Current shading ranges from darkest to light gray white (to distinguish with lakes).
  y.uq=sort(unique(c(y,y.scl)))
  cols<-viridis(length(y.uq),direction=-1)
  shading=y
  for (i in 1:length(y)){
    shading[i]<-cols[which(y.uq==y[i])]
  }
  
  par(mar=c(0,0,2,0))
  if (leg.loc=="beside"){
    layout(matrix(1:2,ncol=2),width=c(.8,.2))
  } else 
    if (leg.loc=="below"){
      layout(matrix(1:2,nrow=2),height=c(.6,.4))
    } else (print("leg.loc options are below or beside"))
  
  plot(map,col=shading,axes=F, lty=map.lty)
  title(main=figtitle,cex.main=main.cex,line=main.line) 
  
  par(mar=c(5, 4, 4, 2) + 0.1)
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  cols.5=cols[seq(1,length(y.uq),length.out=5)]
  lab.5=cols.5
  for (i in 1:5){lab.5[i]=y.uq[which(cols==cols.5[i])[1]]}
  lab.5=round(as.numeric(lab.5),leg.rnd)
  par(mar=c(0,0,0,0))
  if (leg.loc=="beside"){
    legend_image <- as.raster(matrix(cols, ncol=1))
    text(x=1.6, 
         y = seq(0,length(y.uq),length.out=5)/length(y.uq),
         labels = rev(lab.5), cex=leg.cex)
    rasterImage(legend_image, 0, 0, 1,1)
  } else{
    legend_image <- as.raster(matrix(cols, nrow=1))
    text(y=-0.25, 
         x = seq(0,length(y.uq),length.out=5)/(length(y.uq)*.5),
         labels = lab.5, cex=leg.cex)
    rasterImage(legend_image, 0, 0, 2,1)
  }
}

# Packages
library(INLA)
library(rgdal)
library(spdep)
library(maptools)
library(corrplot)
library(visdat)
library(viridis)

# Read shape data
NHtracts=readOGR("data/NHtract_map.shp")

# WPD arrest data
data=read.csv("data/wpd_arrests_race_tract1018_clean.csv")
nrow(data)

data1 = data[!(data$census_tract_code == 990100),]
dim(data1)
summary(data1)
length(unique(data1$census_tract_code))

length(unique(data1$census_tract_code))
length(unique(NHtracts$TRACTCE))

head(unique(data1$census_tract_code))
head(unique(NHtracts$TRACTCE))

length(unique(data1$census_tract_code))
length(unique(NHtracts$TRACTCE))

head(unique(data1$census_tract_code))
head(unique(NHtracts$TRACTCE))

# arrests as a percent of total arrests
data1$total_arrests_prct = data1$arrests_total/(data1$arrests_total + 0.1)
data1$arrests_B_prct = data1$arrests_B/(data1$arrests_total + 0.1)
data1$arrests_W_prct = data1$arrests_W/(data1$arrests_total + 0.1)

# arrests as a percent of the relevant population
data1$total_arrests_pop_prct = (data1$arrests_total)/(data1$ct_pop + 0.1)
data1$arrests_B_pop_prct = (data1$arrests_B)/(data1$ct_black + 0.1)
data1$arrests_W_pop_prct = (data1$arrests_W)/(data1$ct_white + 0.1)

# SIR of total arrests
rvalue_tot = sum(data1$arrests_total)/sum(data1$ct_pop)
pvalue_tot = data1$ct_pop + 0.1
evalue_tot = rvalue_tot * pvalue_tot
SIR_tot = data1$arrests_total/evalue_tot
length(SIR_tot)

# SIR of black arrests
rvalue_b = sum(data1$arrests_B)/sum(data1$ct_black)
pvalue_b = data1$ct_black + 0.1
evalue_b = rvalue_b * pvalue_b
SIR_b = data1$arrests_B/evalue_b
length(SIR_b)

# SIR of white arrest
rvalue_w = sum(data1$arrests_W)/sum(data1$ct_white)
pvalue_w = data1$ct_white + 0.1
evalue_w = rvalue_w * pvalue_w
SIR_w = data1$arrests_W/evalue_w
length(SIR_w)

# create SIR columns
data1$SIR_tot_pop = SIR_tot
data1$SIR_B = SIR_b
data1$SIR_W = SIR_w
summary(data1)


d.inla_1 = read.csv("data/INLAdata.csv")[,-1]

id = d.inla_1$id2
# Poisson Regression - Total arrests
ftot = d.inla_1$arrests_total ~ d.inla_1$black + d.inla_1$poverty + d.inla_1$educBachPlus + d.inla_1$secperctot + d.inla_1$age1824.perc + f(id, model="iid", param=c(2,1))
restot=inla(ftot1, data=d.inla_1, family='poisson', E=eTot, control.compute=list(dic=TRUE, waic=TRUE))
# Poisson Regression - Black arrests
ftot_b = d.inla_1$arrests_B ~ d.inla_1$black + d.inla_1$poverty + d.inla_1$educBachPlus + d.inla_1$secperctot + d.inla_1$age1824.perc + f(id, model="iid", param=c(2,1))
restot_b = inla(ftot_b, data=d.inla_1, family='poisson', E=eBlk, control.compute=list(dic=TRUE, waic=TRUE))
# Poisson Regression - White arrests
ftot_w = d.inla_1$arrests_W ~ d.inla_1$black + d.inla_1$poverty + d.inla_1$educBachPlus + d.inla_1$secperctot + d.inla_1$age1824.perc + f(id, model="iid", param=c(2,1))
restot_w = inla(ftot_w, data=d.inla_1, family='poisson', E=eWht, control.compute=list(dic=TRUE, waic=TRUE))

# R Shiny App
library(shiny)
library(shinythemes)


#Define UI
ui <- shinyUI(fluidPage(
  theme = shinytheme("superhero"),
  titlePanel("2010-2018 Wilmington Police Department Data"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("Year", label="Year", min=2010, max=2018, value=2010, sep="",
                  animate=animationOptions(interval=500, loop=TRUE)),
      radioButtons("data", label="Data:", c("Total Arrests", "White Only Arrests", "Black Only Arrests")),
      radioButtons("adj", label="Data Adjustment:", c("None", "As a Percent of the Population", "As a Percent of Total Arrests",
                                                      "Standardized Incidence Ratio (SIR)", "Poisson Regression")),
      br()
    ),
    mainPanel(
      textOutput("text"),
      plotOutput("map"),
      tableOutput("table")
    )
  )
))

#plot the arrest data
# first output displays the relevant text above each map
server <- shinyServer(function(input,output){
  output$text <- renderText({
    if (input$adj=="None"){
      "No Adjustments specified. The displayed data are arrest counts for the selected data."
    } else
      if (input$adj=="Standardized Incidence Ration (SIR)"){
        "The SIR adjustment was applied. SIR is a method to adjust for tract population. Values greater than one suggest more observed arrests than expected."
      } else
        if (input$adj=="Poisson Regression"){
          "A Poisson regression for adjustment was applied. Higher, darker values indicate areas of increased risk and tract population is indirectly adjusted."
        } else
          if (input$adj=="As a Percent of the Population"){
            "The 'As a Percent of the Population' adjustment displays the selected arrest counts divided by the relevant population."
          } else {
            "The 'As a Percent of Total Arrest' adjust displays the selected arrests counts divided by the total arrest counts."
          }
  })
  # second output displays the relevant maps
  output$map <- renderPlot({
    if (input$data == "Total Arrests" & input$adj == "None"){
      MapData = data1$arrests_total[seq(input$Year-2009, dim(data1)[1],9)]
      MapDataScl = data1$arrests_total
      Caption = paste(input$year, input$data)
    } else
      if (input$data == "White Only Arrests" & input$adj == "None"){
        MapData = data1$arrests_W[seq(input$Year-2009, dim(data1)[1],9)]
        MapDataScl = data1$arrests_W
        Caption = paste(input$year, input$data)
      } else
        if (input$data == "Black Only Arrests" & input$adj == "None"){
          MapData = data1$arrests_B[seq(input$Year-2009, dim(data1)[1],9)]
          MapDataScl = data1$arrests_B
          Caption = paste(input$year, input$data)
        } else
          if (input$data =="Total Arrests" & input$adj == "As a Percent of the Population"){
            MapData = data1$total_arrests_pop_prct[seq(input$Year-2009, dim(data1)[1],9)]
            MapDataScl = data1$total_arrests_pop_prct
            Caption = paste(input$year, input$data)
          } else
            if (input$data == "White Only Arrests" & input$adj == "As a Percent of the Population"){
              MapData = data1$arrests_W_pop_prct[seq(input$Year-2009, dim(data1)[1],9)]
              MapDataScl = data1$arrests_W_pop_prct
              Caption = paste(input$year, input$data)
            } else
              if (input$data == "Black Only Arrests" & input$adj == "As a Percent of the Population"){
                MapData = data1$arrests_B_pop_prct[seq(input$Year-2009, dim(data1)[1],9)]
                MapDataScl = data1$arrests_B_pop_prct
                Caption = paste(input$year, input$data)
              } else
                if (input$data == "Total Arrests" & input$adj == "As a Percent of Total Arrests"){
                  MapData = data1$total_arrests_prct[seq(input$Year-2009, dim(data1)[1],9)]
                  MapDataScl = data1$total_arrests_prct
                  Caption = paste(input$year, input$data)
                } else
                  if (input$data == "White Only Arrests" & input$adj == "As a Percent of Total Arrests"){
                    MapData = data1$arrests_W_prct[seq(input$Year-2009, dim(data1)[1],9)]
                    MapDataScl = data1$arrests_W_prct
                    Caption = paste(input$year, input$data)
                  } else
                    if (input$data == "Black Only Arrests" & input$adj == "As a Percent of Total Arrests"){
                      MapData = data1$arrests_B_prct[seq(input$Year-2009, dim(data1)[1],9)]
                      MapDataScl = data1$arrests_B_prct
                      Caption = paste(input$year, input$data)
                    } else
                      if (input$data == "Total Arrests" & input$adj == "Standardized Incidence Ratio (SIR)"){
                        MapData = data1$SIR_tot_pop[seq(input$Year-2009, dim(data1)[1],9)]
                        MapDataScl = data1$SIR_tot_pop
                        Caption = paste(input$year, input$data)
                      } else
                        if (input$data == "White Only Arrests" & input$adj == "Standardized Incidence Ratio (SIR)"){
                          MapData = data1$SIR_W[seq(input$Year-2009, dim(data1)[1],9)]
                          MapDataScl = data1$SIR_W
                          Caption = paste(input$year, input$data)
                        } else
                          if (input$data == "Black Only Arrests" & input$adj == "Standardized Incidence Ratio (SIR)"){
                            MapData = data1$SIR_B[seq(input$Year-2009, dim(data1)[1],9)]
                            MapDataScl = data1$SIR_B
                            Caption = paste(input$year, input$data)
                          } else
                            if (input$data == "Total Arrests" & input$adj == "Poisson Regression"){
                              MapData = restot1$summary.random$id$mean[seq(input$Year-2009, dim(data1)[1],9)]
                              MapDataScl = restot1$summary.random$id$mean
                              Caption = paste(input$year, input$data)
                            } else
                              if (input$data == "White Only Arrests" & input$adj == "Poisson Regression"){
                                MapData = restot_w$summary.random$id2$mean[seq(input$Year-2009, dim(data1)[1],9)]
                                MapDataScl = restot_w$summary.random$id$mean
                                Caption = paste(input$year, input$data)
                              } else
                                if (input$data == "Black Only Arrests" & input$adj == "Poisson Regression"){
                                  MapData = restot_b$summary.random$id2$mean[seq(input$Year-2009, dim(data1)[1],9)]
                                  MapDataScl = restot_b$summary.random$id2$mean
                                  Caption = paste(input$year, input$data)
                                }
    fillmap2(NHtracts, Caption, MapData, map.lty=0, leg.loc="beside", y.scl=MapDataScl, leg.rnd=3)
  })
})

shinyApp(ui=ui, server=server)