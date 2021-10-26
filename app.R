library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(h4("Algal and environmental parameters")),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input on sliders # default species is F serratus
      sliderInput(inputId = "a.intercept",
                  label = "allometry intercept",
                  min = -5,
                  max = 1,
                  step = 0.01,
                  value = 0.21),
      sliderInput(inputId = "b.slope",
                  label = "allometry slope",
                  min = 1,
                  step = 0.01,
                  max = 4.6,
                  value = 2.85),
      sliderInput(inputId = "Pmax",
                  label = "Max photosyn rate",
                  min = 1,
                  max = 15,
                  step = 0.1,
                  value = 7),
      sliderInput(inputId = "Pmax.slope",
                  label = "rate of change in photosynthesis in air with water content",
                  min = 0,
                  max = 0.02,
                  step = 0.001,
                  value = 0.015),
      sliderInput(inputId = "Pmax.intercept",
                  label = "intecept for photosynthesis in air with water content",
                  min = -1,
                  max = 1,
                  step = 0.01,
                  value = -0.236),
      sliderInput(inputId = "alpha",
                  label = "initial slope PI curve",
                  min = 0.01,
                  max = 1,
                  step = 0.01,
                  value = 0.05),
      sliderInput(inputId = "n.fronds",
                  label = "frond density per m",
                  min = 1,
                  max = 1000,
                  value = 100),
      sliderInput(inputId = "f.mean",
                  label = "mean frond height",
                  min = 0.01,
                  max = 5,
                  step = 0.01,
                  value = 0.25),
      sliderInput(inputId = "f.sd",
                  label = "sd frond height",
                  min = 0,
                  max = 1,
                  step = 0.01,
                  value = 0),
      sliderInput(inputId = "TSM",
                  label = "thallus specific mass",
                  min = 50,
                  max = 1000,
                  step = 1,
                  value = 168),
      sliderInput(inputId = "resp",
                  label = "respiration rate",
                  min = 0,
                  max = 0.05,
                  step = 0.001,
                  value = 0.002),
      sliderInput(inputId = "recovery.threshold",
                  label = "threshold for desiccation recovery",
                  min = 0,
                  max = 60,
                  step = 5,
                  value = 40),
      sliderInput(inputId = "shore.level",
                  label = "height of algae relative to LAT",
                  min = -5,
                  max = 5,
                  step = 0.1,
                  value = 0),
      sliderInput(inputId = "M2.amplitude",
                  label = "M2 tidal amplitude (m)",
                  min = 0,
                  max = 5,
                  step = 0.1,
                  value = 0),
      sliderInput(inputId = "S2.amplitude",
                  label = "S2 tidal amplitude (m)",
                  min = 0,
                  max = 2.5,
                  step = 0.1,
                  value = 0),
      sliderInput(inputId = "S2.phase",
                  label = "S2 tidal phase lag (degrees)",
                  min = 0,
                  max = 360,
                  step = 10,
                  value = 0),
      sliderInput(inputId = "daylength",
                  label = "daylength",
                  min = 1,
                  max = 24,
                  step = 0.5,
                  value = 12),
      sliderInput(inputId = "Isurf",
                  label = "noon irradiance (PAR)",
                  min = 1,
                  max = 2500,
                  step = 100,
                  value = 1600),
      sliderInput(inputId = "lamda",
                  label = "water attenuation coefficient",
                  min = 0.01,
                  max = 3,
                  step = 0.01,
                  value = 0.7),
      sliderInput(inputId = "kappa",
                  label = "canopy attenuation coefficient",
                  min = 0.01,
                  max = 3,
                  step = 0.01,
                  value = 0.5),
      sliderInput(inputId = "D.rate",
                  label = "desiccation rate",
                  min = 0,
                  max = 2*10^-4,
                  step = 10^-5,
                  value = 1.17*10^-4),
      sliderInput(inputId = "Days",
                  label = "days simulated",
                  min = 1,
                  max = 14,
                  step = 1,
                  value = 1)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      fluidRow( 
        column(3,(tableOutput("TAIcaption"))),
     #   h3(textOutput("Int.caption")),
     # Button
     column(3,downloadButton("downloadData", "Download model output"))),
      
      # Output: Timeseries of variables ----
      plotOutput(outputId = "PPlot"),
      plotOutput(outputId = "TPlot"),
      plotOutput(outputId = "PhPlot")
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  library(REdaS)
  
  
  
  
  ##Functions for irradience to vary with time
  Ivary <- function(hours,halfdaylong) {
    Itime <- (1/2*halfdaylong)*(1+cos(pi*(hours-12)/(halfdaylong)))
    Itime<-Itime/halfdaylong
    if (abs(12-hours)>halfdaylong) {Itime<-0}
    return(Itime)
  }
  ###convert continous timestep to hours
  TOD <- function(timestep) {
    hours<-24*(timestep/96-floor(timestep/96))
    return(hours)
  }
 
  output$PPlot <- renderPlot({
    
    ## simulation run duration
    timesteps.max<-96*input$Days
    data.view<-matrix(1:timesteps.max,nrow=timesteps.max,ncol=5) #to generate output for plots (96 is 24 hours)
    #set up frond heights
    frond<-c(rlnorm(1:input$n.fronds, meanlog = log(input$f.mean),sdlog=input$f.sd),rep(0,1000-input$n.fronds)) ##recommended set up at 1000 elements, even if some 0
    frond<-round(frond,2) #round to 1 cm increments

    # set up matrices to help allocate biomass to discrete depth layers
    mat1 <- matrix(rep(0,1000),nrow=1000,ncol=1000) ##maximum 1000 fronds, 1000 depth segments
    ht<-c(999:0)
    ht<-ht/100
    mat1[,1:1000]<-ht # matrix 1000 columns all counting down to 0 from 9.99
    mat2<-matrix(frond,nrow=1000,ncol=1000,byrow=TRUE) #will repeat if length(frond) < 1000
    mat3<-mat2-mat1
    mat3<-1*(mat3>=0) #is 1 at all depths where the frond is there
    
    # equation relatine frond area to length
    a.intercept<-input$a.intercept #F serratus is set up when using the default parameters 
    b.slope<-input$b.slope
    area<-exp(a.intercept+b.slope*log(ht)) #area of fronds
    
    
    area2<-c(area[2:1000],0)
    area.slice<-area-area2
    area.slice.mat<-mat1
    area.slice.mat[,1:1000]<-area.slice
    area.slice.mat<-area.slice.mat*mat3 #gives area of frond in each slice of depth row 0 is 9.99 m of seabed.
    
    canopy<-rowSums(area.slice.mat) #total canopy in each slice above seabed
    TAI<-sum(canopy) #(Thallus area index) area fronds per square metre of seabed
   
    #shore level above LAT
    shore.level<- input$shore.level # defined by user, maximum (range- shore level) currently 9.99 m 
    water.content<-100
    min.water.content<-100
    #timer
    timer<- 1
    
    # Tidal constants
    M2.speed <-	28.984
    S2.speed <-	30
    # N2.speed <-		28.44 #unused tidal constants
    # K1.speed <-		15.041
    # O1.speed <-		13.943
    # M4.speed <-		57.968
    # M6.speed <-		86.952
    # S4.speed <-		60
    # MS4.speed <-		58.984
    M2.phase <-	0
    S2.phase <-	input$S2.phase
    # N2.phase <-		0 #unused tidal constants
    # K1.phase <-		0
    # O1.phase <-		0
    # M4.phase <-		0
    # M6.phase <-		0
    # S4.phase <-		0
    # MS4.phase <-		0
    M2.amplitude <-	input$M2.amplitude #is actually half range as cos function has range -1 to 1
    S2.amplitude <-	input$S2.amplitude # is actually half range as cos function has range -1 to 1
    #N2.amplitude <-		0 #unused tidal constants
    #K1.amplitude <-		0
    #O1.amplitude <-		0
    #M4.amplitude <-		0
    #M6.amplitude <-		0
    #S4.amplitude <-		0
    #MS4.amplitude <-		0
    #range<-2*(M2.amplitude + S2.amplitude + N2.amplitude + K1.amplitude +  O1.amplitude +  M4.amplitude +  M6.amplitude + S4.amplitude +  MS4.amplitude) ## full constants for tidal range
    range<-2*(M2.amplitude + S2.amplitude) ## M2 and S2 components only
    d.rate<-input$D.rate # desiccation rate before any modifications due to layering or weather
    #Light parameters (max surface flux), attenuation rates (canopy and water) and threshold for tolerating desiccation with full recovery
    Isurf<-input$Isurf
    lambda<-input$lamda
    kappa<-input$kappa
    recovery.threshold<-input$recovery.threshold
    #P-I curve parameters
    alpha<-input$alpha #parameters may vary with e.g., temp, user can input values.
    Pmax<-input$Pmax
    Pmax.R<-Pmax
    R<-input$resp #respiration rate (may vary with e.g., temp, user can input values)
    Tsm<-input$TSM #Thallus specific mass, converts area to dry mass as respiration is mass based (per gram dry mass)
    daylength<-input$daylength #(hours) if input is 24, light is set as a constant = Isurf
    d.slope<-input$Pmax.slope
    d.intercept<-input$Pmax.intercept
    em.toggle<-0 #just in case dirty background values for these variables exist
    em.time<-0 #just in case dirty background values for these variables exist
    sub.time<-0 #just in case dirty background values for these variables exist
    
    
    for(timer in 1:timesteps.max) {
      #generate tidal heights above LAT and area of canopy in each layer of water
      tide.ht.M2<-M2.amplitude*cos(deg2rad((timer/4)*M2.speed-M2.phase)) # 0 to range
      tide.ht.S2<-S2.amplitude*cos(deg2rad((timer/4)*S2.speed-S2.phase))
      tide.ht<-tide.ht.M2+tide.ht.S2+range/2
      #light by depth - first work out where surface is and canopy at surface
      wd<-round((tide.ht- shore.level),digits=2)
      if( wd < 0) {wd<-0} else {wd} #water depth is calculated as tidal.ht (0 to max) - (shore.level above minimum tide): set to zero if negative.
      working.ht<-wd-ht #gives depths below surface (as positive values)
      working.ht2<-1*(working.ht<=0) # to gather any surface canopy
      working.ht3<-1*(working.ht>=0) # submerged canopy
      working.canopy<-working.ht2*canopy # surface canopy layers
      surface.canopy<-sum(working.canopy) #scalar of canopy area in surface layer
      working.canopy2<-working.ht3*canopy #area in submerged layers except collation of surface
      working.canopy2[which(working.ht==0)]<-surface.canopy #now includes amount in surface layer
      cummulative.canopy<-cumsum(working.canopy2)
      
      # desiccation of fronds if tide has uncovered them (water depth = 0)
      if(wd==0) {em.toggle<-1} else {em.toggle<-2}
      if(em.toggle==1) {em.time<-em.time+1}
      if(em.toggle==1) {sub.time<-0}
      if(em.toggle==2) {sub.time<-sub.time+1}
      if(em.toggle==2) {em.time<-0}
      if (TAI>1) {dess.rate<-d.rate/TAI} #alters effective desiccation rate to accound for fronds sheltering each other
      if (TAI<=1) {dess.rate<-d.rate}
      if (wd==0) {water.content<-100*exp(-dess.rate*em.time*15*60)} 
      if (water.content<min.water.content) {min.water.content<-water.content}
      
      # light part
      zed<-(working.ht+0.005)*working.ht3 #depth is middle of layer for water attenuation, using working.ht3 means that the modifier only applies to submerged canopy
      cummulative.canopy2<-cummulative.canopy
      cummulative.canopy2[which(working.ht==0)]<-cummulative.canopy2[which(working.ht==0)]-1 #allowing unshaded layer at surface
      if (cummulative.canopy2[which(working.ht==0)]<0) {cummulative.canopy2[which(working.ht==0)]<-0} #avoids negative canopy at surface if thallus area is less than 1 m2
      attenuation<-(-lambda*zed)-(kappa*cummulative.canopy2)
      Isurf.v<-Isurf*Ivary(TOD(timer),daylength/2)
      if (daylength==24) {Isurf.v<-Isurf} #allows user to set continuous light if daylength of 24 is chosen
      Ized<-Isurf.v*exp(attenuation)
      
      #need to fix Ized to surface value for all layers above water surface for multiplication with slice area matrix
      Icanop<-Ized[which(working.ht==0)]*(Ized==Isurf.v)
      Icanop2<-Ized*(Ized<=Ized[which(working.ht==0)])
      I.canop.final<-Icanop+Icanop2 
      
      Pmax.D<-Pmax.R*(d.slope*water.content+d.intercept) #photosynthesis in air, including desiccation. (is Modifier for net psn in Madsen and Maberly 1990, Pmax parameters not established yet)
      if (Pmax.D<0) {Pmax.D<-0} #Pmax in air can't go below zero, not clear if this happens, more info needed, but represents conditions where canopy unlikely to persist anyway.
      Pmax.modifier<-min.water.content*(1/recovery.threshold) #Pmax restored from more than 40% desiccation (default value for the recovery threshold, can vary by species)
      if (Pmax.modifier>1) {Pmax.modifier<-1}
      Pmax.R<-Pmax*Pmax.modifier #altered P.max for recovery 
      if (wd>0) {water.content<-water.content + 15} #90% recovery in 6 steps = 1.5 hours, this command needs to go after Pmax.D modifier if used, plus another equation for how P.max.d responds
      if(water.content>100) {water.content<-100}
      if(wd>0){Pmax.R<-Pmax.R*water.content/100} #simple linear recovery of Pmax with hydration
      if(wd>0) {Pmax.D<-Pmax.R}
      
      # calculation of photosynthesis for each canopy layer
      Psn.numerator<-Pmax.D*alpha*I.canop.final
      Psn.denominator<-sqrt(Pmax.D^2+(alpha*I.canop.final)^2)
      if (sum(Pmax.D*alpha*I.canop.final==0)) {areal.Psn<-c(rep(0,1000))}
      if (sum(Pmax.D*alpha*I.canop.final>0)) {areal.Psn<-Psn.numerator/Psn.denominator} # areal photosynthetic rate 
      gross.psn<-areal.Psn*area.slice.mat #frond by depth matrix, model expansions need to aggregate and disaggregate for mixes of different species if PI or other parameters different
      inst.resp<-R*Tsm*area.slice.mat #respiration rate
      inst.net.psn<-gross.psn-inst.resp #instantaneous net photosynthesis (gross - respiration)
      
      data.view[timer,2]<-wd #gather summary data
      data.view[timer,3]<-water.content #gather summary data
      data.view[timer,4]<-Isurf.v #gather summary data
      data.view[timer,5]<-sum(inst.net.psn) #gather summary data 
    }
    
    data.view<-as.data.frame(data.view)
    names(data.view)<-c("time.step","water.depth","water.content","surface.light","net.photosynthesis.rate")
    plot(data.view[,1]/96,data.view[,5],type = 'l',lwd=3,col="green",xlab = "Days",ylab = "net photosynthesis rate")
    #plot(ht,rowSums(inst.net.psn))
    
    output$TPlot <- renderPlot({
      plot(data.view[,1]/96,data.view[,2],type = 'l',lwd=3,col="navy",xlab = "Days",ylab = "water height above shore level (m)")
    })
    
    output$PhPlot <- renderPlot({
      plot(data.view[,1]/96,data.view[,4]/10,type = 'l',col="orange",lwd=3,xlab = "Days",ylab = "water content (black,%) and surface light/10",ylim = c(0,200))
      lines(data.view[,1]/96,data.view[,3],col="black",lwd=3)
    })
    
    # integrated net photosynthesis
    int.net.psn<-(15*60/2)*(sum(data.view[,5])+sum(data.view[,5])-data.view[1,5]-data.view[nrow(data.view),5]) #converts time step to seconds
    int.net.psn<-int.net.psn *12.011/(1000*1000) #g C m-2 
    
    re <- as.data.frame(cbind(TAI,int.net.psn,Pmax.R))
    output$TAIcaption<-renderTable({ re })
    
    
    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("data.view", ".csv", sep = "")
      },
      content = function(file) {
        write.csv(data.view, file, row.names = FALSE)
      }
    )
  })
  

}
# Create Shiny app ----
shinyApp(ui = ui, server = server)
