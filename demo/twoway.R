  #####################
  ## Fourfold tables ##
  #####################
  
  ### Berkeley Admission Data ###
  ###############################
  data(UCBAdmissions)

  ## unstratified
  ### no margin is standardized
  x <- margin.table(UCBAdmissions, 2:1)
  fourfoldplot(x, std = "i", extended = FALSE)
  ### std. for gender
  fourfoldplot(x, margin = 1, extended = FALSE)
  ### std. for both
  fourfoldplot(x, extended = FALSE)
  
  ## stratified
  fourfoldplot(UCBAdmissions, extended = FALSE)
  fourfoldplot(UCBAdmissions) ## extended plots

  ### Coal Miners Lung Data ###
  #############################
  data(CoalMiners)
  
  ## Fourfold display, both margins equated
  fourfoldplot(CoalMiners, mfcol = c(2,4))

  ## Log Odds Ratio Plot
  summary(l <- oddsratio(CoalMiners))
  g <- seq(25, 60, by = 5)
  plot(l,
       xlab = "Age Group",
       main = "Breathelessness and Wheeze in Coal Miners")
  m <- lm(l ~ g + I(g^2))
  lines(fitted(m), col = "red")
  
  ## Fourfold display, strata equated
  fourfoldplot(CoalMiners, std = "ind.max", mfcol = c(2,4))
  
  ####################
  ## Sieve Diagrams ##
  ####################

  ### Hair Eye Color ###
  ######################
  data(HairEyeColor)

  ## aggregate over `sex':
  (tab <- margin.table(HairEyeColor, 1:2))
  ## plot expected values:
  sieveplot(t(tab), type = "expected", values = "both")

  ## plot sieve diagram:
  sieveplot(t(tab))

  ### Visual Acuity ###
  #####################
  data(VisualAcuity)
  sieveplot(Freq ~ right + left,
            data = VisualAcuity,
            subset = gender == "female",
            reverse.y = FALSE,
            main = "Unaided distant vision data",
            xlab = "Left Eye Grade",
            ylab = "Right Eye Grade")
  
  ### Berkeley Admission ###
  ##########################

  ## -> Larger tables: Cross factors
  ### Cross Gender and Admission
  data(UCBAdmissions)

  (tab <- xtabs(Freq ~ Dept + I(Gender : Admit), data = UCBAdmissions))
  sieveplot(tab, reverse.y = FALSE,
            xlab = "Gender:Admission",
            ylab = "Department",
            main = "Berkeley Admissions Data"
            )

  ######################
  ## Association Plot ##
  ######################
  
  ### Hair Eye Color ###
  ######################
  data(HairEyeColor)
  assocplot(margin.table(HairEyeColor, 1:2),
            col = c("blue","red"),
            xlab = "Hair Color",
            ylab = "Eye Color",
            main = "Association Plot")

  ####################
  ## Agreement Plot ##
  ####################

  ### Sexual Fun ###
  ##################
  data(SexualFun)

  ## Kappa statistics
  Kappa(SexualFun)

  ## Agreement Chart
  agreementplot(t(SexualFun), weights = 1)
  ## Partial Agreement Chart and B-Statistics
  (agreementplot(t(SexualFun),
                 xlab = "Husband's Rating",
                 ylab = "Wife's Rating",
                 main = "Husband's and Wife's Sexual Fun")
   )
  
  ### MS Diagnosis data ###
  #########################
  data(MSPatients)
  ## use e.g., X11(width = 12), or expand graphics device
  par(mfrow = c(1,2))
  agreementplot(t(MSPatients[,,1]), main = "Winnipeg Patients")
  agreementplot(t(MSPatients[,,2]), main = "New Orleans Patients")
  par(mfrow = c(1,1))

  ##################
  ## Ternary Plot ##
  ##################

  ### sample data ###
  ###################
  (x <- rbind(c(A=10,B=10,C=80),
              c(40,30,30),
              c(20,60,20)
              )
   )
  ternaryplot(x,
              cex = 2,
              col = c("black", "blue", "red"),
              coordinates = TRUE
              )

  ### Arthritis Treatment Data ###
  ################################
  data(Arthritis)
  
  ## Build table by crossing Treatment and Sex
  (tab <- as.table(xtabs(~ I(Sex:Treatment) + Improved, data = Arthritis)))
  
  ## Mark groups
  col <- c("red", "red", "blue", "blue")
  pch <- c(1, 19, 1, 19)
  
  ## plot
  ternaryplot(
              tab,
              col = col,
              pch = pch,
              cex = 2,
              bg = "lightgray",
              grid.color = "white",
              labels.color = "white",
              main = "Arthritits Treatment Data"
              )
  ## legend
  legend(0.7, 0.8,
         c("GROUP", rownames(tab)),
         pch = c(NA, pch),
         col = c(NA, col)
         )

  ### Baseball Hitters Data ###
  #############################
  data(Hitters)
  attach(Hitters)

  colors <- c("black","red","green","blue","red","black","blue")
  pch <- substr(levels(Positions), 1, 1)
  ternaryplot(
              Hitters[,2:4],
              pch = as.character(Positions),
              col = colors[codes(Positions)],
              main = "Baseball Hitters Data"
              )
  legend(
         0.8, 0.9,
         legend = c("POSITION(S)", levels(Positions)),
         pch = c("", pch),
         col = c(NA, colors)
         )

  detach(Hitters)

  ### Lifeboats on the Titanic ###
  ################################
  data(Lifeboats)
  attach(Lifeboats)

  ternaryplot(
              Lifeboats[,4:6],
              pch = ifelse(side=="Port", 1, 19),
              col = ifelse(side=="Port", "red", "blue"),
              id  = ifelse(men/total > 0.1, as.character(boat), NA),
              dimnames.position = "edge",
              dimnames = c("Men of Crew", "Men passengers", "Women and Children"),
              main = "Lifeboats on the Titanic"
              )
  legend(
         0.7, 0.8,
         legend = c("SIDE", "Port", "Starboard"),
         pch = c(NA, 1, 19),
         col = c("black", "red", "blue"),
         )

  ## Load against time for Port/Starboard boats
  plot(launch, total,
       pch = ifelse(side == "Port", 1, 19),
       col = ifelse(side == "Port", "red", "darkblue"),
       xlab = "Launch Time",
       ylab = "Total loaded",
       main = "Lifeboats on the Titanic"
       )
  legend(as.POSIXct("1912-04-15 01:48:00"), 70,
         legend = c("SIDE","Port","Starboard"),
         pch = c(NA, 1, 19),
         col = c(NA, "red", "darkblue")
         )
  text(as.POSIXct(launch),
       total,
       labels = as.character(boat),
       pos = 3,
       offset = 0.3
       )
  abline(lm(total ~ as.POSIXct(launch),
            subset = side == "Port"),
         col = "red")     
  abline(lm(total ~ as.POSIXct(launch),
            subset = side == "Starboard"),
         col = "darkblue")     

  detach(Lifeboats)


