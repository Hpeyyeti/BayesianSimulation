#########################################################################################################
#
  ###########PROJECT: A CURVE-FREE BAYESIAN DESICION-THEORETIC DESIGN FOR TWO-AGENT PHASE I TRIALS########
#
#########################################################################################################
#LIST OF VARIABLES:
#n.trials=10000
#cohortsize is 1
#Target DLT theta0 is 0.2
#Strict ordering of dose combinations
#c.par=4
#a=c.par*true.DLT
#b=c.par*(1-true.DLT)
#true.DLT is same as true
#alpha0=1.2 eta0=1
#n.min=10 n.max=50 delata0=0.05 r1=0.5 r2=0.95
#Stopping rules:
#S1: no stopping before n.min is reached
#S2: stop when n.max is reached
#S3: stop when lowest dose is over-toxic
#S4: stop if every higher dose is over-toxic

true.A <- (1/100)*matrix(c(4,8,12,16,
                           10,14,18,22,
                           16,20,24,28,
                           22,26,30,34),nrow=4,ncol=4)

true.B <- (1/100)*matrix(c(2,4,6,8,
                           5,7,9,11,
                           8,10,12,14,
                           11,13,15,17),nrow=4,ncol=4)

true.C <- (1/100)*matrix(c(10,20,30,40,
                           25,35,45,55,
                           40,50,60,70,
                           55,65,75,85),nrow=4,ncol=4)

true.D <- (1/100)*matrix(c(44,48,52,56,
                           50,54,58,62,
                           56,60,64,68,
                           62,66,70,74),nrow=4,ncol=4)

true.E <- (1/100)*matrix(c(8,18,28,29,
                           9,19,29,30,
                           10,20,30,31,
                           11,21,31,41),nrow=4,ncol=4)

true.F <- (1/100)*matrix(c(12,13,14,15,
                           16,18,20,22,
                           44,45,46,47,
                           50,52,54,55),nrow=4,ncol=4)

Run.Clinical.Trials <- function(true.DLT, LFL) {
  ###Initiailse the parameters ###
  n.trials <- 10000
  theta0<-0.20 # target DLT
  cohort<-1
  c.par <- 4
  a <- c.par*true.DLT
  b <- c.par*(1-true.DLT)
  alpha0 <- 1.2
  eta0 <- 1
  nmin <- 10
  nmax <- 50
  delta0 <- 0.05
  r1 <- 0.5
  r2 <- 0.95
  
  MTD.count <- matrix(0,nrow=4,ncol=4) # records selection count for each dose comb
  S2_count <- 0
  S3_count <- 0
  S4_count <- 0
  n.pts <- matrix(0,nrow = 4,ncol = 4) #records count of patients treated at dose combination
  n.tox <- matrix(0,nrow = 4,ncol = 4) # records count of patients who experienced DLT
  k <- numeric(n.trials)
  
  for (trial in 1:n.trials){
    a <- c.par*true.DLT
    b <- c.par*(1-true.DLT)
    i <- 1 #dose A
    j <- 1 #dose B
    k[trial] <- 1 # count of patient
    
    ## Random dose escalation
    repeat{
      toxic.ind <- 0
      n.pts[i,j] = n.pts[i,j] + 1
      ## Check for toxicity     
      if (rbinom(1, 1, true.DLT[i,j]) == 1) {
        toxic.ind <- 1 # indicator for toxicity
        n.tox[i,j] =  n.tox[i,j] + 1
        break
      } else
        if (i == 4 & j == 4) { 
          #reached highest dose combination
          toxic.ind <- 0
          break 
        } else  
          if (i == 4) j <- j+1 else
            if (j == 4) i <- i+1 else
              if (sample(c('H','T'),1) == 'H') { #randomly escalating the dose
                i <- i+1
              } else {
                j <- j+1 
              }
      k[trial] <- k[trial] + 1 #go to next patient
    }
    
    ##Bayesian approach starts
    
    repeat{
      ## Update parameters a,b
      if (toxic.ind == 1){ # if current dose is toxic 
        a[i:4,j:4] = a[i:4,j:4] + 1
      } else if (toxic.ind == 0){ # if current dose is not toxic 
        b[1:i,1:j] = b[1:i,1:j] + 1
      }
      ##Stopping Rules:
      if (k[trial] >= nmin){ # S1
        if ((pbeta(theta0 + delta0, a[1,1], b[1,1], lower.tail=FALSE)) > r1) { 
          # S3:lowest dose toxic? No MTD suggested from S3
          S3_count <- S3_count + 1
          break
        } 
        if (!(i==4 & j==4)){ 
          
          # S4: if all higher doses  (r,s) > (i,j) are over-toxic?
          p <- pbeta(theta0 + delta0, a, b, lower.tail = FALSE)
          if (min(p[which((row(p) >= i) & (col(p) >= j) & ((row(p) + col(p)) > (i+j)))], 
                na.rm = TRUE) > r2) {
            MTD.count[i,j] <- MTD.count[i,j] + 1
            S4_count <- S4_count + 1
            break
          }
        }
      
        if (k[trial] == nmax) { # S2
          MTD.count[i,j] <- MTD.count[i,j] + 1
          S2_count <- S2_count + 1
          break
        } 
      }
      ## Trial continues..stopping criteria not met!!! 
      ## Suggest the dose combination for next patient
      E <- matrix(nrow=4,ncol=4)
      E <- -(alpha0 + eta0) * (theta0*pbeta(theta0, a, b) - 
                                 (a/(a+b))*pbeta(theta0, a+1, b)) - 
        eta0*((a/(a+b)) - theta0)
      
      # Check if LFL: No skipping? or skipping allowed?
      
      if (LFL == "NS") { # i.e. Dose skipping is not allowed
        # ## If (i,j) is current dose,next dose is (r,s) s.t. r ≤ i, j ≤ s, (i+1,j), (i,j+1) 
        y <- which(((row(E) <= i & col(E) <= j)|(row(E)==i+1 & col(E)==j)|(row(E)==i & col(E)==j+1)),arr.ind = T,useNames = F)
        util <- y[which(E[y]==max(E[y]),arr.ind = T)[1],]
        i <- util[1] # suggested dose A
        j <- util[2] # suggested dose B
      } else if (LFL == "WS") { # i.e. Dose skipping is allowed
        util <- which(E == max(E), arr.ind = TRUE, useNames = FALSE)
        ## proposed next dose combination (i,j) that maximaized E
        i <- util[1,1] # suggested dose A
        j <- util[1,2] # suggested dose B
      }
      
      k[trial] <- k[trial] + 1 #go to next patient; 
      n.pts[i,j] <- n.pts[i,j] + 1
      ## Check for toxicity     
      if (rbinom(1, 1, true.DLT[i,j]) == 1) { 
        toxic.ind <- 1
        n.tox[i,j] =  n.tox[i,j] + 1
      } else toxic.ind <- 0
    }
  }
  Display.Results(MTD.count, k, S3_count, true.DLT, LFL)
}

Display.Results <- function(local.MTD.count, local.k, local.S3_count, local.true.DLT, local.LFL) {
  cat ("Design LFL - ", local.LFL, "\n")
  cat ("MTD Selection % ", "\n")
  print ((local.MTD.count/10000)*100)
  cat ("\n\tPercentage of recommendation\n")
  
  if (sum(local.true.DLT) ==  sum(true.A)) {
    cat ("At target - ", (local.MTD.count[10]/10000)*100, "\n")
    
    target_1_10 <-  local.MTD.count[3] + local.MTD.count[4] + local.MTD.count[6] + local.MTD.count[7] +
                     local.MTD.count[8] + local.MTD.count[9] + local.MTD.count[11] + local.MTD.count[12] +
                     local.MTD.count[13] + local.MTD.count[14] + local.MTD.count[15] 
    cat ("1-10 pts of target - ", (target_1_10/10000)*100, "\n")
    
    target_10plus <- local.MTD.count[1] + local.MTD.count[2] + local.MTD.count[5] + local.MTD.count[16] 
    cat ("> or < 10 pts of target - ", (target_10plus/10000)*100, "\n")
    
  } else if (sum(local.true.DLT) == sum(true.B)) {
    cat ("At target - ", 0, "\n")
    
    target_1_10 <-  local.MTD.count[8] + local.MTD.count[11] + local.MTD.count[12] +
                    local.MTD.count[13] + local.MTD.count[14] + local.MTD.count[15] + local.MTD.count[16] 
    cat ("1-10 pts of target - ", (target_1_10/10000)*100, "\n")
    
    target_10plus <- local.MTD.count[1] + local.MTD.count[2] + local.MTD.count[3] + local.MTD.count[4] +
                     local.MTD.count[5] + local.MTD.count[6] + local.MTD.count[7] + local.MTD.count[9] + 
                     local.MTD.count[10]
    cat ("> or < 10 pts of target - ", (target_10plus/10000)*100, "\n")
    
  } else if (sum(local.true.DLT) == sum(true.C)) {
    cat ("At target - ", (local.MTD.count[2]/10000)*100, "\n")
    
    target_1_10 <-  local.MTD.count[3]  + local.MTD.count[5]
    cat ("1-10 pts of target - ", (target_1_10/10000)*100, "\n")
    
    target_10plus <- local.MTD.count[1] + local.MTD.count[4] + local.MTD.count[6] + local.MTD.count[7] +
                     local.MTD.count[8] + local.MTD.count[9] + local.MTD.count[10] + local.MTD.count[11] + 
                     local.MTD.count[12] + local.MTD.count[13] + local.MTD.count[14] + local.MTD.count[15] +  
                     local.MTD.count[16]
    cat ("> or < 10 pts of target - ", (target_10plus/10000)*100, "\n")
    
  } else if (sum(local.true.DLT) == sum(true.D)) {
    cat ("At target - ", 0, "\n")
    
    cat ("1-10 pts of target - ", 0, "\n")
    
    target_10plus <- sum(local.MTD.count)
    cat ("> or < 10 pts of target - ", (target_10plus/10000)*100, "\n")
    
  } else if (sum(local.true.DLT) == sum(true.E)) {
    cat ("At target - ", (local.MTD.count[10]/10000)*100, "\n")
    
    target_1_10 <-  local.MTD.count[2] + local.MTD.count[3] + local.MTD.count[4] + local.MTD.count[6] +
                    local.MTD.count[7] + local.MTD.count[8] + local.MTD.count[11] + local.MTD.count[13] +
                    local.MTD.count[14] 
    cat ("1-10 pts of target - ", (target_1_10/10000)*100, "\n")
    
    target_10plus <- local.MTD.count[1] + local.MTD.count[5] + local.MTD.count[9] +
                     local.MTD.count[12] + local.MTD.count[15] +  local.MTD.count[16] 
    cat ("> or < 10 pts of target - ", (target_10plus/10000)*100, "\n")
    
  } else if (sum(local.true.DLT) == sum(true.F)) {
    cat ("At target - ", (local.MTD.count[7]/10000)*100, "\n")
    
    target_1_10 <-  local.MTD.count[1] +local.MTD.count[2] + local.MTD.count[3] + local.MTD.count[4] + 
                    local.MTD.count[5] +local.MTD.count[6] + local.MTD.count[8] 
      
    cat ("1-10 pts of target - ", (target_1_10/10000)*100, "\n")
    
    target_10plus <- local.MTD.count[9] + local.MTD.count[10] + local.MTD.count[11] +  local.MTD.count[12] +
                     local.MTD.count[13] + local.MTD.count[14] + local.MTD.count[15] +  local.MTD.count[16] 
    cat ("> or < 10 pts of target - ", (target_10plus/10000)*100, "\n")
    
  }
  
  cat("None Recommended ", (local.S3_count/10000)*100, "\n") 
  cat("Average sample size ", mean(local.k), "\n") 

}

######FUNCTION CALLS####

Run.Clinical.Trials(true.A, "NS")
Run.Clinical.Trials(true.A, "WS")
# 
Run.Clinical.Trials(true.B, "NS")
Run.Clinical.Trials(true.B, "WS")
# 
Run.Clinical.Trials(true.C, "NS")
Run.Clinical.Trials(true.C, "WS")
#
Run.Clinical.Trials(true.D, "NS")
Run.Clinical.Trials(true.D, "WS")
#
Run.Clinical.Trials(true.E, "NS")
Run.Clinical.Trials(true.E, "WS")
# 
Run.Clinical.Trials(true.F, "NS")
Run.Clinical.Trials(true.F, "WS")

