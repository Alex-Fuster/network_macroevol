########################################
sim_comp_neutral = function(seed,  nsteps) {
  

    
    set.seed(seed)
    
    # Draw the first species trait
    traits_mat = matrix(nr = Smax, nc = 1)
    traits_mat[1,] = rnorm(1, mean = 0, sd = 1)
    traits_mat = as.data.frame(traits_mat)
    names(traits_mat) = c("z")
    
    # Set the presence/absence matrix
    pres = matrix(0, nr = nsteps, nc = Smax)
    pres[1,1] = 1
    
    # Set the ancestry object
    anc = matrix(NA,nr = Smax, nc = 3)
    
    # Set the extinctions object
    extinct = matrix(NA,nr = Smax, nc = 2)
    
    #  [TABLE DIST_ANC] Set the distance ancestry matrix -----------------------------------
    # (new_spp, ancestor, tip, distance)
    dist_anc = as.data.frame(matrix(NA, nr = Smax, nc = 4))
    colnames(dist_anc) = c("spp", "ancestor", "A/E", "distance")
    dist_anc[1,] <- c(1, 0, "A", 0)
    dist_anc$distance <- as.numeric(dist_anc$distance)
    dist_anc$distance[is.na(dist_anc$distance)] <- 0
    
    #  [TABLE DIST_ANC] record the dist_anc at each timestep
    
    list_dist_anc <- list()
    
    # Record the matrices
    L_list = list()
    
    
    # Species count
    S = 1
    
    Stotact <- matrix(NA, nsteps, 3)
    Stotact <- as.data.frame(Stotact)
    names(Stotact) = c("Step","TotalS","ActualS")
    
    ##########################
    # MAIN LOOP
    for(step in 2:nsteps) {
      ActualS = sum(pres[step-1,])
      #cat("Step = ", step-1," Total S = ", S," Actual S = ", ActualS,'\n')
      Stotact[step-1,] <- c(step-1, S, ActualS)
      
      
      
      
      
      # Test for successful origination probability
      for(esp in 1:Smax) {
        
        if(S >= Smax) break
        
        # Speciation occurs if the species is present
        if(pres[step,esp] == 1) {
          
          
          # [TABLE DIST_ANC] Define distances of Alive spp -------------------
          
          dist_anc[esp,"A/E"] <- "A"
          dist_anc[esp,"distance"] <- as.numeric(dist_anc[esp,"distance"])+1
          
          
          # Species is maintained
          pres[step, esp] = 1
          
          # Test if there is speciation
          test_number = runif(1, 0, 1)
          speciation_prob = u_max/(1 + exp(d * (ActualS - I_max)))
          if(test_number < speciation_prob) {
            
            # Ancestor trait
            traits_anc <- traits_mat[esp,]
            
            # Pick new trait for the first mutant
            traits_mut = rnorm(1, mean = traits_anc, sd = 0.2)
            
           
            # get vector of traits of present species
            
            vec_traits <- traits_mat[which(pres[step,] == 1),]
            
            # get vector of interaction strengths
            
            vec_interaction_strengths <- 1/(vec_traits-traits_mut)
            
            
        
            # Compute the probability of successful establishment #---------------------- P(establishment)
            
            estab_prob = u_0neg + u_1neg*exp(-a_uneg * sum(vec_interaction_strengths))
            
            
  
            # Test if there is establishment of the new species
            
            ## We start to test for one of them
            if(runif(1) < estab_prob) {
              S <- S + 1
              traits_mat[S,] = traits_mut
              pres[step,S] = 1
              anc[S,] = c(step, esp, S) #Time step, ancestor, new_species
              
              # [TABLE DIST_ANC] add new spp to table dist_anc ------------
              
              dist_anc[S, ] <- c(S, esp, "A", 1)
              
              
              if(S >= Smax) break
              
              #pres[step,esp] = 0 #the ancestor disapear because of trait deplacement
            }
            
            # Reasign original traits to the ancestor (not needed I think)
           # traits_mat[esp,] <- traits_anc
            
            
          }
        }
      }
      
      
      if(S >= Smax) break
      #if(ActualS == 0) break
      
      # Compute the weighted interaction matrix among present species
      
   #   traits <- traits_mat[, 1]
      traits <- traits_mat[which(pres[step,] == 1),1]
      
      L <- inv( outer( traits, traits, FUN= "-"))
      
      # crop matrix to present species
      
      #pres_vec = pres[step,]
      
      #L = L[pres_vec,pres_vec]
      
      L_list[[step]] = L                            #---------------------- List of networks/time-step
      
     
      
      # Test for extinction                         #---------------------- P(extinction)
      
      
      # L is a weighted matrix
      
      graph_L <- graph_from_adjacency_matrix(L, weighted=TRUE)
      
      vec_strenghts <- strength(graph_L, mode = "in")
      
      # compute probabilities based on in-degree strenghts
      
      ext_prob = e_0neg * (1 - exp(-a_eneg*0.1*vec_strenghts))
      
      vec_ext_prob <- numeric(Smax)
      
      vec_ext_prob[which(pres[step,] == 1)] <- ext_prob
      
    
      # Perform extinctions
     
      present_spe <- grep(1,pres[step,])
      test_extprob <- rep(0,Smax)
      random_number <- runif(length(present_spe),0,1)
      test_extprob[present_spe] <- random_number
      #cat("S actuel = ", ActualS," longueur random_number = ", length(random_number),'\n')
      pres[step, pres[step-1,] & test_extprob < vec_ext_prob] = 0
      
      if(step > 1){
        
        for(i in 1:Smax) {
          if(pres[step,i] != pres[step-1,i] & pres[step,i] == 0){
            extinct[i,] = c(step, i)
            
            # [TABLE DIST_ANC] Define a new E spp -------------------
            
            dist_anc[i,"A/E"] <- "E"
            
          }
        }
        
        
      }
      
      
      
      # [TABLE DIST_ANC] save table dist_anc at timestep step
      
      list_dist_anc[[step]] <- na.omit(dist_anc)
      
    } # End of main loop
    
    list(pres = pres, 
         traits = traits_mat, 
         anc = anc, 
         extinct = extinct, 
         dynamic = Stotact, 
         L = L, 
         L_list = L_list,
         dist_anc = dist_anc,
         list_dist_anc = list_dist_anc)

}
