march.loop.cov <- function(p){
	population <- list()
	fitness <- array(0,p@populationSize)
	for( i in 1:p@populationSize ){
		population[[i]] <- p@initParameters@fct(p@initParameters)
		fitness[i] <- p@evalParameters@fct(population[[i]],p@evalParameters)
	}
	
	for( g in 1:p@generation ){		
		# Roulette wheel proportional selection
		fitnessProb <- fitness/sum(fitness)
		prop <- cumsum(fitnessProb)-fitnessProb[1]
		
		childrenPopulation <- list()
		childrenFitness <- array(0,p@populationSize)
    
		
		for( i in seq(1,p@populationSize,2) ){
			sel1 <- max(which(prop<runif(1)))
			sel2 <- max(which(prop<runif(1)))
			
			# Do we need to crossover
			if( runif(1)<p@crossoverProb ){	
				children <- p@crossoverParameters@fct(population[[sel1]],population[[sel2]],p@initParameters@AConst)
				child1 <- children$c1
				child2 <- children$c2
			}
			else{
				child1 <- population[[sel1]]
				child2 <- population[[sel2]]
			}
			
			child1 <- p@mutationParameters@fct(child1,p@mutationParameters)
			child2 <- p@mutationParameters@fct(child2,p@mutationParameters)
			
			if( p@optimizing ){
				child1 <- p@optimizingParameters@fct(child1,p@optimizingParameters)
				child2 <- p@optimizingParameters@fct(child2,p@optimizingParameters)
			}
      
  		
			# eval the generated individual
			childrenFitness[i] <- p@evalParameters@fct(child1,p@evalParameters)
			childrenFitness[i+1] <- p@evalParameters@fct(child2,p@evalParameters)	

      childrenPopulation[[i]] <- child1
			childrenPopulation[[i+1]] <- child2
		}
		
		# replacement step
		newPopulation <- list()
		newFitness <- array(0,p@populationSize)
		for( i in 1:p@populationSize){
			# best of parents 
			maxParentId <- which.max(fitness)
			maxChildrenId <- which.max(childrenFitness)
			if( fitness[maxParentId]>childrenFitness[maxChildrenId] ){
				newPopulation[[i]] <- population[[maxParentId]]
				newFitness[i] <- fitness[maxParentId]
				fitness[maxParentId] <- 0
			}
			else{
				newPopulation[[i]] <- childrenPopulation[[maxChildrenId]]
				newFitness[i] <- childrenFitness[maxChildrenId]
				childrenFitness[maxChildrenId] <- 0
			}
		}

		population <- newPopulation
		fitness <- newFitness
		
		meanFitness <- mean(fitness)
		bestFitness <- max(fitness)
		cat(sprintf("Generation %d Mean fitness: %f Best Fitness: %f\n",g,meanFitness,bestFitness))
	}
	
	list(children=childrenPopulation,best=population[[which.max(fitness)]])
}


march.dcmm.cov.ea.initialization <- function(p){
  	NbAMCovar <- sum(p@AMCovar)
  	NbCMCovar <- sum(p@CMCovar)
  	KCovar <- p@y@Kcov
  	AMCovar <- p@AMCovar
  	CMCovar <- p@CMCovar
  	M <- p@M
  	orderVC <- p@orderVC
  	orderHC <- p@orderHC
  	Amodel <- p@Amodel
  	Cmodel <- p@Cmodel
  	K <- p@y@K
  	y <- p@y
  	AConst <- p@AConst
  
  	d <- march.dcmm.cov.constructEmptyDcmm(p@M,p@y,p@orderVC,p@orderHC,p@AMCovar,p@CMCovar,p@Amodel,p@Cmodel)
   
	  placeACovar <- which(AMCovar==1)
  
  	AtmCovar <- 1
  	if(NbAMCovar>0){
  		for (i in 1:NbAMCovar){
      		AtmCovar <- AtmCovar*KCovar[placeACovar[i]]
      	}
  	}
  
  	placeCCovar <- which(CMCovar==1)
  	CtmCovar <- 1
  	if(NbCMCovar>0){
  		for (i in 1:NbCMCovar){
      		CtmCovar <- CtmCovar*KCovar[placeCCovar[i]]
    	}
  	}
    
  	#Initialisation of all the things involved in the hidden process
  	if(AConst==TRUE){
  	  AQ <- array(diag(M),c(1,M,M))
  	  A <- diag(M)
  	  RA <- diag(M)
  	  ATCovar <- list()
  	  APhi <- array(1,c(1,1))
  	  mccov <- march.mccov.sk(M,1,AtmCovar,KCovar,placeACovar,NbAMCovar,matrix(AQ[1,,],M,M),ATCovar)
  	  AProbT <- mccov$ProbT
  	  
  	}else{
  	  if(Amodel=="complete"){
  	
		  #Random initialisation of the full transition matrix between hidden states
      	if(M>1){
    	  	AQ <- 0.1+array(runif(M^max(orderHC,1)*M),c(1,M^max(orderHC,1),M))
    	  	AQ[1,,] <- AQ[1,,]/rowSums(AQ[1,,])
    		
    	  	if(orderHC==0){
    		  	for(m in 2:M){
    			  	AQ[1,m,] <- AQ[1,1,]
    			  }
    		  }
      	}else{
    	  	AQ <- array(1,c(1,1,1))
    	  }
    	
    	  #Random initialisation of the vector of lag Phi, size 1+NbAMCovar
    	  APhi <- 0.1+array(runif(1+NbAMCovar),c(1,1+NbAMCovar))
    	  APhi <- APhi/sum(APhi)
    
    	  #Random initialisation of the matrices of transition between covariates and hidden states
    	  ATCovar <- list()
    	  if(NbAMCovar>0){
    		  for(i in 1:NbAMCovar){
      			ATCovar[[i]] <- 0.1+array(runif(KCovar[placeACovar[i]]*M),c(KCovar[placeACovar[i]],M))
      			ATCovar[[i]] <- ATCovar[[i]]/rowSums(ATCovar[[i]])
    	  	}
    	  }
    
    
        mccov <- march.mccov.sk(M,max(orderHC,1),AtmCovar,KCovar,placeACovar,NbAMCovar,matrix(AQ[1,,],M^max(orderHC,1),M),ATCovar)
    	  AProbT <- mccov$ProbT
    
    	  l <- GMTD_tm_cov(max(orderHC,1),M,APhi,AProbT)
    	  A <- l$HOQ
    	  RA <- l$CRHOQ
    	
  	  }else{
    	  A <- array(0,c(M^orderHC*AtmCovar,M^orderHC))
    
    	  #Random initialisation of AQ
    	  if(Amodel=="mtd"){
      		AQ <- 0.1 + array(runif(M*M),c(1,M,M))
      		AQ[1,,] <- AQ[1,,]/rowSums(AQ[1,,])
    	  }else{
      		AQ <- 0.1+array(runif(M*M*orderHC),c(orderHC,M,M))
      		for(i in 1:orderHC){
        		AQ[i,,] <- AQ[i,,]/rowSums(AQ[i,,])
      		}
    	  }
    
    	  #Random initialisation of APhi
    	  APhi <- 0.1+array(runif(orderHC+NbAMCovar),c(1,orderHC+NbAMCovar))
    	  APhi <- APhi/sum(APhi)
    
    	  ANSS <- matrix(0,M^(orderHC+1)*AtmCovar,1)
    	  AValT <- BuildArrayCombinations(M,orderHC,KCovar[placeACovar],NbAMCovar)
    
    	  #Random initialisation of ATCovar
    	  ATCovar <- list()
    	  if(NbAMCovar>0){
      		for(i in 1:NbAMCovar){
        		ATCovar[[i]] <- 0.1+array(runif(KCovar[placeACovar[i]]*M),c(KCovar[placeACovar[i]],M))
        		ATCovar[[i]] <- ATCovar[[i]]/rowSums(ATCovar[[i]])
      		}
    	  }
    	
    	  #The function BuildArray Q is in march.mtd.R
    	  AProbT <- BuildArrayQ(M,l=orderHC,i0_il=AValT,n_i0_il = ANSS,q=AQ,kcov=KCovar[placeACovar],ncov=NbAMCovar,S=ATCovar)
    	  l <- GMTD_tm_cov(orderHC,M,APhi,matrix(AProbT,M^(orderHC+1)*AtmCovar,orderHC+NbAMCovar))

    	  A <- l$HOQ
    	  RA <- l$CRHOQ  
      }
  	}
  	#Initialisation of all the things involved in the visible process
  	if(Cmodel=="complete"){
  	
  		#Random initialisation of CQ
    	CQ <- 0.1+array(runif(K^orderVC*K*M),c(1,K^orderVC,K,M))
   		if(orderVC>0){
    		for(state in 1:M){
      			CQ[1,,,state] <- CQ[1,,,state]/rowSums(CQ[1,,,state])
    		}
    	}else{
      		for(state in 1:M){
        		CQ[1,,,state] <- CQ[1,,,state]/sum(CQ[1,,,state])
      		}
    	}
    
    	#Random initialisation of CPhi
      	CPhi <- array(1,c(1,1+NbCMCovar,M))
    	
		if(NbCMCovar>0){
			for(state in 1:M){
				CPhi[1,,state] <- 0.1+array(runif(dim(CPhi)[2]),c(1,dim(CPhi)[2]))
				CPhi[1,,state] <- CPhi[1,,state]/sum(CPhi[1,,state])
          	}
     	}
     
    	#Random initialisaton of CTCovar
		CTCovar <- list()
		CTCovartmp <- list()
    	if(NbCMCovar>0){
      		c <- BuildContingencyTableDCMM(y,orderVC,NbCMCovar,placeCCovar)
      		for(i in 1:NbCMCovar){
        		CTCovar[[i]] <- array(0,c(KCovar[placeCCovar[i]],K,M))
        		tmp <- NormalizeTable(c[[orderVC+i]])
        		CTCovartmp[[i]] <- tmp
        		for(j in 1:M){
          			CTCovar[[i]][,,j] <- tmp
        		}
      		}
    	}
    
    	
    	CNSS <- BuildArrayNumberOfSequencesDCMM(y,orderVC,CMCovar,NbCMCovar,placeCCovar)
    	CProbT <- array(0,c(K^(orderVC+1)*CtmCovar,1+NbCMCovar,M))
    	CValT <- BuildArrayCombinations(K,orderVC,KCovar[placeCCovar],NbCMCovar)
    
    	RC <- array(0,c(K^orderVC*CtmCovar,K,M))
    	for(state in 1:M){
        	mccov <- march.mccov.sk(K,orderVC,CtmCovar,KCovar,placeCCovar,NbCMCovar,matrix(CQ[1,,,state],K^orderVC,K),CTCovartmp)
        	CProbT[,,state] <- mccov$ProbT
        	RC[,,state] <- GMTD_tm_cov(orderVC,K,t(CPhi[1,,state]),matrix(CProbT[,,state],K^(orderVC+1)*CtmCovar,1+NbCMCovar))$CRHOQ
        }
	}else{
    	c <- BuildContingencyTableDCMM(y,orderVC,NbCMCovar,placeCCovar)
    	u <- CalculateTheilUDCMM(y,orderVC,c,NbCMCovar,placeCCovar)
   		#Initval corresponds to init=best
   		Initval <- 1
   		CPhi <- array(0,c(1,orderVC+NbCMCovar,M))
   		for(state in 1:M){      			
   			CPhi[,,state] <- u/sum(u)
    	}
    	
    	if(Cmodel=="mtd"){
      		CQ <- array(0,c(1,K,K,M))
      		#Maybe optimize this part since CQ is equal for every state
      		for(state in 1:M){
      			CQ[1,,,state] <- march.dcmm.initq.c(y,1,Initval,orderVC,c,u)
      		}
    	}else{
      		CQ <- array(0,c(orderVC,K,K,M))
      		for(state in 1:M){
      			CQ[,,,state] <- march.dcmm.initq.c(y,2,Initval,orderVC,c,u)
      		}
    	}
    
    	CTCovar <- list()
    	CTCovartmp <- list()
    	#Initialisation of the matrices of transition between covariates and the dependant variable by normalizing the crosstables.
    	if(NbCMCovar>0){
      		c <- BuildContingencyTableDCMM(y,orderVC,NbCMCovar,placeCCovar)
      		for(i in 1:NbCMCovar){
        		CTCovar[[i]] <- array(0,c(KCovar[placeCCovar[i]],K,M))
        		tmp <- NormalizeTable(c[[orderVC+i]])
        		CTCovartmp[[i]] <- tmp
        		for(j in 1:M){
          			CTCovar[[i]][,,j] <- tmp
        		}
      		}
    	}
    
    	CProbT <- array(0,c(K^(orderVC+1)*CtmCovar,orderVC+NbCMCovar,M))
    	CValT <- BuildArrayCombinations(K,orderVC,KCovar[placeCCovar],NbCMCovar)
    	CNSS <- BuildArrayNumberOfSequencesDCMM(y,orderVC,CMCovar,NbCMCovar,placeCCovar)

    	RC <- array(0,c(K^orderVC*CtmCovar,K,M))
    	if(Cmodel=="mtd"){
    		t <- array(CQ[,,,1],c(1,K,K))
    	}else{
      		t <- array(CQ[,,,1],c(orderVC,K,K))
    	}
    	q_i0_il <- array(0,c(K^(orderVC+1)*CtmCovar,orderVC+NbCMCovar))
    	q_i0_il <- BuildArrayQ(K,l=orderVC,i0_il=CValT,n_i0_il = CNSS,q=t,kcov=KCovar[placeCCovar],NbCMCovar,S=CTCovartmp)
    	for(state in 1:M){
      		CProbT[,,state] <- q_i0_il
	      	tm <- GMTD_tm_cov(orderVC,K,as.matrix(t(CPhi[,,state])),array(CProbT[,,state],c(K^(orderVC+1)*CtmCovar,orderVC+NbCMCovar)))
			RC[,,state] <- tm$CRHOQ
   		}
	}
    
    #Initialisation of Pi
    if(M>1){
    	if(orderHC==0){
    		Pi <- array(0,c(AtmCovar,M,1))
    		Pi[,,1] <- RA[1:AtmCovar,]
    	}else{
    		Pi <- array(0,c(M^(orderHC-1)*AtmCovar,M,orderHC))
    		Pi[c(1:(M^(1-1)*AtmCovar)),,1] <- 0.1+array(runif(M^(1-1)*AtmCovar*M),c(M^(1-1)*AtmCovar,M))
    		if(AtmCovar>1){
   				Pi[c(1:(M^(1-1)*AtmCovar)),,1] <- Pi[c(1:(M^(1-1)*AtmCovar)),,1]/rowSums(Pi[c(1:(M^(1-1)*AtmCovar)),,1])
   			}else{
   				Pi[c(1:(M^(1-1)*AtmCovar)),,1] <- Pi[c(1:(M^(1-1)*AtmCovar)),,1]/sum(Pi[c(1:(M^(1-1)*AtmCovar)),,1])
   			}
    		if(orderHC>1){
    			for(i in 2:orderHC){
   					Pi[c(1:(M^(i-1)*AtmCovar)),,i] <- 0.1+array(runif(M^(i-1)*AtmCovar*M),c(M^(i-1)*AtmCovar,M))
      				Pi[c(1:(M^(i-1)*AtmCovar)),,i] <- Pi[c(1:(M^(i-1)*AtmCovar)),,i]/rowSums(Pi[c(1:(M^(i-1)*AtmCovar)),,i])
       			}
       		}
       	}
    }else{
    	Pi <- array(1,c(AtmCovar,1,1))
  	}
  	
    d@A <- A
    d@RB <- RC
    d@Pi <- Pi
    d@ATCovar <- ATCovar
    d@CTCovar <- CTCovar
    d@APhi <- APhi
    d@CPhi <- CPhi
    d@AQ <- AQ
    d@CQ <- CQ
    d@AMCovar <- AMCovar
    d@CMCovar <- CMCovar
    d@Amodel <- Amodel
    d@Cmodel <- Cmodel
    d@AProbT <- AProbT
    d@CProbT <- CProbT
	d    
}

#Initialization of the transition matrix Q for the MTD modeling of the visible process of a DCMM
#Model: 1=MTD, 2=MTDg
#Initval: 1)Q best, 2)Q weighted 3)Q random
march.dcmm.initq.c<-function(y,Model,Initval,order,c,u){
  K=y@K
  N=y@N
  T=y@T
  if(Model==2){
    Q <- array(NA,c(order,K,K))
    for (g in 1:order){
      Q[g,,] <- NormalizeTable(c[[g]])
    }
  }
  else{
  if(Initval==3){
    Q=0.1+array(runif(K*K),c(K,K))
    Q=Q/rowSums(Q)
  }
  if(Initval==1){
    u_tilde<-u[1:order]
    k <- which.max(u_tilde)
    Q <- array(0,c(1,K,K))
    Q[1,,] <- NormalizeTable(c[[k]])
  }
  if(Initval==2){
    Q <- array(0,c(1,K,K))
    q_tilde <- array(0,dim=c(1,K,K))
    for (g in 1:order){
      q_tilde[1,,] <- q_tilde[1,,] + u[g] * c[[g]]
    }
    Q[1,,] <- NormalizeTable(q_tilde[1,,])
  }
  }
  return(Q)
}

march.mccov.sk<-function(K,order,tmCovar,KCovar,placeCovar,NbMCovar,Q,TCovar){
	NSS <- rep(0,c(K^(order+1)*tmCovar))
	ProbT <- array(0,c(K^(order+1)*tmCovar,1+NbMCovar))
	ValT <- array(0,c(K^(order+1)*tmCovar,1+order+NbMCovar))
	
	values<-1:K
	for(i in 1:(order+1)){
		ValT[,i] <- t(kronecker(values,rep(1,K^(order+1)*tmCovar/K^i)))
		values <- kronecker(rep(1,K),values)
	}
	
	if(NbMCovar>0){
		totm <- tmCovar
		totp <- K^(order+1)
		
		for(i in 1:NbMCovar){
			totm <- totm/KCovar[placeCovar[i]]
			values <- 1:KCovar[placeCovar[i]]
			ValT[,order+1+i] <- kronecker(rep(1,totp),kronecker(values,rep(1,totm)))
			totp <- totp*KCovar[placeCovar[i]]
		}
	}
	
	for(i in 1:length(NSS)){
		if(order==0){
			row <- 1
		}else{
			row <- ValT[i,order+1]
			for(g in march.h.seq(order,2,-1)){
				row <- row+K^(order-g+1)*(ValT[i,g]-1)
			}
		}
		ProbT[i,1] <- Q[row,ValT[i,1]]
		if(NbMCovar>0){
			for(j in 1:NbMCovar){
				ProbT[i,j+1] <- TCovar[[j]][ValT[i,order+1+j],ValT[i,1]]
			}
		}
	}
	list(ValT=ValT,ProbT=ProbT)
}

march.mccov.sp<-function(order,Q,ValT,ProbT){
	dq <- dim(Q)[2]
	
	for(i in 1:dim(ValT)[1]){
		if(order==0){
			row <- 1
		}else{
			row <- ValT[i,order+1]
			for(g in march.h.seq(order,2,-1)){
				row <- row+dq^(order-g+1)*(ValT[i,g]-1)
			}
		}
		
		ProbT[i,1] <- Q[row,ValT[i,1]]
	}
	ProbT
}


march.dcmm.cov.ea.evaluation <- function(d,p){
  d@ll <- march.dcmm.h.cov.computeLL(d=d,y=p@ds)
  d@dsL <- sum(p@ds@T)
  
  -1/d@ll
}

march.dcmm.h.cov.computeLL <- function(d,y){
   	placeACovar <- which(d@AMCovar==1)
   	AtmCovar <- 1
   	if(sum(d@AMCovar)>0){
		  for (i in 1:sum(d@AMCovar)){
      		AtmCovar <- AtmCovar*y@Kcov[placeACovar[i]]
		  }
  	}
  
  	placeCCovar <- which(d@CMCovar==1)
  	CtmCovar <- 1
  	if(sum(d@CMCovar)>0){
     	for (i in 1:sum(d@CMCovar)){
      		CtmCovar <- CtmCovar*y@Kcov[placeCCovar[i]]
    	}
  	}

  
  	ll <- 0
  
  	for( n in 1:y@N){ 
    	if( y@T[n]>(d@orderVC+d@orderHC) ){
     		s <- march.dataset.h.extractSequence(y,n)
     		ll <- ll+march.dcmm.fp.cov(d,s,y@cov[n,,],AtmCovar,CtmCovar)$LLAlpha
    	}
  	}
  	ll
}

march.dcmm.cov.ea.mutation <- function(d,p){
  
  NbAMCovar <- sum(d@AMCovar)
  NbCMCovar <- sum(d@CMCovar) 
  	
  placeACovar <- which(d@AMCovar==1)
  placeCCovar <- which(d@CMCovar==1)
  	
  AtmCovar <- 1
 	if(NbAMCovar>0){
    for (i in 1:NbAMCovar){
      AtmCovar <- AtmCovar*d@y@Kcov[placeACovar[i]]
    }
  }
  	
  CtmCovar <-1
  if(NbCMCovar>0){
  	for(i in 1:NbCMCovar){
  		CtmCovar <- CtmCovar*d@y@Kcov[placeCCovar[i]]	
  	}
  }

	#Mutation of APhi
	Adim<-dim(d@APhi)[2]
	d@APhi[1,]<-d@APhi[1,]+(runif(Adim)<p@pMut)*rnorm(Adim)
	d@APhi[1,]<-march.ea.checkBounds(v=d@APhi[1,],lb=0,ub=1)
	d@APhi[1,]<-march.dcmm.ea.h.scale(d@APhi[1,],Adim)
	
	#Mutation of CPhi
	if(NbCMCovar>0 | d@orderVC>0){
		Cdim<-dim(d@CPhi)[2]
		for( i in 1:d@M){
			d@CPhi[1,,i]<-d@CPhi[1,,i]+(runif(Cdim)<p@pMut)*rnorm(Cdim)
			d@CPhi[1,,i]<-march.ea.checkBounds(v=d@CPhi[1,,i],lb=0,ub=1)
			d@CPhi[1,,i]<-march.dcmm.ea.h.scale(d@CPhi[1,,i],Cdim)
		}
	}	
	
  if(sum(d@AMCovar>0)){
  	for (i in 1:sum(d@AMCovar)){
  		for(j in 1:d@y@Kcov[placeACovar[i]]){
  			d@ATCovar[[i]][j,] <- d@ATCovar[[i]][j,]+(runif(d@M)<p@pMut)*rnorm(d@M)
  			d@ATCovar[[i]][j,] <- march.ea.checkBounds(v=d@ATCovar[[i]][j,],lb=0, ub=1)
  			d@ATCovar[[i]][j,] <- march.dcmm.ea.h.scale(d@ATCovar[[i]][j,],d@M)
  		}
  	}
  }
  
  if(sum(d@CMCovar>0)){
  	for (i in 1:sum(d@CMCovar)){
  		for (j in 1:d@y@Kcov[placeCCovar[i]]){
			  for (k in 1:d@M){
				  d@CTCovar[[i]][j,,k]<-d@CTCovar[[i]][j,,k]+(runif(d@y@K)<p@pMut)*rnorm(d@y@K)
				  d@CTCovar[[i]][j,,k]<-march.ea.checkBounds(v=d@CTCovar[[i]][j,,k],lb=0,ub=1)
				  d@CTCovar[[i]][j,,k]<-march.dcmm.ea.h.scale(d@CTCovar[[i]][j,,k],d@y@K)
  		  }
  	  }
  	}
  }
    
  if(p@AConst==FALSE){
    if(d@orderHC>0){
		  for(i in 1:dim(d@AQ)[1]){
  		  for(j in 1:dim(d@AQ)[2]){
  			  d@AQ[i,j,]<-d@AQ[i,j,]+(runif(d@M)<p@pMut)*rnorm(d@M)
  				d@AQ[i,j,]<-march.ea.checkBounds(v=d@AQ[i,j,],lb=0,ub=1)
  				d@AQ[i,j,]<-march.dcmm.ea.h.scale(d@AQ[i,j,],d@M)
        }
      }
    }else{
  	  d@AQ[1,1,] <- d@AQ[1,1,]+(runif(d@M)<p@pMut)*rnorm(d@M)
  		d@AQ[1,1,]<-march.ea.checkBounds(v=d@AQ[1,1,],lb=0,ub=1)
  		d@AQ[1,1,]<-march.dcmm.ea.h.scale(d@AQ[1,1,],d@M)
  		for(i in 1:d@M){
  			d@AQ[1,i,] <-d@AQ[1,1,]
  		}
    }
  }
  
  for( i in 1:dim(d@CQ)[1]){
  	for(j in 1:dim(d@CQ)[2]){
  		for(k in 1:d@M){
  			d@CQ[i,j,,k]<-d@CQ[i,j,,k]+(runif(d@y@K)<p@pMut)*rnorm(d@y@K)
  			d@CQ[i,j,,k]<-march.ea.checkBounds(v=d@CQ[i,j,,k],lb=0,ub=1)
  			d@CQ[i,j,,k]<-march.dcmm.ea.h.scale(d@CQ[i,j,,k],d@y@K)
  		}
  	}
  }
    
  for( i in 1:(d@orderHC) ){
    for( j in 1:(d@M^(i-1)*AtmCovar)){
      d@Pi[j,,i] <- d@Pi[j,,i]+(runif(d@M)<p@pMut)*rnorm(d@M)
      d@Pi[j,,i] <- march.ea.checkBounds(v=d@Pi[j,,i],lb=0,ub=1)
      d@Pi[j,,i] <- march.dcmm.ea.h.scale(d@Pi[j,,i],d@M)
    }
  }
  
  #Build of AProbT, RA, A, CProbT, RC and C
  
	if(d@Amodel=="complete"){
		mccov <- march.mccov.sk(d@M,max(d@orderHC,1),AtmCovar,d@y@Kcov,placeACovar,NbAMCovar,matrix(d@AQ[1,,],d@M^max(d@orderHC,1),d@M),d@ATCovar)
    d@AProbT <- mccov$ProbT
    
    l <- GMTD_tm_cov(max(d@orderHC,1),d@M,d@APhi,d@AProbT)
    d@A <- l$HOQ
	}else{
		ANSS<-matrix(0,d@M^(d@orderHC+1)*AtmCovar,1)
		AValT<-BuildArrayCombinations(d@M,d@orderHC,d@y@Kcov[placeACovar],NbAMCovar)

    d@AProbT <- BuildArrayQ(d@M,l=d@orderHC,i0_il=AValT,n_i0_il = ANSS,q=d@AQ,kcov=d@y@Kcov[placeACovar],ncov=NbAMCovar,S=d@ATCovar)
    l<-GMTD_tm_cov(d@orderHC,d@M,d@APhi,array(d@AProbT,c(d@M^(d@orderHC+1)*AtmCovar,d@orderHC+NbAMCovar)))
    d@A<-l$HOQ
  }
    
  #VISIBLE LEVEL
  if(d@Cmodel=="complete"){
		for(state in 1:d@M){
			CTCovartmp <- list()
			if(NbCMCovar>0){
        for(i in 1:NbCMCovar){
          CTCovartmp[[i]] <- d@CTCovar[[i]][,,state]
        }
      }

      mccov <- march.mccov.sk(d@y@K,d@orderVC,CtmCovar,d@y@Kcov,placeCCovar,NbCMCovar,matrix(d@CQ[1,,,state],d@y@K^d@orderVC,d@y@K),CTCovartmp)
      d@CProbT[,,state] <- mccov$ProbT
      d@RB[,,state] <- GMTD_tm_cov(d@orderVC,d@y@K,t(d@CPhi[1,,state]),matrix(d@CProbT[,,state],d@y@K^(d@orderVC+1)*CtmCovar,1+NbCMCovar))$CRHOQ
    }
	}else{
    d@CProbT <- array(0,c(d@y@K^(d@orderVC+1)*CtmCovar,d@orderVC+NbCMCovar,d@M))
    CValT <- BuildArrayCombinations(d@y@K,d@orderVC,d@y@Kcov[placeCCovar],NbCMCovar)
    CNSS <- BuildArrayNumberOfSequencesDCMM(d@y,d@orderVC,d@CMCovar,NbCMCovar,placeCCovar)
    	
    if(d@Cmodel=="mtd"){
    	t <- array(d@CQ[,,,1],c(1,d@y@K,d@y@K))
    }else{
      t <- array(d@CQ[,,,1],c(d@orderVC,d@y@K,d@y@K))
    }
    q_i0_il <- array(0,c(d@y@K^(d@orderVC+1)*CtmCovar,d@orderVC+NbCMCovar))
    	
    	
    for(state in 1:d@M){
      CTCovartmp <- list()
      if(NbCMCovar>0){
      	for(i in 1:NbCMCovar){
      		CTCovartmp[[i]] <- d@CTCovar[[i]][,,state]
      	}
      }
      d@CProbT[,,state] <- BuildArrayQ(d@y@K,l=d@orderVC,i0_il=CValT,n_i0_il = CNSS,q=t,kcov=d@y@Kcov[placeCCovar],NbCMCovar,S=CTCovartmp)
      tm <- GMTD_tm_cov(d@orderVC,d@y@K,as.matrix(t(d@CPhi[,,state])),array(d@CProbT[,,state],c(d@y@K^(d@orderVC+1)*CtmCovar,d@orderVC+NbCMCovar)))
      d@RB[,,state] <- tm$CRHOQ
    }
  }
  d@ll <- march.dcmm.h.cov.computeLL(d,d@y)
  d
}

march.dcmm.cov.ea.crossover <- function(d1, d2, AConst){
  
	  NbAMCovar <- sum(d1@AMCovar)
 	  NbCMCovar <- sum(d1@CMCovar)
  
  	placeACovar <- which(d1@AMCovar==1)
  	placeCCovar <- which(d1@CMCovar==1)
    
   	AtmCovar <- 1
  	if(NbAMCovar>0){
  		for (i in 1:NbAMCovar){
      		AtmCovar <- AtmCovar*d1@y@Kcov[placeACovar[i]]
      	}
  	}

  
  	c1 <- march.dcmm.cov.constructEmptyDcmm(d1@M,d1@y,d1@orderVC,d1@orderHC,d1@AMCovar,d1@CMCovar,d1@Amodel,d1@Cmodel)			
  	c2 <- march.dcmm.cov.constructEmptyDcmm(d1@M,d1@y,d1@orderVC,d1@orderHC,d1@AMCovar,d1@CMCovar,d1@Amodel,d1@Cmodel)
  	
    if(AConst==FALSE){
  	  if(d1@Amodel=="mtdg"){
  		  for(i in 1:d1@orderHC){
  			  for(j in 1:d1@M){
  				  r <- march.ea.h.sbx(d1@AQ[i,j,],d2@AQ[i,j,],2)
  				  c1@AQ[i,j,] <- r$c1
  				  c2@AQ[i,j,] <- r$c2
  			  }
  		  }
  	  }else{
  		  if(d1@orderHC>0){
  			  for(i in 1:dim(d1@AQ)[2]){
  				  r <- march.ea.h.sbx(d1@AQ[1,i,],d2@AQ[1,i,],2)
  				  c1@AQ[1,i,] <- r$c1
  				  c2@AQ[1,i,] <- r$c2
  			  }
		    }else{
			    r <- march.ea.h.sbx(d1@AQ[1,1,],d2@AQ[1,1,],2)
			    for(i in 1:d1@M){
				    c1@AQ[1,i,] <- r$c1
				    c2@AQ[1,i,] <- r$c2
			    }
		    }
 	    } 
    }else{
      c1@AQ <- array(diag(d1@M),c(1,d1@M,d1@M))
      c2@AQ <- array(diag(d1@M),c(1,d1@M,d1@M))
    }
 	
  	if(d1@Cmodel=="mtdg"){
  		for(i in 1:d1@orderVC){
  			for(j in 1:d1@y@K){
  				for(k in 1:d1@M){
  					r <- march.ea.h.sbx(d1@CQ[i,j,,k],d2@CQ[i,j,,k],2)
  					c1@CQ[i,j,,k] <- r$c1
  					c2@CQ[i,j,,k] <- r$c2
  				}
  			}
  		}
  	}else{
  		for(j in 1:dim(d1@CQ)[2]){
  			for(k in 1:d1@M){
  				r <- march.ea.h.sbx(d1@CQ[1,j,,k],d2@CQ[1,j,,k],2)
  				c1@CQ[1,j,,k] <- r$c1
  				c2@CQ[1,j,,k] <- r$c2
  			}
  		}
	}  
	
	
  	# Crossover the initial probability
  	if(d1@orderHC>0){
  		for( i in 1:d1@orderHC ){
    		for( j in 1:(d1@M^(i-1)*AtmCovar)){
      			r <- march.ea.h.sbx(d1@Pi[j,,i],d2@Pi[j,,i],2)
      			c1@Pi[j,,i] <- r$c1
      			c2@Pi[j,,i] <- r$c2
      		}
    	}
  	}else{
  		for(i in 1:AtmCovar){
  			r <- march.ea.h.sbx(d1@Pi[1,,1],d2@Pi[1,,1],2)
  			c1@Pi[i,,1] <- r$c1
  			c2@Pi[i,,1] <- r$c2
  		}
  	}
  
    r <- march.ea.h.sbx(d1@APhi[1,],d2@APhi[1,],2)
    c1@APhi[1,] <- r$c1
    c2@APhi[1,] <- r$c2
  	
  	if(NbCMCovar>0|d1@orderVC>0){
  		for(i in 1:d1@M){
      		r <- march.ea.h.sbx(d1@CPhi[1,,i],d2@CPhi[1,,i],2)
      		c1@CPhi[1,,i] <- r$c1
      		c2@CPhi[1,,i] <- r$c2
      	}
  	}
	
  
  	if(sum(d1@AMCovar)>0){
    	for(i in 1:sum(d1@AMCovar)){
      		for(j in 1:d1@y@Kcov[placeACovar[i]]){
          		r <- march.ea.h.sbx(d1@ATCovar[[i]][j,],d2@ATCovar[[i]][j,],2)
          		c1@ATCovar[[i]][j,] <- r$c1
          		c2@ATCovar[[i]][j,] <- r$c2
      		}
    	}
  	}
  
  	if(sum(d1@CMCovar)>0){
    	for(i in 1:sum(d1@CMCovar)){
      		for(j in 1:d1@y@Kcov[placeCCovar[i]]){
          		for(l in 1:d1@M){
            		r <- march.ea.h.sbx(d1@CTCovar[[i]][j,,l],d2@CTCovar[[i]][j,,l],2)
            		c1@CTCovar[[i]][j,,l] <- r$c1
            		c2@CTCovar[[i]][j,,l] <- r$c2
          		}
            }
    	}
  	}
  	list(c1=c1,c2=c2)
}


march.dcmm.cov.ea.optimizing<-function(d,p){
	ConstAPhi <- 1
	ConstCPhi <- 1
	placeACovar <- which(d@AMCovar==1)
	placeCCovar <- which(d@CMCovar==1)
	maxAMKCovar <- 1
	maxCMKCovar <- 1
	if(sum(d@AMCovar)>0){
		maxAMKCovar <- max(d@y@Kcov[placeACovar])
	}
	if(sum(d@CMCovar)>0){
		maxCMKCovar <- max(d@y@Kcov[placeCCovar])
	}
	
	#Initialisation of Delta for the hidden level
	if(d@Amodel=="complete"){
		if(sum(d@AMCovar)==0){
			AMDelta <- 0.1
		}else{
			AMDelta <- array(0.1,c(1,d@M^max(d@orderHC,1)+sum(d@AMCovar)*maxAMKCovar+1))
		}
	}else if(d@Amodel=="mtd"){
		AMDelta <- array(0.1,c(1,d@M+sum(d@AMCovar)*maxAMKCovar+1))
	}else{
		AMDelta <- array(0.1,c(d@orderHC,d@M+sum(d@AMCovar)*maxAMKCovar+1))
	}
	
	#Initialisation of Delta for the visible level
	if(d@Cmodel=="complete"){
		if(sum(d@CMCovar)==0){
			CMDelta <- 0.1
		}else{
			CMDelta <- array(0.1,c(d@M,d@y@K^d@orderVC+sum(d@CMCovar)*maxCMKCovar+1))
		}
	}else if(d@Cmodel=="mtd"){
		CMDelta <- array(0.1,c(d@M,d@y@K+sum(d@CMCovar)*maxCMKCovar+1))
	}else{
		CMDelta <- array(0.1,c(d@M,d@y@K*d@orderVC+sum(d@CMCovar)*maxCMKCovar+1))
	}
	
	#Baum-Welch algorithm
	for(i in march.h.seq(1,p@iterBw)){
		referenceLL <- d@ll
		l<- march.dcmm.cov.em(d,AMDelta,CMDelta,ConstAPhi,ConstCPhi)
		d<-l$d
		AMDelta<-l$AMDelta
		CMDelta<-l$CMDelta
		if(abs(d@ll-referenceLL)<=p@stopBw){ break }
		
	}
	d
}

march.dcmm.ea.h.scale <- function(v,j){
  mA <- sum(v)
  
  if( mA==0 ){ v <- matrix(1,ncol=1,nrow=j)*(1/j) }
  else{ v <- v/mA }
  v
}

