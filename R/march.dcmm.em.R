march.dcmm.cov.em <- function(d,AMDelta,CMDelta,ConstAPhi,ConstCPhi){
	
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
  	
	CtmCovar <- 1
  	if(NbCMCovar>0){
    	for (i in 1:NbCMCovar){
      		CtmCovar <- CtmCovar*d@y@Kcov[placeCCovar[i]]
    	}
  	}
  
  
  	SAlpha <- array(0,c(d@M^max(d@orderHC,1), max(d@y@T), d@y@N))
  	SBeta <- array(0,c(d@M^max(d@orderHC,1), max(d@y@T), d@y@N))
  	SAlog <- array(0,c(max(d@y@T),d@y@N))
  	SBlog <- array(0,c(max(d@y@T),d@y@N))
  	LLAlpha <- rep(0,d@y@N)
  	Epsilon <- array(0,c(d@M^max(d@orderHC,1)*AtmCovar,d@M,max(d@y@T)-1,d@y@N))
  	Gamma <- array(0,c(d@M^max(d@orderHC,1)*AtmCovar,max(d@y@T),d@y@N))
  	
	for(n in 1:d@y@N){
		#Extract the current sequence
		s <- march.dataset.h.extractSequence(d@y,n)
		
		#Forward-Backward algorithm
		a <- march.dcmm.fp.cov(d,s,d@y@cov[n,,],AtmCovar,CtmCovar)
		SAlpha[,1:d@y@T[n],n] <- a$SAlpha
		SAlog[1:d@y@T[n],n] <- a$SAlog
		LLAlpha[n] <- a$LLAlpha

		b <- march.dcmm.bp.cov(d,s,d@y@cov[n,,],AtmCovar,CtmCovar)
		SBeta[,1:d@y@T[n],n] <- b$SBeta
		SBlog[1:d@y@T[n],n] <- b$SBlog
		
		#Computation of epsilon and gamma
		Epsilon[,,1:d@y@T[n]-1,n] <- march.dcmm.eps.cov(d,s,d@y@cov[n,,],a$SAlpha,b$SBeta,a$SAlog,b$SBlog,a$LLAlpha,AtmCovar,CtmCovar)
		Gamma[,1:d@y@T[n],n] <- march.dcmm.gam.cov(d,s,d@y@cov[n,,],a$SAlpha,b$SBeta,a$SAlog,b$SBlog,a$LLAlpha,AtmCovar,CtmCovar)
	}

	#Computation of the log-likelihood
	d@ll <- sum(LLAlpha)
	
	##Reestimation of pi, A and C when no mixture modeling and covariates are involved
	#We save the current model to compare if we will have an increase in the log-likelihood
	dtmp <- d
	dtmp <- march.dcmm.bw.cov(d,AtmCovar,CtmCovar,Epsilon,Gamma)
		
	#Computation of the likelihood
	SAlphat <- array(0,c(d@M^max(d@orderHC,1), max(d@y@T), d@y@N))
  	SAlogt <- array(0,c(max(d@y@T),d@y@N))
  	LLAlphat <- rep(0,d@y@N)
  	
  	for(n in 1:d@y@N){
  		s <- march.dataset.h.extractSequence(d@y,n)
		
		a <- march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)
		SAlphat[,1:dtmp@y@T[n],n] <- a$SAlpha
		SAlogt[1:dtmp@y@T[n],n] <- a$SAlog
		LLAlphat[n] <- a$LLAlpha
  	}
	dtmp@ll <- sum(LLAlphat)
	if(dtmp@ll>d@ll){
		d <- dtmp
		SAlpha <- SAlphat
		SAlog <- SAlogt
		LLAlpha <- LLAlphat
	}
	
	##########Reestimation of A when a mixture modeling or some covariates are involved###########
	if(d@M>1 & ((d@Amodel=="complete" & AtmCovar>1) | d@Amodel=="mtd" | d@Amodel=="mtdg")){
		for(n in 1:d@y@N){
			s <- march.dataset.h.extractSequence(d@y,n)
			b <- march.dcmm.bp.cov(d,s,d@y@cov[n,,],AtmCovar,CtmCovar)
			SBeta[,1:d@y@T[n],n] <- b$SBeta
			SBlog[1:d@y@T[n],n] <- b$SBlog
			
			Epsilon[,,1:d@y@T[n]-1,n] <- march.dcmm.eps.cov(d,s,d@y@cov[n,,],SAlpha[,1:s@N,n],b$SBeta,SAlog[1:s@N],b$SBlog,LLAlpha[n],AtmCovar,CtmCovar)

		}
	
	
		#Computation of ANSS (number of appeareace of each combination of lags and covariate in the hidden process)
		ANSS <- rep(0,d@M^(max(d@orderHC,1)+1)*AtmCovar)
		for(n in 1:d@y@N){
			for(p in march.h.seq(d@orderHC+d@orderVC, d@y@T[n]-1)){
				ANSS <- ANSS+as.vector(Epsilon[,,p,n])
			}
		}
		
		AValT <- BuildArrayCombinations(d@M,max(d@orderHC,1),d@y@Kcov[placeACovar],NbAMCovar)
    
    	#Reestimation of APhi
		if(d@Amodel=="complete" & NbAMCovar>0){
			
			#Gradient of the log-likelihood with respect to APhi
			SDip <- GMTD.dp(1+NbAMCovar,ANSS,d@AProbT,d@APhi)
			
			#Optimisation of APhi
			l <- GMTD.op.dcmm.A(d,ANSS,SDip,AMDelta[1,1],0.0001,ConstAPhi,AtmCovar,CtmCovar)
			d <- l$d
			AMDelta[1,1] <- l$MD
			
		}else if(d@Amodel=="mtd" | d@Amodel=="mtdg"){
			
			SDip <- GMTD.dp(max(d@orderHC,1)+NbAMCovar,ANSS,d@AProbT,d@APhi)
		
			l <- GMTD.op.dcmm.A(d,ANSS,SDip,AMDelta[1,1],0.0001,ConstAPhi,AtmCovar,CtmCovar)
			d <- l$d
			AMDelta[1,1] <- l$MD
		}
	
		#Reestimation of AQ
		if(d@Amodel=="complete" & NbAMCovar>0){
			for(g in 1:(d@M^max(d@orderHC,1))){
				
				SDiq <- march.mccov.dq(max(d@orderHC,1),d@M,ANSS,d@AProbT,AValT,d@APhi)
				
				l <- GMTD.oq.dcmm.A.b(d,AtmCovar,CtmCovar,AValT,ANSS,g,SDiq,AMDelta[1,g+1],0.0001)
				d <- l$d
				AMDelta[1,g+1] <- l$MD
				if(d@orderHC==0){
					break
				}
			}
		}else if(d@Amodel=="mtd"){
			for(g in 1:d@M){
				
				SDiq <- GMTD.dq(d@orderHC,d@M,ANSS,d@AProbT,AValT,d@APhi)
				
				l <- GMTD.oq.dcmm.A(d,AtmCovar,CtmCovar,AValT,ANSS,g,SDiq,AMDelta[1,g+1],0.0001)
				d <- l$d
				AMDelta[1,g+1] <- l$MD
				
			}
		}else if(d@Amodel=="mtdg"){
			for(g in 1:d@orderHC){
				if(d@APhi[g]!=0){
					for(r in 1:d@M){
						SDiq <- GMTDg.dq(d@orderHC,d@M,ANSS,d@AProbT,AValT,d@APhi)
						
						l <- GMTDg.oq.dcmm.A(d,AtmCovar,CtmCovar,AValT,ANSS,r,g,SDiq,AMDelta[g,r+1],0.0001)
						
						d <- l$d
						AMDelta[g,r+1] <- l$MD
					}
				}
			}
		}
		
		#Reestimation of ATCovar
		if(sum(d@AMCovar)>0){
			if(d@Amodel=="complete"){
				CCov <- 1
				Delcov <- d@M^d@orderHC
			}else{
				CCov <- d@orderHC
				Delcov <- d@M
			}
			
			for(i in 1:sum(d@AMCovar)){
				SDicov <- GMTD.dcov(CCov+i,max(d@orderHC,1)+1+i,d@M,d@y@Kcov[placeACovar[i]],ANSS,d@AProbT,AValT,d@APhi)
				
				for(g in 1:d@y@Kcov[placeACovar[i]]){
				
					l <- GMTD.ocov.dcmm.A(d,SDicov,AtmCovar,CtmCovar,AValT,ANSS,i,g,CCov+i,AMDelta[1,1+Delcov+(i-1)*max(d@y@Kcov[placeACovar])+g],0.0001)
					d <- l$d
					AMDelta[1,1+Delcov+(i-1)*max(d@y@Kcov)+g] <- l$MD
				}
			}
		}
	}
		
		
	##############Reestimation of the visible process#############
	if((d@Cmodel=="visible"&CtmCovar>1)|d@Cmodel=="mtd"|d@Cmodel=="mtdg"){
		
		for(n in 1:d@y@N){
			#Extract the current sequence
			s <- march.dataset.h.extractSequence(d@y,n)
		
			#Forward-Backward algorithm
			a <- march.dcmm.fp.cov(d,s,d@y@cov[n,,],AtmCovar,CtmCovar)
			SAlpha[,1:d@y@T[n],n] <- a$SAlpha
			SAlog[1:d@y@T[n],n] <- a$SAlog
			LLAlpha[n] <- a$LLAlpha

			b <- march.dcmm.bp.cov(d,s,d@y@cov[n,,],AtmCovar,CtmCovar)
			SBeta[,1:d@y@T[n],n] <- b$SBeta
			SBlog[1:d@y@T[n],n] <- b$SBlog
		
			#Computation of epsilon and gamma
			Epsilon[,,1:d@y@T[n]-1,n] <- march.dcmm.eps.cov(d,s,d@y@cov[n,,],a$SAlpha,b$SBeta,a$SAlog,b$SBlog,a$LLAlpha,AtmCovar,CtmCovar)
			Gamma[,1:d@y@T[n],n] <- march.dcmm.gam.cov(d,s,d@y@cov[n,,],a$SAlpha,b$SBeta,a$SAlog,b$SBlog,a$LLAlpha,AtmCovar,CtmCovar)
		}
		#Computation of the log-likelihood
		d@ll <- sum(LLAlpha)
		
		
		CValT <- BuildArrayCombinations(d@y@K,d@orderVC,d@y@Kcov[placeCCovar],NbCMCovar)
		
		CNSS <- array(0,c(d@y@K^(d@orderVC+1)*CtmCovar,d@M))
		
		if(d@M==1){
			CNSS <- matrix(BuildArrayNumberOfSequencesDCMM(d@y,d@orderVC,d@CMCovar,NbCMCovar,placeCCovar),d@y@K^(d@orderVC+1)*CtmCovar,1)
		}else{
			for(state in 1:d@M){
				for(n in 1:d@y@N){
					for(t in march.h.seq(d@orderVC+1,d@y@T[n]-1)){
						Pstate <- sum(Epsilon[,state,t,n])
						#pos <- d@y@yRaw[n,(t-d@orderVC+1):(t+1)]
						pos <- d@y@y[[n]][(t-d@orderVC+1):(t+1)]
						rowC <- pos[1]
						for(r in march.h.seq(2,d@orderVC+1)){
							rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
						}
						
						rowCovarC <- 1
						if(CtmCovar>1){
							r <- 0
							tKC <- 1
							for(i in d@y@Ncov:1){
								if(d@CMCovar[i]>0){
									r <- r+1
									rowCovarC <- rowCovarC+tKC*(d@y@cov[n,t+1,i]-1)
									tKC <- tKC*d@y@Kcov[i]
								}
							}
						}
						CNSS[(rowC-1)*CtmCovar+rowCovarC,state] <- CNSS[(rowC-1)*CtmCovar+rowCovarC,state]+Pstate
					}
				}
			}
		}
	
		#For each hidden state, we reestimate the corresponding visible process
		for(state in 1:d@M){
			
			#Reestimation of CPhi
			if(d@Cmodel=="complete" & NbCMCovar>0){
				SDip <- GMTD.dp(1+NbCMCovar,CNSS[,state],matrix(d@CProbT[,,state],dim(d@CProbT)[1],dim(d@CProbT)[2]),d@CPhi[,,state])
				l <- GMTD.op.dcmm.C(d,AtmCovar,CtmCovar,CNSS[,state],SDip,state,CMDelta[state,1],0.0001,ConstCPhi)
				d <- l$d
				CMDelta[state,1] <- l$MD
				
			}else if(d@Cmodel=="mtd"|d@Cmodel=="mtdg"){
				SDip <- GMTD.dp(d@orderVC+NbCMCovar,CNSS[,state],matrix(d@CProbT[,,state],dim(d@CProbT)[1],dim(d@CProbT)[2]),d@CPhi[,,state])
				l <- GMTD.op.dcmm.C(d,AtmCovar,CtmCovar,CNSS[,state],SDip,state,CMDelta[state,1],0.0001,ConstCPhi)
				d <- l$d
				CMDelta[state,1] <- l$MD
			}
			
			
			#Reestimation of CQ
			if(d@Cmodel=="complete" & NbCMCovar>0){
				for(g in 1:d@y@K^d@orderVC){
					SDiq <- march.mccov.dq(d@orderVC,d@y@K,CNSS[,state],matrix(d@CProbT[,,state],dim(d@CProbT)[1],dim(d@CProbT)[2]),CValT,d@CPhi[,,state])
					l <- GMTD.oq.dcmm.C.b(d,AtmCovar,CtmCovar,CValT,CNSS[,state],g,state,SDiq,CMDelta[state,g+1],0.0001)
					d <- l$d
					CMDelta[state,g+1] <- l$MD
				}
			}else if(d@Cmodel=="mtd"){
				for(g in 1:d@y@K){
					SDiq <- GMTD.dq(d@orderVC,d@y@K,CNSS[,state],matrix(d@CProbT[,,state],dim(d@CProbT)[1],dim(d@CProbT)[2]),CValT,d@CPhi[,,state])
					l <- GMTD.oq.dcmm.C(d,AtmCovar,CtmCovar,CValT,CNSS[,state],g,state,SDiq,CMDelta[state,g+1],0.0001)
					d <- l$d
					CMDelta[state,g+1] <- l$MD
				}
			}else if(d@Cmodel=="mtdg"){
				for(g in 1:d@orderVC){
					if(d@CPhi[1,g,state]!=0){
						for(r in 1:d@y@K){
							SDiq <- GMTDg.dq(d@orderVC,d@y@K,CNSS[,state],matrix(d@CProbT[,,state],dim(d@CProbT)[1],dim(d@CProbT)[2]),CValT,d@CPhi[,,state])
							
							l <- GMTDg.oq.dcmm.C(d,AtmCovar,CtmCovar,CValT,CNSS[,state],g,r,state,SDiq,CMDelta[state,1+(g-1)*d@y@K+r],0.0001)
							d <- l$d
							CMDelta[state,1+(g-1)*d@y@K+r] <- l$MD
						}
					}
				}
			}
			
			
			#Reestimation of CTCovar
			if(NbCMCovar>0){
				
				#Computation of CCov and Delcov
				if(d@Cmodel=="complete"){
					CCov <- 1
					Delcov <- d@y@K^d@orderVC
				}else if(d@Cmodel=="mtd"){
					if(d@orderVC==0){
						CCov <- 1
					}else{
						CCov <- d@orderVC
					}
					Delcov <- d@y@K
				}else{
					if(d@orderVC==0){
						CCov <- 1
					}else{
						CCov <- d@orderVC
					}
					Delcov <- d@y@K*d@orderVC
				}
				
				for(i in 1:NbCMCovar){
					SDicov <- GMTD.dcov(CCov+i,d@orderVC+1+i,d@y@K,d@y@Kcov[placeCCovar[i]],CNSS[,state],matrix(d@CProbT[,,state],dim(d@CProbT)[1],dim(d@CProbT)[2]),CValT,d@CPhi[1,,state])
					
					for( g in 1:d@y@Kcov[placeCCovar[i]]){
						l <- GMTD.ocov.dcmm.C(d,SDicov,AtmCovar,CtmCovar,CValT,CNSS[,state],i,g,state,CCov+i,CMDelta[state,1+Delcov+(i-1)*max(d@y@Kcov[placeCCovar])+g],0.0001)
						
						d <- l$d
						CMDelta[state,1+Delcov+(i-1)*max(d@y@Kcov[placeCCovar])+g] <- l$MD
					}
				}
			}
			if(d@Cmodel=="mtd"|d@Cmodel=="mtdg"|NbCMCovar>0){
				l <- GMTD_tm_cov(d@orderVC,d@y@K,matrix(dtmp@CPhi[,,state],1,dim(dtmp@CPhi)[2]),matrix(d@CProbT[,,state],dim(d@CProbT)[1],dim(d@CProbT)[2]))
				d@RB[,,state] <- l$CRHOQ
			}
		}
	}
	
	#Computation of the log-likelihood
	tLL <- 0
	for(n in 1:d@y@N){
		s <- march.dataset.h.extractSequence(d@y,n)
		tLL <- tLL+march.dcmm.fp.cov(d,s,d@y@cov[n,1:s@N,],AtmCovar,CtmCovar)$LLAlpha
	}
	d@ll <- tLL
	
	list(d=d,AMDelta=AMDelta,CMDelta=CMDelta)
}


march.dcmm.fp.cov<-function(d,s,Covar,AtmCovar,CtmCovar){
  
	ordHC <- d@orderHC
  	if(ordHC==0){ordHC <- 1}
  
	SAlpha <- matrix(0,d@M^ordHC,s@N)
	SAlog <- rep(0,s@N)
  
  
	# First part t=d@orderVC+1
  	pos <- s@y[1:d@orderVC]
	rowC <- 1
	for(r in march.h.seq(1,d@orderVC)){
		rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
	}
  
	rowCovar <- 1
	if(AtmCovar>1){
		r <- 0
		tKC <- 1
		for(i in d@y@Ncov:1){
			if(d@AMCovar[i]>0){
				r <- r+1
				rowCovar <- rowCovar+tKC*(Covar[d@orderVC+1,i]-1)
				tKC <- tKC*d@y@Kcov[i]
      		}
    	}
	}
  
	rowCovarC <- 1
	if(CtmCovar>1){
		r <- 0
		tKC <- 1
		for(i in d@y@Ncov:1){
			if(d@CMCovar[i]>0){
				r <- r+1
				rowCovarC <- rowCovarC+tKC*(Covar[d@orderVC+1,i]-1)
				tKC <- tKC*d@y@Kcov[i]
			}
		}
	}
  
	for(j in 1:d@M){
		SAlpha[j,d@orderVC+1] <- d@Pi[rowCovar,j,1]*d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[d@orderVC+1],j]
	}
  
	mA <- sum(SAlpha[1:d@M,d@orderVC+1])/d@M
  
	if(mA==0){
		SAlpha[1:d@M,d@orderVC+1] <- rep(1,d@M)
		SAlog[d@orderVC+1] <- log(1/d@M)
	}else{
		SAlpha[1:d@M,d@orderVC+1] <- SAlpha[1:d@M,d@orderVC+1]/mA
		SAlog[d@orderVC+1] <- log(mA)
	}
  
  
	# Second part t takes value between orderVC+2 and orderVC+orderHC
 
	if(ordHC>1){
		for(t in march.h.seq(d@orderVC+2,d@orderVC+ordHC)){
			pos <- s@y[(t-d@orderVC):(t-1)]
			rowC <- 1
			for(r in march.h.seq(1,d@orderVC)){
        		rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
      		}
      		
      		rowCovar <- 1
      		if(AtmCovar>1){
        		r <- 0
        		tKC <- 1
        		for(i in d@y@Ncov:1){
          			if(d@AMCovar[i]>0){
          				r <- r+1
          				rowCovar <- rowCovar+tKC*(Covar[t,i]-1)
          				tKC <- tKC*d@y@Kcov[i]
          			}
        		}
			}
			
    		rowCovarC <- 1
    		if(CtmCovar>1){
      			r <- 0
      			tKC <- 1
      			for(i in d@y@Ncov:1){
        			if(d@CMCovar[i]>0){
        				r <- r+1
        				rowCovarC <- rowCovarC+tKC*(Covar[t,i]-1)
        				tKC <- tKC*d@y@Kcov[i]
        			}
      			}
			}
    
    		for(i in 1:d@M^(t-d@orderVC-1)){
      			for(j in 1:d@M){
					row <- d@M^(t-d@orderVC-1)*(j-1)+i
        			SAlpha[row,t] <- d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t],j]*d@Pi[AtmCovar*(i-1)+rowCovar,j,t-d@orderVC]*SAlpha[i,t-1]
      			}
    		}
    
    		mA <- sum(SAlpha[1:d@M^(t-d@orderVC),t])/(d@M^(t-d@orderVC))
    
    		if(mA==0){
      			SAlpha[1:d@M^(t-d@orderVC),t] <- rep(1,d@M^(t-d@orderVC))
      			SAlog[t] <- log(1/(d@M^(t-d@orderVC)))
    		}else{
				SAlpha[1:d@M^(t-d@orderVC),t] <- SAlpha[1:d@M^(t-d@orderVC),t]/mA
				SAlog[t] <- log(mA)
			}
		}
	}
  
  
	# Third part t>orderVC+orderHC
  
	for(t in march.h.seq(d@orderVC+ordHC,s@N-1)){
		pos <- s@y[(t-d@orderVC+1):t]
		rowC <- 1
		for(r in march.h.seq(1,d@orderVC)){
			rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
		}
		
		rowCovar <- 1
		if(AtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@AMCovar[i]>0){
					r <- r+1
					rowCovar <- rowCovar+tKC*(Covar[t+1,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		rowCovarC <- 1
		if(CtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@CMCovar[i]>0){
					r <- r+1
					rowCovarC <- rowCovarC+tKC*(Covar[t+1,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		for(j in 1:d@M^ordHC){
			j0 <- floor((j-1)/(d@M^(ordHC-1)))+1
			for(k in 1:d@M^ordHC){
				SAlpha[j,t+1] <- SAlpha[j,t+1]+SAlpha[k,t]*d@A[AtmCovar*(k-1)+rowCovar,j]*d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t+1],j0]
			}
		}
    
		mA <- sum(SAlpha[,t+1])/(d@M^ordHC)
    
		if(mA==0){
			SAlpha[,t+1] <- rep(1,d@M^ordHC)
			SAlog[t+1] <- log(1/d@M^ordHC)
		}else{
			SAlpha[,t+1] <- SAlpha[,t+1]/mA
			SAlog[t+1] <- log(mA)
		}
	}
  
  
	# Computation of the log-likelihood
	LLAlpha <- log(sum(SAlpha[,s@N]))+sum(SAlog)
  
  
	list(SAlog=SAlog,LLAlpha=LLAlpha,SAlpha=SAlpha)
}


march.dcmm.bp.cov<-function(d,s,Covar,AtmCovar,CtmCovar){
	
	ordHC <- d@orderHC
	if(ordHC==0){ordHC <- 1}

	RA <- march.dcmm.cov.h.compactA(d)


	SBeta <- array(0,c(d@M^ordHC,s@N))
	SBlog <- rep(0,s@N)

	#Part t=s@N
	SBeta[,s@N] <- rep(1,d@M^ordHC)


	#Part t between orderVC+orderHC and s@N

	if(ordHC<s@N){
		for(t in march.h.seq(s@N-1,d@orderVC+ordHC,-1)){
			pos <- s@y[(t-d@orderVC+1):t]
			rowC <- 1
			for (r in march.h.seq(1,d@orderVC)){
				rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
			}
		
			rowCovar <- 1
			if(AtmCovar>1){
				r <- 0
				tKC <- 1
				for(i in d@y@Ncov:1){
					if(d@AMCovar[i]>0){
						r <- r+1
						rowCovar <- rowCovar+tKC*(Covar[t+1,i]-1)
						tKC <- tKC*d@y@Kcov[i]
					}
				}
			}
		
			rowCovarC <- 1
			if(CtmCovar>1){
				r <- 0
				tKC <- 1
				for(i in d@y@Ncov:1){
					if(d@CMCovar[i]>0){
						r <- r+1
						rowCovarC <- rowCovarC+tKC*(Covar[t+1,i]-1)
						tKC <- tKC*d@y@Kcov[i]
					}
				}
			}
		
			for(i in 1:(d@M^ordHC)){
				for(j in 1:d@M){
					col <- d@M^(ordHC-1)*(j-1)+floor((i-1)/d@M)+1
					SBeta[i,t] <- SBeta[i,t]+RA[AtmCovar*(i-1)+rowCovar,j]*d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t+1],j]*SBeta[col,t+1]
				}
			}
		
			mB <- sum(SBeta[,t])/(d@M^ordHC)
		
			if(mB==0){
				SBeta[,t] <- rep(1,d@M^ordHC)
				SBlog[t] <- log(1/(d@M^ordHC))
			}else{
				SBeta[,t] <- SBeta[,t]/mB
				SBlog[t] <- log(mB)
			}
		}
	}


	# Part t between orderVC and orderVC + orderHC

	if(ordHC>1){
		for(t in march.h.seq(d@orderVC+ordHC-1,d@orderVC+1,-1)){
			pos <- s@y[(t-d@orderVC+1):t]
			rowC <- 1
			for(r in march.h.seq(1,d@orderVC)){
				rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
			}
		
			rowCovar <- 1
			if(AtmCovar>1){
				r <- 0
				tKC <- 1
				for(i in d@y@Ncov:1){
					if(d@AMCovar[i]>0){
						r <- r+1
						rowCovar <- rowCovar+tKC*(Covar[t+1,i]-1)
						tKC <- tKC*d@y@Kcov[i]
					}
				}
			}
		
			rowCovarC <- 1
			if(CtmCovar>1){
				r <- 0
				tKC <- 1
				for(i in d@y@Ncov:1){
					if(d@CMCovar[i]>0){
						r <- r+1
						rowCovarC <- rowCovarC+tKC*(Covar[t+1,i]-1)
						tKC <- tKC*d@y@Kcov[i]
					}
				}
			}
		
			for(i in 1:(d@M^(t-d@orderVC))){
				for(j in 1:d@M){
					col <- d@M^(t-d@orderVC)*(j-1)+i
					SBeta[i,t] <- SBeta[i,t]+d@Pi[AtmCovar*(i-1)+rowCovar,j,t-d@orderVC+1]*d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t+1],j]*SBeta[col,t+1]
				}
			}
		
			mB <- sum(SBeta[,t])/(d@M^(t-d@orderVC))
		
			if(mB==0){
				SBeta[1:(d@M^(t-d@orderVC)),t] <- rep(1,d@M^(t-d@orderVC))
				SBlog[t] <- log(1/(d@M^(t-d@orderVC)))
			}else{
				SBeta[,t] <- SBeta[,t]/mB
				SBlog[t] <- log(mB)
			}	
		}
  	}
  	list(SBeta=SBeta,SBlog=SBlog)
}

march.dcmm.eps.cov<-function(d,s,Covar,SAlpha,SBeta,SAlog,SBlog,LLAlpha,AtmCovar,CtmCovar){
	
	#Init
	tt <- 0
	ordHC <- d@orderHC
	if(ordHC==0){
		ordHC <- 1
		tt <- 1
	}
	Epsilon <- array(0,c(d@M^ordHC*AtmCovar,d@M,s@N-1))
	ZEpsilon <- array(1,c(d@M^ordHC*AtmCovar,d@M,s@N-1))
	
	RA <- march.dcmm.cov.h.compactA(d)
	
	
	#Search for the place of zeros in Epsilon
	for(t in march.h.seq(d@orderVC+1,s@N-1)){
		pos <- s@y[(t-d@orderVC+1):t]
		rowC <- 1
		for( r in march.h.seq(1,d@orderVC)){
			rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
		}
		
		rowCovar <- 1
		if(AtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@AMCovar[i]>0){
					r <- r+1
					rowCovar <- rowCovar+tKC*(Covar[t+1,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		rowCovarC <- 1
		if(CtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@CMCovar[i]>0){
					r <- r+1
					rowCovarC <- rowCovarC+tKC*(Covar[t+1,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		if(t<(d@orderVC+ordHC)){
			for(i in 1:d@M^(t-d@orderVC)){
				for(j in 1:d@M){
					col <- d@M^(t-d@orderVC)*(j-1)+i
					if(SAlpha[i,t]==0 | SBeta[col,t+1]==0 | d@Pi[AtmCovar*(i-1)+rowCovar,j,t-d@orderVC+1]==0 |d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t+1],j]==0){
						ZEpsilon[AtmCovar*(i-1)+rowCovar,j,t] <- 0
					}
				}
			}
		}else{
			for(i in 1:d@M^ordHC){
				for(j in 1:d@M){
					col <- d@M^(ordHC-1)*(j-1)+floor((i-1)/d@M)+1
					if(SAlpha[i,t]==0 | SBeta[col,t+1]==0 | RA[AtmCovar*(i-1)+rowCovar,j]==0 | d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t+1],j]==0){
						ZEpsilon[AtmCovar*(i-1)+rowCovar,j,t] <- 0
					}
				}
			}
		}
	}
	
	
	#Computation of Epsilon
	
	LSAfact <- 0
	for(t in march.h.seq(d@orderVC+1,s@N-1)){
		pos <- s@y[(t-d@orderVC+1):t]
		rowC <- 1
		for(r in march.h.seq(1,d@orderVC)){
			rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
		}
		
		rowCovar <- 1
		if(AtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@AMCovar[i]>0){
					r <- r+1
					rowCovar <- rowCovar+tKC*(Covar[t+1,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		rowCovarC <- 1
		if(CtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@CMCovar[i]>0){
					r <- r+1
					rowCovarC <- rowCovarC+tKC*(Covar[t+1,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		LSAfact <- LSAfact+SAlog[t]
		if(t<(d@orderVC+ordHC)){
			for(i in 1:d@M^(t-d@orderVC)){
				for(j in 1:d@M){
					if(ZEpsilon[AtmCovar*(i-1)+rowCovar,j,t]!=0){
						Epsilon[AtmCovar*(i-1)+rowCovar,j,t] <- LSAfact+log(SAlpha[i,t])+log(d@Pi[AtmCovar*(i-1)+rowCovar,j,t-d@orderVC+1])+log(d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t+1],j])
					}
				}
			}
		}else{
			for(i in 1:d@M^ordHC){
				for(j in 1:d@M){
					col <- d@M^(ordHC-1)*(j-1)+floor((i-1)/d@M)+1
					if(ZEpsilon[AtmCovar*(i-1)+rowCovar,j,t]!=0){
							Epsilon[AtmCovar*(i-1)+rowCovar,j,t] <- LSAfact+log(SAlpha[i,t])+log(RA[AtmCovar*(i-1)+rowCovar,j])+log(d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t+1],j])
					}
				}
			}
		}
	}
	
	LSBfact <- 0
	
	for(t in march.h.seq(s@N-1,d@orderVC+1,-1)){
		rowCovar <- 1
		if(AtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@AMCovar[i]>0){
					r <- r+1
					rowCovar <- rowCovar+tKC*(Covar[t+1,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		LSBfact <- LSBfact+SBlog[t+1]
		if(t>=(d@orderVC+ordHC)){
			for(i in 1:d@M^ordHC){
				for(j in 1:d@M){
					col <- d@M^(ordHC-1)*(j-1)+floor((i-1)/d@M)+1
					if(ZEpsilon[AtmCovar*(i-1)+rowCovar,j,t]!=0){
						Epsilon[AtmCovar*(i-1)+rowCovar,j,t] <- Epsilon[AtmCovar*(i-1)+rowCovar,j,t]+LSBfact+log(SBeta[col,t+1])-LLAlpha
						Epsilon[AtmCovar*(i-1)+rowCovar,j,t] <- exp(Epsilon[AtmCovar*(i-1)+rowCovar,j,t])
					}
				}
			}
		}else{
			for(i in 1:d@M^(t-d@orderVC)){
				for(j in 1:d@M){
					col <- d@M^(t-d@orderVC)*(j-1)+i
					if(ZEpsilon[AtmCovar*(i-1)+rowCovar,j,t]!=0){
						Epsilon[AtmCovar*(i-1)+rowCovar,j,t] <- Epsilon[AtmCovar*(i-1)+rowCovar,j,t]+LSBfact+log(SBeta[col,t+1])-LLAlpha
						Epsilon[AtmCovar*(i-1)+rowCovar,j,t] <- exp(Epsilon[AtmCovar*(i-1)+rowCovar,j,t])
					}
				}
			}
		}
	}
	
	if(tt==1){
		for(t in march.h.seq(d@orderVC+1,s@N-1)){
			if(d@M>1){
				dist <- colSums(Epsilon[,,t,drop=FALSE])
				Epsilon[,,t] <- kronecker(matrix(1,d@M*AtmCovar,1),matrix(dist,1,d@M))/d@M
			}else{
				Epsilon[,,t] <- rep(sum(Epsilon[,,t]),AtmCovar)
			}
		}
	}
	Epsilon
}


march.dcmm.gam.cov<-function(d,s,Covar,SAlpha,SBeta,SAlog,SBlog,LLAlpha,AtmCovar,CtmCovar){
	
	#Init
	ordHC <- d@orderHC
	if(ordHC==0){ordHC <- 1}
	
	Gamma <- array(0,c(d@M^ordHC*AtmCovar,s@N))
	
	#t between orderVC+1 and orderVC+d@orderHC-1
	
	for(t in march.h.seq(d@orderVC+1,d@orderVC+ordHC-1)){
		rowCovar <- 1
		
		if(AtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@AMCovar[i]>0){
					r <- r+1
					rowCovar <- rowCovar+tKC*(Covar[t,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		for(i in 1:d@M^(t-d@orderVC)){
			if(SAlpha[i,t]*SBeta[i,t]!=0){
				Gamma[AtmCovar*(i-1)+rowCovar,t] <- exp(log(SAlpha[i,t]*SBeta[i,t])+sum(SAlog[1:t])+sum(SBlog[t:s@N])-LLAlpha)
			}
		}
	}
	
	
	#t between orderVC+orderHC and T-1
	
	for(t in march.h.seq(d@orderVC+ordHC,s@N-1)){
		
		rowCovar <- 1
		if(AtmCovar>1){
			r <- 0
			tKC <- 1
			for(i in d@y@Ncov:1){
				if(d@AMCovar[i]>0){
					r <- r+1
					rowCovar <- rowCovar+tKC*(Covar[t,i]-1)
					tKC <- tKC*d@y@Kcov[i]
				}
			}
		}
		
		for(i in 1:d@M^ordHC){
			if(SAlpha[i,t]*SBeta[i,t]!=0){
				Gamma[AtmCovar*(i-1)+rowCovar,t] <- exp(log(SAlpha[i,t]*SBeta[i,t])+sum(SAlog[1:t])+sum(SBlog[t:s@N])-LLAlpha)
			}
		}
	}
	
	
	#Part t=T
	
	rowCovar <- 1
	if(AtmCovar>1){
		r <- 0
		tKC <- 1
		for(i in d@y@Ncov:1){
			if(d@AMCovar[i]>0){
				r <- r+1
				rowCovar <- rowCovar+tKC*(Covar[s@N,i]-1)
				tKC <- tKC*d@y@Kcov[i]
			}
		}
	}
	
	for(i in 1:d@M^ordHC){
		if(SAlpha[i,t]!=0){
			Gamma[AtmCovar*(i-1)+rowCovar,s@N] <- exp(log(SAlpha[i,t])+sum(SAlog[1:s@N])+SBlog[s@N]-LLAlpha)
		}
	}
	Gamma
}

march.dcmm.bw.cov<-function(d,AtmCovar,CtmCovar,Epsilon,Gamma){
	
	tt <- 0
	ordHC <- d@orderHC
	if(ordHC==0){
		ordHC <- 1
		tt <- 1
	}
	
	#Reestimation of the transition matrix at the visible level when no mtd modelling and no covariates are involved
	if(d@Amodel=="complete" & AtmCovar==1){
		
		RA <- array(0,c(d@M^ordHC*AtmCovar,d@M))
		
		for(n in 1:d@y@N){
			for(t in march.h.seq(d@orderVC+ordHC,d@y@T[n]-1)){
				RA <- RA+Epsilon[,,t,n]
			}
		}
		
		tot <- rowSums(RA)
		for(i in 1:d@M^ordHC*AtmCovar){
			if(tot[i]!=0){
				RA[i,] <- RA[i,]/tot[i]
			}
		}
	    d@A <- march.dcmm.cov.h.expandRA(d,RA)
	}
	
	
	#Reestimation of the transition matrix at the visible level when no mtd modelling and no covariates are involved
	if(d@Cmodel=="complete" & CtmCovar==1){
		
		RC <- array(0,c(d@y@K^d@orderVC*CtmCovar,d@y@K,d@M))
		
		for(n in 1:d@y@N){
			for(t in march.h.seq(d@orderVC+1,d@y@T[n])){
				
				pos <- d@y@y[[n]][(t-d@orderVC):(t-1)]
				rowC <- 1
				for(r in march.h.seq(1,d@orderVC)){
					rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
				}
				
				rowCovar <- 1
				if(CtmCovar>1){
					r <- 0
					tKC <- 1
					for(i in d@y@Ncov:1){
						if(d@CMCovar[i]>0){
							r <- r+1
							rowCovar <- rowCovar+tKC*(d@y@cov[n,t,i]-1)
							tKC <- tKC*d@y@Kcov[i]
						}
					}
				}
				
				if(t<(d@orderVC+ordHC)){
					for(m in 1:d@M){
						RC[(rowC-1)*CtmCovar+rowCovar,d@y@y[[n]][t],m] <- RC[(rowC-1)*CtmCovar+rowCovar,d@y@y[[n]][t],m]+sum(Gamma[((m-1)*d@M^(t-d@orderVC)*AtmCovar+1):(m*d@M^(t-d@orderVC)*AtmCovar),t,n])
					}
				}else{
					for(m in 1:d@M){
						RC[(rowC-1)*CtmCovar+rowCovar,d@y@y[[n]][t],m] <- RC[(rowC-1)*CtmCovar+rowCovar,d@y@y[[n]][t],m]+sum(Gamma[((m-1)*d@M^(ordHC-1)*AtmCovar+1):(m*d@M^(ordHC-1)*AtmCovar),t,n])
					}
				}
			}
		}
		
		for(m in 1:d@M){
			if(d@orderVC==0){
				tot <- sum(RC[,,m])
				RC[,,m] <- RC[,,m]/tot
			}else{
				tot <- rowSums(RC[,,m])
				for(j in 1:(d@y@K^d@orderVC*CtmCovar)){
					if(tot[j]!=0){
						RC[j,,m] <- RC[j,,m]/tot[j]
					}
				}
			}
		}
		d@RB <- RC
	}


	#Reestimation of the matrix Pi
	if(tt==1){
		if(AtmCovar==1){
			d@Pi[1,,1] <- d@RA[1,]
		}
	}else{
		d@Pi <- array(0,c(d@M^(ordHC-1)*AtmCovar,d@M,ordHC))
		for(n in 1:d@y@N){
			for(m in 1:d@M){
				d@Pi[1:AtmCovar,m,1] <- d@Pi[1:AtmCovar,m,1]+Gamma[((m-1)*AtmCovar+1):(m*AtmCovar),d@orderVC+1,n]
			}
		}
		if(ordHC>1){
			for(t in 2:ordHC){
				for(n in 1:d@y@N){
					for(i in 1:d@M){
						d@Pi[1:(d@M^(t-1)*AtmCovar),i,t] <- d@Pi[1:(d@M^(t-1)*AtmCovar),i,t]+Gamma[((i-1)*d@M^(t-1)*AtmCovar+1):(i*d@M^(t-1)*AtmCovar),t+d@orderVC,n]
					}
				}
			}
		}
		
		for(t in 1:ordHC){
			for(i in 1:(d@M^(t-1)*AtmCovar)){
				tot <- sum(d@Pi[i,,t])
				if(tot>0){
					d@Pi[i,,t] <- d@Pi[i,,t]/tot
				}
			}
		}
	}
	d
}


GMTD.dp<-function(order,NSS,ProbT,Phi){
	SDip <- rep(0,order)
	
	for(k in 1:length(NSS)){
		for(m in 1:order){
			if(NSS[k]>0 & Phi%*%ProbT[k,]!=0){
				SDip[m] <- SDip[m]+NSS[k]*ProbT[k,m]/(Phi%*%ProbT[k,])
			}
		}
	}
	SDip
}

GMTD.op.dcmm.A<-function(d,NSS,SDip,Delta,Stop,Constraint,AtmCovar,CtmCovar){
	
	SDelta <- Delta
	LPhi <- length(d@APhi)
	
	SDip_p <- order(SDip)
	SDip_s <- SDip[SDip_p]
	d_inc <- SDip_p[LPhi]
	d_dec <- SDip_p[1]
		
	if(Constraint==1 & d@APhi[SDip_p[LPhi]]==1){
		MD <- Delta
		return(list(d=d,MD=MD))
	}
	
	if(Constraint==1 & d@APhi[SDip_p[1]]==0){
		for(k in 2:LPhi){
			if(d@APhi[SDip_p[k]]!=0){
				d_dec <- SDip_p[k]
				break
			}
		}
	}
	
	repeat{
		if(Constraint==1){
			maxchange <- min(1-d@APhi[d_inc],d@APhi[d_dec],Delta)
		}else{
			maxchange <- Delta
		}
		dtmp <- d
		tPhi <- d@APhi
		tPhi[d_inc] <- tPhi[d_inc]+maxchange
		tPhi[d_dec] <- tPhi[d_dec]-maxchange
		dtmp@APhi <- tPhi
		
		l <- GMTD_tm_cov(max(dtmp@orderHC,1),dtmp@M,dtmp@APhi,dtmp@AProbT)
		dtmp@A <- l$HOQ
		nbnpCRHOQ <-l$nbnpCRHOQ
		
		#Check for non-degenerescence
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@APhi%*%dtmp@AProbT[k,]
				if(calc==0){
					dtmp@ll <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			SAlpha <- array(0,c(d@M^max(d@orderHC,1),max(d@y@T),d@y@N))
			SAlog <- array(0,c(max(d@y@T),d@y@N))
			
			tLL <- 0
			for(n in 1:d@y@N){
				s <- march.dataset.h.extractSequence(d@y,n)
				tLL <- tLL+march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)$LLAlpha
			}
			dtmp@ll <- tLL
		}
		
		#Test of improvement
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				MD <- 2*Delta
				return(list(d=d,MD=MD))
			}
		}
	}
	list(d=d,MD=MD)
}


march.mccov.dq<-function(order,K,NSS,ProbT,ValT,Phi){
	SDiq <- array(0,c(K^order,K))
	
	for(k in 1:length(NSS)){
		if(NSS[k]>0 & Phi%*%ProbT[k,]!=0){
			
			#Computation of the row in SDiq
			if(order==0){
				row <- 1
			}else{
				row <- ValT[k,order+1]
				for(g in march.h.seq(order,2,-1)){
					row <- row+K^(order-g+1)*(ValT[k,g]-1)
				}
			}
			SDiq[row,ValT[k,1]] <- SDiq[row,ValT[k,1]]+NSS[k]*Phi[1]/(Phi%*%ProbT[k,])
		}
	}
	SDiq
}


GMTD.oq.dcmm.A.b<-function(d,AtmCovar,CtmCovar,ValT,NSS,Qr,SDiq,Delta,Stop){
	
	#Init
	tt <- 0
	ordHC <-d@orderHC
	if(d@orderHC==0){
		ordHC <- 1
		tt <-1
	}
	SDelta <- Delta
	
	#Sort of SDiq
	SDiq <- SDiq[Qr,]
	SDiq_p <- order(SDiq)
	SDiq_s <- SDiq[SDiq_p]
	d_inc <- SDiq_p[d@M]
	d_dec <- SDiq_p[1]
	
	if(d@AQ[1,Qr,SDiq_p[d@M]]==1){
		MD <- Delta
		return(list(d=d,MD=MD))
	}
	
	if(d@AQ[1,Qr,SDiq_p[1]]==0){
		for(k in 2:d@M){
			if(d@AQ[1,Qr,SDiq_p[k]]!=0){
				d_dec <- SDiq_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@AQ[1,Qr,d_inc],d@AQ[1,Qr,d_dec])
		tQ <- d@AQ
		dtmp <- d
		tQ[1,Qr,d_inc] <- tQ[1,Qr,d_inc]+min(Delta,maxchange)
		tQ[1,Qr,d_dec] <- tQ[1,Qr,d_dec]-min(Delta,maxchange)
		
		#Case orderHC=0
		if(tt==1){
			tQ <- kronecker(matrix(1,d@M,1),tQ[1,,drop=FALSE])
		}
		
		dtmp@AQ <- tQ
		
		#New matrix ProbT
		dtmp@AProbT <- march.mccov.sp(ordHC,dtmp@AQ[1,,],ValT,dtmp@AProbT)
		
		#Computation of A 
		l <- GMTD_tm_cov(ordHC,dtmp@M,dtmp@APhi,dtmp@AProbT)
		dtmp@A <- l$HOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		#Check for non-degenerescence
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@APhi%*%dtmp@AProbT[k,]
				if(calc==0){
					dtmp@ll <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			tLL <- 0
			for(n in 1:d@y@N){
				s <- march.dataset.h.extractSequence(d@y,n)
				tLL <- tLL+march.dcmm.fp.cov(d,s,d@y@cov[n,1:s@N,],AtmCovar,CtmCovar)$LLAlpha
			}
			dtmp@ll <- tLL
		}

		#Test of improvement
		test <- 0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		if(test==1){break}
	}
	list(d=d,MD=MD)
}


GMTD.dq<-function(order,K,NSS,ProbT,ValT,Phi){
	SDiq <- array(0,c(K,K))
	
	for(k in 1:length(NSS)){
		for(m in 1:order){
			if(NSS[k]>0 & Phi%*%ProbT[k,]!=0){
				SDiq[ValT[k,m+1],ValT[k,1]] <- SDiq[ValT[k,m+1],ValT[k,1]]+NSS[k]*Phi[m]/(Phi%*%ProbT[k,])
			}
		}
	}
	SDiq
}


GMTD.oq.dcmm.A<-function(d,AtmCovar,CtmCovar,ValT,NSS,Qr,SDiq,Delta,Stop){
	SDelta <- Delta
	
	#Sort of SDiq
	SDiq <- SDiq[Qr,]
	SDiq_p <- order(SDiq)
	SDiq_s <- SDiq[SDiq_p]
	
	d_inc <- SDiq_p[d@M]
	d_dec <- SDiq_p[1]
	
	if (d@AQ[1,Qr,SDiq_p[d@M]]==1){
		MD <- Delta
		return(list(d=d,MD=MD))
	}
	
	if(d@AQ[1,Qr,SDiq_p[1]]==0){
		for(k in 2:d@M){
			if(d@AQ[1,Qr,SDiq_p[k]]!=0){
				d_dec <- SDiq_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@AQ[1,Qr,d_inc],d@AQ[1,Qr,d_dec])
		
		dtmp <- d
		dtmp@AQ[1,Qr,d_inc] <- dtmp@AQ[1,Qr,d_inc]+min(Delta,maxchange)
		dtmp@AQ[1,Qr,d_dec] <- dtmp@AQ[1,Qr,d_dec]-min(Delta,maxchange)
		
		dtmp@AProbT <- GMTD.sp(dtmp@orderHC,dtmp@AQ[1,,],ValT,NSS,dtmp@AProbT)
		
		l <- GMTD_tm_cov(dtmp@orderHC,dtmp@M,dtmp@APhi,dtmp@AProbT)
		dtmp@A <- l$HOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@APhi%*%dtmp@AProbT[k,]
				if(calc==0){
					tLL <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			tLL <- 0
			
			for(n in 1:dtmp@y@N){
				s <- march.dataset.h.extractSequence(dtmp@y,n)
				tLL <- tLL+march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)$LLAlpha
			}
			dtmp@ll <- tLL
		}
		
		test <- 0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		
		if(test==1){break}
	}
	return(list(d=d,MD=MD))	
}


GMTD.sp<-function(order,Q,ValT,NSS,ProbT){
	for(i in 1:length(NSS)){
		for(j in 1:order){
			ProbT[i,j] <- Q[ValT[i,j+1],ValT[i,1]]
		}
	}
	ProbT
}


GMTDg.dq <- function(order,K,NSS,ProbT,ValT,Phi){
	SDiq <- array(0,c(K,K,order))
	
	for(k in 1:length(NSS)){
		for(m in 1:order){
			if(NSS[k]>0 & Phi%*%ProbT[k,]){
				SDiq[ValT[k,m+1],ValT[k,1],m] <- SDiq[ValT[k,m+1],ValT[k,1],m]+NSS[k]*Phi[m]/(Phi%*%ProbT[k,])
			}
		}
	}
	SDiq
}


GMTDg.oq.dcmm.A<-function(d,AtmCovar,CtmCovar,ValT,NSS,Qr,Ql,SDiq,Delta,Stop){
	
	#Init
	SDelta <- Delta
	
	#Sort of SDiq
	SDiq <- SDiq[Qr,,Ql]
	SDiq_p <- order(SDiq)
	SDiq_s <- SDiq[SDiq_p]
	d_inc <- SDiq_p[d@M]
	d_dec <- SDiq_p[1]
	
	if(d@AQ[Ql,Qr,SDiq_p[d@M]]==1){
		MD <- Delta
		return(list(d=d,MD=MD))
	}
	
	if(d@AQ[Ql,Qr,SDiq_p[1]]==0){
		for(k in 2:d@M){
			if(d@AQ[Ql,Qr,SDiq_p[k]]!=0){
				d_dec <- SDiq_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@AQ[Ql,Qr,d_inc],d@AQ[Ql,Qr,d_dec])
		
		dtmp <- d
		dtmp@AQ[Ql,Qr,d_inc] <- dtmp@AQ[Ql,Qr,d_inc]+min(Delta,maxchange)
		dtmp@AQ[Ql,Qr,d_dec] <- dtmp@AQ[Ql,Qr,d_dec]-min(Delta,maxchange)
		
		dtmp@AProbT <- GMTDg.sp(dtmp@orderHC,dtmp@M,dtmp@AQ,ValT,NSS,dtmp@AProbT)
		
		l <- GMTD_tm_cov(dtmp@orderHC,dtmp@M,dtmp@APhi,dtmp@AProbT)
		dtmp@A <- l$HOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@APhi%*%dtmp@AProbT[k,]
				if(calc==0){
					tLL <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			
			tLL <- 0
			
			for(n in 1:dtmp@y@N){
				s <- march.dataset.h.extractSequence(dtmp@y,n)
				tLL <- tLL+march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)$LLAlpha
			}
			dtmp@ll <- tLL
		}
		
		test <- 0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		
		if(test==1){break}
		
	}
	return(list(d=d,MD=MD))
}


GMTD.dcov<-function(CCov,ColVT,K,KCov,NSS,ProbT,ValT,Phi){
	
	SDicov <- matrix(0,KCov,K)
	
	for(k in 1:length(NSS)){
		if(NSS[k]>0 & Phi%*%ProbT[k,]!=0){
			SDicov[ValT[k,ColVT],ValT[k,1]] <- SDicov[ValT[k,ColVT],ValT[k,1]]+NSS[k]*Phi[CCov]/(Phi%*%ProbT[k,])
		}
	}
	SDicov
}

GMTD.ocov.dcmm.A<-function(d,SDicov,AtmCovar,CtmCovar,ValT,NSS,Tl,Tr,Pcol,Delta,Stop){
	
	#Init
	tt <- 0
	ordHC <- d@orderHC
	if(d@orderHC==0){
		ordHC <- 1
		tt <-1
	}
	
	SDelta <- Delta
	
	#Sort of SDcov
	SDicov <- SDicov[Tr,]
	SDicov_p <- order(SDicov)
	SDicov_s <- SDicov[SDicov_p]
	
	d_inc <- SDicov_p[d@M]
	d_dec <- SDicov_p[1]
	
	if(d@ATCovar[[Tl]][Tr,SDicov_p[d@M]]==1){
		MD <- Delta
		return(list(d=d,MD=MD))
	}
	
	if(d@ATCovar[[Tl]][Tr,SDicov_p[1]]==0){
		for(k in 2:d@M){
			if(d@ATCovar[[Tl]][Tr,SDicov_p[k]]!=0){
				d_dec <- SDicov_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@ATCovar[[Tl]][Tr,d_inc],d@ATCovar[[Tl]][Tr,d_dec])
		
		dtmp <- d
		dtmp@ATCovar[[Tl]][Tr,d_inc] <- dtmp@ATCovar[[Tl]][Tr,d_inc]+min(Delta,maxchange)
		dtmp@ATCovar[[Tl]][Tr,d_dec] <- dtmp@ATCovar[[Tl]][Tr,d_dec]-min(Delta,maxchange)
		
		if(d@Amodel=="complete"){
			dtmp@AProbT <- GMTD.mc.spcov(dtmp@ATCovar[[Tl]],ValT,NSS,dtmp@AProbT,Pcol,ordHC)
		}else{
			dtmp@AProbT <- GMTD.spcov(dtmp@ATCovar[[Tl]],ValT,NSS,dtmp@AProbT,Pcol)
		}
		
		l <- GMTD_tm_cov(dtmp@orderHC,dtmp@M,dtmp@APhi,dtmp@AProbT)
		dtmp@A <- l$HOQ
		RA <- l$CRHOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		if(tt==1){
			dtmp@Pi <- RA[1:AtmCovar,]
		}else{
			dtmp@Pi <- d@Pi
		}
		
		tLLtest <- 0
		
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@APhi%*%dtmp@AProbT[k,]
				if(calc==0){
					tLL=-Inf
					tLLtest=1
					break
				}
			}
		}
		if(tLLtest==0){
						
			tLL <- 0
			for(n in 1:d@y@N){
				s <- march.dataset.h.extractSequence(d@y,n)
				tLL <- tLL+march.dcmm.fp.cov(d,s,d@y@cov[n,1:s@N,],AtmCovar,CtmCovar)$LLAlpha
			}
			dtmp@ll <- tLL
		}

		#Test of improvement
		
		test=0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d<-dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		if(test==1){
			break
		}
	}
	list(d=d,MD=MD)

}


GMTD.spcov<-function(TCov,ValT,NSS,ProbT,Pcol){
	for(i in 1:length(NSS)){
		ProbT[i,Pcol] <- TCov[ValT[i,Pcol+1],ValT[i,1]]
	}
	ProbT
}


GMTD.mc.spcov<-function(TCov,ValT,NSS,ProbT,Pcol,order){
	for(i in 1:length(NSS)){
		ProbT[i,Pcol] <- TCov[ValT[i,Pcol+order],ValT[i,1]]
	}
	ProbT
}


GMTD.op.dcmm.C<-function(d,AtmCovar,CtmCovar,NSS,SDip,state,Delta,Stop,Constraint){
	SDelta <- Delta
	LPhi <- dim(d@CPhi)[2]
	
	SDip_p <- order(SDip)
	SDip_s <- SDip[SDip_p]
	d_inc <- SDip_p[LPhi]
	d_dec <- SDip_p[1]
	
	if(Constraint==1 & d@CPhi[1,SDip_p[LPhi],state]==1){
		MD <- Delta
		return(list(d=d,MD=MD))
	}
	
	if(Constraint==1 & d@CPhi[1,SDip_p[1],state]==0){
		for(k in 2:LPhi){
			if(d@CPhi[1,SDip_p[k],state]!=0){
				d_dec <- SDip_p[k]
				break
			}
		}
	}
	
	repeat{
		
		if(Constraint==1){
			maxchange <- min(1-d@CPhi[1,d_inc,state],d@CPhi[1,d_dec,state], Delta)
		}else{
			maxchange <- Delta
		}
		
		dtmp <- d
		
		dtmp@CPhi[1,d_inc,state] <- dtmp@CPhi[1,d_inc,state]+maxchange
		dtmp@CPhi[1,d_dec,state] <- dtmp@CPhi[1,d_dec,state]-maxchange
		
		l <- GMTD_tm_cov(dtmp@orderVC,dtmp@y@K,matrix(dtmp@CPhi[,,state],1,dim(dtmp@CPhi)[2]),dtmp@CProbT[,,state])
		dtmp@RB[,,state] <- l$CRHOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@CPhi[,,state]%*%dtmp@CProbT[k,,state]
				if(calc==0){
					dtmp@ll <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			if(d@M==1){
				dtmp@ll <- as.numeric(GMTD.ll(NSS,dtmp@CProbT[,,state],dtmp@CPhi[,,state]))
			}else{
				tLL <- 0
				for(n in 1:dtmp@y@N){
					s <- march.dataset.h.extractSequence(dtmp@y,n)
					tLL <- tLL+march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)$LLAlpha
				}
			  dtmp@ll <- tLL
			}
		}
		
		test <- 0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		
		if(test==1){break}
	}
	return(list(d=d,MD=MD))
}


GMTD.oq.dcmm.C <- function(d,AtmCovar,CtmCovar,ValT,NSS,Qr,state,SDiq,Delta,Stop){
	SDelta <- Delta
	
	SDiq <- SDiq[Qr,]
	SDiq_p <- order(SDiq)
	SDiq_s <- SDiq[SDiq_p]
	d_inc <- SDiq_p[d@y@K]
	d_dec <- SDiq_p[1]
	
	if(d@CQ[1,Qr,SDiq_p[d@y@K],state]==1){
		return(list(d=d,MD=Delta))
	}
	
	if(d@CQ[1,Qr,SDiq_p[1],state]==0){
		for(k in 2:d@y@K){
			if(d@CQ[1,Qr,SDiq_p[k],state]!=0){
				d_dec <- SDiq_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@CQ[1,Qr,d_inc,state],d@CQ[1,Qr,d_dec,state])
		
		dtmp <- d
		
		dtmp@CQ[1,Qr,d_inc,state] <- dtmp@CQ[1,Qr,d_inc,state]+min(Delta,maxchange)
		dtmp@CQ[1,Qr,d_dec,state] <- dtmp@CQ[1,Qr,d_dec,state]-min(Delta,maxchange)
		
		dtmp@CProbT[,,state] <- GMTD.sp(dtmp@orderVC,dtmp@CQ[1,,,state],ValT,NSS,matrix(dtmp@CProbT[,,state],dim(dtmp@CProbT)[1],dim(dtmp@CProbT)[2]))
		
		l <- GMTD_tm_cov(dtmp@orderVC,dtmp@y@K,matrix(dtmp@CPhi[,,state],1,dim(dtmp@CPhi)[2]),matrix(dtmp@CProbT[,,state],dim(dtmp@CProbT)[1],dim(dtmp@CProbT)[2]))
		dtmp@RB[,,state] <- l$CRHOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@CPhi[,,state]%*%dtmp@CProbT[k,,state]
				if(calc==0){
					dtmp@ll <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			tLL <- 0
			for(n in 1:dtmp@y@N){
				s <- march.dataset.h.extractSequence(dtmp@y,n)
				tLL <- tLL+march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)$LLAlpha
			}
		dtmp@ll <- tLL
		}
		
		test <- 0
		
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		
		if(test==1){break}
		
	}
	return(list(d=d,MD=MD))
}


GMTD.oq.dcmm.C.b<-function(d,AtmCovar,CtmCovar,ValT,NSS,Qr,state,SDiq,Delta,Stop){
	SDelta <- Delta
	
	SDiq <- SDiq[Qr,]
	SDiq_p <- order(SDiq)
	SDiq_s <- SDiq[SDiq_p]
	d_inc <- SDiq_p[d@y@K]
	d_dec <- SDiq_p[1]
	
	if(d@CQ[1,Qr,SDiq_p[d@y@K],state]==1){
		MD <- Delta
		return(list(d=d,MD=MD))
	}
	
	if(d@CQ[1,Qr,SDiq_p[1],state]==0){
		for(k in 2:d@y@K){
			if(d@CQ[1,Qr,SDiq_p[k],state]!=0){
				d_dec <- SDiq_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@CQ[1,Qr,d_inc,state],d@CQ[1,Qr,d_dec,state])
		dtmp <- d
		dtmp@CQ[1,Qr,d_inc,state] <- dtmp@CQ[1,Qr,d_inc,state]+min(Delta,maxchange)
		dtmp@CQ[1,Qr,d_dec,state] <- dtmp@CQ[1,Qr,d_dec,state]-min(Delta,maxchange)
		
		dtmp@CProbT[,,state] <- march.mccov.sp(dtmp@orderVC,dtmp@CQ[1,,,state],ValT,dtmp@CProbT[,,state])
		
		l <- GMTD_tm_cov(dtmp@orderVC,dtmp@y@K,dtmp@CPhi[,,state],dtmp@CProbT[,,state])
		dtmp@RB[,,state] <- l$CRHOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@CPhi[,,state]%*%dtmp@CProbT[k,,state]
				if(calc==0){
					dtmp@ll <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			if(d@M==1){
				dtmp@ll <- GMTD.ll(NSS,dtmp@CProbT[k,,state],dtmp@CPhi[,,state])
			}else{
				tLL <- 0
				for(n in 1:dtmp@y@N){
					s <- march.dataset.h.extractSequence(dtmp@y,n)
					tLL <- tLL+march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)
				}
			dtmp@ll <- tLL
			}
		}
		test <- 0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		
		if(test==1){break}
		
	}
	return(list(d=d,MD=MD))
}


GMTDg.oq.dcmm.C<-function(d,AtmCovar,CtmCovar,ValT,NSS,Ql,Qr,state,SDiq,Delta,Stop){
	SDelta <- Delta
	
	SDiq <- SDiq[Qr,,Ql]
	SDiq_p <- order(SDiq)
	SDiq_s <- SDiq[SDiq_p]
	d_inc <- SDiq_p[d@y@K]
	d_dec <- SDiq_p[1]
	
	if(d@CQ[Ql,Qr,SDiq_p[d@y@K],state]==1){
		return(list(d=d,MD=Delta))
	}
	
	if(d@CQ[Ql,Qr,SDiq_p[1],state]==0){
		for(k in 2:d@y@K){
			if(d@CQ[Ql,Qr,SDiq_p[k],state]!=0){
				d_dec <- SDiq_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@CQ[Ql,Qr,d_inc,state],d@CQ[Ql,Qr,d_dec,state])
		
		dtmp <- d
		
		dtmp@CQ[Ql,Qr,d_inc,state] <- dtmp@CQ[Ql,Qr,d_inc,state]+min(Delta,maxchange)
		dtmp@CQ[Ql,Qr,d_dec,state] <- dtmp@CQ[Ql,Qr,d_dec,state]-min(Delta,maxchange)
		
		dtmp@CProbT[,,state] <- GMTDg.sp(dtmp@orderVC,dtmp@y@K,dtmp@CQ[,,,state],ValT,NSS,dtmp@CProbT[,,state])
		
		l <- GMTD_tm_cov(dtmp@orderVC,dtmp@y@K,matrix(dtmp@CPhi[,,state],1,dim(dtmp@CPhi)[2]),dtmp@CProbT[,,state])
		
		dtmp@RB[,,state] <- l$CRHOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		tLLtest <- 0
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@CPhi[,,state]%*%dtmp@CProbT[k,,state]
				if(calc==0){
					dtmp@ll <- -Inf
					tLLtest <- 1
					break
				}
			}
		}
		
		if(tLLtest==0){
			if(d@M==1){
				dtmp@ll <- as.numeric(GMTD.ll(NSS,dtmp@CProbT[,,state],dtmp@CPhi[,,state]))
			}else{
				tLL<-0
				for(n in 1:dtmp@y@N){
					s <- march.dataset.h.extractSequence(dtmp@y,n)
					tLL <- tLL+march.dcmm.fp.cov(dtmp,s,dtmp@y@cov[n,,],AtmCovar,CtmCovar)$LLAlpha
				}
			dtmp@ll <- tLL
			}
		}
		
		test <- 0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		
		if(test==1){break}
		
	}
	return(list(d=d,MD=Delta))

}


GMTDg.sp<-function(order,K,Q,ValT,NSS,ProbT){
	for(r in 1:length(NSS)){
		for(c in 1:order){
			ProbT[r,c] <- Q[c,ValT[r,c+1],ValT[r,1]]
		}
	}
	ProbT
}

GMTD.ll<-function(NSS,ProbT,Phi){
	LL <- 0
	
	for(k in 1:length(NSS)){
		if(NSS[k]>0){
			calc <- Phi%*%ProbT[k,]
			if(calc==0){
				LL <- -Inf
				break
			}
			LL <- LL+NSS[k]*log(calc)
		}
	}
	LL
}

GMTD.ocov.dcmm.C<-function(d,SDicov,AtmCovar,CtmCovar,ValT,NSS,Tl,Tr,state,Pcol,Delta,Stop){
	SDelta <- Delta
	
	SDicov <- SDicov[Tr,]
	SDicov_p <- order(SDicov)
	SDicov_s <- SDicov[SDicov_p]
	d_inc <- SDicov_p[d@y@K]
	d_dec <- SDicov_p[1]
	
	if(d@CTCovar[[Tl]][Tr,SDicov_p[d@y@K],state]==1){
		return(list(d=d,MD=Delta))
	}
	
	if(d@CTCovar[[Tl]][Tr,SDicov_p[1],state]==0){
		for(k in 2:d@y@K){
			if(d@CTCovar[[Tl]][Tr,SDicov_p[k],state]!=0){
				d_dec <- SDicov_p[k]
				break
			}
		}
	}
	
	repeat{
		maxchange <- min(1-d@CTCovar[[Tl]][Tr,d_inc,state],d@CTCovar[[Tl]][Tr,d_dec,state])
		
		dtmp <- d
		
		dtmp@CTCovar[[Tl]][Tr,d_inc,state] <- dtmp@CTCovar[[Tl]][Tr,d_inc,state]+min(Delta,maxchange)
		dtmp@CTCovar[[Tl]][Tr,d_dec,state] <- dtmp@CTCovar[[Tl]][Tr,d_dec,state]-min(Delta,maxchange)
		
		if(dtmp@Cmodel=="complete"){
			dtmp@CProbT[,,state] <- GMTD.mc.spcov(dtmp@CTCovar[[Tl]][,,state],ValT,NSS,dtmp@CProbT[,,state],Pcol,dtmp@orderVC)
		}else{
			dtmp@CProbT[,,state] <- GMTD.spcov(dtmp@CTCovar[[Tl]][,,state],ValT,NSS,dtmp@CProbT[,,state],Pcol)
		}
		
		l <- GMTD_tm_cov(dtmp@orderVC,dtmp@y@K,matrix(dtmp@CPhi[,,state],1,dim(dtmp@CPhi)[2]),dtmp@CProbT[,,state])
		dtmp@RB[,,state] <- l$CRHOQ
		nbnpCRHOQ <- l$nbnpCRHOQ
		
		tLLtest <- 0
		
		for(k in 1:length(NSS)){
			if(NSS[k]>0){
				calc <- dtmp@CPhi[,,state]%*%dtmp@CProbT[k,,state]
				if(calc==0){
					tLL=-Inf
					tLLtest=1
					break
				}
			}
		}
		if(tLLtest==0){
			if(d@M==1){
				dtmp@ll <- as.numeric(GMTD.ll(NSS,dtmp@CProbT[,,state],dtmp@CPhi[,,state]))
			}else{
				
				tLL <- 0
				for(n in 1:d@y@N){
					s <- march.dataset.h.extractSequence(d@y,n)
					tLL <- tLL+march.dcmm.fp.cov(d,s,d@y@cov[n,1:s@N,],AtmCovar,CtmCovar)$LLAlpha
				}
				dtmp@ll <- tLL
			}
		}

		#Test of improvement
		
		test <- 0
		if(nbnpCRHOQ==0 & dtmp@ll>d@ll){
			d <- dtmp
			if(SDelta==Delta){
				MD <- 2*Delta
			}else{
				MD <- Delta
			}
			break
		}else{
			Delta <- Delta/2
			if(Delta<=Stop){
				test <- 1
				MD <- 2*Delta
				break
			}
		}
		
		if(test==1){break}
		
	}
	return(list(d=d,MD=MD))
}

# Construction of the approximated high-order transition matrix for a MTD model with (if any) covariates
GMTD_tm_cov<-function(Order,K,Phi,ProbT){
  nr=dim(ProbT)[1]
  nc=dim(ProbT)[2]
  CRHOQ=matrix(0,nr/K,K)
  
  #Computation of CRHOQ
  row=0
  col=1
  for(i in 1:nr){
    row=row+1
    CRHOQ[row,col]=ProbT[i,]%*%t(Phi)
    if (row==nr/K){
      row=0
      col=col+1
    }
  }
  for (i in 1:nr/K){
    if(sum(CRHOQ[i,])<0.999){
      CRHOQ[i,]=matrix(0,1,K)
    }
  }
  
  #Computation of HOQ
  if(Order==0){
    HOQ=CRHOQ
  }else{
    HOQ=matrix(0,nr/K,K^Order)
    tt=(nr/K)/K^Order
    offset=0
    for(r1 in 1:K^(Order-1)){
      for(r2 in 1:K){
        for(off in 1:tt){
          offset=offset+1
          for(c in 1:K){
            HOQ[offset,r1+(c-1)*K^(Order-1)]=CRHOQ[offset,c]
          }
        }
      }
    }
  }
  nbnpCRHOQ=length(which(CRHOQ<0|CRHOQ>1))
  return(list(HOQ=HOQ,CRHOQ=CRHOQ,nbnpCRHOQ=nbnpCRHOQ))
}

BuildArrayNumberOfSequencesDCMM<-function(y,order,MCovar,NbMCovar,placeCovar){
	tmCovar <- 1
  	if(NbMCovar>0){
  		for(i in 1:NbMCovar){
      		tmCovar <- tmCovar*y@Kcov[placeCovar[i]]
      	}
  	}
  	
  	#Init
  	NSS <- rep(0,y@K^(order+1)*tmCovar)
  	
  	for(n in 1:y@N){
    	for(t in march.h.seq(1,y@T[n]-order)){
      		pos <- y@y[[n]][t:(t+order)]
      		
      		row1 <- pos[1]
      		for(g in march.h.seq(2,order+1)){
        		row1 <- row1+y@K^(g-1)*(pos[g]-1)
      		}
      		
      		if(NbMCovar>0){
        		row2 <- y@cov[n,t+order,placeCovar[NbMCovar]]
        		totKC <- 1
        		if(NbMCovar>1){
          			for(i in (NbMCovar-1):1 ){
            			totKC <- totKC*y@Kcov[placeCovar[i+1]]
            			row2 <- row2+totKC*(y@cov[n,t+order,placeCovar[i]]-1)
          			}
        		}
      		}
      		
      		if(NbMCovar>0){
        		row1 <- (row1-1)*tmCovar+row2
      		}
      		
      		NSS[row1] <- NSS[row1]+1
    	}
  	}
  	return(NSS)
}

BuildContingencyTableDCMM <- function(y,order,NbCovar,placeCovar){
  n_rows_data <- y@N # number of rows (number of data sequences)
  l<-list() # crosstable (rt, Cg in page 385 of Berchtold, 2001)
  if(order>0){
  for (g in 1:order){
    c=array(0,c(y@K,y@K))
    for (i in 1:n_rows_data){
      for (t in 1:(y@T[i]-g)){ # mc_lag is g in Berchtold, 2001
        past <- y@y[[i]][t]
        present <- y@y[[i]][t+g]
        if(length((which(c(past,present)<1) | (which(c(past,present)>y@K))))==0){
          c[past,present] <- c[past,present] + y@weights[i]
        }
      }
    }
    l[[g]]<-c
  }
  }
  
  if(NbCovar>0){
    for(j in 1:NbCovar){
      kcov=y@Kcov[placeCovar[j]]
      CT=matrix(0,kcov,y@K)
      
      for(n in 1:y@N){
        for (i in 1:y@T[n]){
          row=y@cov[n,i,placeCovar[j]]
          col=y@y[[n]][i]
          CT[row,col]=CT[row,col]+y@weights[n]
        }
      }
      l[[order+j]]<-CT
    }
  }
  return(l)
}

CalculateTheilUDCMM <- function(y,order,c,NbCovar,placeCovar){
  
  u <- array(data=NA,dim=c(1,order+NbCovar))
  for (g in 1:order){
    cg <- c[[g]]
    tc <- sum(cg) # sum of elements (TCg)
    sr <- rowSums(cg) # vector of sums of rows [Cg(.,j)]
    sc <- colSums(cg) # vector of sums of columns [Cg(i,.)]
    
    # the following lines implement equation 14
    num <- 0
    for (i in 1:y@K){
      for (j in 1:y@K){
        if (cg[i,j]!=0){ # if cg[i,j], there is nothing to be added
          num <- num + cg[i,j] * (log2(sc[i]) + log2(sr[j]) - log2(cg[i,j]) - log2(tc))
        }
      }
    }
    
    den <- 0
    
    for (j in 1:y@K){
      if(sc[j]!=0 & tc!=0){
        den <- den + sc[j] * (log2(sc[j]) - log2(tc))
      }
    }
    
    if(den != 0) u[g] = num/den
    else u[g] = NaN
    
  }
  if(NbCovar>0){
    for (g in 1:NbCovar){
      cg <- c[[order+g]]
      tc <- sum(cg) # sum of elements (TCg)
      sr <- rowSums(cg) # vector of sums of rows [Cg(.,j)]
      sc <- colSums(cg) # vector of sums of columns [Cg(i,.)]
      
      # the following lines implement equation 14
      num <- 0
      for (i in 1:y@Kcov[placeCovar[g]]){
        for (j in 1:y@K){
          if (cg[i,j]!=0){ # if cg[i,j], there is nothing to be added
            num <- num + cg[i,j] * (log2(sc[j]) + log2(sr[i]) - log2(cg[i,j]) - log2(tc))
          }
        }
      }
      
      den <- 0
      
      for (j in 1:y@K){
        if(sc[j]!=0 & tc!=0){
          den <- den + sc[j] * (log2(sc[j]) - log2(tc))
        }
      }
      
      if(den != 0) u[order+g] = num/den else u[order+g] = NaN
      
    }
  }
  return(u)
}
