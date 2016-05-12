# TODO: Add comment
# 
# Author: SSP24298
###############################################################################


# The implementation of the main loop of an EA.
# 
# The function executes an evolutionary loop, iterating over classical genetic operations, like crossover, mutation
# evaluation and selection. This algorithm uses a Lamarckian EA, that optimize the children population if 
# p@@optimizing has been set to TRUE, using the function indicated in the slot of the optimizing parameter 
# (\code{\link{march.ea.OptimizingParameters}}).
# 
# @param p a parameter object of type \code{\link{march.ea.Parameters}} describing the option of the EA.
# @return 
#  \item{$children}{Last children population.}
#  \item{$best}{The best individual found during the run.}
#  
# This loop is designed to be generic enough, so every functions that interact really with 
# the genome are written into march.dcmm.ea.R file.
#
march.loop <- function(p){
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
				children <- p@crossoverParameters@fct(population[[sel1]],population[[sel2]])
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

# The implementation of an sbx crossover operator, that emulate bit-flip crossover operation,
# for real numbers. Here the operator is specialized for muting probability [0,1]
march.ea.h.sbx <- function(v1,v2,contiguity){
	c1 <- array(0,length(v1))
	c2 <- array(0,length(v1))
	
	ub <- 1
	lb <- 0
	
	for( i in 1:length(c1) ){		
		if( runif(1)<0.5 ){
			u <- runif(1)
			if( u<0.5 ){ # contracting crossover
				beta <- (2*u)^(1/(contiguity+1))
			}
			else{ # expanding crossover
				beta <- (0.5/(1-u))^(1.0/(contiguity+1))
			}
			
			c1[i] <- ((v1[i]+v2[i])/2)-beta*0.5*abs(v1[i]-v2[i])
			c2[i] <- ((v1[i]+v2[i])/2)+beta*0.5*abs(v1[i]-v2[i])
			
			if( c1[i] > ub )	{ c1[i] = ub }
			else if( c1[i]<lb )	{ c1[i] = lb }
			
			if( c2[i] > ub )	{ c2[i] = ub }
			else if( c2[i]<lb )	{ c2[i] = lb }
			
		}
		else{
			c1[i] <- v1[i]
			c2[i] <- v2[i]
		}
	}
	
	list(c1=c1,c2=c2)
}

march.ea.checkBounds <- function(v,lb,ub){
	slb <- v<lb	
	sub <- v>ub
	
	v <- lb*slb+v*!slb
	v <- ub*sub+v*!sub
	
	v
}