# Copyright (C) 2015  Marco Manfrini, PhD

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

###FUNCTIONS

# GENERAL PURPOSE FUNCTIONS
# SURVIVAL ANALYSIS FUNCTIONS



#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# GENERAL PURPOSE FUNCTIONS

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# openStream(fileName, head)
# Function that opens a stream <fileName> on which R stduot is redirected.

openStream<-function(fileName, head){
  sink(fileName, append=TRUE, split=TRUE)
  if(!missing(head)) { 
    
    cat(paste("--------------------", head, "--------------------", " ", sep="\n"))
    
  }
}

# closeStream()
# Function that closes a stream previously opened on which R stduot is redirected.

closeStream<-function(){
  sink()
}


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SURVIVAL ANALYSIS FUNCTIONS

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

calc_age_at_onset<-function(born_date, onset_date){
  #born_date, onset_date vector (dset[,column])
  bd<-matrix(as.numeric(unlist(strsplit(as.character(born_date), split=c("/"), fixed=T))), ncol=3, byrow=T)
  od<-matrix(as.numeric(unlist(strsplit(as.character(onset_date), split=c("/"), fixed=T))), ncol=3, byrow=T)
  age_at_onset<-matrix(0, dim(bd)[1], 1)
  for(i in 1:dim(bd)[1]){
    a<-od[i,3]-bd[i,3]
    if(a<0){a=-1}
    if(od[i,2]<bd[i,2]) {a<-a-1}
    if(od[i,2]==bd[i,2]) {if(od[i,1]<bd[i,1]) {a<-a-1}}
    age_at_onset[i,1]<-a
  }
  return(age_at_onset)
}

calc_follow_up_monthes<-function(onset_date, last_fu_date, at_monthes){
  #onset_date, last_fu vector (dset[,column])
  od<-matrix(as.numeric(unlist(strsplit(as.character(onset_date), split=c("/"), fixed=T))), ncol=3, byrow=T)
  fd<-matrix(as.numeric(unlist(strsplit(as.character(last_fu_date), split=c("/"), fixed=T))), ncol=3, byrow=T)
  od[,3]<-as.numeric(substr(od[,3], nchar(od[,3])-1, nchar(od[,3])))
  fd[,3]<-as.numeric(substr(fd[,3], nchar(fd[,3])-1, nchar(fd[,3])))
  follow_up_monthes<-matrix(0, dim(od)[1], 1)
  for(i in 1:dim(od)[1]){
    m<-fd[i,3]-od[i,3]
    
    if(m<0) 
      m<-(100+fd[i,3]-od[i,3])
    
    if(m==0) m<-fd[i,2]-od[i,2]
    else {
      if(od[i,2]>fd[i,2]) m<-(m-1)*12+(12-od[i,2])+fd[i,2]
      else m=m*12+(fd[i,2]-od[i,2])
    }
    if(!is.na(at_monthes)) {
      if(m>at_monthes) follow_up_monthes[i,1]<-at_monthes
      else follow_up_monthes[i,1]<-m
    }
    else follow_up_monthes[i,1]<-m
  }
  return(follow_up_monthes)
}

calc_dead_at_monthes<-function(onset_date, death_date, at_monthes){
  #onset_date, last_fu vector (dset[,column])
  dcol<-death_date
  dcol[dcol==""]<-"0/0/0"
  od<-matrix(as.numeric(unlist(strsplit(as.character(onset_date), split=c("/"), fixed=T))), ncol=3, byrow=T)
  dd<-matrix(as.numeric(unlist(strsplit(as.character(dcol), split=c("/"), fixed=T))), ncol=3, byrow=T)
  od[,3]<-as.numeric(substr(od[,3], nchar(od[,3])-1, nchar(od[,3])))
  dd[,3]<-as.numeric(substr(dd[,3], nchar(dd[,3])-1, nchar(dd[,3])))
  dead_at_monthes<-matrix(0, dim(od)[1], 1)
  
  for(i in 1:dim(od)[1]){
    if((dd[i,1]!=0) & (dd[i,2]!=0) & (dd[i,3]!=0)){
      
      m<-dd[i,3]-od[i,3]
      
      if(m<0) m<-(100+dd[i,3]-od[i,3])
      
      if(m==0) m<-dd[i,2]-od[i,2]
      else {
        if(od[i,2]>dd[i,2]) m<-(m-1)*12+(12-od[i,2])+dd[i,2]
        else m=m*12+(dd[i,2]-od[i,2])
      }
    }
  }
  return(dead_at_monthes)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTIONS TO PREPARE BED FILES

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createRepTable<-function(repTab){
  
  chr<-matrix(unlist(strsplit(as.character(repTab[,3]), split=c(":"), fixed=T)), ncol=2, byrow=T)
  
  chr[,1]<-gsub("chr", "", chr[,1])
  
  chr[which(chr[, 1]=="X"), 1]<-23
  
  chr[which(chr[, 1]=="Y"), 1]<-24
  
  pos<-matrix(unlist(strsplit(as.character(chr[,2]), split=c("-"), fixed=T)), ncol=2, byrow=T)
  
  pos[,1]<-gsub(",", "", pos[, 1])
  
  pos[,2]<-gsub(",", "", pos[, 2])
  
  cnv<-repTab[,4]
  
  cnv[cnv=="CN Gain"]<-"gain"
  
  cnv[cnv=="CN Loss"]<-"loss"
  
  cnv[cnv=="High Copy Gain"]<-"duplication"
  
  cnv[cnv=="Homozygous Copy Loss"]<-"deletion"
  
  repTab<-data.frame(as.numeric(chr[,1]), as.numeric(pos[,1]), as.numeric(pos[,2]), cnv, repTab$Sample)
  
  colnames(repTab)<-c("CHR", "START", "END", "CNV_TYPE", "SAMPLE")
  
  return(repTab)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Friedman with post-hoc test

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

friedman.test.with.post.hoc <- function(formu, data, to.print.friedman = T, to.post.hoc.if.signif = T,  to.plot.parallel = T, to.plot.boxplot = T, signif.P = .05, color.blocks.in.cor.plot = T, jitter.Y.in.cor.plot =F)
{
	# formu is a formula of the shape: 	Y ~ X | block
	# data is a long data.frame with three columns:    [[ Y (numeric), X (factor), block (factor) ]]
	
	# Note: This function doesn't handle NA's! In case of NA in Y in one of the blocks, then that entire block should be removed.


	# Loading needed packages
	if(!require(coin))
	{
		print("You are missing the package 'coin', we will now try to install it...")
		install.packages("coin")		
		library(coin)
	}

	if(!require(multcomp))
	{
		print("You are missing the package 'multcomp', we will now try to install it...")
		install.packages("multcomp")
		library(multcomp)
	}
	
	if(!require(colorspace))
	{
		print("You are missing the package 'colorspace', we will now try to install it...")
		install.packages("colorspace")
		library(colorspace)
	}

	
	# get the names out of the formula
	formu.names <- all.vars(formu)
	Y.name <- formu.names[1]
	X.name <- formu.names[2]
	block.name <- formu.names[3]
	
	if(dim(data)[2] >3) data <- data[,c(Y.name,X.name,block.name)]	# In case we have a "data" data frame with more then the three columns we need. This code will clean it from them...

	# Note: the function doesn't handle NA's. In case of NA in one of the block T outcomes, that entire block should be removed.

	# stopping in case there is NA in the Y vector
	if(sum(is.na(data[,Y.name])) > 0) stop("Function stopped: This function doesn't handle NA's. In case of NA in Y in one of the blocks, then that entire block should be removed.")
	
	# make sure that the number of factors goes with the actual values present in the data:
	data[,X.name ] <- factor(data[,X.name ])
	data[,block.name ] <- factor(data[,block.name ])
	number.of.X.levels <- length(levels(data[,X.name ]))
	if(number.of.X.levels == 2) { warning(paste("'",X.name,"'", "has only two levels. Consider using paired wilcox.test instead of friedman test"))}
	
	# making the object that will hold the friedman test and the other.
	the.sym.test <- symmetry_test(formu, data = data,	### all pairwise comparisons	
						   teststat = "max",
						   xtrafo = function(Y.data) { trafo( Y.data, factor_trafo = function(x) { model.matrix(~ x - 1) %*% t(contrMat(table(x), "Tukey")) } ) },
						   ytrafo = function(Y.data){ trafo(Y.data, numeric_trafo = rank, block = data[,block.name] ) }
						)
	# if(to.print.friedman) { print(the.sym.test) }

	
	if(to.post.hoc.if.signif)
		{
			if(pvalue(the.sym.test) < signif.P)
			{
				# the post hoc test
				The.post.hoc.P.values <- pvalue(the.sym.test, method = "single-step")	# this is the post hoc of the friedman test
									

				# plotting
				if(to.plot.parallel & to.plot.boxplot)	par(mfrow = c(1,2)) # if we are plotting two plots, let's make sure we'll be able to see both
									
				if(to.plot.parallel)
				{
					X.names <- levels(data[, X.name])
					X.for.plot <- seq_along(X.names)
					plot.xlim <- c(.7 , length(X.for.plot)+.3)	# adding some spacing from both sides of the plot

					if(color.blocks.in.cor.plot) 
					{
						blocks.col <- rainbow_hcl(length(levels(data[,block.name])))
					} else {
						blocks.col <- 1 # black
					}					

					data2 <- data
					if(jitter.Y.in.cor.plot) {
						data2[,Y.name] <- jitter(data2[,Y.name])
						par.cor.plot.text <- "Parallel coordinates plot (with Jitter)"				
					} else {
						par.cor.plot.text <- "Parallel coordinates plot"
					}				
					
					# adding a Parallel coordinates plot
					matplot(as.matrix(reshape(data2,  idvar=X.name, timevar=block.name,
									 direction="wide")[,-1])  , 
							type = "l",  lty = 1, axes = FALSE, ylab = Y.name, 
							xlim = plot.xlim,
							col = blocks.col,
							main = par.cor.plot.text)
					axis(1, at = X.for.plot , labels = X.names) # plot X axis
					axis(2) # plot Y axis
					points(tapply(data[,Y.name], data[,X.name], median) ~ X.for.plot, col = "red",pch = 4, cex = 2, lwd = 5)
				}
				
				if(to.plot.boxplot)
				{
					# first we create a function to create a new Y, by substracting different combinations of X levels from each other.
					subtract.a.from.b <- function(a.b , the.data)
					{
						the.data[,a.b[2]] - the.data[,a.b[1]]
					}
					
					temp.wide <- reshape(data,  idvar=X.name, timevar=block.name,
									 direction="wide") 	#[,-1]
					wide.data <- as.matrix(t(temp.wide[,-1]))
					colnames(wide.data) <- temp.wide[,1]
						
					Y.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, subtract.a.from.b, the.data =wide.data)
					names.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, function(a.b) {paste(a.b[2],a.b[1],sep=" - ")})
					
					the.ylim <- range(Y.b.minus.a.combos)
					the.ylim[2] <- the.ylim[2] + max(sd(Y.b.minus.a.combos))	# adding some space for the labels
					is.signif.color <- ifelse(The.post.hoc.P.values < .05 , "green", "grey")
					
					boxplot(Y.b.minus.a.combos,
						names = names.b.minus.a.combos ,
						col = is.signif.color,
						main = "Boxplots (of the differences)",
						ylim = the.ylim
						)
					legend("topright", legend = paste(names.b.minus.a.combos, rep(" ; PostHoc P.value:", number.of.X.levels),round(The.post.hoc.P.values,5)) , fill =  is.signif.color )
					abline(h = 0, col = "blue")
					
				}
				
				list.to.return <- list(Friedman.Test = the.sym.test, PostHoc.Test = The.post.hoc.P.values)
				if(to.print.friedman) {print(list.to.return)}				
				return(list.to.return)
				
			}	else {
					print("The results where not significant, There is no need for a post hoc test")
					return(the.sym.test)
				}					
	}
}
