# install these packages first
install.packages("vegan")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("TSA")
install.packages("plyr")

library(vegan)
library(ggplot2)
library(reshape2)
library(TSA)
library(plyr)
# Say you want to simulate a community

# We can use a fisher's series to produce a right-skewed distribution broken into discrete chunks, and then we will apply these relative abundances (probability that an OTU exists in a sample) to a set of OTU id's

fseries<-function(X,C){
	y<-(-1/log(1-C))*(C^X/X)
	return(y)
}

# let's generate a matrix with 10,000 total unique OTUs
otus<-matrix(ncol=1,nrow=0)
for(i in 1:10000){
	otus<-rbind(otus, paste(sample(letters,6),collapse=""))
}

rank<-seq(1,10000,1)
df<-data.frame(otus,rank)
head(df)
df$rel_ab<-fseries(df$rank,.99)
df$rel_ab<-df$rel_ab/sum(df$rel_ab)
#lets check our dataframe
head(df)
tail(df)

# lets look at our distribution of OTUs now
library(ggplot2)

ggplot(df)+geom_point(aes(x=rank,y=rel_ab))+scale_x_log10()+theme_bw()+theme(aspect.ratio=1)

# looks right skewed which is what we should see

# now lets do two more communities that are a bit more even

cs<-c(.99,.999,.9999)
final_df<-data.frame()
for(i in 1:length(cs)){
	df$rel_ab<-fseries(df$rank,cs[i])
	df$C<-cs[i]
	final_df<-rbind(final_df,df)
}

#lets make sure its all there
summary(final_df)	

# ok now plotting it
ggplot(final_df)+geom_point(aes(x=rank,y=rel_ab, colour=factor(C)))+scale_x_log10()+theme_bw()+theme(aspect.ratio=1)+scale_colour_manual(values=c("red","purple","blue"))


# we can see the difference in evenness by changing the scale of the y-axis
ggplot(final_df)+geom_point(aes(x=rank,y=rel_ab, colour=factor(C)))+scale_x_log10()+scale_y_log10()+theme_bw()+theme(aspect.ratio=1)+scale_colour_manual(values=c("red","purple","blue"))


# now we have three communities (or samples) that we can make into an otu table like we like to have
# lets turn it into a nxp matrix where rows are samples and columns are otu id's
library(reshape2)
 dataset<-dcast(final_df, C~otus,value.var="rel_ab")
 
 # let's check it
 head(dataset[,1:10])

# we can run basic community statistics on these communities
library(vegan)
specnumber(dataset[,-1])
diversity(dataset[,-1])
diversity(dataset[,-1])/log(specnumber(dataset[,-1]))
vegdist(dataset[,-1])

#####################################################################

# Challenges #

# Try to simulate 10 communities with a different C value, make a scatterplot of C vs. any diversity measure #

####################################################################

library(TSA)

# based on the Shade paper let's talk about the three different statistics that are generated

func<-function(x){
k=kurtosis(x)
s=skewness(x)
b=(1+(s^2))/(k+3)
return(data.frame(k,s,b))
}

# assume we sample a material 10000 times and are interested in an OTU whose abundance is normal with mu=10 and sd=2
otu<-rnorm(10000,mean=10,sd=2)
hist(otu)
func(otu)
# lets do it with a skewed distribution
otu<-rpois(10000,lambda=1)
hist(otu) 
func(otu)

# The distribution that is used to describe the OTU can change the measurement, which is exactly what the authors are interested in

# Let's sample a single OTU over time t={1,..,100} that peaks in abundance right in the middle t=50

# making time
t<-seq(1,100,1)

# making an abundance of an OTU that is distributed U(0,1)
df<-data.frame(t)
df$abundance<-round(runif(length(t),min=0,max=1))

# something at t=50, we want the otu to increase in abundance by an order of magnitude

df[50,2]<-df[50,2]+1*10

# now lets look at its abundance over time
ggplot(df)+geom_line(aes(x=t,y=abundance))

# lets use that function to see if its bimodal

func(df$abundance)
         # k        s         b
# 1 62.34336 7.221587 0.8134157

# Challenge # 

# Make the OTU increase in abundance at two different points

df<-data.frame(t)
df$abundance<-round(runif(length(t),min=0,max=1))
df[33,2]<-df[33,2]+1*10
df[50,2]<-df[50,2]+1*10
df[66,2]<-df[66,2]+1*10
ggplot(df)+geom_line(aes(x=t,y=abundance))
func(df$abundance)

# this changes b, the amount of bimodality to a level that is worth noting (CRT).


# Now we are going to test whether or not some way of handling the data (in this case, rarefaction) may bias CRT.  We are going to simulate a community that only varies based on rarefaction (in other words, we will resample the same OTU distribution over and over again).  We are simulating a community that is similar to the M3 right palm community as in the paper (8,230 OTUs, 5,031 reads per sample).  NOTE I am working on a publication related to this problem, so if you are interested in more info,  please contact me!

rank<-seq(1,8230,1)
df<-data.frame(rank)
df$rel_ab<-fseries(df$rank,.99)
# you have to scale the rel_ab acordingly (has to be sum = 1)
df$rel_ab<-df$rel_ab/sum(df$rel_ab)
ggplot(df)+geom_point(aes(x=rank,y=rel_ab))+scale_x_log10()
# to sample this data by rarefaction there are several functions, for the sake of this simulation, we are using some base functions in R

# to make one rarefied sample we only have to sample the ranks (i.e. otus) based on their relative abundances
rar_sample<-sample(df$rank,5031,replace=T,prob=df$rel_ab)
df_rar<-data.frame(count(rar_sample))
head(df_rar)
ggplot(df_rar)+geom_point(aes(x=x,y=freq/sum(freq)))+scale_x_log10()

# lets compare the rarefied sample with the true community that we sampled from
# should have richness of 8230
specnumber(df_rar[,1],MARGIN=2) # note that I changed the margin here because the dataframe is shaped differently
diversity(df[,2],MARGIN=2)
diversity(df_rar[,2],MARGIN=2)
diversity(df[,2],MARGIN=2)/log(10000)
diversity(df_rar[,2],MARGIN=2)/log(specnumber(df_rar[,1],MARGIN=2))

# lets look at the differences between the two 
m<-merge(df_rar,df,by.x="x",by.y="rank")
m$freq<-m$freq/sum(m$freq)
m_melt<-melt(m, id="x")
ggplot(m_melt)+geom_point(aes(x=x,y=value,colour=variable))+scale_x_log10()

#seems to be pretty close visually, yet we know there is some loss of information
#now we want to generate this same rarefaction step over a time series.  We are going to assume a time series as we did previously t={1,...,100}
df$rel_ab<-fseries(df$rank,.99)
df$rel_ab<-df$rel_ab/sum(df$rel_ab)
df_time<-data.frame()
for(t in 1:100){
	rar_sample<-sample(df$rank,5031,replace=T,prob=df$rel_ab)
	df_rar<-data.frame(count(rar_sample))
	df_rar$t<-t
	df_time<-rbind(df_time,df_rar)		
}

# changing names for ease
names(df_time)<-c("rank","reads","time")

# now running the analysis
ranks<-as.vector(unique(df_time$rank))
results_final<-data.frame()
for(r in 1:length(ranks)){
	#r<-250
	temp<-subset(df_time, rank==ranks[r])
	
	temp<-arrange(temp, time)
	temp_full<-data.frame(seq(1,100,1))
	names(temp_full)<-c("time")
	temp_merge<-merge(temp_full, temp, by="time",all.x=T)
	#head(temp_merge)
	for(d in 1:dim(temp_merge)[1]){
		if(is.na(temp_merge$reads[d])==T){temp_merge$reads[d]<-0}
	}
	results<-data.frame(func(temp_merge$reads))
	results$rank<-ranks[r]
	results_final<-rbind(results_final,results)
}


# now lets see howm many are CRT

ggplot(results_final)+geom_histogram(aes(b))+geom_vline(xintercept=0.9)
dim(subset(results_final, b > 0.9 ))[1]/dim(results_final)[1]

# what rank are those with high CRT
ggplot(results_final)+geom_point(aes(x=rank,y=b))+geom_hline(yintercept=0.9)

# what if we don't include singletons 

singletons<-ddply(df_time, .(rank),summarise, sample_count=length(time))
df_time_m<-merge(df_time,singletons,by="rank")
head(df_time_m)
df_time_no_singletons<-subset(df_time_m, sample_count > 1)
ranks<-as.vector(unique(df_time_no_singletons$rank))
results_final<-data.frame()
for(r in 1:length(ranks)){
	#r<-250
	temp<-subset(df_time_no_singletons, rank==ranks[r])
	
	temp<-arrange(temp, time)
	temp_full<-data.frame(seq(1,100,1))
	names(temp_full)<-c("time")
	temp_merge<-merge(temp_full, temp, by="time",all.x=T)
	#head(temp_merge)
	for(d in 1:dim(temp_merge)[1]){
		if(is.na(temp_merge$reads[d])==T){temp_merge$reads[d]<-0}
	}
	results<-data.frame(func(temp_merge$reads))
	results$rank<-ranks[r]
	results_final<-rbind(results_final,results)
}

# now lets see howm many are CRT

ggplot(results_final)+geom_histogram(aes(b))+geom_vline(xintercept=0.9)
dim(subset(results_final, b > 0.9 ))[1]/dim(results_final)[1]

# what rank are those with high CRT
ggplot(results_final)+geom_point(aes(x=rank,y=b))+geom_hline(yintercept=0.9)

# Challenge # 

# Simulate one of the other datasets used in the paper

rank<-seq(1,1816,1)
df<-data.frame(rank)
df$rel_ab<-fseries(df$rank,.99)
# you have to scale the rel_ab acordingly (has to be sum = 1)
df$rel_ab<-df$rel_ab/sum(df$rel_ab)
ggplot(df)+geom_point(aes(x=rank,y=rel_ab))+scale_x_log10()

df_time<-data.frame()
for(t in 1:100){
	rar_sample<-sample(df$rank,5134,replace=T,prob=df$rel_ab)
	df_rar<-data.frame(count(rar_sample))
	df_rar$t<-t
	df_time<-rbind(df_time,df_rar)		
}

# changing names for ease
names(df_time)<-c("rank","reads","time")

# now running the analysis
ranks<-as.vector(unique(df_time$rank))
results_final<-data.frame()
for(r in 1:length(ranks)){
	#r<-250
	temp<-subset(df_time, rank==ranks[r])
	
	temp<-arrange(temp, time)
	temp_full<-data.frame(seq(1,100,1))
	names(temp_full)<-c("time")
	temp_merge<-merge(temp_full, temp, by="time",all.x=T)
	#head(temp_merge)
	for(d in 1:dim(temp_merge)[1]){
		if(is.na(temp_merge$reads[d])==T){temp_merge$reads[d]<-0}
	}
	results<-data.frame(func(temp_merge$reads))
	results$rank<-ranks[r]
	results_final<-rbind(results_final,results)
}


# now lets see howm many are CRT

ggplot(results_final)+geom_histogram(aes(b))+geom_vline(xintercept=0.9)
dim(subset(results_final, b > 0.9 ))[1]/dim(results_final)[1]

# what rank are those with high CRT
ggplot(results_final)+geom_point(aes(x=rank,y=b))+geom_hline(yintercept=0.9)

# what if we don't include singletons 

singletons<-ddply(df_time, .(rank),summarise, sample_count=length(time))
df_time_m<-merge(df_time,singletons,by="rank")
head(df_time_m)
df_time_no_singletons<-subset(df_time_m, sample_count > 1)
ranks<-as.vector(unique(df_time_no_singletons$rank))
results_final<-data.frame()
for(r in 1:length(ranks)){
	#r<-250
	temp<-subset(df_time_no_singletons, rank==ranks[r])
	
	temp<-arrange(temp, time)
	temp_full<-data.frame(seq(1,100,1))
	names(temp_full)<-c("time")
	temp_merge<-merge(temp_full, temp, by="time",all.x=T)
	#head(temp_merge)
	for(d in 1:dim(temp_merge)[1]){
		if(is.na(temp_merge$reads[d])==T){temp_merge$reads[d]<-0}
	}
	results<-data.frame(func(temp_merge$reads))
	results$rank<-ranks[r]
	results_final<-rbind(results_final,results)
}

# now lets see howm many are CRT

ggplot(results_final)+geom_histogram(aes(b))+geom_vline(xintercept=0.9)
dim(subset(results_final, b > 0.9 ))[1]/dim(results_final)[1]

# what rank are those with high CRT
ggplot(results_final)+geom_point(aes(x=rank,y=b))+geom_hline(yintercept=0.9)
