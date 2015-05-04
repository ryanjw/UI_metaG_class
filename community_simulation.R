# install these packages first
install.packages("vegan")
install.packages("ggplot2")
install.packages("reshape2")


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

df$rel_ab<-fseries(df$rank,.99)

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