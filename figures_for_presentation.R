# figures created for presentation
library(reshape2)
library(ggplot2)
x<-rnorm(10000)
b0<-rnorm(10000)
b1<-rnorm(10000,2,1)
y<-b0+b1*x

plot(x,y)
df<-data.frame(x,y)
summary(lm(y~x))
ggplot(df)+geom_point(aes(x=x,y=y),alpha=0.5)+theme_bw(base_size=17)+theme(aspect.ratio=1)+geom_smooth(aes(x=x,y=y),method="lm")


x<-rnorm(10000,2)
y<-rpois(10000,2)
z<-rgamma(10000,2)
q<-runif(10000,0,7)
df<-data.frame(x,y,z,q)
df<-melt(df)
ggplot(df)+geom_density(aes(value,fill=variable),alpha=0.75)+theme_bw(base_size=17)+theme(aspect.ratio=1)+scale_fill_discrete(name="Distribution",labels=c("Normal","Poisson","Gamma","Uniform"))+theme(legend.position=c(1,0),legend.justification=c(1,0))