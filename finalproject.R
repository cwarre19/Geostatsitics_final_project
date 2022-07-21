options(noaakey = "FsbwNbQyoIOHLnPdaNSZJQZquVXxakKu")
library(rnoaa)
library(tidyverse)



get_data <- function(x){
  data <- ncdc(datasetid='NORMAL_DLY', limit=100, datatypeid='dly-tavg-normal', startdate = '2010-01-01', enddate = '2010-03-31',  add_units=TRUE, locationid=x$id)
  data$data %>% as_tibble() %>% mutate(city=x$name, cityid=x$id)
}

#collect data on the mean daily temperature for all cities in California with a NOAA temperature monitor. 
datafile <- "~/Desktop/school/geostats/data_final.rds"

if (file.exists(datafile)){
  data <- readRDS(datafile)
} else{
  cities <- ncdc_locs(locationcategoryid = "CITY", sortfield = "name", sortorder='desc')
  ca.cities <- cities$data %>% as_tibble() %>% filter(grepl("CA US", name))
  data <- lapply(1:dim(ca.cities)[1], function(x) get_data(ca.cities[x,]))
  data <- do.call(rbind, data) %>% as_tibble()
}

#collect latitude, longitude data on the stations within the cities
latlong <- readRDS("latlong.rds")

#merge data types
df <- latlong %>% left_join(data, by=c("city"="city")) %>% as_tibble()

df <- df %>% select(city, lat, long, date,value) %>% mutate(lat = as.numeric(lat),long = as.numeric(long))

#let us use a location that is somewhat far from all California places. Let us choose New York City, New York
ref.loc <- c(40.7128, 74.0060)

#let us compute the great circle distance because we know that a 1 latitude change is not always the same amount of distance 
source("greatcircledist.r")

gcd <- unlist(apply(df, 1, function(x) greatcircledist(ref.loc[1],as.numeric(x[2]), ref.loc[2],as.numeric(x[3]))))
df$gcd <- gcd
#these values are really big, let us log base 10 transform them
df <- df %>% select(-c(lat, long)) %>% mutate(gcd = log(gcd, base=10)) %>% mutate(value = as.numeric(value))

library(geoR)

#let us now fit some variograms
#we are interested in one month's worth of data at a time, so let us split by month
#and summarize by taking the mean at each month
df$month <- str_sub(df$date, 6, 7)
df <- df %>% group_by(month, city) %>% mutate(mean=mean(value)) %>% ungroup() 

gdat <- df %>% select(-c(date, value)) %>% distinct()
split.data <- split(gdat, gdat$month)

#let us look at each month separately (change the index on split.data) for now (Think about how to best display data for all 3 months)
summ.df <- split.data[[1]] %>% select(-city) %>% distinct() %>% select(-month) %>% mutate(dummy=0)
gdata <- as.geodata(summ.df, coord.col=c(1,3), data.col = 2)

#raw variogram
vcloud = variog(gdata, option = "cloud")
plot(vcloud, main = "Variogram Cloud")

#let us now look at an experimental variogram
max.dist = vcloud$max.dist/2
#Now, Kitanidis' book recommends 3-6 bins, and Dr. Miller recommends 5-6 bins. 
#But, for sake of being thorough, let us look at the variograms when we iterate from 3 to 6 bins inclusive.   
nbins <- 3:6
par(mfrow=c(2,2))
for (i in nbins){
  uvec = i
  vexp = variog(gdata, option = "bin",max.dist=max.dist, uvec = uvec)
  
  plot(vexp, main=paste0("Bins: ", i), pch=19, cex=1.5, col="blue")
}
#I like 6 bins

#what variogram model to choose? Clearly, this is nonstationary! (what is intrinsic nonstationary?)
#the two variogram models that model intrinsic nonstationarity are the linear and power model 
#I think the power will do better, the variogram is clearly not linear


#let us visually play with the parameters in the power model to find good inital guesses for our parameters of choice
power <- function(x, b, c){
  x^c * b
}
max.dist = vcloud$max.dist/2
vexp = variog(gdata, option = "bin",max.dist=max.dist, uvec = 6)
xvec <- seq(0, max.dist)

plot(vexp)
lines(xvec, power(xvec, b=1.2, c=1.8))

#we can fit the power model parameters b and c 
#(0 < c < 2 by definition)
#Restricted Maximum Likelihood and Weighted Least Squares)

b.guess <- 1.2
c.guess <- 1.8

weightsoption = "npairs"
initialguess = c(1.2, 1.8)
covmodel = "power"
vario.WLS = variofit(vexp, ini.cov.pars = initialguess, weights = weightsoption, cov.model = covmodel)

#the WLS estimate for b is 0.5755  and c is 1.9627

plot(vexp)
lines(xvec, power(xvec, 1.2, 1.8), col="red")
lines(vario.WLS, col="blue")
legend("topleft", 
       legend = c("Visual Fit",  "WLS"), 
       col = c("red", "blue"), 
       pch = c(19,19), 
       bty = "n", 
       pt.cex = 1.1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))



#Now, for estimation we will do a leave one out approach. 
#Namely, we will leave out city j, fit the variogram and do the kriging/OLS to estimate monthly temp at city j

#loc is the value for which we want to estimate the response for
#data is the data frame
idw <- function(loc, data){
  vals <- matrix(data$mean, ncol=1)
  p <- 2
  distance <- loc - data$gcd
  pwr <- distance ^ (-p)
  total.sum <- sum(pwr)
  lambdas <- pwr/total.sum
  lambdas <- matrix(lambdas, nrow=1)
  lambdas %*% vals
}

cov.model <- "power"
get_est <- function(tmp.df, j) {
  temp <- tmp.df[-j, ]
  
  x <- tmp.df[j, 2]
  city.guess <- tmp.df[j, 1]
  z <- tmp.df[j, 4]
  
  summ.df <- temp %>% select(-city) %>% distinct() %>% select(-month) %>% mutate(dummy=0)
  gdata <- as.geodata(summ.df, coord.col=c(1,3), data.col = 2)
  
  vcloud = variog(gdata, option = "cloud")
  max.dist = vcloud$max.dist/2
  vexp = variog(gdata, option = "bin", max.dist=max.dist, uvec = 6)
  vario.WLS = variofit(vexp, ini.cov.pars = initialguess, weights = weightsoption, cov.model = cov.model)
  
  #linear regression
  m <- lm(mean ~ gcd, temp)
  
  #fit, lower, upper 95% CI
  lm.guess <- predict.lm(m, x, interval="confidence") %>% unname()
  
  #idw(x %>% unname() %>% unlist(), temp)
  
  #kriging vars are 0 which makes sense because power model at distance 0 gives covariance of 0; nugget is estimated to be about 0 by variogram fitting (as we can see w our eyes too)
  sk.guess <- krige.conv(gdata, locations = data.frame(gcd = x %>% unname()  %>% unlist(), dummy=0),krige = krige.control(type.krige = "SK", cov.model=cov.model, beta = mean(temp$gcd), obj.model = vario.WLS) )$predict[[1]]
  ok.guess <- krige.conv(gdata, locations = data.frame(gcd = x %>% unname()  %>% unlist(), dummy=0),krige = krige.control(type.krige = "OK", cov.model=cov.model, obj.model = vario.WLS) )$predict[[1]]
  
  c(lm.guess[1], lm.guess[2], lm.guess[3], sk.guess, ok.guess, z %>% unname() %>% unlist())
  
}


rmse <- list()
res <- list()

for (k in 1:2){
  
  tmp.df <- split.data[[k]]
  n.cities <- dim(tmp.df)[1]
  
  if (k==1){
    #for january
    #at position six,Error in solve.default(Vcov, trend.d) : 
    #Lapack routine dgesv: system is exactly singular: U[36,36] = 0
    indexes <- c(1:n.cities)[-6]
  } else{
    indexes <- c(1:n.cities)
  }
  
  results <- do.call(rbind, lapply(indexes, function(j) get_est(tmp.df, j) )) %>% as_tibble()
  
  colnames(results) <- c("lm", "lm.low", "lm.upper", "sk", "ok", "obs")
  
  if (k==1){
    results$city <- tmp.df$city[-6]
    results$month <- tmp.df$month[-6]
  } else{
    results$city <- tmp.df$city
    results$month <- tmp.df$month
  }
  
  r.long <- results %>% pivot_longer(-c(lm.low, lm.upper, month, city), names_to="type", values_to="val")
  r.long <- r.long %>% mutate(city = gsub(", CA US", "", city))
  
  r.long <- r.long %>% mutate(lm.low = ifelse(type=="lm", lm.low, val), lm.upper = ifelse(type=="lm", lm.upper, val))
  
  plt <- ggplot(r.long, aes(x=city, y=val, color=type, group=type)) + geom_point(size=2) + theme_bw()  + theme(axis.text.x = element_text(angle = 45))
  plt <- plt + theme(axis.title.y = element_text(size=15), axis.title.x = element_blank()) + ylab("Monthly Mean Temperature")
  plt <- plt + geom_errorbar(aes(ymin=lm.low, ymax=lm.upper), width=.2, position=position_dodge(0.05)) + labs(color="Method")
  
  
  #RMSE for each method
  n <- (dim(results)[1])
  lm.rmse <- sum ( (results$lm - results$obs)^2 )  / n
  ok.rmse <- sum ( (results$ok - results$obs)^2 )  / n
  sk.rmse <- sum ( (results$sk - results$obs)^2 )  / n
  
  if (k==1){
    plt <- plt + ggtitle("Month: January") + theme(plot.title=element_text(hjust=0.5, size=16)) + theme(legend.position="none")
    month.rmse <- data.frame(lm=lm.rmse, ok=ok.rmse, sk=sk.rmse, month="January")
  }
  
  if (k==2){
    plt <- plt + ggtitle("Month: February") + theme(plot.title=element_text(hjust=0.5, size=16))  + ylab("")
    month.rmse <- data.frame(lm=lm.rmse, ok=ok.rmse, sk=sk.rmse, month="February")
  }
  
  rmse[[k]] <- month.rmse
  res[[k]] <- plt  
}


#rmse comparisons
t <- do.call(rbind, rmse)
t <- t %>% pivot_longer(-c(month), names_to="method", values_to="val")
p <- ggplot(t, aes(x=month, y=val, color=method, group=method)) + geom_point(size=3) + theme_bw() + geom_line()
loocv.rmse <- p + xlab("") + ylab("RMSE in Leave One Out Setup")  + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=16))
ggsave(loocv.rmse, filename="plots/loocv.rmse.pdf", width=8, height=9,units="in", device=cairo_pdf )

#dimensions 8 x 6
#much less data in March! So we won't krig with it! 
par(mfrow=c(1,3))
 hist(split.data[[1]]$gcd, xlab="GCD", main="January")
 hist(split.data[[2]]$gcd, xlab="GCD", main="February")
 hist(split.data[[3]]$gcd, xlab = "GCD", main="March")


 #dimensions 21 x 10
library(cowplot)
fin.res <- plot_grid(res[[1]], res[[2]], nrow=1, rel_widths=c(1.75, 1.25))
ggsave(fin.res, filename="plots/fin.res.pdf", width=15, height=10,units="in", device=cairo_pdf )


#what if we make an approximation and say that the variogram model is stationary (say exponential in nature)
rmse <- list()
res <- list()

cov.model <- "exponential"
for (k in 1:2){
  
  tmp.df <- split.data[[k]]
  n.cities <- dim(tmp.df)[1]
  
  if (k==1){
    #for january
    #at position six,Error in solve.default(Vcov, trend.d) : 
    #Lapack routine dgesv: system is exactly singular: U[36,36] = 0
    indexes <- c(1:n.cities)[-6]
  } else{
    indexes <- c(1:n.cities)
  }
  
  results <- do.call(rbind, lapply(indexes, function(j) get_est(tmp.df, j) )) %>% as_tibble()
  
  colnames(results) <- c("lm", "lm.low", "lm.upper", "sk", "ok", "obs")
  
  if (k==1){
    results$city <- tmp.df$city[-6]
    results$month <- tmp.df$month[-6]
  } else{
    results$city <- tmp.df$city
    results$month <- tmp.df$month
  }
  
  r.long <- results %>% pivot_longer(-c(lm.low, lm.upper, month, city), names_to="type", values_to="val")
  r.long <- r.long %>% mutate(city = gsub(", CA US", "", city))
  
  r.long <- r.long %>% mutate(lm.low = ifelse(type=="lm", lm.low, val), lm.upper = ifelse(type=="lm", lm.upper, val))
  
  plt <- ggplot(r.long, aes(x=city, y=val, color=type, group=type)) + geom_point(size=2) + theme_bw()  + theme(axis.text.x = element_text(angle = 45))
  plt <- plt + theme(axis.title.y = element_text(size=15), axis.title.x = element_blank()) + ylab("Monthly Mean Temperature")
  plt <- plt + geom_errorbar(aes(ymin=lm.low, ymax=lm.upper), width=.2, position=position_dodge(0.05)) + labs(color="Method")
  
  
  #RMSE for each method
  n <- (dim(results)[1])
  lm.rmse <- sum ( (results$lm - results$obs)^2 )  / n
  ok.rmse <- sum ( (results$ok - results$obs)^2 )  / n
  sk.rmse <- sum ( (results$sk - results$obs)^2 )  / n
  
  if (k==1){
    plt <- plt + ggtitle("Month: January") + theme(plot.title=element_text(hjust=0.5, size=16)) + theme(legend.position="none")
    month.rmse <- data.frame(lm=lm.rmse, ok=ok.rmse, sk=sk.rmse, month="January")
  }
  
  if (k==2){
    plt <- plt + ggtitle("Month: February") + theme(plot.title=element_text(hjust=0.5, size=16))  + ylab("")
    month.rmse <- data.frame(lm=lm.rmse, ok=ok.rmse, sk=sk.rmse, month="February")
  }
  
  rmse[[k]] <- month.rmse
  res[[k]] <- plt  
}

library(cowplot)
fin.res <- plot_grid(res[[1]], res[[2]], nrow=1, rel_widths=c(1.75, 1.25))
ggsave(fin.res, filename="plots/exponential.fin.res.pdf", width=15, height=10,units="in", device=cairo_pdf )





