############### MODELING POPULATION DYNAMIC ###############
###### TOPIC: The role of semelparity, delayed reproduction ######
####and density dependence on population dynamics#########

         
##' @aliases Student: ZINZINHEDO Mahoukpégo Luc
##'Prof: Orou GAOUE
##'Date: 14 September 2021
                    

##'Creating projection matrix function called 'projection'. 
#'
#'This function tackes the parameters: sigma1,sigma2,gama and phi. 
#'The output coming from projection are : 
#'(1) projection matrix, (2) lambda0(asymptotic growth rate), 
#'and (4) the elasticities
#'lambda0 and elasticites come from the package @source popbio

########################## The Begining ###########################
rm(list = ls(all=T))


projection <- function(sigma1,sigma2,gama,phi){
 A <- matrix(c(sigma1*(1-gama),phi,sigma1*gama, sigma2), 
              byrow = T, ncol = 2)
 if(sigma1>1|sigma1<0){stop("sigma1 must be bounded in 0 and 1")}
 if(gama>1|gama<0){stop("gama must be bounded in 0 and 1")}
 if(phi<0){stop("phi must be greater or equal to 0")}
library(popbio)
        e.a <- eigen.analysis(A) 
 lambda0 <- e.a$lambda1
Dynamic <- list(projection.matrix = A, lambda0=lambda0,
               Elasticity=e.a$elasticities)
 return(Dynamic)
}

#Try the function projection with B
B <- projection(0.5,0.9,0.1,1.5)


#Install package to store data in excel
install.packages("writexl")
library(writexl)
#######################'TASK 1##########################
##'                                                   ##
#' simulate the effects of delayed reproduction on the asymptotic(lambda0)
#' and transient (lambda1) dynamics of the study species.
Task1 <- function(Gama){
        n <- length(as.vector(Gama))
        g <- list()
        for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5,sigma2 =  0.9,
                                        gama = Gama[i],phi = 1.5)}
    return(g)
} 
maturation.process <- seq(from=0, to=1, length.out = 50) # We can modify it
resG <- Task1(maturation.process)
Lambda0 <- sapply(resG, "[[", "lambda0") 
n <- length(Lambda0)
v <- c()
a <- 1:length(Lambda0)

for (i in 1:length(Lambda0)){
 v[i] <- (1/(maturation.process[a[i+1]]-maturation.process[i]))*log(Lambda0[a[i+1]]/Lambda0[i])}
Lambda1 <- c(Lambda0[1], na.omit(v))
 v0 <- sapply(resG, "[[", 'Elasticity')
Elast_1.1 <- v0[1,] # the 1.1 element of the 2 * 2 elasticity matrix 
Elast_1.2 <- v0[2,] # the 1.2 element of the 2 * 2 elasticity matrix 
Elast_2.1 <- v0[3,] # the 2.1 element of the 2 * 2 elasticity matrix 
Elast_2.2 <- v0[4,] # the 2.2 element of the 2 * 2 elasticity matrix 

Task1.data <- data.frame(maturation.process, Lambda0, Lambda1, 
                         Elast_1.1, Elast_1.2,Elast_2.1, Elast_2.2)

# see the means of elasticities 
cbind(c("Elast_1.1", "Elast_1.2", "Elast_2.1","Elast_2.2"), 
      c(mean(Task1.data$Elast_1.1), mean(Task1.data$Elast_1.2), mean(Task1.data$Elast_2.1),
        mean(Task1.data$Elast_2.2)))
write_xlsx(Task1.data,"~/Maturation.process.xlsx") #Go to Document in your laptop and open the file called Maturation.process.xlsx
par(mfrow=c(1, 2))

plot(Task1.data$maturation.process,Task1.data$Lambda0,type = 'l',
       col='red', xlab ='Maturation rate',lwd=3,
       ylab = 'Asymptotic growth rate lambda0', 
       main = 'Asymptotic dynamic')
plot(Task1.data$maturation.process,Task1.data$Lambda1,
       col='blue', xlab ='Maturation rate', type = 'l', lwd=3,
        ylab = 'Transient growth rate Lambda1', 
        main = 'Transient dynamic')

#######################'TASK 2##########################
##                                                    ##
#'simulate the effects of iteroparity on the asymptotic
#'and transient dynamics of the study species.
#'
 Task2 <- function(Sigma2){
        n <- length(Sigma2)
        g <- list()
        for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5,sigma2=Sigma2[i],
                                        gama=0.1,phi = 1.5)}
        return(g)
} 
Survival.rate <- seq(from=0, to=1, length.out= 50)
resS <- Task2(Survival.rate)

Lambda0 <- sapply(resS, "[[", "lambda0") 
#n <- length(Lambda0)
v <- c()
a <- 1:length(Lambda0)

for (i in 1:length(Lambda0)){
    v[i] <- (1/(Survival.rate[a[i+1]]-Survival.rate[i]))*log(Lambda0[a[i+1]]/Lambda0[i])}
Lambda1 <- c(Lambda0[1], na.omit(v))
v0 <- sapply(resS, "[[", 'Elasticity')
Elast_1.1 <- v0[1,] # the 1.1 element of the 2 * 2 elasticity matrix 
Elast_1.2 <- v0[2,] # the 1.2 element of the 2 * 2 elasticity matrix 
Elast_2.1 <- v0[3,] # the 2.1 element of the 2 * 2 elasticity matrix 
Elast_2.2 <- v0[4,] # the 2.2 element of the 2 * 2 elasticity matrix 

Task2.data <- data.frame(Survival.rate, Lambda0, Lambda1, 
                         Elast_1.1, Elast_1.2,Elast_2.1, Elast_2.2)



# see the means of elasticities 
cbind(c("Elast_1.1", "Elast_1.2", "Elast_2.1","Elast_2.2"), 
      c(mean(Task2.data$Elast_1.1), mean(Task2.data$Elast_1.2), mean(Task2.data$Elast_2.1),
        mean(Task2.data$Elast_2.2)))
write_xlsx(Task2.data,"~/Survival.rate.xlsx") #Go to Document in your laptop and open the file called Survival.rate.xlsx


par(mfrow=c(1, 2))
plot(Task2.data$Survival.rate,Task2.data$Lambda0,type = 'l',
     col='red', xlab ='Survival rate',lwd=3,
     ylab = 'Asymptotic growth rate lambda0', 
     main = 'Asymptotic dynamic')
plot(Task2.data$Survival.rate,Task2.data$Lambda1,
     col='blue', xlab ='Survival rate', type = 'l', lwd=3,
     ylab = 'Transient growth rate Lambda1', 
     main = 'Transient dynamic')





#######################'TASK 3##########################
##                                                    ##
#' simulate the effect of over-compensatory density dependence on the 
#' asymptotic and transient population dynamics for semelparous 
#'(??2 = 0.1) and iteroparous (??2 = 0.9) species.
 
Pop.size <- 10:100 # Set the population size variating from 50 to 100

############## Semelparous
Task3.semelparous <- function(N){
        n <- length(N)
        g <- list()
        for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5, sigma2 = 0.1,
                                        gama= 0.1, phi = 1.5*exp(-N[i]))}
        return(g)
} 
resS.N <- Task3.semelparous(Pop.size)
Lambda0.s <- sapply(resS.N, "[[", "lambda0") 
#n <- length(Lambda0.s)
v <- c()
a <- 1:length(Lambda0.s)

for (i in 1:length(Lambda0.s)){
    v[i] <- (1/(Pop.size[a[i+1]]-Pop.size[i]))*log(Lambda0.s[a[i+1]]/Lambda0.s[i])}
Lambda1.s <- c(Lambda0.s[1], na.omit(v))
v0 <- sapply(resS.N, "[[", 'Elasticity')
Elast_1.1 <- v0[1,] # the 1.1 element of the 2 * 2 elasticity matrix 
Elast_1.2 <- v0[2,] # the 1.2 element of the 2 * 2 elasticity matrix 
Elast_2.1 <- v0[3,] # the 2.1 element of the 2 * 2 elasticity matrix 
Elast_2.2 <- v0[4,] # the 2.2 element of the 2 * 2 elasticity matrix 

Task3.data.s <- data.frame(Pop.size, Lambda0.s, Lambda1.s, 
                         Elast_1.1, Elast_1.2,Elast_2.1, Elast_2.2)
# see the means of elasticities 
cbind(c("Elast_1.1", "Elast_1.2", "Elast_2.1","Elast_2.2"), 
      c(mean(Task3.data.s$Elast_1.1), mean(Task3.data.s$Elast_1.2), mean(Task3.data.s$Elast_2.1),
        mean(Task3.data.s$Elast_2.2)))

write_xlsx(Task3.data.s,"~/Semelparous.xlsx") #Go to Document in your laptop and open the file called Semelparous.xlsx.xlsx
par(mfrow=c(1, 2))

plot(Task3.data.s$Pop.size,Task3.data.s$Lambda0.s,type = 'l',
     col='red', xlab ='Semelparous Poulation size',lwd=3,
     ylab = 'Asymptotic growth rate lambda0.s', 
     main = 'Asymptotic dynamic')
plot(Task3.data.s$Pop.size,Task3.data.s$Lambda1.s,
     col='blue', xlab ='Semelparous Population size', type = 'l', lwd=3,
     ylab = 'Transient growth rate Lambda1', 
     main = 'Transient dynamic')

############## Iteroparous
Task3.iteroparous <- function(N){
    n <- length(N)
    g <- list()
    for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5, sigma2 = 0.9,
                                      gama= 0.1, phi = 1.5*exp(-N[i]))}
    return(g)
}
resI.N <- Task3.iteroparous(Pop.size)
Lambda0.i <- sapply(resI.N, "[[", "lambda0") 
#n <- length(Lambda0.i)
v <- c()
a <- 1:length(Lambda0.i)

for (i in 1:length(Lambda0.i)){
    v[i] <- (1/(Pop.size[a[i+1]]-Pop.size[i]))*log(Lambda0.i[a[i+1]]/Lambda0.i[i])}
Lambda1.i <- c(Lambda0.i[1], na.omit(v))
v0 <- sapply(resI.N, "[[", 'Elasticity')
Elast_1.1 <- v0[1,] # the 1.1 element of the 2 * 2 elasticity matrix 
Elast_1.2 <- v0[2,] # the 1.2 element of the 2 * 2 elasticity matrix 
Elast_2.1 <- v0[3,] # the 2.1 element of the 2 * 2 elasticity matrix 
Elast_2.2 <- v0[4,] # the 2.2 element of the 2 * 2 elasticity matrix 

Task3.data.i <- data.frame(Pop.size, Lambda0.i, Lambda1.i, 
                           Elast_1.1, Elast_1.2,Elast_2.1, Elast_2.2)
# see the means of elasticities 
cbind(c("Elast_1.1", "Elast_1.2", "Elast_2.1","Elast_2.2"), 
      c(mean(Task3.data.i$Elast_1.1), mean(Task3.data.i$Elast_1.2), mean(Task3.data.i$Elast_2.1),
        mean(Task3.data.i$Elast_2.2)))

write_xlsx(Task3.data.i,"~/Iteroparous.xlsx") #Go to Document in your laptop and open the file called Iteroparous.xlsx.xlsx
par(mfrow=c(1, 2))

plot(Task3.data.i$Pop.size,Task3.data.i$Lambda0.i,type = 'l',
     col='red', xlab ='Iteroparous Poulation size',lwd=3,
     ylab = 'Asymptotic growth rate lambda0.s', 
     main = 'Asymptotic dynamic')
plot(Task3.data.i$Pop.size,Task3.data.i$Lambda1.i,
     col='blue', xlab ='Iteroparous Population size', type = 'l', lwd=3,
     ylab = 'Transient growth rate Lambda1', 
     main = 'Transient dynamic')



############################## The End ###################################
