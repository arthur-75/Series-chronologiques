#intrpoudtuin 

#import les donnés 
setwd("~/School/Miashs/MIASHS 6/seris chro/final project")

tmp<- read.csv("Data/tmp.csv")
#################################################################################
####################################################################################

#deviser les donnée 
last_one <- length(tmp$temperature) #de 1 à 2448
tmp[is.na(tmp$temperature),c('X','temperature')] #show all the missing values 
first_part<-c(1,1336) #premier section de donnés 
last_part<- c(1369,last_one) #deuxeme section de donnés 
tset_first<- c(1,(first_part[2]-8)) #test how good is our model for the first part
test_second <- c((last_part[1]+8),last_one)#test how good is our model for the second part
################################################################################## 
####################################################################################

#tout d'abord pour prendre une idee genral de ce qui pass on fait un visulastion des donnés



temperature <- xts::xts(tmp$temperature, order.by=as.POSIXct(tmp$dateheure))
plot(temperature)
#apres a prendre une idée général on doit etre plus presier de chartiser ces donnés
#on peut sentir qu'il y une sesanolité mais doit zommer plus pres pour etre sur 
#et comme dans l'annceé dit que tout les jours on prede 8 fois la dergre de temperteure

   
#cette fonction eveite la repetion de code pour chauqe fois on apple une partie de donnés
#et aussi s'on utlise la differnce ou une tendece ou un varaince  
main_data <- function(inverse=F,log=0,diff=0,lag=0,de_a=c(1,2448),data=tmp){
  P<-data$dateheure[de_a[1]:de_a[2]]
  X<-data$temperature[de_a[1]:de_a[2]]
  
  if (log){ 
    X<-log(X)}
  if (inverse){
    X<-rev(X)} #inverser les donnés
  if (lag){
    X<-diff(X,lag = lag)} #mettre un lag au donnés
  if (diff){
    X<- diff(X,differences = diff)}
  P<-data$dateheure[(de_a[1]+lag+diff):de_a[2]]
  Temperature<-  xts::xts(X, order.by=as.POSIXct(P))
  plot_<<- plot(Temperature)
  return(X)}
  
#voir s'il y a un sesanolite plus pret
#on predre deux les donnés de 2 jour ou 3 jour 
data_2j<- main_data(de_a=c(1,16))
plot_

#on voit une repettion chaque jour avec des petites varabilité 
#et s'on fait pout 4 jour 
data_4j<- main_data(de_a=c(1,32))
plot_X

#et encore on peut voir la autocorlation 
acf(tmp$temperature ,plot = TRUE, na.action = na.pass,lag.max = 32)
#ce graphique confirm la sesanlité puisque on a repertuation tout 8 fois
####################################################################################
####################################################################################

#de plus on voit pas des tendence donc on n a pas besoin de utliser
#la fonction log aux tous les donnés 
####################################################################################
####################################################################################

#test la stastionatité  
library("tseries")
#pour tester la stasionarite on va utliser deux test difrrentes 
#la premper est ADF(Test de Dickey-Fuller) : 

#la deuxeme est KPSS( Test de KPSS) :



options(warn=-1) #Ca aide de ne pas afficher "Warning message:"
test_ <- function(data){
  ADF <- adf.test(data)  
  KPSS <- kpss.test(data, null="Level")
  return(sprintf("le test ADF donne %s pour stastionarité et le test KPSS donne %s pour la stastionarité",
                 ADF$p.value < 0.05 , KPSS$p.value >= 0.05 ))
}
#pour bien affecut les test on doit deviser les donnés en deux parties commme
#il y a un manque de donné au millue ca sera beacoup plus facile à treter les donnés
#sperarment 

X_1 <- main_data(de_a=tset_first,lag=8)
plot_
acf(X_1)
#on vois les L'ACF

#on vois que c'est deje mieux avec lag de 8 
#et on pass les tests :
test_(X_1)
#on vois que c'est bien stastionare'
#on pass rappidement seulment avec moving avarage methode  

#tout d'abord on va essyer de trouver le meilleure d reapterue retard
#construit notre une fonction pour menemiser l'AIC

min_aic<-function(X,P,D,Q) {
  min_p=P[1]
  min_d=D[1]
  min_q=Q[1]
  cpt=0
  for(d in D[1]:D[2]) {
    for(p in P[1]:P[2]) {
      for(q in Q[1]:Q[2] ) {
        tryCatch( {
          Mod<- arima(X,order = c(p,d,q),method='ML')
          Mod____<<-Mod
          aic=Mod$aic
        },error=function(e){})
        if(cpt==0){
          aic_min=aic
          print(aic_min)}
        else{
          if(aic<aic_min){
            aic_min=aic
            min_p=p
            min_d=d
            min_q=q
            print(sprintf('Tour %s et les orders sont : %s , %s , %s et aic est %1f :',cpt,min_p,min_d,min_q,aic_min))
        }}
        cpt = cpt+1
      }
    }
    print("yes")
  }
  aic_<<-aic_min
  return(c(min_p,min_d,min_q))
}
#"0,2,30Tour 30 et les orders sont : 0 , 0 , 30 et aic est 4882.529617 :"
#"10,0,10 Tour 108 et les orders sont : 9 , 0 , 9 et aic est 4875.599745 :"
# "10,2,10"Tour 54 et les orders sont : 4 , 1 , 10 et aic est 4890.488702 :"
#"20,2,20 "Tour 30 et les orders sont : 12 , 0 , 11 et aic est 4874.780577 :"
#pour_model<- min_aic(X_1,c(0,0),c(0,9),c(0,30))

MA_model <-arima(X_1,order=c(0,0,30))


prediction_ <- function(model,pr_=8,data_n=tset_first[2],sta_p=100) {
  MA_predict<- predict(model , pr_ )
  data_<-data_n
  MA_X<-MA_sup<-MA_inf<- tmp$temperature
  for(i in (data_+1):(data_+pr_)){
    MA_X[i]<-MA_predict$pred[i-data_]+MA_X[i-8]
    MA_inf[i]<-MA_predict$pred[i-data_]-qnorm(0.975)*MA_predict$se[i-data_]+MA_X[i-8]
    MA_sup[i]<-MA_predict$pred[i-data_]+qnorm(0.975)*MA_predict$se[i-data_]+MA_X[i-8]
  }
#apres predire linfeire et superurier d'un entrvale de confince de 95%
# on peut les tracer sur la graphi mere et aussi afficher les vrai bon valures aussi 
  end_p<- pr_+data_
  sta_p<- sta_p
  plot(tmp$temperature[sta_p:(end_p)],
      type = 'l',col='blue',
      main="Predction de temperture avec un intervalle de confiance de niveau 95%"
      ,ylab="Temperatures",xlab="Jours",
      ylim = c(min(MA_inf[sta_p:end_p]),
      max(MA_sup[sta_p:end_p])),panel.first = grid(20))
  lines(MA_X[sta_p:end_p],type = 'l',col='orange',lty=3,lwd=3)
  lines(MA_sup[sta_p:end_p],type='l',col='red',lty=3,lwd=2)
  lines(MA_inf[sta_p:end_p],type='l',col='red',lty=3,lwd=2)
  lines(tmp$temperature[sta_p:end_p],lwd=3,
       type = 'l',col='blue')
  MES<- (1/8)*sum((MA_X[(data_+1):end_p]-tmp$temperature[(data_+1):end_p])^2)
  if (!is.na(MES)){
    print(MES)}

}
prediction_(MA_model,8,tset_first[2],1327)
#on vois que nous somme assez proche de la bonne reponse et l'intervalle de confince 
#bien juste 
#maintantit on predre la completment la partie 1 
X_2 <- main_data(de_a=first_part,lag=8)
MA_model_2 <-arima(X_1,order=c(0,0,30))
#avec le modele complet on trouve l'aic =4882.53 qui ne change pas vraiment 

prediction_(MA_model_2,32,first_part[2],1330)


#Alors on attack  autorefersie 















ARIMA_X_8 <- main_data(inverse=T,de_a=test_second)
test_(ARIMA_X_8)

#pour_model_<- min_aic(ARIMA_X,c(0 ,20),c(0,2),c(0,20))
library(MASS)
library(tseries)
library(forecast)
?auto.arima(ARIMA_X_8,d=0,S=8)













  