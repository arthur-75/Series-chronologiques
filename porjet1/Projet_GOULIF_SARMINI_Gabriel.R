#M1 ISEFAR - Series chronologiques TP 

#Auteurs : Gabriel, Goulif, Sarmini 

#WorkyDirectory à changer
setwd("")

#libraries
library(ggpmisc)
library(tseries)




#Sinistres
load("Sinistres.Rdata")

#Question 1
# total des frais vétérinaires engagés par les assurés chaque mois (en million d'euros)
Sinistres
str(Sinistres)

#visualisation
plot(Sinistres)

time(Sinistres) #12 mois de 2005 à 2020

#meilleure visualisation avec ggplot
ggplot(Sinistres, as.numeric = FALSE) + geom_line()+
  geom_line(colour = "blue")+xlab("Année")+ 
  ylab("Cout sinitres en Million")+
  ggtitle("Les Sinistres de 2005 à 2020")

#on remarque une tendance, mais pas de saisonnalité, le modele semble additif
#(amplitude 'constante', en effet la composante tendancielle et le bruit semblent rester constante au cours du temps) 
#et donc pertinent la prediction sera plus simple


#Question 2
source("boxcox.r")
# cette fonction prends en argument une série chronologique "x" un parametre "a" et retourner la série chronologique transformée 
#Essayons quelques transformations pour se rapprocher vraiment d'un modèle additif :
a = 0.2
Sinistres02 = boxcox(Sinistres,a)
plot(Sinistres02,main = "BoxCox Sinistres, a = 0.2")
b = 0.5 
Sinistres05 = boxcox(Sinistres,b)
plot(Sinistres05,main = "BoxCox Sinistres, a = 0.5")

A = seq(0,0.5,0.1)
n = length(Sinistres)
M = matrix(0,n,length(A))
#test également fait pour A = seq(0.5,1,0.1)
for (i in (1:length(A))){
  M[,i] = boxcox(Sinistres,A[i])
}

M = ts(M,start=start(Sinistres),frequency=frequency(Sinistres))

#visualisation de toutes les transformations
ts.plot(M,col=(1:length(A)))
legend('topleft',legend=paste("a=",as.character(A),sep=""),col=(1:length(A)),lty=1)
#On a voulu faire une transformation additive mais il n'y a pas eu de réel changement observé.
#il n'y a pas de difference entre BOXCOX a=0.1 et a=0.9 
#nous ne gardons pas la transformation, et nous travaillerons sur la série telle quelle

#Polynome degres 6
#ajustement de la tendance polynomiale 
lm.out2 = lm(Sinistres ~ I(time(Sinistres)) + I(time(Sinistres)^2) + I(time(Sinistres)^6), data = Sinistres)

#visualisation de la série et de la tendance ajustée
plot(Sinistres,main="Sinistres")
lines(as.vector(time(Sinistres)),lm.out2$fitted.values,col="blue") 


#Question 3
SinistresDetend = Sinistres-lm.out2$fitted.values #il s'agit ici des résidus car nous avons constaté/supposé qu'il n'y a pas de saisonnalité

source("SaisEst.r") 
#cette fonction permet d'extraire la saisonnalité de la série
#nous sommes passé par cette étape pour être sur de ne pas nous tromper sur nos suppositions
p 		= 12
List 	= SaisEst(SinistresDetend,p)
motif 	= List$motif
SEst 	= List$serie

Residus = SinistresDetend-SaisEst(SinistresDetend,p)$serie

par(mfrow=c(1,1))
plot(Residus) 
acf(Residus,lag = 60) #plus le lag est grand plus il y a de la corrélation

h = 50
Box.test(Residus,lag=h,type ="Ljung-Box")
#  p-value > 5%  : On ne rejette pas l'hypothese d'indépendance des résidus

#Question 4

#Choisir les ordres pour l'arima
#le processus MA(q) où q est l’ordre chronologique.
#Tout d’abord le clé principale pour bien modéliser est L’AIC est que plus l’AIC est petit plus la modèle sera bon. Alors oui, on doit minimiser l’AIC dans un modèle pour choisir le bon ordre de q et d où d est l’opérateur différence.
#Bien sur, on peut choisir q à partir de son autocorrélogramme alors c’est bien lié à notre observation et notre goût statistique, cependant j’ai lassé cette méthode de modélisation pour la section suivante autorégressif.
#Et bien, Pour minimiser l’AIC on peut créer deux boucles pour d et q et laisser le code tourner pour tester plusieurs ordres jusqu’à ce qu'elle trouve le plus petit AIC.
#La fonction est :

min_aic<-function(X,P,D,Q) { #X = donnes, P = ordre de la autoregressif, D = differention rendre données stationnaires , Q = ordre de moyenne mobile process arma(p,q) aic=
  min_p=P[1] #vecteur pour autoregressif 
  min_d=D[1] #vecteur pour diff 
  min_q=Q[1] #vecteur pour moyenne mobile  
  cpt=0
  for(d in D[1]:D[2]) { 
    for(p in P[1]:P[2]) {
      for(q in Q[1]:Q[2] ) {
        tryCatch( { #eviter l'errure 
          Mod<- arima(X,order = c(p,d,q),method='ML') #Modiliser avec un combinasion  
          Mod____<<-Mod 
          aic=Mod$aic  #stocker le AIC
        },error=function(e){}) 
        if(cpt==0){ #essyer la premire combinasion  
          aic_min=aic # et stocker le premier AIC
          print(aic_min)
          print(sprintf
                ('Tour %s et les orders sont : %s , %s , %s et aic est %1f :'
                  ,cpt,min_p,min_d,min_q,aic_min))}
        else{
          if(aic<aic_min){ #si la nouvelle combinasion a un AIC plus petite que la precednt 
            aic_min=aic #alors remblase L'AIC 
            min_p=p # et remblase p 
            min_d=d # et remblase d 
            min_q=q # et remblase q
            print(sprintf
                  ('Tour %s et les orders sont : %s , %s , %s et aic est %1f :'
                    ,cpt,min_p,min_d,min_q,aic_min))
          }}
        cpt = cpt+1
      }
    }
    print("Stil running")
  }
  aic_<<-aic_min
  print("Done")
  return(c(min_p,min_d,min_q))
}

pour_model<- min_aic(Sinistres,c(0,12),c(1,1),c(0,12)) #test de p = 0 a 12, q allant également de 0 à 12
pour_model
#d'apres L'AIC on trouve que les ordre 4 1 4 sont les meilleurs

model_sini <-arima(Sinistres,order = c(4,1,4),method='ML')

#On verifie si les residus predits correspendant à notre residus orginal
ggplot(SinistresDetend, as.numeric = FALSE) + geom_line()+
  geom_line(y=model_sini$residuals,colour='red',linetype="dashed")+
  xlab("Année")+ 
  ggtitle("Residus de sinistres en noir + Residus de notre modele en rouge")

#les residus de sinistres en noir
#les residus de notre modèle en rouge
#donc on trouve que c'est presque proche de des residus originaux

#Prediction
pred_sini<-predict(model_sini,12) #utilisation de notre modèle

sini_pred <- pred_sini$pred
ggplot(sini_pred) + geom_line()+
  xlab("Année")+ 
  ylab("Cout sinitres en Million")+
  ggtitle("Prediction de 2021")


#Question 5
#intervalle de variations à 95% 
sini_se <- pred_sini$se

#les 12 dernières valeurs de Sinistres
Sini_last = ts(tail(Sinistres,12),end=end(Sinistres),frequency = frequency(Sinistres))

#valeur predites
sini_val<-ts(c(tail(Sinistres,1),sini_pred),
             end =end(sini_pred),frequency = frequency(sini_pred))

#valeur predites supérieures 
sini_up<-ts(c(tail(Sinistres,1),sini_pred + qnorm(0.975)*sini_se),
            end =end(sini_pred),frequency = frequency(sini_pred))

#valeur predites inférieures
sini_down<-ts(c(tail(Sinistres,1),sini_pred -qnorm(0.975)*sini_se),
              end =end(sini_pred),frequency = frequency(sini_pred))

# Graphe
plot(Sinistres,
     type = 'l',col='blue',
     main="Prediction de sinistres avec un intervalle de confiance de 95%"
     ,ylab="Cout Sinistre",xlab="Année")
lines(sini_val,type = 'l',col='orange',lty=2,lwd=2)
points(sini_up, type = "l", col = 2, lty = 2)
points(sini_down, type = "l", col = 2, lty = 2)






#Primes
#total des primes mensuelles versés par les assurés.
load("Primes.Rdata")

#Question 1
Primes
str(Primes)
#visualisation
plot(Primes)
time(Primes) #12 mois de 2005 à 2020, pareil que sinistres (logique)
#meilleure visualisation avec ggplot
ggplot(Primes, as.numeric = FALSE) + geom_line()+
  geom_line(colour = "blue")+xlab("Année")+ 
  ylab("Primes mensuelles versées")+
  ggtitle("Les Primes de 2005 à 2020") 
#on constate qu'il y a une tendance, une saisonnalité de période 1an (même motif se repetant)
#cependant il n'y a pas la même amplitude, elle varie entre 2005 et 2020 (contrairement à ce qu'on a constaté pour les sinistres)


#Question 2
#on rend le modèle additif avec la fonction log, remarque : on aurait pu utiliser la fonction boxcox avec a = 0
logPrimes=log(Primes)
#Visualisation de la transformation
ggplot(logPrimes, as.numeric = FALSE) + geom_line()+
  geom_line(colour = "blue")+xlab("Année")+ 
  ylab("Modèle additif : Primes mensuelles versées")+
  ggtitle("Les Primes de 2005 à 2020")
#le modèle additif : il semble y avoir la meme amplitude en 2005 comme en 2020
#le modèle additif sera pertinent notamment pour les prédictions


#Question 3
#Ajustons la tendance polynomiale avec un polynome de degres 2 : 
#car la repartition des données semble se rapprocher d'un polynome de degre 2

T = time(logPrimes)
lm.primes = lm(logPrimes ~I(T) + I(T^2), data = logPrimes)

#Visualisation de LogPrimes et de la tendance ajustée
ggplot(logPrimes, as.numeric = FALSE) + geom_line()+
  geom_line(colour = "blue")+xlab("Année")+ 
  ylab("Primes mensuelles versées (log)")+
  ggtitle("Les Primes de 2005 à 2020(log)")+
  geom_line( y=lm.primes$fitted.values,colour = "red") 
#l'ajustement de la tendance semble correct


#Question 4
PrimesDetend = logPrimes-lm.primes$fitted.values #logprime sans la tendance (saison+résidus)
#visualisation de logprimes sans la tendance estimée
ggplot(PrimesDetend) + geom_line() 

p= 12 #car la saisonnalité est de période 1an soit 12 mois
saison= SaisEst(PrimesDetend,p)
motif 	= saison$motif
SEst 	= saison$serie

#Visualisation
plot(motif,type="l",main ="motif")
ts.plot(PrimesDetend,SEst,col=c("black","red")) 
legend('topleft',legend=c("PrimesDetend","SEst"),col=c("black","red"),lty=1)


#Question 5
ResidusPrimes = PrimesDetend-SEst #on enlève la saisonnalité (il ne reste que les résidus)

ggplot(ResidusPrimes) + geom_line()

#autocorrelation : correlogramme
acf(ResidusPrimes,lag = 10)
acf(ResidusPrimes,lag = 20)
acf(ResidusPrimes,lag = 30)

#test de non stationnarité 
#Test de KPSS: H0 = stationnarité 
kpss.test(ResidusPrimes)
#la p-value est égale à 0.0669 > 0.05 : on ne rejette pas l'hypothèse de stationnarité 
#au seuil 5% on peut dire que les residus sont stationnaires


#Question 6
h = 20
Box.test(ResidusPrimes,lag=h,type ="Ljung-Box")
#  p-value > 5%  : On ne rejette pas l'hypothese d'indépendance des résidus 


#Question 7
#Prediction
#Tendance :
newT = seq(2021,2022,(1/12))[-13] #création d'un vecteur contenant les composantes temps des 12 prochains mois
pred_tend=predict(lm.primes,   newdata=data.frame(T=newT)) 
pred_tend_ts= ts(pred_tend,start=c(2021,1),frequency = 12)
ts.plot(logPrimes,pred_tend_ts,col=c("blue","orange"),main="Prédiction tendance")

#Saison
#Utilisons la fonction min_aic définie précédement pour obtenir le meilleur modèle
min_aic(Sinistres,c(0,12),c(1,1),c(0,12)) 
#les ordres donnes par l'algo sont 11,1,5 mais l'arima est meilleure avec les ordres 12,1,5
model_prime <-arima(SEst,order = c(12,1,5), method ='ML')

pred_prime<-predict(model_prime,12)

prime_pred <- pred_prime$pred
ggplot(prime_pred) + geom_line()+
  xlab("Année")+ 
  ggtitle("Prediction saison de 2021")
prime_pred_ts= ts(prime_pred,start=c(2021,1),frequency = 12)
ts.plot(SEst,prime_pred_ts,col=c("blue","orange"),main="Prédiction saison") 

#Tendance + Saison (predits)
#MUt=Mpt + St 
MUt = pred_tend +  prime_pred
ts.plot(logPrimes,MUt  ,col=c("blue","orange"),main="Prédiction MUt (LogPrimes)")


#Question 8
#valeur predites
prim_val<-ts(MUt,end =end(MUt),frequency = frequency(MUt))

#valeur predites supérieures 
prim_up<-ts(MUt + qnorm(0.975)*sd(PrimesDetend),end =end(MUt),frequency = frequency(MUt))

#valeur predites inférieures
prim_down<-ts(MUt-qnorm(0.975)*sd(PrimesDetend),end =end(MUt),frequency = frequency(MUt))
# Graphe
plot(logPrimes,
     type = 'l',col='blue',
     main="Prediction de logPrimes avec un intervalle de confiance de 95%"
     ,ylab="Primes",xlab="Année")
lines(prim_val,type = 'l',col='orange',lty=2,lwd=2)
points(prim_up, type = "l", col = 2, lty = 2)
points(prim_down, type = "l", col = 2, lty = 2)

#Question 9
#Initialement : Primes = Mt*St*Zt avec Mt la tendance, St la saisonnalité, Zt les residus
#Nous avons appliqué le log à la série chronologique Primes afin d'avoir un modèle additifs alors nous avons :
#logPrimes=lm.primes+SEst+PrimesDetend
#E(Primes)=E(exp(logPrimes))
#=E(exp(lm.primes+SEst+PrimesDetend))
#=E(ext(lm.primes))*E(exp(SEst))*E(exp(PrimesDetend)) par indépendance


#Question 10
#a
n=12
m=100
var_res=sd(PrimesDetend)
Mat <- matrix(rnorm(100*12,mean=0,var_res),n,m)

#b
Mat+matrix(MUt,n,m)
sum_mat=ts(Mat+matrix(MUt,n,m),start=c(2021,1),frequency = 12) 
ts.plot(sum_mat)

#c
#on utilise la fonction replicate pour simule 100 fois les 12 mois 
ZSMAT.prime = replicate(100,expr=rnorm(12,0,sd = sd(ResidusPrimes))  )
#ensuite on repete 100 fois les 12 valeurs predites 
rep.prime<- replicate(100,as.vector(MUt)) 

#Ajouter à chacune des trajectoires (à chaque colonne) l’espérance de primes
#et Appliquer exp à ces simulations
trajecto.prime.BB <- exp(ZSMAT.prime+rep.prime) 

#d
#Moyenner chaque ligne k = 1, . . . , 12 pour obtenir un estimateur de l’espérance de chaque valeur prédite.
montoCarlo.prime<- apply(trajecto.prime.BB,1,mean) 
ts.plot(montoCarlo.prime)

#e
#la fonction quantile appliquée à chaque ligne k = 1,.., 12 pour obtenir une estimation du quantile d’ordre 0,975 et 0,025
montoCarlo.prime.all = apply(trajecto.prime.BB, 1, function(yy) c(mean = mean(yy),quantile(yy,c(0.025, 0.975))))

#cette fonction permet de transformer la formule en série chronologique 
put_in_time<- function(data){ 
  return(ts(data,start=c(2021.000,1),frequency = 12))
}

MoCar.prime_m<- put_in_time(montoCarlo.prime.all["mean",])
MoCar.prime_l <- put_in_time(montoCarlo.prime.all["2.5%",])
MoCar.prime_u<- put_in_time(montoCarlo.prime.all["97.5%",])
ts.plot(Primes,MoCar.prime_m,MoCar.prime_l,MoCar.prime_u,col=1:4)


#Monte-Carlo

#Question 1

#Pour sinistres 
#tout d'abord on simule des valeurs de Sinistres par la methode de Monte-Carlo comme
#on a fait pour les Primes
#pour simuler 100 trajectoires des 12 valeurs futures de Sinistres.

#on utilise la fonction replicate pour simuler 100 fois les 12 mois 
ZSMAT.sini = replicate(100,expr=rnorm(12,0,sd = sqrt(var(SinistresDetend )))  )
#ensuit on repete 100 fois les 12 valures predite de sinistres
rep.sini<- replicate(100,as.vector(sini_pred)) 

#Ajouter à chacune des trajectoires (à chaque colonne) l’espérance de sinistres
trajecto.sini.BB <- ZSMAT.sini+rep.sini 

#Moyenner chaque ligne k = 1, . . . , 12 pour obtenir un estimateur de l’espérance de chaque valeur prédite.
montoCarlo.sini<- apply(trajecto.sini.BB,1,mean) 

#la fonction quantile appliquée à chaque ligne k = 1,.., 12 pour obtenir une estimation du quantile d’ordre 0,975 et 0,025
montoCarlo.sini.all = apply(trajecto.sini.BB, 1, function(yy) c(mean = mean(yy),quantile(yy,c(0.025, 0.975))))

MonCar.sini_mean<- put_in_time(montoCarlo.sini.all["mean",])
MonCar.sini_down <- put_in_time(montoCarlo.sini.all["2.5%",])
MonCar.sini_up<- put_in_time(montoCarlo.sini.all["97.5%",])
ts.plot(Sinistres,MonCar.sini_mean,MonCar.sini_up,MonCar.sini_down,col=1:4)


#Question 2
#Pour chacune des 100 trajectoires i = 1, . . . , 100, on calcule le Rvec=S/P

Rvec <- apply(trajecto.sini.BB,2,sum)  / apply(trajecto.prime.BB,2,sum) 

R2021.pred 		= mean(Rvec)
R2021.pred.up 	= quantile(Rvec,0.975)
R2021.pred.low 	= quantile(Rvec,0.025)

# Il faudrait augmenter les primes d'un facteur pour que les deux conditions soient vérifiées. Ramenons le quantile à 2,5% de R à 60% 
# R2021.pred.low/facteur = 0.6 
# facteur = R2021.pred.low/0.6 
facteur = R2021.pred.low/0.6 

# L'augmentation maximale serait de 1.499754 = 149% 
