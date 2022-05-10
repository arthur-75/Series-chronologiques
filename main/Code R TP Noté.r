
# Création données TP noté sinistre primes 
{
	# Simulation primes 
	{
		#Paramètres
		{
			NYears 		= 17
			Origin_P 	= 1
			kappa_s 	= 0.1 
			kappa_m 	= 2
			sigma1      = 0.01
		}
		Motif = c(2.5,1.2,1,1.5,1.1,1,1.5,1.1,1,1.5,1.1,1)
		Motif = Motif - mean(Motif)
		Motif = Motif*kappa_s
		
		plot(Motif,type="l")
		
		S = rep(Motif,NYears)
		S = ts(S,start = c(2005,1),frequency = 12)
		
		n 		= length(S)
	
		
		m 		= (1-(0.5-(1:n)*0.5/n)^2 )*kappa_m
		
		X 		= m + S + rnorm(n,0,sd=sigma1 )
		X 		= ts(X,start = start(S), frequency = frequency(S))
		plot(m)
		dev.new()
		par(mfrow = c(2,1))
		plot(X)
		 
		plot(exp(X +Origin_P ))
		
		Primes =  exp(Origin_P+ X)
	}
	
	# Sinistres 
	{
		#Les sinistres sont proportionnels à la tendance * Facteur d'inflation des couts vétérinaires 
		#On reprend donc le m vu précédemment 
		#Paramètres
		{
			kappa_Sin 	= 1 
			Origin_Sin 	= 6
			Kappa_veto  = 0.0025
			sigma2 		= 0.5
		}
		Sinistres 		= exp(m*kappa_Sin  + (1:n)*Kappa_veto )+ Origin_Sin+ rnorm(n,0,sd=sigma2)
		Sinistres 		= ts(Sinistres,start = start(S), frequency = frequency(S))
		dev.new()
		ts.plot(Primes,Sinistres,col=c("blue","red"))
		
		T 			=  time(Sinistres)
		out			=  lm(Sinistres~I(T) + I(T^2)+I(T^6))
		adjustSin 	= ts(out$fitted,start= start(Sinistres),frequency = frequency(Sinistres))
		
		dev.new()
		ts.plot(Primes,Sinistres,adjustSin,col=c("blue","red","green"))

	}
	
	#Evolution de R=S/T sur 12 mois glissant 
	{
		 h = 12
		 R = c()
		 for(i in (h):n)
		 {
			R[i-h] = sum(Sinistres[(i-h+1):i])/sum(Primes[(i-h+1):i])
		 }
		 
		  R = ts(R ,start=time(Sinistres)[h],frequency=frequency(Sinistres)) 
		  dev.new()
		  par(mfrow=c(2,1))
		  
		  ts.plot(Primes,Sinistres,adjustSin,col=c("blue","red","green"))
		  plot(R)
	}

	Primes0  	= Primes
	Sinistres0 	= Sinistres
	
	Primes      = ts(Primes[1:192],start = start(Primes),frequency = frequency(Primes))
	Sinistres     = ts(Sinistres[1:192],start = start(Primes),frequency = frequency(Primes))
	
}
 
 
# Analyse 
{
	rm(list=ls())
	setwd("C:/Users/tdumo/Dropbox/Documents/DOCUMENTS/ENSEIGNEMENT/NANTERRE/2020-2021/Series Chro/TP/SeriesChroTPNoté")
	load("Primes0.Rdata")
	load("Primes.Rdata")
	load("Sinistres0.Rdata")
	load("Sinistres.Rdata")
	source("SaisEst.r")
	# Modélisation 
	{
		# Primes 
		{
			n 			= length(Primes) 
			
			plot(Primes)
			# Amplitudes des oscillation qui augmente avec la tendance : transformation boxcox 
			# .... on choisit le log 
	
			boxPrimes 	= log(Primes)
			plot(boxPrimes)
			
			T 			= time(boxPrimes)
			out.tend.Primes = lm(boxPrimes~I(T) + I(T^2))
			
			boxPrimes.m = ts(out.tend.Primes$fitted,start = start(boxPrimes),frequency = frequency(boxPrimes))
			ts.plot(boxPrimes,boxPrimes.m ,col=c("black","blue"))
			
			
			boxPrimes.detend = boxPrimes -boxPrimes.m
			
			dev.new()
			plot(boxPrimes.detend )
			
			List.sais 		= SaisEst(boxPrimes.detend,12)
			
			boxPrimes.motif	=  List.sais$motif
			boxPrimes.S		=  List.sais$serie
			
			ts.plot(boxPrimes.detend ,boxPrimes.S	,col=c("black","red"))
			
			boxPrimes.Res 	= boxPrimes.detend - boxPrimes.S
			dev.new()
			par(mfrow = c(2,1))
			plot(boxPrimes.Res)
			acf(boxPrimes.Res,lag=50)
			
			Box.test(boxPrimes.Res, lag= 50,type = "Ljung-Box") 
			
			# Residus indépendants 

			# Prédiction de boxPrimes 
			# Model : boxPrimes = m + S + Z 
			
			#Pour tout t \in[t0 + 1 , .... t0 + q] 
			# On prédit suivant le modèle vu dans l'énoncé : 
			
			time(boxPrimes)
			length(time(boxPrimes))
			# [1] 192
			time(boxPrimes)[(192-11):192]
			time(boxPrimes)[(192-11):192] + 1 
			
			newT 	= time(boxPrimes)[(192-11):192] + 1 


			m.pred 	=   predict(out.tend.Primes , newdata=data.frame(T=newT))
		 
			S.pred 	= boxPrimes.motif
			
		    mu.Pred = m.pred + S.pred 
			 
			sigma2Est 	= var(boxPrimes.Res)
			
			Primes.pred 	= exp(sigma2Est/2) *exp( mu.Pred)
			Primes.pred.var =  (exp(sigma2Est*2) - exp(sigma2Est))* exp(2*mu.Pred)
			
			Primes.pred        = ts(Primes.pred,start=c(2021,1),frequency = 12)
			Primes.pred.up     = ts(Primes.pred + 1.96*sqrt(Primes.pred.var),start=c(2021,1),frequency = 12)
			Primes.pred.low    = ts(Primes.pred - 1.96*sqrt(Primes.pred.var),start=c(2021,1),frequency = 12)
		
			SumPrimes.pred 	= exp(sigma2Est/2) *sum(exp( mu.Pred))
			pred.P.var    = (exp(sigma2Est*2) - exp(sigma2Est))*sum(exp( 2*mu.Pred))
			
			
			dev.new()
			ts.plot(Primes,Primes.pred,Primes.pred.up,Primes.pred.low,xlim = c(2005,2022),col=c("black","red","orange","orange"),lty=c(1,1,2,2))
			
			
			lines(c(2021,2022),c(SumPrimes.pred/12,SumPrimes.pred/12),col="blue")
			lines(c(2021,2022),c((SumPrimes.pred +1.96*sqrt(pred.P.var))/12,(SumPrimes.pred +1.96*sqrt(pred.P.var))/12),col="blue",lty=2)
			lines(c(2021,2022),c((SumPrimes.pred -1.96*sqrt(pred.P.var))/12,(SumPrimes.pred -1.96*sqrt(pred.P.var))/12),col="blue",lty=2)
		
			sum(Primes0[193:length(Primes0)])
			c((SumPrimes.pred -1.96*sqrt(pred.P.var)),SumPrimes.pred,(SumPrimes.pred +1.96*sqrt(pred.P.var)) )
			
		}
	
		#Sinistres
		{
			plot(Sinistres)
			T = time(Sinistres)
			out.tend.Sinistres = lm(Sinistres~I(T) + I(T^2)+I(T^6) )
			Sinistres.m = ts(out.tend.Sinistres$fitted,start = start(Sinistres),frequency = frequency(Sinistres))
			ts.plot(Sinistres,Sinistres.m ,col=c("black","blue"))
		
			Sinistres.detend = Sinistres - Sinistres.m 
			plot(Sinistres.detend )
			
			dev.new()
			par(mfrow = c(2,1))
			plot(Sinistres.detend)
			acf(Sinistres.detend,lag=50)
			
			Box.test(Sinistres.detend, lag= 50,type = "Ljung-Box") 
			
			# Model 	 			: Sinistres = mt  + Zt 
			# Prédiction 			: somme(mt)
			#Variance d'estimation  : 12*var(ZT0)
			
			VarEst = var(Sinistres.detend) 
			
			newT 	= time(Sinistres)[(192-11):192] + 1 

			m.pred 			= predict(out.tend.Sinistres, newdata=data.frame(T=newT))
		 
			Sinistres.pred 	   = m.pred 	
			Sinistres.pred.var = VarEst
			
			Sinistres.pred        = ts(Sinistres.pred,start=c(2021,1),frequency = 12)
			Sinistres.pred.up     = ts(Sinistres.pred + 1.96*sqrt(Sinistres.pred.var),start=c(2021,1),frequency = 12)
			Sinistres.pred.low    = ts(Sinistres.pred - 1.96*sqrt(Sinistres.pred.var),start=c(2021,1),frequency = 12)
		
			SumSinistres.pred 	=  sum(m.pred)
			pred.Sin.var        = 12*VarEst
		 
			dev.new()
			ts.plot(Sinistres,Sinistres.pred,Sinistres.pred.up,Sinistres.pred.low,xlim = c(2005,2022),col=c("black","red","orange","orange"),lty=c(1,1,2,2))
			
			
			lines(c(2021,2022),c(SumSinistres.pred/12,SumSinistres.pred/12),col="blue")
			lines(c(2021,2022),c((SumSinistres.pred +1.96*sqrt(pred.Sin.var))/12,(SumSinistres.pred +1.96*sqrt(pred.Sin.var))/12),col="blue",lty=2)
			lines(c(2021,2022),c((SumSinistres.pred-1.96*sqrt(pred.Sin.var))/12,(SumSinistres.pred -1.96*sqrt(pred.Sin.var))/12),col="blue",lty=2)
		
			sum(Sinistres0[193:length(Sinistres0)])
			c((SumSinistres.pred -1.96*sqrt(pred.Sin.var)),SumSinistres.pred,(SumSinistres.pred +1.96*sqrt(pred.Sin.var)) )
		}
		
		#R = S/T 
		{
			Rplus  = (SumSinistres.pred +1.96*sqrt(pred.Sin.var) )/(SumPrimes.pred -1.96*sqrt(pred.P.var))
			Rmoins = (SumSinistres.pred -1.96*sqrt(pred.Sin.var) )/(SumPrimes.pred +1.96*sqrt(pred.P.var))
			
			c(Rmoins,Rplus)
			sum(Sinistres0[193:length(Sinistres0)])/sum(Primes0[193:length(Primes0)])
		}
	}	
	

}


# Résolution TP 
{
	rm(list=ls())
	setwd("C:/Users/tdumo/Dropbox/Documents/DOCUMENTS/ENSEIGNEMENT/NANTERRE/2020-2021/Series Chro/TP/SeriesChroTPNoté")
	load("Primes0.Rdata")
	load("Primes.Rdata")
	load("Sinistres0.Rdata")
	load("Sinistres.Rdata")
	source("SaisEst.r")
	source("MM.r")


	#Sinistres
	{
		#1 Représentation 
			plot(Sinistres)
			n = length(Sinistres)
			
			# Le modèle additif semble crédible : une tendance + un bruit stationnaire : m + Z
			# Estimons la tendance et, et jugeons de la stationnarité des résidus 
		
		#2 Ajustement
			T = time(Sinistres)
			out.tend.Sinistres = lm(Sinistres~I(T) + I(T^2)+I(T^6) )
			
			# Récupération de l'estimation de m : 
			Sinistres.m = ts(out.tend.Sinistres$fitted,start = start(Sinistres),frequency = frequency(Sinistres))
			ts.plot(Sinistres,Sinistres.m ,col=c("black","blue"))
	
			# On retranche la tendance
			Sinistres.residus = Sinistres - Sinistres.m 
			plot(Sinistres.residus )
			residus.MM = MM(Sinistres.residus,q=12)
			ts.plot(Sinistres.residus,residus.MM,col=c("black","red"))
			# oscille autour de zero : rien qui ne vienne contredire la stationnarité 
			
			# Regardons le carré des résidus 
			residus.MM2 = MM(Sinistres.residus^2,q=12)
			ts.plot(Sinistres.residus^2,residus.MM2,col=c("black","red"))
			# Là non plus, rien qui ne vienne contredire la stationnarité 
			
			# Comparons l'ACF de la première partie de la série avec la deuxième partie de la série 
			residus.1 = ts(Sinistres.residus[1:(n/2)] 	,start= T[1]		, frequency = frequency(Sinistres))
			residus.2 = ts(Sinistres.residus[(n/2+1):n]	,start= T[(n/2+1)]	, frequency = frequency(Sinistres))
			
			dev.new() 
			par(mfrow = c(1,2))
			acf(residus.1,lag=30)
			acf(residus.2,lag=30)
			# Là non plus, rien qui ne vienne contredire la stationnarité 
		
			# colle au modèle Sinistres = mt  + Zt où Zt est stationnaire
		
		#3 Test d'indépendance des résidus
			dev.new()
			par(mfrow = c(2,1))
			plot(Sinistres.residus)
			acf(Sinistres.residus,lag=50)
			
			Box.test(Sinistres.residus, lag= 50,type = "Ljung-Box") 
			#On ne rejette pas l'hypothèse d'indépendance 
		
		#4 Prédiction des 12 prochaines valeurs de Sinistres
			# Dans le modèle  Sinistres = mt  + Zt 
			# Les simistres sont indépendants 
			# la meilleur prédiction linéaire de Sinistres est mt
			# calculons mt pour les 12 prochaines valeurs de T 
			T
			# Les 12 dernières valeurs de T sont : 
			#     2020.000 2020.083 2020.167 2020.250 2020.333 2020.417 2020.500 2020.583  2020.667 2020.750 2020.833 2020.917
			#Les valeurs à prédire correspondent à 
			# t = 2021.000 2021.083 2021.167 2021.250 2021.333 2021.417 2021.500 2021.583  2021.667 2021.750 2021.833 2021.917
			
			newT 					= c(2021.000, 2021.083, 2021.167, 2021.250, 2021.333, 2021.417, 2021.500, 2021.583,  2021.667, 2021.750, 2021.833, 2021.917) 
			m.pred 	    			= predict(out.tend.Sinistres, newdata=data.frame(T=newT))
			Sinistres.pred 			= ts(m.pred,start=c(2021.000,1),frequency = 12) 
			
			ts.plot(Sinistres,Sinistres.pred,col=c("black","red"))
		
		#5 Interval de variation pour la prédiction des Sinistres
			Sinistres.pred.var 	  = var(Sinistres.residus) 
			Sinistres.pred.up     = Sinistres.pred + 1.96*sqrt(Sinistres.pred.var) 
			Sinistres.pred.low    = Sinistres.pred - 1.96*sqrt(Sinistres.pred.var) 
			ts.plot(Sinistres,Sinistres.pred,Sinistres.pred.up ,Sinistres.pred.low ,col=c("black","red","orange","orange"),lty=c(1,1,2,2))
	}
	
	# Primes 
	{
		#1 Représentation 
			n 			= length(Primes) 
			plot(Primes)
			# Amplitudes des oscillation qui augmente avec la tendance : transformation boxcox 
			#L'énoncé nous dit de choisir le log 
		
		#2 Transformation
			logPrimes 	= log(Primes)
			plot(logPrimes)
			
		#3 Ajustement polynomiale de la tendance
			T 				= time(logPrimes)
			out.tend.Primes = lm(logPrimes~I(T) + I(T^2))
				
			logPrimes.m = ts(out.tend.Primes$fitted,start = start(logPrimes),frequency = frequency(logPrimes))
			ts.plot(logPrimes,logPrimes.m ,col=c("black","blue"))
		
		#4 Ajustement de la saisonnalité
			logPrimes.detend = logPrimes -logPrimes.m
		
			dev.new()
			plot(logPrimes.detend )
			
			List.sais 		= SaisEst(logPrimes.detend,12)
			logPrimes.motif	=  List.sais$motif
			logPrimes.S		=  List.sais$serie
		
			ts.plot(logPrimes.detend ,logPrimes.S	,col=c("black","red"))
		
		#5 Stationnarité des résidus 
			logPrimes.residus 	= logPrimes.detend - logPrimes.S
			dev.new()
			par(mfrow = c(2,1))
			plot(logPrimes.residus)
			
			residus.MM = MM(logPrimes.residus,q=12)
			ts.plot(logPrimes.residus,residus.MM,col=c("black","red"))
			# oscille autour de zero : rien qui ne vienne contredire la stationnarité 
			
			# Regardons le carré des résidus 
			residus.MM2 = MM(logPrimes.residus^2,q=12)
			ts.plot(logPrimes.residus^2,residus.MM2,col=c("black","red"))
			# Là non plus, rien qui ne vienne contredire la stationnarité 
			
			# Comparons l'ACF de la première partie de la série avec la deuxième partie de la série 
			residus.1 = ts(logPrimes.residus[1:(n/2)] 	,start= T[1]		, frequency = frequency(logPrimes))
			residus.2 = ts(logPrimes.residus[(n/2+1):n]	,start= T[(n/2+1)]	, frequency = frequency(logPrimes))
			
			dev.new() 
			par(mfrow = c(1,2))
			acf(residus.1,lag=30)
			acf(residus.2,lag=30)
			# Là non plus, rien qui ne vienne contredire la stationnarité 
		
			# colle au modèle logPrimes = mt  + St + Zt où Zt est stationnaire
		
		#6 Test d'indépendance des résidus
			dev.new()
			par(mfrow = c(2,1))
			plot(logPrimes.residus)
			acf(logPrimes.residus,lag=50)
			
			Box.test(logPrimes.residus, lag= 50,type = "Ljung-Box") 
			#On ne rejette pas l'hypothèse d'indépendance 
		 
		#7 	Prédiction des 12 prochaines valeurs de logPrimes 
			#7.(a) Prédiction des 12 prochaines valeurs pour mt
				newT 					= c(2021.000, 2021.083, 2021.167, 2021.250, 2021.333, 2021.417, 2021.500, 2021.583,  2021.667, 2021.750, 2021.833, 2021.917 ) 
				m.pred 	    			= predict(out.tend.Primes, newdata=data.frame(T=newT))
				logPrimes.m.pred 	    = ts(m.pred,start=c(2021.000,1),frequency = 12) 
				
				ts.plot(logPrimes,logPrimes.m.pred ,col=c("black","red"))
			
			#7.(b) Prédiction des 12 prochaines valeurs pour St
				# Très facile : St est 12-périodique le motif est le suivant : 
				ts.plot(logPrimes.motif) 
				
				# Il faut juste savoir, dans la série d'origine, à quel moment du motif on s'arrête. 
				dev.new() 
				par(mfrow=c(3,1))
				ts.plot(logPrimes.motif) 
				ts.plot(logPrimes.S)
				ts.plot(logPrimes.S[(n-11):n])
				
				# La saisonnalité s'est arrêté à la fin du motif 
				# Pour prédire, on commence au début du motif 
				logPrimes.S.pred 	    = ts(logPrimes.motif,start=c(2021.000,1),frequency = 12) 
				ts.plot(logPrimes.S,logPrimes.S.pred ,col=c("black","red"))
			
			#7.(c) Prédiction des 12 prochaines valeurs pour logPrimes = mt + St =MUt
				MUt 				 = logPrimes.m.pred +  logPrimes.S.pred 
				logPrimes.pred       = MUt
				dev.new() 
				ts.plot(logPrimes,logPrimes.pred  ,col=c("black","red"))
		
		#8 Interval de variation pour la prédiction de logPrimes
				logPrimes.pred.var 	  = var(logPrimes.residus) 
				logPrimes.pred.up     = logPrimes.pred + 1.96*sqrt(logPrimes.pred.var) 
				logPrimes.pred.low    = logPrimes.pred - 1.96*sqrt(logPrimes.pred.var) 
				ts.plot(logPrimes,logPrimes.pred,logPrimes.pred.up ,logPrimes.pred.low ,col=c("black","red","orange","orange"),lty=c(1,1,2,2))
				
		#9 Calcul de l'espérance de Primes 
		    #  log(Primes) = logPrimes = mt + St + Zt
			#  Primes = exp(mt + St + Zt)
			#  Primes_t = exp(mt + St)*exp(Zt)
			#  Comme les Zt sont supposés indépendant, les primes_t le sont aussi. 
			# Le meilleur prédicteur linéaire E( Primes_t | Primes_(t-1),Primes_(t-2),Primes_(t-3),...) =  E( Primes_t )
			# Calcul de E(Primes_t) = exp(mt + St)*E(exp(Zt))
			# E(exp(Zt))  = exp(sigma^2/2) 
			# Ainsi, le meilleur prédicteur linéaire n'est pas exp(mt + St)
			# mais  exp(mt + St)*exp(sigma^2/2)
			
			# Cette "subtilité" concerne aussi la construction des intervalles de confiance. Même si les Z sont gaussiens, la loi de e^Zt n'est pas Gaussienne ! 
			# Attention à la construction des intervalles de confiance dans ce cas ! 
		
		#10 Méthode de Monte-Carlo pour la construction des intervalles de confiance : simulations de trajectoires 
			# Dans le modèle  Primes = exp(mt + St + Zt) ce qui est aléatoire : c'est Zt : Simulons plusieurs trajectoires des 12 valeurs futures de Zt
			# Sahant qu'ils sont i.i.d. normales centrées de variance : 			
			logPrimes.pred.var 	 
				
			#10.1 Simuler et stocker 100 trajectoires de Z de longueur 12 dans une matrice 12x100 
				#Alimentons une matrice N = 12 lignes et 100 colonnes 		
				ZPmat = matrix(rnorm(100*12,mean=0,sd = sqrt(logPrimes.pred.var )),12,100 )
				matplot(ZPmat,type="l")
			
			#10.2 Ajouter à chacune des trajectoires l'espérance de logPrimes. Représenter logPrimes avec les 100 trajectoires futures
				# Pour chacune des 100 trajectoires on va calculer les Primes correxpondantes 
				MUmat = matrix(MUt,12,100 )
				MUmat
				matplot(MUmat,type="l")
				
				MClogPrimes = ts( MUmat + ZPmat,start=c(2021.000,1),frequency = 12) 
				ts.plot(MClogPrimes)
				ts.plot(logPrimes,MClogPrimes,col=c("black",rep("red",100)))
				
			#10.3 Appliquer exp à ces simulations pour créer des simulations de trajectoires pour la série d'origine Primes. Représenter
				MCPrimes = exp(MClogPrimes)
				ts.plot(Primes,MCPrimes,col=c("black",rep("red",100)))
				
			#10.4 utilisation de ces simulations pour approcher les intervalles de confiance 
				Primes.pred 	= rep(0,12) 
				Primes.pred.up 	= rep(0,12) 
				Primes.pred.low = rep(0,12) 
				
				for(k in 1:12)
				{
					Primes.pred[k]		= mean(MCPrimes[k,]) 
					Primes.pred.up[k]	= quantile(MCPrimes[k,],0.975) 
					Primes.pred.low[k]	= quantile(MCPrimes[k,],0.025) 
				}
				Primes.pred 	= ts( Primes.pred,start=c(2021.000,1),frequency = 12) 
				Primes.pred.up 	= ts( Primes.pred.up,start=c(2021.000,1),frequency = 12) 
				Primes.pred.low = ts( Primes.pred.low,start=c(2021.000,1),frequency = 12) 
				
				ts.plot(Primes,Primes.pred,Primes.pred.up,Primes.pred.low,col=c("black","red","orange","orange"))
			
	}
	
	#R=S/T
	{
		#1 Simulation Monte-Carlo des valeurs de Sinistres 
			Sinistres.pred.var 	 
			ZSmat = matrix(rnorm(100*12,mean=0,sd = sqrt(Sinistres.pred.var )),12,100 )
			matplot(ZSmat,type="l")
		
			m.mat = matrix(m.pred,12,100 )
			m.mat 
			matplot(m.mat ,type="l")
			
			MCSinistres = ts( m.mat + ZSmat,start=c(2021.000,1),frequency = 12) 
			ts.plot(MCSinistres)
			ts.plot(Sinistres,MCSinistres,col=c("black",rep("red",100)))
		
		#2 Pour chaque trajectoire : calcul de R = S/T 
			
			Rvec = rep(0,100) 
			
			for(i in 1:100)
			{
				Rvec[i] = sum(MCSinistres[,i])/sum(MCPrimes[,i])
			}
			Rvec
			
			R2021.pred 		= mean(Rvec)
			R2021.pred.up 	= quantile(Rvec,0.975)
			R2021.pred.low 	= quantile(Rvec,0.025)
			
			c(R2021.pred.low,R2021.pred,R2021.pred.up )
			
			# Il faudrait augmenter les primes d'un facteur pour que les deux conditions soient vérifiées. Ramenons le quantile à 2,5% de R à 60% 
			# R2021.pred.low/facteur = 0.6 
			# facteur = R2021.pred.low/0.6 
			facteur = R2021.pred.low/0.6 
			
			# L'augmentation maximale  serait de 1.488943 = 148%
			
		
	}


}