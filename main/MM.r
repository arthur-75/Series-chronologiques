MM = function(x,q)
	{
		# x est un argument de la fonction 
		# lorsqu'on rentre dans la fonction on ne connait pas la taille de x 
		# il faut la définir à l'intérieur de la fonction 
		
		n = length(x) 
		
		m = rep(0,n-2*q)
		
		for(t in (q+1):(n-q))
		{
			# décalage de temps par rapport à la définition
			# m[1] correspond au temps t = q+1
			# le pas i correspond au temps t = 
			i = t-q
			m[i] = (1/(2*q+1)) * sum(x[(t-q):(t+q)])
		}
		
		# on converti m en objet de type ts et on précise que m démarre à q+1 
		
		m = ts(m,start = time(x)[q+1],frequency = frequency(x))
		
		return(m)				
	}
