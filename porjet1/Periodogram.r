Periodogram = function(X,range.type = 1){
	#retourne une liste à deux éléments : les omega et les In(omega)
	n = length(X)
	T = 1:n
	
	if(range.type ==1)
	{
		k_min = floor((n-1)/2) 
		k_max = floor(n/2) 
		
		OMEGA = 2*pi*((-k_min):k_max) /n
	}else if(range.type ==2)
	{
		OMEGA = 2*pi*(0:(n-1)) /n
	}
	
	In = rep(0,length(OMEGA))
	
	for(i in 1:length(OMEGA))
	{
		In[i] = 1/n * (Mod(sum(X*complex(real =  cos(-OMEGA[i]*T), imaginary =sin(-OMEGA[i]*T)))))^2
	}
	 
	return( list(OMEGA=OMEGA,In = In))
}
