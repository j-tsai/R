####################################################################################
#Author: Jennifer Tsai 
#main: find roots of nonlinear function
#(Can be used to find local maxima, local minima, and MLE)
#	input:	f, function
#			[a,b], limits to search for brackting intervals
#			root_tol = difference in computed and true roots tolerance
# 			fxn_tol = tolerance of abs value of f computed at the root estimate
# 			itrlim = max number of iterations
#			h = step size
#
#	Example:
#		f = function(x) return(sin(x^2)*exp(-x))
#		main(f,a=-5,b=5,h=1.0e-8,root_tol= 1.0e-9, fxn_tol=1.0e-12, itrlim=100)
#
#####################################################################################
main = function(f,a,b,h=1.0e-8,root_tol= 1.0e-9, fxn_tol=1.0e-12, itrlim=100)
{
	stopifnot( is.function(f))
  	stopifnot( is.numeric(a) & length(a) == 1 )
	stopifnot( is.numeric(b) & length(b) == 1 )
	stopifnot( b>a )
  	stopifnot( is.numeric(itrlim) & length(itrlim) == 1 & (itrlim >0) )
	stopifnot( is.numeric(root_tol) & length(root_tol) == 1 & (root_tol>0))
	stopifnot( is.numeric(fxn_tol) & length(fxn_tol) == 1 & (fxn_tol>0))
	stopifnot( h >= 0)

	###############################################################################
	# bracket (Step 2): Search for bracketing intervals
	#			
	
	bracket <- function(f,a,b){
		
		x = seq(a,b,length.out=2^12) 
		bracket.a = 0
		bracket.b = 0
   	i = 1 		# initialize i
		while (i < (length(x)-1)){
        		if((f(x[i])*f(x[i+1]))<=0){ # root(s) exits within x[i] and x[i+1]
			
			# add the interval to bracketing interval vectors
			bracket.a = c(bracket.a, x[i]) 
			bracket.b = c(bracket.b, x[i+1]) 
			}
			i = i + 1
      	}
##browser()		
		if(length(bracket.a)<= 1 || length(bracket.b) <=1)
			stop("Error: No root in [a,b]")
		else{
			bracket.a = bracket.a[-1]
			bracket.b = bracket.b[-1]
			# Bracketing intervals that contain roots
			return(list(a=bracket.a,b=bracket.b))
 		}
			
	} #end bracket function

	a.list=bracket(f,a,b)$a
	b.list=bracket(f,a,b)$b

##traceback()

	##################################################
	# deriv (step 3): computes accurate numerical derivatives
	#	inputs:	f, function
	#			a, value for evaluating the derivative
	#			h, step size 
	# 
	deriv = function(f,a,h) return((f(a+h) - f(a))/h)

	###############################################################################
	# Step 1: Bisection algorithm to find a root given a braketing interval (a,b)
	# 	input: 	a,b list of bracketing intervals
	# 			
		
	bisect <- function(f, a.list, b.list, root_tol, fxn_tol, itrlim)
	{
		stopifnot( is.function(f))
		stopifnot( b.list > a.list )
		stopifnot( length(a.list) == length(b.list))
		
		#Initialize
		n = length(a.list)
		itr = rep(1,n)
		root=numeric(n)
		error=0.5*abs(b.list-a.list)
		fa = f(a.list)
		fb = f(b.list)

		for ( k in 1:n){
			if(sign(fa[k]) == sign(fb[k]))
				stop("Error: Root not in bracket")
			
			c = (a.list[k]+b.list[k])/2
			##################################################
			# Step 4: Newton's Method
				itr.new = 1
				repeat{
					stopifnot(deriv(f,c,h) != 0)
					temp = c
					c = c - f(c)/deriv(f,c,h)
					if( abs(temp-c) < 1.0e-12 || (c < a.list[k]) || (c > b.list[k]) || itr.new > itrlim) break
				}
			#
			##################################################
			
			while ( error[k] > root_tol & abs(fb[k] - fa[k]) > fxn_tol & itr[k] < 100){
				fc=f(c)
				
				if (fc ==0){
					error[k] = 0
					break #root found
				}
						
				else if (sign(fa[k]) != sign(fc)){ #root lies in [a,c]
					b.list[k] = c
					fb[k] = fc
				}
				else{ #root lines in [c,b]
					a.list[k] = c
					fa[k] = fc
				}
				c = (a.list[k]+b.list[k])/2
				error[k] = 0.5*abs(b.list[k]-a.list[k])
				itr[k] = itr[k] +1
			}
	
			root[k] = c
		}
		return(list(Root=root, Error = error, Iterations = itr))
		#cat(sprintf("Root = %.5f  Error = %e  Iterations =%.0f \n",root,error,itr))
	} # end bisection

	root=bisect(f, a.list, b.list, root_tol, fxn_tol, itrlim)$Root
	error=bisect(f, a.list, b.list, root_tol, fxn_tol, itrlim)$Error
	itr=bisect(f, a.list, b.list, root_tol, fxn_tol, itrlim)$Iterations
	return(list(Root=root, Error = error, Iterations = itr))
}



