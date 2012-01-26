check.null.dim <-
function(vec0){

###########################################################################################################
##
## FUNCTION TO CHECK IF A DIMENSION OF A VECTOR IS NULL
## USED IN MAIN FUNCTION
## RETURNS VECTOR AS A MATRIX
##

	if(is.null(dim(vec0))){
		lov <- length(vec0)
		vec0 <- matrix(vec0,1,lov)
	}
return(vec0)
}

