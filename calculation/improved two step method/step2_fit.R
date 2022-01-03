step2_fit  =  function( abc, y, grs, a ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 3, length( abc ) / 3 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( grs ) )) {
        grs  =  array( grs, dim = c( length( grs ), 1 ) )
    }

    a1  =  a[ 1 ]
    a2  =  a[ 2 ]
    b   =  abc[ 1, ]
    c   =  abc[ 2, ]
    sig =  abc[ 3, ]

    d1  =  array( 1, dim = c( length( y ), 1 ) )
    d2  =  array( 1, dim = c( 1, length( b ) ) )
    g2  =  grs^2
    
    yr  =  y - a1 * grs - a2 * (g2 - 1)
    
    v   =  0.5 * colSums( ( (yr^2) %*% d2 - (d1 %*% b + grs %*% c)^2 - d1 %*% (sig^2) )^2 )
    v
}
