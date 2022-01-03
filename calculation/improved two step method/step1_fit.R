step1_fit  =  function( abc, y, grs ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 2, length( abc ) / 2 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( grs ) )) {
        grs  =  array( grs, dim = c( length( grs ), 1 ) )
    }

    a1  =  abc[ 1, , drop = FALSE ]
    a2  =  abc[ 2, , drop = FALSE ]

    d1  =  array( 1, dim = c( length( y ), 1 ) )
    d2  =  array( 1, dim = c( 1, length( a1 ) ) )
    g2  = grs^2

    v = 0.5 * colSums( (y %*% d2 - grs %*% a1 - (g2-1) %*% a2)^2 )
    v
}

