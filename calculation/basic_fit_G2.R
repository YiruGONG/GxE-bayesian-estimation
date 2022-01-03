basic_fit_G2  =  function( abc, y, grs ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 5, length( abc ) / 5 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( grs ) )) {
        grs  =  array( grs, dim = c( length( grs ), 1 ) )
    }

    a1 =  abc[ 1, ]
    a2 =  abc[ 2, ]
    c  =  abc[ 3, ]
    sigE  =  abc[ 4, ]
    sigN  =  abc[ 5, ]

    d1  =  array( 1, dim = c( length( y ), 1 ) )
    d2  =  array( 1, dim = c( 1, length( a ) ) )
    g2  =  grs^2

    v0  =  0.5 * colSums( ( y %*% d2 -grs %*% a1 - (g2-1) %*% a2 )^2
                          / ( ( 1 + grs %*% c )^2 * (d1%*%(sigE^2) ) + d1 %*% (sigN^2) )
                          + log( ( 1 + grs %*% c )^2 * (d1%*%(sigE^2) ) + d1 %*% (sigN^2) ) ) 
    # v  =  0.5 * colSums( ( y %*% d2 )^2 )
    # 
    # sel  =  s2 > 0
    # v[ sel ]  =  v0[ sel ]
    # v
    v0
}
