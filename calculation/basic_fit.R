basic_fit  =  function( abc, y, grs ){
    if (is.null( dim( abc ) )) {
        abc  =  array( abc, dim = c( 4, length( abc ) / 4 ) )
    }
    if (is.null( dim( y ) )) {
        y  =  array( y, dim = c( length( y ), 1 ) )
    }
    if (is.null( dim( grs ) )) {
        grs  =  array( grs, dim = c( length( grs ), 1 ) )
    }

    a  =   abc[ 1, , drop = FALSE ]
    # b  =   abc[ 2, , drop = FALSE ]
    c  =   abc[ 2, , drop = FALSE ]
    sigE = abc[ 3, , drop = FALSE ]
    sigN = abc[ 4, , drop = FALSE ]

    d1  =  array( 1, dim = c( length( y ), 1 ) )
    d2  =  array( 1, dim = c( 1, length( a ) ) )

    v0  =  0.5 * colSums( ( y %*% d2 - grs %*% a )^2 / ( ( 1 + grs %*% c )^2 * (d1 %*% (sigE^2))  + d1 %*% (sigN^2) )
                      + log( ( 1 + grs %*% c )^2 * (d1 %*% (sigE^2))  + d1 %*% (sigN^2) ) )
    # v   =  0.5 * colSums( ( y %*% d2 )^2 )
# 
#     sel  =  s2 > 0
#     v[ sel ]  =  v0[ sel ]
    v0
}

