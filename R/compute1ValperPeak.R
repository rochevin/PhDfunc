#' Title
#'
#' @param bed A GRanges object
#' @param wig A RLEList object
#' @param w Integer used to extend bed of both side : 2*w
#' @param seqlens A list wich contain the chromsize and their names
#' @param fun A character, both mean and sum are accepted
#'
#' @return A vector of integer for the coverage of each element in bed
#' @export
#'
#' @examples
compute1ValperPeak <- function( bed, wig, w = 20000, seqlens, fun = "sum" ){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    vec = NULL;
    for( i in 1:length( bed ) ){
        if( i %% 100 == 0 ){
            message( i, "/", length( bed ) );
        }
        bedi = bed[i, ];
        chr = as.character( seqnames( bedi ) );
        cov = wig[[chr]];
        center = start( bedi ) + (width(bedi)/2);
        stW = center - w;
        edW = center + w;

        stW[stW < 1] = 1;
        stW[stW > seqlens[chr] | stW > length( cov )] = min( length( cov), seqlens[chr] );
        edW[edW > seqlens[chr] | edW > length( cov )] = min( length( cov), seqlens[chr] );
        v = Views( cov, start = stW, end = edW );
        if( fun == "sum" ){
            vm = sum( v );
        }else if( fun == "mean" ){
            vm =  mean( v );
        }else{
            stop( "ERROR : unknown function fun = ", fun, " - fun must be 'sum' or 'mean'" );
        }
        vec = c( vec, vm );
    }
    return( vec );
}
