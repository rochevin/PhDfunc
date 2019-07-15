#' computeProfile : Get the average/sum coverage of a RleList object for a list of genomic positions
#'
#' @param bed A GRanges object
#' @param wig A RLEList object
#' @param w Integer used to extend bed of both side : 2*w
#' @param span Integr used to determine the size of each bin of our 2*w interval
#' @param seqlens A list wich contain the chromsize and their names
#' @param method A character, both mean and sum are accepted
#'
#' @return A matrix of size length(bed)*(w/span)
#' @export
#'
#' @examples
computeProfile = function( bed, wig, w = 20000, span = 200, seqlens ,method="mean"){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    mat = NULL;
    for( i in 1:length( bed ) ){
        message( i, "/", length( bed ) );
        bedi = bed[i, ];
        chr = as.character( seqnames( bedi ) );
        cov = wig[[chr]];


        center = start( bedi ) + 4;
        stW = center - w;
        edW = center + w;

        if( span == 1 ){
            vm = as.numeric( Views( cov, start = stW, end = edW )[[1]] )
        }else{
            sts = seq( stW, edW - span + 1, span );
            eds = seq( stW + span - 1, edW, span );
            v = Views( cov, start = sts, end = eds );
            if(method =="sum"){
                vm =  sum( v );
            }else {
                vm =  mean( v );
            }

            vm[sts >= seqlens[chr] | sts > length( cov )] = 0;
            vm[sts < 1] = 0;
        }
        mat = rbind( mat, vm );
    }
    #rv = colMeans( mat );
    return( mat );
}
