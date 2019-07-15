#' Title
#'
#' @param bed A GRanges object
#' @param wig A RLEList object
#' @param w Integer used to extend bed of both side : 2*w
#' @param span Integr used to determine the size of each bin of our 2*w interval
#' @param seqlens A list wich contain the chromsize and their names
#' @param fun A character, both mean and sum are accepted
#' @param get either "one" or "all". If get == "one", perform fun, if not, return all values
#'
#' @return If get == "one", return a vector of value for each position w/span, if get =="all", return a matrix with each values.
#' @export
#'
#' @examples
computeProfileGenes = function( bed, wig, w = 3000, span = 200, seqlens , fun="sum",get="one"){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    mat = NULL;
    for( i in 1:length( bed ) ){
        if( i %% 100 == 0 ){
            message( i, "/", length( bed ) );
        }
        bedi = bed[i, ];
        len = end(bedi) - start(bedi);
        if( len < 100 ){
            message( "Warning : Gene too small - n° ", i, " - ", bedi$name );
        }else{
            chr = as.character( seqnames( bedi ) );
            if(!chr %in% names(wig)){
                next;
            }
            cov = wig[[chr]];
            stW_st = start(bedi) - w;
            edW_st = start(bedi);
            stW_ed = end(bedi);
            edW_ed = end(bedi) + w;
            sts_st = seq( stW_st, edW_st - span + 1, span );
            eds_st = seq( stW_st + span - 1, edW_st, span );
            sts_ed = seq( stW_ed, edW_ed - span + 1, span );
            eds_ed = seq( stW_ed + span - 1, edW_ed, span );
            div = len / 100;
            ent = trunc( div );
            dec = div - ent;
            if( dec > 0.75 ){
                span_c = ent + 1;
            }else{
                span_c = ent;
            }
            sts_c = seq( start(bedi), end(bedi) - 1, span_c );
            eds_c = seq( start(bedi) + span_c - 1, end(bedi) + span_c - 1, span_c );
            if( length( sts_c ) > 100 ){
                sts_c = sts_c[-c( 101:length( sts_c ) )];
            }
            if( length( eds_c ) > 100 ){
                eds_c = eds_c[-c( 101:length( eds_c ) )];
            }
            if( length( sts_c ) < 100 | length( eds_c ) < 100 ){
                message( "Warning : Less than 100 regions - n° ", i, " - ", bedi$name);
            }else{
                eds_c[100] = end(bedi);
                v_st = Views( cov, start = sts_st, end = eds_st );
                vm_st =  switch(fun,
                                "sum"=sum(v_st),
                                "mean"=mean( v_st ));
                vm_st[sts_st >= seqlens[chr] | sts_st > length( cov )] = 0;
                vm_st[sts_st < 1] = 0;
                v_ed = Views( cov, start = sts_ed, end = eds_ed );
                vm_ed =  switch(fun,
                                "sum"=sum(v_ed),
                                "mean"=mean( v_ed ));
                vm_ed[sts_ed >= seqlens[chr] | sts_ed > length( cov )] = 0;
                vm_ed[sts_ed < 1] = 0;
                v_c = Views( cov, start = sts_c, end = eds_c );
                vm_c =  switch(fun,
                               "sum"=sum(v_c),
                               "mean"=mean( v_c ));
                vm_c[sts_c >= seqlens[chr] | sts_c > length( cov )] = 0;
                vm_c[sts_c < 1] = 0;
                if( as.character( strand(bedi) ) == "+" ){
                    vec = c( vm_st, vm_c, vm_ed );
                }else{
                    vec = c( rev( vm_ed ), rev( vm_c ), rev( vm_st ) );
                }
                if( length( vec ) == ( 100 + ( 2 * ( w / span ) ) ) ){
                    mat = rbind( mat, vec );
                }else{
                    stop( "ERROR : profiles of different size" );
                }
            }
        }
    }
    if(get == "one"){
        rv = switch(fun,
                    "sum"=colSums( mat ),
                    "mean"=colMeans( mat ));
    }else{
        rv = mat
    }
    return( rv );
}
