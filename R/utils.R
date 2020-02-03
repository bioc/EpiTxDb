
insert_data_into_table <- function(conn, table, data){
    stopifnot(is.list(data))
    placeholders <- paste(rep.int("?", length(data)), collapse = ",")
    SQL <- sprintf("INSERT INTO %s VALUES (%s)", table, placeholders)
    ## dbExecute() emits annoying warnings if 'params' is a named list or if
    ## some of its list elements are factors.
    params <- unname(as.list(data))
    params <- lapply(params,
                     function(x) if (is.factor(x)) as.character(x) else x)
    dbExecute(conn, SQL, params = params)
}

duplicatedIntegerQuads <- S4Vectors:::duplicatedIntegerQuads
orderIntegerQuads <- S4Vectors:::orderIntegerQuads
matchIntegerQuads <- S4Vectors:::matchIntegerQuads

makeFeatureIds <- function(name = NULL, type, start, end,
                           same.id.for.dups = FALSE){
    a <- type
    b <- name
    c <- start
    d <- end
    fac <- c == d
    d[fac] <- ""
    if (!same.id.for.dups) {
        oo <- orderIntegerQuads(a, b, c, d)
        ans <- integer(length(oo))
        ans[oo] <- seq_len(length(oo))
        return(ans)
    }
    ## There should be a better way to do this...
    is_not_dup <- !duplicatedIntegerQuads(a, b, c ,d)
    ua <- a[is_not_dup]
    ub <- b[is_not_dup]
    uc <- c[is_not_dup]
    ud <- d[is_not_dup]
    oo <- orderIntegerQuads(ua, ub, uc, ud)
    ua <- ua[oo]
    ub <- ub[oo]
    uc <- uc[oo]
    ud <- ud[oo]
    matchIntegerQuads(a, b, c, d, ua, ub, uc, ud)
}
