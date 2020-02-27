#' @include AllGenerics.R
NULL

#' @name positionSequence
#'
#' @title Generate integer sequences from position information of \code{Ranges}
#'
#' @description \code{positionSequence} generates sequences of integer values
#' along the range information of \code{x}. This can be used for navigating
#' specific positions on a range information.
#'
#' @param x a \code{Ranges} object, like a
#'   \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#'   \code{\link[IRanges:IRanges-class]{IRanges}}, or a \code{RangesList}
#'   object, like a \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
#'   or \code{\link[IRanges:IRangesList-class]{IRangesList}}
#' @param order \code{TRUE} or \code{FALSE}: Should the position be ordered?
#'   (default: \code{order = FALSE})
#' @param decreasing \code{TRUE} or \code{FALSE}: If \code{order = TRUE} Should
#'   the position be ordered in a decreasing order? (default: \code{order =
#'   FALSE})
#'
#' @return a \code{integer} vector if x is a
#'   \code{\link[GenomicRanges:GRanges-class]{GRanges}} object and a
#'   \code{IntegerList} if x is a
#'   \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
#'
#' @examples
#' library(GenomicRanges)
#' # Returns an integer vector
#' gr <- GRanges("chr1:1-5")
#' positionSequence(gr)
#' # returns an IntegerList
#' grl <- GRangesList("1" = gr,"2" = gr,"3" = gr) # must be named
#' positionSequence(grl)
NULL

# per element of GRangesList unique
.get_strand_u_GRangesList <- function(grl){
    strand_u <- unique(IRanges::CharacterList(strand(grl)))
    ans <- unlist(strand_u)
    if(length(strand_u) != length(ans)){
        stop("Non unqiue strands per GRangesList element.")
    }
    ans
}

# Vectorize version of seq specific for start/ends from a RangesList
.seqs_rl_strand <- function(rl, order = FALSE, decreasing = FALSE){
    strand_u <- .get_strand_u_GRangesList(rl)
    strand_minus <- strand_u == "-"
    ansP <- .seqs_rl_by(rl[!strand_minus])
    ansM <- .seqs_rl_by(rl[strand_minus], by = -1L)
    if(order){
        ansM <- ansM[IRanges::IntegerList(lapply(ansM, order,
                                                 decreasing = decreasing))]
    }
    ans <- c(ansP, ansM)
    ans[match(names(rl),names(ans))]
}

# Vectorize version of seq specific for start/ends from a RangesList with a by
# option
.seqs_rl_by <- function(rl, by = 1L){
    starts <- unlist(start(rl))
    ends <- unlist(end(rl))
    .seqs_l_by(starts, ends, by)
}

#' @importFrom IRanges PartitioningByWidth PartitioningByEnd
#' @importClassesFrom IRanges PartitioningByWidth PartitioningByEnd
# Vectorize version of seq using to input lists
.seqs_l_by <- function(from, to, by = 1L){
    if(is.null(names(from)) || is.null(names(to))){
        stop("Inputs must be named.")
    }
    if(length(from) != length(to)){
        stop("Inputs must have the same length.")
    }
    if(by == 0L){
        stop("by cannot be zero.")
    }
    if(any(names(from) != names(to))){
        stop("Unmatched names.")
    }
    if(by < 0L){ # switch from to around if negative
        tmp <- to
        to <- from
        from <- tmp
        rm(tmp)
    }
    ans <- Map(
        function(f,t){
            seq.int(f,t,by)
        },
        from,
        to)
    ans <- IRanges::IntegerList(ans)
    partitioning <- IRanges::PartitioningByEnd(ans)
    width_x <- IRanges::IntegerList(split(width(partitioning),
                                          names(partitioning)))
    m <- match(unique(names(from)),names(width_x))
    width_x <- width_x[m]
    width_ans <- sum(width_x)
    ans <- relist(unname(unlist(ans)),
                  IRanges::PartitioningByWidth(width_ans,
                                               names = names(width_ans)))
    unique(ans)
}

#' @rdname positionSequence
#' @export
setMethod("positionSequence","Ranges",
          function(x, order = FALSE, decreasing = FALSE){
              class <- paste0(class(x),"List")
              x <- do.call(class,list(x))
              names(x) <- "tmp"
              ans <- .seqs_rl_strand(x, order = order, decreasing = decreasing)
              ans[["tmp"]]
          })

#' @rdname positionSequence
#' @export
setMethod("positionSequence","RangesList",
          function(x, order = FALSE, decreasing = FALSE){
              .seqs_rl_strand(x, order = order, decreasing = decreasing)
          })

# rescale ----------------------------------------------------------------------

#' @name rescale
#'
#' @title Rescaling \code{Ranges} object
#'
#' @description
#' \code{rescale()} rescales \code{IRanges}, \code{GRanges}, \code{IRangesList}
#' and \code{GRangesList} by using \code{to} and \code{from}.
#'
#' @param x a \code{IRanges}, \code{GRanges}, \code{IRangesList} and
#'   \code{GRangesList} object
#' @param to,from an \code{IRanges} object or a numeric vector parallel to
#' \code{x} or with \code{length = 1L}.
#'
#' @return an object of the same type and dimensions as \code{x}
#'
#' @examples
#' x <- IRanges("5-10")
#' # widen the ranges
#' rescale(x, 100, 10)
#' # widen and shift
#' rescale(x, "31-60", "5-14")
NULL

.zoom0 <- function(x, z = 1)
{
    stopifnot(is(x, "Ranges"), is.numeric(z))
    if (length(z) > length(x) && length(z) != 1L)
        stop("'z' is longer than 'x'")
    if (anyNA(z) || min(z) <= -1L)
        stop("'z' contains NAs and/or negative values")
    new_start <- as.integer(start(x) * z)
    new_width <- as.integer(width(x) * z)
    BiocGenerics:::replaceSlots(x, start=new_start, width=new_width)
}

.normarg_scale <- function(scale)
{
    if (is(scale, "IRanges"))
        return(scale)
    if (is.numeric(scale))
        return(IRanges(1L, scale))
    as(scale, "IRanges")
}

#' @rdname rescale
#' @export
setMethod("rescale","IRanges",
    function(x, to = 1L, from = 1L){
        to <- .normarg_scale(to)
        from <- .normarg_scale(from)
        ans <- shift(x, -start(from))
        ans <- .zoom0(ans, width(to) / width(from))
        shift(ans, start(to))
    }
)

#' @rdname rescale
#' @export
setMethod("rescale","IRangesList",
    function(x, to = 1L, from = 1L){
        to <- .normarg_scale(to)
        from <- .normarg_scale(from)
        ans <- shift(x, -start(from))
        ans <- Map(.zoom0,ans,width(to) / width(from))
        ans <- IRanges::IRangesList(ans)
        shift(ans, start(to))
    }
)

#' @rdname rescale
#' @export
setMethod("rescale","GRanges",
    function(x, to = 1L, from = 1L){
        ranges(x) <- rescale(ranges(x), to = to, from = from)
        x
    }
)

#' @rdname rescale
#' @export
setMethod("rescale","GRangesList",
    function(x, to = 1L, from = 1L){
        ranges(x) <- rescale(ranges(x), to = to, from = from)
        x
    }
)
