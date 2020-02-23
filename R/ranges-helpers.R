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
#' @examples
#' library(EpiTxDb)
#' # Returns an integer vector
#' gr <- GRanges("1:5")
#' positionSequence(gr)
#' # returns an IntegerList
#' grl <- GRangesList("1" = gr,"2" = gr,"3" = gr) # must be named
#' positionSequence(gr)
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
              .seqs_rl_strand(x, order = order, decreasing = decreasing)[["tmp"]]
          })

#' @rdname positionSequence
#' @export
setMethod("positionSequence","RangesList",
          function(x, order = FALSE, decreasing = FALSE){
              .seqs_rl_strand(x, order = order, decreasing = decreasing)
          })

