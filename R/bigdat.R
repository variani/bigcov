#' S3 class bigdat.
#' @exportClass bigdat
#

#' S3 class bigdatMat.
#' @exportClass bigdatMat
#

#' S3 class bigdatBigMat.
#' @exportClass bigdatBigMat
#

#' S3 class bigdatBEDMatrix.
#' @exportClass bigdatBEDMatrix
#

#--------------
# Constructors
#--------------

#' @export
bigdat <- function(x, rows = NULL, 
  num_batches = NULL, batch_size = NULL, 
  path = ".", 
  ...)
{
  ### args
  stopifnot(!any(duplicated(rows)))
  
  ### initialize 
  out <- list(path = path)
  
  ### `x` as list
  if(class(x)[1] != "list") {
    x <- list(x)
  }
  
  # define `class`
  x_classes <- sapply(x, function(elem) class(elem)[1])
  stopifnot(all(x_classes == x_classes[1]))
  
  out$class <- x_classes[1]
  
  if(out$class == "matrix") {
    out$data <- x

    oldClass(out) <- c("bigdatMat", "bigdat")
  } else if(out$class == "big.matrix") {
    out$data <- lapply(x, describe)
    
    oldClass(out) <- c("bigdatBigMat", "bigdat")
  } else if(out$class == "big.matrix.descriptor") {
    out$data <- x
    
    oldClass(out) <- c("bigdatBigMat", "bigdat")
  } else if(out$class == "BEDMatrix") {
    out$data <- x
    
    oldClass(out) <- c("bigdatBEDMatrix", "bigdat")
  } else {
   stop("not supported class: ", out$class)
  }

  ### `num_batches`
  if(!is.null(num_batches)) {
    stopifnot(is.null(batch_size))
    
    num_cols_sum <- bigdat_ncols_sum(out)
    batch_size <- floor(num_cols_sum / num_batches)
  }
  
  ### rows
  nrows <- bigdat_nrows(out)
  stopifnot(all(nrows == nrows[1]))
  
  if(is.null(rows)) {
    rows <- seq(1, bigdat_nrow1(out))
  }
  out$rows <- rows
  out$all_rows <- (length(out$rows) == bigdat_nrow1(out))
  
  ### batches
  num_elems <- bigdat_nelem(out)
  
  if(is.null(batch_size)) {
    batches <- data_frame(batch = seq(1, num_elems),
      elem = seq(1, num_elems),
      beg = rep(1, num_elems),
      end = bigdat_ncols(out))
  } else {
    batches <- lapply(seq(1, num_elems), function(elem) {
      num_cols_elem <- bigdat_ncol(out, elem = elem)
      batch_size_elem <- batch_size
      
      if(batch_size_elem > num_cols_elem) {
        batch_size_elem <- num_cols_elem
      }
      
      beg <- seq(1, num_cols_elem, by = batch_size_elem) 
      end <- c(beg[-1] - 1, num_cols_elem)  
 
      stopifnot(all(beg <= num_cols_elem))
      stopifnot(all(end <= num_cols_elem))
      
      data_frame(elem = elem, beg = beg, end = end)  
    })
    batches <- bind_rows(batches)
    
    # add `batch` column
    batches <- mutate(batches, 
      batch = seq(1, nrow(batches))) %>%
    select(batch, everything())
  }
  
  out$batches <- batches
  
  ### return
  return(out)
}

#-----------------
# Export methods
#-----------------

#' @export bigdat_attach
bigdat_attach <- function(x, ...) UseMethod("bigdat_attach")

#' @export bigdat_nelem
bigdat_nelem <- function(x, ...) UseMethod("bigdat_nelem")

#' @export bigdat_nrow
bigdat_nrow <- function(x, ...) UseMethod("bigdat_nrow")

#' @export bigdat_nrow1
bigdat_nrow1 <- function(x, ...) UseMethod("bigdat_nrow1")

#' @export bigdat_nrows
bigdat_nrows <- function(x, ...) UseMethod("bigdat_nrows")

#' @export bigdat_ncol
bigdat_ncol <- function(x, ...) UseMethod("bigdat_ncol")

#' @export bigdat_ncol1
bigdat_ncol1 <- function(x, ...) UseMethod("bigdat_ncol1")

#' @export bigdat_ncols
bigdat_ncols <- function(x, ...) UseMethod("bigdat_ncols")

#' @export bigdat_ncols_sum
bigdat_ncols_sum <- function(x, ...) UseMethod("bigdat_ncols_sum")

#' @export bigdat_nbatch
bigdat_nbatch <- function(x, ...) UseMethod("bigdat_nbatch")

#' @export bigdat_slice
bigdat_slice <- function(x, ...) UseMethod("bigdat_slice")

#' @export bigdat_data
bigdat_data <- function(x, ...) UseMethod("bigdat_data")

#' @export bigdat_batch
bigdat_batch <- function(x, ...) UseMethod("bigdat_batch")


#--------------
# print, show
#--------------

print.bigdat <- function(x, ...)
{
  cat(" ", class(x)[1], " (", class(x$data[[1]]), ")\n", sep = "")

  cat("  - dimensions:", bigdat_nrow1(x), "x", bigdat_ncols_sum(x), "\n")
  
  ### batches
  num_batches <- bigdat_nbatch(x)
  num_top <- 5
  
  if(num_batches < num_top) {
    cat("\n - batches:\n")
    x$batches %>% head(num_top) %>% as.data.frame %>% print
  } else {
    cat("\n - top batches ", num_top, " / ", format(num_batches, big.mark = ",") , ":\n", sep = "")
    x$batches %>% head(num_top) %>% as.data.frame %>% print
  }
}

#--------------------------------------
# pred_attach (`bigdatBigMat` class)
#--------------------------------------

bigdat_attach.bigdatBigMat <- function(x, elem, ...) 
{
  stopifnot(class(x$data[[elem]]) == "big.matrix.descriptor")
  bmat <- attach.big.matrix(x$data[[elem]], path = x$path)
  
  stopifnot(!is.nil(bmat@address))
  
  return(bmat)
}

#-----------------------------------
# bigdat_nrow(s) & bigdat_ncol(s)
#-----------------------------------

bigdat_nelem.bigdat <- function(x, ...) length(x$data)

bigdat_nrow.bigdat <- function(x, elem, ...) nrow(x$data[[elem]])
bigdat_nrow.bigdatBigMat <- function(x, elem, ...) bigdat_attach(x, elem) %>% nrow

bigdat_nrow1.bigdat <- function(x, ...) bigdat_nrow(x, elem = 1, ...)

bigdat_nrows.bigdat <- function(x, ...) sapply(x$data, nrow)
bigdat_nrows.bigdatBigMat <- function(x, ...) sapply(seq(1, length(x$data)), 
  function(elem) bigdat_attach(x, elem) %>% nrow)

bigdat_ncol.bigdat <- function(x, elem, ...) ncol(x$data[[elem]])
bigdat_ncol.bigdatBigMat <- function(x, elem, ...) bigdat_attach(x, elem) %>% ncol

bigdat_ncol1.bigdat <- function(x, ...) bigdat_ncol(x, elem = 1, ...)

bigdat_ncols.bigdat <- function(x, ...) sapply(x$data, ncol)
bigdat_ncols.bigdatBigMat <- function(x, ...) sapply(seq(1, length(x$data)), 
  function(elem) bigdat_attach(x, elem) %>% ncol)
  
bigdat_ncols_sum.bigdat <- function(x, ...) bigdat_ncols(x, ...) %>% sum

bigdat_nbatch.bigdat <- function(x, ...) nrow(x$batches)

#-----------------------------------
# bigdat_slice, bigdat_data, bigdat_batch
#-----------------------------------

bigdat_slice.bigdat <- function(x, elem, rows, cols, ...) x$data[[elem]][rows, cols, drop = FALSE]

bigdat_slice.bigdatBigMat <- function(x, elem, rows, cols, ...) bigdat_attach(x, elem)[rows, cols, drop = FALSE]


bigdat_data.bigdat <- function(x, elem, rows, cols, ...)
{
  ### args
  if(missing(elem)) {
    if(bigdat_nelem(x) == 1) {
      elem <- 1
    } else {
      stop("missing(elem) & bigdat_nelem(x) != 1")
    }
  }
    
  if(missing(rows)) {
    rows <- x$rows
  }
  
  if(missing(cols)) {
    cols <- seq(1, bigdat_ncol(x, elem))
  }
    
  ### get data matrix `mat`
  mat <- bigdat_slice(x, elem, rows, cols)
  
  return(mat)
}

bigdat_batch.bigdat <- function(x, batch, rows, ...)
{
  ### args
  if(missing(batch)) {
    if(bigdat_nbatch(x) == 1) {
      batch <- 1
    } else {
      stop("missing(batch) & bigdat_nbatch(x) != 1")
    }
  }
  
  stopifnot(batch > 0 & batch <= bigdat_nbatch(x))
  
  ### processing
  batch_par <- batch
  batches_row <- filter(x$batches, batch == batch_par)
  stopifnot(nrow(batches_row) == 1)
  
  elem <- batches_row$elem
  cols <- seq(batches_row$beg, batches_row$end)
  
  bigdat_data(x, elem = elem, rows = rows, cols = cols, ...) 
}


