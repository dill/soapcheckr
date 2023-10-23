#' Quickly find knots that will cause issues for soap-film smoother
#' Setting up a soap film smoother is often a hard and frustrating process.
#' This function checks that all knots, not just one knot are within the boundary
#' and futher enough away from the boundary to effectively use a soap-film smoother
#'
#' @param bnd A list with sub-lists that will be the boundary the soap-film smoother
#' will smooth in. Coordinates need to be in a reference system that has meaningful units.
#' For example, coordinates could be in UTMs as UTMs rely on metres.
#' @param knots A dataframe with two columns that are the coordinates of the knots
#' that are to be supplied to the soap-film smoother.
#' Coordinates need to be be in meaniful units. For example, coordinates could be
#' in UTMs as UTMs rely on metres.
#' @param nmax If supplying `nmax` argument to `mgcv::gam()`, provide that
#' same value to `autocruncher()`.
#' @param k If supplying `k` argument to `mgcv::gam()`, provide that
#' same value to `autocruncher()`.
#' @param xname Column name of x coordinate for the knots object
#' @param yname Column name of y coordinate for the knots object
#'
#'
#'
#' @return A vector of the indices of all of the knots that either fall outside
#' the boundary or too close to the boundary and will cause for the soap-film
#' smoother to error.
#'
#' @references
#' This function is based on soap.r in {mgcv}, Copyright Simon Wood 2006-2012.
#' Bugs fixes added by David L Miller, 2019
#'
#' @examples
#' # crunch_ind <- autocruncher(fsb, knots, k=30, xname="v", yname="w"), will add more
#' crunch_ind <- autocruncher(sissabagama_lake_ls, sissabagama_lake_knots)
#'
#' @import mgcv
#' @export
#'


# NB in mgcv nmax is 100 here, where as the default is 200
autocruncher <- function(bnd,
                         knots,
                         k = 10,
                         nmax = 200,
                         xname = "x",
                         yname = "y"){
  ## setup soap film  smooth - nmax is number of grid cells for longest side
  ## it's important that grid cells are square!

  autocrunch.knots <- function(G, knots, x0, y0, dx, dy){
    ## finds indices of knot locations in solution grid
    ## the knot x,y locations are given in the `knots' argument.
    badk <- NULL
    nk <- length(knots$x)
    ki <- rep(0, nk)
    nx <- ncol(G);
    ny <- nrow(G)
    if (nk == 0) return(NULL)
    for (k in 1:nk) {
      i <- round((knots$x[k] - x0) / dx) + 1
      j <- round((knots$y[k] - y0) / dy) + 1
      if (i > 1 && i <= nx && j > 1 && j <= ny) {
        ki[k] <- G[j, i]
        if (ki[k] <= 0) {
          badk <- c(badk, k)
          #str <- paste("knot",k,"is on or outside boundary")
          #stop(str)
        }
      }
    } ## all knots done
    #ki ## ki[k] indexes kth knot in solution grid
    badk
  } ## end crunch.knots

  # ----- THIS NEEDS TO BE FIXED because if you don't supply autocrunch -----
  # with x or ynames it nulls...considering this is used in soap_check it won't
  # produced the plot...I have tried to figure out what this is doing and i'm stuck
  # so instead I added x_name and y_name ot soap_check and added
  # xname and yname arguments to when autocrunch is called
  # boundary names must be x and y!
  bnd <- lapply(bnd, function(x, xname, yname){
    if(!all(names(x) %in% c("x", "y"))){
      x$x <- x[[xname]]
      x[[xname]] <- NULL
      x$y <- x[[yname]]
      x[[yname]] <- NULL
    }
    x
  }, xname = xname, yname = yname)

  # there was a missing if statement that is why it was nulling out
  if(!all(names(knots) %in% c("x", "y"))) {
  knots$x <- knots[[xname]]
  knots[[xname]] <- NULL
  knots$y <- knots[[yname]]
  knots[[yname]] <- NULL
  }

  ## check boundary...
  if (!inherits(bnd,"list")) stop("bnd must be a list.")
  n.loops <- length(bnd)
  if (n.loops!=length(k)) {
    if (length(k)==1) k <- rep(k,n.loops)
    else stop("lengths of k and bnd are not compatible.")
  }
  bnd <- mgcv:::process.boundary(bnd) ## add distances and close any open loops

  ## create grid on which to solve Laplace equation
  ## Obtain grid limits from boundary 'bnd'....
  x0 <- min(bnd[[1]]$x);x1 <- max(bnd[[1]]$x)
  y0 <- min(bnd[[1]]$y);y1 <- max(bnd[[1]]$y)
  if (length(bnd)>1) for (i in 2:length(bnd)) {
    x0 <- min(c(x0,bnd[[i]]$x)); x1 <- max(c(x1,bnd[[i]]$x))
    y0 <- min(c(y0,bnd[[i]]$y)); y1 <- max(c(y1,bnd[[i]]$y))
  } ## now got the grid limits, can set it up

  if (x1-x0>y1-y0) { ## x is longest side
    dy <- dx <- (x1-x0) /(nmax-1)
    nx <- nmax
    ny <- ceiling((y1-y0)/dy)+1
  } else { ## y is longest side
    dy <- dx <- (y1-y0) /(nmax-1)
    ny <- nmax
    nx <- ceiling((x1-x0)/dy)+1
  }
  ## so grid is now nx by ny, cell size is dx by dy (but dx=dy)
  ## x0, y0 is "lower left" cell centre

  ## Create grid index G
  bnc <- mgcv:::bnd2C(bnd) ## convert boundary to form required in C code

  G <- matrix(0, ny, nx)
  nb <- rep(0,bnc$n.loop)

  oo <- .C(mgcv:::C_boundary,G=as.integer(G), d=as.double(G), dto=as.double(G),
           x0=as.double(x0),
           y0 = as.double(y0), dx=as.double(dx), dy = as.double(dy),
           nx=as.integer(nx),as.integer(ny), x=as.double(bnc$x),y=as.double(bnc$y),
           breakCode=as.double(bnc$breakCode),n=as.integer(bnc$n),nb=as.integer(nb))

  ret <- list(G=matrix(oo$G,ny,nx),nb=oo$nb,d=oo$d[oo$d >= 0],
              x0=x0,y0=y0,dx=dx,dy=dy,bnd=bnd)
  rm(oo)
  ## Now create the PDE coefficient matrix
  #  n.inside <- sum(ret$G > - nx*ny)
  #  xx <- rep(0,5*n.inside)
  #  o1 <- .C(C_pde_coeffs,as.integer(ret$G),xx=as.double(xx),ii=as.integer(xx),jj=as.integer(xx),
  #            n=as.integer(0),as.integer(nx),as.integer(ny),as.double(dx),as.double(dy))
  #  ind <- 1:o1$n
  #  X <- sparseMatrix(i=o1$ii[ind]+1,j=o1$jj[ind]+1,x=o1$xx[ind])
  #  er <- expand(lu(X))
  #  ret$Q <- er$Q;ret$U <- er$U;ret$L <- er$L;ret$P <- er$P
  #  ret$ng <- n.inside ## the number of cells to solve for
  #  rm(er);rm(X)
  ## ... so the sparse LU decomposition of X can be used to solve PDE.
  ## X = PLUQ where P and Q are permuation matrices.

  ## now obtain location of knots in solution ...
  ret <- autocrunch.knots(ret$G, knots, x0, y0, dx, dy)

  ret
} ## end of setup.soap

