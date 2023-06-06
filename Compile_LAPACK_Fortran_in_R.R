code = '
subroutine eigvals(n, a, vre, vim, info)
    implicit none
    integer :: n, info
    integer, parameter :: lwork = 65536
    double precision :: a(n, n), vre(n), vim(n)
    double precision, save :: work(lwork)

    call dgeev("n", "n", n, a, n, vre, vim, 0d0, 1, 0d0, 1, work, lwork, info)
end subroutine'

write.table(code,'ev.f90',quote = F,row.names = F,col.names = F)

system('R CMD SHLIB -shared C:/Temp/R/bin/x64/Rlapack.dll ev.f90 -o ev.dll')

dyn.load("ev.dll")

eigvals <- function(a) {
  if (is.matrix(a) && is.double(a) && nrow(a) == ncol(a)) {
    n <- nrow(a)
    s <- .Fortran("eigvals", n = as.integer(n), a = a, vre = double(n), vim = double(n), info = 0L)
    structure(complex(real = s$vre, imaginary = s$vim), info = s$info)
  } else stop("Invalid input")
}

A = crossprod(matrix(rnorm(25),5,5))
eigvals(A)
