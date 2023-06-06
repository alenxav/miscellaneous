
# Fortran code
code = '
subroutine eigvals(n, a, vre, vim, info)
    implicit none
    integer :: n, info
    integer, parameter :: lwork = 65536
    double precision :: a(n, n), vre(n), vim(n)
    double precision, save :: work(lwork)

    call dgeev("n", "n", n, a, n, vre, vim, 0d0, 1, 0d0, 1, work, lwork, info)
end subroutine'

# Write file
write.table(code,'ev.f90',quote = F,row.names = F,col.names = F)

# Compile code
tmp = paste('R CMD SHLIB -shared',gsub('library/base','bin/x64/Rlapack.dll',system.file()),'ev.f90 -o ev.dll')
system(tmp)

# Load dynamic library
dyn.load("ev.dll")

# Function to call
eigvals <- function(a) {
  if (is.matrix(a) && is.double(a) && nrow(a) == ncol(a)) {
    n <- nrow(a)
    s <- .Fortran("eigvals", n = as.integer(n), a = a, vre = double(n), vim = double(n), info = 0L)
    structure(complex(real = s$vre, imaginary = s$vim), info = s$info)
  } else stop("Invalid input")
}

# Test
A = crossprod(matrix(rnorm(25),5,5))
eigvals(A)
