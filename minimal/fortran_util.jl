global const __FLIB__ = "./fortran_util.so"


nint(r::Float64) = ccall((:nint_fortran_,__FLIB__), Cintmax_t, (Ref{Float64},), r)

