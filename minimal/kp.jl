## ---------------------------------------------------
## Copyright (C) 2001-2007 Quantum ESPRESSO group
## This file is distributed under the terms of the
## GNU General Public License. See the file `License'
## in the root directory of the present distribution,
## or http://www.gnu.org/copyleft/gpl.txt .
## ---------------------------------------------------

include("fortran_util.jl")

include("QE.jl")



function errore( calling_routine, message, ierr )
    # ... This is a simple routine which writes an error message to output: 
    # ... if ierr <= 0 it does nothing,
    # ... if ierr  > 0 it stops.
    if ierr <= 0
        return
    else
        throw("Error in routine $calling_routine ($ierr) : $message")
        return
    end
end

##

function kpoint_grid( 
    KPSETTING::Tuple,
    nrot,               # number of bravais lattice symmetries
    s, 
    t_rev, 
    bg,
    npk::Int;           # max number of k-points
    time_reversal = true
    )

    (k1, k2, k3, nk1, nk2, nk3) = KPSETTING
    nkr = nk1*nk2*nk3
    
    #!! coordinates of k points
    #REAL(DP), INTENT(out) :: xk(3,npk)
    xk = zeros(Float64, 3, npk)
    
    #!! weight of k points
    #REAL(DP), INTENT(out) :: wk(npk)
    wk = zeros(Float64, npk)

    # REAL(DP), PARAMETER :: eps=1.0d-5  #N[59]
    eps = 1e-5

    #!!  xkg are the components of the complete grid in crystal axis  #N[70]
    #allocate (xkg( 3,nkr),wkk(nkr))  #N[62]
    wkk = ones(Float64,nkr)
    xkg = ones(Float64,3,nkr)

    # REAL(DP) :: xkr(3), fact, xx, yy, zz #N[54]
    xkr = [0.0, 0.0, 0.0]

    #allocate (equiv( nkr))  #N[63]
    #for nk = 1 : nkr  #N[81]
    #    equiv[nk] = nk  #N[82]
    #end
    equiv = collect(1:nkr) # equivalent to itself

    #!!  this is nothing but consecutive ordering
    @inline n(i,j,k) = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1  
    for i = 1 : nk1
        for j = 1 : nk2
            for k = 1 : nk3
                ##  xkg are the components of the complete grid in crystal axis  #N[70]
                xkg[1,n(i,j,k)] = (i-1+k1//2)//nk1
                xkg[2,n(i,j,k)] = (j-1+k2//2)//nk2
                xkg[3,n(i,j,k)] = (k-1+k3//2)//nk3
            end
        end
    end

    for nk = 1 : nkr
        ##  check if this k-point has already been found equivalent to another
        if ( equiv[nk] ==  nk )
            wkk[nk] = 1.0
            ##  check if there are equivalent k-point to this in the list
            ##  (excepted those previously found to be equivalent to another)
            ##  check both k and -k
            for ns = 1 : nrot
                for i = 1 : 3
                    xkr[i] = s[i,1,ns] * xkg[1,nk] + s[i,2,ns] * xkg[2,nk] + s[i,3,ns] * xkg[3,nk]
                    xkr[i] = xkr[i] - nint(xkr[i])
                end
                if (t_rev[ns] == 1) xkr = -xkr end
                xx = xkr[1]*nk1 - k1//2
                yy = xkr[2]*nk2 - k2//2
                zz = xkr[3]*nk3 - k3//2
                in_the_list = abs(xx-nint(xx))<=eps  &&  abs(yy-nint(yy))<=eps  &&  abs(zz-nint(zz))<=eps
                if ( in_the_list )  #N[110]
                    i = mod( nint( xkr[1]*nk1 - k1//2 + 2nk1 ), nk1 ) + 1
                    j = mod( nint( xkr[2]*nk2 - k2//2 + 2nk2 ), nk2 ) + 1
                    k = mod( nint( xkr[3]*nk3 - k3//2 + 2nk3 ), nk3 ) + 1
                    if ( n(i,j,k)>nk  &&  equiv[n(i,j,k)]==n(i,j,k) )  #N[115]
                        equiv[n(i,j,k)] = nk
                        wkk[nk] = wkk[nk]+1.0
                    else
                        if (equiv[n(i,j,k)]!=nk  ||  n(i,j,k)<nk)
                            errore("kpoint_grid", "something wrong in the checking algorithm", 1)
                        end  #N[120]
                    end
                end
                if ( time_reversal )
                    xx = -xkr[1]*nk1 - k1//2
                    yy = -xkr[2]*nk2 - k2//2
                    zz = -xkr[3]*nk3 - k3//2
                    in_the_list = abs(xx-nint(xx))<=eps  &&  abs(yy-nint(yy))<=eps  &&  abs(zz-nint(zz))<=eps
                    if ( in_the_list )  #N[129]
                        i = mod( nint( xkr[1]*nk1 - k1//2 + 2nk1 ), nk1 ) + 1
                        j = mod( nint( xkr[2]*nk2 - k2//2 + 2nk2 ), nk2 ) + 1
                        k = mod( nint( xkr[3]*nk3 - k3//2 + 2nk3 ), nk3 ) + 1
                        if ( n(i,j,k)>nk  &&  equiv[n(i,j,k)]==n(i,j,k) )   #N[134]
                            equiv[n(i,j,k)] = nk
                            wkk[nk] = wkk[nk]+1.0
                        else  #N[137]
                            if (equiv[n(i,j,k)]!=nk  ||  n(i,j,k)<nk)
                                errore("kpoint_grid", "something wrong in the checking algorithm", 2)
                            end  #N[139]
                        end
                    end
                end
            end
        end
    end

    #!! count irreducible points and order them  #N[148]
    nks = 0
    fact = 0.0
    for nk = 1 : nkr  #N[152]
        if ( equiv[nk] == nk )
            nks = nks+1
            if (nks>npk)
                errore("kpoint_grid","too many k-points",1)
            end
            wk[nks] = wkk[nk]
            fact    = fact + wk[nks]
            #!!!!!  bring back into to the first BZ  #N[158]
            for i = 1 : 3
                xk[i,nks] = xkg[i,nk]-nint(xkg[i,nk])
            end
        end
    end

    #!!  go to cartesian axis (in units 2pi/a0)
    cryst_to_cart( nks, xk, bg, 1 )  #N[165]

    #!!  normalize weights to one
    wk = copy(wk)./fact

    return nks, xk, wk
end
