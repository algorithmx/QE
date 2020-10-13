global const fortran_build_in = ["date_and_time",
"DGEMV",
"DGER",
"dgetrf",
"dgetri",
"zgetrf",
"zgetri",
]


@inline comm(ls)                = [l for l in ls if !startswith(lstrip(string(l[3])),"!")]
@inline trim_comment(l)         = first(split(l,"!",keepempty=false))
@inline process_call1(l)        = last(strip.(split(trim_comment(l[3]), r"(CALL|call|Call)",keepempty=false)))
@inline process_call2(l)        = strip(first(split(l,r"\s*(\(|\&|\,|\!)",keepempty=false)))
@inline process_subroutine1(l)  = last(strip.(split(trim_comment(l[3]), r"(SUBROUTINE|subroutine|Subroutine|INTERFACE)",keepempty=false)))
@inline process_subroutine2(l)  = strip(first(split(l,r"\s*(\(|\&|\,|\!)",keepempty=false)))

lines1 = try readlines(`grep -nrw .  -e "CALL"`) catch _ String[] end
lines2 = try readlines(`grep -nrw .  -e "Call"`) catch _ String[] end
lines3 = try readlines(`grep -nrw .  -e "call"`) catch _ String[] end
clines = [lines1; lines2; lines3]

calls = comm(map(x->strip.(x), split.(clines, ":", keepempty=false)))

demanded = unique(process_call2.(process_call1.(calls)))

lines4 = try readlines(`grep -nrw .  -e "^\s*SUBROUTINE"`) catch _ String[] end
lines5 = try readlines(`grep -nrw .  -e "^\s*subroutine"`) catch _ String[] end
lines6 = try readlines(`grep -nrw .  -e "^\s*Subroutine"`) catch _ String[] end
lines7 = try readlines(`grep -nrw .  -e "^\s*INTERFACE"`) catch _ String[] end
slines = [lines4; lines5; lines6; lines7]

subroutines = comm(map(x->strip.(x), split.(slines, ":", keepempty=false)))

defined = unique(process_subroutine2.(process_subroutine1.(subroutines)))

wrong = [f for f in demanded if occursin(' ',f)]

undefined = [f for f in demanded if (f ∉ defined && f ∉ wrong && f ∉ fortran_build_in)]

println.(undefined)

##

