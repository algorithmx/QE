global const man_nm = Dict(
"A" => "The symbol's value is absolute, and will not be changed by further linking.",
"B" => "The symbol is in the BSS data section.  This section typically contains zero-initialized or uninitialized data, although the exact behavior is system dependent.",
"b" => "The symbol is in the BSS data section.  This section typically contains zero-initialized or uninitialized data, although the exact behavior is system dependent.",
"C" => "The symbol is common.  Common symbols are uninitialized data.  When linking, multiple common symbols may appear with the same name.  If the symbol is defined anywhere, the common symbols are treated as undefined references.",
"D" => "The symbol is in the initialized data section.",
"d" => "The symbol is in the initialized data section.",
"G" => "The symbol is in an initialized data section for small objects.  Some object file formats permit more efficient access to small data objects, such as a global int variable as opposed to a large global array.",
"g" => "The symbol is in an initialized data section for small objects.  Some object file formats permit more efficient access to small data objects, such as a global int variable as opposed to a large global array.",
"i" => "For PE format files this indicates that the symbol is in a section specific to the implementation of DLLs.  For ELF format files this indicates that the symbol is an indirect function.  This is a GNU extension to the standard set of ELF symbol types.  It indicates a symbol which if referenced by a relocation does not evaluate to its address, but instead must be invoked at runtime.  The runtime execution will then return the value to be used in the relocation.",
"I" => "The symbol is an indirect reference to another symbol.",
"N" => "The symbol is a debugging symbol.",
"p" => "The symbols is in a stack unwind section.",
"R" => "The symbol is in a read only data section.",
"r" => "The symbol is in a read only data section.",
"S" => "The symbol is in an uninitialized or zero-initialized data section for small objects.",
"s" => "The symbol is in an uninitialized or zero-initialized data section for small objects.",
"T" => "The symbol is in the text (code) section.",
"t" => "The symbol is in the text (code) section.",
"U" => "The symbol is undefined.",
"u" => "The symbol is a unique global symbol.  This is a GNU extension to the standard set of ELF symbol bindings.  For such a symbol the dynamic linker will make sure that in the entire process there is just one symbol with this name and type in use.",
"V" => "The symbol is a weak object.  When a weak defined symbol is linked with a normal defined symbol, the normal defined symbol is used with no error.  When a weak undefined symbol is linked and the symbol is not defined, the value of the weak symbol becomes zero with no error.  On some systems, uppercase indicates that a default value has been specified.",
"v" => "The symbol is a weak object.  When a weak defined symbol is linked with a normal defined symbol, the normal defined symbol is used with no error.  When a weak undefined symbol is linked and the symbol is not defined, the value of the weak symbol becomes zero with no error.  On some systems, uppercase indicates that a default value has been specified.",
"W" => "The symbol is a weak symbol that has not been specifically tagged as a weak object symbol.  When a weak defined symbol is linked with a normal defined symbol, the normal defined symbol is used with no error.  When a weak undefined symbol is linked and the symbol is not defined, the value of the symbol is determined in a system-specific manner without error.  On some systems, uppercase indicates that a default value has been specified.",
"w" => "The symbol is a weak symbol that has not been specifically tagged as a weak object symbol.  When a weak defined symbol is linked with a normal defined symbol, the normal defined symbol is used with no error.  When a weak undefined symbol is linked and the symbol is not defined, the value of the symbol is determined in a system-specific manner without error.  On some systems, uppercase indicates that a default value has been specified.",
"-" => "The symbol is a stabs symbol in an a.out object file.  In this case, the next values printed are the stabs other field, the stabs desc field, and the stab type.  Stabs symbols are used to hold debugging information.",
"?" => "The symbol type is unknown, or object file format specific.",
)

##

function content_lib_a(lib_a::String)
        _content0 = readlines(`nm $lib_a`)
        _content  = [split(l,' ',keepempty=false) for l in _content0]
        @inline mkpr3(sss) = Symbol(sss[3]) => [sss[2] , Ptr{Nothing}(parse(Int,sss[1],base=16))]
        @inline mkpr2(sss) = Symbol(sss[2]) => [sss[1] , Ptr{Nothing}(0)]
        @inline mkpr(sss)  = length(sss)==3 ? mkpr3(sss) : mkpr2(sss)
        _functions  = Dict(mkpr(l) for l in _content if length(l)>1)
        return _functions
end

##

global const __QE_LIB__ = "/home/dabajabaza/jianguoyun/Workspace/QE/minimal/lib/QE_minimal.so"
QE_functions = content_lib_a(__QE_LIB__)

##

include("pwx_input_file.jl")

ss = zeros(Int32, 3,3,48)
(a_, b_, c_, cosab_, cosac_, cosbc_) = (1.2, 1.2, 1.2, 0.0, 0.0, 0.0)
(alpha, beta, gamma) = (180acos(cosab_)/π, 180acos(cosac_)/π, 180acos(cosbc_)/π)
ibrav_ = Find_ibrav(Int_Tables[221], Find_Lattice(a_,b_,c_,alpha,beta,gamma))
celldm_ = celldm_array(ibrav_, (a_, b_, c_, acos(cosab_), acos(cosac_), acos(cosbc_)))
println(celldm_)
##

ccall(  (:get_symm_base_s_,__QE_LIB__), 
        Cvoid, 
        (Ref{Int32}, Ref{Float64}, Ref{Float64},  Ref{Float64},  Ref{Float64},  Ref{Float64}, Ref{Float64}, Ref{Float64}, 
        Ref{Int32} ), 
        ibrav_,      celldm_,      0.0,           0.0,           0.0,           0.0,          0.0,          0.0,       
        ss  )