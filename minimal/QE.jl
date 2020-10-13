using Libdl

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


global const __QE_LIB__ = "/home/dabajabaza/jianguoyun/Workspace/QE/minimal/lib/QE_minimal.so"

QE_lib = dlopen(__QE_LIB__)

@inline mkpr(sss) = Symbol(sss[3]) => [sss[2] , ( dlsym(QE_lib, Symbol(sss[3]), throw_error=false), 
                                                  Ptr{Nothing}(parse(Int,sss[1],base=16)) )]

QE_content = readlines(`nm $__QE_LIB__`)

QE_functions  = Dict(mkpr(split(l,' ',keepempty=false)) for l in QE_content if startswith(l,"0"))

dlclose(QE_lib)


ss = zeros(Int32, 3,3,48)
ibrav_ = 0
celldm_ = zeros(Float64, 6)
(a_, b_, c_, cosab_, cosac_, cosbc_) = (1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
trd_ht = false
rd_ht = zeros(Float64, 3,3)

ccall(  (:get_symm_base_s_,__QE_LIB__), 
        Cvoid, 
        (Ref{Int32}, Ref{Float64}, Ref{Float64},  Ref{Float64},  Ref{Float64},  Ref{Float64}, Ref{Float64}, Ref{Float64}, 
        Ref{Bool}, Ref{Float64}, Ref{Int32} ), 
        ibrav_,      celldm_,      a_,            b_,            c_,            cosab_,       cosac_,       cosbc_,       
        trd_ht,      rd_ht,        ss )