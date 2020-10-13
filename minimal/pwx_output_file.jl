# pwx_output_file.jl

function verify_QE_result(results, program)
    if program.exec[1]=="pw.x"
        return occursin("JOB DONE.", results[end-1])
    else
        return true
    end
end


## -------------------------------------------------------------
## process pw.x results
## -------------------------------------------------------------

function find_line(res, cont_str)
    findfirst(x->occursin(cont_str,x), res)
end


function extract(res, reg1, reg2)
    pos  = findfirst(x->occursin(reg1,x), res)
    m1   = (pos !== nothing ? match(reg1,res[pos]).match : "")
    m2   = (m1 != "" ? match(reg2,m1).match : "")
    m2
end


function extract_all(res, reg1, reg2)
    pos  = findall(x->occursin(reg1,x), res)
    m1   = [match(reg1,res[p]).match for p in pos]
    m2   = [match(reg2,m).match for m in m1]
    m2
end

parse_float64(x) = (try parse(Float64,strip(x)) catch _ NaN end)

parse_3_float64(x) = (try parse.(Float64, split(strip(x)," ",keepempty=false)) catch _ Vector([NaN,NaN,NaN]) end)

function parse_1_string_3_float64(x)
    x1 = split(x, " ", keepempty=false, limit=2)
    return (string(x1[1]) => parse_3_float64(x1[2]))
end

num_f_rstr = r"-?\d+\.\d+"
total_energy_line_rstr = r"\!\s*total\s+energy\s*\=\s*-?\d+\.\d+\s+Ry"
total_energy_line1_rstr = r"\s*total\s+energy\s*\=\s*-?\d+\.\d+\s+Ry"
Harris_Foulkes_estimate_line_rstr = r"\s*Harris-Foulkes\s+estimate\s*\=\s*-?\d+\.\d+\s+Ry"
estimated_scf_accuracy_line_rstr = r"\s*estimated\s+scf\s+accuracy\s*\<\s*\d+\.\d+\s+Ry"
one_electron_contribution_line_rstr = r"\s*one-electron\s+contribution\s*\=\s*-?\d+\.\d+\s+Ry"
hartree_contribution_line_rstr = r"\s*hartree\s+contribution\s*\=\s*-?\d+\.\d+\s+Ry"
xc_contribution_line_rstr = r"\s*xc\s+contribution\s*\=\s*-?\d+\.\d+\s+Ry"
ewald_contribution_line_rstr = r"\s*ewald\s+contribution\s*\=\s*-?\d+\.\d+\s+Ry"


function extract_energy(result_lines)
    endp = find_line(result_lines, "End of self-consistent calculation")
    res  = result_lines[endp+1:end]
    en = Dict(
        "total_energy"                  => parse_float64(extract(result_lines, total_energy_line_rstr,              num_f_rstr)),
        "Harris_Foulkes_estimate"       => parse_float64(extract(result_lines, Harris_Foulkes_estimate_line_rstr,   num_f_rstr)),
        "estimated_scf_accuracy"        => parse_float64(extract(result_lines, estimated_scf_accuracy_line_rstr,    num_f_rstr)),
        "one_electron_contribution"     => parse_float64(extract(result_lines, one_electron_contribution_line_rstr, num_f_rstr)),
        "hartree_contribution"          => parse_float64(extract(result_lines, hartree_contribution_line_rstr,      num_f_rstr)),
        "xc_contribution"               => parse_float64(extract(result_lines, xc_contribution_line_rstr,           num_f_rstr)),
        "ewald_contribution"            => parse_float64(extract(result_lines, ewald_contribution_line_rstr,        num_f_rstr)),
    )
    return en
end


extract_total_energy(result_lines) = extract_energy(result_lines)["total_energy"]


extract_total_energy_history(result_lines) = parse_float64.(extract_all(result_lines, total_energy_line1_rstr, num_f_rstr))


##

function parse_cell_parameters(lines)
    three_numbers = r"\s*-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+"
    @assert (occursin("CELL_PARAMETERS", lines[1]) 
          && occursin(three_numbers, lines[2]) 
          && occursin(three_numbers, lines[3]) 
          && occursin(three_numbers, lines[4]) )  "Invalid block of CELL_PARAMETERS."
    #
    alat  = parse_float64(extract(lines, r"\(\s*alat\s*\=\s*\d+\.\d+\s*\)", num_f_rstr))
    basis =  parse_3_float64.(lines[2:4])
    return alat, basis
end


function parse_atomic_positions(lines)
    one_str_three_numbers = r"\s*\w+\s*-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+"
    @assert (occursin(r"\s*ATOMIC\_POSITIONS\s*\(\w+\)", lines[1]) 
          && all([occursin(one_str_three_numbers, l) for l in lines[2:end]]))  "Invalid block of ATOMIC_POSITIONS."
    #
    unit = extract(lines, r"\s*ATOMIC\_POSITIONS\s*\(\w+\)", r"\(\w+\)")
    ap   = (parse_1_string_3_float64.(lines[2:end]))
    return unit, ap
end


function extract_bfgs_geom_opt(result_lines)
    BFGS_end      = 
    p             = find_line(result_lines, "End of BFGS Geometry Optimization")
    @assert p !== nothing   "The output file does not contain a line \"End of BFGS Geometry Optimization\"."
    result_lines1 = result_lines[(p+1):end]
    p             = find_line(result_lines1, "Begin final coordinates")
    @assert p !== nothing   "The output file does not contain a line \"Begin final coordinates\"."
    result_lines2 = result_lines1[(p+1):end]
    q             = find_line(result_lines2, "End final coordinates")
    @assert q !== nothing   "The output file does not contain a line \"End final coordinates\"."
    result_lines3 = result_lines2[1:q-1]
    @assert length(result_lines3) > 0 "The output file contains nothing in the final result section."

    cp_line_1   = find_line(result_lines3, "CELL_PARAMETERS")
    alat, basis = (-1e10, Vector{Float64}[])  ## "relax"
    if cp_line_1 !== nothing ##  "vc-relax" 
        alat, basis = parse_cell_parameters(result_lines3[cp_line_1:cp_line_1+3])
    end
    ap_line_1   = find_line(result_lines3, "ATOMIC_POSITIONS")
    unit, ap    = parse_atomic_positions(result_lines3[ap_line_1:end])

    data = Dict("alat"             => alat,
                "ap_unit"          => unit,
                "ATOMIC_POSITIONS" => ap,
                "CELL_PARAMETERS"  => basis)

    return data
end


num_kp_rstr = r"\s*number\s+of\s+k\s+points\s*\=\s*\d+"

extract_num_kpoints(res) = extract(res, num_kp_rstr, r"\d+")

##

