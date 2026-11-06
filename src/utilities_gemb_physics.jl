# GEMB physics members: one tile sweep per unresolved physics configuration, named by GEMB_GlacierSims'
# `physics_member` id `"<grain>-<band decoupling>-<impurity>"` (e.g. `"C-kfit-table"`). Each member is
# calibrated against every altimetry ensemble member, so the outputs of steps 13 and 14 form a
# (run × physics) ensemble. A member's outputs are written to `gemb_physics/<member>/` beside the
# altimetry member's synthesized file, under the altimetry member's own file name, so every file-name
# rule of the altimetry ensemble (`binned_filled_fileparts`) applies unchanged.

# A member id: three or more `-`-separated tokens of letters, digits and `.` naming the GEMB_GlacierSims
# tile sweep (`<grain>-<band decoupling>-<impurity>`, and any further level GEMB_GlacierSims adds),
# optionally followed by `-g<g>`, a precipitation gradient applied here (see
# `gemb_physics_member_parts`). No `_` (file names split on it) and no path separator.
const GEMB_PHYSICS_MEMBER_PATTERN = r"^[A-Za-z0-9.]+(-[A-Za-z0-9.]+){2,}$"

# The token form of a precipitation gradient; never a tile-sweep token.
const _GEMB_GRADIENT_TOKEN = r"^g[0-9]+(\.[0-9]+)?$"

# The `g` token of a precipitation gradient of `g` percent per 100 m: `g2.5`, `g5`.
_gemb_gradient_token(g::Real) = "g" * (isinteger(g) ? string(Int(g)) : string(Float64(g)))

"""
    check_gemb_physics_member(member) -> String

`member` as a `String` when it is a GEMB physics member id: a GEMB_GlacierSims tile-sweep member
(`"<grain>-<band decoupling>-<impurity>"`, `-`-separated tokens of letters, digits and `.`, with any
further level GEMB_GlacierSims adds as a fourth token), optionally followed by `-g<g>`, a
precipitation gradient of `g` > 0 percent per 100 m written in its shortest form (`-g2.5`, `-g5`).
Throws otherwise; the gradient-free member has only the name without a `g` token, so no member has two,
and only the last token may be a gradient.
"""
function check_gemb_physics_member(member)
    member isa AbstractString && occursin(GEMB_PHYSICS_MEMBER_PATTERN, member) ||
        error("$(repr(member)) is not a GEMB physics member id; expected \"<grain>-<band decoupling>-<impurity>\" " *
              "or \"<grain>-<band decoupling>-<impurity>-g<gradient>\", e.g. \"C-kfit-table\" or \"C-kfit-table-g2.5\"")
    parts = split(member, '-')
    any(occursin(_GEMB_GRADIENT_TOKEN, p) for p in parts[1:end-1]) &&
        error("$(repr(member)) is not a GEMB physics member id: a gradient token must be the last")
    if occursin(_GEMB_GRADIENT_TOKEN, last(parts))
        length(parts) >= 4 ||
            error("$(repr(member)) is not a GEMB physics member id; expected \"<grain>-<band decoupling>-<impurity>-g<gradient>\"")
        g = parse(Float64, last(parts)[2:end])
        g > 0 || error("$(repr(member)) is not a GEMB physics member id: the gradient-free member is " *
                       "$(repr(join(parts[1:end-1], '-'))), without a g token")
        last(parts) == _gemb_gradient_token(g) ||
            error("$(repr(member)) is not a GEMB physics member id: write the gradient as " *
                  "$(repr(_gemb_gradient_token(g)))")
    end
    return String(member)
end

"""
    gemb_physics_member_parts(member) -> (; tile_member, precipitation_gradient)

The GEMB_GlacierSims tile-sweep member a GEMB physics member reads (every token but a trailing `g`
token), and the precipitation gradient [% per 100 m] applied to it here (`0.0` without a `g` token).
See [`gemb_band_pscale_multipliers`](@ref).
"""
function gemb_physics_member_parts(member)
    parts = split(check_gemb_physics_member(member), '-')
    occursin(_GEMB_GRADIENT_TOKEN, last(parts)) || return (; tile_member=join(parts, '-'), precipitation_gradient=0.0)
    return (; tile_member=join(parts[1:end-1], '-'), precipitation_gradient=parse(Float64, last(parts)[2:end]))
end

"""
    gemb_physics_member_id(tile_member, precipitation_gradient) -> String

The GEMB physics member that applies `precipitation_gradient` [% per 100 m] to `tile_member`'s sweep;
`tile_member` itself at zero. The inverse of [`gemb_physics_member_parts`](@ref).
"""
function gemb_physics_member_id(tile_member, precipitation_gradient::Real)
    parts = gemb_physics_member_parts(tile_member)
    parts.precipitation_gradient == 0 || error("$(repr(tile_member)) already carries a precipitation gradient")
    precipitation_gradient >= 0 || error("a precipitation gradient is ≥ 0, got $(precipitation_gradient)")
    iszero(precipitation_gradient) && return parts.tile_member
    return check_gemb_physics_member("$(parts.tile_member)-$(_gemb_gradient_token(precipitation_gradient))")
end

"""
    gemb_band_pscale_multipliers(band_centers, w, precipitation_gradient) -> Vector{Float64}

The per-band factors `m_b = 1 + g (z_b − z̄) / 100 m` (g = `precipitation_gradient` / 100) that turn a
geotile's precipitation scaling into each band's, `pscale_b = pscale · m_b`. `z̄` is the mean band
elevation weighted by `w`, the band weights the geotile's series are summed with
(`_gemb_band_weights`), so `Σ w_b m_b = Σ w_b`: the gradient moves precipitation between bands and
leaves the geotile's total unchanged. All ones at `g = 0`.
"""
function gemb_band_pscale_multipliers(band_centers, w, precipitation_gradient::Real)
    length(band_centers) == length(w) ||
        throw(DimensionMismatch("$(length(band_centers)) band centers but $(length(w)) weights"))
    iszero(precipitation_gradient) && return ones(length(band_centers))
    sum(w) > 0 || error("the band weights sum to $(sum(w)); a geotile without ice has no reference elevation")
    zbar = sum(w .* band_centers) / sum(w)
    return 1 .+ (precipitation_gradient / 100) .* (band_centers .- zbar) ./ 100
end

"""
    gemb_gradient_pscale_range(m, w, pscale_axis) -> (lo, hi)

The geotile-mean precipitation scalings at which every band holding ice (`w > 0`) stays on
`pscale_axis` under the multipliers `m` ([`gemb_band_pscale_multipliers`](@ref)):
`pscale · m_b ∈ [minimum(pscale_axis), maximum(pscale_axis)]`. `lo > hi` when none does, which a band
with `m_b ≤ 0` always causes.
"""
function gemb_gradient_pscale_range(m, w, pscale_axis)
    mi = m[w .> 0]
    isempty(mi) && return (minimum(pscale_axis), maximum(pscale_axis))
    minimum(mi) > 0 || return (Inf, -Inf)
    return (minimum(pscale_axis) / minimum(mi), maximum(pscale_axis) / maximum(mi))
end

"""
    gemb_physics_path(path, physics_member) -> String

Where the output `path` of a calibration against `physics_member`'s tile sweep is written:
`<dirname(path)>/gemb_physics/<physics_member>/<basename(path)>`. With `physics_member = nothing`
(a tile sweep run without members), `path` itself.
"""
gemb_physics_path(path, ::Nothing) = path
gemb_physics_path(path, physics_member::AbstractString) =
    joinpath(dirname(path), "gemb_physics", check_gemb_physics_member(physics_member), basename(path))

"""
    gembfit_path(synthesized, physics_member) -> String
    gembfit_dv_path(synthesized, physics_member) -> String

The calibration (step 13, `…_synthesized_gembfit.arrow`) and the calibrated synthesis (step 14,
`…_synthesized_gembfit_dv.jld2`) of the altimetry member `synthesized` (`…_synthesized.jld2`) against
`physics_member`, in that member's folder ([`gemb_physics_path`](@ref)).
"""
gembfit_path(synthesized, physics_member) =
    gemb_physics_path(replace(synthesized, ".jld2" => "_gembfit.arrow"), physics_member)
gembfit_dv_path(synthesized, physics_member) =
    gemb_physics_path(replace(synthesized, ".jld2" => "_gembfit_dv.jld2"), physics_member)

"""
    gemb_physics_member_from_path(path) -> String or nothing

The physics member whose folder (`gemb_physics/<member>/`, [`gemb_physics_path`](@ref)) holds `path`,
or `nothing` when `path` is in no member folder.
"""
function gemb_physics_member_from_path(path)
    parts = splitpath(dirname(path))
    length(parts) >= 2 && parts[end-1] == "gemb_physics" || return nothing
    return check_gemb_physics_member(parts[end])
end

"""
    gemb_physics_source(path) -> (path, physics_member)

The inverse of [`gemb_physics_path`](@ref): `path` with its `gemb_physics/<member>/` folder removed,
and the member, or `(path, nothing)` when it is in no member folder.
"""
function gemb_physics_source(path)
    member = gemb_physics_member_from_path(path)
    isnothing(member) && return (path, nothing)
    return (joinpath(dirname(dirname(dirname(path))), basename(path)), member)
end

"""
    gemb_physics_member(A) -> String or nothing

The physics member the GEMB ensemble `A` (a layer or the DimStack from [`gemb_ensemble_dv`](@ref)) was
built from, read from its metadata; `nothing` for a sweep run without members.
"""
gemb_physics_member(A) = get(DD.metadata(A), "gemb_physics_member", nothing)

"""
    gembfit_dv_files(path2runs_synthesized, physics_members) -> DimArray{String}

The calibrated synthesis files (step 14, `…_synthesized_gembfit_dv.jld2`) of every altimetry ensemble
member in `path2runs_synthesized` (`…_synthesized.jld2`) against every GEMB physics member in
`physics_members`, as a `(run, physics)` DimArray of paths. `run` is labelled with the altimetry
member's synthesized file and `physics` with the member id, so the two ensemble dimensions can be
reduced separately.
"""
function gembfit_dv_files(path2runs_synthesized, physics_members)
    members = check_gemb_physics_member.(collect(physics_members))
    allunique(members) || error("physics members repeat: $(join(members, ", "))")
    runs = collect(String, path2runs_synthesized)
    all(endswith("_synthesized.jld2"), runs) ||
        error("expected `…_synthesized.jld2` altimetry members, got $(join(filter(!endswith("_synthesized.jld2"), runs), ", "))")
    files = [gembfit_dv_path(r, m) for r in runs, m in members]
    return DimArray(files, (Dim{:run}(runs), Dim{:physics}(members)))
end

"""
    runs2rgi_physics(path2runs, physics_members; kwargs...) -> Dict{String,DimArray}

[`runs2rgi`](@ref) for every GEMB physics member in `physics_members`, with the members as a
`physics` dimension beside the altimetry ensemble's `run`: each variable is a DimArray over
`(run, physics, rgi, date)` in the unit `runs2rgi` gives it. `kwargs` are passed to `runs2rgi`.
Throws if two members' regional sums of a variable differ in their dates or units.
"""
function runs2rgi_physics(path2runs, physics_members; kwargs...)
    members = check_gemb_physics_member.(collect(physics_members))
    allunique(members) || error("physics members repeat: $(join(members, ", "))")
    by_member = Dict(m => runs2rgi(path2runs; physics_member=m, kwargs...) for m in members)
    first_sum = by_member[first(members)]
    out = Dict{String,Any}()
    for (varname, A) in first_sum
        drun, drgi, ddate = dims(A, :run), dims(A, :rgi), dims(A, :date)
        units = DD.metadata(A)["units"]
        B = fill(NaN, drun, Dim{:physics}(members), drgi, ddate)
        for m in members
            Am = by_member[m][varname]
            dims(Am, :date) == ddate && DD.metadata(Am)["units"] == units ||
                error("$(varname): physics member $(m) has dates or units other than $(first(members))'s")
            B[physics=At(m)] = parent(Am)
        end
        out[varname] = DD.rebuild(B; metadata=Dict("units" => units, "gemb_physics_members" => members))
    end
    return out
end

# The provenance of a tile sweep read as `physics_member`: a member with a precipitation gradient
# records itself and its gradient in place of the tile sweep's member, so a fit or an ensemble of it is
# never taken for the gradient-free member's (GEMB_GlacierSims compares `physics_*` keys too).
function _gemb_gradient_provenance(provenance, precipitation_gradient)
    iszero(precipitation_gradient) && return provenance
    out = deepcopy(provenance)
    p = out["tile_run_parameters"]
    p["physics_member"] = gemb_physics_member_id(p["physics_member"], precipitation_gradient)
    p["physics_precipitation_gradient"] = Float64(precipitation_gradient)
    return out
end

# Throw unless the `tile_run_parameters` of a tile provenance (`gemb_tile_provenance`) name
# `physics_member`, or name no member when it is `nothing`. `what` names the provenance's source.
function _check_gemb_physics_provenance(parameters, physics_member, what)
    recorded = get(parameters, "physics_member", nothing)
    recorded == physics_member && return nothing
    isnothing(physics_member) &&
        error("$(what) was fitted on, or built from, tiles of GEMB physics member $(repr(recorded)); " *
              "name the member (physics_member = $(repr(recorded))) so its outputs are kept apart")
    error("$(what) was fitted on, or built from, tiles of GEMB physics member $(repr(recorded)), " *
          "not $(repr(physics_member)); members must not be mixed")
end
