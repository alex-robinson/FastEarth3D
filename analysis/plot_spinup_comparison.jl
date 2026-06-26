#!/usr/bin/env julia
# Spin-up comparison plots (branch restart-spinup). For each run: a FOUR-panel map
# row at [init after spin-up, -20 ka, -10 ka, 0 ka] showing
#   - bedrock topography (semi-transparent topographic colors),
#   - ice-sheet surface shaded grey->white with 500 m surface contours,
#   - relative-sea-level contours (labelled) over the ocean,
#   - the barystatic sea level (BSL) printed in the top-right corner.
# Plus one figure of the BSL time series for every run on the same axes.
#
# Usage:  julia analysis/plot_spinup_comparison.jl [label=out.nc ...]
# With no args it uses the RUNS table below (the l64 1-D vs 3-D spin-up experiment
# plus the cross-resolution l128-from-l64 run).

using NCDatasets
using CairoMakie
using Printf

# label => output netCDF. Override on the command line as label=path.
const DEFAULT_RUNS = [
    "3-D spin-up (l64)"        => "runs/cmp_spin3d_l64/out.nc",
    "1-D spin-up (l64)"        => "runs/cmp_spin1d_l64/out.nc",
    "l128 from l64 spin-up"    => "runs/xres_l128_from_l64/out.nc",
]

# Panel 1 is always the first written slice ("init after spin-up"); the rest are the
# nearest slices to these times [ka].
const TIMES_KA = [-20.0, -10.0, 0.0]
const OUTDIR   = "analysis/figs"

const TOPO_CMAP = cgrad([:navy, :dodgerblue, :paleturquoise,
                         :darkgreen, :khaki, :saddlebrown, :white],
                        [0.0, 0.30, 0.48, 0.52, 0.68, 0.88, 1.0])
const TOPO_RANGE = (-6000.0, 6000.0)
const ICE_CMAP   = cgrad([:gray35, :gray70, :white])
const ICE_RANGE  = (0.0, 4000.0)

parse_runs(args) = isempty(args) ? DEFAULT_RUNS :
    [ (p = split(a, "="; limit = 2); String(p[1]) => String(p[2])) for a in args ]

"Roll lon to [-180,180), sort lon and lat ascending; return (lon, lat, perm_lon, perm_lat)."
function grid_axes(ds)
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    lon = map(x -> x > 180 ? x - 360 : x, lon)
    pl = sortperm(lon); pt = sortperm(lat)
    return lon[pl], lat[pt], pl, pt
end

field(ds, name, it, pl, pt) = Float64.(ds[name][:, :, it])[pl, pt]
nearest_index(times_yr, t_ka) = argmin(abs.(times_yr ./ 1000 .- t_ka))

function panel!(ax, ds, it, pl, pt, lon, lat)
    zb  = field(ds, "z_bed",   it, pl, pt)
    hi  = field(ds, "h_ice",   it, pl, pt)
    C   = field(ds, "C_ocean", it, pl, pt)
    surf = zb .+ hi
    has_ice = hi .> 1.0

    heatmap!(ax, lon, lat, zb; colormap = TOPO_CMAP, colorrange = TOPO_RANGE, alpha = 0.65)
    ice_s = map((s, m) -> m ? s : NaN, surf, has_ice)
    heatmap!(ax, lon, lat, ice_s; colormap = ICE_CMAP, colorrange = ICE_RANGE)
    contour!(ax, lon, lat, ice_s; levels = 500:500:4500, color = (:black, 0.45), linewidth = 0.5)

    is_sea = (C .> 0.5) .| (zb .< 0.0)
    rsl = field(ds, "rsl", it, pl, pt)
    rsl_sea = map((r, m) -> m ? r : NaN, rsl, is_sea)
    contour!(ax, lon, lat, rsl_sea; levels = -160:40:160, labels = true,
             labelsize = 8, color = :firebrick, linewidth = 0.8)

    bsl = Float64(ds["bsl"][it])
    text!(ax, 0.98, 0.97; text = @sprintf("BSL = %.0f m", bsl),
          space = :relative, align = (:right, :top), fontsize = 11,
          font = :bold, color = :black, strokecolor = :white, strokewidth = 2)
    return nothing
end

function plot_run(label, path)
    ds = NCDataset(path)
    times = Float64.(ds["time"][:])
    lon, lat, pl, pt = grid_axes(ds)

    # panel index list: first slice (init after spin-up) + nearest to TIMES_KA.
    its    = vcat(1, [nearest_index(times, t) for t in TIMES_KA])
    titles = vcat(@sprintf("init  (t = %.0f yr)", times[1]),
                  [@sprintf("%.0f ka  (t = %.0f yr)", t, times[nearest_index(times, t)]) for t in TIMES_KA])

    fig = Figure(size = (1960, 520))
    Label(fig[0, 1:4], label; fontsize = 18, font = :bold)
    for (j, it) in enumerate(its)
        ax = Axis(fig[1, j]; aspect = DataAspect(), title = titles[j],
                  xlabel = "lon", ylabel = j == 1 ? "lat" : "")
        limits!(ax, -180, 180, -90, 90)
        panel!(ax, ds, it, pl, pt, lon, lat)
    end
    Colorbar(fig[1, 5], colormap = TOPO_CMAP, limits = TOPO_RANGE, label = "bedrock elevation [m]")
    close(ds)
    fn = joinpath(OUTDIR, "maps_" * replace(label, r"[^A-Za-z0-9]" => "_") * ".png")
    save(fn, fig); println("wrote ", fn)
end

function plot_bsl(runs)
    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xlabel = "time [ka]", ylabel = "barystatic sea level [m]",
              title = "Barystatic sea level vs present-day reference")
    for (label, path) in runs
        isfile(path) || continue
        ds = NCDataset(path)
        t = Float64.(ds["time"][:]) ./ 1000
        b = Float64.(ds["bsl"][:])
        lines!(ax, t, b; label = label, linewidth = 2)
        close(ds)
    end
    axislegend(ax; position = :rb)
    fn = joinpath(OUTDIR, "bsl_timeseries.png")
    save(fn, fig); println("wrote ", fn)
end

function main()
    runs = parse_runs(ARGS)
    mkpath(OUTDIR)
    for (label, path) in runs
        isfile(path) || (@warn "missing output, skipping" label path; continue)
        plot_run(label, path)
    end
    plot_bsl(runs)
end

main()
