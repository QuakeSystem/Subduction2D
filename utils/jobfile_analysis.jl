using Printf

"""
Extract basic statistics of the solver from a job.*.*out file.

Finds:
it    itPH  final_err     total_iters  stokes_time(s)  thermal_time(s)
per iteration and stats over all.
"""
function analyse_job(directory=pwd())

    # Find .out files - matches job.<node>.<digits>.out OR job.<node>.<digits>.<version>.out
    out_files = filter(f -> match(r"^job\..+\.out$", basename(f)) !== nothing,
                   readdir(directory, join=true))

    isempty(out_files) && return

    for out_file in out_files
        lines = readlines(out_file)

        # --- Extract version name ---
        # Try from filename first: job.<node>.<digits>.<version>.out
        version_name = let m = match(r"^job\.[^.]+\.\d+\.(.+)\.out$", basename(out_file))
            m !== nothing ? m.captures[1] : nothing
        end

        # Fall back to searching log content
        if version_name === nothing
            for line in lines
                m = match(r"^version is (.+)$", line)
                if m !== nothing
                    version_name = strip(m.captures[1])
                    break
                end
            end
        end

        version_name === nothing && continue

        # --- Parse timesteps ---
        struct_timestep = @NamedTuple{
            it::Int,
            n_phases::Int,
            final_err::Float64,
            stokes_time::Float64,
            thermal_time::Float64,
            total_iters::Int
        }

        timesteps = struct_timestep[]

        current_it   = nothing
        itPH_errors  = Float64[]
        stokes_time  = NaN
        thermal_time = NaN
        total_iters  = 0

        function flush_timestep!()
            current_it === nothing && return
            push!(timesteps, (
                it           = current_it,
                n_phases     = length(itPH_errors),
                final_err    = isempty(itPH_errors) ? NaN : last(itPH_errors),
                stokes_time  = stokes_time,
                thermal_time = thermal_time,
                total_iters  = total_iters,
            ))
        end

        for line in lines
            m = match(r"it \+= 1 = (\d+)", line)
            if m !== nothing
                flush_timestep!()
                current_it   = parse(Int, m.captures[1])
                itPH_errors  = Float64[]
                stokes_time  = NaN
                thermal_time = NaN
                total_iters  = 0
                continue
            end

            m = match(r"itPH = \d+ iter = \s*(\d+) iter/nx.*err = ([0-9eE.+\-]+)", line)
            if m !== nothing
                total_iters = parse(Int, m.captures[1])
                push!(itPH_errors, parse(Float64, m.captures[2]))
                continue
            end

            m = match(r"Total time:\s+([0-9eE.+\-]+) s", line)
            if m !== nothing
                stokes_time = parse(Float64, m.captures[1])
                continue
            end

            m = match(r"solver finished in ([0-9eE.+\-]+) seconds", line)
            if m !== nothing
                thermal_time = parse(Float64, m.captures[1])
                continue
            end
        end
        flush_timestep!()

        isempty(timesteps) && continue

        # --- Write analysis file ---
        outname = joinpath(directory, "analysis_$(version_name).txt")
        open(outname, "w") do f
            println(f, "Analysis of: $(basename(out_file))")
            println(f, "Version:     $version_name")
            println(f, "Timesteps completed: $(length(timesteps))")
            println(f, "")
            println(f, "=" ^ 90)
            println(f, rpad("it", 6),
                       rpad("itPH", 6),
                       rpad("final_err", 14),
                       rpad("total_iters", 13),
                       rpad("stokes_time(s)", 16),
                       rpad("thermal_time(s)", 16))
            println(f, "=" ^ 90)
            for ts in timesteps
                println(f,
                    rpad(ts.it, 6),
                    rpad(ts.n_phases, 6),
                    rpad(@sprintf("%.3e", ts.final_err), 14),
                    rpad(ts.total_iters, 13),
                    rpad(@sprintf("%.3f", ts.stokes_time), 16),
                    rpad(@sprintf("%.4f", ts.thermal_time), 16),
                )
            end
            println(f, "=" ^ 90)
            println(f, "")

            # --- Summary stats ---
            stokes_times  = filter(!isnan, [ts.stokes_time  for ts in timesteps])
            thermal_times = filter(!isnan, [ts.thermal_time for ts in timesteps])
            final_errs    = filter(!isnan, [ts.final_err    for ts in timesteps])
            n_phases      = [ts.n_phases for ts in timesteps]

            println(f, "Summary")
            println(f, "-" ^ 40)
            println(f, "Total timesteps:         $(length(timesteps))")
            isempty(stokes_times)  || println(f, "Stokes time  - mean: $(@sprintf("%.2f", sum(stokes_times)/length(stokes_times))) s  |  min: $(@sprintf("%.2f", minimum(stokes_times))) s  |  max: $(@sprintf("%.2f", maximum(stokes_times))) s")
            isempty(thermal_times) || println(f, "Thermal time - mean: $(@sprintf("%.4f", sum(thermal_times)/length(thermal_times))) s  |  min: $(@sprintf("%.4f", minimum(thermal_times))) s  |  max: $(@sprintf("%.4f", maximum(thermal_times))) s")
            isempty(final_errs)    || println(f, "Final err    - mean: $(@sprintf("%.3e", sum(final_errs)/length(final_errs)))  |  min: $(@sprintf("%.3e", minimum(final_errs)))  |  max: $(@sprintf("%.3e", maximum(final_errs)))")
            isempty(n_phases)      || println(f, "itPH count   - mean: $(@sprintf("%.1f", sum(n_phases)/length(n_phases)))  |  min: $(minimum(n_phases))  |  max: $(maximum(n_phases))")
        end
    end
end


# Run 'jobfile_analysis.jl' on multiple outputs.
# Run from a working directory containing multiple version output folders
# Working directory
#  |
# --- Working directory/v1
# --- Working directory/v2
# --- Working directory/v3
parent = pwd()
for dir in readdir(parent, join=true)
    println(dir)
    isdir(dir) && analyse_job(dir)
end