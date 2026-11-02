# LEAD-chi POST HOC (2026-10-07): second, independent QSM engine (QSM.jl v0.5.4, Kames; Julia 1.10.9).
# Chain (QSM.jl documented default): Laplacian phase unwrapping (mgpcg, Dirichlet on Mask 3) -> magnitude-weighted
# multi-echo linear fit with phase offset -> field in ppm -> V-SHARP (QSM.jl defaults) -> RTS (defaults) and TV
# (defaults; NLTV diverged to NaN on these data). Same coded inputs and the same Mask 3 as the primary run; nothing
# from QSMxT/QSM.rs is used. Written as a top-level loop (a single large function stalled in type inference).
# Usage: julia e2_qsmjl.jl <sessions.tsv> <worker> <nworkers> <workdir>
using QSM, NIfTI, Printf
const GAMMA = 267.52          # 1e6 rad s^-1 T^-1 -> field in ppm = omega / (GAMMA * B0)
println(stderr, "start threads=", Threads.nthreads()); flush(stderr)
write_raw(path, A) = open(io -> write(io, A), path, "w")
logstep(lab, x, t0) = (println(stderr, "  ", lab, " ", x, " ", round(time() - t0, digits = 1), "s"); flush(stderr))

tsv, wk, nw, work = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]), ARGS[4]
lines = filter(!isempty, readlines(tsv))
for (i, l) in enumerate(lines)
    (i - 1) % nw == wk - 1 || continue
    f = split(l, '\t')
    lab = String(f[3])
    out = joinpath(work, lab)
    isfile(joinpath(out, "done")) && continue
    mkpath(out)
    t0 = time()
    try
        B0 = parse(Float64, f[4]); TEs = parse.(Float64, f[5:8])
        bdir = (parse(Float64, f[9]), parse(Float64, f[10]), parse(Float64, f[11]))
        src, maskp = String(f[12]), String(f[13])
        p1 = niread(joinpath(src, "$(lab)_echo-1_part-phase_MEGRE.nii.gz"))
        sz = size(p1.raw)[1:3]
        vsz = (Float64(p1.header.pixdim[2]), Float64(p1.header.pixdim[3]), Float64(p1.header.pixdim[4]))
        phas = Array{Float32}(undef, sz..., 4); mag = Array{Float32}(undef, sz..., 4)
        for e in 1:4
            phas[:, :, :, e] .= Float32.(niread(joinpath(src, "$(lab)_echo-$(e)_part-phase_MEGRE.nii.gz")).raw)
            mag[:, :, :, e] .= Float32.(niread(joinpath(src, "$(lab)_echo-$(e)_part-mag_MEGRE.nii.gz")).raw)
        end
        mask = Array{Bool}(niread(maskp).raw .> 0.5)
        logstep(lab, "loaded", t0)
        uphas = unwrap_laplacian(phas, mask, vsz)
        fit = multi_echo_linear_fit(uphas, TEs; W = mag, mask = mask, phase_offset = true)
        B = Float32.(fit[1] ./ (GAMMA * B0)); B[.!mask] .= 0f0
        logstep(lab, "fit", t0)
        fl, mask2 = vsharp(B, mask, vsz); logstep(lab, "vsharp", t0)
        x_rts = rts(fl, mask2, vsz; bdir = bdir); logstep(lab, "rts", t0)
        x_tv = tv(fl, mask2, vsz; bdir = bdir); logstep(lab, "tv", t0)
        write_raw(joinpath(out, "lf.f32"), Float32.(fl))
        write_raw(joinpath(out, "mask.u8"), UInt8.(mask2))
        write_raw(joinpath(out, "chi_rts.f32"), Float32.(x_rts))
        write_raw(joinpath(out, "chi_tv.f32"), Float32.(x_tv))
        open(joinpath(out, "meta.txt"), "w") do io
            @printf(io, "%d %d %d %.6f %.6f %.6f\n", sz..., vsz...)
        end
        touch(joinpath(out, "done"))
        @printf("%s ses-%s ok %.1fs\n", f[1], f[2], time() - t0)
    catch err
        @printf("%s ses-%s error: %s\n", f[1], f[2], sprint(showerror, err))
    end
    flush(stdout)
end
