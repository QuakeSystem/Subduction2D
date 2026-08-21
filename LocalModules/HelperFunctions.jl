# PREPARE VISUALIZATION SETTINGS
function prepare_visualisation(ni; do_vtk=true, pictures = true,  pvd_name = "Subduction2D", data_folder = "Subduction2D_SZU2019/data", JRversion = "Subduction2D_JRv0.6.1", version=nothing, save_particle_points = false, vtk_every = 25,  particle_vtk_every = 25, picture_every = 25)
    figdir   = joinpath(data_folder, JRversion, version)
    if do_vtk == true
        vtk_dir = joinpath(figdir, "vtk")
        if isfile(joinpath(vtk_dir, "$pvd_name.pvd"))
            rm(joinpath(vtk_dir, "$pvd_name.pvd"))
        end
        take(vtk_dir)
        checkpoint = joinpath(figdir, "checkpoint")
        take(checkpoint)
    end
    vis=(;do_vtk,vtk_dir,pvd_name ,figdir,save_particle_points,vtk_every,particle_vtk_every,pictures,picture_every,checkpoint,Vx_v = @zeros(ni .+ 1...), Vy_v = @zeros(ni .+ 1...),)

    return vis
end