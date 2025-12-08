function Base.show(io::IO, M::JustRelax.Mask)
    inner = M.mask
    sz = size(inner)

    println(io, "$(sz[1])×$(sz[2]) Mask:")

    for i in 1:sz[1]
        for j in 1:sz[2]
            print(io, inner[i, j] != 0 ? "1 " : "0 ")
        end
        println(io)
    end
end


using JustRelax, ParallelStencil
nx, ny = 4, 10

x = rand(nx, ny)
m = x .< 0.5
jra = JustRelax.Mask(nx, ny, 2:3, 2:3)

basevisc = 1e25
viscoblock = jra.mask .* basevisc


# "ConstantDirichletBoundaryCondition" begin
#             ni = 10, 10
#             A = rand(ni...)
#             value = 5.0e0
#             mask = JustRelax.Mask(ni..., 4:7, 4:7)

#             bc = JustRelax.ConstantDirichletBoundaryCondition(value, mask)
