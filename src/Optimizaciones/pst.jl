using JuMP
using Cbc
using PyCall
pfnet = pyimport("pfnet")

# Se obtiene la red
filenames = ["src/Casos/pglib_opf_case30_ieee.m" ,"src/Casos/pglib_opf_case30_ieee__api.m"]
parser = pfnet.PyParserMAT()
nets = [parser.parse(filename) for filename in filenames]

# Pre proceso del caso
for (s, net) in enumerate(nets)
    for br in net.branches
        br.set_as_phase_shifter()
    end
end
δ = Dict((br.index, s) => 0.0 for (s, net) in enumerate(nets) for br in net.branches if br.is_phase_shifter())

# formulacion del problema
optimizer = GLPK.Optimizer
model = Model(with_optimizer(optimizer))

# Vars
α = Dict(br.index => @variable(model,
                               binary=true) for br in nets[1].branches)
β = Dict(br.index => @variable(model) for br in nets[1].branches)
θ = Dict((bus.index, s) => @variable(model) for (s, net) in enumerate(nets)
                                            for bus in net.buses)
pg = Dict((gen.index, s) => @variable(model) for (s, net) in enumerate(nets)
                                             for gen in net.generators if gen.is_slack())
f = Dict((br.index, s) => @variable(model) for (s, net) in enumerate(nets)
                                           for br in net.branches)
f_ = Dict((br.index, s) => @variable(model,
                                     lower_bound = 0.) for (s, net) in enumerate(nets)
                                                       for br in net.branches if br.ratingA>0)
Φ = Dict((br.index, s) => @variable(model) for (s, net) in enumerate(nets)
                                           for br in net.branches if br.is_phase_shifter())
γ = Dict((br.index, s) => @variable(model) for (s, net) in enumerate(nets)
                                           for br in net.branches if br.is_phase_shifter())

# Objective
A, B, C = 1e1, 1e0, 1e3

@objective(model, Min, sum(A*α[br.index]+B*β[br.index] for br in nets[1].branches if br.is_phase_shifter())+
                       sum(C*f_[br.index, s] for (s, net) in enumerate(nets) for br in net.branches if br.ratingA > 0))

# Restricciones
for (s, net) in enumerate(nets)
    for bus in net.buses
        dp = 0
        for gen in bus.generators
            dp += gen.is_slack() ? pg[gen.index, s] : gen.P
        end
        for load in bus.loads
            dp -= load.P
        end
        for br in bus.branches_k
            dp -= f[br.index, s]
        end
        for br in bus.branches_m
            dp += f[br.index, s]
        end
        @constraint(model, dp == 0.)
    end
end

for (s, net) in enumerate(nets)
    for br in net.branches
        ckt = br.index
        k, m = br.bus_k.index, br.bus_m.index

        if br.is_phase_shifter()
            M = 1e1
            @constraint(model, β[ckt] <= α[ckt]*60*float(π)/180)
            @constraint(model, β[ckt] >= α[ckt]*1*float(π)/180)
            @constraint(model, f[ckt, s] == -br.b*(θ[k, s]-θ[m, s]-Φ[ckt, s])-γ[ckt, s])
            @constraint(model, γ[ckt, s]-δ[ckt, s]*(θ[k, s]-θ[m, s]-Φ[ckt, s]) <= (1-α[ckt])*M)
            @constraint(model, γ[ckt, s]-δ[ckt, s]*(θ[k, s]-θ[m, s]-Φ[ckt, s]) >= -(1-α[ckt])*M)
            @constraint(model, γ[ckt, s] <= α[ckt]*M)
            @constraint(model, γ[ckt, s] >= -α[ckt]*M)
            @constraint(model, Φ[ckt, s] <= β[ckt])
            @constraint(model, Φ[ckt, s] >= -β[ckt])
        else
            @constraint(model, f[ckt, s] == -br.b*(θ[k, s]-θ[m, s]))
        end

        if br.ratingA > 0
            @constraint(model, f[ckt, s] <= br.ratingA+f_[ckt, s])
            @constraint(model, f[ckt, s] >= -br.ratingA-f_[ckt, s])
        end
    end
end

for (s, net) in enumerate(nets)
    for bus in net.buses
        if bus.is_slack()
            @constraint(model, θ[bus.index, s] == 0.)
        end
    end
end

optimize!(model)

for (s, net) in enumerate(nets)

    println("Escenario $s: \n")

    for br in net.branches
        ckt = br.index
        carga = round((abs(value(f[ckt, s])))/br.ratingA*100 ,sigdigits=2)
        angle_pst = round(value(Φ[ckt, s])*180/float(π),sigdigits=2)
        pst = value(α[ckt]) > 0 ? " -> PST instalado (ϕ = $angle_pst °)" : ""

        bus_k = br.bus_k.number
        bus_m = br.bus_m.number

        str = "Linea $bus_k - $bus_m [$ckt] : "

        if carga<60
            printstyled(str), printstyled("$carga %", color=:green), printstyled(pst, color=:magenta)
            println()
        elseif 60<=carga<90
            printstyled(str), printstyled("$carga %", color=:yellow), printstyled(pst, color=:magenta)
            println()
        else carga>=90
            printstyled(str), printstyled("$carga %", color=:red), printstyled(pst, color=:magenta)
            println()
        end
    end
end
