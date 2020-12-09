using JuMP
using Ipopt
using PyCall
pfnet = pyimport("pfnet")

# Maximiza el flujo de potencia de una linea

# Obtencion de la red electrica
filename = "src/Casos/pglib_opf_case30_ieee.m"
parser = pfnet.PyParserMAT()
net = parser.parse(filename)

# Indices
bus_indices = [bus.index for bus in net.buses]
gen_indices = [gen.index for gen in net.generators]
br_indices = [br.index for br in net.branches]

# Creacion del modelo
nlp_optimizer = Ipopt.Optimizer
m = Model(with_optimizer(nlp_optimizer, print_level=0))

# Variables
vm = Dict(bus.index => @variable(m, start = 1.0,
                                    lower_bound=bus.v_min,
                                    upper_bound=bus.v_max) for bus in net.buses)
va = Dict(bus.index => @variable(m) for bus in net.buses)
pg = Dict(gen.index => @variable(m, lower_bound=gen.P_min,
                                    upper_bound=gen.P_max) for gen in net.generators)
qg = Dict(gen.index => @variable(m, lower_bound=gen.Q_min,
                                    upper_bound=gen.Q_max) for gen in net.generators)
p_km = Dict(br.index => @variable(m, lower_bound=-br.ratingA,
                                     upper_bound=br.ratingA) for br in net.branches)
q_km = Dict(br.index => @variable(m, lower_bound=-br.ratingA,
                                     upper_bound=br.ratingA) for br in net.branches)
p_mk = Dict(br.index => @variable(m, lower_bound=-br.ratingA,
                                     upper_bound=br.ratingA) for br in net.branches)
q_mk = Dict(br.index => @variable(m, lower_bound=-br.ratingA,
                                     upper_bound=br.ratingA) for br in net.branches)

# Objectivo -> Flujo Maximo en Linea 39
@objective(m, Max, p_km[39]^2 + q_km[39]^2)

# Restricciones
for bus in net.buses # KCL
    dp = 0.
    dq = 0.
    for gen in bus.generators
        dp += pg[gen.index]
        dq += qg[gen.index]
    end
    for load in bus.loads
        dp -= load.P
        dq -= load.Q
    end
    for sh in bus.shunts
        dp -= sh.g * vm[sh.bus.index]^2
        dq += sh.b * vm[sh.bus.index]^2
    end
    for br in bus.branches_k
        dp -= p_km[br.index]
        dq -= q_km[br.index]
    end
    for br in bus.branches_m
        dp -= p_mk[br.index]
        dq -= q_mk[br.index]
    end
    @constraint(m, dp == 0.)
    @constraint(m, dq == 0.)
end

for br in net.branches # KVL
    v_k, v_m = vm[br.bus_k.index], vm[br.bus_m.index]
    θ_k, θ_m = va[br.bus_k.index], va[br.bus_m.index]
    ckt = br.index

    b, g = br.b, br.g
    b_k, g_k = br.b_k, br.g_k
    b_m, g_m = br.b_m, br.g_m
    t = br.ratio
    ϕ = br.phase

    @NLconstraint(m, p_km[ckt] ==  (g_k+g)*t^2*v_k^2 - t*v_k*v_m*(g*cos(θ_k-θ_m-ϕ)+b*sin(θ_k-θ_m-ϕ)))
    @NLconstraint(m, q_km[ckt] == -(b_k+b)*t^2*v_k^2 - t*v_k*v_m*(g*sin(θ_k-θ_m-ϕ)-b*cos(θ_k-θ_m-ϕ)))
    @NLconstraint(m, p_mk[ckt] ==  (g_m+g)*v_m^2 - t*v_k*v_m*(g*cos(θ_k-θ_m-ϕ)-b*sin(θ_k-θ_m-ϕ)))
    @NLconstraint(m, q_mk[ckt] == -(b_m+b)*v_m^2 + t*v_k*v_m*(g*sin(θ_k-θ_m-ϕ)+b*cos(θ_k-θ_m-ϕ)))

    if br.ratingA > 0
        @NLconstraint(m, p_km[ckt]^2 + q_km[ckt]^2 <= br.ratingA^2)
        @NLconstraint(m, p_mk[ckt]^2 + q_mk[ckt]^2 <= br.ratingA^2)
    end

    @constraint(m, θ_k-θ_m <= 30*float(π)/180)
    @constraint(m, θ_k-θ_m >= -30*float(π)/180)
end

for bus in net.buses # Restricciones Barras Slacks
    if bus.is_slack()
        @constraint(m, va[bus.index] == bus.v_ang)
    end
end

optimize!(m)

printstyled("Ramas:", color=:underline)
println()
for br in net.branches
    ckt = br.index
    p_fr, p_to = value(p_km[ckt]), value(p_mk[ckt])
    q_fr, q_to = value(q_km[ckt]), value(q_mk[ckt])
    carga = round(sqrt(max(p_fr^2+q_fr^2, p_to^2+q_to^2)/br.ratingA^2)*100, sigdigits=2)
    bus_k = br.bus_k.number
    bus_m = br.bus_m.number
    str = "Linea $bus_k - $bus_m [$ckt] : "
    if carga<60
        printstyled(str), printstyled("$carga %", color=:green)
        println()
    elseif 60<=carga<90
        printstyled(str), printstyled("$carga %", color=:yellow)
        println()
    else carga>=90
        printstyled(str), printstyled("$carga %", color=:red)
        println()
    end
end