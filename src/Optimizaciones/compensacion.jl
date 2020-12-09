using JuMP
using Ipopt
using Juniper
using Cbc
using PyCall
pfnet = pyimport("pfnet")

# Ubicacion optima de capacitores/reactores para mejorar el perfil de tensiones

# Obtencion de la red electrica
filename = "src/Casos/pglib_opf_case30_ieee.m"
parser = pfnet.PyParserMAT()
net = parser.parse(filename)

# Limites de tension ajustados
for bus in net.buses
    bus.v_max = 1.01
    bus.v_min = 0.99
end

# Indices
bus_indices = [bus.index for bus in net.buses]
gen_indices = [gen.index for gen in net.generators]
br_indices = [br.index for br in net.branches]

# Create model
optimizer = Juniper.Optimizer
nl_solver = with_optimizer(Ipopt.Optimizer, print_level=0)
mip_solver = with_optimizer(Cbc.Optimizer, logLevel=0)
m = Model(with_optimizer(optimizer,
                         nl_solver=nl_solver,
                         mip_solver=mip_solver))

# Vars
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
c_sh = Dict(bus.index => @variable(m, start = 0.,
                                     lower_bound = 0.,
                                     upper_bound = 0.3) for bus in net.buses)
c_sh_bin = Dict(bus.index => @variable(m, binary=true) for bus in net.buses)
r_sh = Dict(bus.index => @variable(m, start=0.,
                                     lower_bound = 0.,
                                     upper_bound = 0.3) for bus in net.buses)
r_sh_bin = Dict(bus.index => @variable(m, binary=true) for bus in net.buses)

# Objective
M = 1
@objective(m, Min, sum(c_sh[bus.index]+r_sh[bus.index] for bus in net.buses)+
                   M * sum(c_sh_bin[bus.index]+r_sh_bin[bus.index] for bus in net.buses))
# Constraints
for bus in net.buses # KCL
    i = bus.index
    dp = 0.
    dq = c_sh[i]*c_sh_bin[i] - r_sh[i]*r_sh_bin[i]
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

for bus in net.buses # Voltage constraints
    if bus.is_slack()
        @constraint(m, va[bus.index] == bus.v_ang)
    end
end

optimize!(m)


for bus in net.buses
    i = bus.index
    k = bus.number
    v_mag = round(value(vm[i]), sigdigits=3)
    v_ang = round(value(va[i]) *180/float(π), sigdigits=2)
    v_max, v_min = bus.v_max, bus.v_min
    Q_sh = round(value(c_sh[i] - r_sh[i]), sigdigits=3)

    println("Barra $k: $v_mag | $v_ang ° | [$v_min:$v_max] | Q_sh: $Q_sh p.u.")
end

printstyled("Barras:", color=:underline)
println()
for bus in net.buses
    i = bus.index
    k = bus.number
    v_mag = round(value(vm[i]), sigdigits=3)
    Q_sh = round(value(c_sh[i] - r_sh[i]), sigdigits=3)
    reactor = abs(value(c_sh[i] - r_sh[i])) >= 1e-4 ?  "reactor: $(round(value(c_sh[i] - r_sh[i]), sigdigits=3))" : ""
    str = "Barra $(bus.number) :"
    if abs(1-v_mag) < 0.999
        printstyled(str), printstyled("$v_mag %  ", color=:yellow), 
        printstyled(reactor, color=:magenta)
        println()
    elseif 0.999 <= abs(1-v_mag) < 1.001
        printstyled(str), printstyled("$v_mag %  ", color=:green),
        printstyled(reactor, color=:magenta)
        println()
    else
        printstyled(str), printstyled("$v_mag %  ", color=:yellow),
        printstyled(reactor, color=:magenta)
        println()
    end
end