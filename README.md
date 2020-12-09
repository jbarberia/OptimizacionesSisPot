# OptimizacionesSisPot

Diversas optimizaciones sobre el caso de prueba IEEE30 de [pglib-opf](https://github.com/power-grid-lib/pglib-opf).
Los problemas se formularon a traves de la libreria de optimización [JuMP](https://github.com/jump-dev/JuMP.jl) y la libreria de sistemas de potencia [PFNET.py](https://github.com/ttinoco/PFNET.py).


### Optimizaciones:

- Despacho Economico (OPF)
- Despacho Economico con Lineas Abiertas (OTS)
- Flujo Maximo
- Ubicación Optima de Reactores/Capacitores
- Ubicación Optima de Transformadores Desfasadores (PST)
- PST + Incorporacion de nuevas lineas

### Requisitos:

Librerias:
- [PyCall](https://github.com/JuliaPy/PyCall.jl)
- [JuMP](https://github.com/jump-dev/JuMP.jl)
- [PFNET.py](https://github.com/ttinoco/PFNET.py)

Solvers:
- [GLPK](https://github.com/jump-dev/GLPK.jl)
- [Cbc](https://github.com/jump-dev/Cbc.jl)
- [Ipopt](https://github.com/jump-dev/Ipopt.jl)
- [Juniper](https://github.com/lanl-ansi/Juniper.jl)




### Uso:

Una vez instalado todos los requisitos se puede llamar a cada script de optimización desde el terminal de julia

```
# Correr el OPF
julia >> include("src/Optimizaciones/opf.jl")
```
