using LinearAlgebra
using DataFrames
using CSV
using Plots
using SparseArrays

# Calcular la matriz de admitancia nodal

function calcular_ybus(lines,nodes)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
    Salida :    Ybus : matriz
    """
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    Ybus = zeros(num_nodes, num_nodes)*1im

    for k = 1:num_lines
        # Nodo de envío
        n1 = lines.FROM[k]
        # Nodo de recibo
        n2 = lines.TO[k]
        # Admitancia de la línea
        yL = 1/(lines.R[k]+lines.X[k]*1im)
        # Susceptancia de la línea
        Bs = lines.B[k]*1im/2
        # Valor del TAP
        t = lines.TAP[k]
        if lines.TAP[k] == 0
            Ybus[n1,n1] += yL + Bs   # Dentro de la diagonal
            Ybus[n1,n2] -= yL        # Fuera de la diagonal
            Ybus[n2,n1] -= yL        # Fuera de la diagonal
            Ybus[n2,n2] += yL + Bs   # Dentro de la diagonal
        else
            Ybus[n1,n1] += (t^2 - t)*yL  # Dentro de la diagonal
            Ybus[n1,n2] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n1] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n2] += (1-t)*yL      # Dentro de la diagonal
        end
    end
    return Ybus
end

## Función principal
lines = DataFrame(CSV.File("Actividad 1/lines.csv"))
nodes = DataFrame(CSV.File("Actividad 1/nodes.csv"))
Ybus = calcular_ybus(lines,nodes)
#Ybus = sparse(Ybus)

# Se calcula la matriz B_bus

function B_bus(lines,nodes)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
    Salida :    Bbus : matriz
    """
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    Bbus = zeros(num_nodes, num_nodes)

    for k = 1:num_lines
        # Nodo de envío
        n1 = lines.FROM[k]
        # Nodo de recibo
        n2 = lines.TO[k]
        # Susceptancia
        BL = 1/(lines.X[k])
        Bbus[n1,n1] += BL        # Dentro de la diagonal
        Bbus[n1,n2] -= BL        # Fuera de la diagonal
        Bbus[n2,n1] -= BL        # Fuera de la diagonal
        Bbus[n2,n2] += BL        # Dentro de la diagonal
    end
    
    return Bbus
end

B = B_bus(lines, nodes)
# display(B)

function flujo_dc(lines,nodes,B_bus)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
    Salida :    P_km (P_lineas) : Vector
                θ : Vector
    """
    # Identificando el nodo slack
    s = nodes[nodes.TYPE .== 3, "NUMBER"]
    B_bus = B_bus[setdiff(1:end, s), setdiff(1:end, s)]

    # Calculo de las potencias nodales
    Pn = nodes.PGEN .- nodes.PLOAD

    # Eliminando el nodo slack (Tipo 3)
    Pn = Pn[1:end .!= s]

    # Angulos nodales
    theta_ = inv(B_bus) * Pn

    # Vector θ completo (slack = 0)
    theta = zeros(nrow(nodes))
    theta[1:end .!= s] = theta_

    # Vector de flujos
    flujos = zeros(nrow(lines))

    # Calculando las pontecias de linea
    for i in 1:nrow(lines)
        from = lines.FROM[i]
        to = lines.TO[i]
        B = 1/lines.X[i]
        flujos[i] = B*(theta[from] - theta[to])
    end

    return theta, flujos

end

function flujo_dc_data(lines,nodes,B_bus)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
    Salida :    P_km (P_lineas) : DataFrame
                θ : DataFrame
                
    """

    # Identificando el nodo slack
    s = nodes[nodes.TYPE .== 3, "NUMBER"]
    B_bus = B_bus[setdiff(1:end, s), setdiff(1:end, s)]

    # Calculo de las potencias nodales
    Pn = nodes.PGEN .- nodes.PLOAD

    # Eliminando el nodo slack (Tipo 3)
    Pn = Pn[1:end .!= s]

    # Angulos nodales
    theta_ = inv(B_bus) * Pn

    # Vector θ completo (slack = 0)
    theta = zeros(nrow(nodes))
    theta[1:end .!= s] = theta_


    # Vector de flujos
    flujos = zeros(nrow(lines))

    # Calculando las pontecias de linea
    for i in 1:nrow(lines)
        from = lines.FROM[i]
        to = lines.TO[i]
        B = 1/lines.X[i]
        flujos[i] = B*(theta[from] - theta[to])
    end

    # Creando un DataFrame para flujos y así visualizarlo mejor
    df_flujos = DataFrame(
            Desde = lines.FROM,
            Hacia = lines.TO,
            Flujo_pu = flujos
    )

    # Creando un DataFrame para los angulos nodales y así visualizarlos mejor
    df_theta = DataFrame(
        Nodo = nodes.NUMBER,
        Angulo_rad = theta
    )

    return df_theta,df_flujos

end

# theta_df, flujos_df = flujo_dc(lines, nodes,B)

# # Para ver los resultados:
# println("Ángulos nodales:")
# display(theta_df)
# println("\nFlujos de potencia:")
# display(flujos_df)

# # Para guardar los ángulos
# CSV.write("Actividad 1/angulos_nodales.csv", theta_df)

# # Para guardar los flujos
# CSV.write("Actividad 1/flujos_potencia.csv", flujos_df)

# B_bus = matriz_B_bus(lines,nodes)
# display(B_bus)

function es_invertible(matriz)

    """
    Entradas: matriz: matriz cuadrada a evaluar

    Salidas: true: si la matriz es invertible (det != 0)
             false: si la matriz no es invertible (det aproximadamente igual a 0 o la operación falla)
    """
    try
        # Intenta calcular el determinante
        det_valor = det(matriz)
        # Verifica si el determinante es cercano a cero
        return !isapprox(det_valor, 0, atol=1e-10)
    catch
        return false
    end
end

function contingencias(lines,nodes,B_bus)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
    Salida :    P_km (P_lineas) : DataFrame
                θ : DataFrame
    Nota:       Estas contingencias son n-1
    """
    # Numero de lineas y nodos
    num_lines = nrow(lines)
    num_nodes = nrow(nodes)

    # Creando DataFrames para los angulo y flujos
    df_theta = DataFrame()
    df_flujos = DataFrame()

    for i in 1:num_lines

        # Creando una copia de la B_bus para modificarla caso a caso
        mBbus = copy(B_bus)

        # Nodo de envío
        n1 = lines.FROM[i]
        # Nodo de recibo
        n2 = lines.TO[i]
        # Quitando el aporte de la línea a la matriz de Susceptancia
        # Susceptancia
        BL = 1/(lines.X[i])
        mBbus[n1,n1] -= BL        # Dentro de la diagonal
        mBbus[n1,n2] += BL        # Fuera de la diagonal
        mBbus[n2,n1] += BL        # Fuera de la diagonal
        mBbus[n2,n2] -= BL        # Dentro de la diagonal

        # Quitando el nodo slack para verificar si es matriz singular
        B_temp = copy(mBbus)
        s = nodes[nodes.TYPE .== 3, "NUMBER"]
        B_temp = B_temp[setdiff(1:end, s), setdiff(1:end, s)]

        # Determinando si la matriz B_bus es singular
        if !es_invertible(B_temp)
            # Creando los vectores para angulos y flujos
            theta = zeros(num_nodes)
            flujos = zeros(num_lines)

            # Asignando un valor para este caso especial
            theta[1:end] .= Inf64
            flujos[1:end] .= Inf64
        else
            theta,flujos = flujo_dc(lines, nodes, mBbus)
        end
        # Haciendo cero el flujo de potencia de la línea que se saca de operación
        flujos[i] = 0
        # Determinando si existen lineas en paralelo
        if i != num_lines && (lines.FROM[i+1] == lines.FROM[i]) && (lines.TO[i+1] == lines.TO[i])

            # Añadiendo las lineas en paralelo
            theta_ = " Angulos en rad caso linea $n1 - $n2 (Caso Paralelo)"
            df_theta[!,theta_] = theta
            flujos_ = "Flujos en p.u. caso linea $n1 - $n2 (Caso Paralelo)" 
            df_flujos[!,flujos_] = flujos
        else
            theta_ = " Angulos en rad caso linea $n1 - $n2"
            df_theta[!,theta_] = theta
            flujos_ = "Flujos en p.u. caso linea $n1 - $n2" 
            df_flujos[!,flujos_] = flujos
        end
    end
    return df_theta,df_flujos

end

theta,flujos = contingencias(lines, nodes, B)

CSV.write("Actividad 1/flujos_potencia_contingencias.csv", flujos)
CSV.write("Actividad 1/angulos_nodales_contingencias.csv", theta)

# Graficando los datos para mayor legibilidad
# Leer datos

data_theta = CSV.read("Actividad 1/angulos_nodales_contingencias.csv", DataFrame)
data_flujos = CSV.read("Actividad 1/flujos_potencia_contingencias.csv",DataFrame)


# Convertir el DataFrame a una matriz
matriz_theta = Matrix(data_theta)
matriz_flujos = Matrix(data_flujos)

# Creando el mapa de calor para los angulos
p1 = heatmap(
    matriz_theta,
    c=:plasma,
    aspect_ratio=:equal,
    xlabel="Casos",
    ylabel="Nodo",
    title="Angulos para las contingencias n-1"
)

# Creando el mapa de calor para los flujos
p2 = heatmap(
    matriz_flujos,
    c=:plasma,
    aspect_ratio=:equal,
    xlabel="Casos",
    ylabel="Lineas",
    title="Flujos de potencia para las lineas"
)

# Combinar los plots
plot(p1, p2, layout=(1,2), size=(1200,500))