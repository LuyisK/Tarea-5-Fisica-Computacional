"""
FÃ­sica Computacional - TS13

Valentina Campos Aguilar
Luis Alfredo Guerrero Camacho 

METRÓPOLIS - MONTE CARLO: Modelo Ising 1-D
"""


import matplotlib.pyplot as plt
import numpy as np

# Se declaran las variables universales.
J = 1
kb = 1 #Valor constante de Boltzmann
B = 0.5 #Valor del campo magnético B
nespines = 100
N = nespines


'''Se crea la función que permite obtener la Energí­a de Ising del sistema de espines determinado donde 
se tiene de parámetros:
    - el arreglo de espines
    - el valor J
La misma retorna el valor de la Energí­a de Ising.'''

def Energía_Ising(arreglo_espines, J):
    valorE = 0
    for i in range(len(arreglo_espines)-1):
        valorE += arreglo_espines[i]* arreglo_espines[i+1]
    return -J *valorE
  
  
''' Se genera una función que permita obtener un arreglo de espines aleatorios con el número de 
espines en el arreglo como parámetro y el arreglo de espines aleatoriosde salida '''

def Espines_Aleatorios(nespines):
    espines = np.random.randint(-1,1,size=nespines)
    for i in range(nespines):
        if espines[i]== 0:
            if np.random.random()>0.5:
                espines[i]= 1
            else:
                espines[i]= -1
    return espines


''' Se crea una función que permita obtener el valor de la magnetización del arreglo de espines que
se introduzca como parámetro de la misma. Esta retorna dicho valor como salida.'''

def Magnetización(arreglo_espines):
    valorM = 0
    for i in arreglo_espines:
        valorM += i
    return valorM 


''' Se crea la función que ejecuta el Modelo de Ising para un arreglo de espines definido y a una 
  temperatura constante, con parámetros:
    - temperatura del sistema
    - pasos a ejecutar en el modelo
Se retorna la lista de valores para la Energías de Ising y magnetización para cada paso ejecutado,
el estado del sistema y el paso donde el sistema alcanza el equilibrio en caso de ser así. '''

def Modelo_Ising(temperatura, npasos):
    
    #Se inicializan las variables a (estado del sistema), b (paso donde se alcanza el equilibrio),
    #y c (variable condicional)
    a = ("No equilibrio")
    b = 1
    c = 0
         
    #Se define el arreglos de espines hacia abajo:
    #arreglo_espines=np.ones([nespines],np.int)*-1

    #Se define el arreglo de espines hacia arriba:
    arreglo_espines=np.ones([nespines],np.int)

    #Se define el arreglo de espines aleatorios:
    #arreglo_espines=CreacionEspines(nespines)

    #Esta lista almacenará el arreglo de espines
    listagráfico=[]
    listagráfico.append(np.array(arreglo_espines))
    
    #Se calculan los valores iniciales de energía y magnetización para el sistema de espines definido
    energía_inicial = Energía_Ising(arreglo_espines, J)
    magnetización_inicial = Magnetización(arreglo_espines)

    #Se crean las listas que almacenarán los valores de energía y magnetización para cada paso
    #ejecutado.
    energía_acum=[]
    energía_acum.append(energía_inicial)
    magnet_acum= []
    magnet_acum.append(magnetización_inicial)
  
    for k in range(npasos):

        #Se calcula la energía del estado i
        ienergía=Energía_Ising(arreglo_espines, J)
        #Se escoge aleatoriamente un espin
        iespin=np.random.randint(nespines)
        #Se cambia el estado del espin
        arreglo_espines[iespin]*=-1
        #Se calcula la energía del estado j
        jenergía=Energía_Ising(arreglo_espines, J)
        #Se calcula la diferencia de energías entre el estado i y j 
        deltaE=jenergía-ienergía
    
        #Se define la probabilidad de acpetación al cambio de estado
        pacept=np.exp(-deltaE/(kb*temperatura))

        #Se comprueba si se realiza el cambio:
        if deltaE>0:
            if np.random.random()<pacept:
                pass
                #Se acepta el cambio de estado del epsin
            else:
                arreglo_espines[iespin]*=-1 #aqui se rechazaria el cambio
                #Se rechaza el cambio de estado del espin
        else:
            pass
        
        #Se calcula la energía y magnetización del estado actual del sistema y se almanecenan los 
        #resultados en la lista correspondiente
        energía_actual = Energía_Ising(arreglo_espines, J)
        energía_acum.append(energía_actual)
        
        magnetización_actual = Magnetización(arreglo_espines)
        magnet_acum.append(magnetización_actual)
        
        #Se actualiza la lista con el arreglo de espines del sistema actual
        listagráfico.append(np.array(arreglo_espines))
    
    #Se genera un método de verificación para determinar si el sistema está en equilibrio
    for i in range(npasos):
        if i % 500 == 0 and c == 0:
            #Se calcula la diferencia de energía cada 500 pasos ejecutados en el modelo 
            delta_eq = abs(energía_acum[i]-energía_acum[i-500])
            #Se define un umbral de acpetación para la diferencia definida anteriormente
            umbral = 5
            #Si la diferencia de energética es menor al umbral establecido entonces se actualizan el
            #estado del sistema, la variable condicional que se actualiza una vez que el equilibrio
            #se alcanza y el paso en el que este equilibrio es alcanzado
            if umbral > delta_eq:
                a = ("Equilibrio")
                b = i
                c += 1
        else:
            pass

    return listagráfico, energía_acum, magnet_acum, a, b


''' Se crea la función que permite obtener la energía y magnetización promedio cuando el Modelo de 
Ising se corre un número establecido de recorridos a una temperatura constante. Los parámetros de la 
función son:
    - el número de recorridos
    - la temperatura del modelo 
    - el número de pasos que ejecuta el modelo
Se obtiene como salida las listas de energía y magnetización promedio obtenidas para cada recorrido'''

def Valores_Promedio(recorridos, temperatura, npasos):
    #Se inicializan las matrices de energía y magnetización.
    matriz_estados_energía = []
    matriz_estados_magnetización = []
    for i in range(recorridos):
        #Se llama al Modelo de Ising para ejecutarse la cantidad de recorridos establecidos.
        listagraf, energiaacum, magnetacum, estado, paso_equilibrio = Modelo_Ising(temperatura, npasos)
        ienergía = energiaacum
        imagnetización = magnetacum
        #Se actualizan las matrices de energía y magnetización con los valores obtenidos de
        #estas variables para cada uno de los recorridos.
        matriz_estados_energía.append(ienergía)
        matriz_estados_magnetización.append(imagnetización)
    #Se trasponen ambas matrices para calcular de mejor forma el promedio energías y magnetización
    #por cada recorrido realizado. 
    matriz_estados_energía = np.array(matriz_estados_energía).T
    matriz_estados_magnetización = np.array(matriz_estados_magnetización).T
    
    #Se inicializan las listas de energía y magnetización promedio.
    energías_promedio = []
    magnetizaciones_promedio = []
    
    #Se calculan los promedios de ambas variables cada npasos y se almacenen en las listas definidas
    for i in range(npasos):
        ienergía_prom = np.mean(matriz_estados_energía[i])
        imagnetización_prom = np.mean(matriz_estados_magnetización[i])
        energías_promedio.append(ienergía_prom)
        magnetizaciones_promedio.append(imagnetización_prom)
    
    return energías_promedio, magnetizaciones_promedio


''' Se crea la función que prmite calcular la energía interna 2 (U2) de un sistema en M recorridos a una 
temperatura definida donde los parámetros son:
    - temperatura del sistema
    - número de pasos a ejecutar en el Modelo de Ising
Se obtiene como salida la lista de energías internas obtenidas por cada recorrido.'''

def Energía_Interna_2(temperatura, npasos):
    #Se inicializa la lista de energía interna 2 y las variables de energía Et y M
    energíaU2 = []
    M = 10
    Et = 0
    for i in range(M):
        #La variable Et se define como el promedio de la lisa de energías obtenida en un recorrido
        #del Modelo de Ising. 
        Et = np.mean(Modelo_Ising(temperatura, npasos)[1])
        #Se actualiza la lista de U2 con los valores de Et al cuadrado obtenidos
        energíaU2.append(Et**2)
    #Para obtener el valor de U2 se promedian los resultados de la lista obtenida anteriormente
    U2 = np.mean(energíaU2)
    return U2


'''Se define la función que permite obtener los cálculos analíticos de la energía interna, la 
magnetización y el calor específico con parámetros:
    - valor J
    - la constante de Boltzmann
    - el número de espines en el sistema
    - la temperatura del sistema
Se obtienen como salidas las tres variables analíticas.'''

def Cálculos_Analíticos(J, kb, N, temperaturas):
    U_An = -N*J*np.tanh(J/(kb*temperaturas))
    M_An = (N*np.exp(J/(kb*temperaturas))*np.sinh(B/(kb*temperaturas))) / np.sqrt(np.exp((2*J)/(kb*temperaturas)) * (np.sinh(B/(kb*temperaturas)))**2 + np.exp((-2*J)/(kb*temperaturas)))
    C_An = (J/(kb*temperaturas))**2 /(np.cosh(J/(kb*temperaturas)))**2
    
    return U_An, M_An, C_An
    
''' Primero se presentan los resultados del Modelo de Ising 1-D diseñado para uno y 10 recorridos'''

#Se definen los parámetros con los que se correrá el Modelo de Ising planteado
temperatura = 20
npasos = 5000
recorridos = 10

                ### Modelo Ising ejecutado una vez ###
                
#Inicialmente, se reporta si el sistema alcanza el equilibrio y en qué paso lo hace
listagraf, energiaacum, magnetacum, estado, paso_equilibrio = (Modelo_Ising(temperatura, npasos))
print("El estado del sistema es:", estado)
print("El equilibrio se alcanza en el paso:", paso_equilibrio)

#Se crean listas de pasos y energías para cuando el modelo es ejecutado una vez
lista_energías = energiaacum
lista_pasos1 = np.arange(len(lista_energías))
   
#Se grafica la energía del sistema en función de los pasos que se ejecutan para el Modelo de Ising 
#1-D recorrido una vez   
fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(lista_pasos1, lista_energías)
ax.set_title("Gráfico de Energía del sistema vrs Pasos")
ax.set_xlabel("Pasos")
ax.set_ylabel("Energía (E)")
plt.show

#Para reportar el valor de la energía del sistema se promedian las energías a partir del paso en el
#que se alcanza el quilibrio 
E = np.asarray(lista_energías)
valorE = np.mean(E[paso_equilibrio:])
print("El valor de la energía en el sistema a una T=", temperatura, "es:", valorE)

#Se crean listas de pasos y magnetizaciones obtenidas para cuando el modelo es ejecutado una vez
lista_magnetizaciones = magnetacum
lista_pasos2 = np.arange(len(lista_magnetizaciones))

#Se grafica la magnetización del sistema en función de los pasos que se ejecutan en el modelo
fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(lista_pasos2, lista_magnetizaciones)
ax.set_title("Gráfico de Magnetización del sistema vrs Pasos")
ax.set_xlabel("Pasos")
ax.set_ylabel("Magnetización (M)")
plt.show

#Para reportar el valor de la magnetización del sistema se promedian las magnetizaciones a partir del
#paso en el que se alcanza el equilibrio
M = np.asarray(lista_magnetizaciones)
valorM = np.mean(M[paso_equilibrio:])
print("El valor de la magnetización en el sistema a una T=", temperatura, "es:", valorM)


            ### Resultados promedio del Modelo de Ising 1-D ejecutado 10 veces ###

#Se definen las listas que contienen la magnetización y energía promedio para cada recorrido ejecutado            
EProm, MProm = Valores_Promedio(recorridos, temperatura, npasos)
lista_eprom = EProm
lista_mprom = MProm
lista_pasos3 = np.arange(len(lista_eprom))
lista_pasos4 = np.arange(len(lista_mprom))    

#Se grafican las energías y magnetizaciones promedio en función de los pasos ejecutados
fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(lista_pasos3, lista_eprom)
ax.set_title("Gráfico de Energía Promedio vrs Pasos")
ax.set_xlabel("Pasos")
ax.set_ylabel("Energía Promedio (<E>)")
plt.show
    
fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(lista_pasos4, lista_mprom)
ax.set_title("Grafico de Magnetización Promedio vrs Pasos")
ax.set_xlabel("Pasos")
ax.set_ylabel("Magnetización Promedio (<M>)")
plt.show   

#Se reportan el valor de energía y magnetización promedio en 10 recorridas sacando la media de 
#las listas de magnetización y energía promedio por cada recorrido
Energía_Promedio = np.asarray(lista_eprom)
valorEP = np.mean(Energía_Promedio)
print("El valor de la energía promedio para", recorridos, "recorridos a una temperatura T=", temperatura, " es:", valorEP)

Magnetización_Promedio = np.asarray(lista_mprom)
valorMP = np.mean(Magnetización_Promedio)
print("El valor de la magnetización promedio para", recorridos, "recorridos a una temperatura T=", temperatura, "es: ", valorMP)


''' Se reportan los resultados obtenidos de energía interna, magnetización y calor específico 
en función de kb*T del Modelo de Ising 1-D desarrollado'''

#Se definen los parámetros con los que se correrá el Modelo de Ising y demás funciones
npasos1 = 2000
recorridos1 = 10

#Se crea el arreglo de temperaturas a estudiar que va de 0.1 a 5.
temperaturas = np.linspace(0.1, 5, 10)

#Se inicializan las listas de magnetización y energía interna
lista_M = []
lista_U = []

#Se calculan la energía interna y magnetización para cada temperatura a analizar
for i in temperaturas:
    #Se llama al Modelo de Ising que se ejecuta 10 recorridos y se promedian los valores de energía
    #interna y magnetización para cada temperatura y se actualizan las respectivas listas
    energíainter, magnetakb = (Valores_Promedio(recorridos1, i, npasos1))
    energíaU = np.mean(energíainter)
    magnetización_kb = np.mean(magnetakb)
    lista_U.append(energíaU)
    lista_M.append(magnetización_kb)


#Se grafica la energía interna y magnetización en función de la variación de kb*T

fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(temperaturas, lista_U)
ax.set_title("Gráfico de U vrs kbT")
ax.set_xlabel("kb*T")
ax.set_ylabel("Energía Interna (U)")
plt.show

fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(temperaturas, lista_M)
ax.set_title("Gráfico de M vrs kbT")
ax.set_xlabel("kb*T")
ax.set_ylabel("Magnetización Promedio (<M>)")
plt.show


#Se inicializa la lista de energía interna 2
listaU2 = []

#Se calcula la energía interna 2 llamando a la función respectiva y anexando los resultados a la lista
for i in temperaturas:
    energía_inter2 = Energía_Interna_2(i, npasos1)
    listaU2.append(energía_inter2)

#Se inicializa la lista de calores específicos y se definen las variables de U y U2 con sus listas
#correspondientes
listaC = []
energía_interna = lista_U
energía_interna2 = listaU2

#El calor específico es calculado para cada temperatura y se anexan los resultados a la lista C
for i in range(len(temperaturas)):
    C = np.abs(((1/N**2)*(energía_interna2[i]-(energía_interna[i])**2))/(kb*(temperaturas[i])**2))
    listaC.append(C)

#Se grafican los calores específicos obtenidos para cada variación de kb*T
fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(temperaturas, listaC)
ax.set_title("Gráfico de C vs kb*T")
ax.set_xlabel("kb*T")
ax.set_ylabel("Calor Específico (C)")
plt.show   


#Se grafican las variables del calor específico, energía interna y magnetización analíticas en
#función de kb*T
fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(temperaturas, Cálculos_Analíticos(J, kb, N, temperaturas)[0])
ax.set_title("Gráfica de U vs kb*T Analítica")
ax.set_xlabel("kb*T")
ax.set_ylabel("Energía Interna (U)")
plt.show     
    
fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(temperaturas, Cálculos_Analíticos(J, kb, N, temperaturas)[1])
ax.set_title("Gráfica de M vs kb*T Analítica")
ax.set_xlabel("kb*T")
ax.set_ylabel("Magnetización (M)")
plt.show      

fig, ax = plt.subplots(figsize = (5,5), dpi=120)
ax.plot(temperaturas, Cálculos_Analíticos(J, kb, N, temperaturas)[2])
ax.set_title("Gráfico de C vs kb*T Analítica")
ax.set_xlabel("kb*T")
ax.set_ylabel("Calor Específico (C)")
plt.show  