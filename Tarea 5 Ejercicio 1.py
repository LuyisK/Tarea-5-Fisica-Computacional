# -*- coding: utf-8 -*-
"""
Física Computacional - TS13

Valentina Campos Aguilar
Luis Alfredo Guerrero Camacho 

MÉTODO MONTE CARLO: Difusión de fotones en el sol
"""

import numpy as np
import matplotlib.pyplot as plt

''' Se crea la función que permite generar el camino aleatorio del fotón con parámetros:
    - el número de iteraciones
    - la lista de posiciones inicial 
    - el paso medio en cada iteracion
La misma retorna la posición del fotón y la lista actualizada de posiciones una vez ejecutado el 
camino aleatorio. '''

def CaminoAleatorio(lista_pos, nIter,pasomedio):
    
    #Se establece el tamaño del paso aleatorias x, y y z
    xPos = (np.random.random()-0.5)*2
    yPos = (np.random.random()-0.5)*2
    zPos = (np.random.random()-0.5)*2
    
    #Se calcula y define L como la distancia que recorre el fotón
    L=np.sqrt((xPos)**2+(yPos)**2+(zPos)**2) 
    
    #Se establecen los cambios de posición en x, y y z que sufrirá el fotón según el camino medio
    #que este recorre
    deltaX= (1/L)*xPos*pasomedio
    deltaY= (1/L)*yPos*pasomedio
    deltaZ= (1/L)*zPos*pasomedio
    
    #Se definen las nuevas posiciones en x, y y z sumando los cambios de posición a las posiciones
    #anteriores
    nuevaxPos = lista_pos[0][nIter-1] + deltaX
    nuevayPos = lista_pos[1][nIter-1] + deltaY
    nuevazPos = lista_pos[2][nIter-1] + deltaZ
    
    #Se actualiza la lista de posiciones con las nuevas posiciones x, y y z
    lista_pos[0].append(nuevaxPos)
    lista_pos[1].append(nuevayPos)
    lista_pos[2].append(nuevazPos)
    

    return lista_pos,L

#Se disminuye el valor del verdadero radio de la región solar
R=0.01 
#Se define el camino medio que realiza el fotón en cada paso
l=5e-5
#Se calcula la cantidad n de pasos promedio que debería dar el fotón para salir de la región solar
Nteórico= (R/l)**2

#Se define la velocidad de la luz    
c=299792458

#Se define el valor arbitrario de veces que se ejecutará el camino aleatorio
ValorArbitrario= 10

#Se inicializan las variables de distancia total recorrida por el fotón (Ltotal), cantidad total 
#de pasos (nTotal) y tiempo total que recorre el fotón (t)
Ltotal=0
nTotal=0
t=0

#Se crea un ciclo que ejecute el algoritmo el número de veces definido
for i in range(0,ValorArbitrario+1):
    #Se define la lista de posiciones inicial, donde el fotón inicia desde el origen
    lista_posiciones= [[0.],[0.],[0.]]
    #Se inicializan las listas de radios y de posiciones
    lista_radios=[]
    lista_pasos=[]
    #Se inicializa el número de pasos en 1
    nPasos=1
    continuar=True
    #Se crea un ciclo que ejecuta el camino aleatorio hasta que la posición del fotón sea mayor 
    #al radio de la región (R)
    while continuar==True:
        #Se ejecuta el camino aleatorio con la lista de posiciones inicial
        lista_posiciones,L= CaminoAleatorio(lista_posiciones, nPasos,l)
        #Se actualiza la distancia total recorrida por el fotón sumando la distancia que este recorre
        #en cada paso
        Ltotal+=L
        #Se calcula el radio de la posición del fotón en cada paso
        Rn = ((np.abs(lista_posiciones[0][nPasos-1]))**2 + (np.abs(lista_posiciones[1][nPasos-1]))**2 + (np.abs(lista_posiciones[2][nPasos-1]))**2)**(1/2)
        #Se actualiza la lista de radios que toma el fotón en cada paso
        lista_radios.append(Rn)
        #Se genera un condicional para determinar si el radio desde el origen hasta  
        #la posicion actual del fotón es mayor que el radio de la región
        if Rn < R:
            continuar=True
            #Se va actualizando la cantidad de pasos que va realizando el fotón
            nPasos+=1
        else:
            continuar=False
        #Se actualiza la lista de pasos ejecutados
        lista_pasos.append(nPasos)
    #Se suman la cantidad de pasos ejecutados por el fotón en todas las iteraciones
    nTotal+=nPasos
    #Se calcula el tiempo total que tarda el fotón en salir de la región 
    t+=Ltotal/c
    #Se calcula el número de pasos promedio que debe ejecutar el fotón para salir de la región
    nPromedio=nTotal/ValorArbitrario
    
#Se calcula el tiempo promedio que tarda el fotón en salir de la región en años
tpromedio=(t/ValorArbitrario)*(1/60)*(1/60)*(1/24)*(1/365.25)

#Lpromedio=Ltotal/ValorArbitrario

print("El N teórico es: ", Nteórico)    
print("El promedio de pasos que deberá ejecutar el fotón será", nPromedio, "para",ValorArbitrario, "iteraciones")
print("El fotón tardará en salir de la región: ",tpromedio, 'años')
#print(lista_posiciones)
#print(lista_radios)


#Se grafica el camino aleatorio del fotón en 3 dimensiones
ax = plt.axes(projection="3d")
plt.figure(figsize=(10,6)) 
ax.plot3D(lista_posiciones[0], lista_posiciones[1],lista_posiciones[2], c='m', lw=0.5)
ax.set_title('Camino aleatorio de un fotón')
plt.show()

#Se hace un gráfico de posición en función de los pasos que ejctuta el fotón para identificar en qué
#rango de pasos este sale de la región solar
plt.figure(figsize=(10,6))
plt.scatter(lista_pasos, lista_radios)
plt.xlabel('pasos')
plt.ylabel('posición (R)')
plt.show()

