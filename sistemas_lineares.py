import numpy as np
# sistema linear
# 1A + 1B + 2C = 9
# 2A + 4B - 3C = 1
# 3A + 6B - 5C = 0

coeficientes = np.array([[1,2,1],[2,3,1],[1,1,2]])
var_independentes = np.array([[10],[15],[9]])
matriz_aumentada = np.concatenate((coeficientes,var_independentes),axis=1)
chute_inicial = np.array([1,1,1])
ciclos = 3

def sistemaValido():
    determinante = np.linalg.det(coeficientes)
    if determinante == 0:
        return False
    else:
        return True

def gaussJordam():
    matriz_aux = matriz_aumentada.astype(float)
    for idxA, valA in enumerate(matriz_aux, 0):
        #tornando pivor igual a 1
        matriz_aux[idxA] = np.multiply(valA, 1/valA[idxA])
        for idxB, valB in enumerate(matriz_aux, 0):
            if idxA != idxB:
                #zerando os fatores abaixo do pivor
                matriz_aux[idxB] = valB - np.multiply(valB[idxA], valA)
    #criando lista com o resultando
    vet_resul = [0]*len(matriz_aux)
    for idx, val in enumerate(matriz_aux, 0):
        vet_resul[idx] = val[len(val)-1]
    return vet_resul

def fatoracaoLU():
    upper = coeficientes.astype(float)
    lower = np.identity(len(coeficientes), dtype = float)
    #determinando a matriz triangular superior e inferior
    for idxA, valA in enumerate(upper, 0):
        for idxB, valB in enumerate(upper, 0):
            if idxA < idxB:
                razao = valB[idxA]/valA[idxA]
                upper[idxB] = valB - np.multiply(razao, valA)
                lower[idxB][idxA] = razao
    sistema_l = var_independentes.astype(float)
    #resolvendo sistema ly = b
    for idxA in range(len(lower)):
        for idxB in range(idxA):
            sistema_l[idxA] = sistema_l[idxA] - sistema_l[idxB]*lower[idxA][idxB]
        sistema_l[idxA] = sistema_l[idxA]/lower[idxA][idxA]
    sistema_u = sistema_l
    #resolvendo sistema ux = y
    for idxA in range(len(upper)-1, -1, -1):
        for idxB in range(idxA+1, len(upper)):
            sistema_u[idxA] = sistema_u[idxA] - sistema_u[idxB]*upper[idxA][idxB]
        sistema_u[idxA] = sistema_u[idxA]/upper[idxA][idxA]
    return sistema_u

def inversa():
    matriz_expandida = np.concatenate((coeficientes,np.identity(len(coeficientes), dtype = float)),axis=1)
    #pivoteamento
    for idxA, valA in enumerate(matriz_expandida, 0):
        matriz_expandida[idxA] = np.multiply(valA, 1/valA[idxA])
        for idxB, valB in enumerate(matriz_expandida, 0):
            if idxA != idxB:
                matriz_expandida[idxB] = valB - np.multiply(valB[idxA], valA)
    inversa = np.split(matriz_expandida, 2, axis = 1)
    return np.matmul(inversa[1], var_independentes)

def cramer():
    resul = [0]*len(var_independentes)
    determinante = np.linalg.det(coeficientes)
    matrix = coeficientes.astype(float)
    for i in range(0, len(coeficientes)):
        matrix = coeficientes.astype(float)
        matrix[:len(coeficientes),i:i+1] = var_independentes
        resul[i] = np.linalg.det(matrix)/determinante
    return resul

def gaussJacobi():
    matriz_anterior = chute_inicial
    matriz_posterior = np.zeros(len(matriz_anterior))
    for i in range(0, ciclos):
        for idx, val in enumerate(coeficientes,0):
            somatorio = np.dot(matriz_anterior,val) - (matriz_anterior[idx]*val[idx])
            matriz_posterior[idx] = (var_independentes[idx] - somatorio)/ val[idx]
        matriz_anterior = np.copy(matriz_posterior)
    return matriz_posterior

def gaussSeidel():
    matriz_anterior = chute_inicial
    matriz_posterior = np.zeros(len(matriz_anterior))
    for i in range(0, ciclos):
        for idx, val in enumerate(coeficientes,0):
            somatorio = np.dot(matriz_anterior,val) - (matriz_anterior[idx]*val[idx])
            matriz_posterior[idx] = float(var_independentes[idx] - somatorio)/ float(val[idx])
            matriz_anterior[idx] = matriz_posterior[idx]
    return matriz_posterior

#verificando se o sistema definido é valido
if sistemaValido() == False:
    quit()

#encontrando a solucao do sistema pelo metodo de Gauss-Jordam
print("\n-----------gaussJordam------------")
gj = gaussJordam()
print(gj)

#encontrando a solucao do sistema pelo metodo de fatoracao LU
print("\n-----------fatoracaoLU------------")
fat = fatoracaoLU()
print(fat.transpose())

#encontrando a solucao do sistema pelo metodo da matriz inversa
print("\n-------------inversa--------------")
inv = inversa()
print(inv.transpose())

#encontrando a solucao do sistema pelo metodo de cramer
print("\n-------------cramer---------------")
cra = cramer()
print(cra)

#encontrando a solucao do sistema pelo metodo de Gauss-Jacobi
print("\n-----------gaussJacobi------------")
gja = gaussJacobi()
print(gja)

#encontrando a solucao do sistema pelo metodo de Gauss-seidel
print("\n-----------gaussSeidel------------")
gs = gaussSeidel()
print(gs,'\n')