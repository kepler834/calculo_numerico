import matplotlib.pyplot as plt
import numpy as np

#recebendo funcao entrada
funcao = input("informe a funcao:")
derivada = input("informe a derivada da funcao:")
iteracao = int(input("informe o numero iteracoes:"))
intervalo_sup = float(input("informe o intervalo superior:"))
intervalo_infe = float(input("informe o intervalo inferior:"))
chut_init = intervalo_sup

def f(x):
    y = funcao
    return float(eval(y))

def d(x):
    y = derivada
    return float(eval(y))

def bisseccao(a,b):
    if f(a)*f(b)<0:
        ciclos = 0
        while ciclos < iteracao:
            ciclos += 1
            x = (a+b)/2
            if f(a)*f(x)<0:
                b = x
            elif f(b)*f(x)<0:
                a = x
        print("\nciclos : ",ciclos)
        return (a+b)/2
    else:
        print("INTERVALO INVALIDO")


def posicaoFalsa(a,b):
    if f(a)*f(b)<0:
        ciclos = 0
        while ciclos < iteracao:
            ciclos += 1
            x = ((a*f(b))-(b*f(a)))/(f(b)-f(a))
            if f(a)*f(x)<0:
                b = x
            elif f(b)*f(x)<0:
                a = x
        print("\nciclos : ",ciclos)
        return (a+b)/2
    else:
        print("INTERVALO INVALIDO")


def newton(a):
    ciclos = 0
    while ciclos < iteracao:
        ciclos += 1
        a = a - (f(a)/d(a))
    print("\nciclos : ",ciclos)
    return a


def secante(a,b):
    ciclos = 0
    while ciclos < iteracao:
        ciclos += 1
        x = a - (f(a)/((f(a)-f(b))/(a-b)))
        b = a
        a = x
    print("\nciclos : ",ciclos)
    return a



#utilizando metodo da bisseccao
print("\n-------------bisseccao-------------")
biss = bisseccao(intervalo_infe, intervalo_sup)
print("x =",biss,", f(x) =",f(biss))

#utilizando metodo da posicao falsa
print("\n-----------posicao falsa-----------")
posf = posicaoFalsa(intervalo_infe, intervalo_sup)
print("x =",posf,", f(x) =",f(posf))

#utilizando metodo de newton
print("\n--------------newton---------------")
newt = newton(intervalo_sup)
print("x =",newt,", f(x) =",f(newt))

#utilizando metodo das secantes
print("\n-------------secante---------------")
sec = secante(intervalo_infe, intervalo_sup)
print("x =",sec,", f(x) =",f(sec),"\n\n")

array = [biss, posf, newt, sec]

#exibindo a funcao graficamente
eixo_x = np.arange(intervalo_infe, intervalo_sup, 0.000001)
fv = np.vectorize(f)
plt.plot(eixo_x, fv(eixo_x), label = "funcao")

#exibindo os metodos numericos graficamente
plt.plot(biss, f(biss), 'ko', label = "bisseccao")
plt.plot(posf, f(posf), 'go', label = "posicaoFalsa")
plt.plot(newt, f(newt), 'ro', label = "newton")
plt.plot(sec, f(sec), 'yo', label = "secante")

#exibindo legendas
plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.legend()
plt.grid(True)

plt.show()