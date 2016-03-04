# coding=utf-8


#Outra coisa a se notar: aquelas divisoes de vetor voltaram... vamos ter que consertar isso no quase-newton.
#Fiz umas funcoes de transposicao que podem facilitar

import sys
import math
import numpy

#--------------Constantes-----------

_GAMMA_ = 0.8
_ETA_ = 0.25
_ITER_ = 1000
_TOLERANCIA_ = 1/10**8

#_TIME_ = 6000 #min

#--------------Funcoes--------------

def fA (x1, x2):
	return - math.exp( (-1)*pow(x1, 2) - pow(pow(x1, 2) - x2, 2) )

def fB (x1, x2):
	return  math.sqrt( pow(x1, 2) + pow(  (pow(x1, 2) - x2), 2) )

def grad_fA(x1, x2):
	deriv_x1 = (-1)*( 2 * x1 + 4 * pow(x1, 3) - 4 * x1 * x2 ) * fA(x1, x2)
	deriv_x2 = 2 * ( pow(x1, 2) - x2 ) * fA(x1, x2)
	return [deriv_x1, deriv_x2]

def grad_fB(x1, x2):
	result = fB(x1, x2)
	if(result == 0.0):
		return [0, 0]
	else:
		deriv_x1 = (x1 + 2 * pow(x1, 3) - 2 * x1 * x2) / result
		deriv_x2 = ( x2 - pow(x1, 2) ) / result
		return [deriv_x1, deriv_x2]

def hess_fA(x1, x2):
	deriv_x1_x1 = (-1)*(2 + 12 * pow(x1, 2) - 4 * x2) * fA(x1, x2) - (2 * x1 + 4 * pow(x1, 3) - 4 * x1 * x2) * grad_fA(x1, x2)[0]
	deriv_x1_x2 = 4 * x1 * fA(x1, x2) - (2 * x1 + 4 * pow(x1, 3) - 4 * x1 * x2) * grad_fA(x1, x2)[1]

	deriv_x2_x1 = 4 * x1 * fA(x1, x2) + 2 * (pow (x1, 2) - x2) * grad_fA(x1, x2)[0]
	deriv_x2_x2 = (-1)*2 * fA(x1, x2) + 2 * (pow(x1, 2) - x2) * grad_fA(x1, x2)[1]

	return [[deriv_x1_x1, deriv_x1_x2], [deriv_x2_x1, deriv_x2_x2]]

def hess_fB(x1, x2):
	deriv_x1_x1 = ( (1 + 6 * pow (x1, 2) - 2 * x2) * fB(x1, x2) - (x1 + 2 * pow(x1, 3) - 2 * x1 * x2) * grad_fB(x1, x2)[0] ) / pow (fB(x1, x2), 2)
	deriv_x1_x2 = (-2 * x1 * fB(x1, x2) - (x1 + 2 * pow(x1, 3) - 2 * x1 * x2) * grad_fB(x1, x2)[1]) / pow(fB(x1, x2) , 2)

	deriv_x2_x1 = (-2 * x1 * fB(x1, x2) - (x2 - pow(x1, 2)) * grad_fB(x1, x2)[0]) / pow(fB(x1, x2), 2)
	deriv_x2_x2 = (fB(x1, x2) - (x2 - pow(x1, 2)) * grad_fB(x1, x2)[1]) / pow(fB(x1, x2), 2)

	return [[deriv_x1_x1, deriv_x1_x2], [deriv_x2_x1, deriv_x2_x2]]


#busca de armijo: argumentos: funcao, gradiente da funcao, ponto x[x1, x2], direcao d[dx1, dx2], taxa de atualizacao de t (gamma)
def armijo(func, grad_func, x1, x2, dx1, dx2, gamma, eta): #retorno: tamanho do passo
	t = 1
	contador = 0
	while func(x1 + t*dx1, x2 + t*dx2) > (func(x1, x2) + eta*t*(grad_func(x1, x2)[0]*dx1 + grad_func(x1, x2)[1]*dx2)):
		t = gamma*t
		contador+=1
	return [t, contador]

#Argumentos: funcao, gradiente da funcao, ponto x[x1, x2]
def gradiente(func, grad_func, x1, x2):
	global _GAMMA_
	x = [x1, x2]
	x_ant = x
	x_0 = [x1, x2]
	contador = 0

	while ( norma(grad_func(x[0], x[1])) > _TOLERANCIA_ and contador < _ITER_ ):
		d = [(-1) * value for value in grad_func(x[0], x[1])]
		a = armijo(func, grad_func, x[0], x[1], d[0], d[1], _GAMMA_, _ETA_)
		t = a[0]
		call_armijo = a[1]

		x_prox = [x[0] + t*d[0], x[1] + t*d[1]]
		x_ant = x
		x = x_prox
		contador += 1
		print(x_prox, func(x_prox[0], x_prox[1]))
		if (x_ant == x):
			break

	if(contador==_ITER_):
		print("Numero maximo de iteracoes")
	return [x_0, contador, call_armijo, x, func(x[0], x[1])]


def newton(func, grad_func, hess_func, x1, x2):
	
	contador = 0
	x0 = [x1, x2]
	x = [x1, x2]
	x_ant = x

	while ( ( norma(grad_func(x[0], x[1]) ) > _TOLERANCIA_ and contador < _ITER_) ):
		aux1 = invMatrix2x2( hess_func(x[0], x[1]) )
		if (aux1 is None): #se a matriz for singular, não ha inversa
			x_prox = x
			break;

		aux2 = [ [(-1) * value for value in aux1[0]], [(-1) * value for value in aux1[1]] ]
		aux3 = grad_func(x[0], x[1])
		d = m_MV(aux2, aux3)
		a = armijo(func, grad_func, x[0], x[1], d[0], d[1], _GAMMA_, _ETA_)
		print(a)
		t = a[0]
		call_armijo = a[1]

		x_prox = [x[0] + t*d[0], x[1] + t*d[1]]
		x_ant = x
		x = x_prox

		contador+=1

		if (x_ant == x):
			break;

	if (contador == _ITER_):
		print ("Numero maximo de iteracoes")

	return [x0, contador, call_armijo, x_prox, func(x_prox[0], x_prox[1])]


def quaseNewton(func, grad_func, hessiana, x1, x2):
	k = 0
	x1_ant = x1
	x2_ant = x2
	x0 = [x1, x2]
	contador = 0

	#Para a primeira iteracao definimos Hk = Identidade
	Hk = [[1, 0], [0, 1]]
	while (((grad_func(x1, x2)[0]!=0) and (grad_func(x1, x2)[1]!=0)) or ((x1 == x1_ant) and (x2 == x2_ant))):
		aux0 = m_MV( Hk, grad_func(x1, x2) )
		print("Hk= ", Hk)
		print("grad= ", grad_func(x1, x2))
		d = [(-1) * value for value in aux0]
		print ("d= ", d)
		a = armijo(func, grad_func, x1, x2, d[0], d[1], _GAMMA_, _ETA_)
		t = a[0]
		call_armijo = a[1]
		print(a)
		x1_prox = x1 + t*d[0]
		x2_prox = x2 + t*d[1]

		p = [x1_prox - x1, x2_prox - x2]
		q = [(grad_func(x1_prox, x2_prox)[0] - grad_func(x1, x2)[0]), (grad_func(x1_prox, x2_prox)[1] - grad_func(x1, x2)[1])]
		print(t, d)
		aux1 = 1 + int( (m_VV((m_VM(q, Hk)), q) / m_VV(p, q)) )
		aux2 = aux1 / m_VV(p,q)
		aux3 = [ [(aux2) * value for value in m_VVT(p,p)[0]], [aux2 * value for value in m_VVT(p,p)[1]] ]  
		aux4 = s_MM(m_MM( m_VVT(p, q), Hk), m_VVT(m_MV(Hk, q), p))
		aux5 = m_VV(p,q)
		aux6 = [ [((-1)*value/aux5) for value in aux4[0]], [ ((-1)*value/aux5) for value in aux4[1] ] ]  

		hess_Est = s_MM(s_MM(Hk, aux3), aux6)

		Hk = hess_Est
		k = k + 1

		x1_ant = x1
		x2_ant = x2
		x1 = x1_prox
		x2 = x2_prox
		contador+=1

	return [x0, contador, call_armijo, [x1, x2], func(x1, x2)]

#Multiplicacao Matriz2x2 x Vetor
def m_MV(m, v):
	a = m[0][0]*v[0] + m[0][1]*v[1]
	b = m[1][0]*v[0] + m[1][1]*v[1]
	return [a, b]

#Multiplicacao Vetor x Matriz2x2
def m_VM(v, m):
	a = v[0]*m[0][0] + v[1]*m[1][0]
	b = v[0]*m[0][1] + v[1]*m[1][1]
	return [a, b]

#Multiplicacao Vetor x Vetor
def m_VV (v1, v2):
	return v1[0]*v2[0] + v1[1]*v2[1]

def m_VVT (v1, v2):
	return [[v1[0]*v2[0], v1[0]*v2[1]], [v1[1]*v2[0], v1[1]*v2[1]]]

def m_MM(m1, m2):
	return [ [ m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0], m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] ], [ m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0], m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] ] ]

def s_MM(m1, m2):
	return [ [ m1[0][0]+m2[0][0], m1[0][1]+m2[0][1] ], [ m1[1][0]+m2[1][0], m1[1][1]+m2[1][1] ] ]

#Somando Vetor + Vetor
def s_VV(v1, v2):
	return [v1[0]+v2[0], v1[1]+v2[1]]

def determinante(M):
	return M[0][0] * M[1][1] - M[0][1] * M[1][0]

#inversa
def invMatrix2x2(M):

	det = determinante(M)

	if det == 0:
		return None

	return [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]

def norma(v):
	return math.sqrt( v[0] ** 2 + v[1] ** 2 )

def vetorTransposto(v):
	return [[v[0]], [v[1]]]

def matrizTransposta(M):
	return [ [M[0][0], M[1][0]], [M[0][1], M[1][1] ] ]

def autovalores(M):
	av1 = ( (M[0][0] + M[1][1]) + math.sqrt( (M[0][0] + M[1][1]) ** 2 - 4 * determinante(M) ) ) / 2
	av2 = ( (M[0][0] + M[1][1]) - math.sqrt( (M[0][0] + M[1][1]) ** 2 - 4 * determinante(M) ) ) / 2

	return [av1, av2]

def ePositivaDefinida(M):
	return (autovalores(M)[0] > 0) and (autovalores(M)[1] > 0)

#-----------Codigo principal------------

# Teste para a funcao de armijo, tudo ok
#def teste(x1, x2):
#	return (0.5)*pow( (x1 - 2) , 2) + pow( (x2 - 1) , 2)
#def grad_teste(x1, x2):
#	return [x1 - 2, 2*x2 - 2]
#print(armijo(teste, grad_teste, 1, 0, 3, 1, 0.8, (0.25)))


#Abaixo, teste para o metodo de Newton: FUNCIONA!

#x aqui é um vetor com x1 e x2
# o otimo deve dar 0 no ponto (0, 0)

def testF(x1, x2):
	return pow(x1, 2) + 2 * pow(x2, 2)

def grad_testF(x1, x2):
	return[2 * x1, 4 * x2]

def hess_testF(x1, x2):
	return [[2, 0], [0, 4]]


#print( newton(testF, grad_testF, hess_testF, 7.5, 9 ) )
#print( gradiente(testF, grad_testF, 7.5, 9 ) )

#print( quaseNewton(testF, grad_testF, hess_testF, 1, 1 ) )
#print( newton(fB, grad_fB, hess_fB, 7.5, 9 ) )

#FIM do teste para Newton

#print( gradiente(fB, grad_fB, 0.1, 1 ) )

print(quaseNewton(fA, grad_fA, hess_fA, 5, 50))


