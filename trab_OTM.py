import sys
import math
import numpy

def fA (x1, x2):
	return - math.exp( -pow(x1, 2) - pow(pow(x1, 2)- x2, 2) )

def fB (x1, x2):
	return  math.sqrt( pow(x1,2) + pow(  (pow(x1,2) - x2), 2) )

def grad_fA(x1, x2):
	deriv_x1 = - ( 2 * x1 + 4 * pow(x1, 3) - 4 * x1 * x2 ) * fA(x1, x2)
	deriv_x2 = 2 * ( pow(x1, 2) - x2 ) * fA(x1, x2)
	return [deriv_x1, deriv_x2]

def grad_fB(x1, x2):
	deriv_x1 = (x1 + 2 * pow(x1, 3) - 2 * x1 * x2) / fB(x1, x2)
	deriv_x2 = ( x2 - pow(x1, 2) ) / fB(x1, x2)
	return [deriv_x1, deriv_x2]

#busca de armijo: arumentos: funcao, gradiente da funcao, ponto x[x1, x2], direcao d[dx1, dx2], taxa de atualizacao de t (gamma)
def armijo(func, grad_func, x1, x2, dx1, dx2, gamma, eta): #retorno: tamanho do passo
	t = 1

	while func(x1 + t*dx1, x2 + t*dx2) > func(x1, x2) + eta*t*(grad_func(x1, x2)[0]*dx1 + grad_fA(x1, x2)[1]*dx2):
		t = gamma*t
	print (x1 + t*dx1, x2 + t*dx2)
	return t


#Argumentos: funcao, gradiente da funcao, ponto x[x1, x2], direcao d[dx1, dx2]
#def gradiente(func, grad_func, x1, x2, dx1, dx2):
#	k = 0
#	x1_ant = x1
#	x2_ant = x2

#	while (((grad_func(x1, x2)[0]!=0) && (grad_func(x1, x2)[1]!=0)) || ((x1 == x1_ant) && (x2 == x2_ant))):
#		d = -(grad_func(x1, x2))
		
#		t = armijo(func, grad_func, x1, x2, d[0], d[1], gamma, eta)
		
#		x1_prox = x1 + t*d[0]
#		x2_prox = x2 + t*d[1]

#		k = k+1

#		x1_ant = x1
#		x2_ant = x2
#		x1 = x1_prox
#		x2 = x2_prox

#	return [x1, x2]


#def newton(func, grad_func, x):
#	k = 0
#	x_ant = x
#	while ( (grad_func != 0) || (x_ant == x) ):
#		d = -(inv_hess(x)) * grad_func(x)
#		t = armijo(func, x, d, grad_func)
#		x_prox = x + t*d
#		k = k + 1
#		x_ant = x
#		x = x_prox
#	return x


def quaseNewton(func, grad_func, x1, x2):
	k = 0
	x1_ant = x1
	x2_ant = X2
	#Para a primeira iteracao definimos Hk = Identidade
	Hk = [[1, 0], [0, 1]]
	while (((grad_func(x1, x2)[0]!=0) and (grad_func(x1, x2)[1]!=0)) or ((x1 == x1_ant) and (x2 == x2_ant))):
		d = - (m_MatrizVetor( (hessiana(func, grad_func, x1, x2)), grad_func(x1, x2) ) )
		t = armijo(func, grad_func, x1, x2, d[0], d[1], gamma, eta)
		x1_prox = x1 + t*d[0]
		x2_prox = x2 + t*d[1]


		p = [x1_prox - x1, x2_prox - x2]
		q = grad_func(x1_prox, x2_prox) - grad_func(x1, x2)
		
#		np.dot(v1, v2) -> multiplicacao de vetores
#		hess_Est = Hk + { 1 + [(m_VetorMatriz(q, Hk))*q ]    }
		k = k + 1
		x_ant = x
		x = x_prox
	return x

#Argumentos: Matrix hessiana, vetor gradiente
def m_MatrizVetor(m, v):
	a = m[0][0]*v[0] + m[0][1]*v[1]
	b = m[1][0]*v[0] + m[1][0]*v[1]
	return [a, b]
	
#Argumentos: Vetor gradiente, matrix hessiana
def m_VetorMatriz(v, m):
	a = v[0]*m[0][0] + v[1]*m[1][0]
	b = v[0]*m[0][1] + v[1]*m[1][1]
	return [a, b]



# Teste para a funcao de armijo, tudo ok
def teste(x1, x2):
	return (0.5)*pow( (x1 - 2) , 2) + pow( (x2 - 1) , 2)

def grad_teste(x1, x2):
	return [x1 - 2, 2*x2 - 2]

armijo(teste, grad_teste, 1, 0, 3, 1, 0.8, (0.25))
