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

#busca de armijo: arumentos: função, ponto (x1, x2), direção (dx1, dx2), taxa de atualização de t (gamma)
# e eta (multiplicador para a direção de descida - ?)
def armijo(func, x1, x2, dx1, dx2, gamma, eta): #retorno: tamanho do passo

	t = 1

	if func == "a":
		while fA(x1 + t * dx1, x2 + t * dx2) > fA(x1, x2) + eta * t * (grad_fA(x1, x2)[0] * dx1 + grad_fA(x1, x2)[1] * dx2):
			t = gamma * t

	if func == "b":
		while fB(x1 + t * dx1, x2 + t * dx2) > fB(x1, x2) + eta * t * (grad_fB(x1, x2)[0] * dx1 + grad_fB(x1, x2)[1] * dx2):
			t = gamma * t

	return t

#Argumentos: funcao, minimizador, ponto a ser verificado, direcao de descida, funcao gradiente
#def gradient (func, minimazer, x, d, grad):
#	t = 1
#	while func(x + t*d) > ( func(x) + NI*t*grad(x)*d ):
#		t = LAMBA*t
#
#O que é NI e LAMBDA???




x1 = 2
x2 = 5

fB(x1, x2)
