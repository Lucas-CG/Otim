import sys
import math
import numpy  

def fA (x1, x3):
	return - math.exp( -pow(x1, 2) - pow(pow(x1, 2)- x2, 2) )

def fB (x1, x2):
	return  math.sqrt( pow(x1,2) + pow(  (pow(x1,2) - x2), 2) )

def grad_fA(x1, x2):
	deriv_x1 =
	deriv_x2 =
	return [deriv_x1, deriv_x2]

def grad_fB(x1, x2):
	deriv_x1 =
	deriv_x2 =
	return [deriv_x1, deriv_x2]

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


