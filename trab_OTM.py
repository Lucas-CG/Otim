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


def armijo():
#Argumentos: funcao, ponto a ser verificado, direcao de descida, funcao gradiente
#def gradient (func, x, d, grad):
#	t = 1
#	while func(x + t*d) > ( func(x) + NI*t*grad(x)*d ):
#		t = LAMBA*t
#	return t
#O que é NI e LAMBDA???

def gradiente(func, grad_func, x):
#	k = 0
#	x_ant = x
#	while ((grad_func(x)!=0) || (x == x_ant)):
#		d = -grad_func(x)
#		t = armijo(func, x, d, grad_func)
#		x_prox = x + t*d
#		k = k+1
#		x_ant = x
#		x = x_prox
#	return x


def newton(func, grad_func, x):
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


def quaseNewton(func, grad_func, x):
#	k = 0
#	x_ant = x
#	Hk = DEFINIR
#	while ( (grad_func(x)!=0) || (x_ant==x) ):
#		d = -(Hk)*grad_func(x)
#		t = armijo(func, x, d, grad_func)
#		x_prox = x + t*d
#		Hk_prox = H_k(x_prox, Hk)
#		k = k + 1
#		x_ant = x
#		x = x_prox

def Hk_prox(x, Hk):
	


x1 = 2
x2 = 5

fB(x1, x2)
