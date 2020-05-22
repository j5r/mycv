# para utilizar, abra um interpretador python como o IDLE
# e importe este arquivo conforme
# >> from exercicio import *
# e então pode usá-lo conforme a descrição abaixo
# pode ser necessário instalar numpy e matplotlib
# para tanto, entre num terminal CMD com
# >> pip install numpy matplotlib

"""
> help()
    abre uma janela de ajuda, que se fecha com a letra Q

> raio(A)
    < retorna o raio espectral de uma matriz quadrada [A]
        (o maior autovalor em módulo)

> run(coeff,iteracoes=1000)
    executa o cálculo da Lyapunov com os parâmetros em memória
    lyapunov(A,B,Q,Sigma0,iteracoes), e plota os gráficos de E||x||²
    > coeff
        coeficiente para multiplicar a matriz [A] antes de executar
        a lyapunov. Inicialmente é calculado como sendo
        coeff = 1/raio(A) - 0.01

> lyapunov(A,B,Q,Sigma0,iteracoes=1000)
    calcula a matriz Sigma_inf, os traços de cada Sigma_k e os erros
    tr(Sigma_k - Sigma_{k-1})
    < retorna  Sigma_inf, traços, erros

> gerarMatriz(m,n)
    < retorna uma matriz m x n de números conforme o modelo normal(0,1)

> gerarDados(n,m)
    gera matrizes com a rotina gerarMatriz(m,n) com as dimensões
    coerentes e as simetrias necessárias
    > n
        é a dimensão de x,
    > m
        é a dimensão do ruído w
    < retorna A,B,Q,Sigma0

> randomizar(n,m)
    gera novos dados em memória, utilizando a rotina gerarDados(n,m)
    > n
        é a dimensão de x,
    > m
        é a dimensão do ruído w
"""

from numpy import matrix, trace
from numpy.linalg import eig
import matplotlib.pyplot as plt
import random



########################################################################
__help = help
def help():
    __help(__name__)
def raio(A):
    return max(abs(eig(A)[0]))

def gerarMatriz(m,n):
    rows = []
    for i in range(m):
        rows.append([random.gauss(0,1) for j in range(n)])
    return matrix(rows)

def gerarDados(n,m):
    A = gerarMatriz(n,n)
    A = A/raio(A)
    B = gerarMatriz(n,m)
    Q = gerarMatriz(m,m)
    Q = Q @ Q.T #simetrizando
    S0 = gerarMatriz(n,n)
    S0 = S0 @ S0.T #simetrizando
    print()
    print("\tRaio Espectral rho(A) = ", raio(A))
    print("\tRaio Espectral rho(BB') = ", raio(B @ B.T))
    print("\tRaio Espectral rho(Q) = ", raio(Q))
    print("\tRaio Espectral rho(Sigma0) = ", raio(S0))
    print()
    coeff = 1/raio(A) - 0.01
    return A,B,Q,S0


def lyapunov(A,B,Q,S0,iteracoes=1000):
    erros = []
    tracos = []
    S = S0.copy()
    S_anterior = S.copy()
    Add = B @ Q @ B.T
    for i in range(iteracoes):
        tracos.append(trace(S))
        S = A @ S @ A.T + Add
        erros.append(trace(S - S_anterior))
        S_anterior = S.copy()

    return S,tracos,erros

# [coeff] esse coeficiente é um fator de escala para a matriz
# [iteracoes] é o número de iterações da Lyapunov
def run(coeff,iteracoes=1000):
    Sigma,tracos,erros = lyapunov(coeff*A,B,Q,Sigma0,iteracoes)
    print()
    print("Raio de (coeff*A) =",raio(coeff*A))
    print()
    print("Sigma_inf = ")
    print(Sigma)

    plt.subplot(2,1,1)
    plt.plot(tracos,"b.-")
    plt.title("Plot de E||x_k||²")

    plt.subplot(2,1,2)
    plt.plot(erros,"b.-")
    plt.title("Plot do erro de cálculo de E||x||²")
    plt.show()



# Gerar dados aleatoriamente
def randomizar(n,m):
    A,B,Q,Sigma0 = gerarDados(n,m)
    raioinversoA = 1/raio(A) - 0.01





########################################################################
A = matrix([# matriz n x n (quadrada)
        [1,2,3],
        [-2,4,1],
        [7,-2,-3]
    ])

raioinversoA = 1/raio(A) - 0.01

B = 100*matrix([#matriz n x m (retangular)
        [1,2],
        [-2,1],
        [3,-2]
    ])
Q = 100*matrix([#matriz m x m (simétrica positiva semidefinida)
        [1,2],
        [2,12]
    ])
Sigma0 = 100*matrix([#matriz n x n (simétrica positiva semidefinida)
            [0.1076667,   0.0068587,   0.0712888],
            [0.0068587,   0.2424449,   0.0498955],
            [0.0712888,   0.0498955,   0.5460524]
        ])


########################## RUN
iteracoes = 100
run(0.1*raioinversoA,iteracoes)




