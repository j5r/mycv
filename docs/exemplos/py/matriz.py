oprint = print;
def print(*args):
    oprint(*args,end="")

######################## INICIO
matriz = [
    [1,2,3,4,5],
    [6,7,8,9,10],
    [11,12,13,14,15]
]

print("Vamos imprimir uma matriz","\n")

for linha in range(len(matriz)):
    for coluna in range(len(matriz[linha])):
        print(matriz[linha][coluna], "\t")
    print("\n")

print("Fim do programa","\n")

######################## FIM
