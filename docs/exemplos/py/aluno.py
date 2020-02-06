oprint = print;
def print(*args):
    oprint(*args,end="")
class Aluno:
    def __init__(self,nome,media):
        self.nome = nome
        self.mediaFinal = media

########################### INICIO

listaDeAlunos = [
    Aluno("Marcos",5),   #aprovado
    Aluno("Ana",3.5),    #recuperação
    Aluno("William",9.5),#aprovado
    Aluno("Márcia",2.5), #reprovado
    ]

for aluno in listaDeAlunos:
    if aluno.mediaFinal >= 5:
        print(aluno.nome,"está Aprovado")
    else:
        if aluno.mediaFinal >= 3:
            print(aluno.nome, "está de Recuperação")
        else:
            print(aluno.nome, "está Reprovado")
    print("\n")

############################ FIM
