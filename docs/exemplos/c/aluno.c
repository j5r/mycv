#include <stdio.h>
//#include <stdlib.h>
#include <string.h>

struct Aluno{
    char nome[30];
    float mediaFinal;
};

struct listaAlunos{
    struct Aluno alunos[50];
    int num_alunos;
};



void main(){
int i;
struct Aluno aluno;

struct listaAlunos listaDeAlunos;
listaDeAlunos.num_alunos = 0;

strcpy(aluno.nome,"Marcos");
aluno.mediaFinal = 5.0;
listaDeAlunos.alunos[0] = aluno;
listaDeAlunos.num_alunos += 1;


strcpy(aluno.nome, "Ana");
aluno.mediaFinal = 3.5;
listaDeAlunos.alunos[1] = aluno;
listaDeAlunos.num_alunos += 1;

strcpy(aluno.nome, "William");
aluno.mediaFinal = 9.5;
listaDeAlunos.alunos[2] = aluno;
listaDeAlunos.num_alunos += 1;


strcpy(aluno.nome, "Marcia");
aluno.mediaFinal = 2.5;
listaDeAlunos.alunos[3] = aluno;
listaDeAlunos.num_alunos += 1;




for(i=0; i<listaDeAlunos.num_alunos; i++){
    aluno = listaDeAlunos.alunos[i];
    if(aluno.mediaFinal >= 5.0)
        printf("%s está aprovado.\n", aluno.nome);
    else if(aluno.mediaFinal >= 3.0)
        printf("%s está de recuperação.\n", aluno.nome);
    else
        printf("%s está reprovado.\n", aluno.nome);
    }

}
