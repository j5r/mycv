#include <stdio.h>


void main(){
int contador = 0;
int soma = 0;
goto secao2;

secao1:
    printf("%d\n",soma);
    printf("fim do programa\n");
goto secaoFim;

secao2:
    contador = contador + 1;
    soma = soma + contador;
    if(contador < 5){
        goto secao2;
    }
    else{
        goto secao1;
    }

secaoFim:
    printf("Fim :)\n");
}
