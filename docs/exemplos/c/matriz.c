#include <stdio.h>


void main(){
    int matriz[3][5] = {
            {1,2,3,4,5},
            {6,7,8,9,10},
            {11,12,13,14,15}
        };
    int linha, coluna;

printf("Vamos imprimir uma matriz\n");
for(linha=0; linha<3; linha++){
    for(coluna=0; coluna<5; coluna++){
        printf("%d\t",matriz[linha][coluna]);
        }
    printf("\n");
    }
printf("Fim do programa\n");
}
