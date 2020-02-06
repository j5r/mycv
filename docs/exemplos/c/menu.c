#include <stdio.h>


void main(){
int opcao;
do{
    printf("MENU--------------\n");
    printf("1 - listar produtos\n");
    printf("2 - cadastrar produtos\n");
    printf("3 - editar produtos\n");
    printf("4 - deletar produtos\n");
    printf("0 - sair\n");
    printf("\t>>> \n");

    scanf("%d", &opcao);

    switch(opcao){
        case 1:
            printf("\t>>> Listando os produtos. . .\n");
            break;
        case 2:
            printf("\t>>> Vamos cadastrar um novo produto. . .\n");
            break;
        case 3:
            printf("\t>>> Vamos escolher um produto para editar. . .\n");
            break;
        case 4:
            printf("\t>>> Vamos escolher um produto para deletar. . .\n");
            break;
        case 0:
            printf("\t>>> Obrigado. Ate breve!\n");
            break;
        default:
            printf("\t>>> Opcao invalida!\n");
        }
    }while(opcao != 0);
}
