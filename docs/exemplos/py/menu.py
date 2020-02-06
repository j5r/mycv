while True:
    print("MENU--------------")
    print("1 - listar produtos")
    print("2 - cadastrar produtos")
    print("3 - editar produtos")
    print("4 - deletar produtos")
    print("0 - sair")
    opção = int(input("\t>>> "))

    if opção == 1:
        print("\t>>> ", "Listando os produtos. . .")
    elif opção == 2:
        print("\t>>> ", "Vamos cadastrar um novo produto. . .")
    elif opção == 3:
        print("\t>>> ", "Vamos escolher um produto para editar. . .")
    elif opção == 4:
        print("\t>>> ", "Vamos escolher um produto para deletar. . .")
    elif opção == 0:
        print("\t>>> ", "Obrigado. Até breve!")
        break
    else:
        print("\t>>> ", "Opção inválida!")
