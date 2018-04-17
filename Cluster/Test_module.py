from module_Flo import directory

directory = directory("/home/charriat/Documents/Script")
print(directory.path)
print(directory.listAll)
print(directory)
print(directory.listExt('directory'))
print(directory.listExt('fasta'))
print(directory.listExt('allfastq'))
