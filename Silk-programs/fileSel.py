import random
FileSelection1 = random.randrange(1,12)
FileSelection2 = 0
while (FileSelection2 == FileSelection1) or (FileSelection2 == 0) :
        FileSelection2 = random.randrange(1,12)
print(FileSelection1, FileSelection2)  
