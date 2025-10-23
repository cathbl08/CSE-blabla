# search for libraray like bla.cpython-312-darwin.so in the build directory:
# import sys
# sys.path.append('/Users/joachim/texjs/lva/ws2324/ScientificComputing/ASC-bla/build')
# from bla import Vector

# import from the installed ASCsoft package:
from ASCsoft.bla import Vector
from ASCsoft.bla import Matrix

x = Vector(3)
y = Vector(3)

for i in range(len(x)):
    x[i] = i
y[:] = 2    

print ("x =", x)
print ("y =", y)
print ("x+3*y =", x+3*y)


x = Vector(10)
x[0:] = 1
print (x)

x[3:7] = 2
print (x)

x[0:10:2] = 3
print (x)

a = Matrix(2,3)
b = Matrix(2,3)

for i in range(2):
    for j in range(3):
        a[i,j] = i+j
        b[i,j] = 10*(i+j)

print ("a =", a)
print ("b =", b)
print ("a + b =", a + b)


print ("3* a =", 3*a) 

# this we dont have implemented yet
# print ("a[:,1] =", a[:,1]) 

a[1,1]= 3
print (a[1,1])

print (a[0,0])

print(a)

