from ASCsoft.bla import *

x = Vector(5)
y = Vector(5)

for i in range(len(x)):
    x[i] = i
y[:] = 3

print ("x+y =",x+y)

v = Vector(3)
v[:] = 7

import pickle 
f = open("file.txt", 'wb')
pickle.dump([2,"hello", v], f)
del f

f2 = open("file.txt", 'rb')
val = pickle.load(f2)
print (val)
print (val[2])