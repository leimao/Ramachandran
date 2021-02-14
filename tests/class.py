# https://stackoverflow.com/questions/366422/what-is-the-pythonic-way-to-avoid-default-parameters-that-are-empty-lists

class A(object):

    def __init__(self, lst = []):
        self.lst = lst

lst = [1,2,3]
a = A()
a.lst.extend(lst)
print(a.lst)

print("-------------")
lst = [1,2,3]
a = A()
print(a.lst)
a.lst.extend(lst)
print(a.lst)

print("-------------")

lst = [1,2,3]
a = A([])
a.lst.extend(lst)
print(a.lst)

print("-------------")
lst = [1,2,3]
a = A()
print(a.lst)
a.lst.extend(lst)
print(a.lst)