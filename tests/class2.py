# https://stackoverflow.com/questions/366422/what-is-the-pythonic-way-to-avoid-default-parameters-that-are-empty-lists
from typing import List, Optional

class A(object):

    def __init__(self, lst: Optional[List[int]] = None):
        
        if lst is None:
            self.lst = [0]
        else:
            self.lst = lst

print("-------------")
lst = [1,2,3]
a = A()
print(a.lst)
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
a = A([0])
print(a.lst)
a.lst.extend(lst)
print(a.lst)

print("-------------")
lst = [1,2,3]
a = A()
print(a.lst)
a.lst.extend(lst)
print(a.lst)

print("-------------")