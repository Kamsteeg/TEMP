#!/usr/bin/env python


list1 = ["aaa", "bbb", "ccc"]
list2 = ["aaa"]

list3 = [["aaa1"],["aaa2"],["aaa3"],["aaa4"]]
list4 = ["aaa2"]

print list1

list1.remove(list2[0])

print list1

print list3

list3.remove(list4)

print list3
