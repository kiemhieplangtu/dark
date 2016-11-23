import operator
l = [15.0, 18, 2, 36, 12, 7.2, 8.8, 5, 6, 10.5, 6.0, 3.0, 2.0, 4.0, 9.0, 7.5, 8.3, 4.7]
print len(l)
print map(operator.sub, l[0:17], l[1:18])
res = map(operator.sub, l[0:17], l[1:18])
print len(res)
print sum(res)/float(len(res))

my_list = range(1, 17)
print my_list[10:16]
print sum(my_list[10:16])