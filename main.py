# -*- coding=utf-8 -*-
__author__ = 'Song'

import matplotlib.pyplot as plt
import math

def test1():
    xarray = [200, 300, 450, 600]
    yarray = [300, 200, 400, 270]
    u = [0.0]
    for i in range(1, len(xarray)):
        temp = 0.0
        temp = math.sqrt((xarray[i]-xarray[i-1])**2 + (yarray[i]-yarray[i-1])**2)
        print temp
        temp += u[i-1]
        u.append(temp)
    print u
    n = len(u)

    for i in range(n):
        u.append(u[i]/u[n-1])
    print u

def test2():

    x = [2, 4, 6, 8]
    for i in range(len(x)):
        if x[i]>5:
            break

    print i
    print x[i]

if __name__ == '__main__':
    # li = []
    # li.append(None)
    # if li[0] is None:
    #     print 'ok'
    # print len(li)
    x = [1.0, 2.5, 3.0, 4.0]
    y = [1.0, 2.5, 2.0, 3.0]

    # x.pop(2)
    # y.pop(2)

    # plt.plot(x, y, 'b-', x, y, 'ro')
    # plt.show()

    test2()