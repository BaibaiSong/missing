# -*- coding=utf-8 -*-
__author__ = 'Song'
import numpy as np
import matplotlib.pyplot as plt
import operator
import math
from matplotlib.patches import Circle

class CubicSpline(object):
    """

    """
    def __init__(self):
        self.data_path = "data"
        self.order = 2
        self.xarray = []
        self.yarray = []
        self.M = []

    def load_data(self):
        with open(self.data_path, 'rU') as txt:
            for line in txt:
                content = line.strip().split(',')
                self.xarray.append(eval(content[0]))
                self.yarray.append(eval(content[1]))
        print self.xarray
        print self.yarray

    def setdata(self, xarray, yarray):
        self.xarray = xarray[:]
        self.yarray = yarray[:]

    def tdma(self, a, b, c, y):

        if b[0] == 0:
            return False
        c[0] = c[0]/b[0]
        y[0] = y[0]/b[0]
        for i in range(1, len(b)):
            temp = float(b[i]-a[i]*c[i-1])
            if temp == 0.0:
                return False
            temp = float(1.0/temp)
            c[i] *= temp
            y[i] = (y[i]-a[i]*y[i-1])*temp

        for i in range(len(b)-2, -1, -1):
            y[i] = y[i]-c[i]*y[i+1]
        return True

    def spline(self):

        d1 = 0.0  # bound var
        d2 = 0.0  # bound var
        boundType = 2  # bound Type 2: nature

        n = len(self.xarray)
        a = np.zeros(n)
        b = np.zeros(n)
        c = np.zeros(n)

        self.M = np.zeros(n)
        b[0] = 2.0
        b[n-1] = 2.0
        if boundType == 1:
            c[0] = 1.0
            a[n-1] = 1.0
            h0 = float(self.xarray[1] - self.xarray[0])
            self.M[0] = 6.0*((self.yarray[1]-self.yarray[0])/h0 - d1)/h0
            hn = float(self.xarray[n-1] - self.xarray[n-2])
            self.M[n-1] = 6.0*(d2 - (self.yarray[n-1] - self.yarray[n-2])/hn)/hn
        else:
            c[0] = 0.0
            a[n-1] = 0.0
            self.M[0] = 2*d1
            self.M[n-1] = 2*d2
        for i in range(1,n-1):
            hi_1 = float(self.xarray[i] - self.xarray[i-1])
            hi = float(self.xarray[i+1] - self.xarray[i])
            a[i] = hi_1/(hi_1+hi)
            b[i] = 2.0
            c[i] = 1-a[i]
            self.M[i] = 6.0*((self.yarray[i+1]-self.yarray[i])/hi - (self.yarray[i]-self.yarray[i-1])/hi_1)/(hi+hi_1)
        # print a
        # print b
        # print c
        # print self.M
        self.tdma(a,b,c,self.M)
        print self.M

    def get_value(self, x):
        # print self.M
        n = len(self.xarray)
        if x<self.xarray[0] or self.xarray[n-1]<x:
            return None
        index = 0
        for i in range(n-1):
            if self.xarray[i]<=x and x<=self.xarray[i+1]:
                index = i
                break
        hi = self.xarray[index+1]-self.xarray[index]
        xi_1 = self.xarray[index+1]-x
        xi = x-self.xarray[index]
        Ai = (self.yarray[index+1]-self.yarray[index])/hi-(self.M[index+1]-self.M[index])*hi/(6.0)
        Bi = self.yarray[index+1]-(self.M[index+1]*hi*hi)/6.0-Ai*self.xarray[index+1]
        y = (xi_1**3) * self.M[index]/(6.0*hi)+(xi**3) * self.M[index+1]/(6.0*hi)+Ai*x+Bi
        dy = -xi_1*xi_1*self.M[index]/(2.0*hi)+xi*xi*self.M[index+1]/(2.0*hi)+Ai
        # temp1 = self.M[index]/(6.0*hi)
        # temp2 = self.M[index+1]/(6.0*hi)
        # temp3 = (self.yarray[index]-self.M[index]*hi*hi/6.0)/hi
        # temp4 = (self.yarray[index+1]-self.M[index+1]*hi*hi/6.0)/hi
        # y = temp1*(self.xarray[index+1]-x)**3+ \
        #     temp2*(x-self.xarray[index])**3+ \
        #     temp3*(self.xarray[index+1]-x)+ \
        #     temp4*(x-self.xarray[index])
        # print temp1,temp2,temp3,temp4
        return y, dy

    def get_newton_value(self, x, k, b):
        n = len(self.xarray)
        if x<self.xarray[0] or self.xarray[n-1]<x:
            return None
        index = 0
        for i in range(n-1):
            if self.xarray[i]<=x and x<=self.xarray[i+1]:
                index = i
                break
        hi = self.xarray[index+1]-self.xarray[index]
        xi_1 = self.xarray[index+1]-x
        xi = x-self.xarray[index]
        Ai = (self.yarray[index+1]-self.yarray[index])/hi-(self.M[index+1]-self.M[index])*hi/(6.0)
        Bi = self.yarray[index+1]-(self.M[index+1]*hi*hi)/6.0-Ai*self.xarray[index+1]
        y = (xi_1**3) * self.M[index]/(6.0*hi)+(xi**3) * self.M[index+1]/(6.0*hi)+Ai*x+Bi-k*x-b
        dy = -xi_1*xi_1*self.M[index]/(2.0*hi)+xi*xi*self.M[index+1]/(2.0*hi)+Ai-k
        return y, dy

    def getfunc(self):
        pass

    def test(self):
        self.load_data()
        self.spline()
        x = np.linspace(27.7, 30, 50)
        # x = np.linspace(0.0,8.0,50)
        # x = np.linspace(0,9,100)
        y = []
        for it in x:
            ty, tdy = self.get_value(it)
            y.append(ty)
        # print x
        # print y
        print self.M
        plt.plot(self.xarray,self.yarray,'ro',x, y, 'b-')
        plt.show()

def sort_by_attr(seq, attr):
    return sorted(seq, key=operator.attrgetter(attr))

class Point(object):
    """
    表示点
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def set(self, x, y):
        self.x = x
        self.y = y

    def get(self):
        return self.x, self.y

def load_data(data_path):

    pointlist = []
    with open(data_path, 'rU') as txt:
        for line in txt:
            content = line.strip().split(',')
            point = Point(eval(content[0]), eval(content[1]))
            pointlist.append(point)
    plist = sort_by_attr(pointlist, 'x')
    xarray = []
    yarray = []
    for item in plist:
        x, y = item.get()
        xarray.append(x)
        yarray.append(y)
    return xarray, yarray

def get_data():
    data_path1 = "data1"
    data_path2 = "data2"
    xarray1, yarray1 = load_data(data_path1)
    xarray2, yarray2 = load_data(data_path2)
    print xarray1
    print yarray1
    print xarray2
    print yarray2
    return xarray1, yarray1, xarray2, yarray2

def newton(cspline, k, b, guess_x, scope_left, scope_right):
    INF = 1.0e-3  # 精度控制
    for i in range(100):  # 迭代次数
        if guess_x<scope_left or guess_x>scope_right:
            break
        # print 'guess_x, k, b:',guess_x, k, b, scope_left, scope_right
        y, dy = cspline.get_newton_value(guess_x, k, b)
        if dy == 0.0:
            return False, None
        else:
            x_next = guess_x - float(y/dy)
            if abs(x_next-guess_x)<INF:
                return True, x_next
            else:
                guess_x = x_next
    return False, None

def get_cross_point(x, yarray_i, cspline1, cspline2, scope_left, scope_right):

    y, dy = cspline1.get_value(x)
    corss_pointx = []
    corss_pointy = []
    if dy == 0.0:  # 如果切线斜率为0，则垂线的竖直
        y2, dy2 = cspline2.get_value(x)
        corss_pointx.append(x)
        corss_pointy.append(y2)
        return corss_pointx, corss_pointy
    else:  # 如果切线斜率不为0，则垂线的斜率为-1/dy
        k = float(-1.0/dy)
        b = yarray_i - k*x   # 得到垂线的斜率和b，即得到垂线方程
        # print 'k,b:', k, b
        tempx = scope_left-1.0
        for guess_x in np.linspace(scope_left, scope_right, 10):
            flag, cpx = newton(cspline2, k, b, guess_x, scope_left, scope_right)  # 可能有多个交点
            if flag is True and (scope_left<=cpx and cpx<=scope_right) and abs(tempx-cpx) > 1.0e-2:
                tempx = cpx
                cpy = k*cpx+b
                corss_pointx.append(cpx)
                corss_pointy.append(cpy)
        # print corss_pointx
        # print corss_pointy
        return corss_pointx, corss_pointy

def test1():
    cs = CubicSpline()
    a = [0, 0.231, 0.5, 1.0]
    b = [2.0, 2.0, 2.0, 2.0]
    c = [1.0, 0.769, 0.5, 0]
    y = [-46.6666, -4.00002, -2.7, -17.4]
    print cs.tdma(a,b,c,y)
    print y
    # print range(len(b)-1, -1, -1)

def newton_root(func, xguess):
    INF = 1.0e-3
    dfunc = func.deriv(m=1)
    for i in range(100):#迭代次数
        y = func(xguess)
        dy = dfunc(xguess)
        if dy == 0.0:
            return False, None
        xnext = xguess - float(y/dy)
        if abs(xnext-xguess)<INF:
            return True, xnext
        else:
            xguess = xnext
    return False, None

def find(cs2, centpointx, centpointy, left, right, scope_left, scope_right):
    """
    在第二条曲线上搜索新点

    :param cs2:
    :param centpointx:
    :param centpointy:
    :param left:
    :param right:
    :param scope_left:
    :param scope_right:
    :return:新半径，第二条曲线上的新点
    """
    if left<scope_left:
        left = scope_left
    if right>scope_right:
        right = scope_right
    mindis = 999.99
    rx = 0.0
    ry = 0.0
    for x in np.linspace(left, right, 100):
        y, dy = cs2.get_value(x)
        dis = distance(x, y, centpointx, centpointy)
        if dis<mindis:
            mindis = dis
            rx = x
            ry = y
    return mindis, rx, ry

def between(spline1, spline2, x1, y1):
    y1line1, dyline1 = spline1.get_value(x1)
    y2line1, dy2line1 = spline2.get_value(x1)
    if y1line1<=y1<=y2line1 or y2line1<=y1<=y1line1:
        return True
    else:
        return False

def getmindis(x1, y1, x2, y2, x, y):
    dis1 = distance(x1, y1, x, y)
    dis2 = distance(x2, y2, x, y)
    if dis1>dis2:
        return x2, y2
    else:
        return x1, y1

def adjust(xarray, yarray, cs1, cs2, cap_x, cap_y, cap_d):
    n = len(xarray)
    scope_left = xarray[0]
    scope_right = xarray[n-1]
    for i in range(n):
        centpx = cap_x[i]
        centpy = cap_y[i]
        y1, dy1 = cs1.get_value(xarray[i])
        current_r = cap_d[i]/2.0
        left = cap_x[i] - current_r
        right = cap_x[i] + current_r
        new_r, rx, ry = find(cs2, cap_x[i], cap_y[i], left, right, scope_left, scope_right)

        if dy1 == 0:  # 垂线竖直的时候
            print 'new_r,(x,y):', new_r, rx, ry
            while abs(current_r - new_r) > 1.0e-2:
                current_r = (current_r + new_r)/2.0   # 搜索半径和原来的半径的一半当做新的半径
                print 'current_r:', current_r
                centpx = xarray[i]   # 垂线竖直，故中心点的X坐标为切点坐标
                temppy1 = yarray[i] - current_r  # Y坐标为切点坐标的y的半径长度
                temppy2 = yarray[i] + current_r
                centpx, centpy = getmindis(centpx, temppy1, centpx, temppy2, rx, ry)
                # print 'centpx, centpy:', centpx, centpy
                left = centpx - current_r  # 再次搜索的左右范围
                right = centpx + current_r
                new_r, rx, ry = find(cs2, centpx, centpy, left, right, scope_left, scope_right)
            cap_x[i] = centpx
            cap_y[i] = centpy
            cap_d[i] = current_r * 2.0  # distance(xarray[i], yarray[i], rx, ry)  # 直径为曲线上两点之间的距离
            continue
        # 当垂线斜率存在
        k = float(-1.0/dy1)
        b = y1 - k*xarray[i]
        # print 'k,b:---------x', k, b, xarray[i]
        # current_r = cap_d[i]/2.0
        # left = cap_x[i] - current_r
        # right = cap_x[i] + current_r
        # new_r, rx, ry = find(cs2, cap_x[i], cap_y[i], left, right, scope_left, scope_right)
        mark = True
        while abs(current_r-new_r)>1.0e-2:
            current_r = (current_r+new_r)/2.0
            a0 = (b-y1)*(b-y1)+xarray[i]*xarray[i]-current_r*current_r
            a1 = 2.0*(k*b-k*y1-xarray[i])
            a2 = 1.0+k*k
            func = np.poly1d(np.array([a2, a1, a0]).astype(float))  # 得到联立的方程，用牛顿迭代求解

            flag = False
            tempx1 = 0.0
            for guess_x in np.linspace(xarray[i]-current_r, xarray[i]+current_r, 10):
                flag, tempx1 = newton_root(func, guess_x)
                if flag is True:
                    break
            if flag is True:  # 得到一个距离切点新半径个距离的解，那么另一个解关于切点对称
                tempx2 = 2.0*xarray[i] - tempx1  # 求得对称点
                tempy2 = k*tempx2 + b  # 两点都在垂线上
                tempy1 = k*tempx1 + b
                centpx, centpy = getmindis(tempx1, tempy1, tempx2, tempy2, rx, ry)
                left = centpx - current_r
                right = centpx + current_r
                new_r, rx, ry = find(cs2, centpx, centpy, left, right, scope_left, scope_right)
            else:
                mark = False
                break
        if mark is True:
            cap_x[i] = centpx
            cap_y[i] = centpy
            cap_d[i] = current_r * 2.0   # distance(xarray[i], yarray[i], rx, ry)
        else:
            # cap_x.pop(i)
            # cap_y.pop(i)
            # cap_d.pop(i)
            print '---i---', i
            cap_x[i] = None
            cap_y[i] = None
            cap_d[i] = None


def show1(xarray1, yarray1, cs1, xarray2, yarray2, cs2, cx, cy, cd):
    n = len(xarray1)
    scope_left = xarray1[0]
    scope_right = xarray1[n-1]
    x = np.linspace(scope_left, scope_right, 50)
    # cs1 = CubicSpline()
    # cs1.setdata(xarray1, yarray1)
    # cs1.spline()
    y1 = []

    # cs2 = CubicSpline()
    # cs2.setdata(xarray2, yarray2)
    # cs2.spline()
    y2 = []

    for i in x:
        y, dy = cs1.get_value(i)
        ty, tdy = cs2.get_value(i)
        y1.append(y)
        y2.append(ty)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    m = len(cx)
    for i in range(m):
        circle = Circle(xy=(cx[i], cy[i]), radius=cd[i]/2.0, alpha=0.2)
        ax.add_patch(circle)

    plt.plot(xarray1, yarray1, 'ro')
    plt.plot(xarray2, yarray2, 'ro')
    plt.plot(x, y1, 'b-')
    plt.plot(x, y2, 'b-')
    plt.plot(cx, cy, 'bo')
    plt.show()

def distance(x1, y1, x2, y2):
    x = float(x1 - x2)
    y = float(y1 - y2)
    return math.sqrt( x*x+y*y )

def main():
    xarray1, yarray1, xarray2, yarray2 = get_data()

    n = len(xarray1)
    scope_left = xarray1[0]
    scope_right = xarray1[n-1]
    # xarr = np.linspace(scope_left, scope_right, 15)
    # yarr = []
    xarr = xarray1[:]
    yarr = yarray1[:]
    cs1 = CubicSpline()
    cs1.setdata(xarray1, yarray1)
    cs1.spline()

    cs2 = CubicSpline()
    cs2.setdata(xarray2, yarray2)
    cs2.spline()

    capabilities_x = []  # 候选的中心点x
    capabilities_y = []  # 候选的中心点y
    capabilities_d = []  # 候选的相应直径
    last_d = abs(float(yarray2[0]-yarray1[0]))  # 上一个直径
    m = len(xarr)
    for i in range(m):
        tempy, dyarr = cs1.get_value(xarr[i])
        yarr.append(tempy)
        cross_x, cross_y = get_cross_point(xarr[i], yarr[i], cs1, cs2, scope_left, scope_right)
        if len(cross_x) == 0:
            capabilities_x.append(xarr[i])
            capabilities_y.append(yarr[i])
            capabilities_d.append(last_d)
        else:   # 选取离当前点最近的交点
            tempdislist = []
            for j in range(len(cross_x)):
                dis = distance(xarr[i], yarr[i], cross_x[j], cross_y[j])
                tempdislist.append(dis)
            disarray = np.array(tempdislist)
            minindex = np.argmin(disarray)

            centpointx = float(xarr[i]+cross_x[minindex])/2.0
            centpointy = float(yarr[i]+cross_y[minindex])/2.0

            bdy1, dbdy1 = cs1.get_value(centpointx)
            bdy2, dbdy2 = cs2.get_value(centpointx)
            if (bdy1 <= centpointy <= bdy2) or (bdy2 <= centpointy <= bdy1):  # 交点在两条曲线之间
                capabilities_x.append(centpointx)
                capabilities_y.append(centpointy)
                capabilities_d.append(tempdislist[minindex])
                last_d = tempdislist[minindex]
            else:
                capabilities_x.append(xarr[i])
                capabilities_y.append(yarr[i])
                capabilities_d.append(last_d)
            # print '---------cross point---------'
            # print xarray1[i], yarray1[i]
            # print centpointx, centpointy
            # print cross_x[minindex], cross_y[minindex]
    print xarr
    print yarr
    print capabilities_x
    print capabilities_y
    print capabilities_d
    adjust(xarr, yarr, cs1, cs2, capabilities_x, capabilities_y, capabilities_d)
    print '---------------------result-----------------------'
    print xarr
    print yarr
    print capabilities_x
    print capabilities_y
    print capabilities_d
    show1(xarray1, yarray1, cs1, xarray2, yarray2, cs2, capabilities_x, capabilities_y, capabilities_d)

if __name__ == '__main__':
    # cs = CubicSpline()
    # cs.test()
    main()