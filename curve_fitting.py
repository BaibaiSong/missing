# -*- coding=utf-8 -*-
__author__ = 'Song'
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.patches import Circle
import time
import math
import operator

class CurveFitting(object):
    """
        最小二乘实现抛物线拟合
    """
    def __init__(self):
        self.data_path = "parabolic"
        self.order = 2  # 阶数
        self.xarray = []
        self.yarray = []

    def load_data(self):
        with open(self.data_path, 'rU') as txt:
            for line in txt:
                content = line.strip().split(',')
                self.xarray.append(eval(content[0]))
                self.yarray.append(eval(content[1]))
        print self.xarray
        print self.yarray

    def create_data(self):
        x = np.arange(-1,1,0.1)
        y = x*x - 2*x
        #生成的曲线上的各个点偏移一下，并放入到xa,ya中去
        xr=[];yr=[];i = 0
        for xx in x:
            yy=y[i]
            d=float(random.randint(90, 110))/100
            i+=1
            self.xarray.append(xx*d)
            self.yarray.append(yy*d)

    def setdata(self, xarray, yarray):
        self.xarray = xarray[:]
        self.yarray = yarray[:]

    def get_coefficient_matrix(self):
        """
        得到系数矩阵
        :return:
        """
        all_power_sum = []
        for i in range(0, 2*self.order+1):
            sum = 0
            for j in range(0, len(self.xarray)):
                sum += (self.xarray[j]**i)
            all_power_sum.append(sum)
        coefficient_matrix = []
        for i in range(0, self.order+1):
            row = all_power_sum[i: i+self.order+1]
            coefficient_matrix.append(row)
        return np.array(coefficient_matrix)

    def get_ymatrix(self):
        """
        得到y值矩阵
        :return:
        """
        ymatrix = []
        for i in range(0, self.order+1):
            ty = 0.0
            for j in range(0, len(self.yarray)):
                ty += self.yarray[j]*(self.xarray[j]**i)
            ymatrix.append(ty)
        return np.array(ymatrix)

    def get_result_matrix(self):
        """
        得到要求的参数矩阵
        :return:
        """
        coefficient_matrix = self.get_coefficient_matrix()
        ymatrix = self.get_ymatrix()
        # print 'coefficient_matrix:'
        # print coefficient_matrix
        # print 'ymatrix:'
        # print ymatrix
        result_matrix = np.linalg.solve(coefficient_matrix, ymatrix)
        return result_matrix

    def show(self):
        x = np.arange(5)
        # plt.plot(self.xarray, self.yarray, color='g', linestyle='', marker='.')
        plt.plot(x, x, color='m', linestyle='-', marker='')
        plt.show()

    def fun(self, result_matrix):
        s = self.xarray[0]
        e = self.xarray[ len(self.xarray)-1 ]
        # step = float(e-s)/20.0
        # if step == 0:
        #     step = 0.1
        # print s,e,step
        step = 20
        xarray = np.linspace(s, e, step)
        yarray = []
        for ix in xarray:
            ty = 0.0
            for j in range(0, self.order+1):
                ta = ix**j
                ta *= result_matrix[j]
                ty += ta
            yarray.append(ty)
        return xarray, yarray

    def get_fun(self):
        """
        得到多项式表达式
        :param result_matrix:
        :return:
        """
        result_matrix = self.get_result_matrix()
        # print result_matrix
        coefficient = np.array(result_matrix[::-1])
        # print 'coefficient:', coefficient
        func = np.poly1d(coefficient.astype(float))
        dfunc1 = func.deriv(m=1)
        return func, dfunc1

    def test(self):
        # 画圆
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # circle = Circle(xy=(1.0,1.0), radius=0.5, alpha=0.5)
        # ax.add_patch(circle)
        # ax.set_xlim(-4,4)
        # ax.set_ylim(-4,4)

        self.load_data()
        # self.create_data()
        func, dfunc1 = self.get_fun()
        func2 = np.poly1d(np.array([-1,0,2]).astype(float))

        x = np.linspace(-1, 1, 20)
        y = func(x)
        y1 = dfunc1(x)
        y2 = func2(x)
        plt.plot(x, y, 'r-', x, y2, 'b-')
        plt.show()
        # time.sleep(2)
        # ax.visible = False
        # plt.show()
        print 'ok'
        # fig.delaxes(ax)
        # plt.draw()
        # self.show()

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
    data_path1 = "para"
    data_path2 = "parabolic"
    xarray1, yarray1 = load_data(data_path1)
    xarray2, yarray2 = load_data(data_path2)
    print xarray1
    print yarray1
    print xarray2
    print yarray2
    return xarray1, yarray1, xarray2, yarray2

def

def Newton_root(func, x_guess): # 牛顿迭代求解
    INF = 1e-4
    dfunc = func.deriv(m=1)
    for i in np.arange(100):
        y = func(x_guess)
        dy = dfunc(x_guess)
        if dy == 0.0:
            return False, None
        else:
            x_next = x_guess - float(y/dy)
            if abs(x_next-x_guess)<INF:
                return True, x_next
            else:
                x_guess = x_next
    return False, None

def vertical_line(current_xpoint, current_ypoint, dfunc1, func2):
    """
    得到垂线函数
    :param current_xpoint:
    :param current_ypoint:
    :return:
    """
    slope = dfunc1(current_xpoint)
    print 'slope: ', slope
    ans = []
    if slope == 0:  # 导数的值为0，曲线切线为水平，垂线为竖直
        ans.append( func2(current_xpoint) )
        return ans
    else:
        k = float(1.0) / slope
        b = current_ypoint + k * current_xpoint
        print 'k and b :', -k, b
        vl = np.poly1d(np.array([-k, b]).astype(float))
        fun = func2 - vl
        # ans = np.roots(fun)
        # return True, ans
        # ------------------备用--------------------------
        ans = []
        for i in np.linspace(-1,1,30):
            flag, tans = Newton_root(fun, i) #牛顿迭代求解
            if flag is True:
                ans.append(tans)
        return flag, np.array(ans)

def get_crossover_point(current_xpoint, current_ypoint, dfunc1, func2):
    flag, solvtion = vertical_line(current_xpoint, current_ypoint, dfunc1, func2)
    print 'solvtion:', solvtion
    if flag is False:
        return False, solvtion
    if type(solvtion[0]) is np.float64:
        return True, solvtion # 注意根的个数问题
    else:
        return False, solvtion

def distance(xp1, yp1, xp2, yp2):
    return math.sqrt( (xp1-xp2)*(xp1-xp2) + (yp1-yp2)*(yp1-yp2) )

def search(func, range_left, range_right, xp, yp):
    xlist = np.linspace(range_left, range_right, 100)
    ylist = func(xlist)
    tdis = []
    for index in range(len(xlist)):
        d = distance(xp, yp, xlist[index], ylist[index])
        tdis.append(d)
    dis = np.array(tdis)
    rindex = np.argmin(dis)
    return xlist[rindex], ylist[rindex], dis[rindex]
def adjust(func1, func2, x_list, y_list, d_list, xrange_left, xrange_right, current_x, current_y):
    INF = 1.0e-5
    good_x = []
    good_y = []
    good_d = []
    for index in range(len(x_list)):
        y1 = func1(x_list[index])
        y2 = func2(x_list[index])
        if (y1<=y_list[index] and y_list[index]<=y2) or (y2<=y_list[index] and y_list[index]<=y1):
            good_x.append(x_list[index])
            good_y.append(y_list[index])
            good_d.append(d_list[index])
    dfunc1 = func1.deriv(m=1)
    for index in range(len(good_x)):
        r = good_d[index]/2.0
        croxp,croyp,tr = search(func2, good_x[index]-r, good_x[index]+r, good_x[index], good_y[index])#tr是直径还是半径？
        next_x = 0.0
        next_y = 0.0
        while abs(r-tr)>INF:
            tr = (r+tr)/2.0
            tangent_k = dfunc1(current_x)
            if tangent_k == 0:
                tempdis = []
                temp_y = [tr, -tr]
                tempdis.append(distance(current_x, temp_y[0], croxp, croyp))
                tempdis.append(distance(current_x, temp_y[1], croxp, croyp))
                tadis = np.array(tempdis)
                minindex = np.argmin(tadis)
                next_x = current_x
                next_y = temp_y[minindex]
            else :
                vertical_k = float(-1.0)/tangent_k
                b = current_y - vertical_k * current_x
                a0 = (b-current_y)*(b-current_y)+current_x*current_x-tr*tr
                a1 = 2.0*(vertical_k*b-vertical_k*current_y-current_x)
                a2 = 1.0 + vertical_k*vertical_k
                vertical_fun = np.poly1d(np.array([a2, a1, a0]).astype(float))

def test():
    cf = CurveFitting()
    cf.load_data()
    func, dfunc1 = cf.get_fun()
    # func2 = np.poly1d(np.array([-1,0,2]).astype(float))
    # func2 = np.poly1d(np.array([-0.44, -0.83, 1.42]).astype(float))
    func2 = np.poly1d(np.array([-0.16, -0.45, 0.16, 2.58]).astype(float))
    # func2 = np.poly1d(np.array([-4.0, 1.35, 4.9, -0.34, 2.02]).astype(float))
    xarray = np.linspace(-1, 1, 40)
    yarray = func(xarray)
    capabilities_x = [] #候选的中心点x
    capabilities_y = [] #候选的中心点y
    capabilities_d = [] #候选的相应直径
    for i in np.arange(len(xarray)):
        flag, cross_xpoint = get_crossover_point(xarray[i], yarray[i], dfunc1, func2) #得到交点x值
        dis = []
        if flag is False:
            continue
        else:
            for index in np.arange(len(cross_xpoint)): #有一个或多个交点
                tempy = func2(cross_xpoint[index]) #得到每个交点的y值
                tempdis = distance(xarray[i], yarray[i], cross_xpoint[index], tempy) #计算当前点与交点的距离
                dis.append(tempdis)
            disarray = np.array(dis)
            minindex = np.argmin(disarray) #选取距离最小的交点，*****最小的点应该只有一个*****
            x = (cross_xpoint[ minindex ] + xarray[i])/2.0
            y = (func2(cross_xpoint[ minindex ]) + yarray[i])/2.0
            print x, y, disarray[minindex]
            capabilities_x.append(x)
            capabilities_y.append(y)
            capabilities_d.append(disarray[minindex])
            # adjust(func, func2, capabilities_x, capabilities_y, capabilities_d, -1, 1, xarray[i], yarray[i])

    plt.plot(xarray, yarray, 'r-')
    xb = np.linspace(-1.5, 1.2, 30)
    plt.plot(xb, func2(xb), 'b-')
    plt.plot(capabilities_x, capabilities_y, 'mo')
    plt.show()


if __name__ == '__main__':
    # cf = CurveFitting()
    # cf.load_data()
    # cf.get_result_matrix()
    # cf.get_fun()
    test()