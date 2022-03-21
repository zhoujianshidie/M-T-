# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 11:35:24 2022

@author: jiange
"""
import numpy as np
import math
from pylab import *
mpl.rcParams['font.sans-serif'] = ['SimHei']
'''
C_11 = int(input("请输入C11"))
C_12 = int(input("请输入C12"))
'''
K_p = 13000000
G_p = 0
C_p = np.array([[K_p, K_p, K_p, 0, 0, 0],
                [K_p, K_p, K_p, 0, 0, 0],
                [K_p, K_p, K_p, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0]],dtype = float)


class MAT:#定义一个矩阵类
    def __init__(self):
        '''
        C_11,C_12是各向同性材料的刚度张量中两个独立的分量
        C定义了一个刚度矩阵

        Parameters
        ----------
        C_11 : TYPE
            DESCRIPTION.
        C_12 : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.v = None
        self.E = None
        self.G = None
        self.k = None
        
        self.C = None
        self.S = None

    
    def _sti(self):
        '''
        获取刚度矩阵

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        v = self.v
        E = self.E
        G = self.G
        self.C = np.array(
            [[1/E,-(v/E),-(v/E),0,0,0],
             [-(v/E),1/E,-(v/E),0,0,0],
             [-(v/E),-(v/E),1/E,0,0,0],
             [0,0,0,1/G,0,0],
             [0,0,0,0,1/G,0],
             [0,0,0,0,0,1/G]],dtype = float)#此为柔度矩阵
        self.C = np.linalg.inv(self.C)#将柔度矩阵求逆得到刚度矩阵
        return self.C#返回刚度矩阵

    def _esl(self):
        '''
        获得Eshelly张量

        Returns
        -------
        None.

        '''
        k = self.k
        u = self.u_0
        if k==1:
            S_2222 = S_3333 = S_1111 = (7-5*u)/(15*(1-u))
            S_1122 = S_1133 = S_2211 = S_3311 = S_2233 = S_3322 = (5*u-1)/(15*(1-u))
            S_1212 = S_1313 = S_2323 = (4-5*u)/(15*(1-u))
        elif k>1:
            g = k/pow((k*k-1),1.5)*(k*math.sqrt(k*k-1)-math.acosh(k))
            S_1111 =(1/(2*(1-u))*(4-2*u-2/(1-k*k))+1/(2*(1-u))*(2*u-4+3/(1-k*k))*g)
            S_2222 =((-3)/(8*(1-u))*k*k/(1-k*k)+1/(4*(1-u))*(1-2*u+9/(4*(1-k*k)))*g)
            S_3333 = S_2222
            S_2233 =(1/(8*(1-u))*(1-1/(1-k*k))+1/(16*(1-u))*(8*u-4+3/(1*(1-k*k)))*g)
            S_3322 = S_2233
            S_2211 =(1/(2*(1-u))*k*k/(1-k*k)-1/(4*(1-u))*(1-2*u+3/(1*(1-k*k)))*g)
            S_3311 = S_2211
            S_1122 =(1/(2*(1-u))*(2*u-1+1/(1-k*k))+1/(4*(1-u))*(2-4*u-3/(1-k*k))*g)
            S_1133 = S_1122
            S_2323 =((-1)/(8*(1-u))*k*k/(1-k*k)+1/(16*(1-u))*(4-8*u+3/(1-k*k))*g)
            S_1212 =(1/(4*(1-u))*(1-2*u+(1+k*k)/(1-k*k))-1/(8*(1-u))*(1-2*u+3*(1+k*k)/(1-k*k))*g)
            S_1313 = S_1212
        elif k<1:
            g = k/pow((1-k*k),1.5)(math.acos(k)-k*math.sqrt(1-k*k))
            S_1111 =(1/(2*(1-u))*(4-2*u-2/(1-k*k))+1/(2*(1-u))*(2*u-4+3/(1-k*k))*g)
            S_2222 =((-3)/(8*(1-u))*k*k/(1-k*k)+1/(4*(1-u))*(1-2*u+9/(4*(1-k*k)))*g)
            S_3333 = S_2222
            S_2233 =(1/(8*(1-u))*(1-1/(1-k*k))+1/(16*(1-u))*(8*u-4+3/(1*(1-k*k)))*g)
            S_3322 = S_2233
            S_2211 =(1/(2*(1-u))*k*k/(1-k*k)-1/(4*(1-u))*(1-2*u+3/(1*(1-k*k)))*g)
            S_3311 = S_2211
            S_1122 =(1/(2*(1-u))*(2*u-1+1/(1-k*k))+1/(4*(1-u))*(2-4*u-3/(1-k*k))*g)
            S_1133 = S_1122
            S_2323 =((-1)/(8*(1-u))*k*k/(1-k*k)+1/(16*(1-u))*(4-8*u+3/(1-k*k))*g)
            S_1212 =(1/(4*(1-u))*(1-2*u+(1+k*k)/(1-k*k))-1/(8*(1-u))*(1-2*u+3*(1+k*k)/(1-k*k))*g)
            S_1313 = S_1212
        
        self.S = np.array(
            [[S_1111,S_1122,S_1133,0,0,0],
             [S_2211,S_2222,S_2233,0,0,0],
             [S_3311,S_3322,S_3333,0,0,0],
             [0,0,0,2*S_2323,0,0],
             [0,0,0,0,2*S_1212,0],
             [0,0,0,0,0,2*S_1212]],dtype = float)#此为Eshelly矩阵
        return self.S#返回Eshelly矩阵

def get_C():
    '''
    获取刚度矩阵

    Returns
    -------
    c : TYPE
        DESCRIPTION.

    '''
    Mat_C = MAT()
    #泊松比，剪切模量，体积模量
    Mat_C.v = 0.35#float(input("请输入v"'\n'))
    Mat_C.E = 3350000000#float(input("请输入E"'\n'))
    Mat_C.G = Mat_C.E/(2*(1+Mat_C.v))#float(input("请输入G"'\n'))
    c = Mat_C._sti()
    return c
    
    
def get_S():
    '''
    获取Eshelly矩阵

    Returns
    -------
    s : TYPE
        DESCRIPTION.

    '''
    Mat_S = MAT()
    Mat_S.k = 1#float(input("请输入k"'\n'))
    Mat_S.u_0 = 0.35#float(input("请输入u0"'\n'))
    s = Mat_S._esl()
    return s
    
Porosity = []
array_E = []
array_G = []
array_v = []
n = int(input("请输入共有几组数据"'\n'))
C_0 = get_C()
S_0 = get_S()
for i in range(0, n):#循环得到模量和孔隙率的关系
    print("请输入第{}组数据".format(i+1),'\n')
    c_1 = float(input("请输入孔隙体积分数"'\n'))
    c_0 = 1-c_1
    
    Porosity.append(c_1)

    #获取有效刚度矩阵
    C_ = C_0 +c_1*np.linalg.inv(np.linalg.inv((C_p-C_0))+c_0*np.dot(S_0,np.linalg.inv(C_0)))
    S_ = np.linalg.inv(C_)
    E_m1 = 1/S_[0,0]
    G_m1 = 1/S_[4,4]
    v_m1 = (0.5*E_m1)/G_m1-1
    array_E.append(E_m1)
    array_G.append(G_m1)
    array_v.append(-v_m1)
#print(array_E, array_G, array_v)
#绘制图像
plt.plot(Porosity, array_E, 'ro-', color='red', alpha=0.8, label='体积模量')
plt.plot(Porosity, array_G, 'bs-', color='blue', alpha=0.8, label='剪切模量')

'''
for x, y in zip(Porosity, array_E):
    plt.text(x, y+0.3, '%.0f' % y, ha='center', va='bottom', fontsize=10.5)
for x, y in zip(Porosity, array_G):
    plt.text(x, y+0.3, '%.0f' % y, ha='center', va='bottom', fontsize=10.5)

plt.plot(Porosity, array_v, 'g^-', color='blue', alpha=0.8, label='泊松比')
for x, y in zip(Porosity, array_v):
    plt.text(x, y+0.3, '%.0f' % y, ha='center', va='bottom', fontsize=10.5)
'''
    
plt.legend(loc="upper right")
plt.xlabel('孔隙率')

plt.show()

