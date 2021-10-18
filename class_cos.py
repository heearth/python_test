#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import datetime


# In[2]:


'''params={
    'sigma':1,          #波动率的波动率
    'gamma':1.5,     #非仿射系数
    'theta':0.09,       #波动率的长期均值
    'kappa':2,          #波动率均值恢复速度
    'v_0':0.09,          #初始波动率
    'r':0.05,           #无风险利率
    'rho':-0.3,          #相关系数    
    'S_0':100,            #股票价格
    'K':100,            #行权价
    'tau':1,            #到期时间
    'L':10,             #RR和RL里面的东西，建议取10
    'N':200             #级数展开项数，200项够了
}'''


# In[20]:


class cos_fourier():
      
    #########————————————————#########
    def __init__(self,params):
        #参数赋值
        self.sigma=params['sigma']          #波动率的波动率
        self.gamma=params['gamma']         #非仿射系数
        self.theta=params['theta']          #波动率的长期均值
        self.kappa=params['kappa']       #波动率均值恢复速度
        self.v_0=params['v_0']          #初始波动率
        self.r=params['r']          #无风险利率
        self.rho=params['rho']          #相关系数
        self.S_0=params['S_0']          #股票价格
        self.K=params['K']            #行权价
        self.tau=params['tau']           #到期时间
        self.N=params['N']             #级数展开项数，200项够了
        self.L=params['L']            #RR和L里面的东西，建议取10

        self.X=np.log(self.S_0/self.K)
        self.delta_d=0.00001 #微分的增量（步长）
        self.array_k=np.arange(1.,self.N) #k=0时在累加时特殊处理；（k是数组时，随k而动的函数返回的都是数组）        
        self.phi=self.array_k*np.pi / (self.RR()-self.LL()) #（随k而动）
        return

    # 特征函数；传入数组或者单个数
    def fun_charact(self,phi):
        '''特征函数用以替代概率密度函数，根据每一次迭代（k取不同的值）相应取不同的值（随k而动）
        #     d=(b**2 -4*a*c)**0.5 #当self.phi=0时的处理方式不同，得到的结果也不同，
        注释掉的方式为：分式上下同乘分母，避免分母=0；
        在使用的是给d添上负号
        #     B=(d-b)/(2*a)*(((1-np.exp(self.tau*d))*(d+b)/(d-b))/((d+b)/(d-b)+np.exp(self.tau*d)))
        '''
        phi=phi
        if isinstance(phi,int) and phi ==0:
            k=0
        else:
            k=np.arange(1.,self.N) 
        a=0.5*self.sigma**2 *self.gamma*self.theta**(self.gamma-1)
        b=self.rho*self.sigma*phi*(self.gamma+1)*0.5*self.theta**((self.gamma-1)/2)*1j - self.kappa
        c=-0.5*phi*1j - 0.5*phi**2
        d=-(b**2 -4*a*c)**0.5    
        h=b + d - b*np.exp(d*self.tau) + d*np.exp(d*self.tau)

        m1=0.5*(1-self.gamma)*self.rho*self.sigma*phi*self.theta**((self.gamma+1)*0.5) *1j + self.theta*self.kappa
        m2=0.5*(1-self.gamma)*self.sigma**2*self.theta**self.gamma
        A_tau=self.r*self.tau*phi*1j +             (a*m1-b*m2)*np.log(2*d/h)/a**2 +             m1*(d-b)*self.tau/(2*a) +             m2*((b-d)**2*self.tau + (4*d/h - 2)*(b+d))/(4*a**2)  
        global B_tau
        B_tau=(d-b)/(2*a)*(((1-np.exp(self.tau*d)))/(1+np.exp(self.tau*d)*(d-b)/(d+b)))
        global A_de_t #对self.tau求导
        global B_de_t#对self.tau求导
        A_de_t=0.5*self.sigma**2*self.theta**self.gamma* (1-self.gamma)*B_tau**2+                (1j*self.rho*self.sigma*phi*self.theta**((self.gamma+1)/2) * ((self.gamma-1)/2) + k*self.theta )*B_tau +1j*self.r*phi
        B_de_t=0.5*self.sigma**2*self.gamma*self.theta**(self.gamma-1)* B_tau**2+                (1j*self.rho*self.sigma*phi*self.theta**((self.gamma-1)/2) * ((self.gamma+1)/2) - k )*B_tau +1j*0.5*phi*(1j*phi-1)

        charact=np.exp( phi*self.X*1j + A_tau + B_tau*self.v_0 )
        return charact
    #########————————————————#########
    #变量n阶累积量和矩-积分上下限，与k无关
    def RR(self):
        RR=np.real(-1j*(self.fun_charact(self.delta_d)-self.fun_charact(0))/self.delta_d            +self.L*np.sqrt(-1*(self.fun_charact(self.delta_d)-2*self.fun_charact(0)            +self.fun_charact(-self.delta_d))/(self.delta_d**2)))
        return RR
    def LL(self):
        RL=np.real(-1j*(self.fun_charact(self.delta_d)-self.fun_charact(0))/self.delta_d            -self.L*np.sqrt(-1*(self.fun_charact(self.delta_d)-2*self.fun_charact(0)            +self.fun_charact(-self.delta_d))/(self.delta_d**2)))
        return RL    
    #########————————————————#########    
    
    #########————————————————#########
    #定积分运算，公式14-15
    def chi_k(self,y,iter_): #（随k而动）
        if iter_ ==0: # k=0
            chi= np.exp(y)  
        else:        
            c0=(self.RR()-self.LL())/(self.array_k*np.pi) #简化一堆
            z=(y- self.LL())/c0 #变量代换
            chi=np.exp(self.LL())*c0*(np.cos(z)*np.exp(z*c0)+np.sin(z)*np.exp(z*c0)/c0)/(c0+1/c0)
        return chi

    def psi_k(self,y,iter_): #（随k而动）
        if  iter_ ==0: # k=0
            psi= y   
        else:        
            psi=(self.RR()-self.LL())/(self.array_k*np.pi)*np.sin(self.array_k*np.pi*(y-self.LL())/(self.RR()-self.LL()))
        return psi
    #########————————————————#########
    #看涨和看跌期权Uk，公式12-13，由公式14-15计算得来
    def Uk_call(self,iter_): #（随k而动）
        uk= 2/(self.RR()-self.LL())*self.K*((self.chi_k(self.RR(),iter_)-self.chi_k(0,iter_))-(self.psi_k(self.RR(),iter_)-self.psi_k(0,iter_)))
        return uk
    def Uk_put(self,iter_): #（随k而动）
        uk= 2/(self.RR()-self.LL())*self.K*((-self.chi_k(self.LL(),iter_)-self.chi_k(0,iter_))+(self.psi_k(self.LL(),iter_)-self.psi_k(0,iter_)))
        return uk
    #########————————————————########
    # 期权价格； 公式10：Uk，特征函数，指数（#k=0时在累加时特殊处理）
    def call_option_price(self):
        #当k不为0，求和项数组
        ex= np.exp(-1j*self.array_k*np.pi* (self.LL()/(self.RR()-self.LL())))
        uk=self.Uk_call(1) #传入0是对k=0时特殊处理，传入1则是对后面的迭代矩阵进行处理
        array_sum= ( ex*self.fun_charact(self.phi) ).real *uk
        #数组求和
        sum1= array_sum.sum()
        #当k=0时，即公式10中展开式的第一项
        iter0=self.fun_charact(0).real*self.Uk_call(0)/2
        sum2=sum1+iter0
        price=np.exp(-self.r*self.tau)*sum2
        return price

    #看涨期权价格－看跌期权价值＝标的资产价格－执行价格现值

    def put_option_price(self):
        #当k不为0，求和项数组
        ex= np.exp(-1j*self.array_k*np.pi* (self.LL()/(self.RR()-self.LL())))
        uk=self.Uk_put(1)
        array_sum= ( ex*self.fun_charact(self.phi) ).real *uk
        #数组求和
        sum1= array_sum.sum()
        #当k=0时，即公式10中展开式的第一项
        iter0=self.fun_charact(0).real*self.Uk_put(0)/2
    #     print(iter0)
        sum2=sum1+iter0
        price=np.exp(-self.r*self.tau)*sum2
        return price
    #########————————————————#########
    #希腊字母值
    def greek_delta(self):
        ex= np.exp(-1j*self.array_k*np.pi* (self.LL()/(self.RR()-self.LL())))
        #当k不等于0时
        sum1= (( 1j*self.phi*ex*self.fun_charact(self.phi) ).real *self.Uk_put(1)).sum()
        #当k=0时
        iter0=0
        sum2=sum1+iter0
        greek=np.exp(-self.r*self.tau)/self.S_0*sum2
        return greek
    def greek_gamma(self):
        ex= np.exp(-1j*self.array_k*np.pi* (self.LL()/(self.RR()-self.LL())))
        #当k不等于0时
        sum1= (( 1j*self.phi*ex*self.fun_charact(self.phi) ).real *self.Uk_put(1)).sum()
        #当k=0时
        iter0=0
        sum2=sum1+iter0
        sum3=(( -self.phi**2*ex*self.fun_charact(self.phi) ).real *self.Uk_put(1)).sum()
        greek=np.exp(-self.r*self.tau)/self.S_0**2 *(-sum2 + sum3)
        return greek
    def greek_pseudo_vega(self):
        ex= np.exp(-1j*self.array_k*np.pi* (self.LL()/(self.RR()-self.LL())))
        #当k不等于0时
        sum1= (( B_tau*ex*self.fun_charact(self.phi) ).real *self.Uk_put(1)).sum()
        #当k=0时
        iter0=0
        sum2=sum1+iter0
        greek=np.exp(-self.r*self.tau)*sum2
        return greek

    def greek_theta(self):
        '''   #使用原始导数定义，价格对到期时长求导
            global self.tau
            c0= self.put_option_price()
            self.tau=self.tau+self.delta_d
            c1= self.put_option_price()
            self.tau=self.tau-self.delta_d 
            greek= (c0-c1)/self.delta_d #一阶导'''
        #特征函数对self.tau求导：
        self.fun_charact(self.phi)
        f_self.tau= (A_de_t+B_de_t*self.v_0)*self.fun_charact(self.phi)

        ex= np.exp(-1j*self.array_k*np.pi* (self.LL()/(self.RR()-self.LL())))
        #当k不等于0时
        sum1= ((ex*self.fun_charact(self.phi) ).real *self.Uk_put(1)).sum()
        iter0=self.fun_charact(0).real*self.Uk_put(0)/2
        sum1=sum1+iter0
        sum2= ((ex*f_self.tau).real*self.Uk_put(1)).sum()
        #k为0的项都为0
        greek=np.exp(-self.r*self.tau)*(-r*sum1+sum2)
        return greek

    def greek_rho(self):
        ex= np.exp(-1j*self.array_k*np.pi* (self.LL()/(self.RR()-self.LL())))
        #当k不等于0时
        sum1= (( ex*self.fun_charact(self.phi) ).real *self.Uk_put(1)).sum()
        iter0=self.fun_charact(0).real*self.Uk_put(0)/2
        sum1=sum1+iter0
        sum2=(( 1j*self.tau*self.phi*ex*self.fun_charact(self.phi) ).real *self.Uk_put(1)).sum()
        greek=np.exp(-self.r*self.tau)/self.S_0**2 *(-self.tau*sum1 + sum2)
        return greek

    #########————————————————#########

    #########————————————————#########

    #########————————————————#########

    #########————————————————#########

