#### 计算物理作业 Ising Model

##### 李子霄 2017301020223

#### problem 1

这个问题分两步。第一步得用蒙特卡洛模拟给出一个二维平面上Ising Model的最终构型，第二步是根据得到的构图通过数据处理的方法得到Tc和$\beta$是多少

代码如下：

```
from __future__ import division
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt


def initialstate(N):   
    ###generate a random spin initial state, with -1 for down and 1 for up
    state = 2*np.random.randint(0,2, size=(N,N))-1
    return state


def mcmove(config, beta):
    

    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N)#chose a random lattice to start the interaction
                b = np.random.randint(0, N)
                s =  config[a, b]
                nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                cost = 2*s*nb
                if cost < 0:
                    s=-s
                elif rand() < np.exp(-cost*beta):
                    s=-s
                config[a, b] = s
    return config


#def calcEnergy(config):

 #   '''Energy of a given configuration'''

  #  energy = 0

   # for i in range(len(config)):

#        for j in range(len(config)):

#            S = config[i,j]

#            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]

#            energy += -nb*S

#    return energy/4.

def calcMag(config):
    #summing up the Magnetization 
    mag = np.sum(config)
    return mag

nt      = 2**5        # 温度点数量
N       = 2**4        # 点阵尺寸, N x N
eqSteps = 2**5       # MC方法平衡步数
mcSteps = 2**9       # MC方法计算步数

n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)
tm = 2.269;    T=np.random.normal(tm, .64, nt)
T  = T[(T>1.2) & (T<3.8)];    nt = np.size(T)

Energy       = np.zeros(nt)
Magnetization  = np.zeros(nt)
Susceptibility = np.zeros(nt)


for m in range(len(T)):
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N)
    iT=1.0/T[m]; iT2=iT*iT;
    

    for i in range(eqSteps):         # equilibrate
        mcmove(config, iT)           # Monte Carlo moves
    
    for i in range(mcSteps):
        mcmove(config, iT)           
        
        Mag = calcMag(config)        # calculate the magnetisation


​        

        M1 = M1 + Mag
        M2 = M2 + Mag*Mag 
        E2 = E2 + Ene*Ene


​        

        Magnetization[m]  = n1*M1
        Susceptibility[m] = (n1*M2 - n2*M1*M1)*iT

f = plt.figure(figsize=(8, 10)); # plot the calculated values    
sp =  f.add_subplot(2, 1, 1 );
plt.plot(T, np.abs(Magnetization), 'o', color='Blue');
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Magnetization ", fontsize=20);
sp =  f.add_subplot(2, 1, 2 );
plt.plot(T, Susceptibility, 'o', color="Green");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Susceptibility", fontsize=20);
```

第一步得到的未经过处理的数据如下：

![image-20201130221536981](C:\Users\李子霄\AppData\Roaming\Typora\typora-user-images\image-20201130221536981.png)

这个图像是吻合相变点附近的行为的，但是不能明确的告诉我们Tc和$\beta$是多少；我们只知道这些数据点满足$M\sim(Tc-T)^\beta$那么这就有两种思路去算了

1.利用$log(M)\sim\beta log(Tc-T)$使得$log(M)\ verse\  log(Tc-T)$是直线来判断哪个Tc正确，再算出相应的$\beta$

代码如下：

```
#### approach 2

log_Magnetization=np.log(np.abs(Magnetization))
f = plt.figure(figsize=(10, 20)); 
for i in np.arange(1,10):
    sp =  f.add_subplot(5, 2 ,i);
    tc=2+0.1*i
    log_tc_t=np.log(tc-T)
    plt.plot(log_tc_t,log_Magnetization,'o')
    plt.xlabel('log(Tc-T)')
    plt.ylabel('log(Magnetization)')
log_tc_t=np.log(np.abs(2.3-T))
log_Magnetization=np.log(np.abs(Magnetization))
z1 = np.polyfit(log_tc_t,log_Magnetization, 1)
# 生成多项式对象
p1 = np.poly1d(z1)
print(z1)
print(p1)
```

![image-20201130214406372](C:\Users\李子霄\AppData\Roaming\Typora\typora-user-images\image-20201130214406372.png)

![image-20201130215249310](C:\Users\李子霄\AppData\Roaming\Typora\typora-user-images\image-20201130215249310.png)

结果表明Tc=2.341，$\beta=0.26$， Tc比较准确，$\beta$差了一倍...

2.利用$M^{\frac{1}{\beta}}\sim Tc-T$使得$M^{\frac{1}{\beta}}$和$Tc$是直线关系来求$Tc$和$\beta$，即先假定一组$\beta$，然后再用这组中$\beta$比较好的数值进一步估计Tc

```
####approach 1 
f = plt.figure(figsize=(10, 30)); #a list of Log(M^(1/beta)) verse Log(T)

for beta_trial in np.arange(1,12):
    sp =  f.add_subplot(8, 2 ,beta_trial);
    beta_trial=1/beta_trial
    M2=Magnetization**beta_trial
    plt.plot(T,M2,'o')
```

![image-20201130223221811](C:\Users\李子霄\AppData\Roaming\Typora\typora-user-images\image-20201130223221811.png)

![image-20201130223304110](C:\Users\李子霄\AppData\Roaming\Typora\typora-user-images\image-20201130223304110.png)

肉眼观察结果表明最接近直线行为的是$\beta=1/9$的情况，Tc=2.3





#### problem 2

时间不够，真的不会做...咋今年出的题都这么难呢，我这半吊子编程水平搞不定啊...

