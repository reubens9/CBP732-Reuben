
# coding: utf-8

# In[1]:


from scipy.integrate import odeint
from scipy.optimize import fsolve
import numpy as np
from matplotlib import pyplot as plt
import pandas
import operator
get_ipython().magic('matplotlib inline')


# # $Succi$

# ##  Mass balance
# 
# $$Xylose + CO_2 + NH_3 = X + Succinic~ acid + Acetic ~acid  + Formic~ acid  + Water$$
# 
# 
# $$CH_2O + CO_2 + NH_3= CH_{1.896}O_{0.603}N_{0.258} + CH_{1.5}O + CH_2O + CH_2O_2 + H_2O$$
# 
# ## Metabolic flux model
# <img src="succi.png" width="400" />

# In[2]:


"""Mass balance"""

Mass_S = np.matrix([[1,1,1    ,1,  1,1,0,0],    #C
                    [2,0,1.896,1.5,2,2,2,3],    #H
                    [1,2,0.603,1,  1,2,1,0],    #O
                    [0,0,0.258,0,  0,0,0,1],    #N   
                   [-1,0,0,    0,  0,0,0,0],    #Xylose
                    [0,0,0,    1,  0,0,0,0],    #Succinic acid
                    [0,0,0,    0,  1,0,0,0],    #Acetic acid 
                    [0,0,0,    0,  0,1,0,0],    #Formic acid
                   ])

MM = 12.11 + 1.0078*1.896 + 16*0.603 + 14.006*0.258
MM 


# In[4]:


DataS = pandas.read_excel('Data.xlsx',sheetname = 'SData')
# print (DataS)
DataS = np.array(DataS)
Slen = range(len(DataS))
Results = pandas.ExcelWriter('Results.xlsx')

def MBs(i): 
    "Component order = Xylose, CO2, X, SA, AA, FA, H20"
    B = np.matrix([0,0,0,0,(40 - DataS[i,1])/(150.13/5),DataS[i,2]/(118.09/4),DataS[i,3]/(60.05/2) ,DataS[i,4]/46.03]).T
    X = np.linalg.solve(Mass_S,B)
    
    return [X[0,0]*(150.13/5), X[1,0]*44.01, X[2,0]*MM, X[3,0]*(118.09/4), X[4,0]*(60.05/2), X[5,0]*46.03]

MB_Succi = [MBs(i) for i in Slen]
DS = list(DataS[:,0])
BS = list(DataS[:,5])

MS = pandas.DataFrame(columns= ["Dilution (1/h)",  "Xylose (g/l)", "$CO_2$ (g/L)","Biomass (g/l)",
                                       "Succinic (g/l)",  "Acetic (g/l)",  "Formic (g/l)"])
Biomass_S = pandas.DataFrame(columns= ["Biomass calculated (g/l)","Biomass (g/l) Ladakis, et al. (2018)"])

for i in Slen:
    MS.loc[i] = [DS[i]] + MB_Succi[i]

LS = list(MS['Biomass (g/l)'])
for i in Slen:
    Biomass_S.loc[i] =  [LS[i]] + [BS[i]]
    
Xe = sum(abs((MS['Biomass (g/l)'] - DataS[:,5])/MS['Biomass (g/l)']))/Slen[-1] * 100


Balance =  MS['Xylose (g/l)']/(150.13/5) + MS['$CO_2$ (g/L)']/44.01 + (MS['Biomass (g/l)']/MM +
                 MS['Succinic (g/l)']/(118.09/4) + MS['Acetic (g/l)']/(60.05/2) + MS['Formic (g/l)']/46.03)
                                            
C_balance_errorS = sum(Balance)/len(Balance)
Biomass_errorS = Xe


# In[5]:


"""Flux model"""

S = np.matrix(pandas.read_excel('Data.xlsx',sheetname = 'Succ'))
alpha, beta = 0.105, 0.165

def rates(i):
    Cxy, Cco2, Cx, Csa, Caa, Cfa = MBs(i)
    C = [Cxy/(150.13/5), Cco2/44.01, Cx/MM, Csa/(118.09/4), Caa/(60.05/2), Cfa/46.03]
    CX = C[2]
    D = DataS[i,0]
    r = np.zeros_like(C)
    for j in range(len(C)):
        r[j] = D*C[j]/CX
        
    return r
Xerror = []
def metafluxs(i):
    r = rates(i)
    rco2 = r[0] + r[2] + r[3] + r[4] + r[5]
    D = DataS[i,0]
    specs = np.matrix([0,0,0,0,0,0,-r[0],0,D,r[3],r[4],r[5]]).T
    F =  np.linalg.solve(S,specs)
    CO2 = F[1,0]*alpha + 1/6*F[4,0] - 1/4*F[8,0] + 1/2*F[10,0]
    Xerror.append(abs((r[1] - CO2)/r[1])*100)
#     print (r[1], CO2, rco2)
#     print (-r[0], F[0,0])
#     print (r[2]- F[1,0])
#     print (r[3]- F[8,0])
#     print (r[4]- F[10,0] - F[11,0])
#     print (r[5]- 0.5*F[11,0])

    return [F[0,0],F[1,0],F[2,0],F[3,0],F[4,0],F[5,0],F[6,0],F[7,0],F[8,0],F[9,0],F[10,0],F[11,0]]

FluxS = [metafluxs(i) for i in Slen]
FS =  pandas.DataFrame(columns= ('µ (1/h)','r_xylose (cmol/cmolX/h)','r_CO_2 (cmol/cmolX/h)',
                                'r_SA (cmol/cmolX/h)', 'r_AA (cmol/cmolX/h)', 'r_FA (cmol/cmolX/h)',))


for i in Slen:
    rCO2 = FluxS[i][1]*alpha + 1/6*FluxS[i][4] - 1/4*FluxS[i][8] + 1/2*FluxS[i][10]

    Bal = FluxS[i][0] - FluxS[i][1] - FluxS[i][8] - FluxS[i][10] - FluxS[i][11] - 0.5*FluxS[i][11]
#     print (Bal, rCO2)

    FS.loc[i] = [FluxS[i][1], FluxS[i][0], rCO2, FluxS[i][8], FluxS[i][10] + FluxS[i][11], 0.5*FluxS[i][11]]

CO2_errorS = sum(Xerror)/len(Xerror)

FS


# #  $Basfia$

# ## Mass balance
# 
# $$Xylose + CO_2 + NH_3 = X + Succinic~ acid + Acetic ~acid  + Formic~ acid + Lactic~ acid$$
# 
# 
# $$CH_2O + CO_2 + NH_3= CH_{1.896}O_{0.603}N_{0.258} + CH_{1.5}O + CH_2O + CH_2O_2 + CH_20 $$
# 
# ## Metabolic flux model
# 
# <img src="Basfia.png" width="400" />

# In[6]:


"""Mass balance"""

Mass_B = np.matrix([[1,1,1,1,1,1,1,0,0],  #C
                    [2,0,1.8,1.5,2,2,2,2,3], #H
                    [1,2,0.5,1,1,2,1,1,0],   #O
                    [0,0,0.2,0,0,0,0,0,1],       #N
                   [-1,0,0,0,0,0,0,0,0],      #Xylose
                    [0,0,0,1,0,0,0,0,0],       #Succinic acid
                    [0,0,0,0,1,0,0,0,0],       #Acetic acid 
                    [0,0,0,0,0,1,0,0,0],       #Formic acid
                    [0,0,0,0,0,0,1,0,0],       #Lactic acid
                  ])

MM = 12 + 1.8 + 16*0.5 + 14*0.2
MM


# In[9]:


DataB = pandas.read_excel('Data.xlsx',sheetname = 'BData')
# print (DataB)
DataB = np.array(DataB)
Blen = range(len(DataB))

def MBb(i): 
    "Component order = Xylose, CO2, X, SA, AA, FA, LA, H20"
    B = np.matrix([0,0,0,0,(40 - DataB[i,1])/(150.13/5),DataB[i,2]/(118.09/4),
                   DataB[i,3]/(60.05/2) ,DataB[i,4]/46.03, DataB[i,5]/(90.08/3)]).T
    X = np.linalg.solve(Mass_B,B)

    return [X[0,0]*(150.13/5), X[1,0]*44.01, X[2,0]*MM, X[3,0]*(118.09/4), X[4,0]*(60.05/2), X[5,0]*46.03, X[6,0]*(90.08/3)]

MB_Basfia = [MBb(i) for i in Blen]
DB = list(DataB[:,0])
BB = list(DataB[:,6])

MB = pandas.DataFrame(columns= ["Dilution (1/h)",  "Xylose (g/l)", "$CO_2$ (g/L)","Biomass (g/l)",
                                       "Succinic (g/l)",  "Acetic (g/l)",  "Formic (g/l)", 'Lactic (g/l)'])

Biomass_B = pandas.DataFrame(columns= ["Biomass calculated (g/l)","Biomass (g/l) Ladakis, et al. (2018)"])

for i in Blen:
    MB.loc[i] = [DB[i]] + MB_Basfia[i]

LB = list(MB['Biomass (g/l)'])
for i in Blen:
    Biomass_B.loc[i] =  [LB[i]] + [BB[i]]

Biomass_errorB = sum(abs((MB['Biomass (g/l)'] - DataB[:,6])/MB['Biomass (g/l)']))/Blen[-1] * 100

Balance =  MB['Xylose (g/l)']/(150.13/5) + MB['$CO_2$ (g/L)']/44.01+ (MB['Biomass (g/l)']/MM +
                   MB['Succinic (g/l)']/(118.09/4) + MB['Acetic (g/l)']/(60.05/2) + 
                   MB['Formic (g/l)']/46.03  + MB['Lactic (g/l)']/(90.08/3))
                                            
C_balance_errorB = sum(Balance)/len(Balance)
MB
Biomass_B


# In[10]:


"""Flux model"""
B = np.matrix(pandas.read_excel('Data.xlsx',sheetname = 'Basf'))
alpha, beta = 0.1,0.1
B[12,1] = beta
B[0,1] = 1 + alpha

def rates(i):
    Cxy, Cco2, Cx, Csa, Caa, Cfa, Cla  = MBb(i)
    C = [Cxy/(150.13/5), Cco2/44.01, Cx/MM, Csa/(118.09/4), Caa/(60.05/2), Cfa/46.03, Cla/(90.08/3)]
    CX = C[2]
    D = DataB[i,0]
    r = np.zeros_like(C)
    for j in range(len(C)):
        r[j] = D/CX*C[j]
        
    return r

Xerror = []

def metafluxs(i):
    r = rates(i)
    rco2 = r[0] +r[2] + r[3] + r[4] + r[5] + r[6]
    D = DataB[i,0]
    specs = np.matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,D,r[3],r[4],r[5],r[6],-r[0]]).T
    F =  np.linalg.solve(B,specs)
    CO2 = F[1,0]*alpha + 1/6*F[4,0] - 1/4*F[8,0] + 1/2*F[10,0] + 1/5*F[15,0] + 1/4*F[16,0]
    Xerror.append((abs(r[1] - CO2)/-r[1])*100)
#     print (r[1], CO2, rco2)
#     Xerror.append(abs((r[0] + F[0,0])/r[0])*100)
#     print (-r[0], F[0,0])
#     print (r[2]- F[1,0])
#     print (r[3]- F[18,0] -  F[16,0])
#     print (r[4]- F[14,0])
#     print (r[5]- 0.5*F[11,0])
#     print (r[6]- F[12,0])
    
    return [F[0,0],F[1,0],F[2,0],F[3,0],F[4,0],F[5,0],F[6,0],F[7,0],F[8,0],F[9,0],
            F[10,0],F[11,0],F[12,0],F[13,0],F[14,0],F[15,0],F[16,0],F[17,0],F[18,0]]

FluxB = [metafluxs(i) for i in Blen]
FB =  pandas.DataFrame(columns= ('µ (1/h)','r_xylose (cmol/cmolX/h)','r_CO_2 (cmol/cmolX/h)',
                                'r_SA (cmol/cmolX/h)', 'r_AA (cmol/cmolX/h)', 'r_FA (cmol/cmolX/h)','r_LA (cmol/cmolX/h)'))


for i in Blen:
    rCO2 = FluxB[i][1]*alpha + 1/6*FluxB[i][4] - 1/4*FluxB[i][8] + 1/2*FluxB[i][10] + 1/5*FluxB[i][15] + 1/4*FluxB[i][16]
    Bal = FluxB[i][0] - FluxB[i][1] - FluxB[i][16] - FluxB[i][18] - FluxB[i][14] - 0.5*FluxB[i][11] - FluxB[i][12]
#     print (Bal,rCO2)
    FB.loc[i] = [FluxB[i][1], FluxB[i][0], rCO2, FluxB[i][16] + FluxB[i][18], FluxB[i][14], 0.5*FluxB[i][11], FluxB[i][12]]

CO2_errorB = sum (Xerror)/len(Xerror)
FB
Biomass_B


# # Analysis

# In[11]:


rATPS = []
v6_v4s = []
rATPB = []
v6_v4b = []

for i in Slen:
    rATPS.append(-1/5*FluxS[i][2] + 1/3*FluxS[i][7] + 1/3*FluxS[i][9] + 1/2*FluxS[i][10] + 1/2*FluxS[i][11] + 5/12*FluxS[i][8]) 
    v6_v4s.append(FluxS[i][6]/FluxS[i][4])
for i in Blen:
    rATPB.append(-1/5*FluxB[i][2] + 1/3*FluxB[i][7] + 1/3*FluxB[i][9] + 1/2*FluxB[i][14] + 
                 1/4*FluxB[i][8] + 1/4*FluxB[i][16] + 1/6*FluxB[i][18]) 
    v6_v4b.append(FluxB[i][6]/FluxB[i][4])
    
plt.figure('ATP')
plt.title('ATP production')
plt.plot(DB,rATPB,'ro', label = 'Basfia')
plt.plot(DS,rATPS,'ko',label = 'Succi')
plt.ylabel('$r_{ATP}$ (mol/cmolX/h)')
plt.xlabel('D (1/h)')
plt.legend(loc = 'best')
plt.savefig('ATP.png')


# In[12]:


plt.figure()
plt.title('Succinic acid production')
plt.plot(FB['µ (1/h)'],FB['r_SA (cmol/cmolX/h)'],'rv', label = "Basfia")
plt.plot(FS['µ (1/h)'],FS['r_SA (cmol/cmolX/h)'],'kv',label = "Succi")
plt.ylabel('Rate (cmol/cmolX/h)')
plt.xlabel('D (1/h)')
plt.legend(loc = 'best')
plt.savefig('SUC.png')
plt.figure()
plt.title('Succinic acid concentration')
plt.plot(FB['µ (1/h)'],MB['Succinic (g/l)'],'rv', label = "Basfia")
plt.plot(FS['µ (1/h)'],MS['Succinic (g/l)'],'kv',label = "Succi")
plt.ylabel('Rate (cmol/cmolX/h)')
plt.xlabel('D (1/h)')
plt.legend(loc = 'best')
plt.savefig('SUCconc.png')


# In[13]:


plt.figure()
plt.title('$CO_2$ consumption')
plt.plot(FB['µ (1/h)'],-FB['r_CO_2 (cmol/cmolX/h)'],'r*', label = "Basfia")
plt.plot(FS['µ (1/h)'],-FS['r_CO_2 (cmol/cmolX/h)'], 'k*',label = "Succi")
plt.ylabel('Rate (cmol/cmolX/h)')
plt.xlabel('D (1/h)')
plt.legend(loc = 'best')
plt.savefig('CO2.png')


# In[14]:


Sv8 = [FluxS[i][8] for i in Slen]
Bv8 = [FluxB[i][8] for i in Blen]

plt.figure()
plt.title('PEP to OXA pathway')
plt.plot(FB['µ (1/h)'],Bv8,'rx', label = "Basfia")
plt.plot(FS['µ (1/h)'],Sv8,'kx',label = "Succi")
plt.ylabel('Rate (cmol/cmolX/h)')
plt.xlabel('D (1/h)')
plt.legend(loc = 'best')
plt.savefig('PEPOXA.png')


# In[15]:


Bv18 = [FluxB[i][18] for i in Blen]
Bv16 = [FluxB[i][16] for i in Blen]

plt.figure()
plt.title('Succinic acid production pathways')
plt.plot(FB['µ (1/h)'],Bv18,'rD', label = "Reductive")
plt.plot(FB['µ (1/h)'],Bv16,'kD',label = "Oxidative")
plt.ylabel('Rate (cmol/cmolX/h)')
plt.xlabel('D (1/h)')
plt.legend(loc = 'best')
plt.savefig('REDOX.png')


# In[16]:


Results = pandas.ExcelWriter('Results.xlsx')
MS.to_excel(Results, sheet_name = 'MBA_Succi')
FS.to_excel(Results, sheet_name = 'FA_Succi')
MB.to_excel(Results, sheet_name = 'MBA_Basfia')
FB.to_excel(Results, sheet_name = 'FA_Basfia')


# In[17]:


MS


# In[18]:


FS


# In[19]:


MB


# In[20]:


FB


# In[21]:


FluxB

