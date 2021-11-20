'''
Copyright <2020> <Chien-Hua Chen>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''



from scipy.stats import gamma
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import comb
import os
import time
import scipy.stats as ss
import sys

#argv
#input file(tsv)
arg1 = sys.argv[1]
#output file(Folder url)
arg2 = sys.argv[2]

#Set Figuresize & dpi
plt.rcParams['figure.figsize'] = [10,10]
plt.rcParams['figure.dpi'] = 100

#filename
filename = '[condition 2] The protein-peptide table for Result from [neo4j] Start neo4j database'

#Read file
start = time.time()
print('Calculating target & decoy LPGF')
f = open('%s'%arg1,'r') 
f.readline() #Columns been not needed!


#Create list
p_accession_ls = []
p_type_ls = []
target_LPGF_ls = []
decoy_LPGF_ls = []

#第一個
f1 = f.readline()
f1 = f1.strip().split('\t')
protein = f1[0]
ptype = f1[1]
if (f1[3] == 'target')&(f1[5] == 'true'):
    score_T_True = [float(f1[4])]
    score_T_False = []
    score_D_True = []
    score_D_False = []       
elif (f1[3] == 'decoy')&(f1[5] == 'true'):
    score_T_True = []
    score_T_False = []    
    score_D_True = [float(f1[4])]
    score_D_False = []    
elif (f1[3] == 'target')&(f1[5] == 'false'):
    score_T_True = []
    score_T_False = [float(f1[4])]
    score_D_True = []
    score_D_False = []   
elif (f1[3] == 'decoy')&(f1[5] == 'false'):
    score_T_True = []
    score_T_False = []     
    score_D_True = []
    score_D_False = [float(f1[4])]    

#第二個以後
while True:
    f1 = f.readline()
    if f1 == '':
        #算最後
        #計算Target
        LPF = sum(score_T_True)
        m = len(score_T_True)
        n = m + len(score_T_False)
        if m > 0:
            LPGF_target = -np.log10((1-gamma.cdf(LPF*np.log(10), a=m, scale=1))*comb(n,m))
        else:
            if n != 0:
                LPM = max(score_T_True+score_T_False)
                LPGF_target = -np.log10(1-(1-10**(-LPM))**n)
            else:
                LPGF_target = 'NA'
        #計算Decoy
        LPF = sum(score_D_True)
        m = len(score_D_True)
        n = m + len(score_D_False)
        if m > 0:
            LPGF_decoy = -np.log10((1-gamma.cdf(LPF*np.log(10), a=m, scale=1))*comb(n,m))
        else:
            if n != 0:
                LPM = max(score_D_True+score_D_False)
                LPGF_decoy = -np.log10(1-(1-10**(-LPM))**n)
            else:
                LPGF_decoy = 'NA'
        p_accession_ls.append(protein)
        p_type_ls.append(ptype)
        target_LPGF_ls.append(LPGF_target)
        decoy_LPGF_ls.append(LPGF_decoy)
        #關閉
        f.close()
        break
        
        
        
    f1 = f1.strip().split('\t')

    
    #如果蛋白一樣就儲存
    if f1[0] == protein:
        if (f1[3] == 'target')&(f1[5] == 'true'):
            score_T_True.append(float(f1[4]))
        elif (f1[3] == 'decoy')&(f1[5] == 'true'):
            score_D_True.append(float(f1[4]))
        elif (f1[3] == 'target')&(f1[5] == 'false'):
            score_T_False.append(float(f1[4]))
        elif (f1[3] == 'decoy')&(f1[5] == 'false'):
            score_D_False.append(float(f1[4]))    
    #如果蛋白不一樣就計算 + 寫 + 儲存新的
    else:
        #計算Target
        LPF = sum(score_T_True)
        m = len(score_T_True)
        n = m + len(score_T_False)
        if m > 0:
            LPGF_target = -np.log10((1-gamma.cdf(LPF*np.log(10), a=m, scale=1))*comb(n,m))
        else:
            if n != 0:
                LPM = max(score_T_True+score_T_False)
                LPGF_target = -np.log10(1-(1-10**(-LPM))**n)
            else:
                LPGF_target = 'NA'
        #計算Decoy
        LPF = sum(score_D_True)
        m = len(score_D_True)
        n = m + len(score_D_False)
        if m > 0:
            LPGF_decoy = -np.log10((1-gamma.cdf(LPF*np.log(10), a=m, scale=1))*comb(n,m))
        else:
            if n != 0:
                LPM = max(score_D_True+score_D_False)
                LPGF_decoy = -np.log10(1-(1-10**(-LPM))**n)
            else:
                LPGF_decoy = 'NA'
        p_accession_ls.append(protein)
        p_type_ls.append(ptype)
        target_LPGF_ls.append(LPGF_target)
        decoy_LPGF_ls.append(LPGF_decoy)        
        
        #儲存新的
        protein = f1[0]
        ptype = f1[1]
        if (f1[3] == 'target')&(f1[5] == 'true'):
            score_T_True = [float(f1[4])]
            score_T_False = []
            score_D_True = []
            score_D_False = []       
        elif (f1[3] == 'decoy')&(f1[5] == 'true'):
            score_T_True = []
            score_T_False = []    
            score_D_True = [float(f1[4])]
            score_D_False = []    
        elif (f1[3] == 'target')&(f1[5] == 'false'):
            score_T_True = []
            score_T_False = [float(f1[4])]
            score_D_True = []
            score_D_False = []   
        elif (f1[3] == 'decoy')&(f1[5] == 'false'):
            score_T_True = []
            score_T_False = []     
            score_D_True = []
            score_D_False = [float(f1[4])]    

end = time.time()
print('used %.2f sec.'%(end-start))


#Inf & NA 轉換
Target_score = target_LPGF_ls.copy()
Decoy_score = decoy_LPGF_ls.copy()
#轉換Target
Target_score_rmInf_rmNA = []
max_score = max([i for i in target_LPGF_ls if (str(i)!='inf')&(i!='NA')])
min_score = min([i for i in target_LPGF_ls if (str(i)!='inf')&(i!='NA')])
for i in Target_score:
    if str(i) == 'inf':
        Target_score_rmInf_rmNA.append(max_score+1)
    elif str(i) == 'NA':
        Target_score_rmInf_rmNA.append(min_score-1)
    else:
        Target_score_rmInf_rmNA.append(i)   
        
#轉換Decoy
Decoy_score_rmInf_rmNA = []
max_score = max([i for i in decoy_LPGF_ls if (str(i)!='inf')&(i!='NA')])
min_score = min([i for i in decoy_LPGF_ls if (str(i)!='inf')&(i!='NA')])
for i in Decoy_score:
    if str(i) == 'inf':
        Decoy_score_rmInf_rmNA.append(max_score+1)
    elif str(i) == 'NA':
        Decoy_score_rmInf_rmNA.append(min_score-1)
    else:
        Decoy_score_rmInf_rmNA.append(i)             

        
#Plot target-LPGF vs decoy-LPGF
plt.scatter(Target_score_rmInf_rmNA,Decoy_score_rmInf_rmNA,marker="v",s=3)
lim_upper = max(Target_score_rmInf_rmNA + Decoy_score_rmInf_rmNA)+1
lim_lower = min(Target_score_rmInf_rmNA + Decoy_score_rmInf_rmNA)-1
plt.plot([lim_lower,lim_upper],[lim_lower,lim_upper],color='red')

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.xlim(lim_lower,lim_upper)
#plt.ylim(lim_lower,lim_upper)
plt.xlabel('Target Score',fontsize=16)
plt.ylabel('Decoy Score',fontsize=16)
plt.title('LPGF',fontsize=20)
plt.savefig('%s/tdLPGF_scatter.png'%arg2)            
plt.clf()



#Find FDRr Threshold
start = time.time()
print('Finding Threshold of FDRr...')

Target_score_sorted = np.array(sorted(Target_score_rmInf_rmNA))
Target_score = np.array(Target_score_rmInf_rmNA) 
Decoy_score = np.array(Decoy_score_rmInf_rmNA) 

FDRr = 0
FDRr_max = 0
for ts in Target_score_sorted:
    do = sum((Target_score < ts) & (Decoy_score >= ts))
    to = sum((Target_score >= ts) & (Decoy_score < ts))
    db = sum((Target_score >= ts) & (Decoy_score >= ts) & (Decoy_score > Target_score))
    tb = sum((Target_score >= ts) & (Decoy_score >= ts) & (Decoy_score < Target_score))
    
    to_tb_db = to+tb+db
    if to_tb_db != 0:
        FDRr = (do+2*db) / to_tb_db
        if FDRr < 0.01:
            if FDRr >= FDRr_max:
                FDRr_max = FDRr
                Threshold = ts
                #print('FDRr=%s, Threshold=%s'%(FDRr,ts))  
            if FDRr == 0:
                break
        else:
            #print('FDRr=%s'%FDRr)
            #print('Threshold=%s'%ts)
            break

end = time.time()
print('used %.2f sec.'%(end-start))


#Find qcPass & qcType
qcPass = (Target_score >= Threshold) & (Target_score > Decoy_score)
do = (Target_score < Threshold) & (Decoy_score >= Threshold)
to = (Target_score >= Threshold) & (Decoy_score < Threshold)
db = (Target_score >= Threshold) & (Decoy_score >= Threshold) & (Decoy_score > Target_score)
tb = (Target_score >= Threshold) & (Decoy_score >= Threshold) & (Decoy_score < Target_score)
Not_gate = {True:False,False:True}
nt = list(map(Not_gate.get, do|to|db|tb))

do_dict = {True:'do',False:''}
to_dict = {True:'to',False:''}
db_dict = {True:'db',False:''}
tb_dict = {True:'tb',False:''}
nt_dict = {True:'nt',False:''}
do_ls = list(map(do_dict.get, do))
to_ls = list(map(to_dict.get, to))
db_ls = list(map(db_dict.get, db))
tb_ls = list(map(tb_dict.get, tb))
nt_ls = list(map(nt_dict.get, nt))
qcPass_ls = list(map(str.lower,list(map(str,qcPass))))
qcType_ls = list(map(''.join,zip(do_ls,to_ls,db_ls,tb_ls,nt_ls)))


df_output = pd.DataFrame(list(zip(p_accession_ls, p_type_ls, target_LPGF_ls, decoy_LPGF_ls, qcPass_ls, qcType_ls)))
df_output.columns = ['p.accession', 'p.type', 'target.LPGF', 'decoy.LPGF', 'qcPass', 'qcType']
df_output.to_csv('%s/qcPass.csv'%arg2,na_rep='NA',index=False)


#Plot Distribution of Decoy
start = time.time()
print('Ploting the distribution of Decoy...')
Decoy_score_rmNA = decoy_LPGF_ls.copy()
Decoy_score_rmNA = np.array([i for i in Decoy_score_rmNA if i != 'NA'])
N = len(Decoy_score_rmNA)
decoy_rank = -np.log10(ss.rankdata(-Decoy_score_rmNA)/N)
    
plt.scatter(decoy_rank,Decoy_score_rmNA,marker="$◯$")
lim_upper = max([max(decoy_rank),max(Decoy_score)])
lim_lower = min([min(decoy_rank),min(Decoy_score)])
plt.plot([lim_lower,lim_upper],[lim_lower,lim_upper],color='red')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.xlim(lim_lower,lim_upper)
#plt.ylim(lim_lower,lim_upper)
plt.xlabel('-log$_{10}$(rank/N)',fontsize=16)
plt.ylabel('Protein Probability (LPGF)',fontsize=16)
plt.title('Distribution of Decoy Protein Scores',fontsize=20)
plt.xticks
plt.savefig('%s/dLPGF_rank.png'%arg2)          
plt.clf()   

end = time.time()
print('used %.2f sec.'%(end-start))
print('Threshold of LPGF=%s'%Threshold)
print('FDRr=%s'%FDRr_max)
print('QC passed peptides = %.2f%%(%s/%s)'%(sum(qcPass)/len(qcPass),sum(qcPass),len(qcPass)))