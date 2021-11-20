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

print('Starting...')
import time
start_total = time.time()
import numpy as np
import pandas as pd
import sqlite3
import sys


#Get Folder List
arg1 = sys.argv[1]
file_url = pd.read_csv(arg1,header=None,dtype={0:str,1:str})
file_url = [i+'_'+j for i,j in zip(list(file_url[0]),list(file_url[1]))]
print('Read Batch List')

#Get Input url
arg2 = sys.argv[2]

#Get DataBase url
arg3 = sys.argv[3]

#Get Output url
arg4 = sys.argv[4]

#Set PSM q-value threshold
arg5 = sys.argv[5]

#Check Target-Decoy and change to T/F:
def TDcheck(protein_ls):
    protein_ls = protein_ls.split(';')
    if any('XXX' in s for s in protein_ls):
        if any('XXX' not in s for s in protein_ls):
            return 0.5
        else:
            return 1
    else:
        return 0
# ------------------ ESTABLISH DATABASE ----------------------------------------------------------   
#Create a Database
print('Creating Database.sqlite')
Database = sqlite3.connect('%s'%arg3)
print('Database.sqlite was Created')

#Read txt File & Save to Database
for i in file_url:
    print('reading %s'%i)
    start = time.time()
    file = pd.read_csv('%s/%s/%s-peptable.txt'%(arg2,i,i),sep='\t',header=0,
                     usecols=['Peptide sequence', 'Protein(s)','percolator svm-score','PSM q-value'],
                     converters={'Protein(s)':TDcheck})
    file = file.loc[file['PSM q-value'] < float(arg5)].loc[:,['Peptide sequence', 'Protein(s)','percolator svm-score']]
    end = time.time()
    print('%s was read. Pulling to Database.sqlite.\nused %.2f sec.'%(i,end-start))
    start = time.time()
    pd.io.sql.to_sql(file,'data',con=Database,if_exists='append',index=False)
    end = time.time()
    print('%s was pulled to Database.sqlite.\nused %.2f sec.'%(i,end-start))
    
# Delete duplicated rows
print('Removing duplicated rows')
start = time.time()
sql_command = 'DELETE FROM data WHERE ROWID NOT IN (SELECT MIN(ROWID) FROM data GROUP BY `Peptide sequence`)'
Database.execute(sql_command)
end = time.time()
print('Duplicated Rows was removed.\nused %.2f sec.'%(end-start))

# ------------------ CALCULATION OF P/Q VALUE ----------------------------------------------------------   
#Database Length
print('Counting Database Length...')
start = time.time()
length = Database.execute("select COUNT(`Peptide sequence`) FROM data WHERE (`Peptide sequence`,`percolator svm-score`) IN (select `Peptide sequence`,MAX(`percolator svm-score`) from data GROUP BY `Peptide sequence`)")
length = list(length)[0][0]
end =time.time()
print('Database Length=%s.used %.2f sec.'%(length,end-start))

#Create a cursor
print('Creating Cursor...')
start = time.time()
c = Database.cursor()
sql_command = 'SELECT `Peptide sequence`,`Protein(s)`,`percolator svm-score` FROM data WHERE (`Peptide sequence`,`percolator svm-score`) IN (select `Peptide sequence`,MAX(`percolator svm-score`) FROM data GROUP BY `Peptide sequence`) ORDER BY `percolator svm-score` DESC'
c.execute(sql_command)
end =time.time()
print('Cursor was created.used %.2f sec.'%(end-start))

#max(SVM-score)+1
svm_score_change = Database.execute('SELECT MAX(`percolator svm-score`) FROM data')
svm_score_change = list(svm_score_change)[0][0]

#Parameters for loop
Peptide = []
Peptide2 = []
d_list = []
t_list = []
k_list = []
d = 0
t = 0
d_dict = {1.0:1,0.5:1,0.0:0}
k_dict = {1.0:-0.5,0.5:-0.5,0.0:0.5}
t_dict = {1.0:0,0.5:1,0.0:1}

#Create txt File
txt = open('%s'%arg4,'w')
txt.write('modpepseq:ID(Peptide),seq,score:float,qValue:float\n')

#Calculate D
print('Calculating D...')
start = time.time()
D = Database.execute("select `Protein(s)`,COUNT(`Protein(s)`) FROM data WHERE (`Peptide sequence`,`percolator svm-score`) IN (select `Peptide sequence`,MAX(`percolator svm-score`) from data GROUP BY `Peptide sequence`) GROUP BY `Protein(s)` ORDER BY `Protein(s)` DESC")
D = list(D)
if len(D) == 3:
    D = D[0][1]+D[1][1]
elif len(D) == 2:
    D = D[0][1]
end =time.time()
print('D=%s.used %.2f sec.'%(D,end-start))

#Start Calculation p/q-values
print('Calculating p/q-values...')
start = time.time()
for i in range(length):
    peptide,protein,score = c.fetchone()
    if score == svm_score_change:
        Peptide.append(peptide)
        Peptide2.append(''.join([j for j in peptide if j.isalpha()]))
        d_list.append(d_dict[protein])
        t_list.append(t_dict[protein])
        k_list.append(k_dict[protein])
    else:
        #Calculation
        d += sum(d_list)
        t += sum(t_list)
        rep_number = len(Peptide)
        LP = -np.log10((np.repeat(d,rep_number) + np.array(k_list)) / D)
        q_value = np.repeat(d/t,rep_number)
        
        #Wrinting
        w = ''.join([m+','+m2+','+str(n)+','+str(o)+'\n' for m,m2,n,o in zip(Peptide,Peptide2,LP,q_value)])
        txt.write(w)
        
        #Delete list
        Peptide = [peptide]
        Peptide2 = [''.join([j for j in peptide if j.isalpha()])]
        d_list = [d_dict[protein]]
        t_list = [t_dict[protein]]
        k_list = [k_dict[protein]]
        svm_score_change = score
        
    #Last Calculation
    if i == length-1:
        #Calculation
        d += sum(d_list)
        t += sum(t_list)
        rep_number = len(Peptide)
        LP = -np.log10((np.repeat(d,rep_number) + np.array(k_list)) / D)
        q_value = np.repeat(d/t,rep_number)
        
        #Wrinting
        w = ''.join([m+','+m2+','+str(n)+','+str(o)+'\n' for m,m2,n,o in zip(Peptide,Peptide2,LP,q_value)])
        txt.write(w)
        txt.close()
end = time.time()
print('used %.2f sec.'%(end-start))
end_total = time.time()
print('Finished.Used %.2f in all'%(end_total-start_total))
# ------------------- Close Database --------------------------
Database.close()
