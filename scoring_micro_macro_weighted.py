

import numpy as np


f=open('scoring.dat')

ypred=[]
yref=[]
for line in f:
    l=line.split()
    ypred.append(l[1])
    yref.append(l[2])

#print(ypred)
#print(yref)

def accuracy(a,b):
    N=len(a)
    count=0
    for i in range(N):
        if a[i]==b[i]:
            count=count+1

    return count/float(N)

def weighted_accuracy(a,b):
    N=len(a)
    label=np.unique(a)
    weight=np.zeros(len(label),dtype=float)
    for j in range(len(label)):
        for i in range(N):
            if a[i]==label[j]:
                weight[j]+=1
    for j in range(len(label)):
        weight[j]=weight[j]/N

    total=[]
    for j in range(len(label)):
        count=0.0
        classwisecount=0
        for i in range(N):
            if a[i]==label[j]:
                classwisecount+=1
                if a[i]==b[i]:
                    count=count+1
        #print(weight[j],count/classwisecount)
        total.append(count/float(classwisecount))

    return [np.mean(total), total]

'''
          Reference
Prediction  A  B  C  D  E
         A 35  0  0  0  2
         B  0  9  5  0  2
         C  0  0 10  2  0
         D  5  1  0 23  0
         E  5  0  0  0  1
'''

#none
mat=[[ 44 , 29 , 33 , 29,  27 , 28 , 26,  33,   2],
 [ 52 ,118,  55,  67,  75,  89,  66,  41,   2],
 [ 21,  28,  21,  15,  15,  25,  21,  14,   0],
 [144, 177, 140, 185, 156, 174, 114, 132,   9],
 [ 21,  38,  24,  24,  30,  33,  23,  15,   3],
 [ 17,  67,  50,  31,  31,  53,  23,  26,   0],
 [193, 292, 184, 213, 223, 250, 199, 164,   9],
 [719, 864, 668, 808, 767, 860, 644, 621,  20],
 [  0,   0,   2,   1,   0,   2,   0,   0,   4]]
'''

#true
mat=[[0.1752988 , 0.11553785 ,0.1314741 , 0.11553785 ,0.10756972 ,0.11155378,
  0.10358566 ,0.1314741 , 0.00796813],
 [0.0920354 , 0.20884956 ,0.09734513, 0.11858407, 0.13274336, 0.15752212,
  0.11681416, 0.07256637, 0.00353982],
 [0.13125 ,   0.175,      0.13125,    0.09375,    0.09375,    0.15625,
  0.13125 ,   0.0875 ,    0.        ], [0.11697807 ,0.14378554 ,0.11372868 ,0.15028432 ,0.12672624, 0.1413485,
  0.09260764, 0.10722989, 0.00731113],
 [0.09952607, 0.18009479, 0.11374408, 0.11374408, 0.14218009, 0.1563981,
  0.10900474, 0.07109005, 0.01421801],
 [0.05704698, 0.22483221, 0.16778523, 0.10402685, 0.10402685, 0.17785235,
  0.07718121, 0.08724832, 0.        ],
 [0.11175449, 0.16907933, 0.10654314, 0.12333526, 0.12912565, 0.1447597,
  0.11522872, 0.09496236, 0.00521135],
 [0.12041534, 0.14469938, 0.11187406, 0.13532072, 0.1284542,  0.14402948,
  0.10785463, 0.10400268, 0.00334952], [0.  ,       0.    ,     0.22222222 ,0.11111111 ,0.   ,      0.22222222,
  0.  ,       0.   ,      0.44444444]]

#pred
mat=[[0.03633361 ,0.01797892, 0.02803738, 0.02112163, 0.02039275, 0.01849406,
  0.02329749, 0.03154876, 0.04081633],
 [0.04293972 ,0.07315561, 0.04672897 ,0.04879825, 0.05664653 ,0.05878468,
  0.05913978, 0.03919694, 0.04081633],
 [0.01734104, 0.01735896, 0.01784197, 0.01092498, 0.01132931, 0.01651255,
  0.0188172,  0.01338432, 0.        ],
 [0.11890999, 0.10973342, 0.11894647, 0.13474144, 0.11782477, 0.11492734,
  0.10215054, 0.12619503, 0.18367347],
 [0.01734104 ,0.02355859, 0.02039082, 0.01747997, 0.02265861, 0.02179657,
  0.02060932, 0.01434034, 0.06122449],
 [0.01403799, 0.04153751, 0.04248088, 0.0225783,  0.0234139,  0.03500661,
  0.02060932, 0.0248566,  0.        ]
, [0.15937242, 0.18102914, 0.15632965 ,0.15513474 ,0.168429  , 0.1651255,
  0.17831541, 0.15678776 ,0.18367347],
 [0.59372419, 0.53564786, 0.5675446 , 0.58849235, 0.57930514, 0.5680317,
  0.57706093, 0.59369025 ,0.40816327],
 [0.  ,       0.  ,       0.00169924 ,0.00072833 ,0. ,        0.001321,
  0.,         0.  ,       0.08163265]]
'''


def macro_avg(a,b):
    N=len(a)
    label=np.unique(a)
    TP=[35,9,10,23,1]
    FP=[2,7,2,6,5]
    FN=[10,1,5,2,4]

    TP=[]
    FP=[]
    FN=[]
    for i in range(len(mat)):
        #print(mat[i][1])
        TP.append(mat[i][i])
        rowsum=0
        colsum=0
        for j in range(len(mat)):
            if j!=i:
                colsum+=mat[i][j]
                rowsum+=mat[j][i]

        FN.append(colsum)
        FP.append(rowsum)

    #print(TP,'\n',FP,'\n',FN)


    Pre=[]
    Rec=[]
    F1=[]
    for k in range(len(TP)):
        v1=TP[k]/float(TP[k]+FP[k])
        v2=TP[k]/float(TP[k]+FN[k])
        Pre.append(v1)
        Rec.append(v2)
        F1.append(2*v1*v2/(v1+v2))

    print("Pre",Pre)
    print("Rec",Rec)

    Pre=sum(Pre)/len(TP)
    Rec=sum(Rec)/len(TP)
    F1= sum(F1)/len(TP)

    return (Pre,Rec,F1)

def micro_avg(a,b):
    N=len(a)
    label=np.unique(a)
    TP=[35,9,10,23,1]
    FP=[2,7,2,6,5]
    FN=[0,0,0,0,0]
    #FN=[10,1,5,2,4]
    TP=[]
    FP=[]
    FN=[]
    for i in range(len(mat)):
        #print(mat[i][1])
        TP.append(mat[i][i])
        rowsum=0
        colsum=0
        for j in range(len(mat)):
            if j!=i:
                colsum+=mat[i][j]
                rowsum+=mat[j][i]

        FN.append(colsum)
        FP.append(rowsum)


    tp=sum(TP)
    fp=sum(FP)
    fn=sum(FN)

    Pre= tp / (tp+fp)
    Rec= tp / (tp+fn)
    F1= 2*Pre*Rec/(Pre+Rec)

    return (Pre,Rec,F1)


def weight_avg(a,b):
    N=len(a)
    weight=[0.024035117177193026, 0.05423493007747847, 0.015327784307611121, 0.11811269159730385, 0.020293122886133032, 0.028544700040778144, 0.1657271702367531, 0.5728609465326585, 0.0008635371440907674]
    TP=[]
    FP=[]
    FN=[]
    for i in range(len(mat)):
        #print(mat[i][1])
        TP.append(mat[i][i])
        rowsum=0
        colsum=0
        for j in range(len(mat)):
            if j!=i:
                colsum+=mat[i][j]
                rowsum+=mat[j][i]

        FN.append(colsum)
        FP.append(rowsum)

    Pre=[]
    Rec=[]
    F1=[]
    for k in range(len(TP)):
        v1=TP[k]/float(TP[k]+FP[k])
        v2=TP[k]/float(TP[k]+FN[k])
        Pre.append(v1)
        Rec.append(v2)
        F1.append(2*v1*v2/(v1+v2))

    for i in range(len(TP)):
        Pre[i]=Pre[i]*weight[i]
        Rec[i]=Rec[i]*weight[i]

        F1[i]=F1[i]*weight[i]

    Pre=sum(Pre)/sum(weight)
    Rec=sum(Rec)/sum(weight)
    F1= sum(F1)/sum(weight)

    return (Pre,Rec,F1)



print("accuracy",accuracy(yref,ypred))
print("weighted accuracy",weighted_accuracy(yref,ypred))
print("macro (P, R, F1) avg",macro_avg(yref,ypred))

print("micro avg",micro_avg(yref,ypred))

print("weight (P, R, F1) avg",weight_avg(yref,ypred))
