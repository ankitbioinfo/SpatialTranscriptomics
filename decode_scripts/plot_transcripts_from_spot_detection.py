

import matplotlib.pyplot as plt
import cv2
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import numpy as np
import cv2 as cv
import pandas as pd
import sys
import os

def makedir(maindir):
    answer=os.path.isdir(maindir)
    if answer==True:
	    pass
    else:
	    os.mkdir(maindir)


def analysis(name,fname):
    data=pd.read_csv(name)
    maindata=data.to_numpy()

    data=maindata[:,[11,12,13,14]].astype(int)

    dim0=np.max(data[:,0])
    dim1=np.max(data[:,1])
    dim2=np.max(data[:,2])
    dim3=np.max(data[:,3])

    #print(dim3,dim1)
    #print("image dimension from tif files  1053 x 795")

    #name='ppCh2_R2'

    dim3=2069 #Width
    dim1=2075 # Height
    dim3=2304
    dim1=2304


    img = np.ones((dim1,dim3,3), np.uint8)

    print('total',len(data))
    for i in range(len(data)):
        ya=[int(data[i,0]),int(data[i,1])]
        xa=[int(data[i,2]),int(data[i,3])]
        # Draw a diagonal blue line with thickness of 5 px cv2.rectangle(img, (x1, y1), (x2, y2), (255,0,0), 2)
        cv.rectangle(img,(xa[0],ya[0]),(xa[1],ya[1]),(0,255,255),1) #blue

    cv2.imwrite(fname+".png",img)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    cv2.imwrite(fname+".tif",img)


def main():
    rounds=[0,1,2,3,4,5]
    channel=[0,1,2]
    datapath='smfish/'
    #datapath=sys.argv[1]
    maindir=datapath+'png_'
    makedir(maindir)
    for r in rounds:
        for c in channel:
            name='IntensityTable_R'+str(r)+'C'+str(c)+'.dat'
            shortname=maindir+'C'+str(c)+'R'+str(r+1)
            analysis(datapath+name,shortname)

def main2():
    rounds=[0,1,2,3,4,5,6]
    channel=[0,1,2]
    #datapath='case5/'
    #datapath=sys.argv[1]

    for r in rounds:
        for c in channel:
            maindir='png_'+'C'+str(c)+'R'+str(r+1)
            #makedir(maindir)
            for i in range(9,10):
                datapath='case'+str(i)
                makedir(datapath)
                name='/IntensityTable_R'+str(r)+'C'+str(c)+'.dat'
                #shortname=maindir+'/'+datapath
                shortname=datapath+'/'+maindir
                analysis(datapath+name,shortname)

main()
#main2()
