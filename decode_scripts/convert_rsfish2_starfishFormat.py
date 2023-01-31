
'''
IntensityTable_allR_allC.dat
IntensityTable_R0C0.dat
IntensityTable_R0C1.dat
IntensityTable_R0C2.dat
IntensityTable_R1C0.dat
IntensityTable_R1C1.dat
IntensityTable_R1C2.dat
IntensityTable_R2C0.dat
IntensityTable_R2C1.dat
IntensityTable_R2C2.dat
IntensityTable_R3C0.dat
IntensityTable_R3C1.dat
IntensityTable_R3C2.dat
'''


fname=['smFISH localizations_R1_z3_G.csv','smFISH localizations_R1_z3_R.csv','smFISH localizations_R1_z3_FR.csv',
'smFISH localizations_R2_z4_G.csv','smFISH localizations_R2_z4_R.csv','smFISH localizations_R2_z4_FR.csv',
'smFISH localizations_R3_z4_G.csv','smFISH localizations_R3_z4_R.csv','smFISH localizations_R3_z4_FR.csv',
'smFISH localizations_R4_z3_G.csv','smFISH localizations_R4_z3_R.csv','smFISH localizations_R4_z3_FR.csv',
'smFISH localizations_R5_z4_G.csv','smFISH localizations_R5_z4_R.csv','smFISH localizations_R5_z4_FR.csv',
'smFISH localizations_R6_z4_G.csv','smFISH localizations_R6_z4_R.csv','smFISH localizations_R6_z4_FR.csv'
]

#x,y,t,c,intensity
header='features,r,c,intensity,z,y,x,radius,spot_id,z_min,z_max,y_min,y_max,x_min,x_max,IntensityTable\n'


count=0
for i in range(6):
    for j in range(3):
        fw=open('./smfish/IntensityTable_R'+str(i)+'C'+str(j)+'.dat','w')
        f=open('./smfishAndy/'+fname[count])
        cont=f.readlines()
        fw.write(header)
        for k in range(1,len(cont)):
            l=cont[k][0:-1].split(',')
            x=l[0]
            y=l[1]
            #print(x,y)
            ymin=float(y)-0.5
            ymax=float(y)+0.5
            xmin=float(x)-0.5
            xmax=float(x)+0.5
            t=l[2]
            c=str(float(l[3]))
            intensity=str(float(l[4]))
            if ((xmin>0) & (ymin>0)):
                fw.write('0,0,'+c+','+intensity+',0,0,0,0,0,0,1,'+str(ymin)+','+str(ymax)+','+str(xmin)+','+str(xmax)+','+intensity+'\n')
        print(len(cont))
        count+=1
