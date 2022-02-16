import pandas as pd

#reference database Martin Guilliams et al https://doi.org/10.1016/j.cell.2021.12.018

df=pd.ExcelFile('1-s2.0-S0092867421014811-mmc1.xlsx')

#print(df.sheet_names)

sheetnames=['Mouse Endothelial DEGs', 'Mouse Endothelial DEPs', 'Mouse Fibroblast DEGs', 'Mouse Fibroblast DEPs', 'Mouse Hepatocyte DEGs', 'Mouse Cholangiocyte DEGs', 'Mouse Cholangiocyte DEPs', 'Mouse HsPCs DEGs', 'Mouse cDC2 DEGs', 'Mouse cDC2 DEPs', 'Mouse Mig. cDCs DEGs', 'Mouse Mig. cDC DEPs', 'Mouse Monocyte DEGs', 'Mouse Monocyte DEPs', 'Mouse Neutrophil DEGs', 'Mouse Neutrophil DEPs', 'Mouse Basophil DEGs', 'Mouse Basophil DEPs', 'Mouse NK cell DEGs', 'Mouse NK cell DEPs', 'Mouse ILC1 DEGs', 'Mouse ILC1 DEPs', 'Mouse T cell DEGs', 'Mouse T cell DEPs', 'Mouse pDCs DEGs', 'Mouse pDC DEPs', 'Mouse B cell DEGs', 'Mouse B cell DEPs', 'Mouse cDC1 DEGs', 'Mouse cDC1 DEPs', 'Mouse KC DEGs', 'Mouse KC DEPs', 'Human Endothelial DEGs', 'Human Endothelial DEPs', 'Human Fibroblast DEGs', 'Human Hepatocyte DEGs', 'Human Cholangiocyte DEGs', 'Human cDC2 DEGs', 'Human cDC2 DEPs', 'Human Mig. cDC DEGs', 'Human Mig. cDC DEPs', 'Human Monocyte DEGs', 'Human Monocyte DEPs', 'Human Neutrophil DEGs', 'Human Neutrophil DEPs', 'Human Basophil DEGs', 'Human Basophil DEPs', 'Human Res. NK cell DEGs', 'Human Res. NK cell DEPs', 'Human Circ. NK NKT cell DEGs', 'Human Circ. NK NKT cell DEPs', 'Human T cell DEGs', 'Human T cell DEPs', 'Human pDC DEGs', 'Human pDC DEPs', 'Human B cell DEGs', 'Human B cell DEPs', 'Human Plasma cell DEGs', 'Human Plasma cell DEPs', 'Human cDC1 DEGs', 'Human cDC1 DEPs', 'Human Macrophage DEGs', 'Human Macrophage DEPs']


def findDEgenes(sheet,mygene):
    data=df.parse(sheet)
    data=data.to_numpy()
    gname=data[:,0]
    check=[]
    for i in range(len(mygene)):
        found=0
        for j in range(len(data)):
            if mygene[i].lower()==gname[j].lower():
                found=1
        check.append(found)
        #print(mygene[i],found)
    return check


mygene=[['alb', 'mup20', 'hpx', 'hsd17b13', 'cyp2f2'],
['akr1c6', 'rgn', 'sult2a8', 'fabp1'],
['fabp1', 'apoa2', 'sord', 'akr1c6'],
['fabp1', 'apoa2', 'clec4g', 'gpr182', 'plpp3', 'kdr'],
['clec4g', 'plpp3', 'kdr', 'maf', 'ehd3'],
['tyrobp', 'ptpn18', 'ccl5', 'laptm5', 'rbm24', 'uckl1os', 'aurkc','dbx2','igkc','hcst','ly6c2'],

['alb', 'apoa2', 'serpina3k', 'mup20', 'hsd17b13', 'cyp2f2','ccl5'],
['clec4g', 'kdr', 'eng', 'plpp3', 'mrc1', 'maf', 'stab2', 'tyrobp'],
['lyz2', 'ctss', 'alox5ap', 'mpeg1', 'spi1', 'h2-ab1', 'h2-aa', 'unc93b1'],
['dcn', 'sparc', 'ecm1', 'angptl2', 'tcf21', 'bgn', 'col14a1', 'colec11'],
['clec4g', 'kdr', 'plpp3', 'maf'],

['igfbp5', 'igfbp6', 'ccl5', 'upk3b', 'krt19', 's100a6','trav5-4','dcn','igkc','olfr1484'],
['krt19', 'upk3b', 'igfbp5', 'igfbp6', 'trav5-4','dcn','mup20','alb','nrgn','grm5','olfr1277'],
['epcam', 'spp1','vwf', 'fmo2', 'fabp4', 'ccl5', 'iqgap1', 'plvap', 'clu', 'prss23','foxp1','krt19','ddit4l','fbln2','igfbp5'],
['dcn', 'ecm1', 'clu','spp1','colec11', 'efemp1', 'lum', 'rgs5','cygb','gdf2'],
['igkc', 'cd79a', 'ighm', 'iglc2', 'iglc1', 'h2-aa', 'h2-ab1', 'h2-eb1'],

['cavin2', 'timp3', 'fabp4', 'ly6a', 'plvap', 'rspo3', 'cd9' ],
['s100a8', 's100a9', 'camp', 'wfdc21', 'pygl', 'lcn2', 'ngp'],
['ptprc', 'trbc2', 'cd3d', 'gimap3', 'ms4a4b','cd3g','ccl5','actb','lck'],
['ly6d', 'igkc', 'ighm', 'h2-aa', 'h2-eb1', 'cd79a', 'cd79b','mef2c','iglc3'],
['spp1', 'clu', 'tm4sf4', 'ly6a', 'cd63', 'plet1','plpp3','ddit4l','app','atp1a1'],

['apoa2', 'fabp1', 'mt1', 'saa1', 'saa2','saa3', 'serpina3k', 'lcn2','wfdc21','stard10','orm2'],
['nkg7', 'gzma', 'prf1', 'ptprcap', 'klrc2', 'anxa6', 'lck','aw112010','dcn','ccl5'],
['lyz2', 'h2-aa', 'h2-ab1', 'h2-eb1', 's100a6', 'actb', 'h2-dmb1'],
['gimap3', 'ltb', 'ms4a4b', 'dapl1', 'cd53', 'dbx2','actb', 'laptm5', 'svbp','ccdc63','tsga10ip','mboat1','cacna2d2','ikbkb'],
['actb', 'cd3g', 'ccl5', 'cd7', 'trbc2', 'adap1', 'elmo2','tmem156','txk','oas1c','tlr12'],

['fabp1', 'adk','fmo5', 'otc', 'cbs', 'fbp1', ],
['mup20', 'alb', 'apoa2', 'hsd17b13','apoa5','serpina3k','pck1','krt18','hpx','c8g','fxyd1'],
['ccl5', 'nkg7', 'trbc2', 'aw112010', 'trbv12-2', 'cd3g', 'nmb', 'cd8a','il2rb','setbp1']]


for gg in range(len(mygene)):
    print('\n\ncluster', gg, mygene[gg])
    for k in range(len(sheetnames)):
        if sheetnames[k].find('DEGs')!=-1:
            if sheetnames[k].find('Mouse')!=-1:
                out=findDEgenes(sheetnames[k],mygene[gg])
                print(sheetnames[k],out,sum(out))
