def extract_cluster2cells(file_in):
    a=open(file_in,'r')

    dict_cluster2cells={}
    for line in a:
        li=line.split()
        cluster=li[0]
        cell=li[1]
        if cluster not in dict_cluster2cells:
            list=[]
            list.append(li[1])
            dict_cluster2cells[cluster] = list
        else:
            dict_cluster2cells[cluster].append(li[1])
    return dict_cluster2cells

dict_normal = extract_cluster2cells('normal_cluster2cells.txt')
dict_sclc = extract_cluster2cells('sclc_cluster2cells.txt')
dict_nsclc = extract_cluster2cells('nsclc_cluster2cells.txt')


a1=open('HMM3_CNV_filtered.txt','r')
w1=open('cells_cnv_sclc.txt','w')
w2=open('cells_cnv_nsclc.txt','w')
w3=open('cells_cnv_normal.txt','w')

normal_cells=[]
sclc_cells=[]
nsclc_cells=[]

for line in a1:
    li=line.split()
    cluster = li[0].split('.')[1]
    if cluster.startswith('sclc'):
        for cell in dict_sclc[cluster]:
            w1.write(cell+'\t'+li[2]+'\t'+li[3]+'\n')
    if cluster.startswith('nsclc'):
        for cell in dict_nsclc[cluster]:
            w2.write(cell+'\t'+li[2]+'\t'+li[3]+'\n')
    if cluster.startswith('normal'):
        for cell in dict_normal[cluster]:
            w3.write(cell+'\t'+li[2]+'\t'+li[3]+'\n')

w1.close()
w2.close()
w3.close()



