##########filter genes
a0=open('shared_genes.txt','r')
a1=open('HMM_CNV_predictions.HMMi3.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat','r')
w=open('HMM3_CNV_filtered.txt','w')

genes_list=[]
for line in a0:
    li=line.split()
    genes_list.append(li[0])

for line in a1:
    li=line.split()
    if li[3] in genes_list:
        w.write(line)
w.close()



