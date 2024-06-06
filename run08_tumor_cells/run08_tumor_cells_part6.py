def output_count(file_in,file_out):
    a1=open(file_in,'r')
    w=open(file_out,'w')
    dict_cell2count={}
    for line in a1:
        li=line.split()
        if li[0] not in dict_cell2count:
            count=0
            if li[1] !='2':
                count=count+1
                dict_cell2count[li[0]]=count
        else:
            if li[1] != '2':
                dict_cell2count[li[0]]=dict_cell2count[li[0]]+1
    for cell in dict_cell2count:
        w.write(cell+'\t'+str(dict_cell2count[cell])+'\n')
    w.close()

output_count('cells_cnv_normal.txt','cells_cnv_percentage_normal.txt')
output_count('cells_cnv_nsclc.txt','cells_cnv_percentage_nsclc.txt')
output_count('cells_cnv_sclc.txt','cells_cnv_percentage_sclc.txt')



