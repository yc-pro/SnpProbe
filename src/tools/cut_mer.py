#coding:gbk
#coding:utf-8
'''
Created on 2019年4月26日

@author: gene
'''
import string

def cut_reads(reads,len_r):
    mm_re = {}
    for ii in range(0,len(reads)-len_r):
        s_temp = reads[ii:ii+2]
        s_sub = reads[ii:ii+len_r]
        if(s_temp == "AG" or s_temp == "AT" or s_temp == "AC" or s_temp == "GA" or s_temp == "GT"  or s_temp == "GC"):
            if(len(s_sub) == len_r):
                mm_re.setdefault(ii,s_sub)
    return mm_re

def dd_reverse(str_in):
    trans_table = string.maketrans('ATCG','TAGC')
    a1 = str_in.translate(trans_table)
    str_r1 = a1[::-1]
    return str_r1
       
def cut(file_in1,file_in2,file_out):
    i_max = 130000000
    #i_max = 13
    cut_len = 60
    i_cc = 0
    i_len = 149
    while(i_cc <= i_max):
        with open(file_in1) as fp1, open(file_in2) as fp2: 
            i_total = 1;
            i_current_len = 0
            for line_w1 in fp1: 
                line_w2= fp2.readline()
                line_r1 = line_w1.strip()
                line_r2 = line_w2.strip()
                if(i_total == 2):
                    i_current_len += 1
                    r_line_r1 = dd_reverse(line_r1)
                    r_line_r2 = dd_reverse(line_r2)
                    if(i_cc <= i_max):
                        if(len(line_r1) >= i_len):
                            i_cc += 2
                            mm_re_1 = cut_reads(line_r1,cut_len)
                            mm_re_2 = cut_reads(r_line_r1,cut_len)
                            for site in mm_re_1:
                                key = mm_re_1[site][0:4]
                                fout = open(file_out+key,'a')
                                fout.write(">"+str(i_current_len)+"_"+str(site)+"_1+\t" + mm_re_1[site] + "\n")
                            for site in mm_re_2:
                                key = mm_re_2[site][0:4]
                                fout = open(file_out+key,'a')
                                fout.write(">"+str(i_current_len)+"_"+str(site)+"_1-\t" + mm_re_2[site] + "\n")
                        if(len(r_line_r2) >= i_len):
                            i_cc += 2
                            mm_re_1 = cut_reads(line_r2,cut_len)
                            mm_re_2 = cut_reads(r_line_r2,cut_len)   
                            for site in mm_re_1:
                                key = mm_re_1[site][0:4]
                                fout = open(file_out+key,'a')
                                fout.write(">"+str(i_current_len)+"_"+str(site)+"_2+\t" + mm_re_1[site] + "\n")
                            for site in mm_re_2:
                                key = mm_re_2[site][0:4]
                                fout = open(file_out+key,'a')
                                fout.write(">"+str(i_current_len)+"_"+str(site)+"_2-\t" + mm_re_2[site] + "\n")
                    else:
                        break;
                elif(i_total == 4):
                    i_total = 0
                i_total += 1
            
    return 0

def create_256():
    mm = {}
    mm.setdefault("A",1)
    mm.setdefault("T",1)
    mm.setdefault("C",1)
    mm.setdefault("G",1)
    for key1 in mm:
        for key2 in mm:
            for key3 in mm:
                for key4 in mm:
                    s1 = key1 + key2+ key3 + key4
                    st = "fout_"+s1+" = open(file_out+\""+s1+"\",'w')"
                    print(st)
                    print("mm.setdefault(\""+s1+"\",fout_"+s1+")")


def test():
    for ii in range(0,100):
        fout = open("D://down1//test.txt",'a')
        fout.write(str(ii)+"\n")
    return 0



if __name__ == '__main__':
    #pp_cs('D://down1//s_k1','D://down1//s_k','D://down1//s_k2')
    #get_snp('D://down1//tt','D://down1//tt_k')
    #reduance_snp('D://down1//tt_k','D://down1//tt_w')
    #det_tj('D://down1//4.w','D://down1//4.d')
    #cut('/Users/yangcheng/Documents/超算数据/t_1.fq','/Users/yangcheng/Documents/超算数据/t_2.fq','/Users/yangcheng/Documents/超算数据/t_cut/')
    '''
    mm = cut_reads("ATGGAAGCTAAAATACAGGAAGTTACGCTTACTCAGTAGTCCTTGAACTGATTTTTTATATTTGAGTGCACATTCAACTTAATCTCCATTTAATATTTGAAAAAAAGAAATAGATTTGTGAAAGGCTGCAAAAAAAGAAATTGCCCCTGAC",60)
    for key in mm:
        print(key,mm[key])
    '''
    cut('D://down1//t1','D://down1//t2','D://down1//t_cut')