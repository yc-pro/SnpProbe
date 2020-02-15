#coding:gbk
'''
Created on 2020��1��15��

@author: Genome
'''

def sample_cl(file_in,file_out):
    fout = open(file_out,'w')
    ic = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
        ii = int(vec[1])
        if(ii >= 4):
            fout.write(">" + str(ic) + "\n" + vec[0]+"\n")
            ic += 1
    return 0

def x1_is_pair(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        if(line[0] != "@"):
            site = vec[2] +":"+vec[3]
            if(mm.has_key(site)):
                mm[site] += vec[0] +":"
            else:
                mm.setdefault(site,vec[0]+":")
    for key in mm:
        vec = mm[key][0:-1].split(":")
        if(len(vec) == 2):
            fout.write(key+"\t"+vec[0] + "\t"+vec[1]+"\n")
    return 0

def x1_pair_del(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[1],"")
        mm.setdefault(vec[2],"")
    id = ""
    for line in open(file_in1):
        line = line.strip()
        if(line[0] == ">"):
            id = line[1:]
        else:
            if(mm.has_key(id)):
                mm[id] += "line"
    for key in mm:
        fout.write(key+"\t"+mm[key]+"\n")
    return 0

'''
����Ѱ��reads���ȶԵ�ref�����ͳ�ƣ��ȶԵ�Ψһλ�õģ����һ����1���ȶԵ����λ�õ���2
������2��0����û����indel
'''
import os
import numpy as np

def is_only(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[2] + ":" + vec[3]    #��λ��referance�ϵıȶ�λ��
        fen = ""
        indel = ""
        
        for ii in range(4,len(vec)):
            ui = vec[ii]
            veckk = ui.split(":")
            if(ui.find("AS:i:") >= 0):  #�ҵ����������
                fen = veckk[2]  
            elif(ui.find("XO:i:") >= 0):    #�ҵ�indel�����
                indel = veckk[2] 
        sk = site+"@"+fen+"@"+ indel
        if(mm.has_key(id)): #idΪreads�ı��
            mm[id] += sk +"#"
        else:
            mm.setdefault(id,sk+"#")
    print(len(mm))
    for key in mm:
        vec = mm[key][0:-1].split("#")
        mm1 = {}
        for ii in range(0,len(vec)):
            vec1 = vec[ii].split("@")
            fen = int(vec1[1])
            if(mm1.has_key(fen)):
                mm1[fen] += 1   #�ԱȶԵ÷ֽ���ͳ�ƣ�����ж����ͬ�÷֣����1
            else:
                mm1.setdefault(fen,1)
        max_key = max(mm1) #�ҵ���С�÷�
        max_value = mm1[max_key]    #�õ���С�÷ֱȶԵ�referance��λ�õ���Ϣ
        if(max_value == 1): #��С�÷��Ƿ�ֻ��ref��һ��λ�ó���
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key): #�ҵ���С�÷ֵ���Ϣ
                    if(vec1[2] =="0"):  #�ж�ʱ���ǲ�����indel�ģ�vec1[0]λ�㣬vec1[1]�÷֣�vec1[2]indel���
                        fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1\t\n")
                    else:
                        fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1\t\n") #����indel
                
        else:
            #���������С�÷�
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key):
                    fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t2\t\n")
    return 0

'''
��ʼѰ����С�÷֣��������¼
'''
def is_minfen(mm):
    mm_re = {}
    for key in mm:
        vec = mm[key][0:-1].split("#")
        mm1 = {}
        for ii in range(0,len(vec)):
            vec1 = vec[ii].split("@")
            fen = int(vec1[1])
            if(mm1.has_key(fen)):
                mm1[fen] += 1   #�ԱȶԵ÷ֽ���ͳ�ƣ�����ж����ͬ�÷֣����1
            else:
                mm1.setdefault(fen,1)
        max_key = max(mm1) #�ҵ���С�÷�
        max_value = mm1[max_key]    #�õ���С�÷ֱȶԵ�referance��λ�õ���Ϣ
        if(max_value == 1): #��С�÷��Ƿ�ֻ��ref��һ��λ�ó���
            str_re = ""
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key): #�ҵ���С�÷ֵ���Ϣ
                    if(vec1[2] =="0"):  #�ж�ʱ���ǲ�����indel�ģ�vec1[0]λ�㣬vec1[1]�÷֣�vec1[2]indel���,���һ�е�1,����U
                        str_re = key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1"
                    else:
                        str_re = key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1" #����indel
            mm_re.setdefault(1,str_re)
                
        else:
            itcc = 1
            #���������С�÷�
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key):
                    str_re = key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t2"  #���һ�е�2����m
                    mm_re.setdefault(itcc,str_re)
                    itcc += 1
    return mm_re

'''
ͳ��ÿ��reads�ıȶ�������ó���С�÷ֵĸ������÷�ֵ��indel���
'''
def is_only_ordely(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    mm_true = {} #�����жϣ��Ƿ�reads�ȶ���˳���
    for line in open(file_in):
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[2] + ":" + vec[3]    #��λ��referance�ϵıȶ�λ��
        fen = ""
        indel = ""
        
        for ii in range(4,len(vec)):
            ui = vec[ii]
            veckk = ui.split(":")
            if(ui.find("AS:i:") >= 0):  #�ҵ����������
                fen = veckk[2]  
            elif(ui.find("XO:i:") >= 0):    #�ҵ�indel�����
                indel = veckk[2] 
        sk = site+"@"+fen+"@"+ indel
        
        #�����ǿ�ʼ�����ж�
        if(mm.has_key(id)): #idΪreads�ı��
            mm[id] += sk +"#"
        elif(len(mm) == 0):
            mm.setdefault(id,sk+"#")
        else:
            
            mm_re = is_minfen(mm)   #�ȰѾ���Ϣд�룬�ҵ���reads����С�÷֣������صȷ���Ϣ��map
            for key in mm_re:
                fout.write(mm_re[key] + "\t" + str(key) + "\n") #�����ļ�д��,�����С�÷֣����������
                
            old_id = min(mm)
            if(mm_true.has_key(old_id)):
                mm_true[old_id] += 1
            else:
                mm_true.setdefault(old_id,1)
                
            mm.clear()
            mm.setdefault(id,sk+"#") #���map����д������Ϣ
    
    #д�����һ��reads�ļ�¼
    mm_re = is_minfen(mm)
    for key in mm_re:
        fout.write(mm_re[key] + "\t" + str(key) + "\n")    
        
    for key in mm_true:
        fout1.write(str(key) + "\t" + str(mm_true[key]) + "\n")
    return 0

'''
Ѱ�������Ϣ
'''
def find_pair(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0] #reads��id
        site = vec[1]   #�ȶԵ�ref��λ�����
        many_site = vec[4]  #�Ƿ��Ƕ���ȶ�λ�㣬1ΪΨһһ��
        indel = vec[3]  #�Ƿ���indel��0��û����indel
        messa = many_site + "\t" + indel+"\t" + id
        if(mm.has_key(site)) :
            mm[site] += messa + "@"
        else:
            mm.setdefault(site,messa+"@")
    print("total site->" + len(mm))
    for key in mm:
        vec = mm[key][0:-1].split("@")
        ss = str(len(vec)) + "\t"
        for ii in range(0,len(vec)):
            ss += vec[ii] + "\t"
        fout.write(key + "\t" + ss + "\n")
    return 0
                
'''
��ȡreads�������ٴεıȶ�
'''
def ready_compair(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    mm1= {}
    for line in open(file_in):  
        #��ȡreads��id
        vec = line.strip().split("\t")
        id1 = vec[4]
        id2 = vec[7]
        mm.setdefault(id1,"")
        mm.setdefault(id2,"")
        mm1.setdefault(line.strip(),1)
    id = ""
    for line in open(file_in1):
        #����reads���е�Ѱ��
        line = line.strip()
        if(line[0] == ">"):
            id = line[1:]
        else:
            if(mm.has_key(id)):
                mm[id] = line
                
    for key in mm1:
        vec = key.split("\t")
        id1 = vec[4]
        id2 = vec[7]
        reads1 = ""
        reads2 = ""
        if(mm.has_key(id1) and mm.has_key(id2)):
            reads1 = mm[id1]
            reads2 = mm[id2]
            fout.write(key + "\t" + reads1 + "\t" + reads2+"\n")    #�����reads����Ϣ
            fout1.write(">" + vec[0]+":a\n" + reads1+"\n")
            fout1.write(">" + vec[0]+":b\n" + reads2+"\n")
        else:
            print(key)
        
    return 0

'''
�ԱȶԺ�������ٴν�����֤���������ȷ�ģ���ô���Ǹ�reads��û����indel����ô��snpλ�ú͸������浽map
����ԭ�е�ref��λ�ý�����֤�������ƥ��ģ��������ȷ����������ӡ����̨
'''
def yanzhen_compair(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    #Ѱ�ұȶ���Ϣ�е�snp��indel
    for line in open(file_in):  
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[2] + ":" + vec[3]    #��λ��referance�ϵıȶ�λ��
        snp = ""
        indel = ""
        ciga = ""
        
        for ii in range(4,len(vec)):
            ui = vec[ii]
            veckk = ui.split(":")
            if(ui.find("XM:i:") >= 0):  #�ҵ����������
                snp = veckk[2]  
            elif(ui.find("XO:i:") >= 0):    #�ҵ�indel�����
                indel = veckk[2] 
            elif(ui.find("MD:Z:") >= 0):    #�ҵ�indel�����
                ciga = veckk[2] 
        sk = site+"@"+snp+"@"+ciga
        if(mm.has_key(id)):
            print(id+"\tmany")
        elif(indel != "0"):
            print(id+"\tindel")
        else:
            mm.setdefault(id,sk)
    #��ԭ�е�λ�����ƥ��
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id1 = vec[0] + ":a"
        id2 = vec[0] + ":b"
        if(mm.has_key(id1) and mm.has_key(id2)):
            vec1 = mm[id1].split("@")
            vec2 = mm[id2].split("@")
            if(vec1[0] == vec2[0] == vec[0]):
                fout.write(line.strip()+"\t"+vec1[1]+"@"+vec1[2]+"\t"+vec2[1]+"@"+vec2[2]+"\n")
            else:
                print(line.strip()+"\t"+vec1[0]+"\t"+vec2[0]+"\tid_err")
        else:
            print(line.strip()+"\tid_notfind") 
    return 0
    
def yz_site(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    #Ѱ�ұȶ���Ϣ�е�snp��indel
    for line in open(file_in):  
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0][0:-2]
        site = vec[2] + ":" + vec[3]    #��λ��referance�ϵıȶ�λ��
        if(id != site):
            print(vec[0]+"\t"+site)
        else:
            reads = vec[9]
            for ii in range(4,len(vec)):
                ui = vec[ii]
                veckk = ui.split(":")
                if(ui.find("XM:i:") >= 0):  #�ҵ����������
                    snp = veckk[2]  
                elif(ui.find("XO:i:") >= 0):    #�ҵ�indel�����
                    indel = veckk[2] 
                elif(ui.find("MD:Z:") >= 0):    #�ҵ�indel�����
                    ciga = veckk[2] 
            sk = site+"@"+snp+"@"+ciga
            
def get_snp_outside(reads,ciga,site):
    vec_site = site.split(":")
    contig = vec_site[0]
    i_begin = int(vec_site[1])
    mm = {}
    for ii in range(0,len(ciga)):   #�õ�cigaÿһ��snp����λ��
        st_ci = ciga[ii]
        if(st_ci == "A" or st_ci == "T" or st_ci == "C" or st_ci == "G"):
            mm.setdefault(ii,st_ci)
    mm = sorted(mm.items(),key=lambda x:x[0],reverse=False) #��λ�ý�������
    mm_re = {}
    i_old = 0
    i_total = 0
    for key in mm:
        site = key[0]
        base = key[1]
        if(site == 0):
            st_ci = "0"
        else:
            st_ci = ciga[i_old:site]
        i_old = site + 1
        i_total += int(st_ci) + 1
        str_key = contig + ":" + str(i_begin+i_total-1)
        reads_base = reads[i_total-1]
        mm_re.setdefault(str_key,base+"@" + reads_base)
        #print(st_ci,base,i_total)
    #mm_re = sorted(mm_re.items(),key=lambda x:x[0],reverse=False)
    return mm_re

'''
����������
'''
def snp_outside(file_in,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    for line in open(file_in):  
        vec = line.strip().split("\t")
        if(line.find("N") > 0):
            fout1.write(line.strip() +  "\tN\n")    #����N�ļ�¼��ֱ�����������
            continue
        mm = {}
        reads1 = vec[8]
        reads2 = vec[9]
        vec1 = vec[10].split("@")
        vec2 = vec[11].split("@")
        i1 = int(vec1[0])
        i2 = int(vec2[0])
        ciga1 = vec1[1]
        ciga2 = vec2[1]
        str_snp_outside = ""
        if(i1 > 0):
            mm_temp = get_snp_outside(reads1,ciga1,vec[0])  #����snp��λ�ã��Լ�������
            for key in mm_temp:
                if(mm.has_key(key)):    #����map����һ��λ���£��ж��ٸ�������
                    mm[key] += mm_temp[key] + ":"   
                else:
                    mm.setdefault(key,mm_temp[key] + ":")
        if(i2 > 0):
            mm_temp = get_snp_outside(reads2,ciga2,vec[0])
            for key in mm_temp:
                if(mm.has_key(key)):
                    mm[key] += mm_temp[key] + ":"
                else:
                    mm.setdefault(key,mm_temp[key] + ":")
        for key in mm:
            vecio = mm[key][0:-1].split(":")
            mm_kk = {}
            str_kk = ""
            for ik in range(0,len(vecio)):
                message = vecio[ik]
                str_kk = message    #��������һ�ֻ�����
                if(mm_kk.has_key(message)):
                    mm_kk[message] += 1
                else:
                    mm_kk.setdefault(message,1)
            
            if(len(mm_kk) == 1):
                if(str_kk[-1] == "N"):
                    fout1.write(key + "\t" + str_kk + "\tN\n") 
                else:
                    fout.write(key + "\t" + str_kk + "\t1\n")
                    str_snp_outside += key + "#" + str_kk + "#1~"   #���һ��1�Ǳ�ʾ����ط�����1�ֻ�����
            else:
                vecpp = mm[key].split(":")
                str_snp_outside += key + "#" + mm[key] + "#" + str(len(vecpp)) +"~"#��һ����ʾλ�㣬�ڶ�����ʾ�����ͣ���3����ʾ���л�����
                fout.write(key + "\t" + mm[key] + "\t"+str(len(vecpp))+"\n")
        if(len(str_snp_outside) > 0):
            vecsnp = str_snp_outside[0:-1].split("~")
            fout2.write(vec[0] + "\t" + vec[8] + ":" + vec[9] + "\t" + vec[10] + ":" + vec[11] + "\t" + str_snp_outside +"\t" + str(len(vecsnp))+ "\n")
    return 0
       
'''
����������
'''
def snp_inside(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[1].split(":")
        vec2 = vec[0].split(":")
        contig = vec2[0]
        reads_begin = int(vec2[1])
        mm = {}
        for ii in range(0,len(vec1)):
            reads = vec1[ii]
            #print(len(reads),reads)
            for kk in range(0,len(reads)):
                base = reads[kk]
                if(kk == 59):
                    print(vec[0],base)
                
                if(mm.has_key(kk)):
                    mm[kk] += base
                else:
                    mm.setdefault(kk,base)
        ssk = ""
        i_snp = 0
        for key in mm:
            vv = mm[key]
            mm_temp = {}
            for kk in range(0,len(vv)):
                v_base = vv[kk]
                mm_temp.setdefault(v_base,1)
            str_snp = ""
            if(len(mm_temp) > 1):
                for kk in mm_temp:
                    str_snp += kk + "@"
                i_snp += 1
                str_snp = str_snp[0:-1]
                snp_site = reads_begin + key
                ssk += str(snp_site) + "#" + str(key) + "#" + str_snp+ "~"
        fout.write(line.strip() + "\t" + ssk + "\t" + str(i_snp) + "\n")
    return 0
 
'''
������,�޸�Ϊ������
'''
def snp_inside_ref(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[5][0:-1].split("~")
        vec2 = vec[0].split(":")
        contig = vec2[0]
        for ii in range(0,len(vec1)):
            vec3 = vec1[ii].split("#")
            site = vec3[0]
            id = contig+":"+site
            snp = vec3[2]
            if(mm.has_key(id)):
                mm[id] += snp + ":"
            else:
                mm.setdefault(id,snp+":")
    for key in mm:  #ͳ��ÿ��λ���ж���SNP
        vec = mm[key][0:-1].split(":")
        mm_temp = {}
        for ii in range(0,len(vec)):
            mm_temp.setdefault(vec[ii],1)
        str_snp = ""
        for kk in mm_temp:
            str_snp += kk+":"
        fout.write(key+"\t" +str_snp +"\t" + str(len(mm_temp)) + "\n")
    return 0

'''
ȥ��1������ıȶԼ�¼����sam�ļ�����
'''
def n_d0(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],0)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0][0:-2]
        if(mm.has_key(id)):
            fout.write(line.strip()+"\n")
        else:
            fout1.write(line.strip()+"\n")
    return 0

'''
ȥ������N159��������9��û��ԵĶ��Ӽ�¼����sam�ļ�����
'''
def n_d1(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        mm.setdefault(line.strip(),0)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0][0:-2]
        if(mm.has_key(id)):
            mm[id] += 1
        else:
            fout.write(line.strip()+"\n")
    for key in mm:
        if(mm[key] != 1):
            print(key)
    return 0

'''
��֤�£�������ת����λ�㣬�Ƿ񶼰����������������
'''
def yz_outside_inside(file_in,file_in1):
    mm = {}
    for line in open(file_in):
        vec = line.strip().split(" ")
        mm.setdefault(vec[1],"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            pass
        else:
            print(line.strip())
    return 0

'''
�õ�λ�������е�λ����Ϣ
'''
def step1(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split(":")
        print(line.strip());
        site = vec[1]
        val = os.popen("samtools view ../12_29/d.sort.bam " + line.strip() + "-" + site)
        for line in val.readlines():
            fout.write(line.strip()+"\n")
    return 0

'''
��ʼȥ��2�ۼ����ϵ�λ��
'''
def n_d2(file_in,file_in1,file_in2,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[0].split(":")
        id = vec1[0]+":"+vec1[1]
        mm.setdefault(id,0)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0][0:-2]
        if(mm.has_key(id)):
            mm[id] += 1
        else:
            fout.write(line.strip()+"\n")
    for line in open(file_in2):
        vec = line.strip().split("\t")
        id = vec[0]
        if(mm.has_key(id)):
            mm[id] += 1
        else:
            fout1.write(line.strip()+"\n")
    return 0

'''
ȥ������������͵�
@file_name,Ϊ�ȶԵ�bam�ļ�
'''
def step2(file_in,file_out,file_name):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split(":")
        site = vec[1]
        val = os.popen("samtools view " + file_name+" " + line.strip() + "-" + site)
        for line in val.readlines():
            fout.write(line.strip()+"\n")
    return 0


'''
�õ�ÿ��snpλ���£��������ٳɶԵ�Ƭ�Σ���һ�����Ѿ�ȷ����ÿ��snp�£�ֻ��һ�ֻ����͵�����µ�
'''
def get_snp_dp(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[5][0:-1].split("~")
        vec2 = vec[0].split(":")
        contig = vec2[0]
        for ii in range(0,len(vec1)):
            vec3 = vec1[ii].split("#")
            site = vec3[0]
            id = contig+":"+site
            snp = vec3[2]
            if(mm.has_key(id)):
                mm[id] += snp + ":"
            else:
                mm.setdefault(id,snp+":")
    for key in mm:
        vec = mm[key][0:-1].split(":")
        mm_temp = {}
        for ii in range(0,len(vec)):
            mm_temp.setdefault(vec[ii],1)
        str_snp = ""
        for kk in mm_temp:
            str_snp += kk+":"
        fout.write(key+"\t" +str_snp +"\t" + str(len(vec)) + "\n")
    return 0

'''
���Ҹ������
@file_name,Ϊ�ȶԵ�bam�ļ�
'''
def step3(file_in,file_out,file_out1,file_out2,file_name):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    mm = {}
    for line in open(file_in):
        vecaa = line.strip().split("\t")
        vec = vecaa[0].split(":")
        site = vec[1]
        icc = int(vecaa[2]) * 2
        acc=  0
        val = os.popen("samtools view " + file_name+".sort.bam " + vecaa[0] + "-" + site)
        mmk = {}
        for line_k in val.readlines():
            mmk.setdefault(line_k.strip(),1)
            acc+= 1
        fout.write(line.strip() + "\t" + str(acc) + "\t" + str(icc) + "\n")
        if(acc != icc):
            for kk in mmk:
                fout2.write(kk+"\n")
            fout1.write(line.strip() + "\t" + str(acc) + "\t" + str(icc) + "\n")
    return 0

'''
��ʼ��ȥ��Щ���д��ϵ�λ��
'''
def n_d3(file_in,file_in1,file_in2,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],0)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        veckk = vec[0].split(":")
        id = veckk[0] + ":" + veckk[1]
        if(mm.has_key(id)):
            mm[id] += 1
        else:
            fout.write(line.strip()+"\n")
    for line in open(file_in2):
        vec = line.strip().split("\t")
        id = vec[0]
        if(mm.has_key(id)):
            mm[id] += 1
        else:
            fout1.write(line.strip()+"\n")
    return 0

'''
ͳ��λ��֮��Ĳ�ֵ
'''
def tj_site_cha(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    mm3 = {}
    old_scaff = ""
    old_i_site = 0
    for line in open(file_in):
        vec = line.strip().split(":")
        scaff = vec[0]
        i_site = int(vec[1])
        id = line.strip()
        if(mm1.has_key(id)):
            mm1[id] += 1
        else:
            mm1.setdefault(id,1)
        
        if(old_scaff != scaff):
            old_scaff = scaff
            old_i_site = i_site
        elif(old_scaff == scaff):
            temp_vv = i_site - old_i_site 
            if(temp_vv > 0):
                ii_det = temp_vv / 1000;
                if(ii_det > 100):
                    print(line)
                if(mm.has_key(ii_det)):
                    mm[ii_det] += 1
                else:
                    mm.setdefault(ii_det,1)
                    
                if(mm3.has_key(ii_det)):
                    mm3[ii_det] += temp_vv
                else:
                    mm3.setdefault(ii_det,temp_vv)
                    
                old_i_site = i_site
    
    for key in mm:
        i_begin = key * 1000
        i_end = (key + 1) * 1000
        fout.write(str(i_begin) + "\t" + str(i_end) + "\t" + str(mm[key]) + "\t" + str(mm3[key])  + "\n")
    for key in mm1:
        ii = mm1[key]
        if(ii != 2):
            print(key+"\t" + str(ii))
    return 0

'''
ȥ��1������ıȶԼ�¼����sam�ļ�����
'''
def get_jia_outside(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],0)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        if(mm.has_key(id)):
            fout.write(line.strip()+"\n")
    return 0

'''
��ʼ��ѡ����
'''
def select_reads(file_in,file_out,file_name):
    fout = open(file_out,'w')
    for line in open(file_in):
        vecaa = line.strip().split("\t")
        vec = vecaa[0].split(":")
        site = vec[1]
        val = os.popen("samtools view " + file_name+".sort.bam " + vecaa[0] + "-" + site)
        mmk = {}
        mmreads = {}
        contig = vec[0]
        for line_k in val.readlines():
            vec_k = line_k.strip().split("\t")
            vec_k1 = vec_k[0].split(":")
            reads = vec_k[9]
            id = vec_k[0][0:-2]
            if(vec_k1[0] == vec[0]):
                site_k = int(vec_k1[1])
                if(mmk.has_key(site_k)):
                    mmk[site_k] += 1
                else:
                    mmk.setdefault(site_k,1)
                if(mmreads.has_key(id)):
                    mmreads[id] += reads + ":"
                else:
                    mmreads.setdefault(id,reads+":")
            else:
                print(line.strip());
        
        mmk = sorted(mmk.items(),key=lambda x:x[0],reverse=False) #��λ�ý�������
        i_len = len(mmk)
        i_t = i_len /2
        i_i = 0
        str_out = ""
        i_select = ""
        for key in mmk:
            site = key[0]
            base = key[1]
            if(base == 2):
                str_out += str(site) + ":"
                if(i_i == i_t):
                    i_select = contig + ":" + str(site)
            else:
                print(line.strip())
            i_i += 1
        select_reads = mmreads[i_select]
        fout.write(vecaa[0] + "\t" + vecaa[1] + "\t" + str(i_len) + "\t" + str_out + "\t" + i_select + "\t" + select_reads + "\n")
    return 0

'''
ͳ���Ӻϻ��Ǵ���
'''
def tj_1(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[1][0:-1])
        vec1 = vec[1][0:-1].split(":")
        mm1.setdefault(vec1[1],0)
        mm1.setdefault(vec1[0],0)
        
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        if(mm1.has_key(id)):
            mm1[id] = int(vec[1])
    
    for key in mm:
        vec = mm[key].split(":")
        i1 = mm1[vec[0]]
        i2 = mm1[vec[1]]
        str_snp = "N"
        if(i1 > 0 and i2 >0 ):
            str_snp = "AB"
        elif(i1 > 0 and i2 == 0):
            str_snp = "AA"
        elif(i1 == 0 and i2 > 0):
            str_snp = "BB"
        
        fout.write(key + "\t" + mm[key] + "\t" + str(i1) + ":" + str(i2) + "\t" + str_snp + "\n")
    return 0

'''
��50����ͳ����Ϣ����
file_in��50���������
file_name��·������
'''
def tj_2(file_in,file_in1,file_name,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm1.setdefault(vec[0],"")
        mm.setdefault(vec[0],"")
        
    for linekk in open(file_in1):
        linekk = linekk.strip()
        file_aa = file_name + "/" + linekk + ".tj"
        for line in open(file_aa):
            vec = line.strip().split("\t")
            vec1 = vec[2].split(":")
            i_total = int(vec1[0]) + int(vec1[1])
            id = vec[0]
            if(mm1.has_key(id)):
                mm[id] += str(i_total) + ","
                mm1[id] += vec[2] + ","
    
    for key in mm:
        ik = mm1[key]
        fout.write(key + "\t" + mm[key] + "\t" + ik + "\n")
    return 0

'''
�����ֵ,����ͱ�׼��
'''
def reads_mean(file_in,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    mm = {}
    mmstd = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[1][0:-1].split(",")
        ll = []
        for key in range(0,len(vec1)):
            kk = int(vec1[key])
            ll.append(kk)
        while(len(ll) < 50):
            ll.append(0)
        i_mean = np.mean(ll)#��ֵ
        i_var = np.var(ll,ddof=1)#����
        i_std = np.std(ll,ddof=1)#��׼��
        ii = round(i_mean,1)    #��������
        ii_std = round(i_std,1)
        if(mm.has_key(ii)):
            mm[ii] += 1
        else:
            mm.setdefault(ii,1)
        if(mmstd.has_key(ii_std)):
            mmstd[ii_std] += 1
        else:
            mmstd.setdefault(ii_std,1)            
            
        fout.write(line.strip()+"\t"+str(i_mean)+"\t"+str(i_var)+"\t"+str(ii_std)+"\n")
    for key in mm:
        fout1.write(str(key)+"\t"+str(mm[key])+"\n");
    print("===")
    for key in mmstd:
        fout2.write(str(key)+"\t"+str(mmstd[key])+"\n");    
    return 0

def final_pair(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,"")
        
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec1 = vec[3].split("#")
        id = vec1[0]
        if(mm.has_key(id)):
            fout.write(line.strip()+"\t")
        else:
            fout1.write(line.strip()+"\t")
    return 0

def tj_overlap(file_in,file_out):
    fout = open(file_out,'w')
    i_sum = 0
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split("#")
        mm = {}
        for ii in range(0,len(vec1)):
            vec2 = vec1[ii].split(":")
            i_site = int(vec2[1])
            mm.setdefault(i_site,0)
        max_key = max(mm) #�ҵ���С�÷�
        min_key = min(mm)
        i_det = max_key - min_key + 60
        i_sum += i_det
        fout.write(line.strip() + "\t" + str(min_key) + "\t" + str(max_key) + "\t" + str(i_det) + "\n")
    print(i_sum)
    return 0

'''
�õ�ÿ��snpλ���£��������ٳɶԵ�Ƭ�Σ���һ�����Ѿ�ȷ����ÿ��snp�£�ֻ��һ�ֻ����͵�����µ�
'''
def snp_dp(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[5][0:-1].split("~")
        vec2 = vec[0].split(":")
        contig = vec2[0]
        
        mm1 = {} #�жϸ�reads���Ա��Ƿ�����ref��һ�µĴ�����
        veck1 = vec[3][0:-1].split("~")
        for ik in range(0,len(veck1)):
            ik_id = veck1[ik].split("#")[0]
            mm1.setdefault(ik_id,"N")
            
        for ii in range(0,len(vec1)):
            vec3 = vec1[ii].split("#")
            site = vec3[0]
            id = contig+":"+site
            snp = vec3[2]
            if(mm1.has_key(id)):
                mm1[id] = snp
            else:
                mm1.setdefault(id,snp)
        
        for key in mm1:
            snp = mm1[key]
            if(mm.has_key(key)):
                mm[key] += snp + ":"
            else:
                mm.setdefault(key,snp+":")

                    
                    
    for key in mm:
        vec = mm[key][0:-1].split(":")
        mm_temp = {}
        for ii in range(0,len(vec)):
            mm_temp.setdefault(vec[ii],1)
        str_snp = ""
        iK = len(vec)
        for kk in mm_temp:
            str_snp += kk+":"
            if(kk == "N"):
                ik = 0 #����Ǻ����Ա��Ƿ�����ref��һ�µĴ����ͣ���ô�Ѹ�λ�����Ϣȥ��������λ��ɨ���ʱ��ͨ���ж���ȣ���ΪΪ0�����Ը�λ���ȥ��
        fout.write(key+"\t" +str_snp +"\t" + str(ik) + "\n")
    return 0

def x5_fa(file_in,file_out):
    fout = open(file_out,'w')
    id = ""
    for line in open(file_in):
        line = line.strip()
        if(line[0] == ">"):
            id = line[1:]
        else:
            fout.write(">"+id+"_a\n" + line+"\n")
            fout.write(">"+id+"_b\n" + line+"\n")
            fout.write(">"+id+"_c\n" + line+"\n")
            fout.write(">"+id+"_d\n" + line+"\n")
            fout.write(">"+id+"_e\n" + line+"\n")
    return 0

if __name__ == '__main__':
    #print("Lachesis_group18__9_contigs__length_31774926:23659499_a"[0:-2])
    #snp_outside('E://super_down//tt','E://super_down//tt.he','E://super_down//tt.err','E://super_down//kk.err')
    #snp_inside('E://super_down//kt','E://super_down//kt.he')
    #tj_site_cha('E://super_down//lo','E://super_down//lo.he')
    i = 20.3421417798 
    ii = round(i,1)
    ss = ">12"
    print(ss[1:])
    #tj_overlap('/Users/yangcheng/Documents/kk','/Users/yangcheng/Documents/kk.a')
    #snp_dp('/Users/yangcheng/Documents/kl','/Users/yangcheng/Documents/kl.a')
    #reads_mean('/Users/yangcheng/Documents/kk','/Users/yangcheng/Documents/kk.a')
    #mean_sd('/Users/yangcheng/Documents/log_mean')
    #get_snp_outside('ACAATTGTTTCACAATGGCAATCTTTTGTCTGGGACAGCCAAAAATCGTGCAGTGTATCC','A46G2C3A4T','Lachesis_group6__25_contigs__length_38613331:26189578')