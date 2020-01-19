#coding:gbk
'''
Created on 2020��1��15��

@author: Genome
'''

'''
����Ѱ��reads���ȶԵ�ref�����ͳ�ƣ��ȶԵ�Ψһλ�õģ����һ����1���ȶԵ����λ�õ���2
������2��0����û����indel
'''
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
                fout.write(mm_re[key] + "\t" + str(key) + "\n") #�����ļ�д��
                
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
            fout1.write(">" + id1+"\n" + reads1+"\n")
            fout1.write(">" + id2+"\n" + reads2+"\n")
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
        id1 = vec[4]
        id2 = vec[7]
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
    
if __name__ == '__main__':
    print(">hello word"[1:])