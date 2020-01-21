#coding:gbk
'''
Created on 2020年1月15日

@author: Genome
'''

'''
进行寻找reads，比对到ref情况的统计，比对到唯一位置的，最后一列是1，比对到多个位置的是2
倒数第2列0代表没发生indel
'''
def is_only(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[2] + ":" + vec[3]    #定位到referance上的比对位置
        fen = ""
        indel = ""
        
        for ii in range(4,len(vec)):
            ui = vec[ii]
            veckk = ui.split(":")
            if(ui.find("AS:i:") >= 0):  #找到分数的情况
                fen = veckk[2]  
            elif(ui.find("XO:i:") >= 0):    #找到indel的情况
                indel = veckk[2] 
        sk = site+"@"+fen+"@"+ indel
        if(mm.has_key(id)): #id为reads的编号
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
                mm1[fen] += 1   #对比对得分进行统计，如果有多个相同得分，则加1
            else:
                mm1.setdefault(fen,1)
        max_key = max(mm1) #找到最小得分
        max_value = mm1[max_key]    #得到最小得分比对到referance的位置的信息
        if(max_value == 1): #最小得分是否只在ref的一个位置出现
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key): #找到最小得分的信息
                    if(vec1[2] =="0"):  #判断时候是不包含indel的，vec1[0]位点，vec1[1]得分，vec1[2]indel情况
                        fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1\t\n")
                    else:
                        fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1\t\n") #包含indel
                
        else:
            #包含多个最小得分
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key):
                    fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t2\t\n")
    return 0

'''
开始寻找最小得分，并输出记录
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
                mm1[fen] += 1   #对比对得分进行统计，如果有多个相同得分，则加1
            else:
                mm1.setdefault(fen,1)
        max_key = max(mm1) #找到最小得分
        max_value = mm1[max_key]    #得到最小得分比对到referance的位置的信息
        if(max_value == 1): #最小得分是否只在ref的一个位置出现
            str_re = ""
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key): #找到最小得分的信息
                    if(vec1[2] =="0"):  #判断时候是不包含indel的，vec1[0]位点，vec1[1]得分，vec1[2]indel情况,最后一列的1,就是U
                        str_re = key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1"
                    else:
                        str_re = key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t1" #包含indel
            mm_re.setdefault(1,str_re)
                
        else:
            itcc = 1
            #包含多个最小得分
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key):
                    str_re = key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\t2"  #最后一列的2就是m
                    mm_re.setdefault(itcc,str_re)
                    itcc += 1
    return mm_re

'''
统计每条reads的比对情况，得出最小得分的个数，得分值，indel情况
'''
def is_only_ordely(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    mm_true = {} #用于判断，是否reads比对是顺序的
    for line in open(file_in):
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[2] + ":" + vec[3]    #定位到referance上的比对位置
        fen = ""
        indel = ""
        
        for ii in range(4,len(vec)):
            ui = vec[ii]
            veckk = ui.split(":")
            if(ui.find("AS:i:") >= 0):  #找到分数的情况
                fen = veckk[2]  
            elif(ui.find("XO:i:") >= 0):    #找到indel的情况
                indel = veckk[2] 
        sk = site+"@"+fen+"@"+ indel
        
        #以下是开始进行判断
        if(mm.has_key(id)): #id为reads的编号
            mm[id] += sk +"#"
        elif(len(mm) == 0):
            mm.setdefault(id,sk+"#")
        else:
            
            mm_re = is_minfen(mm)   #先把旧信息写入，找到该reads的最小得分，并返回等分信息的map
            for key in mm_re:
                fout.write(mm_re[key] + "\t" + str(key) + "\n") #进行文件写入
                
            old_id = min(mm)
            if(mm_true.has_key(old_id)):
                mm_true[old_id] += 1
            else:
                mm_true.setdefault(old_id,1)
                
            mm.clear()
            mm.setdefault(id,sk+"#") #清空map后再写入新信息
    
    #写人最后一个reads的记录
    mm_re = is_minfen(mm)
    for key in mm_re:
        fout.write(mm_re[key] + "\t" + str(key) + "\n")    
        
    for key in mm_true:
        fout1.write(str(key) + "\t" + str(mm_true[key]) + "\n")
    return 0

'''
寻找配对信息
'''
def find_pair(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0] #reads的id
        site = vec[1]   #比对到ref的位点情况
        many_site = vec[4]  #是否是多个比对位点，1为唯一一个
        indel = vec[3]  #是否发生indel，0是没发生indel
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
提取reads，进行再次的比对
'''
def ready_compair(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    mm1= {}
    for line in open(file_in):  
        #提取reads的id
        vec = line.strip().split("\t")
        id1 = vec[4]
        id2 = vec[7]
        mm.setdefault(id1,"")
        mm.setdefault(id2,"")
        mm1.setdefault(line.strip(),1)
    id = ""
    for line in open(file_in1):
        #进行reads序列的寻找
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
            fout.write(key + "\t" + reads1 + "\t" + reads2+"\n")    #输出含reads的信息
            fout1.write(">" + vec[0]+":a\n" + reads1+"\n")
            fout1.write(">" + vec[0]+":b\n" + reads2+"\n")
        else:
            print(key)
        
    return 0

'''
对比对后的数据再次进行验证，如果是正确的，那么就是该reads并没发生indel，那么将snp位置和个数保存到map
再与原有的ref的位置进行验证，如果是匹配的，这进行正确输出，否则打印到后台
'''
def yanzhen_compair(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    #寻找比对信息中的snp，indel
    for line in open(file_in):  
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[2] + ":" + vec[3]    #定位到referance上的比对位置
        snp = ""
        indel = ""
        ciga = ""
        
        for ii in range(4,len(vec)):
            ui = vec[ii]
            veckk = ui.split(":")
            if(ui.find("XM:i:") >= 0):  #找到分数的情况
                snp = veckk[2]  
            elif(ui.find("XO:i:") >= 0):    #找到indel的情况
                indel = veckk[2] 
            elif(ui.find("MD:Z:") >= 0):    #找到indel的情况
                ciga = veckk[2] 
        sk = site+"@"+snp+"@"+ciga
        if(mm.has_key(id)):
            print(id+"\tmany")
        elif(indel != "0"):
            print(id+"\tindel")
        else:
            mm.setdefault(id,sk)
    #与原有的位点进行匹配
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
    #寻找比对信息中的snp，indel
    for line in open(file_in):  
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0][0:-2]
        site = vec[2] + ":" + vec[3]    #定位到referance上的比对位置
        if(id != site):
            print(vec[0]+"\t"+site)
        else:
            reads = vec[9]
            for ii in range(4,len(vec)):
                ui = vec[ii]
                veckk = ui.split(":")
                if(ui.find("XM:i:") >= 0):  #找到分数的情况
                    snp = veckk[2]  
                elif(ui.find("XO:i:") >= 0):    #找到indel的情况
                    indel = veckk[2] 
                elif(ui.find("MD:Z:") >= 0):    #找到indel的情况
                    ciga = veckk[2] 
            sk = site+"@"+snp+"@"+ciga
            
def get_snp_outside(reads,ciga,site):
    vec_site = site.split(":")
    contig = vec_site[0]
    i_begin = int(vec_site[1])
    mm = {}
    for ii in range(0,len(ciga)):   #得到ciga每一个snp所在位置
        st_ci = ciga[ii]
        if(st_ci == "A" or st_ci == "T" or st_ci == "C" or st_ci == "G"):
            mm.setdefault(ii,st_ci)
    mm = sorted(mm.items(),key=lambda x:x[0],reverse=False) #对位置进行排序
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
计算外坐标
'''
def snp_outside(file_in,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    for line in open(file_in):  
        vec = line.strip().split("\t")
        if(line.find("N") > 0):
            fout1.write(line.strip() +  "\tN\n")    #含有N的记录，直接输出到错误
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
            mm_temp = get_snp_outside(reads1,ciga1,vec[0])  #计算snp的位置，以及基因型
            for key in mm_temp:
                if(mm.has_key(key)):    #放入map，看一个位点下，有多少个基因型
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
                str_kk = message    #用来带出一种基因型
                if(mm_kk.has_key(message)):
                    mm_kk[message] += 1
                else:
                    mm_kk.setdefault(message,1)
            
            if(len(mm_kk) == 1):
                if(str_kk[-1] == "N"):
                    fout1.write(key + "\t" + str_kk + "\tN\n") 
                else:
                    fout.write(key + "\t" + str_kk + "\t1\n")
                    str_snp_outside += key + "#" + str_kk + "#1~"   #最后一个1是表示这个地方，有1种基因型
            else:
                vecpp = mm[key].split(":")
                str_snp_outside += key + "#" + mm[key] + "#" + str(len(vecpp)) +"~"#第一个表示位点，第二个表示基因型，第3个表示集中基因型
                fout.write(key + "\t" + mm[key] + "\t"+str(len(vecpp))+"\n")
        if(len(str_snp_outside) > 0):
            vecsnp = str_snp_outside[0:-1].split("~")
            fout2.write(vec[0] + "\t" + vec[8] + ":" + vec[9] + "\t" + vec[10] + ":" + vec[11] + "\t" + str_snp_outside +"\t" + str(len(vecsnp))+ "\n")
    return 0
       
'''
计算内坐标
'''
def snp_inside(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[1][0:-1].split(":")
        vec2 = vec[0].split(":")
        contig = vec2[0]
        reads_begin = int(vec2[1])
        mm = {}
        for ii in range(0,len(vec1)):
            reads = vec1[ii]
            for kk in range(0,len(reads)):
                base = reads[kk]
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
 
if __name__ == '__main__':
    #print("Lachesis_group18__9_contigs__length_31774926:23659499_a"[0:-2])
    #snp_outside('E://super_down//tt','E://super_down//tt.he','E://super_down//tt.err','E://super_down//kk.err')
    snp_inside('E://super_down//kt','E://super_down//kt.he')
    #get_snp_outside('ACAATTGTTTCACAATGGCAATCTTTTGTCTGGGACAGCCAAAAATCGTGCAGTGTATCC','A46G2C3A4T','Lachesis_group6__25_contigs__length_38613331:26189578')