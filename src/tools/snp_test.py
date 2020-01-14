#coding:UTF-8
'''
Created on 

@author: hankou
'''
from collections import Counter
import operator
#import scipy.stats as st
import numpy as np
from _socket import AI_ADDRCONFIG
import math
import os
import string

def pp(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split(" ")
        mm.setdefault(vec1[1],1)
        
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        id = vec1[2] + ":" + vec1[3]
        if(mm.has_key(id)):
            q_id = vec1[0]
            q_reads = vec1[9]
            #fout.write(">" + q_id + "\n")
            #fout.write(q_reads + "\n")
            fout.write(id + "\t" + q_id + "\t" + q_reads + "\n")
    return 0

def double_ab(file_in,file_out,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out2,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split(" ")
        id = vec1[0]
        s1 = vec1[1] + ":" + vec1[2] + "\t"
        if(mm.has_key(id)):
            mm[id] += s1
        else:
            mm.setdefault(id,s1)
        
    for key in mm:
        vec1 = mm.get(key)[0:-1].split("\t")
        if(len(vec1) == 2):
            id1 = vec1[0].split(":")[0]
            reads1 = vec1[0].split(":")[1]
            id2 = vec1[1].split(":")[0]
            reads2 = vec1[1].split(":")[1]
            fout.write(">" + id1 + "_a\n" + reads1 + "\n")
            fout.write(">" + id1 + "_b\n" + reads1 + "\n")
            fout.write(">" + id2 + "_a\n" + reads2 + "\n")
            fout.write(">" + id2 + "_b\n" + reads2 + "\n")
            fout1.write(key + "\t" + id1 + "\t" + id1 + "_a\t" + id1 + "_b\t" + id2  + "\t" + id2+ "_a\t" + id2 + "_b\t\n")
        else:
            print(key)
    return 0

def pp_cs(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
    #print(vec1[1])
        mm.setdefault(line,1)
        
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
    #print(vec1[0],vec1[1]);
        id = vec1[2] + ":" + vec1[3]
        if(mm.has_key(id)):
            q_id = vec1[0]
            q_reads = vec1[9]
            fout.write(id + "\t" + q_id + "\t" + q_reads + "\n")
    return 0

def tiqu_4_sam(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
        
    for line in open(file_in1):
        if(line[0] != "@"):
            vec1 = line.strip().split("\t")
            id = vec1[2] + ":" + vec1[3]
            if(mm.has_key(id)):
                fout.write(line.strip() + "\n")
    return 0

def get_snp(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        line = line.strip()
        vec1 = line.split("\t")
        i0 = line.find("NM:i:0")
        i1 = line.find("NM:i:1")
        if(i1 > 0):
            vec2 = vec1[len(vec1) - 3].split(":")
            pd = vec2[2]
            w_site = 0
            if(pd.find("A") > 0):
                w_site = pd.find("A")
            elif(pd.find("G") > 0):
                w_site = pd.find("G")
            elif(pd.find("C") > 0):
                w_site = pd.find("C")
            elif(pd.find("T") > 0):
                w_site = pd.find("T")
            s_site = int(pd[0:w_site])
            i_snp = int(vec1[3]) + s_site
            
            s_ref_jj = pd[w_site]
            s_snp_jj = vec1[9][s_site]
            fout.write(vec1[2]+":" +str(i_snp)+"\t" +s_ref_jj+":" +s_snp_jj + "\n")
        elif(i0 > 0):
            pass
        else:
            print(line)
        
    return 0

def get_snp_reads(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        line = line.strip()
        vec1 = line.split("\t")
        reads = vec1[9]
        i0 = line.find("NM:i:0")
        i1 = line.find("NM:i:1")
        if(i1 > 0):
            vec2 = vec1[len(vec1) - 3].split(":")
            pd = vec2[2]
            w_site = 0
            if(pd.find("A") > 0):
                w_site = pd.find("A")
            elif(pd.find("G") > 0):
                w_site = pd.find("G")
            elif(pd.find("C") > 0):
                w_site = pd.find("C")
            elif(pd.find("T") > 0):
                w_site = pd.find("T")
            s_site = int(pd[0:w_site])
            i_snp = int(vec1[3]) + s_site
            
            s_ref_jj = pd[w_site]
            s_snp_jj = vec1[9][s_site]
            
            s_1 = reads[0:s_site]
            s_2 = reads[s_site + 1:len(reads)]
            reads_ref = s_1 + s_ref_jj + s_2
            reads_snp = reads
            id = vec1[2]+":" +str(i_snp);
            fout.write(id + "\t" +reads_ref + ":" + reads_snp + ":" + s_ref_jj+":" +s_snp_jj +":" + pd +":" +vec1[3]+ "\n")
        elif(i0 > 0):
            pass
        else:
            print(line)
        
    return 0

def reduance_snp(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        id = vec1[0]
        snp = vec1[1] + "\t"
        if(not mm.has_key(id)):
            mm.setdefault(id,snp)
        else:
            mm[id] += snp
    for key in mm:
        vec1 = mm.get(key)[0:-1].split("\t")
        mm1 = {}
        for key1 in vec1:
            mm1.setdefault(key1,"1")
        fout.write(key + "\t" + str(len(vec1)) + "\t" +str(len(mm1))+ "\t"+ mm.get(key) + "\n")
    return 0

def det(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        vec2 = vec1[0].split(":")
        id = vec2[0]
        site = vec2[1]+":"
        if(not mm.has_key(id)):
            mm.setdefault(id,site)
        else:
            mm[id] += site
    for key in mm:
        vec1 = mm.get(key)[0:-1].split(":")
        mm1 = {}
        for key1 in vec1:
            mm1.setdefault(int(key1),int(key1))
        mm_temp = sorted(mm1.items(), key = operator.itemgetter(1), reverse = False)
        #print(mm_temp)
        old_vv = 0
        max_vv = 0
        min_vv = 1000000000000
        ii = 1
        for iw in range(0,len(mm_temp)):
            vv = mm_temp[iw][1]
            if(ii == 1):
                ii += 1
                old_vv = vv
            else:
                temp_vv = vv - old_vv
                old_vv = vv
                if(temp_vv >= max_vv):
                    #print(vv)
                    max_vv = temp_vv
                if(temp_vv <= min_vv):
                    print(vv)
                    min_vv = temp_vv
        fout.write(key+"\t"+str(len(vec1))+"\t"+str(min_vv)+"\t"+str(max_vv)+"\n")
    return 0


def det_tj(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        vec2 = vec1[0].split(":")
        id = vec2[0]
        site = vec2[1]+":"
        if(not mm.has_key(id)):
            mm.setdefault(id,site)
        else:
            mm[id] += site
    mm1 = {}
    mm2 = {}
    for key in mm:
        vec1 = mm.get(key)[0:-1].split(":")
        for key1 in vec1:
            mm1.setdefault(int(key1),int(key1))
        mm_temp = sorted(mm1.items(), key = operator.itemgetter(1), reverse = False)
        #print(mm_temp)
        old_vv = 0
        ii = 1
        for iw in range(0,len(mm_temp)):
            vv = mm_temp[iw][1]
            if(ii == 1):
                ii += 1
                old_vv = vv
            else:
                temp_vv = vv - old_vv
                old_vv = vv
                ii_det = temp_vv / 1000;
                print(temp_vv,ii_det)
                if(mm2.has_key(ii_det)):
                    mm2[ii_det] += 1
                else:
                    mm2.setdefault(ii_det,1)
    for key in mm2:
        fout.write(str(key)+"\t"+str(mm2[key]) + "\n")
    return 0

def tiqu_snp(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec1 = line.strip().split(":")
        reads = vec1[2]
        pd = vec1[len(vec1) - 1]
        w_site = 0
        if(pd.find("A") > 0):
            w_site = pd.find("A")
        elif(pd.find("G") > 0):
            w_site = pd.find("G")
        elif(pd.find("C") > 0):
            w_site = pd.find("C")
        elif(pd.find("T") > 0):
            w_site = pd.find("T")
        s_site = int(pd[0:w_site])
        s_ref_jj = pd[w_site]
        s_snp_jj = reads[s_site]
        s_out_snp = vec1[0] + ":" + vec1[1] + ":" + reads +  ":" + str(s_site) +":" + s_ref_jj + ":" + s_snp_jj;
        fout.write(s_out_snp+"\n")
    return 0

def zl_snp(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0],"")
    for line in open(file_in1):
        vec1 = line.strip().split(":")
        id = vec1[0] + ":" + vec1[1]
        if(mm.has_key(id)):
            mm[id] = vec1[2]
            
    for key in mm:
        value = mm[key]
        if(value != "+"):
            fout.write(key+"\t" + value + "\n")
        else:
            print(key + "\n")
    return 0

def get_snp_ab_reads(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0],"+")
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        id = vec1[0]
        if(mm.has_key(id)):
            mm[id] = vec1[1]
            
    for key in mm:
        value = mm[key]
        if(value != "+"):
            fout.write(key+"\t" + value + "\n")
        else:
            print(key + "\n")
    return 0

def create_ab(file_in,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    id = ""
    str_a = ""
    str_b = ""
    i_50 = 0;
    inp = 0;
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        vec2 = vec1[2].split(":")
        if(len(id) == 0):
            id = vec1[0]
            str_a = vec2[0] + "N"
            str_b = vec2[1] + "N"
            i_50 = len(str_a) - 1
            fout2.write(line.strip() + "\t" + str(i_50) + "\n")
        elif(id == vec1[0]):
            str_a += vec2[0] + "N"
            str_b += vec2[1] + "N"
            i_50 = len(str_a) - 1
            fout2.write(line.strip() + "\t" + str(i_50) + "\n")
        else:
            print(id,inp);
            inp = 1
            fout.write(">" +id + "_a\n")
            fout.write(str_a+"\n")
            fout1.write(">" +id + "_b\n")
            fout1.write(str_b+"\n")
            id = vec1[0]
            str_a = vec2[0] + "N"
            str_b = vec2[1] + "N"
            i_50 = len(str_a) - 1
            fout2.write(line.strip() + "\t" + str(i_50) + "\n")
    print(id,inp);
    fout.write(">" +id + "_a\n")
    fout.write(str_a+"\n")
    fout1.write(">" +id + "_b\n")
    fout1.write(str_b+"\n")
    i_50 = len(str_a) - 1
    fout2.write(line.strip() + "\t" + str(i_50) + "\n")
    return 0

def create_del_fa(file_name):
    file_out = "data_fa/" + file_name + ".fa"
    file_in = "del/" + file_name
    fout = open(file_out,'w')
    ii = 1
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        fout.write(">" + str(ii) + "_" + vec1[1] + "\n")
        fout.write(vec1[0] + "\n")
    return 0

def summary_dp(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        if(line[0] != "@"):
            vec1 = line.strip().split("\t")
            id_pv = vec1[0].split("_")[1]
            print(vec1[0],id_pv)
            if(mm.has_key(id_pv)):
                mm[id_pv] += 1
            else:
                mm.setdefault(id_pv,1)
    for key in mm:
        fout.write(key+"\t"+mm[key]+"\n")
    return 0

def summary_dp_new(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}

    for line in open(file_in1):
        if(line[0] != "@"):
            vec1 = line.strip().split("\t")
            id = vec1[2]+ ":" + vec1[3]
            id_pv = vec1[0].split("_")[1]
            if(mm.has_key(id)):
                print(line)
            else:
                mm.setdefault(id,id_pv)
                
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        id_a = vec1[0] + "_a:" + vec1[3]
        id_b = vec1[0] + "_b:" + vec1[3]
        i_a = 0
        i_b = 0
        if(mm.has_key(id_a)):
            i_a = int(mm[id_a])
        if(mm.has_key(id_b)):
            i_b = int(mm[id_b]) 
        str_snp = ""
        if(i_a > 0 and i_b > 0):
            str_snp = "AB"
        elif(i_a > 0 and i_b == 0):
            str_snp = "AA"
        else:
            str_snp = "BB"
        str_dp = str(i_a) + ":" + str(i_b)
        total_dp = i_a + i_b
        fout.write(line.strip() + "\t" + str_snp + "\t" + str_dp + "\t" + str(total_dp) + "\n")
    return 0

def fx_tj(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       my_map = Counter(vec1)
       i_a = my_map["AA"];
       i_b = my_map["BB"];
       i_ab = my_map["AB"];
       i_total = i_a + i_b + i_ab
       i_c = 0;
       i_z = 0;
       if(i_ab > 0):
           i_z = 1
       if(i_a > 0):
           i_c += 1
       if(i_b > 0):
           i_c += 1
       fout.write(id + "\t" + str(i_c) + "\t" + str(i_z) + "\t" + str(i_total) + "\t" + str(i_a) + "\t" + str(i_b) + "\t" + str(i_ab) + "\n")
    return 0

def kf_1_2_1(oa,ob,oab):
    total = oa + ob + oab
    da = 0.25
    db = 0.25
    dab = 0.5
    
    ea = da * total
    eb = db * total
    eab = dab * total
    
    oe_a = oa - ea
    oe_b = ob - eb
    oe_ab = oab - eab
    
    double_oe_a = oe_a * oe_a
    double_oe_b = oe_b * oe_b
    double_oe_ab = oe_ab * oe_ab
    
    r_a = double_oe_a / ea
    r_b = double_oe_b / eb
    r_ab = double_oe_ab / eab
    
    total_re = r_a + r_b + r_ab
    return total_re

def kf_2_1(oa,oab):
    total = oa + oab
    da = 0.66
    dab = 1-da
    
    ea = da * total
    eab = dab * total
    
    oe_a = oa - ea
    oe_ab = oab - eab
    
    double_oe_a = oe_a * oe_a
    double_oe_ab = oe_ab * oe_ab
    
    r_a = double_oe_a / ea
    r_ab = double_oe_ab / eab
    
    total_re = r_a + r_ab
    return total_re

def kf_3_1(oa,oab):
    total = oa + oab
    da = 0.75
    dab = 1-da
    
    ea = da * total
    eab = dab * total
    
    oe_a = oa - ea
    oe_ab = oab - eab
    
    double_oe_a = oe_a * oe_a
    double_oe_ab = oe_ab * oe_ab
    
    r_a = double_oe_a / ea
    r_ab = double_oe_ab / eab
    
    total_re = r_a + r_ab
    return total_re


def kf_1_1(oa,ob):
    total = oa + ob
    da = 0.5
    db = 0.5
    
    ea = da * total
    eb = db * total
    
    oe_a = oa - ea
    oe_b = ob - eb
    
    double_oe_a = oe_a * oe_a
    double_oe_b = oe_b * oe_b

    r_a = double_oe_a / ea
    r_b = double_oe_b / eb
    
    total_re = r_a + r_b
    
    return total_re

def kf_33(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       ia = int(vec1[4])
       ib = int(vec1[5])
       iab = int(vec1[6])
       #kf = kf_1_2(ia,ib,iab)
       #fout.write(line.strip() + "\t" + str(kf) + "\n")
    return 0

def kf_33_1_2(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       ia = int(vec1[4])
       ib = int(vec1[5])
       iab = int(vec1[6])
       o_in = 0
       if(ia > 0):
           o_in = ia
       elif(ib > 0):
           o_in = ib
       else:
           print(line);
       kf = kf_1_2(o_in,iab)
       fout.write(line.strip() + "\t" + str(kf) + "\n")
    return 0

def create_mhd_data(file_in,file_out):
    fout = open(file_out,'w')
    ii = 1
    ic = 1
    old_group = ""
    fout.write("SNP\tCHR\tBP\tP\n")
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       group = vec1[0].split("_")[1]
       kf = vec1[len(vec1) - 1]
       snpi = "rs" + str(ii)
       ii += 1
       if(old_group == ""):
           old_group = group
       elif(old_group == group):
           ic += 1
       else:
           ic = 1
           old_group = group
       fout.write(snpi + "\t" + old_group[5:]  + "\t" + str(ic) + "\t" + kf + "\n")
    return 0

def return_site(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm={}
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       mm.setdefault(vec1[0],1)
    for line in open(file_in1):
       vec1 = line.strip().split("\t")
       id =  vec1[0] + ":" + vec1[3]
       if(mm.has_key(id)):
           fout.write(id + "\t" + vec1[1] + "\t" + vec1[2] + "\n")
    return 0

def qj_tj(file_in,file_out,ioio):
    fout = open(file_out,'w')
    mm={}
    mm1={}
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       site = vec1[1]
       id = vec1[0].split(":")[0]
       #print(id.split("_"))
       I_length = id.split("_")[7]
       mm1.setdefault(id,I_length)
       if(mm.has_key(id)):
           mm[id] += site + ":"
       else:
           mm.setdefault(id,site + ":")
    fout.write("SNP\tCHR\tBP\tP\n")
    i_cc = 1
    for key in mm1:
        print(key+"\n")
        i_len = int(mm1[key])
        mm_temp = {}
        for i_site in range(0,i_len/ioio + 1):
            mm_temp.setdefault(i_site,0)
        i_data = mm.get(key)[0:-1]
        #print(i_data)
        vec1 = i_data.split(":")
        for kk in vec1:
            i_kk = int(kk) / ioio
            if(mm_temp.has_key(i_kk)):
                mm_temp[i_kk] += 1
            else:
                print(kk,i_kk)
        o_id = key.split("_")[1][5:]
        for kk in mm_temp:
            #fout.write(str(key)+"\t"+str(kk)+"\t"+str(mm_temp[kk]) + "\n")
            s_1 = "rs" + str(i_cc)
            i_cc += 1
            fout.write(s_1 + "\t" + o_id  + "\t" + str(kk) + "\t" + str(mm_temp[kk]) + "\n")
            
    return 0

def qj_tj_2(file_in,file_out,ioio):
    fout = open(file_out,'w')
    mm={}
    mm1={}
    mm2 = {}
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       site = vec1[1]
       id = vec1[0].split(":")[0]
       #print(id.split("_"))
       I_length = id.split("_")[7]
       mm1.setdefault(id,I_length)
       mm2.setdefault(id,0)
       if(mm.has_key(id)):
           mm[id] += site + ":"
       else:
           mm.setdefault(id,site + ":")
    
    for key in mm1:
        #print(key+"\n")
        i_data = mm.get(key)[0:-1]
        vec1 = i_data.split(":")
        
        current_ss = 0
        current_site = int(vec1[current_ss])
        next_ss = current_ss + 1
        begin_site = current_site
        end_site = current_site
        link_site = str(begin_site) + ":"
        while(next_ss <= len(vec1) - 1):
            #print(next_ss,len(vec1))
            size_aa = str(len(link_site.split(":"))- 1)
            if(next_ss == 34275448):
                pass
            next_site = int(vec1[next_ss])
            diff = next_site - current_site
            if(diff <= ioio):
                end_site = next_site
                link_site += str(next_site) + ":"
            else:
                if(int(size_aa) > 2):
                    mm2[key] += 1
                fout.write(key + "\t" +link_site + "\t" + size_aa+ "\n")
                begin_site = next_site
                link_site = str(begin_site) + ":"
            current_site = next_site
            next_ss += 1

        size_aa = str(len(link_site.split(":"))- 1)
        if(int(size_aa) > 2):
            mm2[key] += 1
        fout.write(key + "\t" +link_site + "\t" + size_aa + "\n")
    i_scaff = 0
    i_qj = 0
    
    for key in mm2:
        is_x = False
        if(mm2[key] > 0):
           #print(key,mm2[key])
           i_qj += mm2[key]
           is_x = True
        if(is_x == True):
            i_scaff += 1
            
    print(ioio,i_scaff,i_qj)
           
    return 0

def huanyuan_1(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm={}
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       site = vec1[0]
       mm.setdefault(site,"")
    
    for line in open(file_in1):
       vec1 = line.strip().split("\t")
       site = vec1[0] + ":" + vec1[3]
       if(mm.has_key(site)):
           mm[site] = vec1[5]
           
    for key in mm:
        fout.write(key+"\t" + mm[key] + "\n")
    return 0

def count_tj(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       s_out = id + "\t"
       i_a = 0
       i_b = 0
       for ii in range(1,len(vec1)):
           tt = vec1[ii]
           if(tt != id):
               vec2 = tt.split(":")
               i_a += int(vec2[0])
               i_b += int(vec2[1])
               s_out += tt + "\t"
       if(id == "Lachesis_group17__22_contigs__length_31834974:3714076"):
           pass
       kf11 = kf_1_1(i_a,i_b)
       kf12 = 0
       ip = 0
       if(i_a <= i_b):
           ip = i_b
           i_b = i_a
           i_a = ip
       kf12 = kf_2_1(i_a, i_b)
       fout.write(s_out+ str(i_a) + "_" + str(i_b) + "\t" + str(kf11) + "\t" + str(kf12)  + "\n")
    return 0

def count_tj_all(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       s_out = id + "\t"
       i_0_0  = 0
       i_0_1 = 0 
       for ii in range(1,len(vec1)):
           tt = vec1[ii]
           if(tt != id):
               s_out += tt + "\t"
               if(tt == "N"):
                   i_0_0 += 1
               else:
                   i_0_1 += 1
       fout.write(s_out+ "\t" + str(i_0_1)  + "\n")
    return 0

def w_2_2line(r11,r12,r21,r22):
    r1_all = float(r11)+r12
    r2_all = float(r21)+r22
    rall = r1_all + r2_all
    r_1 = float(r11) + r21
    r_2 = float(r12) + r22
    e11 = r1_all * r_1 /rall
    e12 = r1_all * r_2 /rall
    e21 = r2_all * r_1 /rall
    e22 = r2_all * r_2 /rall
    k11 = (r11 - e11) * (r11 - e11) /e11
    k12 = (r12 - e12) * (r12 - e12) /e12
    k21 = (r21 - e21) * (r21 - e21) /e21
    k22 = (r22 - e22) * (r22 - e22) /e22
    kk = round(k11 + k12 + k21 + k22,3)
    return kk

def llb_83(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       a_33 = int(vec1[34])
       b_33 = int(vec1[35])
       all_33 = float(a_33 + b_33)
       a_33_d = a_33 / all_33 * 100
       b_33_d = b_33 / all_33 * 100
       
       a_50 = int(vec1[88])
       b_50 = int(vec1[89])
       all_50 = float(a_50 + b_50)
       
       a_50_d = a_50 / all_33 * 100
       b_50_d = b_50 / all_33 * 100
       if(all_50 >= 100 and all_33 >= 100):
           kf = w_2_2line(a_33_d, b_33_d, a_50_d, b_50_d)
           fout.write(id + "\t" + vec1[34] + "\t" + vec1[35] + "\t" + vec1[88] + "\t" + vec1[89] + "\t" + str(kf) + "\n")
       
    return 0

def tiqu_id(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
       vec1 = line.strip()
       mm.setdefault(vec1)
    for line in open(file_in1):
       vec1 = line.strip().split("\t")
       if(mm.has_key(vec1[0])):
           fout.write(line)
    return 0

def llb_123(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    i_cc = 0
    kk = 0
    oo = ""
    fout1.write("SNP\tCHR\tBP\tP\n")
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       a_33 = int(vec1[1])
       b_33 = int(vec1[2])
       all_33 = float(a_33 + b_33)
       p1 = a_33 / all_33 * 100
       p2 = b_33 / all_33 * 100
       a_33_d = float('%.2f' %p1)
       b_33_d = float('%.2f' %p2)
       
       a_50 = int(vec1[3])
       b_50 = int(vec1[4])
       all_50 = float(a_50 + b_50)
       p3 = a_50 / all_33 * 100
       p4 = b_50 / all_33 * 100
       a_50_d = float('%.2f' %p3)
       b_50_d = float('%.2f' %p4)
       if(all_50 >= 100 and all_33 >= 100):
           kf = w_2_2line(a_33_d, b_33_d, a_50_d, b_50_d)
           o_id = vec1[0].split("_")[1][5:]
           fout.write(line.strip() + "\t" + str(kf) + "\n")
           s_1 = "rs" + str(i_cc)
           i_cc += 1
           if(o_id != "0" and o_id != oo):
               kk = 0
               oo = o_id
           else:
               oo = o_id
           kk += 1
           fout1.write(s_1 + "\t" + o_id  + "\t" + str(kk) + "\t" + str(kf) + "\n")
    return 0

def llb_33_11(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    i_cc = 0
    kk = 0
    oo = ""
    fout1.write("SNP\tCHR\tBP\tP\n")
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       a_33 = int(vec1[1])
       b_33 = int(vec1[2])
       all_33 = float(a_33 + b_33)
       p1 = a_33 / all_33 * 100
       p2 = b_33 / all_33 * 100
       a_33_d = float('%.2f' %p1)
       b_33_d = float('%.2f' %p2)
       
       kf= 1000
       
       if(a_33_d >= b_33_d):
           kf = kf_2_1(a_33_d, b_33_d)
       else:
           kf = kf_2_1(b_33_d, a_33_d)

       a_50 = int(vec1[3])
       b_50 = int(vec1[4])
       all_50 = float(a_50 + b_50)
       p3 = a_50 / all_33 * 100
       p4 = b_50 / all_33 * 100
       a_50_d = float('%.2f' %p3)
       b_50_d = float('%.2f' %p4)
       kf11 = kf_1_1(a_50_d, b_50_d)
       
       o_id = vec1[0].split("_")[1][5:]
       fout.write(line.strip() + "\t" + str(kf) + "\t" + str(kf11) + "\n")
       s_1 = "rs" + str(i_cc)
       i_cc += 1
       if(o_id != "0" and o_id != oo):
           kk = 0
           oo = o_id
       else:
           oo = o_id
       kk += 1
       fout1.write(s_1 + "\t" + o_id  + "\t" + str(kk) + "\t" + str(kf) + "\n")
    return 0


def site_329(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0].split(":")[0]
       i_site = int(vec1[1]) - 150
       s1 = "echo '>" + vec1[0] + "' >> 1w.temp"
       ss = "grep -A1 "+id+" v12.fa|awk '{print(substr($0,"+str(i_site)+",301));}' >> 1w.temp"
       fout.write(s1 +"\n")
       fout.write(ss + "\n")
    return 0

def test_llb(bb):
    a1 = 30 * bb
    a2 = 70 * bb
    b1 = 50 * bb
    b2 = 50 * bb
    
    
    aa = float(a1 + a2)
    p1 = a1 / aa
    p2 = a2 / aa
  
    bb = float(b1 + b2)
    i1 = b1 / bb
    i2 = b2 / bb
         
    kf1 = w_2_2line(a1,a2,b1,b2)
    kf2 = kf_1_1(a1,a2)
    
    kfp1 = w_2_2line(p1,p2,i1,i2)
    kfi2 = kf_1_1(p1,p2)
    
    print(kf1,kf2,kfp1,kfi2)

def kf_1_1_2_1(file_in,file_out):
    fout = open(file_out,'w')
    i_cc = 0
    kk = 0
    oo = ""
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       a_33 = int(vec1[1])
       b_33 = int(vec1[2])
       all_33 = float(a_33 + b_33)
       p1 = a_33 / all_33 * 100
       p2 = b_33 / all_33 * 100
       a_33_d = float('%.2f' %p1)
       b_33_d = float('%.2f' %p2)
       
       a_50 = int(vec1[3])
       b_50 = int(vec1[4])
       all_50 = float(a_50 + b_50)
       p3 = a_50 / all_33 * 100
       p4 = b_50 / all_33 * 100
       a_50_d = float('%.2f' %p3)
       b_50_d = float('%.2f' %p4)
       
       kf = 0
       kf1 = kf_1_1(a_50_d, b_50_d)
       if(a_33_d >= b_33_d):
           kf = kf_2_1(a_33_d, b_33_d)
       else:
           kf = kf_2_1(b_33_d, a_33_d)
           
       fout.write(line.strip() + "\t" + str(kf) + "\t" + str(kf1) + "\n")
            
    return 0

def kf_3_1_1_1(file_in,file_out):
    fout = open(file_out,'w')
    i_cc = 0
    kk = 0
    oo = ""
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       a_33 = int(vec1[1])
       b_33 = int(vec1[2])
       all_33 = float(a_33 + b_33)
       p1 = a_33 / all_33 * 100
       p2 = b_33 / all_33 * 100
       a_33_d = float('%.2f' %p1)
       b_33_d = float('%.2f' %p2)
       
       a_50 = int(vec1[3])
       b_50 = int(vec1[4])
       all_50 = float(a_50 + b_50)
       p3 = a_50 / all_33 * 100
       p4 = b_50 / all_33 * 100
       a_50_d = float('%.2f' %p3)
       b_50_d = float('%.2f' %p4)
       
       kf = kf_1_1(a_33_d, b_33_d)
       kf1 = 0
       if(a_50_d >= b_50_d):
           kf1 = kf_3_1(a_50_d, b_50_d)
       else:
           kf1 = kf_3_1(b_50_d, a_50_d)
           
       fout.write(line.strip() + "\t" + str(kf) + "\t" + str(kf1) + "\n")
            
    return 0


def create_fa_50(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       reads = vec1[2].split(":")[0]
       for i in range(1,50):
           fout.write(">" + id + "_" + str(i) + "\n")
           fout.write(reads)
    return 0

def gene_3490(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       gene = vec1[2]
       if(mm.has_key(gene)):
           mm[gene] += id + "@"
       else:
           mm.setdefault(gene,id + "@")
           
    for key in mm:
        vec = mm[key][0:-1].split("@")
        ss = ""
        for kk in vec:
            ss += kk + "\t"
        fout.write(key + "\t" + str(len(vec)) + "\t" + ss + "\n")
    return 0


def tj_30_cx(file_in,fle_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       cc = vec1[1]
       mm.setdefault(id,cc)
       
    ir = 0
    ip = 0
    mm1 = {}
    istr = ""
    for line in open(fle_in1):
        vec1 = line.strip().split("\t")
        id = vec1[0]
        ii = mm[id]
        
        for ic in range(1,len(vec1)):
            ww = vec1[ic]
            istr = ""
                
            if(ww != id):
                vec2 = ww.split(":")
                ia = int(vec2[0])
                ib = int(vec2[1])
                if(ii == "B"):
                    if(ib > 0):
                        istr = "1"
                    else:
                        istr = "0"
                else:
                    if(ia >0):
                        istr = "1"
                    else:
                        istr = "0"
                if(mm1.has_key(ic)):
                   mm1[ic] += istr
                else:
                    mm1.setdefault(ic,istr) 
    for key in mm1:
        uu = mm1[key]
        k1 = 0
        for k in range(0,len(uu)):
            if(uu[k] == "1"):
                k1 += 1
        fout.write(str(key)+"\t"+uu+"\t"+str(k1) + "\n")
    return 0

def tiqu_3490(file_in,fle_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       mm.setdefault(id,1)
       
    for line in open(fle_in1):
        vec1 = line.strip().split("\t")
        id = vec1[0]
        if(mm.has_key(id)):
            fout.write(line.strip()+"\n")
    return 0

def fx_3490(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       str_fx = id + "\t"
       for key in range(1,len(vec1)):
           fx =vec1[key]
           if(fx != id):
                vec2 = fx.split(":")
                if(len(vec2) > 1):
                    ia = int(vec2[0])
                    ib = int(vec2[1])
                    
                    if(ia + ib >= 9):
                        if(ia > 0 and ib > 0):
                            str_fx += "AB\t"
                        elif(ia > 0 and ib == 0):
                            str_fx += "AA\t"
                        else:
                            str_fx += "BB\t"
                    else:
                        str_fx += fx+"\t"
                else:
                    str_fx += "=\t"
       fout.write(str_fx + "\n")
    return 0

def zhuanhuan_3490(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}

    for line in open(file_in):
       vec1 = line.strip().split("\t")
       id = vec1[0]
       for key in range(0,len(vec1)):
           fx =vec1[key]
           if(mm.has_key(key)):
               mm[key] += fx +"\t"
           else:
               mm.setdefault(key,fx + "\t")
    for key in mm:
        vec1 = mm[key].split("\t")
        my_map = Counter(vec1)
        i_a = my_map["AA"];
        i_b = my_map["BB"];
        i_ab = my_map["AB"];
        i_total = i_a + i_b + i_ab
        str_00 = str(i_a) + "\t" + str(i_b) + "\t" + str(i_ab) + "\t" + str(i_total)
        fout.write(str(key) + "\t" + mm[key] + "\t" + str_00 + "\n")
    return 0

def tiqu_5(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm_scaff = {}
    mm_scaff.setdefault("NC_007112.7","")
    mm_scaff.setdefault("NC_007113.7","")
    mm_scaff.setdefault("NC_007114.7","")
    mm_scaff.setdefault("NC_007115.7","")
    mm_scaff.setdefault("NC_007116.7","")
    mm_scaff.setdefault("NC_007117.7","")
    mm_scaff.setdefault("NC_007118.7","")
    mm_scaff.setdefault("NC_007119.7","")
    mm_scaff.setdefault("NC_007120.7","")
    mm_scaff.setdefault("NC_007121.7","")
    mm_scaff.setdefault("NC_007122.7","")
    mm_scaff.setdefault("NC_007123.7","")
    mm_scaff.setdefault("NC_007124.7","")
    mm_scaff.setdefault("NC_007125.7","")
    mm_scaff.setdefault("NC_007126.7","")
    mm_scaff.setdefault("NC_007127.7","")
    mm_scaff.setdefault("NC_007128.7","")
    mm_scaff.setdefault("NC_007129.7","")
    mm_scaff.setdefault("NC_007130.7","")
    mm_scaff.setdefault("NC_007131.7","")
    mm_scaff.setdefault("NC_007132.7","")
    mm_scaff.setdefault("NC_007133.7","")
    mm_scaff.setdefault("NC_007134.7","")
    mm_scaff.setdefault("NC_007135.7","")
    mm_scaff.setdefault("NC_007136.7","")
    for line in open(file_in):
       if(line[0:1] != "@"):
           vec1 = line.strip().split("\t")
           id = vec1[0]
           scaff = vec1[2]
           site = vec1[3]
           mm_id = scaff + ":" + site
           if(mm.has_key(mm_id)):
               mm[mm_id] += "," + id
           else:
               mm.setdefault(mm_id,id)
    for key in mm:
        scaff = key.split(":")[0]
        if(mm_scaff.has_key(scaff)):
            vec1 = mm[key].split(",")
            fout.write(key+"\t"+str(len(vec1))+"\t"+mm[key]+"\n")
    return 0

def is_map(vec1):
    mm = {}
    for key_temp in range(0,len(vec1) - 1):
        key = vec1[key_temp]
        id = key.split("_")[0]
        if(mm.has_key(id)):
            mm[id] += 1
        else:
            mm.setdefault(id,1)
    ss = str(len(mm)) + "\t"
    for key in mm:
        ss += key + "\t" + str(mm[key]) + "\t"
    return ss

def del_re(file_in,file_out):
    fout = open(file_out,'w')
    site = ""
    ss = ""
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        if(vec1[0] != site):
            if(len(site) > 2):
                vec3 = ss.split(":")
                #sw = is_map(vec1)
                fout.write(site +"\t"  + str(len(vec3)  -1) + "\t" + ss  + "\n")
            site = vec1[0]
            ss = vec1[1] + ":"
        else:
            ss += vec1[1] + ":"
                
    return 0

def is_pair(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        vec2 = vec1[2].split(":")
        mm = {}
        for ii in range(0,len(vec2) - 1):
            ss = vec2[ii][0:-3]
            if(mm.has_key(ss)):
                mm[ss] += 1
            else:
                mm.setdefault(ss,1)
        if(len(mm) == 2):
            it = 0
            sid = ""
            for key1 in mm:
                if(mm[key1] == 5):
                    it += 1
                    sid += key1 + "\n"
            if(it == 2):
                fout.write(line.strip() + "\n")
                print(sid)
    return 0

def is_2(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        if(mm.has_key(vec1[0])):
            mm[vec1[0]] += vec1[1] + ":"
        else:
            mm.setdefault(vec1[0],vec1[1] + ":")
    for key in mm:
        reads_all = mm[key].split(":")
        fout.write(key+"\t"+str(len(reads_all)) + "\t" + mm[key] + "\n")
    return 0

def find_pair_reads(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        vec2 = vec1[2].split(":")
        mm_temp = {}
        for ii in range(0,len(vec2) - 1):
            sw = vec2[ii][0:-2]
            mm_temp.setdefault(sw,1)
        for key in mm_temp:
            mm.setdefault(key,vec1[0])
    ii = 1
    id = ""
    for line in open(file_in):
        reads = line.strip()
        if(ii == 1):
            id = reads[1:-1]
        else:
            ii = 0
            if(mm.has_key(id)):
                fout.write(id+"\t"+reads+"\t"+mm[key]+"\n")
        ii += 1
    return 0

def line_reads(file_in,file_outa,file_outb):
    fouta = open(file_outa,'w')
    foutb = open(file_outb,'w')
    scaff = ""
    is_1 = 1
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        reads = vec1[1]
        s_temp = vec1[2].split(":")[0]
        if(is_1 == 1):
            if(scaff != s_temp):
                if(len(scaff) < 2):
                    fouta.write(">"+s_temp+"\n")
                else:
                    fouta.write("\n>"+s_temp+"\n")
            fouta.write(reads+"N")
        else:
            if(scaff != s_temp):
                if(len(scaff) < 2):
                    foutb.write(">"+s_temp+"\n")
                else:
                    foutb.write("\n>"+s_temp+"\n")
                scaff = s_temp
            foutb.write(reads+"N")
            is_1 = 0
        is_1 += 1
    return 0

def create_fa(file_in):
    file_out = file_in + ".fa"
    fouta = open(file_out,'w')
    i_line = 1
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        fouta.write(">"+str(i_line) + "_"+vec1[1]+"\n")
        fouta.write(vec1[0]+"\n")
    return 0
        
    
def find_all(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    scaff = ""
    is_1 = 1
    mm = {}
    for line in open(file_in):
        line = line.strip()
        if(is_1 == 1):
            scaff = line[1:]
        else:
            is_1 = 0
            mm.setdefault(scaff,line)
        is_1 += 1
    print("over\t"+file_in1);
    id = ""
    is_1 = 1
    for line in open(file_in1):
        line = line.strip()
        if(is_1 == 1):
            id = line[1:]
        else:
            is_1 = 0
            for key in mm:
                reads = mm[key]
                site = reads.find(line)
                if(site >= 0):
                    fout.write(id + "\t" + key + "\t" + str(site) + "\n")
        is_1 += 1

def reads_tj(file_name):
    file_ina = "ref_a/" + file_name + ".sam"
    file_inb = "ref_b/" + file_name + ".sam"
    file_out = file_name + ".tj"
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_ina):
        vec1 = line.strip().split("\t")
        if(len(vec1) > 10):
            id = vec1[0] + "_A:"
            scaff = vec1[2] + ":" + vec1[3]
            mm.setdefault(scaff,id)
    for line in open(file_inb):
        vec1 = line.strip().split("\t")
        if(len(vec1) > 10):
            id = vec1[0] + "_B:"
            scaff = vec1[2] + ":" + vec1[3]
            if(mm.has_key(scaff)):
                mm[scaff] += id
            else:
                mm.setdefault(scaff,id)
    for key in mm:
        kk = mm[key][0:-1].split(":")
        i_total = 0
        i_a= 0
        i_b= 0
        for k1 in range(0,len(kk)-1):
            vv= kk[k1].split("_")
            a_type = vv[2]
            if(a_type == "A"):
                i_a = int(vv[1])
            else:
                i_b = int(vv[1])
        i_total = i_a + i_b
        fout.write(key+"\t"+mm[key] + "\t"+str(len(kk)) + "\t"+str(i_a) + ":" + str(i_b) + "\t"+str(i_total) + "\n")
    return 0

def reads_he(file_site,file_total,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_site):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0],"")
    for file_name in open(file_total):
        file_name = file_name.strip()
        print(file_name)
        i_file = file_name + ".tj"
        for line in open(i_file):
            vec1 = line.strip().split("\t")
            id = vec1[0]
            a_b = vec1[3]
            total = vec1[4]
            ikk = file_name + "," + a_b +"," + total + ";"
            mm[id] += ikk
    for key in mm:
        vec1 = mm[key][0:-1].split(";")
        ff = ""
        a_b = ""
        i_total = ""
        for kk in range(0,len(vec1)):
            vec2 = vec1[kk].split(",")
            ff += vec2[0] + ","
            a_b += vec2[1] + ","
            i_total += vec2[2] + ","
        fout.write(key + "\t" + i_total + "\t" + a_b + "\t" + ff + "\n")
    return 0

def del_reaa(file_site,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_site):
        line = line.strip()
        if(mm.has_key(line)):
            mm[line] += 1
        else:
            mm.setdefault(line,1)
    for key in mm:
        fout.write(key+"\t"+str(mm[key]))
    return 0

def reads_mean(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[1][0:-1].split(",")
        ll = []
        
        for key in range(0,len(vec1)):
            kk = int(vec1[key])
            ll.append(kk)
        i_tt = len(ll)
        while(len(ll) < 52):
            ll.append(0)
            
        i_mean = np.mean(ll)
        i_var = np.var(ll,ddof=1)
        ii = round(i_mean,1)
        if(mm.has_key(ii)):
            mm[ii] += 1
        else:
            mm.setdefault(ii,1)
        #fout.write(vec[0]+"\t"+str(i_mean)+"\t"+str(math.sqrt(i_var))+"\n")
        fout.write(vec[0]+"\t"+str(i_mean)+"\t"+str(math.sqrt(i_var))+"\t"+str(i_tt)+"\t"+vec[1]+"\t"+vec[2]+"\n")
    for key in mm:
        print(str(key)+"\t"+str(mm[key]));
    return 0

def reads_sd(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        ss = np.sqrt(float(vec[2]))
        ii = round(ss,1)
        if(mm.has_key(ii)):
            mm[ii] += 1
        else:
            mm.setdefault(ii,1)
    for key in mm:
        print(str(key)+"\t"+str(mm[key]));
    return 0

def tiqu_ref_50(file_site,file_out):
    fout = open(file_out,'w')
    id = ""
    for line in open(file_site):
        line = line.strip()
        if(len(line) < 100):
            id = line[1:]
        else:
            i_site= 1
            while(i_site < len(line)):
                ss = line[i_site-1:i_site-1 +50]
                fout.write(id+":"+str(i_site)+"\t"+ss + "\n")
                i_site += 51
    return 0

def change_mhd(file_site,file_out):
    fout = open(file_out,'w')
    fout.write("SNP\tCHR\tBP\tP\m")
    ii = 0
    it = 0
    scaff = ""
    for line in open(file_site):
        vec = line.strip().split("\t")
        if(scaff == ""):
            scaff = vec[0]
        elif(scaff != vec[0]):
            scaff = vec[0]
            it  = 0
        fout.write("rs" + str(ii) + "\t" + scaff[5:-2] + "\t" + str(it) + "\t" + vec[6] + "\n")
        ii += 1
        it += 1   
    return 0

def gene_type(file_site,file_out):
    fout = open(file_out,'w')
    for line in open(file_site):
        vec = line.strip().split("\t")
        vec1 = vec[3][0:-1].split(",")
        vec2 = vec[2][0:-1].split(",")
        aa_1 = 0
        ab_1 = 0
        bb_1 = 0
        aa_2 = 0
        ab_2 = 0
        bb_2 = 0
        
        for ii in range(0,len(vec1)):
            gt = vec1[ii][0]
            a = int(vec2[ii].split(":")[0])
            b = int(vec2[ii].split(":")[1])
            aa = 0
            bb = 0
            ab = 0
            if(a > 0 and b == 0):
                aa = 1
            elif(a == 0 and b > 0):
                bb = 1
            elif(a > 0 and b > 0):
                ab = 1
            if(gt == "c"):
                aa_1 += aa
                ab_1 += ab
                bb_1 += bb
            else:
                aa_2 += aa
                ab_2 += ab
                bb_2 += bb                
        fout.write(line.strip() + "\t" + str(aa_1) + "\t" + str(bb_1) + "\t" + str(ab_1) + "\t"  + str(aa_2) + "\t" + str(bb_2) + "\t" + str(ab_2) + "\n")
    
    return 0

def shai_gene_type(file_site,file_out):
    fout = open(file_out,'w')
    for line in open(file_site):
        vec = line.strip().split("\t")
        a_1 = int(vec[4])
        b_1 = int(vec[5])
        ab_1 = int(vec[6])
        a_2 = int(vec[7])
        b_2 = int(vec[8])
        ab_2 = int(vec[9])
        type1 = 0
        type2 = 0
        if(a_1 > 0):
            type1 += 1
        if(b_1 > 0):
            type1 += 1
        if(ab_1 > 0):
            type1 += 1
        if(a_2 > 0):
            type2 += 1
        if(b_2 > 0):
            type2 += 1
        if(ab_2 > 0):
            type2 += 1
        iswrite = 0
        if(type1 == 2 and type2 == 2):
            if(a_1 > 0 and ab_1 > 0 and b_2 > 0 and ab_2 > 0):
                iswrite += 1
            elif(b_1 > 0 and ab_1 > 0 and a_2 > 0 and ab_2 > 0):
                iswrite += 1
        if(iswrite > 0):
            fout.write(line.strip() + "\n")
    return 0

def chang_gtf(file_site,file_out):
    fout = open(file_out,'w')
    for line in open(file_site):
        vec = line.strip().split("\t")
        vec1 = vec[8].split(";")
        ss = vec1[1] + ";" + vec1[0] + ";" + vec1[2]
        fout.write(vec[0] + "\t" + vec[1] + "\t" + vec[2] + "\t" + vec[3] + "\t" + vec[4] + "\t" + vec[5] + "\t" + vec[6] + "\t" + vec[7] + "\t"  + ss + "\n")
    return 0

def del_re_kk(file_site,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_site):
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[1]
        if(mm.has_key(id)):
            mm[id] += site + ":"
        else:
            mm.setdefault(id,site + ":")
    for key in mm:
        vec = mm[key].split(":")
        fout.write(key + "\t" + str(len(vec)) + "\t" + mm[key] + "\n")
    return 0


def shai_21_type(file_site,file_out):
    fout = open(file_out,'w')
    for line in open(file_site):
        vec = line.strip().split("\t")
        a_1 = int(vec[4]) + int(vec[5])
        ab_1 = int(vec[6])
        b_2 = int(vec[7]) + int(vec[8])
        ab_2 = int(vec[9])
        kf11 = kf_1_1(a_1, ab_1)
        kf13 = kf_1_1(b_2, ab_2)
        fout.write(line.strip() + "\t" + str(kf11) + "\t" + str(kf13) + "\n")
    return 0

def find_pair_reads11(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0],vec1[1])
    print("load over");
    for line in open(file_in1):
        vec1 = line.strip().split(":")
        id = vec1[0]
        if(mm.has_key(id)):
            fout.write(id+"\t"+vec1[1]+"\t"+mm[id]+"\n")
    return 0

def get_a_reads(file_in,file_out):
    fout = open(file_out,'w')
    ii = 1
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        if(ii == 2):
            fout.write(">"+vec1[0]+"\n"+vec1[1] + "\n")
        else:
            ii = 0
        ii += 1
    return 0

def pute_dd(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    id_old = ""
    site_old = 0
    for line in open(file_in):
        vec1 = line.strip().split(":")
        id = vec1[0]
        site = int(vec1[1])
        if(id_old != id):
            id_old = id
            site_old = site
        else:
            dd = site - site_old
            it = dd / 500
            fout.write(id + "\t" + str(site_old) + "\t" + str(site) + "\t" + str(it) + "\n")
            if(mm.has_key(it)):
                mm[it] += 1
            else:
                mm.setdefault(it,1)
            site_old = site
    for key in mm:
        fout.write(str(key) + "\t" + str(mm[key]) + "\n")
    return 0

def zebre_in_gene(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        ii = int(vec1[1])
        scaff = vec1[0]
        mm.setdefault(ii,scaff)
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        i_begin = int(vec1[3])
        i_end = int(vec1[4])
        ssw = ""
        for key in mm:
            scaff = mm[key]
            if(scaff == vec1[0] and i_begin <= key and key <= i_end):
                ssw+=str(key)+":"
        if(len(ssw)> 2):
            fout.write(ssw + "\t" + line.strip() + "\n")
    return 0

def get_message(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        id = vec1[0]
        ss = vec1[1] + "\t" + vec1[2]
        mm.setdefault(id,ss)
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        id = vec1[1]
        if(mm.has_key(id)):
            fout.write(line.strip() + "\t" + mm[id] + "\n")
    return 0

def tj_site(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split(":")
        id = vec1[0]
        site = vec1[1] + ":"
        if(mm.has_key(id)):
            mm[id] += site
        else:
            mm.setdefault(id,site)
    print(len(mm))
    for key in mm:
        vec1 = mm[key][0:-1].split(":")
        site_old = int(vec1[0])
        mm_k = {}
        mm_k.setdefault(0,"")
        mm_k.setdefault(1,"")
        mm_k.setdefault(2,"")
        mm_k.setdefault(3,"")
        mm_k.setdefault(4,"")
        mm_k.setdefault(5,"")
        mm_k.setdefault(6,"")
        mm_k.setdefault(7,"")
        mm_k.setdefault(8,"")
        mm_k.setdefault(9,"")
        mm_k.setdefault(10,"")
        mm_k.setdefault(-1,"")
        for ii in range(1,len(vec1)):
            site_current = int(vec1[ii])
            cha = site_current - site_old
            site_old = site_current
            ik = cha / 1000
            if(mm_k.has_key(ik)):
                mm_k[ik] += str(cha) + ":"
            else:
                mm_k[-1] += str(cha) + ":"
        s_out = key + "\t"
        aa = 0
        bb = 0
        vv = key.split("_")
        size = int(vv[7])
        cc = float(0)
        for kk in mm_k:
            vec2 = mm_k[kk][0:-1].split(":")
            if(kk == -1):
                print(key,vec2)
            else:
                i_total = 0
                for gg in vec2:
                    if(len(gg) > 0):
                        i_total += int(gg)
                aa += i_total
                bb += len(vec2)
                cc = float(aa / size)
                #s_out += str(aa) + "\t" + str(bb)  + "\t"
                s_out += str(i_total) + "\t" + str(len(vec2))  + "\t"
        fout.write(s_out + "\n")
            
    return 0


def get_ab(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        ss = line.strip()
        mm.setdefault(ss,"")
        mm1.setdefault(ss,"")
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        id = vec1[2] + ":" + vec1[3]
        if(mm.has_key(id)):
            mm[id] += vec1[1] + ":"
            mm1[id] += vec1[0] + ":"
    for key in mm:
        fout.write(key + "\t" + mm[key] + "\t" + mm1[key] + "\n")
    return 0

'''
void find_reads_ab(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        ss = line.strip()
        mm.setdefault(ss,"")
        mm1.setdefault(ss,"")
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        id = vec1[2] + ":" + vec1[3]
        if(mm.has_key(id)):
            mm[id] += vec1[1] + ":"
            mm1[id] += vec1[0] + ":"
    for key in mm:
        fout.write(key + "\t" + mm[key] + "\t" + mm1[key] + "\n")
    return 0
'''

def find_c_gene(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        ss = line.strip()
        mm.setdefault(ss,1)
    for line in open(file_in1):
        ss = line.strip()
        if(not mm.has_key(ss)):
            fout.write(line)
    
def aa(file_in):
    old = 0
    for line in open(file_in):
        vec = line.strip().split(":")
        cu = int(vec[1])
        if(old == 0):
            old = cu
        else:
            cc = cu -old
            if(cc >= 15000):
                print(line)
            old = cu
            
def del_re_bb(file_in,file_in1,out_all,out_1,out_2):
    out_all = open(out_all,'w')
    out_1 = open(out_1,'w')
    out_2 = open(out_2,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0],vec1[1])
    print("load over")
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        id = vec1[0]
        if(mm.has_key(id)):
            out_all.write(id + "\t" + mm[id] + "\t" + vec1[1] + "\n")
            mm.pop(id)
        else:
            out_2.write(id + "\t" + vec1[1] + "\n")
    for key in mm:
        out_1.write(key + "\t" + mm[key] + "\n")
    return 0

def gff_mrna_static(file_in,file_out,file_js):
    fout = open(file_out,'w')
    foutjs = open(file_js,'w')
    mm = {}
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        if(len(vec1) >= 8):
            type = vec1[2]
            vec2 = vec1[8].split(";")
            if(type == "gene"):
               id = vec2[0][3:]
               mm.setdefault(id,"")
            elif(type == "mRNA"):
                id = vec2[0][3:]
                pid = vec2[1][7:]
                mm[pid] += id + ";"
    ii = 0
    mmk = {}
    for key in mm:
        vv = mm[key]
        vec1 = vv.split(";")
        ll = len(vec1) - 1
        if(ll > 0):
            ii += 1
            ik = vec1[0][4:]
            vec2 = ik.split("-")
            mmk.setdefault(vec2[0],1)
        fout.write(key + "\t" + str(ll) + "\t" + vv + "\n")
        
    print(ii)
    for key in mmk:
        ss = "grep -A1 '" + key + "' rna.fa"
        foutjs.write(ss + "\n")
    
def kk(file_in,out_2):
    out_2 = open(out_2,'w')
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        ik = int(vec1[1])
        if(ik > 0):
            vec2 = vec1[2].split(";")
            ik = vec2[0][4:]
            vec3 = ik.split("-")
            ss = "grep -A1 '" + vec3[0] + "' rna.fa"
            out_2.write(ss + "\n")
    return 0

def ddd(file_in,file_in1,file_1,file_2,file_3):
    fout1 = open(file_1,'w')
    fout2 = open(file_2,'w')
    fout3 = open(file_3,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        if(mm.has_key(line)):
            mm[line] += 1
        else:
            mm.setdefault(line,1)
    print(len(mm))
    mm1 = {}
    for line in open(file_in1):
        line = line.strip()
        if(mm.has_key(line)):
            fout3.write(line+"\n")
            mm.pop(line)
        else:
            mm1.setdefault(line,1)
    print(len(mm1))
    for key in mm:
        fout1.write(key+"\t" + str(mm[key])+ "\n")
        
def ddd1(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vv = line.strip().split("\t")
        mm.setdefault(vv[0],vv[1])
    for line in open(file_in1):
        vv = line.strip().split("\t")
        if(mm.has_key(vv[0])):
            mm[vv[0]] += "\t" + vv[1]
    for key in mm:
        fout.write(key + "\t" + mm[key] + "\n")

def tj_bl_ty(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vv = line.strip().split("\t")
        ii = int(vv[1])
        if(mm.has_key(ii)):
            mm[ii] += 1
        else:
            mm.setdefault(ii,1)
    for key in mm:
        fout.write(str(key) + "\t" + str(mm[key]) + "\n")
        
def tj_bl_gy(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vv = line.strip().split("\t")
        i1 = int(vv[1])
        i2 = int(vv[2])
        i_total = i1 + i2
        ik = float(i1) / i2
        if(mm.has_key(i_total)):
            mm[i_total] += 1
        else:
            mm.setdefault(i_total,1)
        if(mm1.has_key(ik)):
            mm1[ik] += 1
        else:
            mm1.setdefault(ik,1)
    for key in mm:
        fout.write(str(key) + "\t" + str(mm[key]) + "\n")
    for key in mm1:
        fout1.write(str(key) + "\t" + str(mm1[key]) + "\n")
             
def sort_aa(fine_in,file_name):
    in_aa = "sort/" + file_name + "/"
    out_aa = "del/" + file_name + "/"
    mm_file = {}
    for line in open(fine_in):
        aa = line.strip()
        mm_file.setdefault(aa,1)
    mm = {}
    for key in mm_file:
        in_bb = in_aa + key +".sort"
        out_bb = open(out_aa + key,'w')
        for line in open(in_bb):
            vec1 = line.strip().split("\t")
            reads = vec1[1]
            if(mm.has_key(reads)):
                mm[reads] += 1
            else:
                mm.setdefault(reads,1)
        for key in mm:
            out_bb.write(key + "\t" + str(mm[key]) + "\n")
    return 0
   
def del_re_128(file_name):
    in_m = "del/bl_m/" + file_name + ".sort"
    in_f = "del/bl_f/" + file_name + ".sort"
    out_all = open("del/bl/gy." + file_name,'w')
    out_m = open("del/bl/m." + file_name,'w')
    out_f = open("del/bl/f." + file_name,'w')
    mm = {}
    for line in open(in_m):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0],vec1[1])
    for line in open(in_f):
        vec1 = line.strip().split("\t")
        id = vec1[0]
        if(mm.has_key(id)):
            out_all.write(id + "\t" + mm[id] + "\t" + vec1[1] + "\n")
            mm.pop(id)
        else:
            out_f.write(id + "\t" + vec1[1] + "\n")
    for key in mm:
        out_m.write(key + "\t" + mm[key] + "\n")
    return 0

def del_kk(file_in,file_out):
    fout = open(file_out,'w')
    ss = ""
    ic = 0
    for line in open(file_in):
        vv = line.strip().split("\t")
        s1 = vv[1]
        if(ss == ""):
            ss = s1
            ic = 1
        elif(ss == s1):
            ic += 1
        else:
            fout.write(ss + "\t" + str(ic) + "\n")
            ss = s1
            ic = 1
        
def tj_aa(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vv = line.strip().split("\t")
        ii = int(vv[1])
        if(mm.has_key(ii)):
            mm[ii] += 1
        else:
            mm.setdefault(ii,1)
    for key in mm:
        fout.write(str(key) + "\t" + str(mm[key]) + "\n")
        
''' 
瀛楃涓叉埅鍙栫▼搴�
'''
def cut_100(file_in,file_out):
    fout = open(file_out,'w')
    ii=1
    id = ""
    for line in open(file_in):
        line = line.strip()
        if(ii == 1):
            id = line
        else:
            ii = 0
            ik = len(line) / 100 + 1
            for i in range(0,ik):
                i_re = line[i*100:(i+1)*100]
                fout.write(id+":" + str(i) + "\n")
                fout.write(i_re+"\n")
        ii += 1
        
def create_fa_aa(file_in,file_out):
    fout = open(file_out,'w')
    ii=1
    for line in open(file_in):
        vec = line.strip().split("\t")
        fout.write(">"+str(ii)+"_"+vec[1]+"\t"+vec[2]+"\n")
        fout.write(vec[0]+"\n")
        ii += 1
        
def get_n_site(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        if(len(line) > 100):
            for ii in range(0,len(line)):
                ss = line[ii]
                stst = line[ii + 1:ii + 61]
                if(ss == "N"):
                    fout.write(str(ii + 2) + "\n")    
                    print(str(ii + 2) + "\t" + stst)   
                    
def oo(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in1):
        line = line.strip()
        mm.setdefault(line,1)
    print(len(mm))
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        if(len(vec1) > 10):
            q_id = vec1[0]
            r_id = vec1[3]
            if(mm.has_key(r_id)):
                fout.write(q_id+"\n")
              
def dd(file_in,file_in1,file_in2,file_out):
    fout = open(file_out,'w')
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm1.setdefault(vec[0],"")
        mm2.setdefault(vec[1],vec[0])
    print(len(mm1))
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        if(len(vec1) > 10):
            q_id = vec1[0]
            r_id = vec1[3]
            if(mm1.has_key(r_id)):
                mm1[r_id] += q_id
    for line in open(file_in2):
        vec1 = line.strip().split("\t")
        if(mm2.has_key(vec1[0])):
            cc = vec[1]
            id = mm2[vec1[0]]
            dd = mm1[id]
            vec2 = dd.split("_")
            tt = vec2[1]
            fout.write(id + "\t" + dd + "\t" + cc + "\t" + dd +"\n")
    return 0
                  
def cat_fa(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    s_k = ""
    for line in open(file_in):
         vec1 = line.strip().split("\t")
         ii = len(s_k) + 1
         fout1.write(str(ii) + "\t" + line.strip() + "\n")
         s_k += vec1[0] + "N"
    fout.write(">gy\n")
    fout.write(s_k)
    
def is_xy(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
         line = line.strip()
         mm.setdefault(line,"1")
    for line in open(file_in1):
         line = line.strip()
         if(mm.has_key(line)):
             pass
         else:
             fout.write(line+"\n")
             
             
def find_snp(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        scaff = vec[0]
        site = vec[1]
        id = scaff + ":" + site
        snp = vec[3] + ":" +vec[4]
        mm.setdefault(id,snp)
        if(mm1.has_key(scaff)):
            mm1[scaff] += site + ","
        else:
            mm1.setdefault(scaff,site+",")
            
    for line in open(file_in1):
        vec = line.strip().split(":")
        scaff = vec[0]
        site = int(vec[1])
        vec1 = mm1.get(scaff)[0:-1].split(",")
        s_out = ""
        i_t = 0
        for ss in range(0,len(vec1)):
            i_ss = int(vec1[ss])
            if(site <= i_ss and site + 60 >= i_ss):
                id = scaff + ":" + str(i_ss)
                snp = mm.get(id)
                s_out += str(i_ss) + "~" + snp + "\t"
                i_t += 1
        if(len(s_out) > 2):
            fout.write(line.strip() + "\t" + str(i_t) + "\t" + s_out + "\n")
    return 0
            
def is_one(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
         line = line.strip()
         mm.setdefault(line,"1")
    for line in open(file_in1):
         line = line.strip()
         if(mm.has_key(line)):
             pass
         else:
             fout.write(line+"\n")
                     
def is_onesite(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
         line = line.strip()
         mm.setdefault(line,"1")
    
    for line in open(file_in1):
         vec = line.strip().split("\t")
         id = vec[0]
         
         if(mm.has_key(id)):
            ii = vec[2] + ":" + vec[3]
            if(mm1.has_key(ii)):
                mm1[ii] += id + "\t"
            else:
                mm1.setdefault(ii,id+"\t")
    for key in mm1:
        vec = key.split("\t")
        fout.write(key + "\t"+str(len(vec)) + "\t"+mm1[key] + "\n")
        
def is_onesite_new(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
         vec = line.strip().split(":")
         line = line.strip()
         id = vec[0]
         ii = int(vec[1])
         for kk in range(0,59):
             i1 = ii + kk
             s_id = id + ":" + str(i1)
             mm.setdefault(s_id,line)
    print("load over")
    for line in open(file_in1):
         vec = line.strip().split("\t")
         id = vec[0] + ":" + vec[1]
         if(mm.has_key(id)):
             snp = vec[3] + ":" + vec[4]
             fout.write(mm[id] + "\t" + snp + "\t" + vec[0] + "\t" + vec[1] + "\n")
       
def is_onesite_del(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
         vec = line.strip().split(":")
         id = vec[0]
         ii = vec[3] + "~" + vec[1] + "\t"
         if(mm.has_key(id)):
             mm[id] += ii
         else:
             mm.setdefault(id,ii)
    for key in mm:
        vec = mm[key].split("\t")
        fout.write(key + "\t" + str(len(vec) -1) + "\t" + mm[key] + "\n")
 
def print_awk(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        snp = vec[0]
        site = vec[1]
        ii = int(site) - 180
        fout.write("echo '>" + site + "_" + snp + "'\n")
        fout.write("awk '{print(substr($0,"+str(ii)+",361));}' k2\n")  
    return 0

def get_only_one_fa(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[2],1)
    id = ""
    for line in open(file_in1):
        line = line.strip()
        ik = line[1:]
        if(mm.has_key(ik)):
            id = ik
        elif(len(id) > 0):
            fout.write(">" + id + "\n")
            fout.write(line + "\n")
            id = ""
    return 0

def get_num(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        s_ref = vec[0][-3]
        snp_ref = vec[0][-1]
        i_ref = 0
        i_snp = 0
        dp = int(vec[3])
        pp = vec[4]
        for ii in range(0,len(pp)):
            kk = pp[ii]
            if(kk == "." or kk == ","):
                i_ref += 1
            elif(kk.upper() == snp_ref):
                i_snp += 1
        fout.write(vec[0] + "\t" + vec[3] + "\t" + s_ref + "\t" + str(i_ref) + "\t" + snp_ref  + "\t" + str(i_snp) + "\n")
    return 0


def chai_3_181(file_in,file_out):
    fout = open(file_out,'w')
    ii = 1
    id = ""
    for line in open(file_in):
        line = line.strip()
        if(ii == 1):
            id = line
        else:
            ii = 0
            s1 = line[0:150]
            s2 = line[211:]
            s3 = line[106:256]
            fout.write(id+"_1\n")
            fout.write(s1+"\n")
            fout.write(id+"_2\n")
            fout.write(s3+"\n")
            fout.write(id+"_3\n")
            fout.write(s2+"\n")
        ii += 1
    return 0


def duiyin_aa(file_in,file_in1,file_in2,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    ma = {}
    mb = {}
    for line in open(file_in):
        line = line.strip()
        ma.setdefault(line,1)
    for line in open(file_in1):
        line = line.strip()
        mb.setdefault(line,1)
    for line in open(file_in2):
        vec = line.strip().split("\t")
        v1_id = vec[1]
        v1_reads = vec[2]
        v2_id = vec[3]
        v2_reads = vec[4]
        if(ma.has_key(v1_id)):
            if(mb.has_key(v2_id)):
                fout.write(line.strip() + "\t")
                fout1.write(">" + v1_id + "\n")
                fout1.write(v1_reads+"\n")
                fout2.write(">" + v2_id + "\n")
                fout2.write(v2_reads+"\n")
    return 0

def snp_aa(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[0].split(":")
        site_star = int(vec1[1])
        a1 = vec[2]
        a2 = vec[4]
        no_pair = 0
        no_site = 0
        s_snp = ""
        for kk in range(0,60):
            k1 = a1[kk]
            k2 = a2[kk]
            if(k1 != k2):
                no_pair += 1
                no_site = kk
                s_snp = k1 + ":" + k2
        if(no_pair == 1):
            ref_site = site_star + no_site
            ref_id =vec1[0] + ":"+str(ref_site)
            if(mm.has_key(ref_id)):
                mm[ref_id] += str(site_star) + ":"
            else:
                mm.setdefault(ref_id,str(site_star) + ":")
            fout.write(line.strip() + "\t" + str(no_site) + "\t" +ref_id +"\t"+ s_snp + "\n")
    for key in mm:
        kk = mm[key][0:-1]
        vec = kk.split(":")
        vec1 = key.split(":")
        i_l = len(vec)
        ik = vec[i_l/2]
        print(vec1[0] + ":" + ik)
        fout1.write(key + "\t" + str(len(vec)) + "\t" + kk + "\n")
    return 0


def is_one_site(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        if(len(vec) > 10):
            site = vec[2] + ":" + vec[3]
            if(mm.has_key(site)):
                mm[site] += vec[0] + ":"
            else:
                mm.setdefault(site,vec[0] + ":")
    for key in mm:
        ii = mm[key]
        vec = ii[0:-1].split(":")
        mm1 = {}
        ss = ""
        if(len(vec) == 5):
            for jj in range(0,5):
                ik = vec[jj]
                ia = ik[0:-2]
                ss = ia
                mm1.setdefault(ia,1)
            if(len(mm1) == 1):
                fout.write(key + "\t" + ss + "\t" + ss[0] + "\n")
    return 0    
        
def find_snp18(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            fout.write(line.strip() + "\n")
            
def is_3bp(file_in,file_out):
    fout = open(file_out,'w')
    id = ""
    i_begin = ""
    mm = {}
    for line in open(file_in):
        vec = line.strip().split(":")
        i_cc = int(vec[1])
        if(id == ""):
            id = vec[0]
            i_begin = i_cc
        elif(id == vec[0]):
            if(i_cc - i_begin <= 3):
                print(vec[0]+"\t"+str(i_begin)+"\t" + str(i_cc))
                ii1 = vec[0] + ":" + str(i_begin)
                ii2 = vec[0] + ":" + str(i_cc)
                if(mm.has_key(ii1)):
                    mm[ii1] += 1
                else:
                    mm.setdefault(ii1,1)
                if(mm.has_key(ii2)):
                    mm[ii2] += 1
                    
                else:
                    mm.setdefault(ii2,1)
        i_begin = i_cc
        id = vec[0]
    for key in mm:
        fout.write(key + "\t" + str(mm[key]) + "\n")
    return 0

def match_1(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[1])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            fout.write(line.strip() + "\t" + mm[vec[0]] + "\n")
            vec1 = vec[0].split(":")
            st = "samtools mpileup -A -B -C 0 -q 0 -Q 0 -t DP,AD --region "+vec1[0]+":"+vec1[1]+"-"+vec1[1]+" -f v12.fa del.sort.bam >> snp.detail"
            fout1.write(st+"\n")
        else:
            print(vec[0])
    return 0

def match_11(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        snp_site = int(vec[0].split(":")[1])
        i_begin = snp_site - 61
        i_end = snp_site
        vec1 = vec[2].split(":")
        str_i = ""
        for ii in range(0,len(vec1)):
            i_site = int(vec1[ii])
            if(i_site >= i_begin and i_site <= i_end):
                str_i += vec1[ii] + ":"
            else:
                print(line)
        fout.write(vec[0] + "\t" + vec[1] + "\t" + str_i[0:-1] + "\t" + vec[3] + "\n")
    return 0

def match_12(file_in,file_in1):
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[1])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[6]
        snp = vec[7]
        if(mm.has_key(id)):
            snpw1 = mm[id]
            snpw2 = snpw1[2] + ":" + snpw1[0]
            if(snp != snpw1 and snp != snpw2):
                print("-\t" + line.strip())
    return 0

def match_2(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],line.strip())
    print("load over")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        scaff = vec[0].split(":")[0]
        vec1 = vec[2].split(":")
        f_id_l = scaff + ":" + vec1[0]
        k_l = mm.get(f_id_l)
        
        if(len(vec1) == 1):
            kk_1 = k_l.split("\t")
            fout.write(line.strip() + "\t" + kk_1[2] + "\t" +kk_1[5] + "\t0\t0\t60\n")
        else:
            f_id_r = scaff + ":" + vec1[len(vec1) - 1]
            k_r = mm.get(f_id_r)
            kk_1 = k_l.split("\t")
            kk_2 = k_r.split("\t")
            site_l = int(kk_1[5])
            site_r = int(kk_2[5])
            str_l = kk_1[2][0:site_l]
            str_r = kk_2[2][site_r:]
            str1 = str_l + str_r
            fout.write(line.strip() + "\t" + str1 + "\t" + kk_1[5]  + "\t" + str(site_l) + "\t" + str(59-site_r) + "\t" + str(len(str1)) + "\n")
            
        is_write = 0
        for ik in range(0,len(vec1)):
            f_id = scaff + ":" + vec1[ik]
            ka_l = mm.get(f_id)
            kk_1 = ka_l.split("\t")
            snp_site = int(kk_1[5])
            if(snp_site >= 10 and snp_site <= 50 and is_write == 0):
                fout1.write("1\t" + line.strip() +"\t" + ka_l+ "\n")
                is_write = 1
        if(is_write == 0):
            fout1.write("2\t" + line.strip() +"\t" + k_l+ "\n")
    return 0

def cut_aa(file_out):
    ss = "CAATGAAGCAATTTAACTTTCCTTCATTACTAGTTCTAAAGAGAAATTTTAGAATTAGTAACGGCTTGGGCTCAAGTTTTAAATTGTCAACTTGAGATTTGAATCCGTGGAGGAAGGAAGTAGTTCCTTACAAACGAGGGGTTTTAAGACACTCTGTTGTTGTTTCTTATTTATTTACTCACTTGTGCCATCAAACTGTTGTATAAATGCAATATCACGTGGCTTGTGGCGTGCTGATATACAGCCATATCACATGGCTACAAGTGTGATACTGCTCATAAATTCAGTGAAATAATTATGCTTCAGTATAATATTAAAAAGAAAACCATCACCTCAATCAATCAATCAATCAATCAATAATATTATTATTAATAATAATAATTTGTCAAAATGATAATTTG"
    fout = open(file_out,'w')
    ii = 0
    while(ii < len(ss)):
        si = ss[ii:ii+60]
        ii += 30
        fout.write(">z2_"+ str(ii) + "\n")
        fout.write(si+"\n")
    ss = "CAAATTATCATTTTGACAAATTATTATTATTAATAATAATATTATTGATTGATTGATTGATTGATTGAGGTGATGGTTTTCTTTTTAATATTATACTGAAGCATAATTATTTCACTGAATTTATGAGCAGTATCACACTTGTAGCCATGTGATATGGCTGTATATCAGCACGCCACAAGCCACGTGATATTGCATTTATACAACAGTTTGATGGCACAAGTGAGTAAATAAATAAGAAACAACAACAGAGTGTCTTAAAACCCCTCGTTTGTAAGGAACTACTTCCTTCCTCCACGGATTCAAATCTCAAGTTGACAATTTAAAACTTGAGCCCAAGCCGTTACTAATTCTAAAATTTCTCTTTAGAACTAGTAATGAAGGAAAGTTAAATTGCTTCATTG"
    ii = 0
    while(ii < len(ss)):
        si = ss[ii:ii+60]
        ii += 30
        fout.write(">f2_"+ str(ii) + "\n")
        fout.write(si+"\n")
    '''
    s1 = ss[0:60]
    s2 = ss[40:100]
    s3 = ss[80:140]
    s4 = ss[120:180]
    s5 = ss[160:220]
    s6 = ss[200:260]
    s7 = ss[240:300]
    s8 = ss[280:340]
    s9 = ss[320:380]
    s10 = ss[360:420]
    s11 = ss[400:460]
    s12 = ss[420:480]
    s13 = ss[460:-1]
    s1 = ss[0:60]
    s2 = ss[30:90]
    s3 = ss[60:-1]
    s4 = ss[90:150]
    s5 = ss[120:180]
    s6 = ss[150:210]
    s7 = ss[180:240]
    s8 = ss[210:270]
    s9 = ss[240:300]
    s10 = ss[270:330]
    s11 = ss[300:360]
    s12 = ss[330:390]
    s13 = ss[360:420]
    s14 = ss[390:440]
    s15 = ss[410:470]
    s16 = ss[440:500]
    for ii in range(0,5):
        fout.write(">11_"+str(ii)+"\n" +s1 + "\n")
        fout.write(">12_"+str(ii)+"\n" +s2 + "\n")
        fout.write(">13_"+str(ii)+"\n" +s3 + "\n")
        fout.write(">4_"+str(ii)+"\n" +s4 + "\n")
        fout.write(">5_"+str(ii)+"\n" +s5 + "\n")
        fout.write(">6_"+str(ii)+"\n" +s6 + "\n")
        fout.write(">7_"+str(ii)+"\n" +s7 + "\n")
        fout.write(">8_"+str(ii)+"\n" +s8 + "\n")
        fout.write(">9_"+str(ii)+"\n" +s9 + "\n")
        fout.write(">10_"+str(ii)+"\n" +s10 + "\n")
        fout.write(">11_"+str(ii)+"\n" +s11 + "\n")
        fout.write(">12_"+str(ii)+"\n" +s12 + "\n")
        fout.write(">13_"+str(ii)+"\n" +s13 + "\n")
        fout.write(">14_"+str(ii)+"\n" +s13 + "\n")
        fout.write(">15_"+str(ii)+"\n" +s13 + "\n")
        fout.write(">16_"+str(ii)+"\n" +s13 + "\n")
        '''
        
def ti_aa(file_in,file_out1,file_out2):
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        if(vec[1] == "6"):
            vec1 = vec[2].split(":")
            mm = {}
            ik = ""
            for ii in range(0,5):
                kk = vec1[ii][0:-2]
                ik = kk
                mm.setdefault(kk,1)
            if(len(mm) == 1):
                fout1.write(vec[0] + "\t" + ik + "\n")
        elif(vec[1] == "11"):
            vec1 = vec[2].split(":")
            mm = {}
            ik = ""
            ik1 = ""
            for ii in range(0,10):
                kk = vec1[ii][0:-2]
                if(ik == ""):
                    ik = kk
                elif(ik != kk):
                    ik1 = kk
                mm.setdefault(kk,1)
            if(len(mm) == 2):
                fout2.write(vec[0] + "\t" + ik + "\t" + ik1 + "\n")
    return 0
          
def fa(file_in,file_out):
    fout = open(file_out,'w')
    ii = 0
    for line in open(file_in):
        fout.write(">"+str(ii)+ "\n")
        fout.write(line.strip() + "\n")
        ii += 1
       
def xy_tiqu(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[3]
        id = vec[0]
        if(mm.has_key(site)):
            mm[site] += id + ":"
        else:
            mm.setdefault(site,id+":")
    mmcc = {}
    for key in mm:
        vv = mm[key][0:-1]
        vec = vv.split(":")
        if(len(vec) == 5):
            mmk = {}
            ik = ""
            for ii in range(0,len(vec)):
                vec1 = vec[ii].split("_")
                mmk.setdefault(vec1[0],1)
                ik = vec1[0] + "_" + vec1[2]
            if(len(mmk) == 1):
                vecw = ik.split("_")
                if(mmcc.has_key(vecw[1])):
                    mmcc[vecw[1]] += 1
                else:
                    mmcc.setdefault(vecw[1],1)
                fout.write(key + "\t" + mm[key] + "\n")
    for key in mmcc:
        fout.write(key + "\t" + str(mmcc[key]) + "\n")
         
         
def map_xy(file_in,file_in1,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[1])
    print("load over")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if vec[0] in mm:
            kk = mm.pop(vec[0])
            fout2.write(line.strip() + "\t" + kk + "\n")
        else:
            fout1.write(line.strip() + "\n")
    for key in mm:
        fout.write(key + "\t" + mm[key] + "\n")
    return 0

def find_x5_x10(file_in,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        if(len(vec) > 6):
            id = vec[0]
            scaff = vec[2]
            site = vec[3]
            ii = scaff + ":" + site
            if(mm.has_key(ii)):
                mm[ii] += id + ":"
            else:
                mm.setdefault(ii,id+":")
    for key in mm:
        vec = mm[key][0:-1].split(":")
        i_site = len(vec)
        mmk = {}
        for ii in range(0,i_site):
            ik = vec[ii][0:-2]
            mmk.setdefault(ik,1)
        i_mm = len(mmk)
        s_oo = ""
        for key1 in mmk:
            s_oo += key1 + ":"
        if(i_site == 10 and i_mm == 2):
            fout.write(key + "\t" + s_oo[0:-1] + "\n")
        elif(i_site == 5 and i_mm == 1):
            fout1.write(key + "\t" + s_oo[0:-1] + "\n")
        elif(i_site == 15 and i_mm == 3):
            fout2.write(key + "\t" + s_oo[0:-1] + "\n")
    return 0
        
def get_pv(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[1].split(":")
        mm.setdefault(vec[0],vec[1])
        for ii in range(0,len(vec1)):
            kk = vec1[ii]
            mm1.setdefault(kk,"")
            mm2.setdefault(kk,"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        reads = vec[1]
        id = vec[0]
        cc = vec[2]
        if(mm1.has_key(id)):
            mm1[id] = reads
            mm2[id] = cc
    for key in mm:
        value = mm[key]
        vec1 = value.split(":")
        s1 = ""
        s2 = ""
        for ii in range(0,len(vec1)):
            kk = vec1[ii]
            reads = mm1[kk]
            cc = mm2[kk]
            s1 += reads + ":"
            s2 += cc + ":"
        fout.write(key +"\t" + value + "\t" + s2[0:-1] + "\t" + s1[0:-1] + "\n")
    return 0

def chang_snp(old_snp):
    new_snp = old_snp
    if(old_snp == "T:A"):
        new_snp = "A:T"
    elif(old_snp == "C:A"):
        new_snp = "A:C"
    elif(old_snp == "G:A"):
        new_snp = "A:G"
    elif(old_snp == "C:T"):
        new_snp = "T:C"
    elif(old_snp == "G:T"):
        new_snp = "T:G"
    elif(old_snp == "G:C"):
        new_snp = "C:G"
    return new_snp

def x10_snp(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[3].split(":")
        vec2 = vec[1].split(":")
        vec3 = vec[0].split(":")
        scaff = vec3[0]
        be_site = int(vec3[1])
        reads1 = vec1[0]
        reads2 = vec1[1]
        is_diff = 0
        snp_site = 0
        snp_type = ""
        ii_kk = 0
        for ii in range(0,len(reads1)):
            ik = reads1[ii]
            ak = reads2[ii]
            if(ik != ak):
                is_diff += 1
                snp_site = scaff + ":" + str(be_site + ii)
                snp_type = ik + ":" + ak
                ii_kk = ii
        snp_type = chang_snp(snp_type)
        if(is_diff == 1):
            fout.write(line.strip() + "\t" + snp_site + "\t" + snp_type + "\t" + str(ii_kk) + "\n")
            fout1.write(">" + vec2[0] + "_a\n" + vec1[0] + "\n" + ">" + vec2[0] + "_b\n" + vec1[0] + "\n" +">" + vec2[0] + "_c\n" + vec1[0] + "\n" +">" + vec2[0] + "_d\n" + vec1[0] + "\n" +">" + vec2[0] + "_e\n" + vec1[0] + "\n")
            fout1.write(">" + vec2[1] + "_a\n" + vec1[1] + "\n" + ">" + vec2[1] + "_b\n" + vec1[1] + "\n" +">" + vec2[1] + "_c\n" + vec1[1] + "\n" +">" + vec2[1] + "_d\n" + vec1[1] + "\n" +">" + vec2[1] + "_e\n" + vec1[1] + "\n")
    return 0

def x15_snp(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[3].split(":")
        vec2 = vec[2].split(":")
        vec3 = vec[1].split(":")
        vec4 = vec[0].split(":")
        scaff = vec4[0]
        be_site = int(vec4[1])
        
        mm = {}
        mm.setdefault(int(vec2[0]),0)
        mm.setdefault(int(vec2[1]),1)
        mm.setdefault(int(vec2[2]),2)
        sk = ""
        if(len(mm) == 3):
            mk = sorted(mm.items(),key=lambda x:x[0],reverse=False)
            b_1 = mk[1][1]
            b_2 = mk[2][1]
            b_3 = mk[0][1]
            reads1 = vec1[b_1]
            reads2 = vec1[b_2]
            reads3 = vec1[b_3]
            
            is_diff = 0
            snp_site = 0
            snp_type = ""
            ii_kk = 0
            for ii in range(0,len(reads1)):
                ik = reads1[ii]
                ak = reads2[ii]
                bk = reads3[ii]
                if(ik != ak or ik != bk or ak != bk):
                    is_diff += 1
                    snp_site = scaff + ":" + str(be_site + ii)
                    snp_type = ik + ":" + ak
                    ii_kk = ii
            snp_type = chang_snp(snp_type)
            if(is_diff == 1):
                sk = vec3[b_1] + ":" + vec3[b_2] + "\t" + vec2[b_1] + ":" + vec2[b_2] + "\t" + vec1[b_1] + ":" + vec1[b_2]
                fout.write(vec[0] + "\t" + sk + "\t" + snp_site + "\t" + snp_type + "\t" + str(ii_kk) + "\n")
                fout1.write(">" + vec3[b_1] + "_a\n" + vec1[b_1] + "\n" + ">" + vec3[b_1] + "_b\n" + vec1[b_1] + "\n" +">" + vec3[b_1] + "_c\n" + vec1[b_1] + "\n" +">" + vec3[b_1] + "_d\n" + vec1[b_1] + "\n" +">" + vec3[b_1] + "_e\n" + vec1[b_1] + "\n")
                fout1.write(">" + vec3[b_2] + "_a\n" + vec1[b_2] + "\n" + ">" + vec3[b_2] + "_b\n" + vec1[b_2] + "\n" +">" + vec3[b_2] + "_c\n" + vec1[b_2] + "\n" +">" + vec3[b_2] + "_d\n" + vec1[b_2] + "\n" +">" + vec3[b_2] + "_e\n" + vec1[b_2] + "\n")
    return 0


def x15_snp_aa(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[3].split(":")
        vec2 = vec[2].split(":")
        vec3 = vec[1].split(":")
        vec4 = vec[0].split(":")
        scaff = vec4[0]
        be_site = int(vec4[1])
        
        mm = {}
        mm.setdefault(int(vec2[0]),0)
        mm.setdefault(int(vec2[1]),1)
        mm.setdefault(int(vec2[2]),2)
        sk = ""
        if(len(mm) == 3):
            mk = sorted(mm.items(),key=lambda x:x[0],reverse=False)
            b_1 = mk[1][1]
            b_2 = mk[2][1]
            reads1 = vec1[b_1]
            reads2 = vec1[b_2]
            
            is_diff = 0
            snp_site = 0
            snp_type = ""
            ii_kk = 0
            for ii in range(0,len(reads1)):
                ik = reads1[ii]
                ak = reads2[ii]
                if(ik != ak):
                    is_diff += 1
                    snp_site = scaff + ":" + str(be_site + ii)
                    snp_type = ik + ":" + ak
                    ii_kk = ii
            snp_type = chang_snp(snp_type)
            if(is_diff == 1):
                sk = vec3[b_1] + ":" + vec3[b_2] + "\t" + vec2[b_1] + ":" + vec2[b_2] + "\t" + vec1[b_1] + ":" + vec1[b_2]
                fout.write(vec[0] + "\t" + sk + "\t" + snp_site + "\t" + snp_type + "\t" + str(ii_kk) + "\n")
                fout1.write(">" + vec3[b_1] + "_a\n" + vec1[b_1] + "\n" + ">" + vec3[b_1] + "_b\n" + vec1[b_1] + "\n" +">" + vec3[b_1] + "_c\n" + vec1[b_1] + "\n" +">" + vec3[b_1] + "_d\n" + vec1[b_1] + "\n" +">" + vec3[b_1] + "_e\n" + vec1[b_1] + "\n")
                fout1.write(">" + vec3[b_2] + "_a\n" + vec1[b_2] + "\n" + ">" + vec3[b_2] + "_b\n" + vec1[b_2] + "\n" +">" + vec3[b_2] + "_c\n" + vec1[b_2] + "\n" +">" + vec3[b_2] + "_d\n" + vec1[b_2] + "\n" +">" + vec3[b_2] + "_e\n" + vec1[b_2] + "\n")
    return 0

def Get_vcf_snp(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0] + ":" + vec[1]
        ref_site = vec[3]
        snp_site = vec[4]
        snp = ""
        if(len(snp_site) == 1):
            snp = ref_site + ":" + snp_site
        elif(len(snp_site) == 3):
            snp = snp_site[0] + ":" + snp_site[2]
        else:
            print(line)
        snp = chang_snp(snp)
        fout.write(site + "\t" + snp + "\t" + ref_site+"\n")
        
def select_one(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[1])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[4]
        snp = vec[5]
        if(mm.has_key(id)):
            y_snp = mm[id]
            if(y_snp == snp):
                fout.write(line.strip() + "\n")
                mm.pop(id)
    for key in mm:
        print(key+"\t"+mm[key]);
        
def test_aa(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[3]
        read1 = vec1[0]
        read2 = vec1[1]
        mm1.setdefault(read1,id)
        mm2.setdefault(read2,id)
        mm2.setdefault(id,line.strip())
    for line in open(file_in1):
        vec = line.strip().split("\t")
        rr = vec[0]
        if(mm1.has_key(rr)):
            id = mm1[rr]
            mm2.pop(rr)
        elif(mm2.has_key(rr)):
            id = mm2[rr]
            mm2.pop(rr)
            
    for key in mm2:
        fout.write(mm2[key] + "\n");
        
def test_bb(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2].split(":")
        i1 = int(vec1[0])
        i2 = int(vec1[1])
        if(i1 > 10 and i2 > 10):
            fout.write(line.strip() + "\n")
        

def count_a_b(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[3].split(":")
        reads1 = vec1[0]
        reads2 = vec1[1]
        mm1.setdefault(reads1,0)
        mm1.setdefault(reads2,0)
        mm2.setdefault(id,line.strip())
        
    for line in open(file_in1):
        vec = line.strip().split("\t")
        reads = vec[0]
        if(mm1.has_key(reads)):
            ii = int(vec[1])
            mm1[reads] += ii
            
    for key in mm2:
        line = mm2[key]
        vec = line.split("\t")
        vec1 = vec[3].split(":")
        reads1 = vec1[0]
        reads2 = vec1[1]
        i1 = mm1[reads1]
        i2 = mm1[reads2]
        ito = i1 + i2
        snp_id = vec[6]
        snp = vec[4]
        message = vec[0] + "~" +vec[1] + "~" + vec[2] + "~" + vec[3] + "~" + vec[5] + "~" + vec[8]
        fout.write(snp_id + ":" + snp + "\t" + message + "\t" + str(i1) + ":" + str(i2) + str(ito) + "\n")
        
def test_os_before(file_in,file_out,file_name):
    fout = open(file_out,'w')
    for line in open(file_in):
        llk = line.strip()
        vec = line.strip().split("\t")
        site = vec[0].split(":")[1]
        i_site = int(site)
        ss = vec[0] + "~" + vec[1] + "~" + vec[2]
        #val = os.popen("samtools view  k5.sort.bam Lachesis_group0__131_contigs__length_56955660:505-505")
        val = os.popen("samtools view " + file_name + " " + vec[0] + "-" + site)
        mm = {}
        for line in val.readlines():
            vecw = line.strip().split("\t")
            i_begin = int(vecw[3])
            i_d = i_site - i_begin
            if(mm.has_key(id)):
                mm[id] += vecw[0][0:-2] + ":"
            else:
                mm.setdefault(id,vecw[0][0:-2] + ":")    
        mk = sorted(mm.items(),key=lambda x:x[0],reverse=False)
        value = list(mk)[0][1]
        mm1 = {}
        vec1 = value[0:-1].split(":")
        for ii in range(0,len(vec1)):
            mm1.setdefault(vec1[ii],1)
        mw = sorted(mm1.items(),key=lambda x:x[0],reverse=False)
        small = list(mw)[0][0]
        fout.write(llk + "\t" + small + "\n")
    return 0

def test_os(file_in,file_out,file_out1,file_out2,file_name):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    for line in open(file_in):
        line_aa = line.strip()
        vec = line.strip().split("\t")
        site = vec[0].split(":")[1]
        ss = vec[0] + "~" + vec[1] + "~" + vec[2]
        #val = os.popen("samtools view  k5.sort.bam Lachesis_group17__21_contigs__length_31864045:28775586-28775586")
        val = os.popen("samtools view " + file_name + " " + vec[0] + "-" + site)
        
        mm = {}
        mm1 = {}
        mm2 = {}
        for line in val.readlines():
            vecw = line.strip().split("\t")
            fout.write(ss  + "\t" + line.strip() + "\n")
            id = vecw[0][0:-2]
            site = int(vecw[3])
            reads = vecw[9]
            mm.setdefault(site,id)
            mm1.setdefault(id,line.strip())
            if(reads[0] == "A"):
                mm2.setdefault(site,id)
        if(len(mm) >= 1):
                mk = sorted(mm.items(),key=lambda x:x[0],reverse=False)
                mkw = sorted(mm2.items(),key=lambda x:x[0],reverse=False)
                i_len = len(mkw) / 2
                site = list(mkw)[i_len][0]
                id = list(mkw)[i_len][1]
                str_map = mm1[id]
                #print(ss,site,id,str_map);
                fout1.write(ss + "\t" + str(site) + "\t" + str(id) + "\t" + str_map + "\n")
                str_min = mm1[list(mk)[0][1]]
                str_max = mm1[list(mk)[len(mk) - 1][1]]
                fout2.write(ss + "~" + str_min + "~" + str_max + "\n")
        else:
            print(line_aa);
    return 0

def get_reads_ww(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[1].split(":")
        reads1 = vec1[0]
        reads2 = vec1[1]
        mm1.setdefault(reads1,"")
        mm1.setdefault(reads2,"")
        mm2.setdefault(id,line.strip())
        
    ii = 1
    is_find = 0
    ik = ""
    for line in open(file_in1):
        line = line.strip()
        if(ii == 1):
            id = line[1:]
            if(mm1.has_key(id)):
                is_find = 1
                ik = id
        else:
            ii = 0
            if(is_find == 1):
                mm1[ik] = line
            ik = ""
            is_find = 0
        ii += 1
    for key in mm2:
        line = mm2[key]
        vec = line.split("\t")
        vec1 = vec[1].split(":")
        id1 = vec1[0]
        id2 = vec1[1]
        reads1 = mm1[id1]
        reads2 = mm1[id2]
        j1 = id1.split("_")[1]
        j2 = id2.split("_")[1]
        
        is_diff = 0
        snp_type = ""
        ii_kk = 0
        for ii in range(0,len(reads1)):
            ik = reads1[ii]
            ak = reads2[ii]
            if(ik != ak):
                is_diff += 1
                snp_type = ik + ":" + ak
                ii_kk = ii
        
        snp_type = chang_snp(snp_type)
        if(is_diff == 1):
            fout.write(line + "\t" + j1 + ":" + j2 + "\t" + reads1 + ":" + reads2 + "\t" + snp_type + "\t" + str(ii_kk) + "\n")
            fout1.write(">" + id1 + "_a\n" + reads1 + "\n" + ">" + id1 + "_b\n" + reads1 + "\n" +">" + id1 + "_c\n" + reads1 + "\n" +">" + id1 + "_d\n" + reads1 + "\n" +">" + id1 + "_e\n" + reads1 + "\n")
            fout1.write(">" + id2 + "_a\n" + reads2 + "\n" + ">" + id2 + "_b\n" + reads2 + "\n" +">" + id2 + "_c\n" + reads2 + "\n" +">" + id2 + "_d\n" + reads2 + "\n" +">" + id2 + "_e\n" + reads2 + "\n")
    return 0

def pair_os(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[2]
        if(not mm.has_key(id)):
            mm.setdefault(id,vec[0])
        else:
            print(id,mm[id],vec[0])
    print(len(mm))
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec1 = vec[1].split(":")
        is_find = 0
        value = ""
        if(mm.has_key(vec1[0])):
            is_find += 1
            value = mm[vec1[0]]
        elif(mm.has_key(vec1[1])):
            is_find += 1
            value = mm[vec1[1]]
        if(is_find > 0):
            fout.write(line.strip() + "\t" + value)
    return 0
        
def tj_kk(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip();
        if(mm.has_key(line)):
            mm[line] += 1
        else:
            mm.setdefault(line,1)
    for key in mm:
        fout.write(mm+"\t" + str(mm[key]))

'''
python tiqu_os.py k5.type 270w.tiqu  k5.sort.bam  寰楀埌杩欎簺snp鐨剅eads锛屽湪270w.pbs鑴氭湰褰撲腑
python tj_pp.py temp/all_no.id k5.map k5.hege   2448291涓�
awk -F '~' '{print($1"\t"$2"\t"$3);}' k5.hege > k5_shai.type
'''
def tj_pp(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        id = line.strip().split("\t")[0]
        mm.setdefault(id,1)
    id = ""
    is_map = 0
    for line in open(file_in1):
        vec = line.strip().split("\t")
        snp = vec[0]
        if(id == ""):
            id = snp
        elif(id != snp):
            if(is_map == 0):
                fout.write(id + "\n")
            id = snp
            is_map = 0
        reads_id = vec[1][0:-2]
        if(mm.has_key(reads_id)):
            is_map = 1
    
def dd_reverse(str_in):
    trans_table = string.maketrans('ATCG','TAGC')
    a1 = str_in.translate(trans_table)
    str_r1 = a1[::-1]
    return str_r1


        
def pj_101(file_in,file_in1,file_out,file_out1,file_out2):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    mm = {}
    mmw = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[1].split(":")
        vec2 = vec[3].split(":")
        i_snp = int(vec[5])
        id1 = vec1[0]
        id2 = vec1[1]
        snp = vec[4]
        mm.setdefault(id1,snp)
        mm.setdefault(id2,snp)
        mmw.setdefault(id1,i_snp)
        mmw.setdefault(id2,i_snp)
    print(len(mm))
    for line in open(file_in1):
        vec = line.strip().split("~")
        print(vec[0])
        #i_beg = int(vec[0].split(":")[1])
        snp = vec[1]
        vec1 = vec[3].split("\t")
        vec2 = vec[4].split("\t")
        snp1 = mm[vec1[0][0:-2]]
        snp2 = mm[vec2[0][0:-2]]
        if(snp1 == snp and snp2 == snp):
            reads1 = vec1[9]
            reads2 = vec2[9]
            i1 = mmw[vec1[0][0:-2]]
            i2 = mmw[vec2[0][0:-2]]
            str_lef = reads1[0:i1]
            str_rig = reads2[i2:]
            str_pj = str_lef + str_rig
            snp_site = i1 + 1
            fout.write(vec[0] + "\t" + vec[1] + "\t" + str(snp_site) + "\t" + str(len(str_lef)) + "\t" + str(len(str_rig) - 1) + "\t" + str(len(str_pj)) + "\t" + str_pj + "\n")
            re_vv = ""
            ku = vec[1][0]
            if(ku == "A"):
                re_vv = "T:"
            elif(ku == "T"):
                re_vv = "A:"
            elif(ku == "G"):
                re_vv = "C:"
            else:
                re_vv = "G:"
            ku = vec[1][2]
            if(ku == "A"):
                re_vv += "T"
            elif(ku == "T"):
                re_vv += "A"
            elif(ku == "G"):
                re_vv += "C"
            else:
                re_vv += "G"
            s1 = str_lef + snp[0] + str_rig[1:]
            s2 = str_lef + snp[2] + str_rig[1:]
            fout2.write(vec[0] + "\t" + s1 + "\t" + s2 +  "\n")
            #print(len(str_rig))
            if(len(str_lef) >= 50 and str(len(str_rig) - 1) >= 50):
                sub_lef = str_lef[-50:]
                sub_rig = str_rig[0:51]
                sks = sub_lef + sub_rig
                fout1.write(vec[0] + "\t+\t" + vec[1] + "\t51\t"+str(len(sks))+"\t" + sks + "\t" + re_vv + "\t" + dd_reverse(sks) + "\n")
            elif(len(str_lef) >= 50):
                sub_lef = str_lef[-50:]
                sks = sub_lef + str_rig
                fout1.write(vec[0] + "\t+\t" + vec[1] + "\t51\t"+str(len(sks))+"\t" + sks + "\t" + re_vv + "\t" + dd_reverse(sks) + "\n")
            elif((len(str_rig) - 1) >= 50):
                sks = str_lef + str_rig[0:51]
                #print(str_lef)
                #print(str_rig)
                #print(str_rig[0:51])
                fout1.write(vec[0] + "\t-\t" + vec[1] + "\t"+str(len(str_lef) + 1)+"\t"+str(len(sks))+"\t" + sks + "\t" + re_vv + "\t"  + dd_reverse(sks) + "\n")
        else:
            print(line.strip())
    return 0

def yz_ff(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm1.setdefault(vec[0],vec[1])
        mm2.setdefault(vec[0],vec[2])
    is_err = 0
    iid = ""
    for line in open(file_in1):
        vec = line.strip().split("~")
        id = vec[0]
        if(mm1.has_key(id)):
            if(iid == ""):
                iid = id
            elif(iid != id):
                if(is_err == 0):
                    fout.write(iid + "\n")
                else:
                    print(iid)
                iid = id
                is_err = 0
                
            reads = vec[2].split("\t")[10]
        
            k1 = mm1[id]
            k2 = mm2[id]
            i1 = k1.find(reads)
            i2 = k2.find(reads)
            if(i1 >= 0 or i2 >= 0):
                pass
            else:
                is_err += 1
    if(is_err == 0):
        fout.write(iid + "\n")
    return 0

def tiqu_os(file_in,file_out,file_name):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0].split(":")[1]
        ss = vec[0] + "~" + vec[1] + "~" + vec[2]
        val = os.popen("samtools view " + file_name + " " + vec[0] + "-" + site)
        
        for line in val.readlines():
            fout.write(ss  + "\t" + line.strip() + "\n")
            
def find_db(file_in,file_in1,file_in2,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,"")
    for line in open(file_in1):
        vec = line.strip().split("~")
        id = vec[0]
        if(mm.has_key(id)):
            vec1 = vec[2].split("\t")
            value = mm[id]
            reads = vec1[10]
            rid = vec1[1][0:-2]
            print(reads)
            if(value == "" and reads[0] == "A"):
                mm[id]= rid
                mm1.setdefault(rid,line.strip())
    for key in mm1:
        fout.write(key + "\t" + mm1[key] + "\n")
    return 0

def select_one_210w(file_in,file_out,file_name):
    fout = open(file_out,'w')
    for line in open(file_in):
        linew = line.strip()
        vec = line.strip().split(":")
        site = vec[1]
        val = os.popen("samtools view " + file_name + " " + vec[0] + ":" + site + "-" + site)
        mm = {}
        for line in val.readlines():
            vecw = line.strip().split("\t")
            id = vecw[0][0:-2]
            site = int(vecw[3])
            reads = vecw[9]
            mm.setdefault(site,id)
        if(len(mm) >= 1):
                mk = sorted(mm.items(),key=lambda x:x[0],reverse=False)
                i_len = len(mk) / 2
                site = list(mk)[i_len][0]
                id = list(mk)[i_len][1]
                fout.write(linew + "\t" + str(site) + "\t" + str(id)  + "\n")
        else:
            print(line);
    return 0

def find_db_pair(file_in,file_in1,file_in2,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    mm2 = {}
    
    for line in open(file_in2):
        vec = line.strip().split("\t")
        mm2.setdefault(vec[0],vec[1])
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[2],vec[0])
        mm1.setdefault(vec[2],vec[1])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec1 = vec[1].split(":")
        id1 = vec1[0]
        id2 = vec1[1]
        id = ""
        if(mm.has_key(id1)):
            id = id1
        elif(mm.has_key(id2)):
            id = id2
        if(id != ""):
            snp_site = mm[id]
            site = mm1[id]
            i1 = int(vec[0].split(":")[1])
            i2 = int(snp_site.split(":")[1])
            if(i2 - i1 <= 70):
                if(mm2.has_key(snp_site)):
                    gatk_snp = mm2[snp_site]
                    if(gatk_snp == vec[4]):
                        fout.write(line.strip() + "\t" + snp_site + "\t" + site + "\t" + id + "\n")
                    else:
                        print("no match\t" + snp_site)
                else:
                    print("no gatk\t" + snp_site)
            else:
                print("> 70\t"+snp_site)     
    return 0       

def hebing_52(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        a_b = ""
        i_a = ""
        for ii in range(1,len(vec)):
            ss = vec[ii]
            if(ss != id):
                if(ss.find(":") > 0):
                    a_b += ss + ","
                else:
                    i_a += ss + ","
        fout.write(id + "\t" + i_a + "\t" + a_b + "\n")
    return 0

def fx52_1(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[1][0:-1].split(",")
        vec2 = vec[2][0:-1].split(",")
        str_fx = ""
        i_succ = 0
        for ii in range(0,len(vec1)):
            i_total = int(vec1[ii])
            vec3 = vec2[ii].split(":")
            i_a = int(vec3[0])
            i_b = int(vec3[1])
            c_a = 0
            c_b = 0
            if(i_total >= 7):
                i_succ += 1 
                if(mm.has_key(ii)):
                    mm[ii] += 1
                else:
                    mm.setdefault(ii,1)
                if(i_a > 0 and i_b > 0):
                    str_fx += "AB,"
                    c_a += 1
                    c_b += 1
                elif(i_a > 0 and i_b == 0):
                    str_fx += "AA,"
                    c_a += 2
                elif(i_a == 0 and i_b > 0):
                    str_fx += "BB,"
                    c_b += 2
            else:
                str_fx += "N,"
        fout.write(line.strip() + "\t" + str_fx + "\t" + str(i_succ) + "\t" + str(c_a) + "\t" + str(c_b)+"\n")
        if(mm1.has_key(i_succ)):
            mm1[i_succ] += 1
        else:
            mm1.setdefault(i_succ,1)
    for key in mm:
        print(str(key) + "\t" + str(mm[key]))
    print("======")
    for key in mm1:
        print(str(key) + "\t" + str(mm1[key]))
    return 0

def kf_52(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in): 
        vec = line.strip().split("\t")
        i_total = int(vec[6])
        i_a = float(vec[9])
        i_b = float(vec[10])
        l_aa = i_a * i_a *i_total
        l_ab = 2 * i_a * i_b * i_total
        l_bb = i_b * i_b * i_total
        s_aa = 0
        s_ab = 0
        s_bb = 0
        vec1 = vec[5][0:-1].split(",")
        for ii in range(0,len(vec1)):
            si = vec1[ii]
            if(si == "AA"):
                s_aa += 1
            elif(si == "BB"):
                s_bb += 1
            elif(si == "AB"):
                s_ab += 1
        k_aa = (s_aa - l_aa)*(s_aa - l_aa) / l_aa
        k_ab = (s_ab - l_ab)*(s_ab - l_ab) / l_ab
        k_bb = (s_bb - l_bb)*(s_bb - l_bb) / l_bb
        k_total = k_aa + k_ab + k_bb
        #print(i_total,i_a,i_b,l_aa,l_ab,l_bb,s_aa,s_ab,s_bb,k_aa,k_ab,k_bb,k_total)
        so = str(i_total)+"\t" + str(i_a)+"\t" + str(i_b)+"\t" + str(l_aa)+"\t" + str(l_ab)+"\t" + str(l_bb)+"\t" + str(s_aa)+"\t" + str(s_ab)+"\t" + str(s_bb)+"\t" + str(k_aa)+"\t" + str(k_ab)+"\t" + str(k_bb)+"\t" + str(k_total)

        fout.write(line.strip() + "\t" + so + "\n")
    return 0

def kf_tj(file_in,file_cc):
    kf = float(file_cc)
    fout = open(file_cc+".out",'w')
    fout1 = open(file_cc+".log",'w')
    mm = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        i_total = float(vec[23])
        if(i_total < kf):
            fout.write(line.strip() + "\n")
            mean = round(float(vec[3]),1)
            if(mm.has_key(mean)):
                mm[mean] += 1
            else:
                mm.setdefault(mean,1)
    for key in mm:
        fout1.write(str(key)+"\t" + str(mm[key]) + "\n")
    return 0

def kf_tj1(file_in,file_out,file_out1):
    fout1 = open(file_out,'w')
    fout2 = open(file_out1,'w')
    for line in open(file_in): 
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(",")
        i_384 = 0
        i_664 = 0
        for ii in range(0,len(vec1)):
            vec2 = vec1[ii].split(":")
            i_a = float(vec2[0])
            i_b = float(vec2[1])
            if(i_a > 0 and i_b > 0):
                total = i_a + i_b
                da = 0.5
                db = 0.5
                
                ea = da * total
                eb = db * total
                
                oe_a = i_a - ea
                oe_b = i_b - eb
                
                double_oe_a = oe_a * oe_a
                double_oe_b = oe_b * oe_b
            
                r_a = double_oe_a / ea
                r_b = double_oe_b / eb
                
                total_re = r_a + r_b
                if(total_re >= 3.84):
                    i_384 +=1
                if(total_re >= 6.64):
                    i_664 +=1
        if(i_384 == 0):
            fout1.write(line.strip()+"\n")
        if(i_664 == 0):
            fout2.write(line.strip()+"\n")
    return 0

def mean_map(file_in):
    mm = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        mean = round(float(vec[3]),1)
        if(mm.has_key(mean)):
            mm[mean] += 1
        else:
            mm.setdefault(mean,1)
    for key in mm:
        print(str(key)+"\t" + str(mm[key]))
    return 0

def mess(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],line.strip())
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec1 = vec[0].split("~")
        id = vec1[0]
        if(mm.has_key(id)):
            fout.write(mm[id] + "\t" + line.strip() + "\n")
    return 0

def depth_static(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        mm.setdefault(id,int(vec[1]))
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec1 = vec[3].split(":")
        i1 = mm[vec1[0]]
        i2 = mm[vec1[1]]
        i_total = i1 + i2
        fout.write(vec[0] + "\t" + str(i1) + "\t" + str(i2) + "\t" + str(i_total) + "\n")
    return 0

def cx_pd(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm_p = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        i_total = 0
        ss = ""
        for ii in range(1,len(vec)):
            kk = vec[ii]
            print(kk)
            if(kk != id and kk != "0"):
                i_total += 1
            if(kk != id):
                ss += kk + "\t"
        mm.setdefault(id,i_total)
        mm_p.setdefault(id,ss)
    mm1 = {}
    mm1_p = {}
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        i_total = 0
        ss = ""
        for ii in range(1,len(vec)):
            kk = vec[ii]
            if(kk != id and kk != "0"):
                i_total += 1
            if(kk != id):
                ss += kk + "\t"
        mm1.setdefault(id,i_total)
        mm1_p.setdefault(id,ss)
    for key in mm:
        ik = mm1[key]
        ss1 = mm_p[key]
        ss2 = mm1_p[key]
        fout.write(key + "\t" + str(mm[key]) + "\t" + str(ik) + "\t" + ss1 + ss2 + "\n")
    return 0

def sex_1(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[0].split("_")
        i_m = int(vec1[1])
        i_f = int(vec1[2])
        i_c = int(vec[2])
        i_t = i_m  + i_f
        #if(i_c == -36 and i_t >= 40):
        i_kf = kf_2_1(i_f,i_m)
        fout.write(line.strip() + "\t" + str(i_kf) + "\t" + str(i_t)+ "\n")
    return 0

def sex_2(file_in):
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[3]
        print("echo '>" + vec[0]  + "'")
        print("awk '{print(substr($0,"+id+",60));}' ref.fa")
    return 0

def sex_3(file_in):
    for line in open(file_in):
        vec = line.strip().split("\t")
        reads = vec[1]
        ir = reads[0:4]
        print("grep '" + reads + "' " + ir+".sort");
    return 0

def sex_4(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    ii = 0
    for line in open(file_in):
        vec = line.strip().split("\t")
        reads = vec[1]
        id = vec[0]
        if(mm.has_key(reads)):
            mm[reads] += id
        else:
            mm.setdefault(reads,id)
    for key in mm:
        vv = mm[key].split(">")
        ss = len(vv) - 1
        ii += ss
        fout.write("0\t" + key + "\t"+str(ss)+"\t" + mm[key] + "\n")
    print(ii)
    return 0

def sex_5(file_in,file_out):
    fout = open(file_out,'w')
    ii = 1
    for line in open(file_in):
        if(ii == 1):
            fout.write(line.strip() + "_-\n")
        if(ii == 2):
            ik = dd_reverse(line.strip())
            fout.write(ik + "\n")
        elif(ii == 4):
            fout.write("+\n")
            fout.write(line.strip() + "\n")
            ii = 0
        ii += 1
    return 0

def qj_sex_mark(file_in,file_in1,file_gy,file_1,file_2):
    mm = {}
    fout_gy = open(file_gy,'w')
    fout_1 = open(file_1,'w')
    fout_2 = open(file_2,'w')
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0],vec1[1])
    igy = 1
    it1 = 1
    it2 = 1
    print("load over")
    for line in open(file_in1):
        vec1 = line.strip().split("\t")
        reads = vec1[0]
        if(mm.has_key(reads)):
            fout_gy.write(">" + str(igy) + "_" + str(mm[reads]) + "_" + str(vec1[1]) + "\n")
            fout_gy.write(reads + "\n")
            mm.pop(reads)
            igy += 1
        else:
            fout_2.write(">" + str(it2) + "_" + str(vec1[1]) + "\n")
            fout_2.write(reads + "\n")
            it2 += 1
    for key in mm:
        fout_1.write(">" + str(it1) + "_" + str(mm[key]) + "\n")
        fout_1.write(key + "\n")
        it1 += 1
    return 0


def pj_md(file_in,file_in1,file_1):
    mm = {}
    fout_1 = open(file_1,'w')
    for line in open(file_in):
        vec1 = line.strip().split("\t")
        mm.setdefault(vec1[0][1:],int(vec1[1]))
    mm1 = {}
    mm2 = {}
    ii = 1
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(ii > 1):
            id = vec[1]
            ist = id[0]
            site = vec[2]
            if(ist == "L"):
                if(mm1.has_key(id)):
                    mm1[id] += site + ":"
                else:
                    mm1.setdefault(id,site + ":")
            else:
                if(mm2.has_key(id)):
                    mm2[id] += site + ":"
                else:
                    mm2.setdefault(id,site + ":")
        ii += 1
    #fout_1.write("snp\tchr\tpos\n")
    for key in mm1:
        id = key.split("_")[1]
        veck1 = mm1[key][0:-1].split(":")
        xlen = len(veck1)
        mm4 = {}
        for kk in range(0,xlen):
            mm4.setdefault(kk,int(veck1[kk]))
        mm_temp = sorted(mm4.items(),key=lambda x:x[0],reverse=False)
        for kk in range(0,len(mm_temp)):
            fout_1.write(id + ":" + str(mm_temp[kk][1]) + "\t" + id + "\t" + str(mm_temp[kk][1]) + "\n")
    i_length = 0
    for key in mm2:
        vec1 = mm2[key][0:-1].split(":")
        lenkk = mm[key]
        mm5 = {}
        for kk in range(0,len(vec1)):
            mm5.setdefault(kk,int(vec1[kk]))
        mm_temp = sorted(mm5.items(),key=lambda x:x[0],reverse=False)
        for kk in range(0,len(mm_temp)):
            i_site = int(mm_temp[kk][1])
            i_len = i_length + i_site
            fout_1.write("mix:" + str(i_len) + "\tmix\t" + str(i_len) + "\n")
        i_length += lenkk
    
    return 0


def get_00(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    ilen = 0
    i_ref = ""
    for line in open(file_in):
        i_ref = line.strip()
        ilen = len(i_ref)
    ik = 1
    while(ik < ilen):
        mm.setdefault(ik,0)
        iw = i_ref[ik-1:ik + 59]
        mm1.setdefault(ik,iw)
        ik += 100

    for line in open(file_in1):
        ii = int(line.strip().split(":")[1])
        if(mm.has_key(ii)):
            mm[ii] += 1
    for key in mm:
        if(mm[key] == 0):
            siu = mm1[key]
            fout.write(str(key) + "\t" + siu + "\n")        
    print(len(mm))
    return 0

def get_11(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        i_site = int(vec[0])
        mm.setdefault(i_site,vec[1])
    print(len(mm))
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(line[0] != "@"):
            i_site = int(vec[3])
            for key in mm:
                i_begin = key - 10 
                i_end = key + 90
                if(i_begin <= i_site and i_end >= i_site):
                    mm1.setdefault(key,mm[key])
    print(len(mm1))
    for key in mm:
        if(not mm1.has_key(key)):
            fout.write(str(key) + "\t" + mm[key] + "\n")
            
def find_aaaa(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        mm1.setdefault(vec[0],1)
    for kk in mm:
        file_in2 = "/newlustre/home/liaolanjie/yc/cut_mix_3/sort/bl_m/" + kk + ".sort"
        for line in open(file_in2):
            vec = line.strip().split("\t")
            if(mm1.has_key(vec[1])):
                fout.write(line.strip() + "\n")
    return 0

def find_bb(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm1.setdefault(vec[0],1)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm1.has_key(vec[0])):
            fout.write(line.strip() + "\n")
    return 0


def zhan():
    i = 1
    while(i == 1):
        pass
    return 0


def indel_only(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0] + "_" + vec[1]
        indel = vec[3] + "_" + vec[4]
        mm.setdefault(id,indel)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0] + "_" + vec[1]
        indel = vec[3] + "_" + vec[4]
        if(mm.has_key(id)):
            fout.write("gy\t" + id + "\t" + mm[id] + "\t" + indel+"\n")
            mm.pop(id)
        else:
            fout.write("ty2\t" + id + "\t" + indel+"\n")
    for key in mm:
        fout.write("ty1\t" + key + "\t" + mm[key] +"\n")
    return 0

def k1(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        mm.setdefault(line.strip(),1)
    mm1 = {}
    for key in mm:
        ff = "D://down1//" + key
        for line in open(ff):
            vec = line.strip().split("\t")
            id = vec[0] + ":" + vec[1]
            snp = vec[3] + ":" + vec[4]
            if(mm1.has_key(id)):
                mm1[id] += snp + "~"
            else:
                mm1.setdefault(id,snp+"~")
    for key in mm1:
        vec = mm1[key][0:-1].split("~")
        mm2 = {}
        ik = ""
        for ii in range(0,len(vec)):
            mm2.setdefault(vec[ii],1)
            ik = vec[ii]
        if(len(vec) == 8):
            fout.write(key+"\t" +ik + "\n" )
    return 0

def k2(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0].split(":")[0]
        if(mm.has_key(site)):
            mm[site] += vec[1] + "~"
        else:
            mm.setdefault(site,vec[1]+ "~")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        site = vec[0].split(":")[0]
        if(mm.has_key(site)):
            mm.pop(site)
    for key in mm:
        fout.write(key+"\t"+mm[key] + "\n")
def dd_reverse_aa(str_in):
    trans_table = string.maketrans('ATCG','TAGC')
    a1 = str_in.translate(trans_table)
    str_r1 = a1[::-1]
    return str_r1

def tiqu_aa(file_in,file_in1,file_in2,file_out):
    fout = open(file_out,'w')
    mm = {}
    i = 1
    s1 = ""
    s2 = ""

    for line in open(file_in):
        line = line.strip()
        if(i == 1):
            s1 = line
        elif(i == 2):
            s2 = line
        elif(i == 4):
            mm.setdefault(s2,s1+"\t"+line)
            i = 0;
        i += 1
    i = 1
    s1 = ""
    s2 = ""
    for line in open(file_in1):
        line = line.strip()
        if(i == 1):
            s1 = line
        elif(i == 2):
            s2 = line
        elif(i == 4):
            mm.setdefault(s2,s1+"\t"+line)
            i = 0;
        i += 1
    for line in open(file_in2):
        vec = line.strip().split("\t")
        r1 = vec[0]
        r2 = dd_reverse_aa(r1)
        for key in mm:
            i1 = key.find(r1)
            i2 = key.find(r2)
            if(i1 >= 0 or i2 >= 0):
                vec1 = mm[key].split("\t")
                fout.write(vec1[0]+"\n")
                fout.write(key+"\n+\n")
                fout.write(vec1[1]+"\n")
    return 0

def kf_tt(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}

    for line in open(file_in):
        ik = line.strip()[0:4]
        if(mm.has_key(ik)):
            mm[ik] += 1
        else:
            mm.setdefault(ik,1)
    for ii in mm:
        fout.write(str(ii)+"\t" + str(mm[ii]) + "\n")
        
    return 0

def fx52_1_new(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[1][0:-1].split(",")
        vec2 = vec[2][0:-1].split(",")
        str_fx = ""
        i_succ = 0
        for ii in range(0,len(vec1)):
            i_total = int(vec1[ii])
            vec3 = vec2[ii].split(":")
            i_a = int(vec3[0])
            i_b = int(vec3[1])
            c_a = 0
            c_b = 0
            if(i_total >= 7 and i_a <= 23 and i_b <=23):
                i_succ += 1 
                if(mm.has_key(ii)):
                    mm[ii] += 1
                else:
                    mm.setdefault(ii,1)
                if(i_a > 0 and i_b > 0):
                    str_fx += "AB,"
                    c_a += 1
                    c_b += 1
                elif(i_a > 0 and i_b == 0):
                    str_fx += "AA,"
                    c_a += 2
                elif(i_a == 0 and i_b > 0):
                    str_fx += "BB,"
                    c_b += 2
            else:
                str_fx += "N,"
        fout.write(line.strip() + "\t" + str_fx + "\t" + str(i_succ) + "\t" + str(c_a) + "\t" + str(c_b)+"\n")
        if(mm1.has_key(i_succ)):
            mm1[i_succ] += 1
        else:
            mm1.setdefault(i_succ,1)
    for key in mm:
        print(str(key) + "\t" + str(mm[key]))
    print("======")
    for key in mm1:
        print(str(key) + "\t" + str(mm1[key]))
    return 0

def shai_kk(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in): 
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[1][0:-1].split(",")
        vec2 = vec[2][0:-1].split(",")
        is_write = 1
        for ii in range(0,len(vec1)):
            i_total = int(vec1[ii])
            vec3 = vec2[ii].split(":")
            i_a = int(vec3[0])
            i_b = int(vec3[1])
            c_a = 0
            c_b = 0
            if(i_a <= 2 or i_a >= 16 or i_b <=2 or i_b >= 16):
                is_write = 0
        if(is_write == 1):
            fout.write(line.strip() + "\n")
            
def fan1(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        if(line[0] != "@"):
            id = vec[2] + "@" + vec[3]
            reads = vec[9]
            mm.setdefault(id,reads)
    i_total = 0
    i_stat = 0
    i_scaff = ""
    for key in mm:
        vec = key.split("@")
        scaff = vec[0]
        stat = int(vec[1])
        if(i_scaff == scaff):
            i_dd = stat - i_stat 
            if(i_dd <= 150):
                i_total += i_dd
            else:
                i_total += 150
        else:
            i_scaff = scaff
        i_stat = stat
        fout.write(key+"\t" + mm[key] + "\n")
    print(i_total)
    return 0

def fan2(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in): 
        vec = line.strip().split("\t")
        id = vec[0]
        reads = vec[1]
        fout.write(">" + id + "_a\n" + reads + "\n")
        fout.write(">" + id + "_b\n" + reads + "\n")
        fout.write(">" + id + "_c\n" + reads + "\n")
        fout.write(">" + id + "_d\n" + reads + "\n")
        fout.write(">" + id + "_e\n" + reads + "\n")
    return 0

def fan3(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        if(line[0] != "@"):
            id = vec[2] + "@" + vec[3]
            id_reads = vec[0]
            if(mm.has_key(id)):
                mm[id] += id_reads +":"
            else:
                mm.setdefault(id,id_reads +":")
    for key in mm:
        fout.write(key+"\t" + mm[key] + "\n")
    return 0

def find_x10_2type(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        if(len(vec) > 6):
            id = vec[0]
            scaff = vec[2]
            site = vec[3]
            ii = scaff + ":" + site
            if(mm.has_key(ii)):
                mm[ii] += id + ":"
            else:
                mm.setdefault(ii,id+":")
    for key in mm:
        vec = mm[key][0:-1].split(":")
        i_site = len(vec)
        mmk = {}
        for ii in range(0,i_site):
            ik = vec[ii][0:-2]
            if(mmk.has_key(ik)):
                mmk[ik] += 1
            else:
                mmk.setdefault(ik,1)
        s_oo = ""
        for key1 in mmk:
            s_oo += key1 + "@" + str(mmk[key1]) + "\t"
        fout.write(key + "\t" + s_oo[0:-1] + "\n")
    return 0

def shai_x10_2type(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        biaozhun = 0
        ss_k = ""
        for ii in range(1,len(vec)):
            cc = vec[ii]
            vec1 = cc.split("@")
            id = vec1[0]
            num_cc = int(vec1[1])
            if(num_cc == 5):
                fout.write(id+"\t" + vec[0] + "\n")
                biaozhun += 1
                ss_k += id + "\t"
        if(len(vec) == 3 and biaozhun == 2):
            fout1.write(vec[0] + "\t" +  ss_k + "\n")
    return 0

def find_x5fa(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    ic = 1
    is_find = 0
    ip = ""
    for line in open(file_in1):
        line = line.strip()
        if(ic == 1):
            id = line[1:-2]
            if(mm.has_key(id)):
                is_find = 1
                ip = line
        else:
            ic = 0
            if(is_find == 1):
                fout.write(ip+"\n" + line+"\n")
            is_find = 0
            ip = ""
        ic += 1
    return 0

def find_x5fa_new(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    ic = 1
    is_find = 0
    ip = ""
    for line in open(file_in1):
        line = line.strip()
        if(ic == 1):
            id = line[1:]
            if(mm.has_key(id)):
                is_find = 1
                ip = line
        else:
            ic = 0
            if(is_find == 1):
                fout.write(ip+"_a\n" + line+"\n")
                fout.write(ip+"_b\n" + line+"\n")
                fout.write(ip+"_c\n" + line+"\n")
                fout.write(ip+"_d\n" + line+"\n")
                fout.write(ip+"_e\n" + line+"\n")
            is_find = 0
            ip = ""
        ic += 1
    return 0


def jz1(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in): 
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[1].split(":")
        i_same = 0
        for ii in range(0,len(vec1)):
            ik = vec1[ii]
            if(ik.find(id) >= 0):
                i_same += 1
        if(i_same == 5):
            fout.write(id+"\n")
        else:
            print(id)
    return 0

def jj(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in): 
        line = line.strip()
        ss = ""
        for ii in range(0,len(line)):
            sk = line[ii]
            if(sk == "A"):
                ss += "1"
            elif(sk == "T"):
                ss += "2"
            elif(sk == "C"):
                ss += "3"
            elif(sk == "G"):
                ss += "4"
            else:
                ss += "5"
        fout.write(ss+"\n")
    return 0


def in_gene(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        vec1 = vec[0].split(":")
        scaff = vec1[0]
        site = vec1[1].split("~")[0]
        if(mm.has_key(scaff)):
            mm[scaff] += site + ":"
        else:
            mm.setdefault(scaff,site+":")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        beg = int(vec[3])
        end = int(vec[4])
        vec1 = vec[8].split(";")
        geneid = vec1[0]
        if(mm.has_key(id)):
            vec2 = mm[id][0:-1].split(":")
            mm1 = {}
            for ii in range(0,len(vec2)):
                ik = int(vec2[ii])
                if(beg <= ik and end >= ik):
                    mm1.setdefault(ik,geneid)
            for key in mm1:
                fout.write(id + ":" +str(key)+"\t"+mm1[key]+"\n")
            
    return 0

def bl_male_female(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in): 
        vec = line.strip().split("\t")
        site = vec[0]
        indel = vec[1]
        if(mm.has_key(site)):
            mm[site] += indel + ":"
        else:
            mm.setdefault(site,indel+":")
    mm1 = {}
    print(len(mm))
    for key in mm:
        vec = mm[key][0:-1].split(":")
        if(len(vec) == 8):
            mm1.setdefault(key,mm[key][0:-1])
    print(len(mm1))
    for line in open(file_in1):
        vec = line.strip().split("\t")
        site = vec[0]
        indel = vec[1]
        if(mm1.has_key(site)):
            mm1.pop(site)
    for key in mm1:
        fout.write(key+"\t" + mm1[key] + "\n")
            
    return 0

def del_chai_128(file_name):
    in_file = "sort/all/" + file_name + ".sort"
    out_file = "del/" + file_name
    fout = open(out_file,'w')
    str_l = ""
    int_l = 0
    for line in open(in_file):
        line = line.strip()
        if(str_l == ""):
            str_l = line
        elif(str_l != line):
            fout.write(str_l +"\t" + str(int_l) + "\n")
            str_l = line
            int_l = 0
        int_l += 1
    fout.write(str_l +"\t" + str(int_l) + "\n")
    return 0
    
def del_1(file_name):
    in_file = "del/all/" + file_name
    in1_file = "all_50.1"
    out_file = "del/no_1/" + file_name
    fout = open(out_file,'w')
    mm = {}
    for line in open(in_file):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[1])
    print(len(mm))
    for line in open(in1_file):
        line = line.strip()
        if(mm.has_key(line)):
            mm.pop(line)
    print(len(mm))
    for key in mm:
        fout.write(key+"\t" + mm[key] + "\n")
    return 0

def is_aakk(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("_")
        mm.setdefault(vec[1],vec[0])
    for key in mm:
        ikk = mm[key] + "_" + key
        s1 = "samtools view  1_k2.sort.bam "+ mm[key] + ":" + key + "-" + key
        s2 = "samtools view  100_k2.sort.bam "+ mm[key] + ":" + key + "-" + key
        val1 = os.popen(s1)
        mm1 = {}
        for line in val1.readlines():
            vec1 = line.strip().split("\t")
            fec = vec1[11].split(":")[2]
            mm1.setdefault(int(fec),vec1[0])
        mk = sorted(mm1.items(),key=lambda x:x[0],reverse=False)
        value = list(mk)[0][0]
        val2 = os.popen(s2)
        mm2 = {}
        for line in val2.readlines():
            vec1 = line.strip().split("\t")
            fec = vec1[11].split(":")[2]
            mm2.setdefault(int(fec),vec1[0])
        mk1 = sorted(mm2.items(),key=lambda x:x[0],reverse=False)
        value1 = list(mk1)[0][0]
        fout.write(ikk + "\t" + str(value) + "\t" + str(value1) + "\t" + str(value - value1) + "\n")
    return 0

def tiqu_gat(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    ik = 1
    id = ""
    for line in open(file_in):
        line = line.strip()
        if(ik == 1):
            id = line[1:]
        else:
            ik = 0
            mm.setdefault(id,line)
        ik += 1
    for line in open(file_in1):
        vec = line.strip().split("\t")
        veci = vec[0].split(":")
        reads = mm[veci[0]]
        site = int(veci[1])
        i_begin = site - 40
        i_end = site + 40
        sr = reads[i_begin:i_end]
        fout.write(line.strip()+"\t" + sr + "\n")
    return 0

def tj_file_128_pv(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    mm1 = {}
    for key in mm:
        file_in1 = "del/all/" + key
        for line in open(file_in1):
            vec = line.strip().split("\t")
            cc = int(vec[1])
            if(mm1.has_key(cc)):
                mm1[cc] += 1
            else:
                mm1.setdefault(cc,1)
    for key in mm1:
        fout.write(str(cc) + "\t" + str(mm1[key]) + "\n")
    return 0

def only_pair(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip().split("\t")
        mm.setdefault(line,1)
    mm1 = {}
    for key in mm:
        file_in1 = "del/all/" + key
        for line in open(file_in1):
            vec = line.strip().split("\t")
            cc = int(vec[1])
            if(mm1.has_key(cc)):
                mm1[cc] += 1
            else:
                mm1.setdefault(cc,1)
    for key in mm1:
        fout.write(str(cc) + "\t" + str(mm1[key]) + "\n")
    return 0

def tiqu_site_kk(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    ii = 1
    reads = ""
    for line in open(file_in1):
        if(ii == 2):
            reads = line.strip()
        ii += 1
    for key in mm:
        ii = reads.find(key) +1
        fout.write(str(key) + "\t" + str(ii) + "\n")
    return 0

def tiqu_sam_kk(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[1],vec[0])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(len(vec) > 5):
            site = vec[3]
            if(mm.has_key(site)):
                fout.write(vec[2]+":" + vec[3]+"\t" + vec[0] + "\t" + vec[9] + "\t" + vec[11] + "\t" + mm[site]+"\n")
    return 0

def tiqu_sam_mess(file_in,file_in1,file_in2,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,"")
    ii = ""
    i1 = 1
    for line in open(file_in1):
        line = line.strip()
        if(i1 == 1):
            ii = line[1:]
        else:
            i1 = 0
            if(mm.has_key(ii)):
                mm[ii] = line
        i1 += 1   
    for line in open(file_in2):
        vec = line.strip().split("\t")
        id = vec[2]
        qid = vec[0][0:-2]
        if(mm.has_key(id) and id != qid):
            fout.write(vec[2] + "\t" + vec[0] + "\t" + vec[9] + "\t" + vec[11] + "\t" + mm[id]+"\n")
    return 0

def pd_link_snp(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[1]
        if(mm.has_key(id)):
            mm[id] += site + ":"
        else:
            mm.setdefault(id,site + ":")
    for key in mm:
        yi = mm[key]
        vec = yi[0:-1].split(":")
        mk1 = sorted(vec,key=lambda x:x[0],reverse=False)
        i_old = 0
        is_find = 0
        i_k = 0
        i_max = 0
        s_o9 = ""
        for ii in range(0,len(mk1)):
            i_new = int(mk1[ii])
            s_o9 += mk1[ii] + "_"
            if(i_old == 0):
                i_old = i_new
            elif(i_new - i_old == 1):
                is_find = 1
                i_k += 1
                if(i_k >= i_max):
                    i_max = i_k
            else:
                i_k = 0
            i_old = i_new
        if(is_find == 1):
            fout.write(key+"\t" + s_o9[0:-1] + "\t" + str(i_max +1) + "\n")
            print(key,mk1,i_max+1)
    return 0

def tiqu(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        reads1 = vec[1]
        id2 = vec[2][0:-2]
        reads2 = vec[3]
        mm.setdefault(id,reads1)
        mm.setdefault(id2,reads2)
        mm1.setdefault(id,id2)
    for key in mm1:
        r1 = mm[key]
        r2 = mm[mm1[key]]
        fout.write(r1+"\t"+key+"\t"+mm1[key]+"\n")
        fout.write(r2+"\t"+mm1[key] + "\t"+key+"\n")

def tiqu_sam_mess(file_in,file_in1,file_in2,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,"")
    ii = ""
    i1 = 1
    for line in open(file_in1):
        line = line.strip()
        if(i1 == 1):
            ii = line[1:]
        else:
            i1 = 0
            if(mm.has_key(ii)):
                mm[ii] = line
        i1 += 1   
    for line in open(file_in2):
        vec = line.strip().split("\t")
        id = vec[2]
        if(mm.has_key(id)):
            refii = mm.pop(id)
            fout.write(vec[2] + "\t" + refii + "\t" + vec[0][0:-2] + "\t" + vec[9] + "\t" + vec[11] + "\n")
    return 0

def pd_link_snp_pp(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[1]
        if(mm.has_key(id)):
            mm[id] += site + ":"
        else:
            mm.setdefault(id,site + ":")
    for key in mm:
        yi = mm[key]
        vec = yi[0:-1].split(":")
        mk1 = sorted(vec,key=lambda x:x[0],reverse=False)
        i_old = 0
        is_find = 0
        i_k = 0
        i_max = ""
        s_o9 = ""
        for ii in range(0,len(mk1)):
            i_new = int(mk1[ii])
            s_o9 += mk1[ii] + "_"
            if(i_old == 0):
                i_old = i_new
                i_max = i_old + ":"
            elif(i_new - i_old == 1):
                is_find = 1
                i_k += 1
                i_max += i_new + ":"
            else:
                if(is_find == 1):
                    print(key)
                i_max =""
                is_find = 0
                i_k = 0
            i_old = i_new
    return 0

def sex_mark_pair_1(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        st1 = ""
        for ii in range(1,len(vec)):
            ik = vec[ii]
            if(str.isdigit(ik)):
                st1 += ik +":"
        vec1 = st1[0:-1].split(":")
        st_male = ""
        st_female = ""
        i_male = 0
        i_female = 0
        for ii in range(0,8):
            ik = vec1[ii]
            if(ik != "0"):
                i_male += 1
            st_male += ik + ":"
        for ii in range(0,len(vec1)):
            ik = vec1[ii]
            if(ik != "0"):
                i_female += 1
            st_female += ik + ":"
        fout.write(vec[0]+"\t" + str(i_male) + "\t" + str(i_female) + "\t" + st_male + "\t" + st_female + "\n")      
    return 0

def sex_mark_pair_2(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm_male = {}
    mm_female = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm_male.setdefault(vec[0],int(vec[1]))
        mm_female.setdefault(vec[0],int(vec[2]))   
    i_line = 1
    st1 = ""
    for line in open(file_in1):
        vec = line.strip().split("\t")
        sss = vec[0]
        if(i_line == 1):
            st1 = sss
        else:
            male_1 = mm_male[st1]
            female_1 = mm_female[sss]
            fout.write(st1 + "\t" + str(male_1) + "\t" + sss + "\t" + str(female_1) + "\n")
            i_line = 0
        i_line += 1
    return 0

def sex_mark_pair_3(file_in,file_out):
    fout = open(file_out,'w')
    ik = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[3][0:-1].split(":")
        for ii in range(0,len(vec1)):
            reads = vec1[ii]
            r_reads = dd_reverse(reads)
            fout.write(">z_"+str(ik) + "\n" + reads+"\n")
            fout.write(">f_"+str(ik) + "\n" + r_reads+"\n")
        ik += 1
    return 0

def reads_x5(file_in,file_out,f_name):
    fout = open(file_out,'w')
    ik = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
        fout.write(">"+str(ik) + "_"+f_name + "_" +vec[1] +"_a\n"+vec[0]+"\n")
        fout.write(">"+str(ik) + "_"+f_name + "_" +vec[1] +"_b\n"+vec[0]+"\n")
        fout.write(">"+str(ik) + "_"+f_name + "_" +vec[1] +"_c\n"+vec[0]+"\n")
        fout.write(">"+str(ik) + "_"+f_name + "_" +vec[1] +"_d\n"+vec[0]+"\n")
        fout.write(">"+str(ik) + "_"+f_name + "_" +vec[1] +"_e\n"+vec[0]+"\n")
        ik += 1
    return 0

def pd_gy(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            fout.write(line.strip() + "\n")
        else:
            fout1.write(line.strip() + "\n")
    return 0

def pd_gy_zf(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[2])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        kk = mm.pop(vec[0])
        vv = vec[2][0:-1].split(":")
        vv1 = kk[0:-1].split(":")
        mmk = {}
        for ii in range(0,len(vv)):
            mmk.setdefault(vv[ii],1)
        for ii in range(0,len(vv1)):
            ik = vv1[ii]
            ik_r = dd_reverse(ik)
            if(mmk.has_key(ik_r)):
                pass
            else:
                mmk.setdefault(ik,2)
        ss = ""
        for key in mmk:
            ss += key + ":"
        vv2 = ss[0:-1].split(":")
        fout.write(vec[0] + "\t" + str(len(vv2)) + "\t" + ss + "\n")
    return 0

def find_x2_reads(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    i_tt = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        id = vec[0]
        if(vec[1] == "2"):
            for ii in range(0,len(vec1)):
                reads = vec1[ii]
                fout1.write(reads + "\t" + id + "\n")
                fout.write(">" + str(i_tt) +  "_a\n" + reads + "\n")
                fout.write(">" + str(i_tt) +  "_b\n" + reads + "\n")
                fout.write(">" + str(i_tt) +  "_c\n" + reads + "\n")
                fout.write(">" + str(i_tt) +  "_d\n" + reads + "\n")
                fout.write(">" + str(i_tt) +  "_e\n" + reads + "\n")
            i_tt += 1
    return 0

def find_1178_bl_gap(file_in,file_out):
    fout = open(file_out,'w')
    id = ""
    old_site = 0
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        r_id = vec[2]
        r_site = int(vec[3])
        if(id == ""):
            id = r_id
            old_site = r_site
        elif(id == r_id):
            dd = r_site - old_site
            if(dd > 60):
                print(r_site,old_site,dd)
                ik = dd/10
                if(mm.has_key(ik)):
                    mm[ik] += id + ":"
                else:
                    mm.setdefault(ik,id+":") 
        id = r_id
        old_site = r_site
    for key in mm:
        vec = mm[key][0:-1].split(":")
        fout.write(str(key*10)+"\t" + str(len(vec)) + "\t" + mm[key] + "\n")
    return 0

def xyjj_c_t(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    for line in open(file_in1):
        line = line.strip()
        if(mm.has_key(line)):
            mm[line] = 2
        else:
            print(line)
    for key in mm:
        fout.write(key + "\t" + str(mm[key]) + "\n")
    #python c_t.py c.id t.id c.tj > t.ty
    return 0 


def kk_pp(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split(" ")
        fout.write(vec[0] + "\t" + vec[1] + "\n")
    return 0 

def qj_iii(file_in,file_in1):
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        for ii in range(0,len(vec1)):
            mm.setdefault(vec1[ii],1)
    ii = 0
    ik = 0
    for line in open(file_in1):
        line = line.strip()
        if(mm.has_key(line)):
            print(line)
            ii += 1
        ik += 1
    print(ii,ik)
    return 0 


def find_x234_reads_chai(file_in,file_out,file_out2,file_out3,file_out4,is_z):
    fout = open(file_out,'w')
    fout2 = open(file_out2,'w')
    fout3 = open(file_out3,'w')
    fout4 = open(file_out4,'w')
    i_tt = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        id = vec[0]
        if(vec[1] == "2" or vec[1] == "3" or vec[1] == "4"):
            for ii in range(0,len(vec1)):
                reads = vec1[ii]
                id_kk = ">" + str(i_tt) +  "_"+vec[1] + "_" + is_z
                fout.write(id_kk + "\t" + reads + "\t" + id  + "\t" + vec[1] +"\n")
                if(vec[1] == "2"):
                    fout2.write(id_kk + "_a\n" + reads + "\n")
                    fout2.write(id_kk + "_b\n" + reads + "\n")
                    fout2.write(id_kk + "_c\n" + reads + "\n")
                    fout2.write(id_kk + "_d\n" + reads + "\n")
                    fout2.write(id_kk + "_e\n" + reads + "\n")
                elif(vec[1] == "3"):
                    fout3.write(id_kk + "_a\n" + reads + "\n")
                    fout3.write(id_kk + "_b\n" + reads + "\n")
                    fout3.write(id_kk + "_c\n" + reads + "\n")
                    fout3.write(id_kk + "_d\n" + reads + "\n")
                    fout3.write(id_kk + "_e\n" + reads + "\n")
                elif(vec[1] == "4"):
                    fout4.write(id_kk + "_a\n" + reads + "\n")
                    fout4.write(id_kk + "_b\n" + reads + "\n")
                    fout4.write(id_kk + "_c\n" + reads + "\n")
                    fout4.write(id_kk + "_d\n" + reads + "\n")
                    fout4.write(id_kk + "_e\n" + reads + "\n")
                i_tt += 1
    return 0

def hu_bu(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        sp = ""
        for ii in range(0,len(vec1)):
            f_re = vec1[ii]
            z_re = dd_reverse(f_re)
            sp += z_re + ":"
        ss =vec[0] + "\t" + vec[1]  +"\t" + sp;
        fout.write(ss+"\n")
    return 0 

def kk_pp_ll(file_in,file_out):
    fout = open(file_out,'w')
    i_tt = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        nn = vec[1]
        vec1 = vec[2][0:-1].split(":")
        for ii in range(0,len(vec1)):
            reads = vec1[ii]
            id_kk = ">" + id + "_" + nn + "_"  + str(i_tt)
            fout.write(id_kk + "a\n" + reads + "\n")
            fout.write(id_kk + "b\n" + reads + "\n")
            fout.write(id_kk + "c\n" + reads + "\n")
            fout.write(id_kk + "d\n" + reads + "\n")
            fout.write(id_kk + "e\n" + reads + "\n")
            i_tt += 1
    return 0 

def mess1(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        reads1 = vec[1]
        reads2 = vec[3]
        mm2.setdefault(id,reads1+":"+reads2)
        mm1.setdefault(reads1,0)
        mm1.setdefault(reads2,0)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        reads = vec[0]
        i_cc = int(vec[1])
        mm1[reads] = i_cc
    for key in mm2:
        vec = mm2[key].split(":")
        reads1 = vec[0]
        reads2 = vec[1]
        cc1 = mm1[reads1]
        cc2 = mm1[reads2]
        gene_type = "N"
        if(cc1 > 0):
            if(cc2 > 0):
                gene_type = "AB"
            else:
                gene_type = "AA"
        elif(cc2 > 0):
            gene_type = "BB"
        sk = str(cc1) + ":" + str(cc2)
        fout.write(key + "\t" + gene_type + "\t" + sk + "\t" + mm2[key] + "\n")
    return 0 

def mess2(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,0)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        mk1 = {}
        mk1.setdefault("AB",0)
        mk1.setdefault("BB",0)
        mk1.setdefault("N",0)
        mk1.setdefault("AA",0)
        mk2 = {}
        mk2.setdefault("AB",0)
        mk2.setdefault("BB",0)
        mk2.setdefault("N",0)
        mk2.setdefault("AA",0)
        is_male = 1
        if(mm.has_key(id)):
            for ii in range(1,len(vec)):
                ssv = vec[ii]
                if(ssv == id):
                    is_male = 0
                    continue
                if(is_male == 1):
                    mk1[ssv] += 1
                else:
                    mk2[ssv] += 1
            fout.write(id+"\tAB\t"+str(mk1["AB"]) + "\tBB\t"+ str(mk1["BB"]) +  "\tN\t" + str(mk1["N"]) + "\tAB\t"+str(mk2["AB"]) + "\tBB\t"+ str(mk2["BB"]) +  "\tN\t" + str(mk2["N"]) + "\n")
    return 0

def mess3(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    i_total = 0
    i_err = 0
    i_ture = 0
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        is_find = 0
        for ii in range(0,len(vec1)):
            rr = vec1[ii]
            if(mm.has_key(rr)):
                is_find = 1
        if(is_find == 0):
            fout.write(line.strip() + "\n")
            i_ture += 1
        else:
            i_err += 1
        i_total += 1
    print(i_total,i_err,i_ture)
    return 0 

def mess3_jhs(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],1)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[1])):
            fout.write(vec[0] + "\n")
        elif(mm.has_key(vec[2])):
            fout1.write(vec[0] + "\n")
    return 0 

def mess4(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        if(mm.has_key(id)):
            vec1 = vec[3][0:-1].split(":")
            for ii in range(0,len(vec1)):
                fout.write(vec1[ii] + "\n")
    return 0 

def snp_1(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        vec1 = vec[2][0:-1].split(":")
        for ii in range(0,len(vec1)):
            mm.setdefault(vec1[ii],"N")
        mm1.setdefault(id,vec[2][0:-1])
    print("load over")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(len(vec) > 13):
            reads = vec[9]
            site = vec[2] + ":" + vec[3]
            if(mm.has_key(reads)):
                mm[reads]  = site
    print("write begin")
    for key in mm1:
        vec = mm1[key].split(":")
        ss = ""
        mmk = {}
        is_find  = 1
        for ii in range(0,len(vec)):
            reads = vec[ii]
            site = mm[reads]
            if(site == "N"):
                is_find = 0 #鏈変竴涓猺eads鍦ㄥ熀鍥犵粍涓婃病鍖归厤鍒�
            mmk.setdefault(site,1)
        if(len(mmk) > 1):
            is_find = 2 #涓嶅悓reads锛屽尮閰嶅埌浜嗗涓綅缃�
        for kk in mmk:
            ss += mmk[kk] + ":"
        ss = ss[0:-1]
        if(ss != key):
            is_find = 3 #鍖归厤鐨勪綅缃拰鍘熷厛5鍊嶇殑鍘熷浣嶇疆涓嶇浉鍚岋紝杩欐槸灏嗚鐩栨帀is_find=2鐨勫��
        fout.write(key + "\t" +mm1[key] + "\t" + ss + "\t" +str(is_find) + "\n")
    return 0 


def snp_2_aa(file_in,file_out,file_out1,file_outerr):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fouterr = open(file_outerr,'w')
    i_err = 0
    i_total = 0
    for line in open(file_in):
        vec = line.strip().split("\t")
        i_total += 1
        if(vec[3] != "1"):
            fouterr.write(line.strip() + "\n")
            i_err += 1
            continue
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
                fout1.write(contig +":" + str(snp_site) + "\t" + vec[0] + "\t" + str_snp  + "\t" + str(len(mm_temp))+ "\n")
                ssk += str(snp_site) + ":" + str(key) + ":" + str_snp + ":" + str(len(mm_temp)) + "#"
        fout.write(line.strip() + "\t" + str(i_snp) + "\t" + ssk + "\n")
    print(file_in,i_err,i_total)
    return 0
        
def snp_del_re(file_in,file_out,file_outerr):
    fout = open(file_out,'w')
    fouterr = open(file_outerr,'w')
    i_err2 = 0
    i_err = 0
    i_total = 0
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0]
        id = vec[1]
        snp = vec[2]
        if(vec[3] != "2"):
            fouterr.write(line.strip() + "\tmanysnp\n")
            i_err2 += 1
        else:
            if(mm.has_key(site)):
                mm[site] += snp + ":"
                mm1[site] += id + ":"
            else:
                mm.setdefault(site,snp+":")
                mm1.setdefault(site,id+":")
    for key in mm:
        i_total += 1
        vec = mm[key][0:-1].split(":")
        mm_temp = {}
        str_snp = ""
        snp_only = ""
        for kk in range(0,len(vec)):
            mm_temp.setdefault(vec[kk])
            str_snp += vec[kk] + ":"
            snp_only = vec[kk]
        str_snp = str_snp[0:-1]
        if(len(mm_temp) > 1):
            i_err += 1
            fouterr.write(key+"\t"+str_snp+"\t" + mm1[key] + "\tsnp_no_map\n")
        else:
            fout.write(key+"\t"+snp_only+"\t" + mm1[key] + "\n")
    print(i_err,i_err2,i_total)
    return 0

def snp_os(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0].split(":")[1]
        val = os.popen("samtools view 2_snp.sort.bam " + vec[0] + "-" + site)
        mm = {}
        for line in val.readlines():
            vec1 = line.strip().split("\t")
            s0 = vec1[3]
            sr = vec1[9]
            if(mm.has_key(s0)):
                mm[s0] += sr+":"
            else:
                mm.setdefault(s0,sr+":")
        s0_out = ""
        s1_out = ""
        for key in mm:
            s0_out += key + ":"
            s1_out += mm[key] + "#"
        fout.write(vec[0] + "\t" + str(len(mm)) + "\t" + s0_out + "\t" + s1_out + "\n")
    return 0
            
def reads_pv(file_name):
    file_in = file_name
    file_out = file_name + ".tj"
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        icc = int(vec[1])
        if(mm.has_key(icc)):
            mm[icc] += 1
        else:
            mm.setdefault(icc,1)
    for key in mm:
        fout.write(str(key)+"\t"+str(mm[key]) + "\n")
    return 0
            
def snp_234_pv(file_in,file_out):
    mm = {}
    fout = open(file_out,'w')
    for line in open(file_in):
        line = line.strip()
        if(mm.has_key(line)):
            mm[line] += 1
        else:
            mm.setdefault(line,1)
    for key in mm:
        fout.write(key+"\t" + str(mm[key]) + "\n")
    return 0

def snp_reads_pv(file_in,file_out):
    mm = {}
    mm1 = {}
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        i_cc = int(vec[1])
        if(mm1.has_key(i_cc)):
            mm1[i_cc] += 1
        else:
            mm1.setdefault(i_cc,1)
        if(i_cc > 0):
            vec1 = vec[3][0:-2].split(":")
            for ii in range(0,len(vec1)):
                read = vec1[ii]
                if(read[0] == "#"):
                    read = read[1:]
                if(mm.has_key(read)):
                    mm[read] += 1
                else:
                    mm.setdefault(read,1)
    for key in mm:
        fout.write(key+"\t" + str(mm[key]) + "\n")
    for key in mm1:
        print(key,mm1[key])
    return 0

def snp_del_reads(file_in,file_in1,file_out):
    mm = {}
    fout = open(file_out,'w')
    for line in open(file_in1):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],1)
    for line in open(file_in):
        vec = line.strip().split("\t")
        i_cc = int(vec[1])
        is_write = 1
        if(i_cc > 0):
            vec1 = vec[3][0:-2].split(":")
            for ii in range(0,len(vec1)):
                read = vec1[ii]
                if(read[0] == "#"):
                    read = read[1:]
                if(mm.has_key(read)):
                    is_write = 0
            vec1 = vec[3][0:-1].split("#")
            is_ma = 0
            for ii in range(0,len(vec1)):
                vec2 = vec1[ii][0:-1].split(":")
                if(len(vec2) == 2):
                    is_ma = 1
        if(is_write == 1 and is_ma == 1):
            fout.write(line.strip() + "\n")
    return 0

def snp_del_cl(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[3][0:-1].split("#")
        is_ma = 0
        for ii in range(0,len(vec1)):
            vec2 = vec1[ii][0:-1].split(":")
            if(len(vec2) == 2):
                is_ma = 1
        if(is_ma == 1):
            fout.write(line.strip() + "\n")
        else:
            print(line.strip())
    return 0

def get_234_reads(file_in,file_out2,file_out3,file_out4):
    fout2 = open(file_out2,'w')
    fout3 = open(file_out3,'w')
    fout4 = open(file_out4,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        id = vec[0]
        if(vec[1] == "2" or vec[1] == "3" or vec[1] == "4"):
            i_kk = 1
            for ii in range(0,len(vec1)):
                reads = vec1[ii]
                id_kk = ">" + id+  ":"+vec[1]+"_"+str(i_kk)
                if(vec[1] == "2"):
                    fout2.write(id_kk + "\n" + reads + "\n")
                elif(vec[1] == "3"):
                    fout3.write(id_kk + "\n" + reads + "\n")
                elif(vec[1] == "4"):
                    fout4.write(id_kk + "\n" + reads + "\n")
                i_kk += 1
    return 0

def yz_234_reads(file_in,file_out,xx):
    i_x = int(xx)
    fout = open(file_out,'w')
    mm = {}
    mmk = {}
    for line in open(file_in):
        if(line[0] != "@"):
            vec = line.strip().split("\t")
            id = vec[0]
            site = vec[2] + ":" + vec[3]
            reads = vec[9]
            if(mm.has_key(site)):
                mm[site] += id + "#"
                mmk[site] += reads+":"
            else:
                mm.setdefault(site,id+"#")
                mmk.setdefault(site,reads+":")
    print(len(mm))
    for key in mm:
        vec = mm[key][0:-1].split("#")
        is_err = 0
        mm1 = {}
        for ii in range(0,len(vec)):
            str_temp = vec[ii]
            vec1 = str_temp.split(":")
            site = vec1[0] + ":" + vec1[1]
            mm1.setdefault(str_temp,1)
            if(site != key):
                is_err = 1
        if(len(mm1) != i_x):
            is_err = 2
        reads = mmk[key]
        fout.write(key + "\t" + mm[key] + "\t" + str(is_err) +"\t"+str(xx)+"\t" + reads +"\n")
    return 0

def tiqu_234_reads(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        if(vec[2] == "0"):
            mm.setdefault(vec[0],1)
            fout1.write(vec[0] + "\t" + vec[3]+"\t"+vec[4]+"\n")
    print(len(mm))
    for line in open(file_in1):
        if(line[0] != "@"):
            vec = line.strip().split("\t")
            vec1 = vec[0].split(":")
            site = vec1[0] + ":" + vec1[1]
            if(mm.has_key(site)):
                fout.write(line.strip()+"\n")
    return 0

def tiqu_all_reads(file_in,file_out):
    fout = open(file_out,'w')
    icc = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        for ii in range(0,len(vec1)):
            fout.write(">" + str(icc) + "\n" + vec1[ii] + "\n")
            icc +=1
        
    return 0

def get_base(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[0].split(":")
        i_site = int(vec1[1])
        e_site = vec[0]+"-"+str(i_site+1)
        mm.setdefault(vec[0],e_site)
        
    for key in mm:
        vv = mm[key]
        s_k = "samtools mpileup -A -B -C 0 -q 0 -Q 0 -t DP,AD --region "+ vv+" -f v13.fa 12_29/all.sort.bam"
        val = os.popen(s_k)
        for line in val.readlines():
            fout.write(key+"\t"+line)
    return 0

def is_only(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        site = vec[2] + ":" + vec[3]
        fen = ""
        indel = ""
        
        for ii in range(4,len(vec)):
            ui = vec[ii]
            #print(ui,ui.find("AS:i:"))
            veckk = ui.split(":")
            if(ui.find("AS:i:") >= 0):
                fen = veckk[2]
            elif(ui.find("XO:i:") >= 0):
                indel = veckk[2]
        sk = site+"@"+fen+"@"+ indel
        if(mm.has_key(id)):
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
                mm1[fen] += 1
            else:
                mm1.setdefault(fen,1)
        max_key = max(mm1) #鎵惧埌鏈�灏忓緱鍒�
        max_value = mm1[max_key]
        if(max_value == 1): #鏄惁鏈�灏忓緱鍒嗗�兼槸鍞竴鐨�
            for ii in range(0,len(vec)): #鎻愬彇鏈�灏忓緱鍒嗚褰�
                vec1 = vec[ii].split("@")
                fen = int(vec1[1])
                if(fen == max_key):
                    if(vec1[2] =="0"): #鍒ゆ柇鏈�灏忓緱鍒嗘槸鍚︿笉鏄痠ndel
                        fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\n")
                    else:
                        fout1.write(key+"\t"+mm[key]+"\tsmall_indel\n")
                
        else:
            fout1.write(key+"\t"+mm[key]+"\tmany_small_value\n")
    return 0

def is_pairaa(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        ref_site = vec[1]
        if(mm.has_key(ref_site)):
            mm[ref_site] += vec[0] + "#"
        else:
            mm.setdefault(ref_site,vec[0] + "#")
    for key in mm:
        vec = mm[key][0:-1].split("#")
        if(len(vec) != 2):
            print(key+"\t"+mm[key]+"\t1")
        else:
            id1 = vec[0]
            id2 = vec[1]
            if(id1.find(key) < 0 or id2.find(key) < 0):
                print(key+"\t"+mm[key]+"\t2")
            else:
                fout.write(key+"\t"+id1 + "\t" + id2+"\n")
    return 0

def find_pair(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[1],1)
        mm.setdefault(vec[2],1)
    
    ik = ""
    
    for line in open(file_in1):
        line = line.strip()
        if(line[0] == ">"):
            ik = line[1:]
        else:
            if(mm.has_key(ik)):
                fout.write(">" + ik + "a\n" + line + "\n")
                fout.write(">" + ik + "b\n" + line + "\n")
                fout.write(">" + ik + "c\n" + line + "\n")
                fout.write(">" + ik + "d\n" + line + "\n")
                fout.write(">" + ik + "e\n" + line + "\n")
    return 0 

def find_pair_kk(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[1],"")
        mm.setdefault(vec[2],"")
        mm1.setdefault(vec[0],vec[1] +"#" +vec[2])
    ik = ""
    
    for line in open(file_in1):
        line = line.strip()
        if(line[0] == ">"):
            ik = line[1:]
        else:
            if(mm.has_key(ik)):
                mm[ik] = line
    for key in mm1:
        vec = mm1[key].split("#")
        read1 = mm[vec[0]]
        read2 = mm[vec[1]]
        fout.write(key + "\t2\t"+read1+":"+read2+":\t"+vec[0]+":"+vec[1]+":\n")
    
    return 0 

def snp_2(file_in,file_out,file_out1,file_outerr):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fouterr = open(file_outerr,'w')
    i_err = 0
    i_total = 0
    for line in open(file_in):
        vec = line.strip().split("\t")
        i_total += 1
        if(vec[3] != "1"):
            fouterr.write(line.strip() + "\n")
            i_err += 1
            continue
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
                fout1.write(contig +":" + str(snp_site) + "\t" + vec[0] + "\t" + str_snp  + "\t" + str(len(mm_temp))+ "\n")
                ssk += str(snp_site) + ":" + str(key) + ":" + str_snp + ":" + str(len(mm_temp)) + "#"
        fout.write(line.strip() + "\t" + str(i_snp) + "\t" + ssk + "\n")
    print(file_in,i_err,i_total)
    return 0

def get_snp_2(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    i_err = 0
    i_total = 0
    for line in open(file_in):
        vec = line.strip().split("\t")
        i_total += 1
        vec1 = vec[2][0:-1].split(":")
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
                fout1.write(contig +":" + str(snp_site) + "\t" + vec[0] + "\t" + str_snp  + "\t" + str(len(mm_temp))+ "\n")
                ssk += str(snp_site) + ":" + str(key) + ":" + str_snp + ":" + str(len(mm_temp)) + "#"
        fout.write(line.strip() + "\t" + str(i_snp) + "\t" + ssk + "\n")
    print(file_in,i_err,i_total)
    return 0

def is_match_snp(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0]
        snp = vec[2]
        if(mm.has_key(site)):
            mm[site] += snp + ":"
        else:
            mm.setdefault(site,snp+":")
    print(len(mm))
    for key in mm:
        mm1 = {}
        vec = mm[key][0:-1].split(":")
        snp = ""
        for ii in range(0,len(vec)):
            mm1.setdefault(vec[ii],1)
            snp = vec[ii]
        if(len(mm1) == 1):
            fout.write(key+"\t"+str(len(vec))+"\t"+snp+"\n")
        else:
            fout1.write(key+"\t"+str(len(vec))+"\n")
    return 0

def create_fa_kk(file_in,file_out):
    fout = open(file_out,'w')
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        vec2 = vec[3][0:-1].split(":")
        id1 = vec2[0] + ":" + vec2[1] + ":" + vec2[2]
        id2 = vec2[3] + ":" + vec2[4] + ":" + vec2[5]
        fout.write(">"+id1+"\n"+vec1[0]+"\n")
        fout.write(">"+id2+"\n"+vec1[1]+"\n")
    return 0

def del_reduance(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0]
        snp = vec[2]
        if(mm.has_key(site)):
            mm[site] += snp + ":"
        else:
            mm.setdefault(site,snp+":")
    print(len(mm))
    for key in mm:
        mm1 = {}
        vec = mm[key][0:-1].split(":")
        snp = ""
        for ii in range(0,len(vec)):
            mm1.setdefault(vec[ii],1)
            snp = vec[ii]
        fout.write(key+"\t"+str(len(vec))+"\n")
    return 0

def del_one(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec_id = vec[0].split(":")
        contig = vec_id[0]
        site = int(vec_id[1])
        site_1 = site + 1
        s1 = contig + ":"+str(site)
        mm.setdefault(s1,"")
    
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec_id = vec[0].split(":")
        contig = vec_id[0]
        site = int(vec_id[1])
        site_1 = site + 1
        
        s1 = contig + ":"+str(site)
        s2 = contig + ":"+str(site_1)
        
        is_find = 0
        if(mm.has_key(s1)):
            is_find = 1
        elif(mm.has_key(s2)):
            is_find = 1
        if(is_find == 1):
            #print(line.strip())
            vec1 = vec[5][0:-1].split("#")
            i_find = 0 
            for ii in range(0,len(vec1)):
                i_site = vec1[ii].split(":")[0]
                if(i_site == site or i_site == site_1):
                    i_find += 1
                    print(line)
            if(i_find == 2):
                fout.write(line.strip()+"\n")
            else:
                fout1.write(line.strip()+"\n")
        else:
            fout.write(line.strip()+"\n")
    return 0

def del_out_id(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[1])):
            pass
        else:
            fout.write(line.strip()+"\n")
    return 0


def del_sam_id(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec2 = vec[3][0:-1].split(":")
        id1 = vec2[0] + ":" + vec2[1] + ":" + vec2[2]
        id2 = vec2[3] + ":" + vec2[4] + ":" + vec2[5]
        mm.setdefault(id1,"")
        mm.setdefault(id2,"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            pass
        else:
            fout.write(line.strip()+"\n")
    return 0

def del_in_site(file_in,file_in1,file_out,file_out1,file_out2,file_out3):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    fout3 = open(file_out3,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec2 = vec[1].split(":")
        i_dep = int(vec2[0])
        i_snp = int(vec2[1])
        if(i_dep > 0):
            vec1 = vec[2][0:-1].split("@")
            for ii in range(0,len(vec1)):
                mm.setdefault(vec1[ii],"")
        if(i_dep < i_snp):
            print(line.strip())
    print(len(mm))
    for line in open(file_in1):
        vec = line.strip().split("\t")
        vec2 = vec[3][0:-1].split(":")
        id1 = vec2[0] + ":" + vec2[1] + ":" + vec2[2]
        id2 = vec2[3] + ":" + vec2[4] + ":" + vec2[5]
        is_find = 0
        if(mm.has_key(id1)):
            is_find += 1
        if(mm.has_key(id2)):
            is_find += 1
        if(is_find == 0):
            fout.write(line.strip()+"\n")
            
        elif(is_find == 2):
            fout1.write(line.strip()+"\n")
            fout3.write(id1+"\n"+id2+"\n")
        else:
            fout3.write(id1+"\n"+id2+"\n")
            fout2.write(line.strip()+"\n")
        
    return 0

def del_sam_id_3(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            pass
        else:
            fout.write(line.strip()+"\n")
    return 0

def create_5_file(file_in,file_in1,file_out1,file_in2,file_out2,file_out3):
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    fout3 = open(file_out3,'w')
    mm = {}
    mm1 = {}
    mm2 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            fout1.write(line.strip()+"\n")
            mm1.setdefault(vec[1],"")
    for line in open(file_in2):
        vec = line.strip().split("\t")
        if(mm1.has_key(vec[0])):
            fout2.write(line.strip()+"\n")
            vec1 = vec[0].split(":")
            contig = vec1[0]
            site = vec1[1]
            if(mm2.has_key(contig)):
                mm2[contig] += site +":"
            else:
                mm2.setdefault(contig,site+":")
    for key in mm2:
        vec = mm2[key][0:-1].split(":")
        contig = key
        old_site = ""
        aList = list()
        for ii in range(0,len(vec)):
            aList.append(int(vec[ii]))
        aList.sort()
        for ii in range(0,len(aList)):
            old_site += str(aList[ii]) + ":"
        fout3.write(key+"\t"+str(len(aList)) + "\t"+ old_site+"\n")
        
        
def del_5_new(file_in,file_in1,file_out1,file_in2,file_out2,file_out3):
    fout1 = open(file_out1,'w')
    fout2 = open(file_out2,'w')
    fout3 = open(file_out3,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        if(mm.has_key(vec[0])):
            mm1.setdefault(vec[1],"")
        else:
            fout1.write(line.strip()+"\n")
            
    for line in open(file_in2):
        vec = line.strip().split("\t")
        if(mm1.has_key(vec[0])):
            vec2 = vec[3][0:-1].split(":")
            id1 = vec2[0] + ":" + vec2[1] + ":" + vec2[2]
            id2 = vec2[3] + ":" + vec2[4] + ":" + vec2[5]
            fout3.write(id1+"\n"+id2+"\n")
        else:
            fout2.write(line.strip()+"\n")
        
def link_group_kk(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    p_contig = ""
    p_site = 0
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        contig = vec[2]
        site = int(vec[3])
        id = vec[0]
        if(p_contig == ""):
            p_contig = contig
            mm.setdefault(site,id +":")
        elif(p_contig != contig):
            m_key = min(mm)
            vec1 = mm.pop(m_key).split(":")
            i_del = ""
            if(len(vec1) == 7):
                id1 = vec1[0] + ":" + vec1[1] + ":" + vec1[2]
                id2 = vec1[3] + ":" + vec1[4] + ":" + vec1[5]
                i_del = vec1[0] + ":" + vec1[1]
                fout.write(id1+"\n"+id2+"\n")
                s_out = ""
                for key in mm:
                    s_out += mm[key] + "#"
                fout1.write(contig + "\t" +i_del+"\t" + str(len(mm)) + "\t" + s_out+"\n")
            else:
                print(contig+":"+str(m_key)+"\t")
                
            s_out = ""
            for key in mm:
                s_out += mm[key] + "#"
            fout1.write(contig + "\t" +i_del+"\t" + str(len(mm)) + "\t" + s_out+"\n")
            mm.clear()
            mm.setdefault(site,id +":")
            p_contig = contig
        else:
            if(site - p_site <= 60):
                if(mm.has_key(site)):
                    mm[site] += id +":"
                else:
                    mm.setdefault(site,id +":")
            else:
                m_key = min(mm)
                vec1 = mm.pop(m_key).split(":")
                i_del = ""
                if(len(vec1) == 7):
                    id1 = vec1[0] + ":" + vec1[1] + ":" + vec1[2]
                    id2 = vec1[3] + ":" + vec1[4] + ":" + vec1[5]
                    i_del = vec1[0] + ":" + vec1[1]
                    fout.write(id1+"\n"+id2+"\n")
                    s_out = ""
                    for key in mm:
                        s_out += mm[key] + "#"
                    fout1.write(contig + "\t" +i_del+"\t" + str(len(mm)) + "\t" + s_out+"\n")
                else:
                    print(contig+":"+str(m_key)+"\t")
                mm.clear()
                mm.setdefault(site,id +":")
            
        p_site = site
        
def link_group(file_in,file_out1):
    fout1 = open(file_out1,'w')
    p_contig = ""
    p_site = 0
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        contig = vec[2]
        site = int(vec[3])
        id = vec[0]
        if(p_contig == ""):
            p_contig = contig
            mm.setdefault(site,id +":")
        elif(p_contig != contig):
            m_key = min(mm)
            vec1 = mm[m_key].split(":")
            if(len(vec1) == 7):
                s_out = ""
                for key in mm:
                    s_out += mm[key] + "#"
                fout1.write(p_contig +"\t" + str(len(mm)) + "\t" + s_out+"\n")
            else:
                print(p_contig+":"+str(m_key)+"\t")
                
            mm.clear()
            mm.setdefault(site,id +":")
            p_contig = contig
        else:
            if(site - p_site <= 60):
                if(mm.has_key(site)):
                    mm[site] += id +":"
                else:
                    mm.setdefault(site,id +":")
            else:
                m_key = min(mm)
                vec1 = mm[m_key].split(":")
                if(len(vec1) == 7):
                    s_out = ""
                    for key in mm:
                        s_out += mm[key] + "#"
                    fout1.write(contig +"\t" + str(len(mm)) + "\t" + s_out+"\n")
                else:
                    print(contig+":"+str(m_key)+"\t")
                mm.clear()
                mm.setdefault(site,id +":")
            
        p_site = site


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
        fout.write(key+"\t"+str(len(vec)) + "\t"+mm[key]+"\n")
        
def snp_kp(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0]
        cc_snp = int(vec[4])
        mm.setdefault(id,cc_snp)
        mm1.setdefault(id,vec[5])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        contig = vec[0]
        vec1 = vec[2][0:-1].split("#")
        
        mmk = {}
        s_jia = ""
        is_aa = 1
        for ii in range(0,len(vec1)):
            vec2 = vec1[ii][0:-1].split(":")
            if(len(vec2) == 6):
                id1 = vec2[0] + ":" + vec2[1] + ":" + vec2[2]
                id2 = vec2[3] + ":" + vec2[4] + ":" + vec2[5]
                site1 = int(vec2[1])
                site2 = int(vec2[4])
                if(site1 != site2):
                    is_aa = 0
                    print(id1)
                elif(vec2[0] != contig or vec2[3] != contig):
                    is_aa = 4
                else:
                    iop = contig+":"+str(site1)
                    if(mm.has_key(iop)):
                        cc = mm[iop]
                        s_jia += str(cc)+"#"
                        if(mmk.has_key(cc)):
                            mmk[cc] += str(site1) + "#"
                        else:
                            mmk.setdefault(cc,str(site1) + "#")
                    else:
                        print(iop,line)
                        is_aa = 2
            else:
                is_aa = 0
        if(is_aa == 1):
            min_key = min(mmk)
            vec_min = mmk[min_key][0:-1].split("#")
            aList = list()
            for ii in range(0,len(vec_min)):
                aList.append(vec_min[ii])
            aList.sort()
            dai_site = aList[len(aList)/2]
            cid = contig+":"+str(dai_site)
            snp = mm1[cid]
            fout.write(line.strip() + "\t"+s_jia[0:-1] + "\t"+str(min_key) + "\t"+cid+"\t"+snp+"\n")
        else:
            print(line.strip()+"\t"+str(is_aa))
            
        
def db_reads(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[5],"")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        contig = vec[0]
        if(mm.has_key(contig)):
            fout.write(line.strip())

def mdgbb(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        vec1 = vec[2][0:-1].split(":")
        mm.setdefault(">" + vec1[0],"")
        mm.setdefault(">" + vec1[1],"")
    id = ""
    for line in open(file_in1):
        line = line.strip()
        if(line[0] == ">"):
            id = line
        else:
            if(mm.has_key(id)):
                fout.write(id+"\n"+line+"\n")
    return mm
            
def tiqu_x2_det(file_in):
    for line in open(file_in):
        vec = line.strip().split("\t")
        site = vec[0].split(":")[1]
        val = os.popen("samtools view x2.sort.bam " + vec[0] + "-" + site)
        ii = 0
        for line in val.readlines():
            ii += 1
        val = os.popen("samtools view 12_29/x3.sort.bam " + vec[0] + "-" + site)
        it = 0
        for line in val.readlines():
            it += 1
        if(it != ii):
            print(line.strip())
                    
if __name__ == '__main__':
    #pp_cs('D://down1//s_k1','D://down1//s_k','D://down1//s_k2')
    #get_snp_reads('D://down1//log_aa','D://down1//log_bb')
    #reduance_snp('D://down1//tt_k','D://down1//tt_w')
    #det_tj('D://down1//4.w','D://down1//4.d')
    #create_ab('/Users/yangcheng/Documents/瓒呯畻鏁版嵁/r1','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/r1_a','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/r1_b')
    #summary_dp_new('D://down1//op1','D://down1//op2','D://down1//op3')
    #fx_tj('D://down1//i1','D://down1//i1.tj')
    #kf_33('D://down1//3_2_dp28','D://down1//3_2_dp28.kf')
    #create_mhd_data('D://down1//3_2_dp28.kf','D://down1//3_2_dp28.mhd')
    #kf_33_1_2('D://down1//3_2','D://down1//3_2.kf')
    #create_mhd_data('D://down1//3_2.kf','D://down1//3_2.mhd')
    #qj_tj('D://down1//r_site','D://down1//r_site.txt')
    #qj_tj('D://down1//1438_site','D://down1//1438_site_15k.txt',15000)
    #qj_tj_2('D://down1//1737_site','D://down1//1737_site_10k.txt',10000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_1k.txt',1000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_2k.txt',2000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_3k.txt',3000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_4k.txt',4000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_5k.txt',5000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_6k.txt',6000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_7k.txt',7000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_8k.txt',8000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_9k.txt',9000)
    #qj_tj_2('D://down1//1438_site','D://down1//1438_site_10k.txt',10000)
    #count_tj('D://down1//33.firsthand','D://down1//33.secondhand')
    #llb_123('D://down1//69_forKF.tq','D://down1//69_forKF.aa','D://down1//69_forKF.mhd')
    #llb_33_11('D://down1//69_forKF_6.63.aa','D://down1//69_forKF_6.63.kf','D://down1//69_forKF_6.63.mhd')
    #site_329('D://down1//329_site','D://down1//329_site.js')
    #test_llb(1.00)
    #test_llb(2.00)
    #test_llb(3.00)
    '''
    qj_tj_2('D://down1//83_site','D://down1//83_site_5k',5000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_10k',10000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_15k',15000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_20k',20000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_25k',25000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_30k',30000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_35k',35000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_40k',40000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_45k',45000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_50k',50000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_55k',55000)
    qj_tj_2('D://down1//83_site','D://down1//83_site_60k',60000)
    '''
    #kf_1_1_2_1('D://down1//83_forkf.tq','D://down1//83_forkf.3')
    #kf_3_1_1_1('D://down1//83_forkf.tq','D://down1//83_forkf.2')
    #site_329('D://down1//1w_site','D://down1//1w_site.js')
    #gene_3490('D://down1//3490.tiqu','D://down1//3490.gene')
    #tj_30_cx('D://down1//30_dead.txt','D://down1//16_83.fx','D://down1//83_16.re')
    #fx_3490('D://down1//83_3490.fx','D://down1//83_3490.ab')
    #zhuanhuan_3490('D://down1//83_3490.ab','D://down1//83_3490.83')
    #del_re('D://down1//ui','D://down1//ui.txt')
    '''
    ss = "368525986_6_A:"
    print(ss[0:-1].split(":"));
    a = [4,5,2,3,1]
    
    print("this->",np.var(a))
    '''
    #reads_mean('D://down1//sample','D://down1//sample.txt')
    #tiqu_ref_50('D://down1//2','D://down1//2.txt')
    #gene_type('D://down1//w1','D://down1//w2.txt')
    #pute_dd('D://down1//i_id.sort','D://down1//i_id.dd')
    #zebre_in_gene('D://down1//bmy_all_site.txt','D://down1//z11.gene','D://down1//kk.txt')
    #tj_site('D://down1//t_i','D://dowwn1//match_sort.tj')
    #aa('D://down1//kk')
    #gff_mrna_static('D://down1//z11_small.gff','D://down1//z11_mrna.tj','D://down1//id_1.js')
    #kk('D://down1//z11_mrna.tj','D://down1//kk.tj')
    #ddd1('D://down1//c.ty','D://down1//c_log','D://down1//cc.txt')
    #cut_100('D://down1//w1','D://down1//w1.tj')
    #cat_fa('D://down1//tt','D://down1//tt.tj','D://down1//tt1.tj')
    #is_onesite_new('/Users/yangcheng/Documents/瓒呯畻鏁版嵁/t2','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/t1','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/t3')
    #find_snp('/Users/yangcheng/Documents/瓒呯畻鏁版嵁/t1','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/t2','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/t3')
    #get_num('/Users/yangcheng/Documents/瓒呯畻鏁版嵁/181.tiqu','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/181.tj')
    #print_awk('/Users/yangcheng/Documents/瓒呯畻鏁版嵁/k1','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/k2')
    #chai_3_181('D://down1//361.fa','D://down1//361_3.fa')
    #snp_aa('D://down1//k1','D://down1//k1.fa')
    #is_one_site('/Users/yangcheng/Documents/瓒呯畻鏁版嵁/kk','/Users/yangcheng/Documents/瓒呯畻鏁版嵁/kk1')
    #match_2('D://down1//t1.18','D://down1//t2.18','D://down1//t3.18','D://down1//t4.18')
    #find_x5_x10('D://down1//tt','D://down1//tt1','D://down1//tt2')
    #x15_snp('D://down1//kt','D://down1//wt','D://down1//it')
    #tj_pp('E://super_down//k1','E://super_down//k2','E://super_down//it')
    #pj_101('D://down1//lo.pv','D://down1//lu','D://down1//it','D://down1//it1','D://down1//ki.w')
    #crea_fff('D://down1//k.101','D://down1//wt.101')
    #yz_ff('D://down1//op','D://down1//ko','D://down1//it')
    #fx52_1('D://down1//ki','D://down1//ki.txt')
    #kf_52('E://super_down//iu','E://super_down//kk.pj')
    #kf_tj('E://super_down//kk.pj','5.99')
    #sex_1('E://super_down//gy_1','E://super_down//gy_1.kf')
    #sex_4('E://super_down//kk.o','E://super_down//kk_d.male')
    #pj_md('E://super_down//contig.length','E://super_down//kt','E://super_down//kt.kk')
    #get_00('E://super_down//ref.fa','E://super_down//site_sort.id','E://super_down//kk')
    #get_11('E://super_down//kk','E://super_down//gy.sam','E://super_down//ww')
    #k1('D://down1//file_male','D://down1//male.snp')
    #k2('D://down1//male.snp','D://down1//female.snp','D://down1//kk')
    #kf_tt('D://down1//52_value.kf','D://down1//52_value.kf.tj')
    #cx_pd('D://down1//kf','D://down1//km','D://down1//it')
    #jz1('E://super_down//k1','E://super_down//k1.jz')
    #shai_x10_2type('D://down1//ko','D://down1//ko.snp','D://down1//ko.tt')
    ss = "1234"
    #print(ss[1:])
    #tiqu_aa('D://down1//find_r1.fq','D://down1//find_r2.fq','D://down1//kk.c','D://down1//kk.w')
    #in_gene('D://down1//k.id','D://down1//gene_id','D://down1//kk.w')
    #tiqu_gat('E://super_down//bl//2.fa','E://super_down//bl//223_message.txt','E://super_down//bl//80.mess')
    #tiqu_site_kk('D://down1//male.reads','D://down1//ref.fa','D://down1//ref.site')
    #pd_link_snp('E://super_down//12.ac2','E://super_down//12.ui')
    #sex_mark_pair_1('E://super_down//he.tj','E://super_down//he.popo')
    #sex_mark_pair_2('E://super_down//he.popo','E://super_down//731.reads','E://super_down//731.tj')
    #pd_link_snp('D://down1//1.ac1','D://down1//1.su')
    #find_1178_bl_gap('E://super_down//kk','E://super_down//kk.popo')
    #kk_pp('D://down1//site_sort.cc','D://down1//site_sort_f.cc')
    #kk_pp('D://down1//shai_sort.cc','D://down1//12_27.cc')
    #kk_pp('D://down1//shai_sort_old.cc','D://down1//12_27_old.cc')
    #qj_iii('G://yc//jhs//100_2.gap','D://down1//ref_ch.id')
    #hubu('G://yc//jhs//100_2.gap','D://down1//ref_ch.id')
    #kk_pp_ll('D://down1//test','D://down1//test.fa')
    #mess2('E://super_down//y.id','E://super_down//he.zl','E://super_down//y.site')
    #mess2('D://down1//x.id','D://down1//he.zl','D://down1//x.site')
    #snp_2('E://super_down//test.t1','E://super_down//test.snp','E://super_down//snp.site','E://super_down//test.err')
    #snp_reads_pv('D://down1//kt','D://down1//kt.cc')
    #snp_del_cl('D://down1//pv1','D://down1//pv1.cc')
    #is_only('E://super_down//test.sam','E://super_down//test.cc','E://super_down//test.err')
    #is_pairaa('E://super_down//test.cc','E://super_down//test.kk')
    ss = ">Lachesis_group18__9_contigs__length_31774926:5957029:2_1"
    #print(ss[1:])
    #create_fa_kk('E://super_down//tt','E://super_down//tt.kk')
    #del_one('D://down1//kk','D://down1//kk.cc','D://down1//kk.err')
    
    snp_kp('D://down1//g0','D://down1//g1','D://down1//g2')
    