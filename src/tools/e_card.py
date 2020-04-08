#coding:gbk
''' 
Created on 2019-1-29

@author: Genome
'''
from collections import Counter
import operator
import numpy as np
from _socket import AI_ADDRCONFIG
import math
import os
import string

def find_x5(file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    mm.setdefault("AC",1)
    mm.setdefault("AG",1)
    mm.setdefault("AT",1)
    mm.setdefault("GA",1)
    mm.setdefault("GC",1)
    mm.setdefault("GT",1)
    mmk = {}
    for key in mm:
        file_in = key + ".sort"
        #file_in = "E://super_down//kk"
        for line in open(file_in):
            vec1 = line.strip().split(" ")
            ic = vec1[0]
            vec = vec1[1].split(":")
            id = vec[0] + ":" + vec[1]
            read = vec[2]
            if(ic == "5"):
                if(mmk.has_key(id)):
                    mmk[id] += read + ":"
                else:
                    mmk.setdefault(id,read + ":")
            else:
                fout1.write(id+"\n")
    for key in mmk:
        vv = mmk[key][0:-1].split(":")
        fout.write(key + "\t" + str(len(vv)) + "\t" + mmk[key] + "\n")
    return 0

def cut_kk(file_in):
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("/")
        vec1 = vec[1].split("-")
        vec2 = vec[len(vec)-1].split(".")
        #print("mv " + line.strip() + " gz/"+vec1[1]+"_"+vec2[1]+".fastq.gz")
        ss = "sh bi.js "+vec1[1]
        mm.setdefault(ss,1)
    for key in mm:
        print(key)

def cut_kk(file_in):
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("/")
        vec1 = vec[1].split("-")
        vec2 = vec[len(vec)-1].split(".")
        #print("mv " + line.strip() + " gz/"+vec1[1]+"_"+vec2[1]+".fastq.gz")
        ss = "sh bi.js "+vec1[1]
        mm.setdefault(ss,1)
    for key in mm:
        print(key)
   
def ya_1(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        reads = vec[9]
        if(mm.has_key(id)):
            mm[id] += reads + ":"
        else:
            mm.setdefault(id,reads+":")
    for key in mm:
        vec = mm[key][0:-1].split(":")
        fout.write(key+"\t"+str(len(vec))+"\t"+mm[key][0:-1]+"\n")
             

def ya_1(file_in,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        if(line[0] == "@"):
            continue
        vec = line.strip().split("\t")
        id = vec[0]
        reads = vec[9]
        if(mm.has_key(id)):
            mm[id] += reads + ":"
        else:
            mm.setdefault(id,reads+":")
    for key in mm:
        vec = mm[key][0:-1].split(":")
        fout.write(key+"\t"+str(len(vec))+"\t"+mm[key][0:-1]+"\n")
        
def ya_2(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mmk = {}
    id = ""
    for line in open(file_in1):
        line = line.strip()
        if(line[0] == ">"):
            id = line[1:]
        else:
            mmk.setdefault(line,id)
        
    for line in open(file_in):
        vec = line.strip().split("\t")
        id = vec[0][0:-2]
        reads = vec[2]
        if(mm.has_key(id)):
            mm[id] += reads+":"
        else:
            mm.setdefault(id,reads+":")
    for key in mm:
        vec = mm[key][0:-1].split(":")
        is_find = 0
        if(len(mmk) == 2):
            #print(vec)
            if(mmk.has_key(vec[0]) and mmk.has_key(vec[1])):
               fout.write(key+"\t"+str(len(vec))+"\t"+mm[key][0:-1]+"\n")
               
def ya_3(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    mmk = {}
    for line in open(file_in):
        line = line.strip()
        mmk.setdefault(line,1)
        
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[1]
        if(mmk.has_key(id)):
            i3 = float(vec[2])
            op = math.modf(i3)
            ii = op[1]
            if(op[0] > 0.500000000):
                ii += 0.5
            if(mm.has_key(ii)):
                mm[ii] += 1
            else:
                mm.setdefault(ii,1)
    for key in mm:
        fout.write(str(key)+"\t"+str(mm[key])+"\n")
        
def ya_4(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        line = line.strip()
        mm.setdefault(line,1)
    mmk = {}
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        if(mm.has_key(id)):
            i3 = float(vec[1])/60
            op = round(float(i3),1)
            if(mmk.has_key(op)):
                mmk[op] += 1
            else:
                mmk.setdefault(op,1)
            fout.write(line.strip()+"\t"+str(op)+"\n")
    for key in mmk:
        print(str(key)+"\t"+str(mmk[key])+"\n")
        
def ya_5(file_in,file_in1,file_in2,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[0],vec[1])
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[1]
        if(mm.has_key(id)):
            fout.write(id+"\t"+mm[id]+"\t"+vec[2]+"\t3\n")
            mm.pop(id)
    for line in open(file_in2):
        vec = line.strip().split("\t")
        id = vec[1]
        if(mm.has_key(id)):
            fout.write(id+"\t"+mm[id]+"\t"+vec[2]+"\t3\n")
            mm.pop(id)
    print(len(mm))
        
def ya_6(file_in,file_in1,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
    mm = {}
    mm1 = {}
    for line in open(file_in):
        vec = line.strip().split("\t")
        mm.setdefault(vec[3],"")
        mm.setdefault(vec[4],"")
        mm1.setdefault(vec[0],vec[3]+":"+vec[4])
        fout1.write(">"+vec[0]+"_a\n"+vec[3]+"\n>"+vec[0]+"_b\n"+vec[4]+"\n")
    for line in open(file_in1):
        vec = line.strip().split("\t")
        id = vec[0]
        if(mm.has_key(id)):
            mm[id] += vec[2]+":"
    i_total = 0
    it_p = 0
    for key in mm1:
        vec = mm1[key].split(":")
        re1 = mm[vec[0]][0:-1].split(":")
        re2 = mm[vec[1]][0:-1].split(":")
        mmk = {}
        for kk in range(0,len(re1)):
            mmk.setdefault(re1[kk],1)
        mmw = {}
        for kk in range(0,len(re2)):
            mmw.setdefault(re2[kk],1)
        #for kk in range(0,5):
        it = 1
        for kk in mmk:
            if(it_p < 5):
                fout.write(">"+key+"_a_"+str(it)+"\n"+kk+"\n")
            it_p += 1
            it += 1
        it = 1
        it_p = 0
        for kk in mmw:
            if(it_p < 5):
                fout.write(">"+key+"_b_"+str(it)+"\n"+kk+"\n")
            it_p += 1
            it += 1
        i_total += len(mmk)+len(mmw)
        print(key+"\t"+str(len(mmk))+"\t"+str(len(mmw))+"\t"+str(len(mmk)+len(mmw)))
    print(i_total)
if __name__ == "__main__":
    #test()
    #unzip('E://super_down//461_x5.sort.bam')
    #cut_kk('/Users/yangcheng/Documents/supper_down/ff')
    #pysam_tt()
    ya_6('/Users/yangcheng/Documents/supper_down/site_15','/Users/yangcheng/Documents/supper_down/ku_noN.cut','/Users/yangcheng/Documents/supper_down/ku.fa','/Users/yangcheng/Documents/supper_down/ku_ref.fa')