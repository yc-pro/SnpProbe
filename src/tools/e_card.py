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

if __name__ == "__main__":
    #test()
    #unzip('E://super_down//461_x5.sort.bam')
    find_x5('E://super_down//k1','E://super_down//k2')
    #pysam_tt()