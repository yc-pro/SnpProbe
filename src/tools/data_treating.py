#coding:gbk
''' 
Created on 2020年1月14日
苹果数据测试
@author: YC
'''

from collections import Counter
import operator
#import scipy.stats as st
import numpy as np
from _socket import AI_ADDRCONFIG
import math
import os
import string


def x5_reads(file_in,file_out):
    fout = open(file_out,'w')
    id = ""
    for line in open(file_in):
        line = line.strip()
        if(line[0] == ">"):
            id = line
        else:
            fout.write(id+"_a\n"+line+"\n")
            fout.write(id+"_b\n"+line+"\n")
            fout.write(id+"_c\n"+line+"\n")
            fout.write(id+"_d\n"+line+"\n")
            fout.write(id+"_e\n"+line+"\n")
            
def x5_id(file_in,file_in1,file_out):
    fout = open(file_out,'w')
    mm = {}
    for line in open(file_in):
        vec = line.strip().split(":")
        mm.setdefault(vec[2],"")
    id = ""
    for line in open(file_in1):
        line = line.strip()
        if(line[0] == ">"):
            id = line
        else:
            if(mm.has_key(line)):
                mm[line] = id
    for key in mm:
        fout.write(mm[key] + "\n" + key+"\n")
        
                    
if __name__ == '__main__':
    print("Hello word!!!git")