#coding:gbk
'''
Created on 2019-1-29

@author: Genome
'''
import Bio
#import pysam_windows as pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import bgzf

def test():
    #create a sequence object
    my_seq = Seq('CATGTAGACTAG')
     
    #print out some details about it
    print ('seq %s is %i bases long' % (my_seq, len(my_seq)))
    print ('reverse complement is %s' % my_seq.reverse_complement())
    print ('protein translation is %s' % my_seq.translate())

    
def unzip(file_in):
    #orchid_dict = SeqIO.index(file_in, "genbank")
    #print(len(orchid_dict))
    handle = bgzf.BgzfReader(file_in, "r")
    aa = handle.read(100000)
    print(aa)
    #print(handle.readline())
    
def pysam_tt():
    #samfile = pysam.csamtools()
    #print(samfile)
    pass
    
def find_a_ab(file_in):
    for line in open(file_in):
        vec1 = line.strip().split(" ")
        ii = 0
        for iw in range(2,len(vec1)):
            ik = vec1[iw]
            if(vec1[0] == "*tig00000005_pilon_pilon_6481"):
                print(ik)
            if(ik != "ab"  and ik != "a" and ik != "-"):
                ii += 1
                break
        if(ii == 0):
            print(line.strip());
    return 0

def cl_smj_1(file_in):
    for line in open(file_in):
        vec1 = line.strip().split(" ")
        ii = 0
        for iw in range(2,len(vec1)):
            ik = vec1[iw]
            if(vec1[0] == "*tig00000005_pilon_pilon_6481"):
                print(ik)
            if(ik != "ab"  and ik != "a" and ik != "-"):
                ii += 1
                break
        if(ii == 0):
            print(line.strip());
    return 0

def cl_smj_2(file_in,file_out):
    fout = open(file_out,'w')
    mm1 = {}
    mm2 = {}
    icc = 1
    for line in open(file_in):
        vec1 = line.strip().split(" ")
        sk = ""
        if(icc >= 4):
            for ii in range(2,len(vec1)):
                sk += vec1[ii] + ":"
            mm1.setdefault(vec1[0],sk)
            mm2.setdefault(vec1[0],sk)
        icc += 1
    print(len(mm1))
    for key in mm1:
        vec = mm1[key].split(":")
        mm2.pop(key)
        str_same = ""
        ii_same = 0
        print(key)
        for key222 in mm2:
            is_same = 1
            vec1 = mm2[key222].split(":")
            for ii in range(0,len(vec1)):
                i1 = vec[ii]
                i2 = vec1[ii]
                if(i1 == i2):
                    is_same = 0
            if(is_same == 1):
                str_same += key222 + ":"
                ii_same += 1
        print(key + "\t" + str(ii_same) + "\t" + str_same)
    return 0

if __name__ == "__main__":
    #test()
    #unzip('E://super_down//461_x5.sort.bam')
    #find_a_ab('D://down1//kk')
    cl_smj_2('D://down1//kk','D://down1//kk.aa')
    #pysam_tt()