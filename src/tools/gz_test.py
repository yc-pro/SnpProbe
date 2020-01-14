#coding:gbk
'''
Created on 2019-1-29

@author: Genome
'''
import os
import gzip
import tarfile

#库包解压缩与压缩
def untar(fname, dirs):
    t = tarfile.open(fname)
    t.extractall(path = dirs)
     
def tar(fname):
    t = tarfile.open(fname + ".tar.gz", "w:gz")
    for root, dir, files in os.walk(fname):
        print(root, dir, files)
        for file in files:
            fullpath = os.path.join(root, file)
            t.add(fullpath)
    t.close()
    
def read_gz_file(file_in,file_out):
    f_out = gzip.open(file_out, "wb")
    if os.path.exists(file_in):
        for line in gzip.open(file_in, 'rb'):
            f_out.write(line)
    else:
        print('the path [{}] is not exist!'.format(file_in))

if __name__ == "__main__":
    read_gz_file('E://super_down//GCF_000005845.2_ASM584v2_genomic.fna.gz','E://super_down//gy.gz')
    