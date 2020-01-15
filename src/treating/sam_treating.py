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

if __name__ == '__main__':
    print("hello word")