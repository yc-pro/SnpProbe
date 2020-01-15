'''
Created on 2020��1��15��

@author: Genome
'''
def is_only(file_in,file_out,file_out1):
    fout = open(file_out,'w')
    fout1 = open(file_out1,'w')
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
                        fout.write(key+"\t"+ vec1[0] + "\t" + vec1[1] + "\t"  +vec1[2]+"\n")
                    else:
                        fout1.write(key+"\t"+mm[key]+"\tsmall_indel\n") #����indel
                
        else:
            #���������С�÷�
            fout1.write(key+"\t"+mm[key]+"\tmany_small_value\n")
    return 0

if __name__ == '__main__':
    pass