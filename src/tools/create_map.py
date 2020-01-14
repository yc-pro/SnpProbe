'''
Created on 2019Äê1ÔÂ29ÈÕ

@author: Genome
'''
import pandas as pd
import matplotlib.pyplot as plt

def load_data(file_in):
    mm = {}
    ii = 1
    for line in open(file_in):
        vec = line.strip().split("\t")
def map_aa(unrate):
    first_twele = unrate[0:12]
    plt.plot(first_twele['DATE'],first_twele['VALUE'])
    plt.xticks(rotation = 45)
    plt.xlabel('Month')
    plt.ylabel('Unemployment')
    plt.title('ssss')
    plt.show()
    
if __name__ == "__main__":