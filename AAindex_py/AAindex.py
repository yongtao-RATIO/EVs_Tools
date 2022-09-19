'''
    2022.09.17 20:01
    此代码旨在用爬虫的手段获得 AAindex数据库 的信息
    用法：提供存有 index 的 <indexlist.txt>，运行获得集合的氨基酸特征CSV文件

'''

import requests
import re
import pandas as pd
import json

#   <1.得到 AAindex索引目录>
url = 'https://www.genome.jp/aaindex/AAindex/list_of_indices'
response = requests.get(url)
#print(response)
html_txt = response.content.decode('utf-8')
#print(html_txt)

list = html_txt.split(sep='\n')
print(list[0])                                                                  ### 一共有566个氨基酸特征索引 List of 566 Amino Acid Indices in AAindex ver.9.2
list_index = list[5:571]                                                        ### 从第 6 行 到 第 571 行，是index信息，而索引值分别对应 5 和 570，取切片
#print(list_index)


#   <2.根据目录 创建 {'索引'：'解释'} 的字典>
index_dict = { i.split(' ',1)[0]:i.split(' ',1)[1] for i in list_index }
#print(index_dict)
# for l in index_dict.keys():
#     print(l)



#   <3.爬取氨基酸特征>


def pdindex(interpret,index):
    url_index = 'https://www.genome.jp/entry/aaindex:' + index
    tree_txt = requests.get(url_index).content.decode('utf-8')
    #print(tree_txt)
    a = re.findall(r'\sA/L *R/K *N/M *D/F *C/P *Q/S *E/T *G/W *H/Y *I/V\n'
                   r'\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\n'
                   r'\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\s.*?\n',tree_txt)
    b = a[0].split()
    #print(b)

    df = pd.DataFrame({ '%s'%index : ['%s'%interpret,'','','','','','','','','','','','','','','','','','',''],
                        'index':["A", "C", "D", "E", "F",
                                 "G", "H", "I", "K", "L",
                                 "M", "N", "P", "Q", "R",
                                 "S", "T", "V", "W", "Y"],
                        'data':[b[10], b[14], b[13], b[16], b[23],
                                b[17], b[18], b[19], b[21], b[20],
                                b[22], b[12], b[24], b[15], b[11],
                                b[25], b[26], b[29], b[27], b[28]]
                       })

    dict = {"A":b[10], "C":b[14], "D":b[13], "E":b[16], "F":b[23],
            "G":b[17], "H":b[18], "I":b[19], "K":b[21], "L":b[20],
            "M":b[22], "N":b[12], "P":b[24], "Q":b[15], "R":b[11],
            "S":b[25], "T":b[26], "V":b[29], "W":b[27], "Y":b[28]
            }
    with open('datadict.txt','a') as h:
        h.write(str(dict) + '\n')

    return df

#   <4.合并为一个 df 输出为 CSV>
def csvindex(indexlist):
    frames = []
    for j in indexlist:                                                         ### 例子：indexlist = ['BIGC670101','BULH740101','CHOC760102','CHOP780209','GEIM800106']
        index = j
        interpret = index_dict[j]
        frames.append(pdindex(interpret, index))

    df_All = pd.concat(frames,axis=1)
    print(df_All.head())
    df_All.to_csv('AAindex.csv',sep=',',index=True,header=True)


with open('indexlist.txt','r') as f:
    indexlist = [k.split()[0] for k in f.readlines()]
    csvindex(indexlist)