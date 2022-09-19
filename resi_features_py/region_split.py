'''
    2022.09.15 21:39
    根据 CVB1、CVB3、E6、E30，我们获得 EVBs受体结合位置 在 VP1-BC, VP1-CD, VP1-EF, VP1-GH, VP1-C-, VP2-EF, VP3-CD, VP3-GH, VP3-C- 上，截取这些片层的残基信息，
    进行分列操作，便于之后开展 replace 特征分析

'''


import pandas as pd
import numpy as np

df = pd.read_csv('B_loop_summary_demo1.csv')
df.set_index('serotype', inplace=True, drop=False)
# print(df)
# print(df.columns)


frames = []
for i in df.columns[2:]:
    N = df[i].str.split('',expand=True)
    frames.append(N)

df_all = pd.concat(frames,axis=1)
print(df_all)

df_all.to_csv('B_loop_summary_demo2.csv',sep=',',index=True,header=True)