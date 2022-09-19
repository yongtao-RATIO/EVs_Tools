'''
    2022.07.29 22:03
    读取字符串中特定字符的位置，统计出现次数
    我们写出来用于统计各环区预测表位位置(统一标准后)和数量
'''


def index_get(spattern,sstring):                               #读取位置
    import re

    index_list = ''

    for i in re.finditer(spattern,sstring):                    #与findall不同，finditer生成一个iterator（迭代器），需通过迭代展示
        index_list+='%d'%(i.end()) + ','                       #i.end():返回结束的位置
    print(index_list)
    return index_list + '\n'

def pattern_count(spattern,sstring):

    scount = '%d'%(sstring.count(spattern))
    return scount + '\n'


spattern1 = '1'
txt0 = open('chushi01.txt','r')
txt1 = open('index_site.txt','a+',encoding='utf-8')
txt2 = open('count_site.txt','a+',encoding='utf-8')

for sstring1 in txt0.readlines():
    txt1.write(index_get(spattern1,sstring1))
    txt2.write(pattern_count(spattern1,sstring1))

txt0.close()
txt1.close()
txt2.close()