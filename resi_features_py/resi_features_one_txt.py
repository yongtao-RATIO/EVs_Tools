'''
    2022.08.18 19:20
    我们对受体结合位置处的 氨基酸残基的理化特征 统计
    这里主要总结了氨基酸的 分子量、质谱分子量、体积、范德华体积、等电点、极性、疏水性、酸碱性、芳香族、亲水性评估值、可及性、可塑性和抗原概率
    得到多个TXT文件，每个TXT对应一种特征

'''


#1氨基酸的分子量
def mole_weight(resi):

    dict_weight = {
        "A\n":"89.09", "C\n":"121.16", "D\n":"133.1", "E\n":"147.13", "F\n":"165.19",
        "G\n":"75.07", "H\n":"155.16", "I\n":"131.17", "K\n":"146.19", "L\n":"131.17",
        "M\n":"149.21", "N\n":"132.12", "P\n":"115.13", "Q\n":"146.15", "R\n":"174.2",
        "S\n":"105.09", "T\n":"119.12", "V\n":"117.15", "W\n":"204.23", "Y\n":"181.19",
        "\n":""}

    return dict_weight[resi] + "\n"

#2氨基酸的质谱分子量
def msm_weight(resi):

    dict_msmweight = {
        "A\n":"71.07", "C\n":"103.1", "D\n":"115.08", "E\n":"129.11", "F\n":"147.17",
        "G\n":"57.05", "H\n":"137.14", "I\n":"113.15", "K\n":"128.17", "L\n":"113.15",
        "M\n":"131.19", "N\n":"114.1", "P\n":"97.11", "Q\n":"128.13", "R\n":"156.18",
        "S\n":"87.07", "T\n":"101.14", "V\n":"99.13", "W\n":"186.2", "Y\n":"163.17",
        "\n":""}

    return dict_msmweight[resi] + "\n"

#3氨基酸的体积
def volume(resi):

    dict_volume = {
        "A\n":"88.6", "C\n":"108.5", "D\n":"111.1", "E\n":"138.4", "F\n":"189.9",
        "G\n":"60.1", "H\n":"153.2", "I\n":"166.7", "K\n":"168.6", "L\n":"166.7",
        "M\n":"162.9", "N\n":"114.1", "P\n":"112.7", "Q\n":"143.8", "R\n":"173.4",
        "S\n":"89", "T\n":"116.1", "V\n":"140", "W\n":"227.8", "Y\n":"193.6",
        "\n":""}

    return dict_volume[resi] + "\n"

#4氨基酸的范德华体积
def vdw_volume(resi):

    dict_vdwvolume = {
        "A\n":"67", "C\n":"86", "D\n":"91", "E\n":"109", "F\n":"135",
        "G\n":"48", "H\n":"118", "I\n":"124", "K\n":"135", "L\n":"124",
        "M\n":"124", "N\n":"96", "P\n":"90", "Q\n":"114", "R\n":"148",
        "S\n":"73", "T\n":"93", "V\n":"105", "W\n":"163", "Y\n":"141",
        "\n":""}

    return dict_vdwvolume[resi] + "\n"

#5氨基酸的等电点
def PI(resi):

    dict_PI = {
        "A\n":"6", "C\n":"5.02", "D\n":"2.77", "E\n":"3.22", "F\n":"5.48",
        "G\n":"5.97", "H\n":"7.47", "I\n":"5.94", "K\n":"9.59", "L\n":"5.98",
        "M\n":"5.74", "N\n":"5.41", "P\n":"6.3", "Q\n":"5.65", "R\n":"11.15",
        "S\n":"5.68", "T\n":"5.64", "V\n":"5.96", "W\n":"5.89", "Y\n":"5.66",
        "\n":""}

    return dict_PI[resi] + "\n"

#6氨基酸的极性
def polar(resi):

    '''
    N:非极性
    P:极性不带电
    +:极性带正电
    -:极性带负电

    N   AFILMPVW    8
    P   CGNQSTY     7
    +   HKR         3
    -   DE          2
    '''

    dict_polar = {
        "A\n":"N", "C\n":"P", "D\n":"-", "E\n":"-", "F\n":"N",
        "G\n":"P", "H\n":"+", "I\n":"N", "K\n":"+", "L\n":"N",
        "M\n":"N", "N\n":"P", "P\n":"N", "Q\n":"P", "R\n":"+",
        "S\n":"P", "T\n":"P", "V\n":"N", "W\n":"N", "Y\n":"P",
        "\n":""}

    return dict_polar[resi] + "\n"

#7氨基酸的疏水性
def phobic(resi):

    '''
    B:疏水
    L:亲水

    B   AFILMPVW        8
    L   CDEGHKNQRSTY    12
    '''

    dict_phobic = {
        "A\n":"B", "C\n":"L", "D\n":"L", "E\n":"L", "F\n":"B",
        "G\n":"L", "H\n":"L", "I\n":"B", "K\n":"L", "L\n":"B",
        "M\n":"B", "N\n":"L", "P\n":"B", "Q\n":"L", "R\n":"L",
        "S\n":"L", "T\n":"L", "V\n":"B", "W\n":"B", "Y\n":"L",
        "\n":""}

    return dict_phobic[resi] + "\n"

#8氨基酸的酸碱性
def acid_base(resi):

    '''
    N:中性
    A:酸性
    B:碱性

    N   ACFGILMNPQSTVWY     15
    A   DE                  2
    B   HKR                 3
    '''

    dict_acid_base = {
        "A\n":"N", "C\n":"N", "D\n":"A", "E\n":"A", "F\n":"N",
        "G\n":"N", "H\n":"B", "I\n":"N", "K\n":"B", "L\n":"N",
        "M\n":"N", "N\n":"N", "P\n":"N", "Q\n":"N", "R\n":"B",
        "S\n":"N", "T\n":"N", "V\n":"N", "W\n":"N", "Y\n":"N",
        "\n":""}

    return dict_acid_base[resi] + "\n"

#9芳香族氨基酸
def aromatic(resi):

    '''
    A:是
    N:否

    A   FHWY    4
    N   ...     16
    '''

    dict_aromatic = {
        "A\n":"N", "C\n":"N", "D\n":"N", "E\n":"N", "F\n":"A",
        "G\n":"N", "H\n":"A", "I\n":"N", "K\n":"N", "L\n":"N",
        "M\n":"N", "N\n":"N", "P\n":"N", "Q\n":"N", "R\n":"N",
        "S\n":"N", "T\n":"N", "V\n":"N", "W\n":"A", "Y\n":"A",
        "\n":""}

    return dict_aromatic[resi] + "\n"

#10亲水性
def hydrovalue(resi):

    '''
    氨基酸指数表（AAindex）: https://www.genome.jp/aaindex/
    查询: HOPT810101
    亲水性定量
    '''

    dict_hydrovalue = {
        "A\n":"-0.5", "C\n":"-1", "D\n":"3", "E\n":"3", "F\n":"-2.5",
        "G\n":"0", "H\n":"-0.5", "I\n":"-1.8", "K\n":"3", "L\n":"-1.8",
        "M\n":"-1.3", "N\n":"0.2", "P\n":"0", "Q\n":"0.2", "R\n":"3",
        "S\n":"0.3", "T\n":"-0.4", "V\n":"-1.5", "W\n":"-3.4", "Y\n":"-2.3",
        "\n":""}

    return dict_hydrovalue[resi]

#11表面可及性
def accessuface(resi):

    '''
    氨基酸指数表（AAindex）: https://www.genome.jp/aaindex/
    查询: JANJ780101
    平均可接触表面面积
    '''

    dict_accessuface = {
        "A\n": "27.8", "C\n": "15.5", "D\n": "60.6", "E\n": "68.2", "F\n": "25.5",
        "G\n": "24.5", "H\n": "50.7", "I\n": "22.8", "K\n": "103", "L\n": "27.6",
        "M\n": "33.5", "N\n": "60.1", "P\n": "51.5", "Q\n": "68.7", "R\n": "94.7",
        "S\n": "42", "T\n": "45", "V\n": "23.7", "W\n": "34.7", "Y\n": "55.2",
        "\n": ""}

    return dict_accessuface[resi]

#12可塑性
def flexibility(resi):

    '''
    氨基酸指数表（AAindex）: https://www.genome.jp/aaindex/
    查询: KARP850103
    可塑性
    '''

    dict_flexibility = {
        "A\n": "0.892", "C\n": "0.925", "D\n": "0.932", "E\n": "0.933", "F\n": "0.914",
        "G\n": "0.923", "H\n": "0.894", "I\n": "0.827", "K\n": "1.057", "L\n": "0.921",
        "M\n": "0.804", "N\n": "0.93", "P\n": "0.932", "Q\n": "0.885", "R\n": "0.901",
        "S\n": "0.923", "T\n": "0.934", "V\n": "0.913", "W\n": "0.803", "Y\n": "0.837",
        "\n": ""}

    return dict_flexibility[resi]

#13抗原概率
def proantigen(resi):

    '''
    抗原概率
    '''

    dict_proantigen = {
        "A\n": "0.115", "C\n": "-0.12", "D\n": "0.065", "E\n": "-0.071", "F\n": "-0.141",
        "G\n": "-0.184", "H\n": "0.312", "I\n": "-0.292", "K\n": "0.206", "L\n": "0.075",
        "M\n": "-0.385", "N\n": "-0.077", "P\n": "-0.053", "Q\n": "-0.011", "R\n": "0.058",
        "S\n": "-0.026", "T\n": "-0.045", "V\n": "-0.013", "W\n": "-0.114", "Y\n": "0.013",
        "\n": ""}

    return dict_proantigen[resi]


resitxt = open('residues.txt','r')

moletxt = open('mole_weight.txt','w+',encoding='utf-8')
msmtxt = open('msm_weight.txt','w+',encoding='utf-8')
voltxt = open('volume.txt','w+',encoding='utf-8')
vdwtxt = open('vdw_volume.txt','w+',encoding='utf-8')
PItxt = open('PI.txt','w+',encoding='utf-8')
pltxt = open('polar.txt','w+',encoding='utf-8')
phtxt = open('phobic.txt','w+',encoding='utf-8')
abtxt = open('acid_base.txt','w+',encoding='utf-8')
arotxt = open('aromatic.txt','w+',encoding='utf-8')
hydtxt = open('hydrovalue.txt','w+',encoding='utf-8')
acctxt = open('accessuface.txt','w+',encoding='utf-8')
fletxt = open('flexibility.txt','w+',encoding='utf-8')
epitxt = open('proantigen.txt','w+',encoding='utf-8')

resilist = resitxt.readlines()


a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13 = '','','','','','','','','','','','',''

for resi in resilist:
    a1 += mole_weight(resi)
    a2 += msm_weight(resi)
    a3 += volume(resi)
    a4 += vdw_volume(resi)
    a5 += PI(resi)
    a6 += polar(resi)
    a7 += phobic(resi)
    a8 += acid_base(resi)
    a9 += aromatic(resi)
    a10 += hydrovalue(resi)
    a11 += accessuface(resi)
    a12 += flexibility(resi)
    a13 += proantigen(resi)

moletxt.write(a1)
msmtxt.write(a2)
voltxt.write(a3)
vdwtxt.write(a4)
PItxt.write(a5)
pltxt.write(a6)
phtxt.write(a7)
abtxt.write(a8)
arotxt.write(a9)
hydtxt.write(a10)
acctxt.write(a11)
fletxt.write(a12)
epitxt.write(a13)

#关闭文件
resitxt.close()
moletxt.close()
msmtxt.close()
voltxt.close()
vdwtxt.close()
PItxt.close()
pltxt.close()
phtxt.close()
abtxt.close()
arotxt.close()
hydtxt.close()
acctxt.close()
fletxt.close()
epitxt.close()