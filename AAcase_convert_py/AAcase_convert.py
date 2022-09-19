# 氨基酸英文全称转缩写

aafulname=open('aafulname.txt','r')
abbrname=open('abbreviate.txt','a+',encoding='utf-8')
str1=aafulname.readlines()
print(str1)

def aa_abbreviate(aa0):
    if aa0=='GLY':
        aa1='G'
    elif aa0=='ALA':
        aa1='A'
    elif aa0=='VAL':
        aa1='V'
    elif aa0=='LEU':
        aa1='L'
    elif aa0=='ILE':
        aa1='I'
    elif aa0=='PHE':
        aa1='F'
    elif aa0=='PRO':
        aa1='P'
    elif aa0=='TRP':
        aa1='W'
    elif aa0=='SER':
        aa1='S'
    elif aa0=='TYR':
        aa1='Y'
    elif aa0=='CYS':
        aa1='C'
    elif aa0=='MET':
        aa1='M'
    elif aa0=='ASP':
        aa1='D'
    elif aa0=='ASN':
        aa1='N'
    elif aa0=='GLN':
        aa1='Q'
    elif aa0=='GLU':
        aa1='E'
    elif aa0=='THR':
        aa1='T'
    elif aa0=='LYS':
        aa1='K'
    elif aa0=='ARG':
        aa1='R'
    elif aa0=='HIS':
        aa1='H'
    else:
        aa1='-0-'

    return (aa1+'\n')

for i in str1:
    j=i.splitlines()[0]
    k=aa_abbreviate(j)
    abbrname.write(k)

aafulname.close()
abbrname.close()