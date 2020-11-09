# Copyright 2020 zicheng Zhang(18551701375@163.com)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords3"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2018"]#pubmed中的主题词相关文献列表
mypapers=mydb["papers"]#pubmed中文献信息表


mytopicdb=myclient["cs2018_19"]
mydata=mytopicdb["cs2018_score_19"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2018_score_19_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_melanoma = log((29138919 - 88340 + 0.5) / (88340 + 0.5), 10)
    idf_tumor = log((29138919 - 1343849 + 0.5) / (1343849 + 0.5), 10)
    idf_cell = log((29138919 - 5400577 + 0.5) / (5400577 + 0.5), 10)
    idf_negative = log((29138919 - 1343849 + 0.5) / (1343849 + 0.5), 10)
    idf_pdl1 = log((29138919 - 6813 + 0.5) / (6813 + 0.5), 10)
    idf_expression = log((29138919 - 1909989 + 0.5) / (1909989 + 0.5), 10)


    idf_ele_1 = log((13670358 - 2096 + 0.5) / (2096 + 0.5), 10)
    idf_ele_2 = log((13670358 - 3355 + 0.5) / (3355 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 3355 + 0.5) / (3355 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 113103 + 0.5) / (113103 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 63159 + 0.5) / (63159 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 18848 + 0.5) / (18848 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 11900 + 0.5) / (11900 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 1818 + 0.5) / (1818 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        gx4=0
        len_freq=0
        melanoma_score=0
        tumor_score = 0
        cell_score = 0
        negative_score = 0
        pdl1_score = 0
        expression_score = 0

        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                melanoma = [True for x in wordfreq.items() if 'melanoma' in x]


                # ---------------摘要统计-------------------#


                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'melanoma' in key1:
                        melanoma_score = melanoma_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'tumor' in key1:
                        tumor_score = tumor_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'cell' in key1:
                        cell_score = cell_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'negative' in key1:
                        negative_score = negative_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'pdl1' in key1:
                        pdl1_score = pdl1_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'expression' in key1:
                        expression_score = expression_score + wordfreq[key]


                bm25_melanoma_score = (((k1+1)*melanoma_score)/((k1*(b1+(1-b1)*(len_freq/85)))+melanoma_score))
                bm25_tumor_score =(((k1+1)*tumor_score)/((k1*(b1+(1-b1)*(len_freq/85)))+tumor_score))
                bm25_cell_score = (((k1+1)*cell_score)/((k1*(b1+(1-b1)*(len_freq/85)))+cell_score))
                bm25_negative_score = (((k1 + 1) * negative_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + negative_score))
                bm25_pdl1_score = (((k1 + 1) * pdl1_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + pdl1_score))
                bm25_expression_score = (((k1 + 1) * expression_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + expression_score))


                bm25_ab_score =idf_melanoma*bm25_melanoma_score+idf_tumor*bm25_tumor_score+idf_cell*bm25_cell_score+idf_negative*bm25_negative_score+idf_pdl1*bm25_pdl1_score+idf_expression*bm25_expression_score

                idf_para=[{str(melanoma_score):idf_melanoma},{str(tumor_score):idf_tumor},{str(cell_score):idf_cell},{str(negative_score):idf_negative},{str(pdl1_score):idf_pdl1},{str(expression_score):idf_expression}]

                # ---------------共现分析摘要-------------------#
                if len(melanoma)!=0 and melanoma[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'tumor' in key:
                            gx = idf_tumor
                        if 'cell' in key:
                            gx1 = idf_cell
                        if 'negative' in key:
                            gx2 = idf_negative
                        if 'pdl1' in key:
                            gx3 = idf_pdl1
                        if 'expression' in key:
                            gx4 = idf_expression

            # ---------------共现分析化学-------------------#

                if len(melanoma) != 0 and melanoma[0]:
                    for ele in ChemicalNameList:
                        if 'PD-L1' in ele['NameOfSubstance']:
                            gx = idf_pdl1
                            break

            # ---------------共现分析关键字-------------------#

                if len(melanoma) != 0 and melanoma[0]:
                    for eleK in KeywordsList:
                        if 'pd-l1' in str(eleK).lower():
                            gx = idf_pdl1
                            break


             # ---------------共现分析医学主题词-------------------#

                if len(melanoma) != 0 and melanoma[0]:
                    for eleM in MeshHeadingNameList:
                        if 'PD-L1' in eleM['MeshHeadingName']:
                            gx = idf_pdl1
                            break


                for ele in ChemicalNameList:
                    if 'CD274 protein, human' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'B7-H1 Antigen' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break




                for eleM in MeshHeadingNameList:
                    if 'CD274 protein, human' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'B7-H1 Antigen' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Skin Neoplasms' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break


                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    if 'Male' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    if 'Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break


                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_9
                        break


                for eleK in KeywordsList:
                    if 'melanoma' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'pd-l1' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_2
                        break


                total_gx=gx1+gx2+gx3+gx+gx4
                cmk_len=len(ChemicalNameList) + len(MeshHeadingNameList) + len(KeywordsList)
                bm25_cmk_len=ss1 + ss2 + ss4
                bm25_cmk_score = (((k2 + 1) * bm25_cmk_len) / ((k2 * (b2 + (1 - b2) * (cmk_len / 13))) + bm25_cmk_len))
                bm25_score=bm25_ab_score+bm25_cmk_score+total_gx
                if(bm25_score>yuzhi):
                    mydict = {"PMID": x['PMID'],"ab_score":bm25_ab_score,"idf_para":idf_para,
                              "cmk_len":cmk_len,"cmk_freq":bm25_cmk_len,"bm25_cmk_score":bm25_cmk_score,"gx":total_gx,"bm25_score":bm25_score,
                               "ChemicalNameList":x['ChemicalNameList'],"MeshHeadingNameList":x['MeshHeadingNameList'],"KeywordsList":x['KeywordsList']}
                    y = mydata.insert_one(mydict)
                    k=k+1
                    print(str(y) + '---------' + str(k))

def count(mysort,mycount,topic):
    for x in mysort.find({}, {'PMID', 'ab_score','idf_para', 'cmk_len', 'cmk_freq', 'bm25_cmk_score','gx','bm25_score',
                              'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}):
        kk = 0
        for y in mytopic.find({"topic": topic}, {'PMID', 'relate'}):
            if x['PMID'] == y['PMID']:
                mydict = {"PMID": x['PMID'], "related": y['relate'], "ab_score":x["ab_score"],"idf_para":x['idf_para'],
                          "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'],'bm25_cmk_score':x['bm25_cmk_score'],'gx':x['gx'],
                          "bm25_score": x['bm25_score'],
                          "ChemicalNameList": x['ChemicalNameList'], "MeshHeadingNameList": x['MeshHeadingNameList'],
                          "KeywordsList": x['KeywordsList']}
                ss = mycount.insert_one(mydict)
                print(ss)
                kk = kk + 1
        if (kk == 0):
            mydict = {"PMID": x['PMID'], "related": -1, "ab_score": x["ab_score"], "idf_para": x['idf_para'],
                      "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'], 'bm25_cmk_score': x['bm25_cmk_score'],
                      'gx': x['gx'],
                      "bm25_score": x['bm25_score'],
                      "ChemicalNameList": x['ChemicalNameList'], "MeshHeadingNameList": x['MeshHeadingNameList'],
                      "KeywordsList": x['KeywordsList']}
            ss = mycount.insert_one(mydict)
            print(ss)

if __name__ == '__main__':
    sortsecond(mywords,mydata,11)
    count(mydata,mycount,"19")



