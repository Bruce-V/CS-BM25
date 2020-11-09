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


mytopicdb=myclient["cs2018_40"]
mydata=mytopicdb["cs2018_score_40"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2018_score_40_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_breast = log((29138919 - 222402 + 0.5) / (222402 + 0.5), 10)
    idf_erbb2 = log((29138919 - 3359 + 0.5) / (3359 + 0.5), 10)



    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0+ 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 8455 + 0.5) / (8455 + 0.5), 10)
    idf_ele_5 = log((13670358 - 22004 + 0.5) / (22004 + 0.5), 10)
    idf_ele_6 = log((13670358 - 22 + 0.5) / (22 + 0.5), 10)
    idf_ele_7 = log((13670358 - 24866 + 0.5) / (24866 + 0.5), 10)


    idf_eleM_1 = log((25389659 - 264316 + 0.5) / (264316 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 9290 + 0.5) / (9290+ 0.5), 10)
    idf_eleM_3 = log((25389659 - 7755 + 0.5) / (7755 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 24866 + 0.5) / (24866 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 22004 + 0.5) / (22004 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 56304 + 0.5) / (56304 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 330 + 0.5) / (330 + 0.5), 10)
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
        breast_score=0
        erbb2_score = 0

        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                breast = [True for x in wordfreq.items() if 'breast' in x]


                # ---------------摘要统计-------------------#


                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'breast' in key1:
                        breast_score = breast_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'erbb2' == key1:
                        erbb2_score = erbb2_score + wordfreq[key]




                bm25_breast_score = (((k1+1)*breast_score)/((k1*(b1+(1-b1)*(len_freq/85)))+breast_score))
                bm25_erbb2_score = (((k1 + 1) * erbb2_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + erbb2_score))


                bm25_ab_score =idf_breast*bm25_breast_score+idf_erbb2*bm25_erbb2_score

                idf_para=[{str(breast_score):idf_breast},{str(erbb2_score):idf_erbb2}]

                # ---------------共现分析摘要-------------------#
                if len(breast) != 0 and breast[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'erbb2' == key:
                            gx = idf_erbb2

            # ---------------共现分析化学-------------------#
                if len(breast) != 0 and breast[0]:
                    for ele in ChemicalNameList:
                        if 'ERBB2' in ele['NameOfSubstance']:
                            gx = idf_erbb2
                            break

            # ---------------共现分析关键字-------------------#
                if len(breast) != 0 and breast[0]:
                    for eleK in KeywordsList:
                        if 'erbb2' in str(eleK).lower():
                            gx = idf_erbb2
                            break

             # ---------------共现分析医学主题词-------------------#
                if len(breast) != 0 and breast[0]:
                    for eleM in MeshHeadingNameList:
                        if 'ERBB2' in eleM['MeshHeadingName']:
                            gx = idf_erbb2
                            break


                for ele in ChemicalNameList:
                    if 'Breast Neoplasms' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'Rare Diseases' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'Mastectomy, Segmental' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break
                for ele in ChemicalNameList:
                    if 'ERBB2 protein, human' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_4
                        break
                for ele in ChemicalNameList:
                    if 'Receptor, ErbB-2' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_5
                        break
                for ele in ChemicalNameList:
                    if 'herstatin' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_6
                        break
                for ele in ChemicalNameList:
                    if 'Intercellular Signaling Peptides and Proteins' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_7
                        break

                for eleM in MeshHeadingNameList:
                    if 'Breast Neoplasms' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'Rare Diseases' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Mastectomy, Segmental' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    if 'herstatin' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_4
                        break

                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    if 'Female' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    if 'Middle Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    if 'Intercellular Signaling Peptides and Proteins' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_8
                        break

                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_9
                        break
                for eleM in MeshHeadingNameList:
                    if 'ERBB2 protein, human' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_10
                        break
                for eleM in MeshHeadingNameList:
                    if 'Receptor, ErbB-2' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_11
                        break

                for eleK in KeywordsList:
                    if 'breast' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'erbb2' in str(eleK).lower():
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
    sortsecond(mywords,mydata,8)
    count(mydata,mycount,"40")



