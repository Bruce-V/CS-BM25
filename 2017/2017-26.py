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
from math import  log

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords3"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2017"]#pubmed中的主题词相关文献列表
mypapers=mydb["papers"]#pubmed中文献信息表


mytopicdb=myclient["cs2017_26"]
mydata=mytopicdb["cs2017_score_26"]#按词表长度改进过后的2次排序表

mycount = mytopicdb["cs2017_score_26_related"]#聚类后对应与主题相关联的文献






def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_breast = log((29138919 - 357041 + 0.5) / (357041 + 0.5), 10)
    idf_nras = log((29138919 - 1657 + 0.5) / (1657 + 0.5), 10)
    idf_amplification = log((29138919 - 110630 + 0.5) / (110630 + 0.5), 10)

    idf_ele_1 = log((13670358 - 629 + 0.5) / (629 + 0.5), 10)
    idf_ele_2 = log((13670358 - 8363 + 0.5) / (8363 + 0.5), 10)
    idf_ele_3 = log((13670358 - 7222 + 0.5) / (7222 + 0.5), 10)
    idf_ele_4 = log((13670358 - 2398 + 0.5) / (2398 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 268832 + 0.5) / (268832 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 8363 + 0.5) / (8363 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 7221 + 0.5) / (7221 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 88975 + 0.5) / (88975 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 702001 + 0.5) / (702001 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 31674 + 0.5) / (31674 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 249 + 0.5) / (249 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq=0
        breast_score = 0
        nras_score = 0
        amplification_score = 0
        gx=0
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
                    if 'breast' in key:
                        # ss3 = ss3 + 3
                        breast_score = breast_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'nras' in key1:
                        # ss3 = ss3 + 3
                        nras_score = nras_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'amplification' in key1:
                        # ss3 = ss3 + 3
                        amplification_score = amplification_score + wordfreq[key]
                # ---------------共现分析摘要-------------------#
                if len(breast) != 0 and breast[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'nras' in key:
                                gx=idf_nras
                                break
            # ---------------共现分析化学-------------------#
                if len(breast) != 0 and breast[0]:
                        for ele in ChemicalNameList:
                            if 'NRAS' in ele['NameOfSubstance']:
                                gx = idf_nras
                                break
            # ---------------共现分析关键字-------------------#
                if len(breast) != 0 and breast[0]:
                        for eleK in KeywordsList:
                            if 'nras' in str(eleK).lower():
                                gx = idf_nras
                                break
                # ---------------共现分析医学主题词-------------------#
                if len(breast) != 0 and breast[0]:
                        for eleM in MeshHeadingNameList:
                            if 'NRAS' in eleM['MeshHeadingName']:
                                gx = idf_nras
                                break

                bm25_breast_score = (((k1 + 1) * breast_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + breast_score))
                bm25_nras_score = (((k1 + 1) * nras_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + nras_score))
                bm25_amplification_score = (((k1 + 1) * amplification_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + amplification_score))


                bm25_ab_score = idf_breast * bm25_breast_score + idf_nras * bm25_nras_score +idf_amplification * bm25_amplification_score

                idf_para = [{str(breast_score): idf_breast}, {str(nras_score): idf_nras},{str(amplification_score): idf_amplification}]
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'NRAS' in ele['NameOfSubstance']:
                        ss1 = ss1 +idf_ele_1
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'GTP Phosphohydrolases' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'Proto-Oncogene Proteins B-raf' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'Melanoma' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_4
                        break

                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'NRAS' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Breast Neoplasms' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'GTP Phosphohydrolases' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Proto-Oncogene Proteins B-raf' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_4
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Melanoma' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 +idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Female' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Middle Aged' in eleM['MeshHeadingName']:
                        ss2 = ss2 +idf_eleM_8
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Young' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_9
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 +idf_eleM_10
                        break
                for eleK in KeywordsList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'breast cancer' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'nras' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_2
                        break
                total_gx = gx
                cmk_len = len(ChemicalNameList) + len(MeshHeadingNameList) + len(KeywordsList)
                bm25_cmk_len = ss1 + ss2 + ss4
                bm25_cmk_score = (((k2 + 1) * bm25_cmk_len) / ((k2 * (b2 + (1 - b2) * (cmk_len / 13))) + bm25_cmk_len))
                bm25_score = bm25_ab_score + bm25_cmk_score + total_gx

                if (bm25_score > yuzhi):
                    mydict = {"PMID": x['PMID'], "ab_score": bm25_ab_score, "idf_para": idf_para,
                              "cmk_len": cmk_len, "cmk_freq": bm25_cmk_len, "bm25_cmk_score": bm25_cmk_score,
                              "gx": total_gx,
                              "bm25_score": bm25_score,
                              "ChemicalNameList": x['ChemicalNameList'],
                              "MeshHeadingNameList": x['MeshHeadingNameList'], "KeywordsList": x['KeywordsList']}
                    y = mydata.insert_one(mydict)
                    k = k + 1
                    print(str(y) + '---------' + str(k))

def count(mysort,mycount,topic):
    for x in mysort.find({},
                         {'PMID', 'ab_score', 'idf_para', 'cmk_len', 'cmk_freq', 'bm25_cmk_score', 'gx', 'bm25_score',
                          'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}):
        kk = 0
        for y in mytopic.find({"topic": topic}, {'PMID', 'relate'}):
            if x['PMID'] == y['PMID']:
                mydict = {"PMID": x['PMID'], "related": y['relate'], "ab_score": x["ab_score"],
                          "idf_para": x['idf_para'],
                          "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'], 'bm25_cmk_score': x['bm25_cmk_score'],
                          'gx': x['gx'],
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
    sortsecond(mywords,mydata,7)
    count(mydata,mycount,"26")
