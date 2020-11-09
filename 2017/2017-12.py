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
mytopic=mydb["topics2017"]#pubmed中的主题词相关文献列表
mypapers=mydb["papers"]#pubmed中文献信息表


mytopicdb=myclient["cs2017_12"]
mydata=mytopicdb["cs2017_score_12"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2017_score_12_related"]#聚类后对应与主题相关联的文献






def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_colon = log((29138919 - 133274 + 0.5) / (133274 + 0.5), 10)
    idf_v600e = log((29138919 - 1887 + 0.5) / (1887 + 0.5), 10)
    idf_braf = log((29138919 - 8527 + 0.5) / (8527 + 0.5), 10)

    idf_ele_1 = log((13670358 - 5411 + 0.5) / (5411 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 7222 + 0.5) / (7222 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 67366 + 0.5) / (67366 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 16004324 + 0.5) / (16004324 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 702001 + 0.5) / (702001 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)


    idf_eleK_1 = log((5435471 - 3978 + 0.5) / (3978 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 7526 + 0.5) / (7526 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 415 + 0.5) / (415 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq = 0
        colon_score = 0
        v600e_score = 0
        braf_score = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                colon = [True for x in wordfreq.items() if 'colon' in x]
                # ---------------摘要统计-------------------#
                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]
                for key in wordfreq:
                    if 'colon' in key:
                        colon_score = colon_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'braf' in key1:
                        braf_score = braf_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'v600e' in key1:
                        v600e_score = v600e_score + wordfreq[key]


                #---------------共现分析摘要-------------------#
                if len(colon) != 0 and colon[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'v600e' in key:
                                gx=idf_v600e
                                break
                if len(colon) != 0 and colon[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'braf' in key:
                                gx2 = idf_braf
                                break

                bm25_colon_score = (((k1 + 1) * colon_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + colon_score))
                bm25_v600e_score = (((k1 + 1) * v600e_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + v600e_score))
                bm25_braf_score = (((k1 + 1) * braf_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + braf_score))


                bm25_ab_score = idf_colon * bm25_colon_score + idf_v600e * bm25_v600e_score +  idf_braf * bm25_braf_score

                idf_para = [{str(colon_score): idf_colon},
                            {str(v600e_score): idf_v600e}, {str(braf_score): idf_braf}]

                # ---------------共现分析化学-------------------#
                if len(colon) != 0 and colon[0]:
                        for ele in ChemicalNameList:
                            if 'V600E' in ele['NameOfSubstance']:
                                gx = idf_v600e
                                break
                if len(colon) != 0 and colon[0]:
                        for ele in ChemicalNameList:
                            if 'BRAF' in ele['NameOfSubstance']:
                                gx2 = idf_braf
                                break

                # ---------------共现分析医学主题词-------------------#
                if len(colon)!=0 and colon[0]:
                        for eleM in MeshHeadingNameList:
                            if 'V600E' in eleM['MeshHeadingName']:
                                gx = idf_v600e
                                break

                if len(colon) != 0 and colon[0]:
                        for eleM in MeshHeadingNameList:
                            if 'BRAF' in eleM['MeshHeadingName']:
                                gx2 = idf_braf
                                break


                # ---------------共现分析关键字-------------------#
                if len(colon) != 0 and colon[0]:
                        for eleK in KeywordsList:
                            if 'v600e' in str(eleK).lower():
                                gx = idf_v600e
                                break

                if len(colon) != 0 and colon[0]:
                        for eleK in KeywordsList:
                            if 'braf' in str(eleK).lower():
                                gx2 = idf_braf
                                break

                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'BRAF' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'V600E' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'Proto-Oncogene Proteins B-raf' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break

                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Colonic Neoplasms' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'V600E' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'BRAF' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_4
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Male' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Middle Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Young' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_8
                        break

                for eleK in KeywordsList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'colon cancer' in str(eleK).lower():
                        ss4 = ss4 +idf_eleK_1
                        break
                for eleK in KeywordsList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'braf' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_2
                        break
                for eleK in KeywordsList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'v600e' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_3
                        break

                total_gx = (gx + gx1 + gx2 + gx3)
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
    sortsecond(mywords,mydata,5.8)
    count(mydata,mycount,"12")


