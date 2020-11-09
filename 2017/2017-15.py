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


mytopicdb=myclient["cs2017_15"]
mydata=mytopicdb["cs2017_score_15"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2017_score_15_related"]#聚类后对应与主题相关联的文献






def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_cervical = log((29138919 - 195296 + 0.5) / (195296 + 0.5), 10)
    idf_stk11 = log((29138919 - 310 + 0.5) / (310 + 0.5), 10)

    idf_ele_1 = log((13670358 - 858 + 0.5) / (858 + 0.5), 10)
    idf_ele_2 = log((13670358 - 56231 + 0.5) / (56231 + 0.5), 10)
    idf_ele_3 = log((13670358 - 393 + 0.5) / (393 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 70343 + 0.5) / (70343 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 33210 + 0.5) / (33210 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 379394 + 0.5) / (379394 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 56231 + 0.5) / (56231 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 181173 + 0.5) / (181173 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 13064 + 0.5) / (13064 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 702001 + 0.5) / (702001 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 17047 + 0.5) / (17047 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 63 + 0.5) / (63 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq=0
        cervical_score = 0
        stk11_score = 0
        gx=0
        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                cervical = [True for x in wordfreq.items() if 'cervical' in x]
                # ---------------摘要统计-------------------#
                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]


                for key in wordfreq:
                    if 'cervical' in key:
                        # ss3 = ss3 + 3
                        cervical_score = cervical_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'stk11' in key1:
                        # ss3 = ss3 + 3
                        stk11_score = stk11_score + wordfreq[key]
                # ---------------共现分析摘要-------------------#
                if len(cervical) != 0 and cervical[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'stk11' in key:
                                gx=idf_stk11
                                break
            # ---------------共现分析化学-------------------#
                if len(cervical) != 0 and cervical[0]:
                        for ele in ChemicalNameList:
                            if 'STK11' in ele['NameOfSubstance']:
                                gx = idf_stk11
                                break
            # ---------------共现分析关键字-------------------#
                if len(cervical) != 0 and cervical[0]:
                        for eleK in KeywordsList:
                            if 'stk11' in str(eleK).lower():
                                gx = idf_stk11
                                break
                # ---------------共现分析医学主题词-------------------#
                if len(cervical) != 0 and cervical[0]:
                        for eleM in MeshHeadingNameList:
                            if 'STK11' in eleM['MeshHeadingName']:
                                gx = idf_stk11
                                break

                bm25_cervical_score = (((k1 + 1) * cervical_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + cervical_score))
                bm25_stk11_score = (((k1 + 1) * stk11_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + stk11_score))

                bm25_ab_score = idf_cervical * bm25_cervical_score + idf_stk11 * bm25_stk11_score

                idf_para = [{str(cervical_score): idf_cervical}, {str(stk11_score): idf_stk11}]

                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'STK11' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'Protein-Serine-Threonine Kinases' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'Phosphorylation' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break
                for ele in ChemicalNameList:
                    # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
                    if 'Proteolysis' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_4
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'STK11' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Uterine Cervical Neoplasms' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Cervical Vertebrae' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Radiography' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_4
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Protein-Serine-Threonine Kinases' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Phosphorylation' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Proteolysis' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_8
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Female' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_9
                        break
                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'Young' in eleM['MeshHeadingName']:
                        ss2 = ss2 +idf_eleM_10
                        break

                for eleM in MeshHeadingNameList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_11
                        break

                for eleK in KeywordsList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'cervical' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
                    if 'stk11' in str(eleK).lower():
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
    sortsecond(mywords,mydata,5.5)
    count(mydata,mycount,"15")


