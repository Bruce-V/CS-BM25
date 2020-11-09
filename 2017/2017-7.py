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


mytopicdb=myclient["cs2017_7"]
mydata=mytopicdb["cs2017_score_7"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2017_score_7_related"]#聚类后对应与主题相关联的文献






def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_lung = log((29138919 - 494482 + 0.5) / (494482 + 0.5), 10)
    idf_egfr = log((29138919 - 36925 + 0.5) / (36925 + 0.5), 10)
    idf_l858r = log((29138919 - 805 + 0.5) / (805 + 0.5), 10)

    idf_ele_1 = log((13670358 - 7893 + 0.5) / (7893 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 37086 + 0.5) / (37086 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 200349 + 0.5) / (200349 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 43476 + 0.5) / (43476 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 491273 + 0.5) / (491273 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 46148 + 0.5) / (46148 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 5875 + 0.5) / (5875 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_12 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 67747 + 0.5) / (67747 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq = 0
        lung_score = 0
        egfr_score = 0
        l858r_score = 0
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

                lung = [True for x in wordfreq.items() if 'lung' in x]
                # ---------------摘要统计-------------------#
                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]

                for key in wordfreq:
                    if 'lung' in key:
                        lung_score = lung_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'egfr' in key1:
                        egfr_score = egfr_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'l858r' in key1:
                        l858r_score = l858r_score + wordfreq[key]


                #---------------共现分析摘要-------------------#
                if len(lung) != 0 and lung[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'l858r' in key:
                                gx=idf_l858r
                                break
                if len(lung) != 0 and lung[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'egfr' in key:
                                gx2 = idf_egfr
                                break


                # ---------------共现分析化学-------------------#
                if len(lung) != 0 and lung[0]:
                        for ele in ChemicalNameList:
                            if 'L858R' in ele['NameOfSubstance']:
                                gx = idf_l858r
                                break
                if len(lung) != 0 and lung[0]:
                        for ele in ChemicalNameList:
                            if 'EGFR' in ele['NameOfSubstance']:
                                gx2 = idf_egfr
                                break

                # ---------------共现分析医学主题词-------------------#
                if len(lung) != 0 and lung[0]:
                        for eleM in MeshHeadingNameList:
                            if 'L858R' in eleM['MeshHeadingName']:
                                gx = idf_l858r
                                break

                if len(lung) != 0 and lung[0]:
                        for eleM in MeshHeadingNameList:
                            if 'EGFR' in eleM['MeshHeadingName']:
                                gx2 = idf_egfr
                                break


                # ---------------共现分析关键字-------------------#
                if len(lung) != 0 and lung[0]:
                        for eleK in KeywordsList:
                            if 'l858r' in str(eleK).lower():
                                gx = idf_l858r
                                break

                if len(lung) != 0 and lung[0]:
                        for eleK in KeywordsList:
                            if 'egfr' in str(eleK).lower():
                                gx2 = idf_egfr
                                break

                bm25_lung_score = (((k1 + 1) * lung_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + lung_score))
                bm25_egfr_score = (((k1 + 1) * egfr_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + egfr_score))
                bm25_l858r_score = (((k1 + 1) * l858r_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + l858r_score))

                bm25_ab_score = idf_lung * bm25_lung_score + idf_egfr * bm25_egfr_score + idf_l858r * bm25_l858r_score

                idf_para = [{str(lung_score): idf_lung}, {str(egfr_score): idf_egfr},
                            {str(l858r_score): idf_l858r}]

                for ele in ChemicalNameList:
                    if 'EGFR' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'L858R' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'ErbB Receptors' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break
                for ele in ChemicalNameList:
                    if 'Non-Small-Cell Lung' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_4
                        break
                for ele in ChemicalNameList:
                    if 'Carcinoma' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_5
                        break

                for eleM in MeshHeadingNameList:
                    if 'Lung Neoplasms' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'ErbB Receptors' in eleM['MeshHeadingName']:
                        ss2 = ss2 +idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Carcinoma' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    if 'Non-Small-Cell Lun' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_4
                        break
                for eleM in MeshHeadingNameList:
                    if 'EGFR' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    if 'L858R' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    if 'Lupus Nephritis' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_8
                        break
                for eleM in MeshHeadingNameList:
                    if 'Female' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_9
                        break
                for eleM in MeshHeadingNameList:
                    if 'Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_10
                        break
                for eleM in MeshHeadingNameList:
                    if 'Middle Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_11
                        break
                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_12
                        break
                for eleK in KeywordsList:
                    if 'lung' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'egfr' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_2
                        break
                for eleK in KeywordsList:
                    if 'l858r' in str(eleK).lower():
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
     sortsecond(mywords,mydata,7)
     count(mydata,mycount,"7")


