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


mytopicdb=myclient["cs2017_3"]
mydata=mytopicdb["cs2017_score_3"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2017_score_3_related"]#聚类后对应与主题相关联的文献






def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_meningioma = log((29138919 - 11153 + 0.5) / (11153 + 0.5), 10)
    idf_nf2 = log((29138919 - 1431 + 0.5) / (1431 + 0.5), 10)
    idf_akt1 = log((29138919 - 2258 + 0.5) / (2258 + 0.5), 10)
    idf_e17k = log((29138919 - 51 + 0.5) / (51 + 0.5), 10)
    idf_k322 = log((29138919 - 5 + 0.5) / (5 + 0.5), 10)

    idf_ele_1 = log((13670358 - 906 + 0.5) / (906 + 0.5), 10)
    idf_ele_2 = log((13670358 - 3691 + 0.5) / (3691 + 0.5), 10)
    idf_ele_3 = log((13670358 - 36676 + 0.5) / (36676 + 0.5), 10)
    idf_ele_4 = log((13670358 - 664 + 0.5) / (664 + 0.5), 10)
    idf_ele_5 = log((13670358 - 2056 + 0.5) / (2056 + 0.5), 10)
    idf_ele_6 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 664 + 0.5) / (664 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 15403 + 0.5) / (15403 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 18215 + 0.5) / (18215 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 36676 + 0.5) / (36676 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 2056 + 0.5) / (2056 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_12 = log((25389659 - 1582 + 0.5) / (1582 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 3125 + 0.5) / (3125 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 426 + 0.5) / (426 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 208 + 0.5) / (208 + 0.5), 10)
    idf_eleK_4 = log((5435471 - 3 + 0.5) / (3 + 0.5), 10)
    idf_eleK_5 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        meningioma_score = 0
        nf2_score = 0
        akt1_score = 0
        k322_score = 0
        e17k_score = 0
        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']

                meningioma = [True for x in wordfreq.items() if 'meningioma' in x]
                # ---------------摘要统计-------------------#
                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]

                for key in wordfreq:
                    if 'meningioma' in key:
                        meningioma_score = meningioma_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'nf2' in key1:
                        nf2_score = nf2_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'akt1' in key1:
                        akt1_score = akt1_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'k322' in key1:
                        k322_score = k322_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'e17k' in key1:
                        e17k_score = e17k_score + wordfreq[key]

                # ---------------共现分析摘要-------------------#
                if len(meningioma) != 0 and meningioma[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'k322' in key:
                                gx=idf_k322
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'nf2' in key:
                                gx2 =idf_nf2
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'e17k' in key:
                                gx1=idf_e17k
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for key in wordfreq:
                            key = cop.sub('', key)
                            if 'akt1' in key:
                                gx3 = idf_akt1
                                break

                # ---------------共现分析化学-------------------#
                if len(meningioma) != 0 and meningioma[0]:
                        for ele in ChemicalNameList:
                            if 'e17k' in ele['NameOfSubstance']:
                                gx1 = idf_e17k
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for ele in ChemicalNameList:
                            if 'k322' in ele['NameOfSubstance']:
                                gx = idf_k322
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for ele in ChemicalNameList:
                            if 'nf2' in ele['NameOfSubstance']:
                                gx2 = idf_nf2
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for ele in ChemicalNameList:
                            if 'akt1' in ele['NameOfSubstance']:
                                gx3 = idf_akt1
                                break


                # ---------------共现分析医学主题词-------------------#
                if len(meningioma) != 0 and meningioma[0]:
                        for eleM in MeshHeadingNameList:
                            if 'e17k' in eleM['MeshHeadingName']:
                                gx1 = idf_e17k
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for eleM in MeshHeadingNameList:
                            if 'k322' in eleM['MeshHeadingName']:
                                gx = idf_k322
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for eleM in MeshHeadingNameList:
                            if 'nf2' in eleM['MeshHeadingName']:
                                gx2 = idf_nf2
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for eleM in MeshHeadingNameList:
                            if 'akt1' in eleM['MeshHeadingName']:
                                gx3 = idf_akt1
                                break


                # ---------------共现分析关键字-------------------#
                if len(meningioma) != 0 and meningioma[0]:
                        for eleK in KeywordsList:
                            if 'e17k' in str(eleK).lower():
                                gx1 = idf_e17k
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for eleK in KeywordsList:
                            if 'k322' in str(eleK).lower():
                                gx = idf_k322
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for eleK in KeywordsList:
                            if 'nf2' in str(eleK).lower():
                                gx2 = idf_nf2
                                break
                if len(meningioma) != 0 and meningioma[0]:
                        for eleK in KeywordsList:
                            if 'akt1' in str(eleK).lower():
                                gx3 = idf_akt1
                                break

                bm25_meningioma_score = (((k1 + 1) * meningioma_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + meningioma_score))
                bm25_nf2_score = (((k1 + 1) * nf2_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + nf2_score))
                bm25_akt1_score = (((k1 + 1) * akt1_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + akt1_score))
                bm25_k322_score = (((k1 + 1) * k322_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + k322_score))
                bm25_e17k_score = (((k1 + 1) * e17k_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + e17k_score))

                bm25_ab_score = idf_meningioma * bm25_meningioma_score + idf_nf2 * bm25_nf2_score + idf_akt1 * bm25_akt1_score + idf_k322 * bm25_k322_score + idf_e17k * bm25_e17k_score

                idf_para = [{str(meningioma_score): idf_meningioma}, {str(nf2_score): idf_nf2},
                            {str(akt1_score): idf_akt1}, {str(k322_score): idf_k322}, {str(e17k_score): idf_e17k}]


                for ele in ChemicalNameList:
                    if 'NF2' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'AKT1' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'Proto-Oncogene Proteins c-akt' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break
                for ele in ChemicalNameList:
                    if 'Neurofibromin 2' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_4
                        break
                for ele in ChemicalNameList:
                    if 'Class I Phosphatidylinositol 3-Kinases' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_5
                        break
                for ele in ChemicalNameList:
                    if 'Neurofibromatosis 2' in ele['NameOfSubstance']:
                        ss2 = ss2 + idf_ele_6
                        break

                for eleM in MeshHeadingNameList:
                    if 'Neurofibromin 2' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'Meningeal Neoplasms' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'AKT1' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    if 'Meningioma' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_4
                        break
                for eleM in MeshHeadingNameList:
                    if 'Proto-Oncogene Proteins c-akt' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    if 'Class I Phosphatidylinositol 3-Kinases' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    if 'Female' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_8
                        break
                for eleM in MeshHeadingNameList:
                    if 'Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_9
                        break
                for eleM in MeshHeadingNameList:
                    if 'Middle Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_10
                        break

                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_11
                        break
                for eleM in MeshHeadingNameList:
                    if 'Neurofibromatosis 2' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_12
                        break

                for eleK in KeywordsList:
                    if 'meningioma' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'nf2' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_2
                        break
                for eleK in KeywordsList:
                    if 'akt1' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_3
                        break
                for eleK in KeywordsList:
                    if 'e17k' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_4
                        break
                for eleK in KeywordsList:
                    if 'k322' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_5
                        break

                total_gx = (gx + gx1 + gx2 + gx3)
                cmk_len = len(ChemicalNameList) + len(MeshHeadingNameList) + len(KeywordsList)
                bm25_cmk_len = ss1 + ss2 + ss4
                bm25_cmk_score=(((k2+1)*bm25_cmk_len)/((k2*(b2+(1-b2)*(cmk_len/13)))+bm25_cmk_len))
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
    sortsecond(mywords,mydata,8)
    count(mydata,mycount,"3")



