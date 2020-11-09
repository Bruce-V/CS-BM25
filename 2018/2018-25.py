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


mytopicdb=myclient["cs2018_25"]
mydata=mytopicdb["cs2018_score_25"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2018_score_25_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_melanoma = log((29138919 - 88340 + 0.5) / (88340 + 0.5), 10)
    idf_high = log((29138919 - 6347306 + 0.5) / (6347306 + 0.5), 10)
    idf_serum = log((29138919 - 986648 + 0.5) / (986648 + 0.5), 10)
    idf_ldh = log((29138919 - 27163 + 0.5) / (27163 + 0.5), 10)
    idf_level = log((29138919 - 1646688 + 0.5) / (1646688 + 0.5), 10)

    idf_ele_1 = log((13670358 - 40230 + 0.5) / (40230 + 0.5), 10)
    idf_ele_2 = log((13670358 - 80567 + 0.5) / (80567 + 0.5), 10)
    idf_ele_3 = log((13670358 - 1673 + 0.5) / (1673 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 88975 + 0.5) / (88975 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 40230 + 0.5) / (40230 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 113103 + 0.5) / (113103 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 80567 + 0.5) / (80567 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 0 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 11900 + 0.5) / (11900 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 24043 + 0.5) / (24043 + 0.5), 10)
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
        serum_score = 0
        ldh_score = 0
        level_score = 0
        high_score=0
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
                    if 'high' in key1:
                        high_score = high_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'serum' in key1:
                        serum_score = serum_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'ldh' in key1:
                        ldh_score = ldh_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'level' in key1:
                        level_score = level_score + wordfreq[key]

                bm25_melanoma_score = (((k1+1)*melanoma_score)/((k1*(b1+(1-b1)*(len_freq/85)))+melanoma_score))
                bm25_serum_score = (((k1 + 1) * serum_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + serum_score))
                bm25_ldh_score = (((k1 + 1) * ldh_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + ldh_score))
                bm25_level_score = (((k1 + 1) * level_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + level_score))
                bm25_high_score = (((k1 + 1) * high_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + high_score))

                bm25_ab_score =idf_melanoma*bm25_melanoma_score+idf_serum*bm25_serum_score+idf_ldh*bm25_ldh_score+idf_level*bm25_level_score+idf_high*bm25_high_score

                idf_para=[{str(melanoma_score):idf_melanoma},{str(serum_score):idf_serum},{str(ldh_score):idf_ldh},{str(level_score):idf_level},{str(high_score):idf_high}]

                # ---------------共现分析摘要-------------------#
                if len(melanoma)!=0 and melanoma[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'high' in key:
                            gx1 = idf_high
                        if 'serum' in key:
                            gx2 = idf_serum
                        if 'ldh' in key:
                            gx3 = idf_ldh
                        if 'level' in key:
                            gx4 = idf_level

            # ---------------共现分析化学-------------------#
                if len(melanoma) != 0 and melanoma[0]:
                    for ele in ChemicalNameList:
                        if 'LDH' in ele['NameOfSubstance']:
                            gx3 = idf_ldh
                            break
            # ---------------共现分析关键字-------------------#
                if len(melanoma) != 0 and melanoma[0]:
                    for eleK in KeywordsList:
                        if 'ldh' in str(eleK).lower():
                            gx3 = idf_ldh
                            break
             # ---------------共现分析医学主题词-------------------#
                if len(melanoma) != 0 and melanoma[0]:
                    for eleM in MeshHeadingNameList:
                        if 'LDH' in eleM['MeshHeadingName']:
                            gx3 = idf_ldh
                            break

                for ele in ChemicalNameList:
                    if 'L-Lactate Dehydrogenase' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'Isoenzymes' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'LDH' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break


                for eleM in MeshHeadingNameList:
                    if 'Melanoma' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'L-Lactate Dehydrogenase' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Skin Neoplasms' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    if 'Isoenzymes' in eleM['MeshHeadingName']:
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
                    if 'Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    if 'LDH' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_8
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
                    if 'ldh' in str(eleK).lower():
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
    sortsecond(mywords,mydata,12)
    count(mydata,mycount,"25")



