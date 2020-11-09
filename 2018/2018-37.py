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


mytopicdb=myclient["cs2018_37"]
mydata=mytopicdb["cs2018_score_37"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2018_score_37_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_nonsmall = log((29138919 - 39891 + 0.5) / (39891 + 0.5), 10)
    idf_cell = log((29138919 - 5400577 + 0.5) / (5400577 + 0.5), 10)
    idf_carcinoma = log((29138919 - 494907 + 0.5) / (494907 + 0.5), 10)
    idf_nsclc = log((29138919 - 39891 + 0.5) / (39891 + 0.5), 10)
    idf_met = log((29138919 - 3359 + 0.5) / (3359 + 0.5), 10)



    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0+ 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 25413 + 0.5) / (25413 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 87973 + 0.5) / (87973 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 282359 + 0.5) / (282359+ 0.5), 10)
    idf_eleM_3 = log((25389659 - 25413 + 0.5) / (25413 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 24866 + 0.5) / (24866 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 19 + 0.5) / (19 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 519082 + 0.5) / (519082 + 0.5), 10)
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
        nonsmall_score=0
        cell_score = 0
        carcinoma_score = 0
        nsclc_score = 0
        met_score = 0
        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                nonsmall = [True for x in wordfreq.items() if 'non-small' in x]
                cell = [True for x in wordfreq.items() if 'cell' in x]
                carcinoma = [True for x in wordfreq.items() if 'carcinoma' in x]
                nsclc = [True for x in wordfreq.items() if 'nsclc' in x]
                # ---------------摘要统计-------------------#


                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]
                for key in wordfreq:
                    if 'non-small' in key:
                        nonsmall_score = nonsmall_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'cell' in key1:
                        cell_score = cell_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'carcinoma' in key1:
                        carcinoma_score = carcinoma_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'nsclc' in key1:
                        nsclc_score = nsclc_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'met' == key1:
                        met_score = met_score + wordfreq[key]



                bm25_nonsmall_score = (((k1+1)*nonsmall_score)/((k1*(b1+(1-b1)*(len_freq/85)))+nonsmall_score))
                bm25_cell_score = (((k1 + 1) * cell_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + cell_score))
                bm25_carcinoma_score = (((k1 + 1) * carcinoma_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + carcinoma_score))
                bm25_nsclc_score = (((k1 + 1) * nsclc_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + nsclc_score))
                bm25_met_score = (((k1 + 1) * met_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + met_score))

                bm25_ab_score =idf_nonsmall*bm25_nonsmall_score+idf_cell*bm25_cell_score+idf_carcinoma*bm25_carcinoma_score+idf_nsclc*bm25_nsclc_score+idf_met*bm25_met_score

                idf_para=[{str(nonsmall_score):idf_nonsmall},{str(cell_score):idf_cell},{str(carcinoma_score):idf_carcinoma},{str(nsclc_score):idf_nsclc},{str(met_score):idf_met}]

                # ---------------共现分析摘要-------------------#
                if len(nonsmall) != 0 and nonsmall[0] and len(cell) != 0 and cell[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'met' == key:
                            gx = idf_met
                if len(nsclc) != 0 and nsclc[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'met' == key:
                            gx = idf_met
            # ---------------共现分析化学-------------------#
                if len(nonsmall) != 0 and nonsmall[0] and len(cell) != 0 and cell[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for ele in ChemicalNameList:
                        if 'MET' in ele['NameOfSubstance']:
                            gx = idf_met
                            break
                if len(nsclc) != 0 and nsclc[0]:
                    for ele in ChemicalNameList:
                        if 'MET' in ele['NameOfSubstance']:
                            gx = idf_met
                            break
            # ---------------共现分析关键字-------------------#
                if len(nonsmall) != 0 and nonsmall[0] and len(cell) != 0 and cell[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for eleK in KeywordsList:
                        if 'met' in str(eleK).lower():
                            gx = idf_met
                            break
                if len(nsclc) != 0 and nsclc[0]:
                    for eleK in KeywordsList:
                        if 'met' in str(eleK).lower():
                            gx = idf_met
                            break
             # ---------------共现分析医学主题词-------------------#
                if len(nonsmall) != 0 and nonsmall[0] and len(cell) != 0 and cell[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for eleM in MeshHeadingNameList:
                        if 'MET' in eleM['MeshHeadingName']:
                            gx = idf_met
                            break
                if len(nsclc) != 0 and nsclc[0]:
                    for eleM in MeshHeadingNameList:
                        if 'MET' in eleM['MeshHeadingName']:
                            gx = idf_met
                            break

                for ele in ChemicalNameList:
                    if 'Carcinoma' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'Immunohistochemistry' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'Methionine' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break



                for eleM in MeshHeadingNameList:
                    if 'Carcinoma' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'Immunohistochemistry' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Methionine' == eleM['MeshHeadingName']:
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
                    if 'non-small cell carcinoma' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'met' in str(eleK).lower():
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
    sortsecond(mywords,mydata,10)
    count(mydata,mycount,"37")



