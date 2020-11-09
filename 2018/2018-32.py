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


mytopicdb=myclient["cs2018_32"]
mydata=mytopicdb["cs2018_score_32"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2018_score_32_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_leukemia = log((29138919 - 193081 + 0.5) / (193081 + 0.5), 10)
    idf_abl1 = log((29138919 - 660 + 0.5) / (660 + 0.5), 10)



    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 6082 + 0.5) / (6082 + 0.5), 10)
    idf_ele_4 = log((13670358 - 649 + 0.5) / (649 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 56273 + 0.5) / (56273 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 1195 + 0.5) / (1195+ 0.5), 10)
    idf_eleM_3 = log((25389659 - 6082 + 0.5) / (6082 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 865241 + 0.5) / (865241 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 23285 + 0.5) / (23285 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 228 + 0.5) / (228 + 0.5), 10)
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
        leukemia_score=0
        abl1_score = 0

        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                leukemia = [True for x in wordfreq.items() if 'leukemia' in x]


                # ---------------摘要统计-------------------#


                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'leukemia' in key1:
                        leukemia_score = leukemia_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'abl1' == key1:
                        abl1_score = abl1_score + wordfreq[key]




                bm25_leukemia_score = (((k1+1)*leukemia_score)/((k1*(b1+(1-b1)*(len_freq/85)))+leukemia_score))
                bm25_abl1_score = (((k1 + 1) * abl1_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + abl1_score))


                bm25_ab_score =idf_leukemia*bm25_leukemia_score+idf_abl1*bm25_abl1_score

                idf_para=[{str(leukemia_score):idf_leukemia},{str(abl1_score):idf_abl1}]

                # ---------------共现分析摘要-------------------#
                if len(leukemia) != 0 and leukemia[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'abl1' == key:
                            gx = idf_abl1

            # ---------------共现分析化学-------------------#
                if len(leukemia) != 0 and leukemia[0]:
                    for ele in ChemicalNameList:
                        if 'abl' in ele['NameOfSubstance']:
                            gx = idf_abl1
                            break

            # ---------------共现分析关键字-------------------#
                if len(leukemia) != 0 and leukemia[0]:
                    for eleK in KeywordsList:
                        if 'abl1' in str(eleK).lower():
                            gx = idf_abl1
                            break

             # ---------------共现分析医学主题词-------------------#
                if len(leukemia) != 0 and leukemia[0]:
                    for eleM in MeshHeadingNameList:
                        if 'abl' in eleM['MeshHeadingName']:
                            gx = idf_abl1
                            break


                for ele in ChemicalNameList:
                    if 'Leukemia' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'Genes, abl' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'Fusion Proteins, bcr-abl' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break




                for eleM in MeshHeadingNameList:
                    if 'Leukemia' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'Genes, abl' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Fusion Proteins, bcr-abl' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
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
                    if 'Preschool' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break


                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_9
                        break


                for eleK in KeywordsList:
                    if 'leukemia' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'abl1' in str(eleK).lower():
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
    count(mydata,mycount,"32")



