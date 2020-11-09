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


mytopicdb=myclient["cs2018_39"]
mydata=mytopicdb["cs2018_score_39"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2018_score_39_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_anaplastic = log((29138919 - 9173 + 0.5) / (9173 + 0.5), 10)
    idf_large = log((29138919 - 2311 + 0.5) / (2311 + 0.5), 10)
    idf_cell = log((29138919 - 9173 + 0.5) / (9173 + 0.5), 10)
    idf_lymphoma = log((29138919 - 2311 + 0.5) / (2311 + 0.5), 10)
    idf_alcl=log((29138919 - 2311 + 0.5) / (2311 + 0.5), 10)
    idf_alk=log((29138919 - 2311 + 0.5) / (2311 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 16271+ 0.5) / (16271 + 0.5), 10)
    idf_ele_3 = log((13670358 - 2615 + 0.5) / (2615 + 0.5), 10)
    idf_ele_4 = log((13670358 - 33948 + 0.5) / (33948 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 1786 + 0.5) / (1786 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 16271 + 0.5) / (16271+ 0.5), 10)
    idf_eleM_3 = log((25389659 - 2615 + 0.5) / (2615 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 33948 + 0.5) / (33948 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 1862883 + 0.5) / (1862883 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 35956 + 0.5) / (35956 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 226 + 0.5) / (226 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 29226 + 0.5) / (29226 + 0.5), 10)
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
        anaplastic_score=0
        large_score = 0
        cell_score = 0
        lymphoma_score = 0
        alcl_score = 0
        alk_score=0
        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                anaplastic = [True for x in wordfreq.items() if 'anaplastic' in x]
                large = [True for x in wordfreq.items() if 'large' in x]
                cell = [True for x in wordfreq.items() if 'cell' in x]
                lymphoma = [True for x in wordfreq.items() if 'lymphoma' in x]
                alcl = [True for x in wordfreq.items() if 'cell' in x]

                # ---------------摘要统计-------------------#


                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'anaplastic' in key1:
                        anaplastic_score = anaplastic_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'large' in key1:
                        large_score = large_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'cell' in key1:
                        cell_score = cell_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'lymphoma' in key1:
                        lymphoma_score = lymphoma_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'alcl' in key1:
                        alcl_score = alcl_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'alk' == key1:
                        alk_score = alk_score + wordfreq[key]




                bm25_anaplastic_score = (((k1+1)*anaplastic_score)/((k1*(b1+(1-b1)*(len_freq/85)))+anaplastic_score))
                bm25_large_score = (((k1 + 1) * large_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + large_score))
                bm25_cell_score = (((k1 + 1) * cell_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + cell_score))
                bm25_lymphoma_score = (((k1 + 1) * lymphoma_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + lymphoma_score))
                bm25_alcl_score = (((k1 + 1) * alcl_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + alcl_score))
                bm25_alk_score = (((k1 + 1) * alk_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + alk_score))

                bm25_ab_score =idf_anaplastic*bm25_anaplastic_score+idf_large*bm25_large_score+idf_cell*bm25_cell_score+idf_lymphoma*bm25_lymphoma_score+idf_alcl*bm25_alcl_score+idf_alk*bm25_alk_score

                idf_para=[{str(anaplastic_score):idf_anaplastic},{str(large_score):idf_large},{str(cell_score):idf_cell},{str(lymphoma_score):idf_lymphoma},{str(alcl_score):idf_alcl},{str(alk_score):idf_alk}]

                # ---------------共现分析摘要-------------------#
                if len(anaplastic) != 0 and anaplastic[0] and len(large) != 0 and large[0] and len(cell) != 0 and cell[0] and len(lymphoma) != 0 and lymphoma[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'alk' == key:
                            gx = idf_alk
                if len(alcl) != 0 and alcl[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'alk' == key:
                            gx = idf_alk
            # ---------------共现分析化学-------------------#
                if len(anaplastic) != 0 and anaplastic[0] and len(large) != 0 and large[0] and len(cell) != 0 and cell[0] and len(lymphoma) != 0 and lymphoma[0]:
                    for ele in ChemicalNameList:
                        if 'ALK' in ele['NameOfSubstance']:
                            gx = idf_alk
                            break
                if len(alcl) != 0 and alcl[0]:
                    for ele in ChemicalNameList:
                        if 'ALK' in ele['NameOfSubstance']:
                            gx = idf_alk
                            break
            # ---------------共现分析关键字-------------------#
                if len(anaplastic) != 0 and anaplastic[0] and len(large) != 0 and large[0] and len(cell) != 0 and cell[0] and len(lymphoma) != 0 and lymphoma[0]:
                    for eleK in KeywordsList:
                        if 'alk' in str(eleK).lower():
                            gx = idf_alk
                            break
                if len(alcl) != 0 and alcl[0]:
                    for eleK in KeywordsList:
                        if 'alk' in str(eleK).lower():
                            gx = idf_alk
                            break
             # ---------------共现分析医学主题词-------------------#
                if len(anaplastic) != 0 and anaplastic[0] and len(large) != 0 and large[0] and len(cell) != 0 and cell[0] and len(lymphoma) != 0 and lymphoma[0]:
                    for eleM in MeshHeadingNameList:
                        if 'ALK' in eleM['MeshHeadingName']:
                            gx = idf_alk
                            break
                if len(alcl) != 0 and alcl[0]:
                    for eleM in MeshHeadingNameList:
                        if 'ALK' in eleM['MeshHeadingName']:
                            gx = idf_alk
                            break

                for ele in ChemicalNameList:
                    if 'Lymphoma, Large-Cell, Anaplastic' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'Receptor Protein-Tyrosine Kinases' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'Anaplastic Lymphoma Kinase' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break
                for ele in ChemicalNameList:
                    if 'Protein-Tyrosine Kinases' == ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_4
                        break



                for eleM in MeshHeadingNameList:
                    if 'Lymphoma, Large-Cell, Anaplastic' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'Receptor Protein-Tyrosine Kinases' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Anaplastic Lymphoma Kinase' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    if 'Protein-Tyrosine Kinases' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_4
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
                    if 'Child' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break


                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_9
                        break


                for eleK in KeywordsList:
                    if 'anaplastic large cell lymphoma' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'alk' in str(eleK).lower():
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
    sortsecond(mywords,mydata,20)
    count(mydata,mycount,"39")



