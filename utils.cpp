﻿#include <QCoreApplication>
#include <QDir>
#include <iostream>
#include "include/hiredis/hiredis.h"
#include "neuron_editing/neuron_format_converter.h"
#include "utils.h"
#include "include/hiredis/hiredis.h"
void dirCheck(QString dirBaseName)
{
    //QCoreApplication::applicationDirPath()获取程序所在目录
    if(!QDir(QCoreApplication::applicationDirPath()+"/"+dirBaseName).exists())
    {
        QDir(QCoreApplication::applicationDirPath()).mkdir(dirBaseName);
    }
}

QStringList getSwcInBlock(const QString msg,const V_NeuronSWC_list& testVNL)
{
    /*
     * p1:brain_id;res;x;y;z;size;socket.descriptor
     * p2:Neuron_id/name
     * 返回：文件名，文件路径
     */
    QStringList paraList=msg.split(";",Qt::SkipEmptyParts);
    QString brain_id=paraList.at(0).trimmed();//1. tf name/RES  2. .v3draw// test:17302;RES;x;y;z;b
    //0: 18465/RESx18000x13000x5150
    //1: 12520
    //2: 7000
    //3: 2916
    int res=paraList.at(1).toInt();//>0
    int xpos=paraList.at(2).toInt();
    int ypos=paraList.at(3).toInt();
    int zpos=paraList.at(4).toInt();
    int blocksize=paraList.at(5).toInt();

    {
        QString name=brain_id;
        int x1=xpos-blocksize;
        int x2=xpos+blocksize;
        int y1=ypos-blocksize;
        int y2=ypos+blocksize;
        int z1=zpos-blocksize;
        int z2=zpos+blocksize;
        int cnt=pow(2,res-1);

        dirCheck("tmp");
        QString BBSWCNAME=QCoreApplication::applicationDirPath()+"/tmp/blockGet__"+name+QString("__%1__%2__%3__%4__%5__%6__%7.ano.eswc")
                .arg(x1).arg(x2).arg(y1).arg(y2).arg(z1).arg(z2).arg(cnt);
        x1*=cnt;x2*=cnt;y1*=cnt;y2*=cnt;z1*=cnt;z2*=cnt;
        V_NeuronSWC_list tosave;
        for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<testVNL.seg.size();i++)
        {
            NeuronTree SS;
            V_NeuronSWC seg_temp =  testVNL.seg.at(i);
            seg_temp.reverse();
            for(std::vector<V_NeuronSWC_unit>::size_type j=0;j<seg_temp.row.size();j++)
            {
                if(seg_temp.row.at(j).x>=x1&&seg_temp.row.at(j).x<=x2
                        &&seg_temp.row.at(j).y>=y1&&seg_temp.row.at(j).y<=y2
                        &&seg_temp.row.at(j).z>=z1&&seg_temp.row.at(j).z<=z2)
                {
                    tosave.seg.push_back(seg_temp);
                    break;
                }
            }
        }
        qDebug()<<"get nt size:"<<tosave.seg.size();
        auto nt=V_NeuronSWC_list__2__NeuronTree(tosave);
        writeESWC_file(BBSWCNAME,nt);
        return {BBSWCNAME.right(BBSWCNAME.size()-BBSWCNAME.lastIndexOf('/')),BBSWCNAME};
    }
}

QStringList getApoInBlock(const QString msg,const QList <CellAPO>& wholePoint)
{
    /*
     * p1:brain_id;res;x;y;z;size;socket.descriptor
     * p2:Neuron_id/name
     * 返回：文件名，文件路径
     */
    QStringList paraList=msg.split(";",Qt::SkipEmptyParts);
    QString brain_id=paraList.at(0).trimmed();//1. tf name/RES  2. .v3draw// test:17302;RES;x;y;z;b
    //0: 18465/RESx18000x13000x5150
    //1: 12520
    //2: 7000
    //3: 2916
    int res=paraList.at(1).toInt();//>0
    int xpos=paraList.at(2).toInt();
    int ypos=paraList.at(3).toInt();
    int zpos=paraList.at(4).toInt();
    int blocksize=paraList.at(5).toInt();

    {
        QString name=brain_id;
        int x1=xpos-blocksize;
        int x2=xpos+blocksize;
        int y1=ypos-blocksize;
        int y2=ypos+blocksize;
        int z1=zpos-blocksize;
        int z2=zpos+blocksize;
        int cnt=pow(2,res-1);


        dirCheck("tmp");
        QString BBAPONAME=QCoreApplication::applicationDirPath()+"/tmp/blockGet__"+name+QString("__%1__%2__%3__%4__%5__%6__%7.ano.apo")
                .arg(x1).arg(x2).arg(y1).arg(y2).arg(z1).arg(z2).arg(cnt);
        x1*=cnt;x2*=cnt;y1*=cnt;y2*=cnt;z1*=cnt;z2*=cnt;
        qDebug()<<"x1,x2,y1,y2,z1,z2"<<x1<<x2<<y1<<y2<<z1<<z2;
        QList <CellAPO> to_save;
        for(auto marker:wholePoint)
        {
            if(marker.x>=x1&&marker.x<=x2
              &&marker.y>=y1&&marker.y<=y2
              &&marker.z>=z1&&marker.z<=z2)
            {
                to_save.append(marker);
            }
        }
        qDebug()<<"to_save.size()="<<to_save.size();
        writeAPO_file(BBAPONAME,to_save);
        return {BBAPONAME.right(BBAPONAME.size()-BBAPONAME.lastIndexOf('/')),BBAPONAME};
    }
}

NeuronTree convertMsg2NT(QStringList pointlist,int client,int user, int isMany, int mode)
{
    NeuronTree newTempNT;
    newTempNT.listNeuron.clear();
    newTempNT.hashNeuron.clear();
    int cnt=pointlist.size();
    int timestamp=QDateTime::currentMSecsSinceEpoch();
    if(isMany==0){
        for(int i=0;i<cnt;i++)
        {
            NeuronSWC S;
            QStringList nodelist=pointlist[i].split(' ',Qt::SkipEmptyParts);
            if(nodelist.size()<4) return NeuronTree();
            S.n=i+1;
            S.type=nodelist[0].toUInt();

            S.x=nodelist[1].toFloat();
            S.y=nodelist[2].toFloat();
            S.z=nodelist[3].toFloat();
            switch (mode) {
            case 0:S.r=user*10+client;break;
            case 1:S.r=user;break;
            case 2:S.r=client;break;
            }

            if(i==0) S.pn=-1;
            else S.pn=i;
            S.timestamp=timestamp;
            newTempNT.listNeuron.push_back(S);
            newTempNT.hashNeuron.insert(S.n,newTempNT.listNeuron.size());
        }
    }
    else if(isMany==1){
        int index=0;
        for(int i=0;i<cnt;i++)
        {
            if(pointlist[i]!="$"){
                NeuronSWC S;
                QStringList nodelist=pointlist[i].split(' ',Qt::SkipEmptyParts);
                if(nodelist.size()<4) return NeuronTree();
                S.n=i+1;
                S.type=nodelist[0].toUInt();

                S.x=nodelist[1].toFloat();
                S.y=nodelist[2].toFloat();
                S.z=nodelist[3].toFloat();
                switch (mode) {
                case 0:S.r=user*10+client;break;
                case 1:S.r=user;break;
                case 2:S.r=client;break;
                }

                if(index==0) S.pn=-1;
                else S.pn=i;
                S.timestamp=timestamp;
                newTempNT.listNeuron.push_back(S);
                newTempNT.hashNeuron.insert(S.n,newTempNT.listNeuron.size());
                index++;
            }
            else{
                index=0;
            }

        }
    }

    return newTempNT;
}

vector<V_NeuronSWC>::iterator findseg(vector<V_NeuronSWC>::iterator begin,vector<V_NeuronSWC>::iterator end,const V_NeuronSWC seg)
{
    vector<V_NeuronSWC>::iterator result=end;
    double mindist=1.5;
    const std::vector<V_NeuronSWC_unit>::size_type cnt=seg.row.size();

    while(begin!=end)
    {
        if(begin->row.size()==cnt)
        {
            double dist=0;
            for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<cnt;i++)
            {
                auto node=begin->row.at(i);
                dist+=sqrt(
                           pow(node.x-seg.row[i].x,2)
                          +pow(node.y-seg.row[i].y,2)
                          +pow(node.z-seg.row[i].z,2)
                           );
            }
            if(dist/cnt<mindist)
            {
                mindist=dist/cnt;
                result=begin;
            }

            dist=0;
            for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<cnt;i++)
            {
                auto node=begin->row.at(i);
                dist+=sqrt(
                           pow(node.x-seg.row[cnt-i-1].x,2)
                          +pow(node.y-seg.row[cnt-i-1].y,2)
                          +pow(node.z-seg.row[cnt-i-1].z,2)
                           );
            }
            if(dist/cnt<mindist)
            {
                mindist=dist/cnt;
                result=begin;
            }
        }
        begin++;
    }
    return result;
}

double distance(const CellAPO &m1,const CellAPO &m2)
{
    return sqrt(
                (m1.x-m2.x)*(m1.x-m2.x)+
                (m1.y-m2.y)*(m1.y-m2.y)+
                (m1.z-m2.z)*(m1.z-m2.z)
            );
}

double distance(const double x1, const double x2, const double y1, const double y2, const double z1, const double z2){
    return sqrt(
                (x1-x2)*(x1-x2)+
                (y1-y2)*(y1-y2)+
                (z1-z2)*(z1-z2)
        );
}

double getSegLength(V_NeuronSWC &seg){
    int size=seg.row.size();
    double sum=0;
    for(int i=0;i<size-1;i++)
    {
        sum+=distance(seg.row[i].x,seg.row[i+1].x,seg.row[i].y,seg.row[i+1].y,seg.row[i].z,seg.row[i+1].z);
    }
    return sum;
}

double getPartOfSegLength(V_NeuronSWC &seg, int index1, int index2){
    double sum=0;
    for(int i=index1; i<index2; i++){
        sum+=distance(seg.row[i].x,seg.row[i+1].x,seg.row[i].y,seg.row[i+1].y,seg.row[i].z,seg.row[i+1].z);
    }
    return sum;
}

int findnearest(const CellAPO &m,const QList<CellAPO> &markers)
{
    int index=-1;
    double thres=1;
    for(int i=0;i<markers.size();i++){
        if(distance(m,markers[i])<thres){
            index=i;
        }
    }
    return index;
}

void init()
{
    for(int i=0;i<10;i++){
        for(int usercnt=100;usercnt<=100;usercnt+=10){
            QDir dir("/");
            for(int msgcnt=100;msgcnt<=100;msgcnt+=10){
                auto anoname=QString("/Users/huanglei/Desktop/dataserver/data/18454/18454_00049/pressure_18454_00049_%1_%2_%3/pressure_18454_00049_%1_%2_%3.ano")
                        .arg(usercnt).arg(msgcnt).arg(4000+i);
                QFile anof(anoname),apof(anoname+".apo"),swcf(anoname+".eswc");
                QString data=QString("APOFILE=pressure_18454_00049_%1_%2_%3.ano.apo\nSWCFILE=pressure_18454_00049_%1_%2_%3.ano.eswc").arg(usercnt).arg(msgcnt).arg(4000+i);
                dir.mkdir(QString("/Users/huanglei/Desktop/dataserver/data/18454/18454_00049/pressure_18454_00049_%1_%2_%3").arg(usercnt).arg(msgcnt).arg(4000+i));
                if(anof.open(QIODevice::WriteOnly)){
                    anof.write(data.toStdString().c_str(),data.size());
                    anof.close();
                }
                if(apof.open(QIODevice::WriteOnly)){
                    apof.close();
                }
                if(swcf.open(QIODevice::WriteOnly)){
                    swcf.close();
                }
            }
        }
    }

}

void stringToXYZ(string xyz, float& x, float& y, float& z){
    QString xyzstr = QString::fromStdString(xyz);
    QStringList xyzstrs = xyzstr.split("_");
    x = xyzstrs[0].toFloat();
    y = xyzstrs[1].toFloat();
    z = xyzstrs[2].toFloat();
}

void getSegmentsForOthersDetect(V_NeuronSWC_list& last1MinSegments, V_NeuronSWC_list& segmentsForOthersDetect, V_NeuronSWC_list segments){
    map<string, set<size_t>> wholeGrid2SegIDMap = getWholeGrid2SegIDMap(segments);

    set<size_t> tobeInsertSegIds;
    for(size_t i=0; i<last1MinSegments.seg.size(); i++){
        V_NeuronSWC seg = last1MinSegments.seg[i];

        for(size_t j=0; j<seg.row.size(); j++){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(wholeGrid2SegIDMap[gridKey].size()>=2){
                for(auto it=wholeGrid2SegIDMap[gridKey].begin(); it!=wholeGrid2SegIDMap[gridKey].end(); it++){
                    tobeInsertSegIds.insert(*it);
                }
            }
        }
    }
    set<size_t> tmpSegIds;

    for(auto it=tobeInsertSegIds.begin(); it!=tobeInsertSegIds.end(); it++){
        V_NeuronSWC seg = segments.seg[*it];

        for(size_t j=0; j<seg.row.size(); j++){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(wholeGrid2SegIDMap[gridKey].size()>=2){
                for(auto it2=wholeGrid2SegIDMap[gridKey].begin(); it2!=wholeGrid2SegIDMap[gridKey].end(); it2++){
                    tmpSegIds.insert(*it2);
                }
            }
        }
    }

    for(auto it=tmpSegIds.begin(); it!=tmpSegIds.end(); it++){
        tobeInsertSegIds.insert(*it);
    }

    for(auto it=tobeInsertSegIds.begin(); it!=tobeInsertSegIds.end(); it++){
        segmentsForOthersDetect.append(segments.seg[*it]);
    }
}

void getSegmentsForMissingDetect(V_NeuronSWC_list& last3MinSegments, V_NeuronSWC_list& segmentsForMissingDetect, V_NeuronSWC_list segments){
    map<string, set<size_t>> wholeGrid2SegIDMap = getWholeGrid2SegIDMap(segments);

    set<size_t> tobeInsertSegIds;
    for(size_t i=0; i<last3MinSegments.seg.size(); i++){
        V_NeuronSWC seg = last3MinSegments.seg[i];

        for(size_t j=0; j<seg.row.size(); j++){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(wholeGrid2SegIDMap[gridKey].size()>=2){
                for(auto it=wholeGrid2SegIDMap[gridKey].begin(); it!=wholeGrid2SegIDMap[gridKey].end(); it++){
                    tobeInsertSegIds.insert(*it);
                }
            }
        }
    }
    set<size_t> tmpSegIds;

    for(auto it=tobeInsertSegIds.begin(); it!=tobeInsertSegIds.end(); it++){
        V_NeuronSWC seg = segments.seg[*it];

        for(size_t j=0; j<seg.row.size(); j++){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(wholeGrid2SegIDMap[gridKey].size()>=2){
                for(auto it2=wholeGrid2SegIDMap[gridKey].begin(); it2!=wholeGrid2SegIDMap[gridKey].end(); it2++){
                    tmpSegIds.insert(*it2);
                }
            }
        }
    }

    for(auto it=tmpSegIds.begin(); it!=tmpSegIds.end(); it++){
        tobeInsertSegIds.insert(*it);
    }

    for(auto it=tobeInsertSegIds.begin(); it!=tobeInsertSegIds.end(); it++){
        segmentsForMissingDetect.append(segments.seg[*it]);
    }

}

void reverseSeg(V_NeuronSWC& seg){
    reverse(seg.row.begin(), seg.row.end());
    for(int i=0; i<seg.row.size(); i++){
        seg.row[i].n=i+1;
        seg.row[i].parent=i+2;
    }
    seg.row[seg.row.size()-1].parent=-1;
}

int getPointInSegIndex(string point, V_NeuronSWC& seg){
    int index=-1;
    for(int i=0; i<seg.row.size(); i++){
        float xLabel = seg.row[i].x;
        float yLabel = seg.row[i].y;
        float zLabel = seg.row[i].z;
        QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
        string gridKey = gridKeyQ.toStdString();
        if(point==gridKey){
            index=i;
            break;
        }
    }
    return index;
}

map<string, set<size_t>> getWholeGrid2SegIDMap(V_NeuronSWC_list inputSegments){
    map<string, set<size_t>> wholeGrid2SegIDMap;

    for(size_t i=0; i<inputSegments.seg.size(); ++i){
        V_NeuronSWC seg = inputSegments.seg[i];

        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            wholeGrid2SegIDMap[gridKey].insert(size_t(i));
        }
    }

    return wholeGrid2SegIDMap;
}

int isOverlapOfTwoSegs(QTextStream& logOut, V_NeuronSWC& seg1, V_NeuronSWC& seg2){
    if(seg1.row.size() == 1 || seg2.row.size() == 1){
        return 0;
    }

    double length1 = getSegLength(seg1);
    double length2 = getSegLength(seg2);
    double minDensity = min(length1/seg1.row.size(), length2/seg2.row.size());
    double minLength = min(length1, length2);
    //    double mindist = 3.5 - 3.5 * 1.0 / minLength;
    double mindist = 0.1 + 5 * sqrt(minLength)/(sqrt(minLength) + 11);
    //    double mindist_thres = 3.5 - 3.5 * 1.0 / minLength;
    double mindist_thres = 0.1 + 5 * sqrt(minLength)/(sqrt(minLength) + 11);

    if(minDensity < 5){
        mindist = 0.03 + 1.8 * sqrt(minLength)/(sqrt(minLength) + 11);
        mindist_thres = 0.03 + 1.8 * sqrt(minLength)/(sqrt(minLength) + 11);
    }

    if(seg1.row.size() == seg2.row.size()){
//        qDebug()<<"seg1.row.size() == seg2.row.size()";
        double dist=0;
        const std::vector<V_NeuronSWC_unit>::size_type cnt=seg1.row.size();
        for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<cnt;i++)
        {
            auto node=seg1.row.at(i);
            dist+=sqrt(
                pow(node.x-seg2.row[i].x,2)
                +pow(node.y-seg2.row[i].y,2)
                +pow(node.z-seg2.row[i].z,2)
            );
        }

        if(dist/cnt<mindist){
            qDebug()<<"1: "<<dist/cnt;
            seg1.printInfo(logOut);
            seg2.printInfo(logOut);
            return 1;
        }
        dist=0;
        for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<cnt;i++)
        {
            auto node=seg1.row.at(i);
            dist+=sqrt(
                pow(node.x-seg2.row[cnt-i-1].x,2)
                +pow(node.y-seg2.row[cnt-i-1].y,2)
                +pow(node.z-seg2.row[cnt-i-1].z,2)
            );
        }

        if(dist/cnt<mindist){
            logOut <<"2: "<<dist/cnt;
            seg1.printInfo(logOut);
            seg2.printInfo(logOut);
            return 1;
        }
        return 0;
    }

    bool isReverse = false;
    V_NeuronSWC seg_short = seg1;
    V_NeuronSWC seg_long = seg2;

    double length_short = length1;
    double length_long = length2;

    if(seg1.row.size() > seg2.row.size()){
        seg_short = seg2;
        seg_long = seg1;
        length_short = length2;
        length_long = length1;
        isReverse = true;
    }
//    qDebug()<<"seg_short"<<seg_short.row.size();
//    qDebug()<<"seg_long"<<seg_long.row.size();

    int long_index1 = -1;
    int long_index2 = -1;
    int start_index = -1;
    double mindist1 = 100;
    double mindist2 = 100;
    double mindist_final = 100;
    int seg_short_index = 0;
    for(auto it=seg_long.row.begin(); it!=seg_long.row.end(); it++){
        double dist = distance(it->x, seg_short.row[seg_short.row.size()-1].x, it->y, seg_short.row[seg_short.row.size()-1].y, it->z, seg_short.row[seg_short.row.size()-1].z);
        if(dist<mindist1){
            mindist1 = dist;
            long_index1 = it - seg_long.row.begin();
        }
    }

    for(auto it=seg_long.row.begin(); it!=seg_long.row.end(); it++){
        double dist = distance(it->x, seg_short.row[0].x, it->y, seg_short.row[0].y, it->z, seg_short.row[0].z);
        if(dist<mindist2){
            mindist2 = dist;
            long_index2 = it - seg_long.row.begin();
        }
    }

    if(long_index1 == -1 && long_index2 == -1)
        return 0;

    if(mindist1 < mindist2){
        start_index = long_index1;
        mindist_final = mindist1;
        seg_short_index = seg_short.row.size()-1;
    }else{
        start_index = long_index2;
        mindist_final = mindist2;
        seg_short_index = 0;
    }

    if(start_index == -1 || mindist_final >= mindist_thres){
        return 0;
    }else if(seg_long.row.size()-start_index >= seg_short.row.size()){
        double dist=0;
        const std::vector<V_NeuronSWC_unit>::size_type cnt=seg_short.row.size();
        for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<cnt;i++)
        {
            V_NeuronSWC_unit node;
            if(seg_short_index == seg_short.row.size()-1)
                node=seg_short.row.at(cnt-i-1);
            if(seg_short_index == 0)
                node=seg_short.row.at(i);

            dist+=sqrt(
                pow(node.x-seg_long.row[start_index+i].x,2)
                +pow(node.y-seg_long.row[start_index+i].y,2)
                +pow(node.z-seg_long.row[start_index+i].z,2)
            );
        }

        if(dist/cnt<mindist_thres){
            logOut <<"1:  "<<dist/cnt << "\n";
//            seg1.printInfo(logOut);
//            seg2.printInfo(logOut);
            if(!isReverse){
                if(length_short <= length_long)
                    return 1;
                else
                    return 2;
            }
            else{
                if(length_short <= length_long)
                    return 2;
                else
                    return 1;
            }
        }

    }else if(start_index+1 >= seg_short.row.size()){
        double dist=0;
        const std::vector<V_NeuronSWC_unit>::size_type cnt=seg_short.row.size();
        for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<cnt;i++)
        {
            V_NeuronSWC_unit node;
            if(seg_short_index == seg_short.row.size()-1)
                node=seg_short.row.at(cnt-i-1);
            if(seg_short_index == 0)
                node=seg_short.row.at(i);

            dist+=sqrt(
                pow(node.x-seg_long.row[start_index-i].x,2)
                +pow(node.y-seg_long.row[start_index-i].y,2)
                +pow(node.z-seg_long.row[start_index-i].z,2)
            );
        }

        if(dist/cnt<mindist_thres){
//            seg1.printInfo(logOut);
//            seg2.printInfo(logOut);
            logOut <<"2:  "<<dist/cnt << "\n";
            if(!isReverse){
                if(length_short <= length_long)
                    return 1;
                else
                    return 2;
            }
            else{
                if(length_short <= length_long)
                    return 2;
                else
                    return 1;
            }
        }

    }else{
        return 0;
    }

    return 0;
}

QStringList V_NeuronSWCToSendMSG(V_NeuronSWC seg)
{
    QStringList result;
    for(int i=0;i<seg.row.size();i++)   //why  i need  < 120, does msg has length limitation? liqi 2019/10/7
    {
        V_NeuronSWC_unit curSWCunit = seg.row[i];
        result.push_back(QString("%1 %2 %3 %4").arg(curSWCunit.type).arg(curSWCunit.x).arg(curSWCunit.y).arg(curSWCunit.z));
        //        if(i==seg.row.size()-1)
        //            AutoTraceNode=XYZ(GlobalCroods.x,GlobalCroods.y,GlobalCroods.z);
    }
    return result;
}

RGB8 getColorFromType(int type){
    RGB8 color;
    color.r=0;
    color.g=0;
    color.b=0;
    if(type<0 || type>20)
        return color;
    else{
        color.r=neuron_type_color[type][0];
        color.g=neuron_type_color[type][1];
        color.b=neuron_type_color[type][2];
        return color;
    }
}

set<int> getQCMarkerNearBy(vector<V_NeuronSWC> &segs, const QList<CellAPO> &markers){
    set<int> indexs;
    for(auto segIt=segs.begin(); segIt!=segs.end(); segIt++){
        auto seg = *segIt;
        bool getFirst = false;
        bool getSecond = false;
        float x1 = seg.row[0].x;
        float y1 = seg.row[0].y;
        float z1 = seg.row[0].z;
        float x2 = seg.row[seg.nrows()-1].x;
        float y2 = seg.row[seg.nrows()-1].y;
        float z2 = seg.row[seg.nrows()-1].z;
        for(auto it=markers.begin(); it!=markers.end(); it++){
            if(getFirst && getSecond)
                break;

            if(distance(x1, it->x, y1, it->y, z1, it->z) < 1 && it->comment=="quality_control"){
                getFirst = true;
                indexs.insert(it-markers.begin());
            }
            if(distance(x2, it->x, y2, it->y, z2, it->z) < 1 && it->comment=="quality_control"){
                getSecond = true;
                indexs.insert(it-markers.begin());
            }
        }
    }

    return indexs;
}

void mergeresultFiles(const string& folderPath){
    ostringstream outputContent;

    //遍历文件夹
    for(const auto& entry : filesystem::directory_iterator(folderPath)){
        if(entry.is_regular_file() && entry.path().extension() == ".txt"){
            ifstream inputFile(entry.path());
            if(!inputFile){
                cerr << "Failed to open input file: " << entry.path();
                continue;
            }

            //将当前文本文件的内容逐行写入到内存缓冲区
            ostringstream fileContent;
            fileContent << inputFile.rdbuf();
            outputContent << entry.path().filename();
            outputContent << fileContent.str();
            outputContent << "\n\n";

            inputFile.close();
        }
    }

    // 将内存缓冲区中的内容写入到目标文件中
    ofstream outputFile(folderPath + "/summary_results.txt");
    if(!outputFile){
        cerr << "Failed to open output file: " << folderPath + "/summary_results.txt" << std::endl;
        return;
    }
    outputFile << outputContent.str();
    outputFile.close();

    cout << "Text files merged successfully!" << endl;
}

void mergeResultCSVFiles(const string& folderPath, const string& outputPath){
    ofstream outputFile(outputPath);
    if(!outputFile){
        cerr << "Failed to open output File: " << outputPath << endl;
        return;
    }
    bool isFirst = true;
    for(const auto& entry : filesystem::directory_iterator(folderPath)){
        if(entry.is_regular_file() && entry.path().extension() == ".csv"){
            ifstream inputFile(entry.path());
            if(!inputFile){
                 cerr << "Failed to open input file: " << entry.path();
                continue;
            }

            string header;
            if(getline(inputFile, header)){
                if(isFirst){
                    outputFile << header << endl;
                    isFirst = false;
                }
            }

            string line;
            while(getline(inputFile, line)){
                outputFile << line << endl;
            }

            inputFile.close();
        }
    }

    outputFile.close();
    cout << "Csv files merged successfully!" << endl;
}

void mergeResultCSVFilesAll(const QStringList& allPaths, const string& outputPath){
    ofstream outputFile(outputPath);
    if(!outputFile){
        cerr << "Failed to open output File: " << outputPath << endl;
        return;
    }
    bool isFirst = true;
    for(auto path:allPaths){
        ifstream inputFile(path.toStdString());
        if(!inputFile){
            cerr << "Failed to open input file: " << path.toStdString();
            continue;
        }

        string header;
        if(getline(inputFile, header)){
            if(isFirst){
                outputFile << header << endl;
                isFirst = false;
            }
        }

        string line;
        while(getline(inputFile, line)){
            outputFile << line << endl;
        }

        inputFile.close();
    }

    outputFile.close();
    cout << "Csv files merged successfully!" << endl;
}

void mergeErrorResultCSVFilesAll(const QStringList& allPaths, const string& outputPath){
    ofstream outputFile(outputPath);
    if(!outputFile){
        cerr << "Failed to open output File: " << outputPath << endl;
        return;
    }
    bool isFirst = true;
    for(auto path:allPaths){
        ifstream inputFile(path.toStdString());
        if(!inputFile){
            cerr << "Failed to open input file: " << path.toStdString();
            continue;
        }

        string header;
        if(getline(inputFile, header)){
            if(isFirst){
                outputFile << header << endl;
                isFirst = false;
            }
        }

        string line;
        while(getline(inputFile, line)){
            vector<string> tokens = splitCSVLine(line);
            string somaDetectedCondition = tokens[1];
            int allNum = stoi(tokens[14]);
            if(somaDetectedCondition == "已检测出胞体位置"){
                if(allNum != 0){
                    outputFile << line << endl;
                }
            }
            else if(somaDetectedCondition == "自动设置多分叉为胞体位置"){
                if(allNum != 0){
                    outputFile << line << endl;
                }
            }
            else{
                outputFile << line << endl;
            }
        }

        inputFile.close();
    }

    outputFile.close();
    cout << "Csv files merged successfully!" << endl;
}

std::vector<std::string> splitCSVLine(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool insideQuotes = false;

    for (char c : line) {
        if (c == '\"') {
            insideQuotes = !insideQuotes; // Toggle the insideQuotes flag
        } else if (c == ',' && !insideQuotes) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(field); // Add the last field

    // Optionally, you can trim quotes from fields
    for (auto& f : fields) {
        if (!f.empty() && f.front() == '\"' && f.back() == '\"') {
            f = f.substr(1, f.size() - 2);
        }
    }

    return fields;
}

QStringList getAllTargetPaths(QString rootPath, QStringList filters){
    QStringList swcFiles;
    QDirIterator it(rootPath, filters, QDir::Files, QDirIterator::Subdirectories);
    while(it.hasNext()){
        swcFiles << it.next();
    }
    return swcFiles;
}
