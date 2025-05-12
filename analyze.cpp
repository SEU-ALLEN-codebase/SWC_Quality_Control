﻿#include "analyze.h"
#include "sort_swc.h"

vector<int> getMulfurcationsCountNearSoma(float dist_thre, XYZ somaCoordinate, V_NeuronSWC_list segments){
    vector<int> counts;

    map<string, set<size_t> > wholeGrid2SegIDMap;
    map<string, bool> isEndPointMap;

    set<string> allPoints;

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];

        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            wholeGrid2SegIDMap[gridKey].insert(size_t(i));
            allPoints.insert(gridKey);

            if(j == 0 || j == seg.row.size() - 1){
                isEndPointMap[gridKey] = true;
            }
        }
    }

    //末端点和分叉点
    vector<string> points;
    vector<set<int>> linksIndex;

    map<string,int> pointsIndexMap;

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];
        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(j==0 || j==seg.row.size()-1){
                //在pointsIndexMap中找不到某个线的末端点
                if(pointsIndexMap.find(gridKey) == pointsIndexMap.end()){
                    points.push_back(gridKey);
                    linksIndex.push_back(set<int>());
                    //                    linksIndexVec.push_back(vector<int>());
                    pointsIndexMap[gridKey] = points.size() - 1;
                }
            }else{
                if(wholeGrid2SegIDMap[gridKey].size()>1 &&
                    isEndPointMap.find(gridKey) != isEndPointMap.end() &&
                    pointsIndexMap.find(gridKey) == pointsIndexMap.end()){
                    points.push_back(gridKey);
                    linksIndex.push_back(set<int>());
                    //                    linksIndexVec.push_back(vector<int>());
                    pointsIndexMap[gridKey] = points.size() - 1;
                }
            }
        }
    }
    qDebug()<<"points size: "<<points.size();

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];
        vector<int> segIndexs;
        set<int> segIndexsSet;
        segIndexs.clear();
        segIndexsSet.clear();
        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(pointsIndexMap.find(gridKey) != pointsIndexMap.end()){
                int index = pointsIndexMap[gridKey];
                if(segIndexsSet.find(index) == segIndexsSet.end()){
                    segIndexs.push_back(index);
                    segIndexsSet.insert(index);
                }
            }
        }
        //        qDebug()<<"i : "<<i<<"seg size: "<<seg.row.size()<<" segIndexsSize: "<<segIndexs.size();
        for(size_t j=0; j<segIndexs.size()-1; ++j){
            if(segIndexs[j] == 1 || segIndexs[j+1] == 1){
                qDebug()<<segIndexs[j]<<" "<<segIndexs[j+1];
            }
            linksIndex[segIndexs[j]].insert(segIndexs[j+1]);
            //            linksIndexVec[segIndexs[j]].push_back(segIndexs[j+1]);
            linksIndex[segIndexs[j+1]].insert(segIndexs[j]);
            //            linksIndexVec[segIndexs[j+1]].push_back(segIndexs[j]);
        }
    }

    qDebug()<<"link map end";

    int biCount=0;
    int mulCount=0;
    int othersCount=0;

    for(size_t i=0; i<points.size(); ++i){
        if(linksIndex[i].size() > 3){
            NeuronSWC s;
            stringToXYZ(points[i],s.x,s.y,s.z);
            if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y, s.z, somaCoordinate.z)<dist_thre)
                mulCount++;
        }else if(linksIndex[i].size() == 3){
            NeuronSWC s;
            stringToXYZ(points[i],s.x,s.y,s.z);
            if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y, s.z, somaCoordinate.z)<dist_thre)
                biCount++;
        }else if(linksIndex[i].size()==2){
            NeuronSWC s;
            stringToXYZ(points[i],s.x,s.y,s.z);
            if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y, s.z, somaCoordinate.z)<dist_thre)
                othersCount++;
        }
    }

    counts.push_back(othersCount);
    counts.push_back(biCount);
    counts.push_back(mulCount);
    qDebug()<<"counts: "<<counts[0]<<" "<<counts[1]<<" "<<counts[2];
    return counts;

}

map<string, set<int>> getColorChangedPoints(V_NeuronSWC_list segments){
    map<string, set<int>> point2TypeMap;

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];

        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            point2TypeMap[gridKey].insert(seg.row[j].type);
        }
    }

    map<string, set<int>> specPointsMap;
    for(auto it=point2TypeMap.begin(); it!=point2TypeMap.end(); it++){
        if(it->second.size()<=1)
            continue;
        specPointsMap[it->first]=it->second;
    }

    return specPointsMap;
}

set<string> getDissociativeSegEndPoints(V_NeuronSWC_list segments){
    map<string, set<size_t> > wholeGrid2SegIDMap;

    set<string> dissociativePoints;

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];

        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            wholeGrid2SegIDMap[gridKey].insert(size_t(i));

        }
    }

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];
        bool flag=true;
        string savedgridKey;

        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            int size=wholeGrid2SegIDMap[gridKey].size();

            if(j == 0)
                savedgridKey=gridKey;
            if(size > 1){
                flag = false;
            }
        }

        if(flag){
            dissociativePoints.insert(savedgridKey);
        }
    }

    return dissociativePoints;
}

set<string> getAngleErrPoints(QTextStream& logOut, float dist_thre, XYZ somaCoordinate, V_NeuronSWC_list& segments, bool isSomaExists, bool needConsiderType){
    set<string> angleErrPoints;

    if(!isSomaExists){
        return angleErrPoints;
    }
    map<string, set<int>> point2TypeMap;
    map<string, set<string>> parentMap;
    map<string, set<string>> childMap;
    map<string, set<size_t> > wholeGrid2SegIDMap;
    map<string, bool> isEndPointMap;
    set<string> allPoints;

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];

        vector<int> rowN2Index(seg.row.size()+1);

        for(size_t j=0; j<seg.row.size(); ++j){
            rowN2Index[seg.row[j].n]=j;
        }

        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            point2TypeMap[gridKey].insert(seg.row[j].type);
            wholeGrid2SegIDMap[gridKey].insert(size_t(i));
            allPoints.insert(gridKey);

            if(j == 0 || j == seg.row.size() - 1){
                isEndPointMap[gridKey] = true;
            }

            if(seg.row[j].parent!=-1){
                float x2Label=seg.row[rowN2Index[seg.row[j].parent]].x;
                float y2Label=seg.row[rowN2Index[seg.row[j].parent]].y;
                float z2Label=seg.row[rowN2Index[seg.row[j].parent]].z;
                QString parentKeyQ=QString::number(x2Label) + "_" + QString::number(y2Label) + "_" + QString::number(z2Label);
                string parentKey=parentKeyQ.toStdString();
                parentMap[gridKey].insert(parentKey);
                childMap[parentKey].insert(gridKey);
            }
        }
    }

    //末端点和分叉点
    vector<string> points;
    vector<set<int>> linksIndex;

    map<string,int> pointsIndexMap;

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];
        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(j==0 || j==seg.row.size()-1){
                //在pointsIndexMap中找不到某个线的末端点
                if(pointsIndexMap.find(gridKey) == pointsIndexMap.end()){
                    points.push_back(gridKey);
                    linksIndex.push_back(set<int>());
                    //                    linksIndexVec.push_back(vector<int>());
                    pointsIndexMap[gridKey] = points.size() - 1;
                }
            }else{
                if(wholeGrid2SegIDMap[gridKey].size()>1 &&
                    isEndPointMap.find(gridKey) != isEndPointMap.end() &&
                    pointsIndexMap.find(gridKey) == pointsIndexMap.end()){
                    points.push_back(gridKey);
                    linksIndex.push_back(set<int>());
                    //                    linksIndexVec.push_back(vector<int>());
                    pointsIndexMap[gridKey] = points.size() - 1;
                }
            }
        }
    }
    qDebug()<<"points size: "<<points.size();

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];
        vector<int> segIndexs;
        set<int> segIndexsSet;
        segIndexs.clear();
        segIndexsSet.clear();
        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(pointsIndexMap.find(gridKey) != pointsIndexMap.end()){
                int index = pointsIndexMap[gridKey];
                if(segIndexsSet.find(index) == segIndexsSet.end()){
                    segIndexs.push_back(index);
                    segIndexsSet.insert(index);
                }
            }
        }
        //        qDebug()<<"i : "<<i<<"seg size: "<<seg.row.size()<<" segIndexsSize: "<<segIndexs.size();
        for(size_t j=0; j<segIndexs.size()-1; ++j){
            if(segIndexs[j] == 1 || segIndexs[j+1] == 1){
                qDebug()<<segIndexs[j]<<" "<<segIndexs[j+1];
            }
            linksIndex[segIndexs[j]].insert(segIndexs[j+1]);
            //            linksIndexVec[segIndexs[j]].push_back(segIndexs[j+1]);
            linksIndex[segIndexs[j+1]].insert(segIndexs[j]);
            //            linksIndexVec[segIndexs[j+1]].push_back(segIndexs[j]);
        }
    }

    qDebug()<<"link map end";

    set<string> bifurcationPoints;
    for(size_t i=0; i<points.size(); ++i){
        if(linksIndex[i].size() == 3){
            NeuronSWC s;
            stringToXYZ(points[i],s.x,s.y,s.z);
            if(!isSomaExists)
                bifurcationPoints.insert(points[i]);
            else if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y, s.z, somaCoordinate.z)>dist_thre)
                bifurcationPoints.insert(points[i]);
        }
    }

    map<string, vector<int>> twoSegsMap;
    map<string, vector<int>> threeSegsMap;
    for(auto it=bifurcationPoints.begin(); it!=bifurcationPoints.end();){
        set<int> types = point2TypeMap[*it];
//        qDebug()<<"types: "<<point2TypeMap[*it].size();
        bool flag = true;
        if(needConsiderType){
            if (types.size()!=1) {
                it = bifurcationPoints.erase(it); // 通过迭代器去除元素，并返回下一个有效迭代器
                flag = false;
            }
            else if(*types.begin()!=3 && *types.begin()!=4){
                it = bifurcationPoints.erase(it);
                flag = false;
            }
        }
        if(flag){
            set<size_t> segIds = wholeGrid2SegIDMap[*it];
            bool isVaild = true;
            vector<int> segIndexs;
            if(segIds.size() == 3){
                map<int, int> segId2CoorIndex;
                int parentSegId;
                double minDist = 100000;
                //先找到离soma最近的点
                if(isSomaExists){
                    for(auto segIt=segIds.begin(); segIt!=segIds.end(); segIt++){
                        V_NeuronSWC seg = segments.seg[*segIt];
                        int index = getPointInSegIndex(*it, seg);
                        segId2CoorIndex[*segIt] = index;
                        if(index==0){
                            XYZ coor = seg.row[1];
                            if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z) < minDist){
                                minDist = distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z);
                                parentSegId = *segIt;
                            }
                        }
                        else if(index==seg.row.size()-1){
                            XYZ coor = seg.row[seg.row.size()-2];
                            if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z) < minDist){
                                minDist = distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z);
                                parentSegId = *segIt;
                            }
                        }
                    }

                    //调整线段走向
                    for(auto mapIt = segId2CoorIndex.begin(); mapIt != segId2CoorIndex.end(); mapIt++){
                        if(parentSegId == mapIt->first && mapIt->second != 0){
                            reverseSeg(segments.seg[mapIt->first]);
                        }
                        if(parentSegId != mapIt->first && mapIt->second != segments.seg[mapIt->first].row.size() - 1){
                            reverseSeg(segments.seg[mapIt->first]);
                        }
                    }
                }

                for(auto segIt=segIds.begin(); segIt!=segIds.end(); segIt++){
                    double length = getSegLength(segments.seg[*segIt]);
                    if(length < 15){
    //                            it = bifurcationPoints.erase(it);
                        isVaild=false;
                        break;
                    }
                    int index = getPointInSegIndex(*it, segments.seg[*segIt]);
                    if(index == 0){
                        segIndexs.insert(segIndexs.begin(), *segIt);
                    }
                    else{
                        segIndexs.push_back(*segIt);
                    }
                }
            }
            else if(segIds.size() == 2){
                //先放parent, 再放chi
                int countNoParent=0;
                map<int, int> segId2CoorIndex;
                int parentSegId;
                XYZ minDistCoor;
                double minDist = 100000;
                //先找到离soma最近的点
                for(auto segIt=segIds.begin(); segIt!=segIds.end(); segIt++){
                    V_NeuronSWC seg = segments.seg[*segIt];
                    int index = getPointInSegIndex(*it, seg);
                    segId2CoorIndex[*segIt] = index;
                    if(index==0){
                        XYZ coor = seg.row[1];
                        if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z) < minDist){
                            minDist = distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z);
                            parentSegId = *segIt;
                            minDistCoor = coor;
                        }
                    }
                    else if(index==seg.row.size()-1){
                        XYZ coor = seg.row[seg.row.size()-2];
                        if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z) < minDist){
                            minDist = distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y, coor.z, somaCoordinate.z);
                            parentSegId = *segIt;
                            minDistCoor = coor;
                        }
                    }
                    else{
                        XYZ coor1 = seg.row[index + 1];
                        XYZ coor2 = seg.row[index - 1];
                        if(distance(coor1.x, somaCoordinate.x, coor1.y, somaCoordinate.y, coor1.z, somaCoordinate.z) < minDist){
                            minDist = distance(coor1.x, somaCoordinate.x, coor1.y, somaCoordinate.y, coor1.z, somaCoordinate.z);
                            parentSegId = *segIt;
                            minDistCoor = coor1;
                        }
                        if(distance(coor2.x, somaCoordinate.x, coor2.y, somaCoordinate.y, coor2.z, somaCoordinate.z) < minDist){
                            minDist = distance(coor2.x, somaCoordinate.x, coor2.y, somaCoordinate.y, coor2.z, somaCoordinate.z);
                            parentSegId = *segIt;
                            minDistCoor = coor2;
                        }
                    }
                }

                //调整线段走向
                for(auto mapIt = segId2CoorIndex.begin(); mapIt != segId2CoorIndex.end(); mapIt++){
                    if(parentSegId == mapIt->first && minDistCoor != segments.seg[mapIt->first].row[mapIt->second + 1]){
                        reverseSeg(segments.seg[mapIt->first]);
                    }
                    if(parentSegId != mapIt->first && mapIt->second != segments.seg[mapIt->first].row.size() - 1){
                        reverseSeg(segments.seg[mapIt->first]);
                    }
                }

                for(auto segIt=segIds.begin(); segIt!=segIds.end(); segIt++){
                    double length = getSegLength(segments.seg[*segIt]);
                    if(length < 15){
                        //                            it = bifurcationPoints.erase(it);
                        isVaild=false;
                        break;
                    }
                    int index = getPointInSegIndex(*it, segments.seg[*segIt]);
                    if(index != segments.seg[*segIt].row.size() - 1){
                        segIndexs.insert(segIndexs.begin(), *segIt);
                    }
                    else{
                        segIndexs.push_back(*segIt);
                    }
                }
            }

            if(isVaild){
                if(segIndexs.size() == 2)
                    twoSegsMap[*it] = segIndexs;
                else if(segIndexs.size() == 3){
                    threeSegsMap[*it] = segIndexs;
                }
                it++;
            }
            else{
                it = bifurcationPoints.erase(it);
            }
        }
    }

    for(auto it = threeSegsMap.begin(); it != threeSegsMap.end(); it++){
        XYZ curCoor;
        stringToXYZ(it->first, curCoor.x, curCoor.y, curCoor.z);
        XYZ parCoor = segments.seg[it->second[0]].row[1];
        vector<XYZ> chiCoors;
        for(auto chiSegIt=it->second.begin() + 1; chiSegIt!=it->second.end(); chiSegIt++){
            XYZ chiCoor;
            chiCoor.x = segments.seg[*chiSegIt].row[segments.seg[*chiSegIt].row.size() - 2].x;
            chiCoor.y = segments.seg[*chiSegIt].row[segments.seg[*chiSegIt].row.size() - 2].y;
            chiCoor.z = segments.seg[*chiSegIt].row[segments.seg[*chiSegIt].row.size() - 2].z;
            chiCoors.push_back(chiCoor);
        }
        XYZ parCoor_real = parCoor, chiCoor1_real = chiCoors[0], chiCoor2_real = chiCoors[1];

        float angle1 = 0.0, angle2 = 0.0;

        QVector3D vector1(parCoor_real.x-curCoor.x, parCoor_real.y-curCoor.y, parCoor_real.z-curCoor.z);
        double parent_dist = getPartOfSegLength(segments.seg[it->second[0]], 0, 1);
        int parent_vector_count = 1;
        for(int j = 2; j < segments.seg[it->second[0]].row.size(); j++){
            QVector3D tmp_vector(segments.seg[it->second[0]].row[j].x - segments.seg[it->second[0]].row[j - 1].x,
                                 segments.seg[it->second[0]].row[j].y - segments.seg[it->second[0]].row[j - 1].y,
                                 segments.seg[it->second[0]].row[j].z - segments.seg[it->second[0]].row[j - 1].z);
            vector1 += tmp_vector;
            parent_vector_count++;
            parent_dist += getPartOfSegLength(segments.seg[it->second[0]], j - 1, j);
            if(parent_dist >= 26){
                break;
            }
        }

        V_NeuronSWC_unit pre_unit;
        double child1_dist = 0.0;
        int child1_vector_count = 0;
        for(auto coorIt = segments.seg[it->second[1]].row.rbegin(); coorIt != segments.seg[it->second[1]].row.rend(); coorIt++)
        {
            if(coorIt == segments.seg[it->second[1]].row.rbegin()){
                pre_unit = *coorIt;
                continue;
            }
            QVector3D vector2(coorIt->x-pre_unit.x, coorIt->y-pre_unit.y, coorIt->z-pre_unit.z);
            child1_dist += distance(coorIt->x, pre_unit.x, coorIt->y, pre_unit.y, coorIt->z, pre_unit.z);
            pre_unit = *coorIt;
            angle1 += calculateAngleofVecs(vector1, vector2);
            child1_vector_count++;
            if(child1_dist >= 26){
                break;
            }
        }
        angle1 /= child1_vector_count;

        double child2_dist = 0.0;
        int child2_vector_count = 0;
        for(auto coorIt = segments.seg[it->second[2]].row.rbegin(); coorIt != segments.seg[it->second[2]].row.rend(); coorIt++)
        {
            if(coorIt == segments.seg[it->second[2]].row.rbegin()){
                pre_unit = *coorIt;
                continue;
            }
            QVector3D vector3(coorIt->x-pre_unit.x, coorIt->y-pre_unit.y, coorIt->z-pre_unit.z);
            child2_dist += distance(coorIt->x, pre_unit.x, coorIt->y, pre_unit.y, coorIt->z, pre_unit.z);
            pre_unit = *coorIt;
            angle2 += calculateAngleofVecs(vector1, vector3);
            child2_vector_count++;
            if(child2_dist >= 26){
                break;
            }
        }
        angle2 /= child2_vector_count;

        if((angle1>0 && angle1<60) || (angle2>0 && angle2<60)){
            angleErrPoints.insert(it->first);
        }

    }

    for(auto it = twoSegsMap.begin(); it != twoSegsMap.end(); it++){
        XYZ curCoor;
        stringToXYZ(it->first, curCoor.x, curCoor.y, curCoor.z);
        int index = getPointInSegIndex(it->first, segments.seg[it->second[0]]);

        XYZ parCoor = segments.seg[it->second[0]].row[index + 1];
        vector<XYZ> chiCoors;
        chiCoors.push_back(segments.seg[it->second[0]].row[index - 1]);
        chiCoors.push_back(segments.seg[it->second[1]].row[segments.seg[it->second[1]].row.size() - 2]);

        XYZ parCoor_real = parCoor, chiCoor1_real = chiCoors[0], chiCoor2_real = chiCoors[1];

        float angle1 = 0.0, angle2 = 0.0;

        QVector3D vector1(parCoor_real.x-curCoor.x, parCoor_real.y-curCoor.y, parCoor_real.z-curCoor.z);
        double parent_dist = getPartOfSegLength(segments.seg[it->second[0]], index, index + 1);
        int parent_vector_count = 1;
        for(int j = index + 2; j < segments.seg[it->second[0]].row.size(); j++){
            QVector3D tmp_vector(segments.seg[it->second[0]].row[j].x - segments.seg[it->second[0]].row[j - 1].x,
                                 segments.seg[it->second[0]].row[j].y - segments.seg[it->second[0]].row[j - 1].y,
                                 segments.seg[it->second[0]].row[j].z - segments.seg[it->second[0]].row[j - 1].z);
            vector1 += tmp_vector;
            parent_vector_count++;
            parent_dist += getPartOfSegLength(segments.seg[it->second[0]], j - 1, j);
            if(parent_dist >= 26){
                break;
            }
        }

        V_NeuronSWC_unit pre_unit;
        double child1_dist = 0.0;
        int child1_vector_count = 0;
        for(auto coorIt = segments.seg[it->second[1]].row.rbegin(); coorIt != segments.seg[it->second[1]].row.rend(); coorIt++)
        {
            if(coorIt == segments.seg[it->second[1]].row.rbegin()){
                pre_unit = *coorIt;
                continue;
            }
            QVector3D vector2(coorIt->x-pre_unit.x, coorIt->y-pre_unit.y, coorIt->z-pre_unit.z);
            child1_dist += distance(coorIt->x, pre_unit.x, coorIt->y, pre_unit.y, coorIt->z, pre_unit.z);
            pre_unit = *coorIt;
            angle1 += calculateAngleofVecs(vector1, vector2);
            child1_vector_count++;
            if(child1_dist >= 26){
                break;
            }
        }
        angle1 /= child1_vector_count;

        pre_unit.x = curCoor.x;
        pre_unit.y = curCoor.y;
        pre_unit.z = curCoor.z;
        double child2_dist = 0.0;
        int child2_vector_count = 0;
        for(int i = index - 1; i >= 0; i--)
        {
            XYZ coor = segments.seg[it->second[0]].row[i];
            QVector3D vector3(coor.x-pre_unit.x, coor.y-pre_unit.y, coor.z-pre_unit.z);
            child2_dist += distance(coor.x, pre_unit.x, coor.y, pre_unit.y, coor.z, pre_unit.z);
            pre_unit.x = coor.x;
            pre_unit.y = coor.y;
            pre_unit.z = coor.z;

            angle2 += calculateAngleofVecs(vector1, vector3);
            child2_vector_count++;
            if(child2_dist >= 26){
                break;
            }
        }
        angle2 /= child2_vector_count;

        if((angle1>0 && angle1<60) || (angle2>0 && angle2<60)){
            angleErrPoints.insert(it->first);
        }
    }
    return angleErrPoints;
}

float calculateAngleofVecs(QVector3D vector1, QVector3D vector2){
    // 计算两个向量的夹角
    qreal dotProduct = QVector3D::dotProduct(vector1, vector2);
    qreal magnitude1 = vector1.length();
    qreal magnitude2 = vector2.length();

    qreal cosineAngle = dotProduct / (magnitude1 * magnitude2);
    // 计算弧度值
    qreal angleRadians = qAcos(cosineAngle);

    // 将弧度转换为角度
    qreal angleDegrees = qRadiansToDegrees(angleRadians);

    return angleDegrees;
}

//将soma点的半径设置为1.234
void setSomaPointRadius(QString fileSaveName, V_NeuronSWC_list segments, XYZ somaCoordinate){
    map<string, set<size_t>> wholeGrid2SegIDMap = getWholeGrid2SegIDMap(segments);
    QString somaGridQ = QString::number(somaCoordinate.x) + "_" + QString::number(somaCoordinate.y) + "_" + QString::number(somaCoordinate.z);
    string somaGrid = somaGridQ.toStdString();

    V_NeuronSWC seg = segments.seg[*wholeGrid2SegIDMap[somaGrid].begin()];
    for(int i=0; i<seg.row.size(); i++){
        float xLabel = seg.row[i].x;
        float yLabel = seg.row[i].y;
        float zLabel = seg.row[i].z;
        QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
        string gridKey = gridKeyQ.toStdString();
        if(gridKey == somaGrid){
            segments.seg[*wholeGrid2SegIDMap[somaGrid].begin()].row[i].r = 1.234;
            break;
        }
    }

    writeESWC_file(fileSaveName, V_NeuronSWC_list__2__NeuronTree(segments));
}

int getSomaNumberFromSwcFile(QString filePath, float r){
    int number = -1;
    if (filePath.endsWith(".swc") || filePath.endsWith(".SWC") || filePath.endsWith(".eswc") || filePath.endsWith(".ESWC"))
    {
        QFile qf(filePath);
        QString arryRead;
        if(!qf.open(QIODevice::ReadOnly|QIODevice::Text)){
            return -2;
        }
        arryRead=qf.readAll();

        qf.close();

        QStringList arryListWrite= arryRead.split("\n");
        //        for(int i=0;i<arryListWrite.size();i++){
        //            qDebug()<<arryListWrite.at(i);
        //        }

        // QIODevice::Text:以文本方式打开文件，读取时“\n”被自动翻译为换行符，写入时字符串结束符会自动翻译为系统平台的编码，如 Windows 平台下是“\r\n”
        if (!qf.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            return -2;
        }
        QTextStream streamWrite(&qf);
        for(int i=0;i<arryListWrite.size()-1;i++){      //这里到arryListWrite.size()-1是因为arryListWrite数组按照\n分 段时，最后一行尾部有个\n，所以数组最后一个值为空，需要将它去掉
            if(arryListWrite.at(i).contains("#")){
                streamWrite<<arryListWrite.at(i)<<"\n";
            }else{
                QString contentWrite= arryListWrite.at(i);
                QStringList swcInfo=contentWrite.split(' ',Qt::SkipEmptyParts);
                if(swcInfo[5] == QString::number(r))
                {
                    number = swcInfo[0].toInt();
                    swcInfo[5] = "1.000";
                }
                contentWrite=swcInfo.join(' ');
                streamWrite<<contentWrite<<"\n";
            }
        }
        qf.close();
        return number;
    }

}

void analyzeSomaNearBy(QTextStream& logOut, XYZ somaCoordinate, V_NeuronSWC_list segments, float dist_thre){
    logOut << "begin analyzeSomaNearBy...\n";
//    if(!myServer->isSomaExists){
//        qDebug()<<"soma not detected!";
//    }
    vector<int> counts=getMulfurcationsCountNearSoma(dist_thre, somaCoordinate, segments);
    if((counts[1] + counts[2])==1){
        logOut << "the soma has been connected to one point\n";
    }else if((counts[1] + counts[2])==0 && counts[0]==1){
        logOut << "the soma has been connected to one point\n";
    }else{
        logOut << "the soma is not connected to one point!";
    }
    logOut << "analyzeSomaNearBy end\n";
}

vector<CellAPO> analyzeColorMutation(QTextStream& logOut, bool isSomaExists, XYZ somaCoordinate, V_NeuronSWC_list segments, float dist_thre){
    logOut << "begin analyzeColorMutation...\n";
    vector<CellAPO> errorMarkers;
    if(!isSomaExists){
        return errorMarkers;
    }
    bool result=true;
//    if(!myServer->isSomaExists){
//        qDebug()<<"soma not detected!";
//        QString tobeSendMsg="/FEEDBACK_ANALYZE_ColorMutation:";
//        tobeSendMsg += QString("server %1 %2").arg(useridx).arg(-1);
//        sendmsgs({tobeSendMsg});
//        return;
//    }
    map<string,set<int>> specPointsMap=getColorChangedPoints(segments);
    set<string> resultSet;
    for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
        int size = 0;
        for(auto it2=it->second.begin(); it2!=it->second.end(); it2++){
            if(*it2!=2 && *it2!=3 && *it2!=4 && *it2!=1)
                size++;
        }
        if(size!=0){
            resultSet.insert(it->first);
        }
    }
//    for(auto it=resultSet.begin();it!=resultSet.end();){
//        NeuronSWC s;
//        stringToXYZ(*it, s.x, s.y, s.z);
//        if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y,
//                     s.z, somaCoordinate.z)<dist_thre)
//        {
//            it=resultSet.erase(it);
//        }else{
//            it++;
//        }
//    }

    if(resultSet.size()!=0){
        logOut << "color mutation exists\n";
        for(auto it=resultSet.begin(); it!=resultSet.end(); it++){
            NeuronSWC s;
            stringToXYZ(*it, s.x, s.y, s.z);
            CellAPO marker;
            marker.name="";
            marker.comment="Color mutation";
            marker.orderinfo="";
            marker.color.r=200;
            marker.color.g=20;
            marker.color.b=0;
            marker.x=s.x;
            marker.y=s.y;
            marker.z=s.z;
            errorMarkers.push_back(marker);
        }
//        return errorMarkers;
    }

    set<string> otherTypesCoor = resultSet;
    resultSet.clear();

    int case_type=0;
    for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
        if(it->second.find(2)!=it->second.end() && it->second.find(3)!=it->second.end()){
            resultSet.insert(it->first);
        }
    }
    if(resultSet.size()!=1)
        result=false;
    else{
        string gridKey=*resultSet.begin();
        XYZ coor;
        stringToXYZ(gridKey, coor.x, coor.y, coor.z);
        if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y,
                     coor.z, somaCoordinate.z)<dist_thre)
            case_type=1;
        else
            case_type=2;
    }

    if(case_type==0){
        logOut << "axons not connected or the number of axons is more than 1\n";
    }
    else if(case_type==1){
        for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
            if(it->second.find(2)!=it->second.end() && it->second.find(4)!=it->second.end()){
                resultSet.insert(it->first);
            }
            if(it->second.find(4)!=it->second.end() && it->second.find(3)!=it->second.end()){
                resultSet.insert(it->first);
            }
        }
        if(resultSet.size()!=1)
            result=false;
    }
    else if(case_type==2){
        string gridKey;
        for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
            if(it->second.find(4)!=it->second.end() && it->second.find(3)!=it->second.end()){
                resultSet.insert(it->first);
                gridKey=it->first;
            }
        }
        if(resultSet.size()>2)
        {
            result=false;
        }
        else if(resultSet.size()==2){
            XYZ coor;
            stringToXYZ(gridKey, coor.x, coor.y, coor.z);
            if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y,
                         coor.z, somaCoordinate.z)>dist_thre)
                result=false;
        }

        int curCount=resultSet.size();
        for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
            if(it->second.find(4)!=it->second.end() && it->second.find(2)!=it->second.end()){
                resultSet.insert(it->first);
            }
        }
        if(curCount!=resultSet.size())
            result=false;

    }

//    for(auto it=resultSet.begin();it!=resultSet.end();){
//        NeuronSWC s;
//        stringToXYZ(*it, s.x, s.y, s.z);
//        if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y,
//                     s.z, somaCoordinate.z)<dist_thre)
//        {
//            it=resultSet.erase(it);
//        }else{
//            it++;
//        }
//    }

    if(result && errorMarkers.size() == 0){
        logOut << "no color mutation\n";
    }else{
        logOut << "color mutation exists\n";
        for(auto it=resultSet.begin(); it!=resultSet.end(); it++){
            if(otherTypesCoor.find(*it) != otherTypesCoor.end()){
                continue;
            }

            NeuronSWC s;
            stringToXYZ(*it, s.x, s.y, s.z);

            CellAPO marker;
            marker.name="";
            marker.comment="Color mutation";
            marker.orderinfo="";
            marker.color.r=200;
            marker.color.g=20;
            marker.color.b=0;
            marker.x=s.x;
            marker.y=s.y;
            marker.z=s.z;
            errorMarkers.push_back(marker);
        }
    }
    logOut << "analyzeColorMutation end\n";
    return errorMarkers;
}

vector<CellAPO> analyzeColorMutationForHB(QTextStream& logOut, bool isSomaExists, XYZ somaCoordinate, V_NeuronSWC_list segments, float dist_thre){
    logOut << "begin analyzeColorMutation...\n";
    vector<CellAPO> errorMarkers;
    if(!isSomaExists){
        return errorMarkers;
    }

    bool result=true;
    //    if(!myServer->isSomaExists){
    //        qDebug()<<"soma not detected!";
    //        QString tobeSendMsg="/FEEDBACK_ANALYZE_ColorMutation:";
    //        tobeSendMsg += QString("server %1 %2").arg(useridx).arg(-1);
    //        sendmsgs({tobeSendMsg});
    //        return;
    //    }
    map<string,set<int>> specPointsMap=getColorChangedPoints(segments);
    set<string> resultSet;
    for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
        int size = 0;
        for(auto it2=it->second.begin(); it2!=it->second.end(); it2++){
            if(*it2!=2 && *it2!=3 && *it2!=4 && *it2!=1)
                size++;
        }
        if(size!=0){
            resultSet.insert(it->first);
        }
    }
//    for(auto it=resultSet.begin();it!=resultSet.end();){
//        NeuronSWC s;
//        stringToXYZ(*it, s.x, s.y, s.z);
//        if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y,
//                     s.z, somaCoordinate.z)<dist_thre)
//        {
//            it=resultSet.erase(it);
//        }else{
//            it++;
//        }
//    }
    if(resultSet.size()!=0){
        logOut << "color mutation exists\n";
        for(auto it=resultSet.begin(); it!=resultSet.end(); it++){
            NeuronSWC s;
            stringToXYZ(*it, s.x, s.y, s.z);
            CellAPO marker;
            marker.name="";
            marker.comment="Color mutation";
            marker.orderinfo="";
            marker.color.r=255;
            marker.color.g=128;
            marker.color.b=0;
            marker.x=s.x;
            marker.y=s.y;
            marker.z=s.z;
            errorMarkers.push_back(marker);
        }
//        return errorMarkers;
    }
    set<string> otherTypesCoor = resultSet;
    resultSet.clear();

    int case_type=0;
    for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
        if(it->second.find(2)!=it->second.end() && it->second.find(3)!=it->second.end()){
            resultSet.insert(it->first);
        }
    }
    if(resultSet.size() > 1)
        result=false;
    else if(resultSet.size() == 1){
        string gridKey=*resultSet.begin();
        XYZ coor;
        stringToXYZ(gridKey, coor.x, coor.y, coor.z);
        if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y,
                     coor.z, somaCoordinate.z)<dist_thre)
            case_type=1;
        else
            case_type=2;
    }
    else{
        case_type = 4;
    }

    if(case_type==0){
        logOut << "axons not connected or the number of axons is more than 1\n";
    }
    else if(case_type==1){
        for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
            if(it->second.find(2)!=it->second.end() && it->second.find(4)!=it->second.end()){
                resultSet.insert(it->first);
            }
            if(it->second.find(4)!=it->second.end() && it->second.find(3)!=it->second.end()){
                resultSet.insert(it->first);
            }
        }
        if(resultSet.size()!=1)
            result=false;
    }
    else if(case_type==2){
        string gridKey;
        for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
            if(it->second.find(4)!=it->second.end() && it->second.find(3)!=it->second.end()){
                resultSet.insert(it->first);
                gridKey=it->first;
            }
        }
        if(resultSet.size()>2)
        {
            result=false;
        }
        else if(resultSet.size()==2){
            XYZ coor;
            stringToXYZ(gridKey, coor.x, coor.y, coor.z);
            if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y,
                         coor.z, somaCoordinate.z)>dist_thre)
                result=false;
        }

        int curCount=resultSet.size();
        for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
            if(it->second.find(4)!=it->second.end() && it->second.find(2)!=it->second.end()){
                resultSet.insert(it->first);
            }
        }
        if(curCount!=resultSet.size())
            result=false;

    }
    else if(case_type == 4){
        string gridKey;
        for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
            if(it->second.find(2)!=it->second.end() && it->second.find(4)!=it->second.end()){
                resultSet.insert(it->first);
            }
        }
        if(resultSet.size()!=0)
            result=false;
        else{
            for(auto it=specPointsMap.begin();it!=specPointsMap.end();it++){
                if(it->second.find(3)!=it->second.end() && it->second.find(4)!=it->second.end()){
                    resultSet.insert(it->first);
                    gridKey=it->first;
                }
            }
            if(resultSet.size()==0){
                result = true;
            }
            else{
                if(resultSet.size() == 1){
                    XYZ coor;
                    stringToXYZ(gridKey, coor.x, coor.y, coor.z);
                    if(distance(coor.x, somaCoordinate.x, coor.y, somaCoordinate.y,
                                 coor.z, somaCoordinate.z)>dist_thre){
                        result = false;
                    }
                }
                else{
                    result = false;
                }
            }
        }
    }

//    for(auto it=resultSet.begin();it!=resultSet.end();){
//        NeuronSWC s;
//        stringToXYZ(*it, s.x, s.y, s.z);
//        if(distance(s.x, somaCoordinate.x, s.y, somaCoordinate.y,
//                     s.z, somaCoordinate.z)<dist_thre)
//        {
//            it=resultSet.erase(it);
//        }else{
//            it++;
//        }
//    }

    if(result && errorMarkers.size() == 0){
        logOut << "no color mutation\n";
    }else{
        logOut << "color mutation exists\n";
        for(auto it=resultSet.begin(); it!=resultSet.end(); it++){
            if(otherTypesCoor.find(*it) != otherTypesCoor.end()){
                continue;
            }
            NeuronSWC s;
            stringToXYZ(*it, s.x, s.y, s.z);

            CellAPO marker;
            marker.name="";
            marker.comment="Color mutation";
            marker.orderinfo="";
            marker.color.r=255;
            marker.color.g=128;
            marker.color.b=0;
            marker.x=s.x;
            marker.y=s.y;
            marker.z=s.z;
            errorMarkers.push_back(marker);
        }
    }
    logOut << "analyzeColorMutation end\n";
    return errorMarkers;
}

vector<CellAPO> analyzeDissociativeSegs(QTextStream& logOut, V_NeuronSWC_list segments, XYZ somaCoordinate){
    logOut << "begin analyzeDissociativeSegs...\n";
    vector<CellAPO> errorMarkers;
    set<string> dissociativePoints=getDissociativeSegsMarkerPoints(V_NeuronSWC_list__2__NeuronTree(segments).listNeuron, somaCoordinate);

    if(dissociativePoints.size()==0){
        logOut << "no dissociative segs\n";
    }else{
        logOut << "dissociative seg exists\n";
        for(auto it=dissociativePoints.begin(); it!=dissociativePoints.end(); it++){
            NeuronSWC s;
            stringToXYZ(*it, s.x, s.y, s.z);

            CellAPO marker;
            marker.name="";
            marker.comment="Isolated branch";
            marker.orderinfo="";
            marker.color.r=200;
            marker.color.g=20;
            marker.color.b=0;
            marker.x=s.x;
            marker.y=s.y;
            marker.z=s.z;
            errorMarkers.push_back(marker);
        }
    }
    logOut << "analyzeDissociativeSegs end\n";
    return errorMarkers;
}

set<string> getDissociativeSegsMarkerPoints(QList<NeuronSWC> neuron, XYZ somaCoordinate){
    return getTreeMarkerPoints(neuron, somaCoordinate);
}

vector<CellAPO> analyzeAngles(QTextStream& logOut, XYZ somaCoordinate, V_NeuronSWC_list& segments, float dist_thre, bool isSomaExists){
    logOut << "begin analyzeAngles...\n";
    vector<CellAPO> errorMarkers;

    set<string> angleErrPoints=getAngleErrPoints(logOut, dist_thre, somaCoordinate, segments, isSomaExists, false);

    if(angleErrPoints.size()==0){
        logOut << "no angle-error dendrite bifurcations\n";
    }else{
        logOut << "angle-error dendrite bifurcation exists\n";
        for(auto it=angleErrPoints.begin(); it!=angleErrPoints.end(); it++){
            NeuronSWC s;
            stringToXYZ(*it, s.x, s.y, s.z);

            CellAPO marker;
            marker.name="";
            marker.comment="Angle error";
            marker.orderinfo="";
            marker.color.r=0;
            marker.color.g=200;
            marker.color.b=20;
            marker.x=s.x;
            marker.y=s.y;
            marker.z=s.z;
            errorMarkers.push_back(marker);
        }

    }
    logOut << "analyzeAngles end\n";
    return errorMarkers;
}
