#include "colldetection.h"
#include "analyze.h"
#include "sort_swc.h"
#include <iostream>
#include <vector>
#include "detect_crossing/utilities.h"
#include "detect_crossing/SwcReader.h"
#include <mutex>
#include <filesystem>
#include "detect_crossing/CrossingDetect.h"
#include "detect_crossing/ResultWriter.h"
#include <iomanip>

QString CollDetection::swcInputDirPath="";

CollDetection::CollDetection(QString infilepath, QString infilename, QString logpath, QString apopath, QString anopath, QString swcoutputpath, QString tmpswcpath, QString resultpath,
                             QString croppedApomulfurPath, QString croppedApobifurPath, QString croppedApoloopPath, QString croppedApomissingPath,
                             QString croppedApocrossingPath, QString croppedApooverlapPath, QString croppedApodissoPath, QString croppedSwcmulfurPath,
                             QString croppedSwcbifurPath, QString croppedSwcloopPath, QString croppedSwcmissingPath,
                             QString croppedSwccrossingPath, QString croppedSwcoverlapPath, QString croppedSwcdissoPath, QObject* parent){
    accessManager=new QNetworkAccessManager(this);
    SuperUserHostAddress="http://114.117.165.134:26020/SuperUser";
    BrainTellHostAddress="http://114.117.165.134:26000/release";
    maxRes.x = maxRes.y = maxRes.z = 0;
    subMaxRes.x = subMaxRes.y = subMaxRes.z = 0;
    inFile=infilepath;
    inFilename=infilename;
    logPath=logpath;
    apoPath=apopath;
    anoPath=anopath;
    swcOutputPath=swcoutputpath;

    croppedApoMulfurPath = croppedApomulfurPath;
    croppedApoBifurPath = croppedApobifurPath;
    croppedApoLoopPath = croppedApoloopPath;
    croppedApoMissingPath = croppedApomissingPath;
    croppedApoCrossingPath = croppedApocrossingPath;
    croppedApoOverlapPath = croppedApooverlapPath;
    croppedApoDissoPath = croppedApodissoPath;
    croppedSwcMulfurPath = croppedSwcmulfurPath;
    croppedSwcBifurPath = croppedSwcbifurPath;
    croppedSwcLoopPath = croppedSwcloopPath;
    croppedSwcMissingPath = croppedSwcmissingPath;
    croppedSwcCrossingPath = croppedSwccrossingPath;
    croppedSwcOverlapPath = croppedSwcoverlapPath;
    croppedSwcDissoPath = croppedSwcdissoPath;

    tmpInFile = tmpswcpath;
    resultPath=resultpath;
    logFile=new QFile(logPath);
    if(!logFile->open(QIODevice::Append | QIODevice::Text)){
        qDebug() << "cannot open logfile";
    }
    logOut.setDevice(logFile);
    auto nt=readSWC_file(inFile);
    segments=NeuronTree__2__V_NeuronSWC_list(nt);
    isSomaExists = false;
    somaCoordinate = XYZ(0.0, 0.0, 0.0);
}

XYZ CollDetection::getSomaCoordinate(QString apoPath){
    logOut << "begin getSomaCoordinate...\n";
    isSomaExists=false;
    XYZ coordinate;
    coordinate.x=-1;
    coordinate.y=-1;
    coordinate.z=-1;
    QFile qf(apoPath);
    if (!qf.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug()<<"apofile open error";
        return coordinate;
    }
    char *buf;
    char _buf[1000];
    qf.readLine(_buf, sizeof(_buf));
    for (buf=_buf; (*buf && *buf==' '); buf++);
    if (buf[0]=='#' ||buf[0]=='\0')
    {
        if(!qf.atEnd())
        {
            qf.readLine(_buf, sizeof(_buf));
            for (buf=_buf; (*buf && *buf==' '); buf++);
        }
        else{
            qDebug()<<"apofile format error";
            return coordinate;
        }
    }
    else{
        qDebug()<<"apofile format error";
        return coordinate;
    }
    QStringList qsl = QString(buf).split(",");
    if (qsl.size()==0){
        qDebug()<<"apofile format error";
        return coordinate;
    }
    else{
        for (int i=4; i<qsl.size(); i++)
        {
            qsl[i].truncate(200); //change from 99 to 200, 20121212, by PHC
            if (i==4) coordinate.z = qsl[i].toFloat();
            if (i==5) coordinate.x = qsl[i].toFloat();
            if (i==6)
            {
                coordinate.y = qsl[i].toFloat();
                isSomaExists=true;
                break;
            }
        }
    }
    logOut << "getSomaCoordinate end\n";
    return coordinate;
}

void CollDetection::detectAll(){
    logOut << "------------------------------------\n";
    logOut << "Start Detection Time: ";
    // 获取当前时间
    QDateTime currentTime = QDateTime::currentDateTime();
    qint64 currentTimeStamp = currentTime.toSecsSinceEpoch();
    // 将时间格式化为字符串
    QString formattedTime = currentTime.toString("yyyy-MM-dd HH:mm:ss");
    logOut << formattedTime << "\n";
    logOut << "------------------------------------\n";

    logOut << "begin detectAll...\n";
    if(segments.name == "invalid_swc"){
        logOut << "invalid swc\n";
        logOut << "detectAll end\n";
        logOut << "------------------------------------\n";
        return;
    }
    if(segments.seg.size() == 0)
    {
        logOut << "swc is empty!\n";
        logOut << "detectAll end\n";
        logOut << "------------------------------------\n";
        return;
    }
    //尝试从swc中获取soma坐标
    QStringList swcNodeList = getSWCSpecNInfo(inFile, -1);
    if(swcNodeList.size() != 1){
        isSortedSwc = false;
        isSomaExists = false;
    }
    else{
        auto swcNodeInfo = swcNodeList[0].split(' ',Qt::SkipEmptyParts);
        isSortedSwc = true;
        isSomaExists = true;
        somaCoordinate.x = swcNodeInfo[2].toFloat();
        somaCoordinate.y = swcNodeInfo[3].toFloat();
        somaCoordinate.z = swcNodeInfo[4].toFloat();
        if(fabs(somaCoordinate.x) < 1e-5 || fabs(somaCoordinate.y) < 1e-5 || fabs(somaCoordinate.z) < 1e-5){
            logOut << "invalid swc\n";
            logOut << "detectAll end\n";
            logOut << "------------------------------------\n";
            return;
        }
    }

    getImageRES();
//    removeErrorSegs(segments);
    tuneErrorSegs(segments);

    writeESWC_file(swcOutputPath, V_NeuronSWC_list__2__NeuronTree(segments));

    auto nt=readSWC_file(swcOutputPath);
    segments=NeuronTree__2__V_NeuronSWC_list(nt);

    detectOthers();
    detectOverlapSegs(segments);

    colorMutationMarkers = analyzeColorMutationForHB(logOut, isSomaExists, somaCoordinate, segments, 8);
    addMarkers(colorMutationMarkers);
    dissociativeSegsMarkers = analyzeDissociativeSegs(logOut, segments, somaCoordinate);
    addMarkers(dissociativeSegsMarkers);
    anglesMarkers = analyzeAngles(logOut, somaCoordinate, segments, 8, isSomaExists);
    addMarkers(anglesMarkers);

    bool isSuccess = sortSWCAndDetectLoop(swcOutputPath, tmpInFile);
    if(!isSuccess){
        sortSWC(swcOutputPath,tmpInFile,0);
    }
    setSWCRadius(tmpInFile,1);

    detectTips();
//    detectCrossings();

    QFile::remove(tmpInFile);
    logOut << "detectAll end\n";
    logOut << "------------------------------------\n";
    logOut << "End Detection Time: ";
    // 获取当前时间
    QDateTime currentTime2 = QDateTime::currentDateTime();
    qint64 currentTimeStamp2 = currentTime2.toSecsSinceEpoch();
    // 将时间格式化为字符串
    QString formattedTime2 = currentTime2.toString("yyyy-MM-dd HH:mm:ss");
    logOut << formattedTime2 << "\n";
    logOut << "Consumed Time: " << currentTimeStamp2 - currentTimeStamp << "s\n";
    logOut << "------------------------------------\n";
}

void CollDetection::generateResult(){
    if(segments.seg.size() == 0 || segments.name == "invalid_swc")
    {
        return;
    }

    logOut << "begin generateResult...\n";

    std::ofstream resultOut(resultPath.toStdString());

    resultOut << "胞体位置检测情况: ";
    if(isSortedSwc){
        resultOut << "已检测出胞体位置";
    }
    else if(isSomaExists){
        resultOut << "自动设置多分叉为胞体位置";
    }
    else{
        resultOut << "未检测出胞体位置";
    }
    resultOut << endl;
    resultOut << "错误类型\t\t";
    resultOut << "数量\t\t";
    resultOut << "marker颜色\t\t";
    resultOut << "状态\t" << endl;
    resultOut << "多分叉\t\t\t";
    resultOut << mulFurcationMarkers.size() << "\t\t";
    resultOut << "棕色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "临近二分叉\t\t";
    resultOut << nearBifurcationMarkers.size()/2 << "\t\t";
    resultOut << "黄色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "环\t\t\t";
    resultOut << (loopMarkers.size()+1)/2 << "\t\t";
    resultOut << "白色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "角度异常\t\t";
    resultOut << anglesMarkers.size() << "\t\t";
    resultOut << "绿色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "末端缺失\t\t";
    resultOut << tipUndoneMarkers.size() << "\t\t";
    resultOut << "粉色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "交叉错误\t\t";
    resultOut << crossingMarkers.size() << "\t\t";
    resultOut << "淡紫色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "重叠线段\t\t";
    resultOut << overlapSegsMarkers.size()/2 << "对\t\t";
    resultOut << "蓝色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "漂浮分支\t\t";
    resultOut << dissociativeSegsMarkers.size() << "\t\t";
    resultOut << "红色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "颜色突变\t\t";
    resultOut << colorMutationMarkers.size() << "\t\t";
    resultOut << "橙色\t\t\t";
    resultOut << "已标记\t" << endl;
    resultOut << "异常分支\t\t";
    resultOut << errorSegsNum << "\t\t";
    resultOut << "无\t\t\t";
    resultOut << "已去除\t" << endl;
    resultOut.close();

    QString csvPath = resultPath.left(resultPath.size() - 4) + ".csv";
    std::ofstream resultCSVOut(csvPath.toStdString());
    resultCSVOut << "文件目录," << "胞体位置检测情况," << "多分叉(棕色)," << "邻近二分叉(黄色)," << "环(白色)," << "角度异常(绿色)," << "末端缺失(粉色)," << "交叉错误(淡紫色)," << "重叠线段(蓝色)," << "漂浮分支(红色)," << "颜色突变(橙色)," << "异常分支(无)" << endl;
    //获取相对路径
    QDir swcInputDir(swcInputDirPath);
    QString relativePath = swcInputDir.relativeFilePath(inFile);
    resultCSVOut << "\"" <<relativePath.toStdString() << "\"" << ",";
//    if(inFile.endsWith(".swc"))
//        resultCSVOut << inFile.toStdString() << ".swc,";
//    if(inFile.endsWith(".eswc"))
//        resultCSVOut << inFilename.toStdString() << ".eswc,";
    if(isSortedSwc){
        resultCSVOut << "已检测出胞体位置";
    }
    else if(isSomaExists){
        resultCSVOut << "自动设置多分叉为胞体位置";
    }
    else{
        resultCSVOut << "未检测出胞体位置";
    }
    resultCSVOut << ",";
    resultCSVOut << mulFurcationMarkers.size() << ",";
    resultCSVOut << nearBifurcationMarkers.size()/2 << ",";
    resultCSVOut << (loopMarkers.size()+1)/2 << ",";
    resultCSVOut << anglesMarkers.size() << ",";
    resultCSVOut << tipUndoneMarkers.size() << ",";
    resultCSVOut << crossingMarkers.size() << ",";
    resultCSVOut << overlapSegsMarkers.size()/2<< "对,";
    resultCSVOut << dissociativeSegsMarkers.size() << ",";
    resultCSVOut << colorMutationMarkers.size() << ",";
    resultCSVOut << errorSegsNum << endl;
    resultCSVOut.close();

    writeAPO_file(apoPath, markers);
    writeESWC_file(swcOutputPath, V_NeuronSWC_list__2__NeuronTree(segments));
    QFile* anoFile = new QFile(anoPath);
    anoFile->open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream anoOut(anoFile);
    anoOut << "APOFILE=" << inFilename << ".apo\n";
    if(inFile.endsWith(".swc"))
        anoOut << "SWCFILE=" << inFilename << ".swc\n";
    if(inFile.endsWith(".eswc"))
        anoOut << "SWCFILE=" << inFilename << ".eswc\n";
    anoFile->flush();
    anoFile->close();
    delete anoFile;

    logOut << "generateResult end\n";
    logOut << "\n\n";
}

void CollDetection::getApoForCrop(QString fileSaveName, vector<CellAPO> points){
    QList <CellAPO> markers;
    for(int i=0;i<points.size();i++){
        CellAPO marker;
        marker.color.r=points[i].color.r;
        marker.color.g=points[i].color.g;
        marker.color.b=points[i].color.b;
        marker.x=points[i].x;
        marker.y=points[i].y;
        marker.z=points[i].z;
        markers.append(marker);
    }

    writeAPO_file(fileSaveName, markers);
}

void CollDetection::getCropedSwc(QString fileInputPath, QString fileSavePath, XYZ coor1, XYZ coor2){

    NeuronTree nt=readSWC_file(fileInputPath);
    V_NeuronSWC_list segments=NeuronTree__2__V_NeuronSWC_list(nt);

    V_NeuronSWC_list tosave;
    for(std::vector<V_NeuronSWC_unit>::size_type i=0;i<segments.seg.size();i++)
    {
        NeuronTree SS;
        const V_NeuronSWC &seg_temp =  segments.seg.at(i);
        for(std::vector<V_NeuronSWC_unit>::size_type j=0;j<seg_temp.row.size();j++)
        {
            if(seg_temp.row.at(j).x>=coor1.x&&seg_temp.row.at(j).x<=coor2.x
                &&seg_temp.row.at(j).y>=coor1.y&&seg_temp.row.at(j).y<=coor2.y
                &&seg_temp.row.at(j).z>=coor1.z&&seg_temp.row.at(j).z<=coor2.z)
            {
                tosave.seg.push_back(seg_temp);
                break;
            }
        }
    }

    NeuronTree savent=V_NeuronSWC_list__2__NeuronTree(tosave);
    writeESWC_file(fileSavePath,savent);
}

void CollDetection::convertCoordInCropedSwc(QString filePath, XYZ coor1){
    if (filePath.endsWith(".swc") || filePath.endsWith(".SWC") || filePath.endsWith(".eswc") || filePath.endsWith(".ESWC"))
    {
        QFile qf(filePath);
        QString arryRead;
        if(!qf.open(QIODevice::ReadOnly|QIODevice::Text)){
            return;
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
            qDebug()<<"swcfile cannot be opened!";
            return;
        }
        QTextStream streamWrite(&qf);
        for(int i=0;i<arryListWrite.size()-1;i++){      //这里到arryListWrite.size()-1是因为arryListWrite数组按照\n分段时，最后一行尾部有个\n，所以数组最后一个值为空，需要将它去掉
            if(arryListWrite.at(i).contains("#")){
                streamWrite<<arryListWrite.at(i)<<"\n";
            }else{
                QString contentWrite= arryListWrite.at(i);
                QStringList swcInfo=contentWrite.split(' ',Qt::SkipEmptyParts);
                double x=swcInfo[2].toDouble();
                double y=swcInfo[3].toDouble();
                double z=swcInfo[4].toDouble();
                x-=coor1.x;
                y-=coor1.y;
                z-=coor1.z;

                swcInfo[2]=QString::number(x, 'f', 3);
                swcInfo[3]=QString::number(y, 'f', 3);
                swcInfo[4]=QString::number(z, 'f', 3);
                contentWrite=swcInfo.join(' ');
                streamWrite<<contentWrite<<"\n";
            }
        }
        qf.close();
    }
}

void CollDetection::getApoAndCroppedSwc(){
    setSWCRadius(swcOutputPath, 1);

    getApoForCrop(croppedApoMulfurPath, mulFurcationMarkers);
    getApoForCrop(croppedApoBifurPath, nearBifurcationMarkers);
    getApoForCrop(croppedApoLoopPath, loopMarkers);
    getApoForCrop(croppedApoMissingPath, tipAllMarkers);
    getApoForCrop(croppedApoCrossingPath, crossingAllMarkers);
    getApoForCrop(croppedApoOverlapPath, overlapSegsMarkers);
    getApoForCrop(croppedApoDissoPath, dissociativeSegsMarkers);

    int size = 32;
    int z_size = 10;

    for(auto s=mulFurcationMarkers.begin(); s!=mulFurcationMarkers.end(); s++){
        int x=int(s->x);
        int y=int(s->y);
        int z=int(s->z);
       // XYZ coor1=XYZ(max(x-size, 0), max(y-size, 0), max(z-size, 0));
       // XYZ coor2=XYZ(min(x+size, int(maxRes.x - 1)), min(y+size, int(maxRes.y - 1)), min(z+size, int(maxRes.z - 1)));

        XYZ coor1=XYZ(x-size, y-size, z-z_size);
        XYZ coor2=XYZ(x+size, y+size, z+z_size);

        QString coorName=QString::number(x)+"_"+QString::number(y)+"_"+QString::number(z);
        QString cropSwcDirPath = croppedSwcMulfurPath + "/" + coorName;
        QDir dir;
        if(!dir.exists(cropSwcDirPath))
            dir.mkpath(cropSwcDirPath);
        getCropedSwc(swcOutputPath, cropSwcDirPath+"/optical.eswc", coor1, coor2);
        convertCoordInCropedSwc(cropSwcDirPath+"/optical.eswc", coor1);
    }

    for(auto s=nearBifurcationMarkers.begin(); s!=nearBifurcationMarkers.end(); s++){
        int x=int(s->x);
        int y=int(s->y);
        int z=int(s->z);
        // XYZ coor1=XYZ(max(x-size, 0), max(y-size, 0), max(z-size, 0));
        // XYZ coor2=XYZ(min(x+size, int(maxRes.x - 1)), min(y+size, int(maxRes.y - 1)), min(z+size, int(maxRes.z - 1)));

        XYZ coor1=XYZ(x-size, y-size, z-z_size);
        XYZ coor2=XYZ(x+size, y+size, z+z_size);

        QString coorName=QString::number(x)+"_"+QString::number(y)+"_"+QString::number(z);
        QString cropSwcDirPath = croppedSwcBifurPath + "/" + coorName;
        QDir dir;
        if(!dir.exists(cropSwcDirPath))
            dir.mkpath(cropSwcDirPath);
        getCropedSwc(swcOutputPath, cropSwcDirPath+"/optical.eswc", coor1, coor2);
        convertCoordInCropedSwc(cropSwcDirPath+"/optical.eswc", coor1);
    }

    for(auto s=loopMarkers.begin(); s!=loopMarkers.end(); s++){
        int x=int(s->x);
        int y=int(s->y);
        int z=int(s->z);
        // XYZ coor1=XYZ(max(x-size, 0), max(y-size, 0), max(z-size, 0));
        // XYZ coor2=XYZ(min(x+size, int(maxRes.x - 1)), min(y+size, int(maxRes.y - 1)), min(z+size, int(maxRes.z - 1)));

        XYZ coor1=XYZ(x-size, y-size, z-z_size);
        XYZ coor2=XYZ(x+size, y+size, z+z_size);

        QString coorName=QString::number(x)+"_"+QString::number(y)+"_"+QString::number(z);
        QString cropSwcDirPath = croppedSwcLoopPath + "/" + coorName;
        QDir dir;
        if(!dir.exists(cropSwcDirPath))
            dir.mkpath(cropSwcDirPath);
        getCropedSwc(swcOutputPath, cropSwcDirPath+"/optical.eswc", coor1, coor2);
        convertCoordInCropedSwc(cropSwcDirPath+"/optical.eswc", coor1);
    }

    for(auto s=tipAllMarkers.begin(); s!=tipAllMarkers.end(); s++){
        int x=int(s->x);
        int y=int(s->y);
        int z=int(s->z);
        // XYZ coor1=XYZ(max(x-size, 0), max(y-size, 0), max(z-size, 0));
        // XYZ coor2=XYZ(min(x+size, int(maxRes.x - 1)), min(y+size, int(maxRes.y - 1)), min(z+size, int(maxRes.z - 1)));

        XYZ coor1=XYZ(x-size, y-size, z-z_size);
        XYZ coor2=XYZ(x+size, y+size, z+z_size);

        QString coorName=QString::number(x)+"_"+QString::number(y)+"_"+QString::number(z);
        QString cropSwcDirPath = croppedSwcMissingPath + "/" + coorName;
        QDir dir;
        if(!dir.exists(cropSwcDirPath))
            dir.mkpath(cropSwcDirPath);
        getCropedSwc(swcOutputPath, cropSwcDirPath+"/optical.eswc", coor1, coor2);
        convertCoordInCropedSwc(cropSwcDirPath+"/optical.eswc", coor1);
    }

    for(auto s=crossingAllMarkers.begin(); s!=crossingAllMarkers.end(); s++){
        int x=int(s->x);
        int y=int(s->y);
        int z=int(s->z);
        // XYZ coor1=XYZ(max(x-size, 0), max(y-size, 0), max(z-size, 0));
        // XYZ coor2=XYZ(min(x+size, int(maxRes.x - 1)), min(y+size, int(maxRes.y - 1)), min(z+size, int(maxRes.z - 1)));

        XYZ coor1=XYZ(x-size, y-size, z-z_size);
        XYZ coor2=XYZ(x+size, y+size, z+z_size);

        QString coorName=QString::number(x)+"_"+QString::number(y)+"_"+QString::number(z);
        QString cropSwcDirPath = croppedSwcCrossingPath + "/" + coorName;
        QDir dir;
        if(!dir.exists(cropSwcDirPath))
            dir.mkpath(cropSwcDirPath);
        getCropedSwc(swcOutputPath, cropSwcDirPath+"/optical.eswc", coor1, coor2);
        convertCoordInCropedSwc(cropSwcDirPath+"/optical.eswc", coor1);
    }

    for(auto s=overlapSegsMarkers.begin(); s!=overlapSegsMarkers.end(); s++){
        int x=int(s->x);
        int y=int(s->y);
        int z=int(s->z);
        // XYZ coor1=XYZ(max(x-size, 0), max(y-size, 0), max(z-size, 0));
        // XYZ coor2=XYZ(min(x+size, int(maxRes.x - 1)), min(y+size, int(maxRes.y - 1)), min(z+size, int(maxRes.z - 1)));

        XYZ coor1=XYZ(x-size, y-size, z-z_size);
        XYZ coor2=XYZ(x+size, y+size, z+z_size);

        QString coorName=QString::number(x)+"_"+QString::number(y)+"_"+QString::number(z);
        QString cropSwcDirPath = croppedSwcOverlapPath + "/" + coorName;
        QDir dir;
        if(!dir.exists(cropSwcDirPath))
            dir.mkpath(cropSwcDirPath);
        getCropedSwc(swcOutputPath, cropSwcDirPath+"/optical.eswc", coor1, coor2);
        convertCoordInCropedSwc(cropSwcDirPath+"/optical.eswc", coor1);
    }

    for(auto s=dissociativeSegsMarkers.begin(); s!=dissociativeSegsMarkers.end(); s++){
        int x=int(s->x);
        int y=int(s->y);
        int z=int(s->z);
        // XYZ coor1=XYZ(max(x-size, 0), max(y-size, 0), max(z-size, 0));
        // XYZ coor2=XYZ(min(x+size, int(maxRes.x - 1)), min(y+size, int(maxRes.y - 1)), min(z+size, int(maxRes.z - 1)));

        XYZ coor1=XYZ(x-size, y-size, z-z_size);
        XYZ coor2=XYZ(x+size, y+size, z+z_size);

        QString coorName=QString::number(x)+"_"+QString::number(y)+"_"+QString::number(z);
        QString cropSwcDirPath = croppedSwcDissoPath + "/" + coorName;
        QDir dir;
        if(!dir.exists(cropSwcDirPath))
            dir.mkpath(cropSwcDirPath);
        getCropedSwc(swcOutputPath, cropSwcDirPath+"/optical.eswc", coor1, coor2);
        convertCoordInCropedSwc(cropSwcDirPath+"/optical.eswc", coor1);
    }
}

void CollDetection::detectTips(){
    logOut << "begin detectTips...\n";
    map<string, set<size_t>> allPoint2SegIdMap = getWholeGrid2SegIDMap(segments);
    tipPoints=tipDetection(segments, true, allPoint2SegIdMap, 30);

    allPoint2SegIdMap = getWholeGrid2SegIDMap(segments);
    tipPoints=tipDetection(segments, false, allPoint2SegIdMap, 30);

    writeESWC_file(tmpInFile, V_NeuronSWC_list__2__NeuronTree(segments));
    handleTip(tipPoints);
    logOut << "detectTips end\n";
    qDebug() << "detectTips end\n";
}

void CollDetection::detectCrossings(){
    logOut<<"begin detectCrossings...\n";
    QJsonArray infos = crossingDetection();
    handleCrossing(infos);
    logOut<<"detectCrossings end\n";
    qDebug() << "detectCrossings end\n";
}

void CollDetection::detectOthers(){
    logOut<<"begin detectOthers...\n";
    vector<NeuronSWC> outputSpecialPoints = specStructsDetection(segments);
    vector<NeuronSWC> bifurPoints;
    vector<NeuronSWC> mulfurPoints;

    for(int i=0;i<outputSpecialPoints.size();i++){
        if(outputSpecialPoints[i].type == 6)
            bifurPoints.push_back(outputSpecialPoints[i]);
        else if(outputSpecialPoints[i].type == 8)
            mulfurPoints.push_back(outputSpecialPoints[i]);
    }

    handleMulFurcation(mulfurPoints, 8);
    handleNearBifurcation(bifurPoints);
    logOut<<"detectOthers end\n";
}

void CollDetection::detectLoops(){
    logOut<<"begin detectLoops...\n";
    vector<NeuronSWC> outputSpecialPoints = loopDetection(segments);
    handleLoop(outputSpecialPoints);
    logOut<<"detectLoops end\n";
}

vector<NeuronSWC> CollDetection::specStructsDetection(V_NeuronSWC_list& inputSegList, double adjacent_dist_thre, double soma_dist_thre){
    logOut << "begin specStructsDetection...\n";
    vector<NeuronSWC> outputSpecialPoints;
    if(inputSegList.seg.size() == 0)
        return outputSpecialPoints;

    map<string, set<size_t> > wholeGrid2segIDmap;
    map<string, bool> isEndPointMap;
    map<string, set<string>> parentMap;

    set<string> allPoints;
    map<string, set<string>> childMap;

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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
            wholeGrid2segIDmap[gridKey].insert(size_t(i));
            allPoints.insert(gridKey);

            if(seg.row[j].parent!=-1){
                float x2Label=seg.row[rowN2Index[seg.row[j].parent]].x;
                float y2Label=seg.row[rowN2Index[seg.row[j].parent]].y;
                float z2Label=seg.row[rowN2Index[seg.row[j].parent]].z;
                QString parentKeyQ=QString::number(x2Label) + "_" + QString::number(y2Label) + "_" + QString::number(z2Label);
                string parentKey=parentKeyQ.toStdString();
                parentMap[gridKey].insert(parentKey);
                childMap[parentKey].insert(gridKey);
            }

            if(j == 0 || j == seg.row.size() - 1){
                isEndPointMap[gridKey] = true;
            }
        }
    }

    //末端点和分叉点
    vector<string> points;
    vector<set<int>> linksIndex;
    //    vector<vector<int>> linksIndexVec;
    map<string,int> pointsIndexMap;

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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
                if(wholeGrid2segIDmap[gridKey].size()>1 &&
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

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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
        for(size_t j=0; j<segIndexs.size()-1; ++j){
            linksIndex[segIndexs[j]].insert(segIndexs[j+1]);
            //            linksIndexVec[segIndexs[j]].push_back(segIndexs[j+1]);
            linksIndex[segIndexs[j+1]].insert(segIndexs[j]);
            //            linksIndexVec[segIndexs[j+1]].push_back(segIndexs[j]);
        }
    }

    int maxFurcationsNum = 0;
    int maxFurcationIndex = -1;
    for(size_t i=0; i<points.size(); ++i){
//        qDebug()<<i<<" link size: "<<linksIndex[i].size();
        if(linksIndex[i].size() > 3){
            if(maxFurcationsNum < linksIndex[i].size()){
                maxFurcationsNum = linksIndex[i].size();
                maxFurcationIndex = i;
            }
            logOut << i << " link size: " << linksIndex[i].size() << "\n";
            NeuronSWC s;
            stringToXYZ(points[i],s.x,s.y,s.z);
            s.type = 8;
            if(isSomaExists){
                if(distance(s.x,somaCoordinate.x,s.y,somaCoordinate.y,s.z,somaCoordinate.z)>soma_dist_thre)
                    outputSpecialPoints.push_back(s);
            }
            else
                outputSpecialPoints.push_back(s);
        }
    }

    if(!isSomaExists && maxFurcationsNum >= 5){
        NeuronSWC s;
        stringToXYZ(points[maxFurcationIndex],s.x,s.y,s.z);
        isSomaExists = true;
        somaCoordinate.x = s.x;
        somaCoordinate.y = s.y;
        somaCoordinate.z = s.z;
        for(auto it = outputSpecialPoints.begin(); it != outputSpecialPoints.end(); it++){
            if(fabs(s.x - it->x) < 1e-5 && fabs(s.y - it->y) < 1e-5 && fabs(s.z - it->z) < 1e-5){
                outputSpecialPoints.erase(it);
                break;
            }
        }
    }

    vector<vector<size_t>> pairs;
    set<size_t> pset;

    size_t pre_tip_id=-1;
    size_t cur_tip_id=-1;

    double soma_radius=30;
    for(size_t i=0; i<points.size(); i++){
        if(linksIndex[i].size() == 3){
            pre_tip_id=cur_tip_id;
            cur_tip_id=i;
            if(pre_tip_id!=-1){
                NeuronSWC n1;
                stringToXYZ(points[pre_tip_id],n1.x,n1.y,n1.z);
                n1.type=6;
                NeuronSWC n2;
                stringToXYZ(points[cur_tip_id],n2.x,n2.y,n2.z);
                n2.type=6;
                set<size_t> n1Segs=wholeGrid2segIDmap[points[pre_tip_id]];
                set<size_t> n2Segs=wholeGrid2segIDmap[points[cur_tip_id]];
                int count1=0,count2=0;
                for(auto it1=n1Segs.begin();it1!=n1Segs.end();it1++)
                {
                    if(getSegLength(inputSegList.seg[*it1])>40)
                        count1++;
                }

                for(auto it2=n2Segs.begin();it2!=n2Segs.end();it2++)
                {
                    if(getSegLength(inputSegList.seg[*it2])>40)
                        count2++;
                }
                if(!(count1>=2&&count2>=2)){
                    continue;
                }
                if(isSomaExists){
                    if(distance(n1.x,somaCoordinate.x,n1.y,somaCoordinate.y,n1.z,somaCoordinate.z)>soma_radius
                        &&distance(n2.x,somaCoordinate.x,n2.y,somaCoordinate.y,n2.z,somaCoordinate.z)>soma_radius){
                        double dist=distance(n1.x,n2.x,n1.y,n2.y,n1.z,n2.z);
                        if(distance((n1.x+n2.x)/2,somaCoordinate.x,(n1.y+n2.y)/2,somaCoordinate.y,(n1.z+n2.z)/2,somaCoordinate.z)>1e-7&&dist<adjacent_dist_thre){
                            vector<size_t> v={pre_tip_id,cur_tip_id};
                            pairs.push_back(v);
                            pset.insert(pre_tip_id);
                            pset.insert(cur_tip_id);
                        }
                    }
                }
                else{
                    double dist=distance(n1.x,n2.x,n1.y,n2.y,n1.z,n2.z);
                    if(dist<adjacent_dist_thre){
                        vector<size_t> v={pre_tip_id,cur_tip_id};
                        pairs.push_back(v);
                        pset.insert(pre_tip_id);
                        pset.insert(cur_tip_id);
                    }
                }

            }
        }
    }

    for(auto it=pset.begin(); it!=pset.end(); it++){
        logOut << QString::fromStdString(points[*it]) << "\n";
        NeuronSWC n;
        stringToXYZ(points[*it],n.x,n.y,n.z);
        n.type=6;
        outputSpecialPoints.push_back(n);
    }

    logOut << "specStructsDetection end\n";
    return outputSpecialPoints;
}

vector<NeuronSWC> CollDetection::loopDetection(V_NeuronSWC_list& inputSegList){
    logOut << "begin loopDetection...\n";
    vector<NeuronSWC> outputSpecialPoints;

    map<string, set<size_t> > wholeGrid2segIDmap;
    map<string, bool> isEndPointMap;
    map<string, set<string>> parentMap;

    set<string> allPoints;
    map<string, set<string>> childMap;

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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
            wholeGrid2segIDmap[gridKey].insert(size_t(i));
            allPoints.insert(gridKey);

            if(seg.row[j].parent!=-1){
                float x2Label=seg.row[rowN2Index[seg.row[j].parent]].x;
                float y2Label=seg.row[rowN2Index[seg.row[j].parent]].y;
                float z2Label=seg.row[rowN2Index[seg.row[j].parent]].z;
                QString parentKeyQ=QString::number(x2Label) + "_" + QString::number(y2Label) + "_" + QString::number(z2Label);
                string parentKey=parentKeyQ.toStdString();
                parentMap[gridKey].insert(parentKey);
                childMap[parentKey].insert(gridKey);
            }

            if(j == 0 || j == seg.row.size() - 1){
                isEndPointMap[gridKey] = true;
            }
        }
    }

    //末端点和分叉点
    vector<string> points;
    vector<set<int>> linksIndex;
    //    vector<vector<int>> linksIndexVec;
    map<string,int> pointsIndexMap;

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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
                if(wholeGrid2segIDmap[gridKey].size()>1 &&
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

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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
            linksIndex[segIndexs[j]].insert(segIndexs[j+1]);
            //            linksIndexVec[segIndexs[j]].push_back(segIndexs[j+1]);
            linksIndex[segIndexs[j+1]].insert(segIndexs[j]);
            //            linksIndexVec[segIndexs[j+1]].push_back(segIndexs[j]);
        }
    }

    bool isDeleteEnd = false;
    while(!isDeleteEnd){
        isDeleteEnd = true;
        for(int i=0; i<points.size(); ++i){
            if(linksIndex[i].size() == 1){
                int linkIndex = *(linksIndex[i].begin());
                linksIndex[i].clear();
                linksIndex[linkIndex].erase(std::find(linksIndex[linkIndex].begin(),linksIndex[linkIndex].end(),i));
                isDeleteEnd = false;
            }
        }
    }

    //检测3条及3条以上的边构成的环

    vector<string> newpoints;

    for(size_t i=0; i<points.size(); ++i){
        if(linksIndex[i].size()>=2)
            newpoints.push_back(points[i]);
    }

    set<string> specPoints;
    size_t start=0;
    for(size_t i=0; i<newpoints.size(); ++i){
//        qDebug()<<QString::fromStdString(newpoints[i])<<" "<<parentMap[newpoints[i]].size();
        /*if(newLinksIndexVec[i].size()>=2&&counts[i]>=3&&newLinksIndexVec[i].size()!=counts[i])*/
        if(parentMap[newpoints[i]].size()>=2){
            specPoints.insert(newpoints[i]);
            loopNum++;
            logOut << QString::fromStdString(newpoints[i]) << "\n";
//            start=i+1;
//            qDebug()<<"loop exists";
        }

    }

    //检测2条边构成的环
    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
//        if(seg.row.size()<4)
//            continue;
        float xLabel1 = seg.row[0].x;
        float yLabel1 = seg.row[0].y;
        float zLabel1 = seg.row[0].z;
        float xLabel2=seg.row[seg.row.size()-1].x;
        float yLabel2=seg.row[seg.row.size()-1].y;
        float zLabel2=seg.row[seg.row.size()-1].z;
        QString gridKeyQ1 = QString::number(xLabel1) + "_" + QString::number(yLabel1) + "_" + QString::number(zLabel1);
        string gridKey1 = gridKeyQ1.toStdString();
        QString gridKeyQ2 = QString::number(xLabel2) + "_" + QString::number(yLabel2) + "_" + QString::number(zLabel2);
        string gridKey2 = gridKeyQ2.toStdString();
        set<size_t> segSet1=wholeGrid2segIDmap[gridKey1];
        set<size_t> segSet2=wholeGrid2segIDmap[gridKey2];
        set<size_t> intersectionSet;
        set_intersection(segSet1.begin(),segSet1.end(),segSet2.begin(),segSet2.end(),inserter( intersectionSet , intersectionSet.begin()));

        if(intersectionSet.size()>=2){
            logOut << "exists two-edge loop\n";
            logOut << "gridKey1: " << QString::fromStdString(gridKey1) << "\n";
            logOut << "gridKey2: " << QString::fromStdString(gridKey2) << "\n";
            loopNum++;
            specPoints.insert(gridKey1);
            specPoints.insert(gridKey2);
        }
    }

    for(auto it=specPoints.begin(); it!=specPoints.end(); it++){
        NeuronSWC s;
        stringToXYZ(*it,s.x,s.y,s.z);
        s.type=0;
        outputSpecialPoints.push_back(s);
    }

    logOut << "loopDetection end\n";
    return outputSpecialPoints;
}

vector<NeuronSWC> CollDetection::tipDetection(V_NeuronSWC_list& inputSegList, bool removeFlag, map<string, set<size_t>> allPoint2SegIdMap, double dist_thresh){
    logOut << "begin tipDetection...\n";

    vector<NeuronSWC> outputSpecialPoints;
    if(inputSegList.seg.size()==0 || maxRes.x == 0)
        return outputSpecialPoints;

    map<string, bool> isEndPointMap;
    map<string, set<size_t> > wholeGrid2segIDmap;

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];

        for(size_t j=0; j<seg.row.size(); ++j){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            wholeGrid2segIDmap[gridKey].insert(size_t(i));

            if(j == 0 || j == seg.row.size() - 1){
                isEndPointMap[gridKey] = true;
            }
        }
    }

    //末端点和分叉点
    vector<string> points;
    vector<set<int>> linksIndex;
    //    vector<vector<int>> linksIndexVec;
    map<string,int> pointsIndexMap;

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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
                if(wholeGrid2segIDmap[gridKey].size()>1 &&
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

    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        V_NeuronSWC seg = inputSegList.seg[i];
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

    //detect tips
    set<string> tips;

    for(int i=0;i<inputSegList.seg.size();i++){
        V_NeuronSWC seg = inputSegList.seg[i];
        float xLabel1 = seg.row[0].x;
        float yLabel1 = seg.row[0].y;
        float zLabel1 = seg.row[0].z;
        QString gridKeyQ1 = QString::number(xLabel1) + "_" + QString::number(yLabel1) + "_" + QString::number(zLabel1);
        string gridKey1 = gridKeyQ1.toStdString();
        float xLabel2 = seg.row[seg.row.size()-1].x;
        float yLabel2 = seg.row[seg.row.size()-1].y;
        float zLabel2 = seg.row[seg.row.size()-1].z;
        QString gridKeyQ2 = QString::number(xLabel2) + "_" + QString::number(yLabel2) + "_" + QString::number(zLabel2);
        string gridKey2 = gridKeyQ2.toStdString();
        if(wholeGrid2segIDmap[gridKey1].size()==1 && allPoint2SegIdMap[gridKey1].size()==1 && wholeGrid2segIDmap[gridKey2].size()>1)
        {
            if(isSomaExists&&sqrt((xLabel1-somaCoordinate.x)*(xLabel1-somaCoordinate.x)+
                (yLabel1-somaCoordinate.y)*(yLabel1-somaCoordinate.y)+(zLabel1-somaCoordinate.z)*(zLabel1-somaCoordinate.z))>50)
                tips.insert(gridKey1);
            else if(!isSomaExists)
                tips.insert(gridKey1);
        }
        if(wholeGrid2segIDmap[gridKey2].size()==1 && allPoint2SegIdMap[gridKey2].size()==1 && wholeGrid2segIDmap[gridKey1].size()>1)
        {
            if(isSomaExists&&sqrt((xLabel2-somaCoordinate.x)*(xLabel2-somaCoordinate.x)+
                (yLabel2-somaCoordinate.y)*(yLabel2-somaCoordinate.y)+(zLabel2-somaCoordinate.z)*(zLabel2-somaCoordinate.z))>50)
                tips.insert(gridKey2);
            else if(!isSomaExists)
                tips.insert(gridKey2);
        }
    }

    for(auto it=tips.begin();it!=tips.end();it++){
        vector<size_t> visitedSegIds;
        size_t segId=*wholeGrid2segIDmap[*it].begin();
        visitedSegIds.push_back(segId);
        V_NeuronSWC seg = inputSegList.seg[segId];
        float xLabel0 = seg.row[0].x;
        float yLabel0 = seg.row[0].y;
        float zLabel0 = seg.row[0].z;
        QString gridKeyQ0 = QString::number(xLabel0) + "_" + QString::number(yLabel0) + "_" + QString::number(zLabel0);
        string gridKey0 = gridKeyQ0.toStdString();
        float tipBranchLength=0;
        bool isReverse=false;
        if(wholeGrid2segIDmap[gridKey0].size()!=1)
        {
            isReverse=true;
        }
        bool flag=true;
        while(true){
            int size=seg.row.size();
            vector<int> indexs(size);
            for(int m=0;m<size;m++)
                indexs[m]=m;
            if(isReverse)
                reverse(indexs.begin(),indexs.end());
            for(int i=0;i<size;i++){
                int index=indexs[i];
                float xLabel = seg.row[index].x;
                float yLabel = seg.row[index].y;
                float zLabel = seg.row[index].z;
                QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
                string gridKey = gridKeyQ.toStdString();
                vector<string>::iterator it2=find(points.begin(),points.end(),gridKey);
                if(it2!=points.end()){
                    int index2=it2-points.begin();
                    if(linksIndex[index2].size()>=3){
                        flag=false;
                        break;
                    }
                    else{
                        if(index==seg.row.size()-1)
                            break;
                        tipBranchLength+=distance(xLabel,seg.row[index+1].x,
                                                    yLabel,seg.row[index+1].y,
                                                    zLabel,seg.row[index+1].z);
                        if(tipBranchLength>=dist_thresh)
                            break;
                        continue;
                    }
                }
                tipBranchLength+=distance(xLabel,seg.row[index+1].x,
                                            yLabel,seg.row[index+1].y,
                                            zLabel,seg.row[index+1].z);
                if(tipBranchLength>=dist_thresh)
                    break;
            }

            if(tipBranchLength>=dist_thresh||!flag)
                break;
            float xLabel = seg.row[indexs[size-1]].x;
            float yLabel = seg.row[indexs[size-1]].y;
            float zLabel = seg.row[indexs[size-1]].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            if(wholeGrid2segIDmap[gridKey].size()!=2)
            {
                tipBranchLength=0;
                break;
            }
            for(auto segIt=wholeGrid2segIDmap[gridKey].begin(); segIt!=wholeGrid2segIDmap[gridKey].end(); segIt++){
                if(segId != *segIt)
                {
                    segId = *segIt;
                    break;
                }
            }

            if(find(visitedSegIds.begin(),visitedSegIds.end(),segId)==visitedSegIds.end())
                visitedSegIds.push_back(segId);
            else
            {
                tipBranchLength=0;
                break;
            }
            seg = inputSegList.seg[segId];
            float xLabel2 = seg.row[0].x;
            float yLabel2 = seg.row[0].y;
            float zLabel2 = seg.row[0].z;
            QString gridKeyQ2 = QString::number(xLabel2) + "_" + QString::number(yLabel2) + "_" + QString::number(zLabel2);
            string gridKey2 = gridKeyQ2.toStdString();
            if(gridKey2!=gridKey)
                isReverse=true;
            else
                isReverse=false;
        }

        if(tipBranchLength>=dist_thresh){
            NeuronSWC s;
            stringToXYZ(*it,s.x,s.y,s.z);
            s.type = 10;
            if(s.x>33&&s.x+33<maxRes.x&&s.y>33&&s.y+33<maxRes.y&&s.z>33&&s.z+33<maxRes.z)
            {
                QString qKey = QString::number(s.x) + "_" + QString::number(s.y) + "_" + QString::number(s.z);
                string key = qKey.toStdString();
                outputSpecialPoints.push_back(s);
            }

            CellAPO marker;
            marker.x = s.x;
            marker.y = s.y;
            marker.z = s.z;
            marker.color.r = 250;
            marker.color.g = 100;
            marker.color.b = 120;
            tipAllMarkers.push_back(marker);
        }
    }

    logOut << "tipDetection end\n";
    if(!removeFlag)
        return outputSpecialPoints;

    vector<V_NeuronSWC> delSegVec;
    for(auto it=outputSpecialPoints.begin(); it!=outputSpecialPoints.end(); it++){
        float xLabel = it->x;
        float yLabel = it->y;
        float zLabel = it->z;
        QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
        string gridKey = gridKeyQ.toStdString();

        size_t segId=*wholeGrid2segIDmap[gridKey].begin();

        V_NeuronSWC& seg = inputSegList.seg[segId];

        float xLabel0 = seg.row[0].x;
        float yLabel0 = seg.row[0].y;
        float zLabel0 = seg.row[0].z;
        QString gridKeyQ0 = QString::number(xLabel0) + "_" + QString::number(yLabel0) + "_" + QString::number(zLabel0);
        string gridKey0 = gridKeyQ0.toStdString();

        bool isReverse=false;
        if(wholeGrid2segIDmap[gridKey0].size()!=1)
        {
            isReverse=true;
        }

        int size=seg.row.size();
        vector<int> indexs(size);
        for(int m=0;m<size;m++)
            indexs[m]=m;
        if(isReverse)
            reverse(indexs.begin(),indexs.end());

//        if(size<=7){
//            auto segIt=findseg(inputSegList.seg.begin(),inputSegList.seg.end(),seg);
//            if(segIt!=inputSegList.seg.end())
//            {
//                delSegVec.push_back(seg);
//            }
//            continue;
//        }

        qDebug()<<"before: ";
        seg.printInfo();

        int count=1;
        while(count<=1){
            auto tmpIt = seg.row.end();
            if(!isReverse){
                tmpIt=seg.row.begin();
            }
            else{
                tmpIt=seg.row.end()-1;
            }
            seg.row.erase(tmpIt);
            count++;
        }

        int nodeNo = 1;
        for (vector<V_NeuronSWC_unit>::iterator it_unit = seg.row.begin();
             it_unit != seg.row.end(); it_unit++)
        {
            it_unit->data[0] = nodeNo;
            it_unit->data[6] = nodeNo + 1;
            ++nodeNo;
        }
        (seg.row.end()-1)->data[6]=-1;
        qDebug()<<"after: ";
        seg.printInfo();
    }

//    for(auto delIt=delSegVec.begin(); delIt!=delSegVec.end(); delIt++){
//        auto segIt=findseg(inputSegList.seg.begin(),inputSegList.seg.end(),*delIt);
//        if(segIt!=inputSegList.seg.end())
//        {
//            inputSegList.seg.erase(segIt);
//        }
//    }

    return outputSpecialPoints;
}

QJsonArray CollDetection::crossingDetection(){
    logOut << "begin crossingDetection...";
    if(maxRes.x == 0){
        return QJsonArray();
    }
    QString swcFileName=inFilename + ".ano.eswc";

    std::filesystem::path swcPath(tmpInFile.toStdString());

    std::vector<NeuronUnit> neurons;
    std::filesystem::path fullPath = swcPath;
    ESwc swc(fullPath.string());
    neurons = swc.getNeuron();

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<util::Node> nodes;
    std::vector<int> rootNodeIds;

    std::vector<CrossingDetect::KeyPointsType> keypointsList;
    std::vector<CrossingDetect::BranchesType> selectedBranchesList;
    std::vector<util::ImageResolutionInfo> imageResolutionInfoList;

    int curNode = 0;
    for (auto &neuron: neurons) {
        nodes.push_back({.n=neuron.n, .parent=neuron.parent, .x=neuron.x, .y=neuron.y, .z=neuron.z});
        if (neuron.parent == -1) {
            rootNodeIds.push_back(neuron.n);
        }
        curNode++;

        if ((rootNodeIds.size() == 2) || (rootNodeIds.size() == 1 && neurons.size() == curNode)) {

            util::Node lastNode{};
            int lastRootNodeId = -1;
            bool multiRoot = false;
            if (rootNodeIds.size() == 2) {
                multiRoot = true;
            }

            if (multiRoot) {
                lastNode = nodes.back();
                lastRootNodeId = rootNodeIds.back();
                nodes.pop_back();
                rootNodeIds.pop_back();
            }

            CrossingDetect instance;
            instance.initializeNodeData(nodes, rootNodeIds);
            instance.generateBranches();
            instance.selectBranches();
            instance.generateNearestKeyPoint();
            instance.removeFromKeyPoint();

            auto &keyPoints = instance.getKeyPoint();
            auto &selectedBranch = instance.getSelectedBranch();

            auto edgeThreshold = ConfigManager::getInstance().edgePointIgnoreThreshold;
            for (auto &point: keyPoints) {
                if (maxRes.x > edgeThreshold * 2
                    && maxRes.y > edgeThreshold * 2
                    && maxRes.z > edgeThreshold * 2) {
                    auto isXConstrictionOK = point.first.first.first.x >= edgeThreshold &&
                                             point.first.first.first.x <= maxRes.x - edgeThreshold;
                    auto isYConstrictionOK = point.first.first.first.y >= edgeThreshold &&
                                             point.first.first.first.y <= maxRes.y - edgeThreshold;
                    auto isZConstrictionOK = point.first.first.first.z >= edgeThreshold &&
                                             point.first.first.first.z <= maxRes.z - edgeThreshold;

                    if (!(isXConstrictionOK && isYConstrictionOK && isZConstrictionOK)) {
                        point.second = false;
                        continue;
                    }
                }
            }

            keypointsList.push_back(keyPoints);
            qDebug()<<"keyPoints.size(): " << keyPoints.size();
            selectedBranchesList.push_back(selectedBranch);

            // clear data before processing more root nodes
            nodes.clear();
            rootNodeIds.clear();

            if (multiRoot) {
                nodes.push_back(lastNode);
                rootNodeIds.push_back(lastRootNodeId);
            }
        }
    }

    ResultWriter writer("", "");
    QJsonArray infos = writer.getData(keypointsList, selectedBranchesList, crossingAllMarkers);
    qDebug()<<infos;
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = endTime - startTime;

    logOut << "\n";
    logOut << "Using time:"<<std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()<<"ms\n";
    logOut << "--------------------------\n";
    logOut << "crossingDetection end\n";
    return infos;
}

void CollDetection::handleMulFurcation(vector<NeuronSWC>& outputSpecialPoints, double dist_thre){
    logOut << "handleMulFurcation begin...\n";
    vector<CellAPO> errorMarkers;
    for(int i=0;i<outputSpecialPoints.size();i++){
        RGB8 color = getColorFromType(outputSpecialPoints[i].type);
//        if(isSomaExists)
//        {
//            if(distance(outputSpecialPoints[i].x, somaCoordinate.x, outputSpecialPoints[i].y, somaCoordinate.y,
//                         outputSpecialPoints[i].z, somaCoordinate.z) > dist_thre)
//            {
//                CellAPO marker;
//                marker.name="";
//                marker.comment="quality_control";
//                marker.orderinfo="";
//                marker.color.r=188;
//                marker.color.g=94;
//                marker.color.b=37;
//                marker.x=outputSpecialPoints[i].x;
//                marker.y=outputSpecialPoints[i].y;
//                marker.z=outputSpecialPoints[i].z;
//                errorMarkers.push_back(marker);
//            }
//        }
//        else{
        CellAPO marker;
        marker.name="";
        marker.comment="quality_control";
        marker.orderinfo="";
        marker.color.r=188;
        marker.color.g=94;
        marker.color.b=37;
        marker.x=outputSpecialPoints[i].x;
        marker.y=outputSpecialPoints[i].y;
        marker.z=outputSpecialPoints[i].z;
        errorMarkers.push_back(marker);
//        }
    }
    addMarkers(errorMarkers);
    mulFurcationMarkers = errorMarkers;
    logOut << "handleMulFurcation end\n";
}

void CollDetection::handleLoop(vector<NeuronSWC>& outputSpecialPoints){
    logOut << "handleLoop begin...\n";
    vector<CellAPO> errorMarkers;
    for(int i=0;i<outputSpecialPoints.size();i++){
        CellAPO marker;
        marker.name="";
        marker.comment="quality_control";
        marker.orderinfo="";
        marker.color.r=255;
        marker.color.g=255;
        marker.color.b=255;
        marker.x=outputSpecialPoints[i].x;
        marker.y=outputSpecialPoints[i].y;
        marker.z=outputSpecialPoints[i].z;
        errorMarkers.push_back(marker);
    }
    addMarkers(errorMarkers);
    loopMarkers = errorMarkers;
    logOut << "handleLoop end\n";
}

void CollDetection::handleNearBifurcation(vector<NeuronSWC>& bifurPoints){
    logOut << "handleNearBifurcation begin...\n";
    vector<CellAPO> errorMarkers;
    for(int i=0;i<bifurPoints.size();i++){
        CellAPO marker;
        marker.name="";
        marker.comment="quality_control";
        marker.orderinfo="";
        marker.color.r=220;
        marker.color.g=200;
        marker.color.b=0;
        marker.x=bifurPoints[i].x;
        marker.y=bifurPoints[i].y;
        marker.z=bifurPoints[i].z;
        errorMarkers.push_back(marker);
    }
    addMarkers(errorMarkers);
    nearBifurcationMarkers = errorMarkers;
    logOut << "handleNearBifurcation end\n";
}

void CollDetection::handleTip(vector<NeuronSWC>& tipPoints){
    logOut << "handleTip begin...\n";

    if(tipPoints.size()!=0){
        QHttpMultiPart *multiPart = new QHttpMultiPart(QHttpMultiPart::FormDataType);
        QHttpPart filePart;
        QString swcFileName=inFilename + ".ano.eswc";
        QString fileSavePath=tmpInFile;

        // 创建一个QFile对象，用于读取要上传的文件
        QFile *file = new QFile(fileSavePath);
        if(!file->open(QIODevice::Text|QIODevice::ReadWrite)){
            qDebug() << "cannot open file in handleTip";
        }
        file->setPermissions(QFileDevice::ReadOwner|QFileDevice::WriteOwner);

        filePart.setHeader(QNetworkRequest::ContentTypeHeader, QVariant("text/plain"));
        filePart.setHeader(QNetworkRequest::ContentDispositionHeader, QVariant("form-data; name=\"swcFile\"; filename=\""+ swcFileName + "\"")); // file为后端定义的key，filename即为excel文件名
        filePart.setBodyDevice(file);
        file->setParent(multiPart); // 文件将由multiPart对象进行管理

        multiPart->append(filePart);
        QNetworkRequest fileRequest(QUrl(SuperUserHostAddress+"/detect/file/for-missing"));
        // 发送HTTP POST请求
        qDebug()<<"123";
        QNetworkReply* fileReply = accessManager->post(fileRequest, multiPart);
        multiPart->setParent(fileReply); // reply对象将负责删除multiPart对象
        QEventLoop tmpeventLoop;
        connect(fileReply, &QNetworkReply::finished, &tmpeventLoop, &QEventLoop::quit);
        tmpeventLoop.exec(QEventLoop::ExcludeUserInputEvents);

        qDebug()<<"123";
        if(fileReply->error())
        {
            logOut << "SENDFILEERROR!\n";
            logOut << fileReply->errorString() << "\n";
        }
        int fileResCode=fileReply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
        logOut << "sendFile "<< fileResCode << "\n";
        file->close();
        fileReply->deleteLater();

        QJsonObject json;
        QString obj=image;
        json.insert("obj",obj);
//        json.insert("res", myServer->RES);
        json.insert("swcNameWithNoSuffix", inFilename);
        QJsonArray coorList;
        for(int i=0; i<tipPoints.size();i++){
            QJsonObject coor;
            coor.insert("x", tipPoints[i].x);
            coor.insert("y", tipPoints[i].y);
            coor.insert("z", tipPoints[i].z);
            coorList.append(coor);
        }
        json.insert("coors",coorList);

        QJsonDocument document;
        document.setObject(json);
        QString str=QString(document.toJson());
        QByteArray byteArray=str.toUtf8();

        // 创建一个QNetworkRequest对象，设置URL和请求方法
        QNetworkRequest request(QUrl(SuperUserHostAddress+"/detect/missing"));
        request.setHeader(QNetworkRequest::ContentTypeHeader,"application/json");
        //        request.setRawHeader("Content-Type", "multipart/form-data; boundary=" + multiPart->boundary());

        // 发送HTTP POST请求
        QNetworkReply* reply = accessManager->post(request, byteArray);

        QEventLoop eventLoop;
        connect(reply, &QNetworkReply::finished, &eventLoop, &QEventLoop::quit);
        eventLoop.exec(QEventLoop::ExcludeUserInputEvents);

        if(reply->error())
        {
            logOut << "ERROR!\n";
            logOut << reply->errorString() << "\n";
        }
        int code=reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
        logOut<<"handleTip "<<code<<"\n";
        QByteArray responseData = reply->readAll();
        vector<NeuronSWC> markPoints;
        if(code==200)
        {
            //解析json
            QJsonParseError json_error;
            QJsonDocument doucment = QJsonDocument::fromJson(responseData, &json_error);
            if (json_error.error == QJsonParseError::NoError) {
                if (doucment.isObject()) {
                    const QJsonObject obj = doucment.object();
                    QString objCode;
                    QString objMsg;
                    if (obj.contains("code")) {
                        QJsonValue value = obj.value("code");
                        if (value.isString()) {
                            objCode= value.toString();
                            logOut << "code : " << objCode<<"\n";
                        }
                    }
                    if (obj.contains("msg")) {
                        QJsonValue value = obj.value("msg");
                        if (value.isString()) {
                            QString objMsg = value.toString();
                            logOut << "msg : " << objMsg<<"\n";
                        }
                    }
                    if (obj.contains("data")&&objCode=="200") {
                        QJsonValue value = obj.value("data");
                        if (value.isArray()) {  // Version 的 value 是数组
                            QJsonArray array = value.toArray();
                            int nSize = array.size();
                            for (int i = 0; i < nSize; ++i) {
                                QJsonValue mapValue = array.at(i);
                                if (mapValue.isObject()) {
                                    QJsonObject info = mapValue.toObject();
                                    float x,y,z;
                                    int y_pred;
                                    if (info.contains("coors")) {
                                        QJsonValue listValue = info.value("coors");
                                        if (listValue.isArray()) {
                                            QJsonArray listArray = listValue.toArray();
                                            QJsonValue xValue = listArray.at(0);
                                            QJsonValue yValue = listArray.at(1);
                                            QJsonValue zValue = listArray.at(2);
                                            x=xValue.toDouble();
                                            y=yValue.toDouble();
                                            z=zValue.toDouble();
                                            logOut << "x: " << x/2 << " y: " << y/2 << " z: " << z/2 <<" ";
                                        }
                                    }
                                    if (info.contains("y_pred")) {
                                        QJsonValue predValue = info.value("y_pred");
                                        y_pred = predValue.toInt();
                                        if(y_pred == 1)
                                            logOut << "missing\n";
                                        else if(y_pred == 0){
                                            logOut << "done\n";
                                        }
                                    }
                                    if(y_pred==1){
                                        NeuronSWC s;
                                        s.x=x;
                                        s.y=y;
                                        s.z=z;
                                        s.type=10;
                                        markPoints.push_back(s);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            logOut<<"handle tip error!\n";
        }

        vector<CellAPO> errorMarkers;
        for(int i=0;i<markPoints.size();i++){
            CellAPO marker;
            marker.name="";
            marker.comment="quality_control";
            marker.orderinfo="";
            marker.color.r=250;
            marker.color.g=100;
            marker.color.b=120;
            marker.x=markPoints[i].x;
            marker.y=markPoints[i].y;
            marker.z=markPoints[i].z;
            errorMarkers.push_back(marker);
        }
        addMarkers(errorMarkers);
        tipUndoneMarkers = errorMarkers;

        //清理资源
        reply->deleteLater();
    }
    logOut << "handleTip end\n";
}

void CollDetection::handleCrossing(QJsonArray& infos){
    logOut << "begin handleCrossing...\n";
    QString swcFileName=inFilename+".ano.eswc";
    QString fileSavePath=tmpInFile;
    // 创建一个QFile对象，用于读取要上传的文件
    QFile *file = new QFile(fileSavePath);
    file->open(QIODevice::Text|QIODevice::ReadWrite);

    if(!infos.isEmpty()){
        QHttpMultiPart *multiPart = new QHttpMultiPart(QHttpMultiPart::FormDataType);
        QHttpPart filePart;

        filePart.setHeader(QNetworkRequest::ContentTypeHeader, QVariant("text/plain"));
        filePart.setHeader(QNetworkRequest::ContentDispositionHeader, QVariant("form-data; name=\"swcFile\"; filename=\""+ swcFileName + "\"")); // file为后端定义的key，filename即为excel文件名
        filePart.setBodyDevice(file);

        multiPart->append(filePart);
        QNetworkRequest fileRequest(QUrl(SuperUserHostAddress+"/detect/file/for-crossing"));
        // 发送HTTP POST请求
        QNetworkReply* fileReply = accessManager->post(fileRequest, multiPart);
        multiPart->setParent(fileReply); // reply对象将负责删除multiPart对象
        QEventLoop tmpeventLoop;
        connect(fileReply, &QNetworkReply::finished, &tmpeventLoop, &QEventLoop::quit);
        tmpeventLoop.exec(QEventLoop::ExcludeUserInputEvents);

        if(fileReply->error())
        {
            logOut << "SENDFILEERROR!\n";
            logOut << fileReply->errorString() << "\n";
        }
        int fileResCode=fileReply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
        logOut<<"sendFile "<<fileResCode<<"\n";

        QJsonObject json;
        QString obj=image;
        json.insert("obj",obj);
//        json.insert("res", myServer->RES);
        json.insert("swcNameWithNoSuffix", inFilename);
        json.insert("infos",infos);

        QJsonDocument document;
        document.setObject(json);
        QString str=QString(document.toJson());
        QByteArray byteArray=str.toUtf8();

        // 创建一个QNetworkRequest对象，设置URL和请求方法
        QNetworkRequest request(QUrl(SuperUserHostAddress+"/detect/crossing"));
        request.setHeader(QNetworkRequest::ContentTypeHeader,"application/json");
        //        request.setRawHeader("Content-Type", "multipart/form-data; boundary=" + multiPart->boundary());

        // 发送HTTP POST请求
        QNetworkReply* reply = accessManager->post(request, byteArray);

        QEventLoop eventLoop;
        connect(reply, &QNetworkReply::finished, &eventLoop, &QEventLoop::quit);
        eventLoop.exec(QEventLoop::ExcludeUserInputEvents);

        if(reply->error())
        {
            logOut << "ERROR!\n";
            logOut << reply->errorString() <<"\n";
        }
        int code=reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
        logOut<<"handleCrossing "<<code<<"\n";
        QByteArray responseData = reply->readAll();
        vector<NeuronSWC> markPoints;
        if(code==200)
        {
            //解析json
            QJsonParseError json_error;
            QJsonDocument doucment = QJsonDocument::fromJson(responseData, &json_error);
            if (json_error.error == QJsonParseError::NoError) {
                if (doucment.isObject()) {
                    const QJsonObject obj = doucment.object();
                    QString objCode;
                    QString objMsg;
                    if (obj.contains("code")) {
                        QJsonValue value = obj.value("code");
                        if (value.isString()) {
                            objCode= value.toString();
                            logOut << "code : " << objCode <<"\n";
                        }
                    }
                    if (obj.contains("msg")) {
                        QJsonValue value = obj.value("msg");
                        if (value.isString()) {
                            QString objMsg = value.toString();
                            logOut << "msg : " << objMsg <<"\n";
                        }
                    }
                    if (obj.contains("data")&&objCode=="200") {
                        QJsonValue value = obj.value("data");
                        if (value.isArray()) {  // Version 的 value 是数组
                            QJsonArray array = value.toArray();
                            int nSize = array.size();
                            for (int i = 0; i < nSize; ++i) {
                                QJsonValue mapValue = array.at(i);
                                if (mapValue.isObject()) {
                                    QJsonObject info = mapValue.toObject();
                                    float x,y,z;
                                    int y_pred;
                                    if (info.contains("coors")) {
                                        QJsonValue listValue = info.value("coors");
                                        if (listValue.isArray()) {
                                            QJsonArray listArray = listValue.toArray();
                                            QJsonValue xValue = listArray.at(0);
                                            QJsonValue yValue = listArray.at(1);
                                            QJsonValue zValue = listArray.at(2);
                                            x=xValue.toDouble();
                                            y=yValue.toDouble();
                                            z=zValue.toDouble();
                                            logOut << "x: " << x/2 << " y: " << y/2 <<" z: " << z/2 << " ";
                                        }
                                    }
                                    if (info.contains("y_pred")) {
                                        QJsonValue predValue = info.value("y_pred");
                                        y_pred = predValue.toInt();
                                        if(y_pred == 0)
                                            logOut << "crossing error\n";
                                        else if(y_pred == 1){
                                            logOut << "crossing right\n";
                                        }

                                    }
                                    if(y_pred==0){
                                        NeuronSWC s;
                                        s.x=x;
                                        s.y=y;
                                        s.z=z;
                                        s.type=18;
                                        markPoints.push_back(s);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            logOut<<"handle crossing error!\n";
        }

        vector<CellAPO> errorMarkers;
        for(int i=0;i<markPoints.size();i++){
            CellAPO marker;
            marker.name="";
            marker.comment="quality_control";
            marker.orderinfo="";
            marker.color.r=168;
            marker.color.g=128;
            marker.color.b=255;
            marker.x=markPoints[i].x;
            marker.y=markPoints[i].y;
            marker.z=markPoints[i].z;
            errorMarkers.push_back(marker);
        }
        addMarkers(errorMarkers);
        crossingMarkers = errorMarkers;

        //清理资源
        reply->deleteLater();
    }

    //清理资源
    file->close();
    file->deleteLater();
    logOut << "handleCrossing end\n";
}

void CollDetection::handleOverlapSegs(set<string>& outputSpecialPoints){
    logOut << "handleOverlapSegs begin...\n";
    vector<CellAPO> errorMarkers;
    for(auto it=outputSpecialPoints.begin(); it!=outputSpecialPoints.end(); it++){
        NeuronSWC s;
        stringToXYZ(*it, s.x, s.y, s.z);
        CellAPO marker;
        marker.name="";
        marker.comment="quality_control";
        marker.orderinfo="";
        marker.color.r=0;
        marker.color.g=20;
        marker.color.b=200;
        marker.x=s.x;
        marker.y=s.y;
        marker.z=s.z;
        errorMarkers.push_back(marker);
    }
    addMarkers(errorMarkers);
    overlapSegsMarkers = errorMarkers;
    logOut << "handleOverlapSegs end\n";
}

void CollDetection::sortSWC(QString fileOpenName, QString fileSaveName, double thres, V3DLONG rootid){
    logOut << "begin sortSWC...\n";
    QList<NeuronSWC> neuron, result;
    if (fileOpenName.endsWith(".swc") || fileOpenName.endsWith(".SWC") || fileOpenName.endsWith(".eswc") || fileOpenName.endsWith(".ESWC"))
        neuron = readSWC_file(fileOpenName).listNeuron;
    if (!SortSWC(neuron, result , rootid, thres))
    {
        logOut<<"Error in sorting swc"<<"\n";
    }
    if (!export_list2file(result, fileSaveName, fileOpenName))
    {
        logOut<<"Error in writing swc to file"<<"\n";
    }
    logOut << "sortSWC end\n";
}

bool CollDetection::sortSWCAndDetectLoop(QString fileOpenName, QString fileSaveName, V3DLONG rootid){
    logOut << "begin sortSWCAndDetectLoop...\n";
    QList<NeuronSWC> neuron, result;
    V_NeuronSWC_list segments;
    if (fileOpenName.endsWith(".swc") || fileOpenName.endsWith(".SWC") || fileOpenName.endsWith(".eswc") || fileOpenName.endsWith(".ESWC")){
        auto nt = readSWC_file(fileOpenName);
        neuron = nt.listNeuron;
        segments = NeuronTree__2__V_NeuronSWC_list(nt);
    }

    if (!SortSWCAndDetectLoop(neuron, segments, result, loopMarkers, rootid))
    {
        int count = 0;
        if(loopMarkers.size()!=0){
            logOut<<"Error in sortSWCAndDetectLoop"<<"\n";
            logOut<<"have detected loops"<<"\n";
            addMarkers(loopMarkers);
        }
        logOut << "sortAndDetectLoop end\n";
        return false;
    }
    if (!export_list2file(result, fileSaveName, fileOpenName))
    {
        logOut<<"Error in writing swc to file"<<"\n";
        logOut << "sortAndDetectLoop end\n";
        return false;
    }
    logOut << "sortAndDetectLoop end\n";
    return true;
}

void CollDetection::setSWCRadius(QString filePath, int r){
    logOut << "begin setSWCRadius...\n";
    if (filePath.endsWith(".swc") || filePath.endsWith(".SWC") || filePath.endsWith(".eswc") || filePath.endsWith(".ESWC"))
    {
        QFile qf(filePath);
        QString arryRead;
        if(!qf.open(QIODevice::ReadOnly|QIODevice::Text)){
            return;
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
            qDebug()<<"swcfile cannot be opened!";
            return;
        }
        QTextStream streamWrite(&qf);
        for(int i=0;i<arryListWrite.size()-1;i++){      //这里到arryListWrite.size()-1是因为arryListWrite数组按照\n分 段时，最后一行尾部有个\n，所以数组最后一个值为空，需要将它去掉
            if(arryListWrite.at(i).contains("#")){
                streamWrite<<arryListWrite.at(i)<<"\n";
            }else{
                QString contentWrite= arryListWrite.at(i);
                QStringList swcInfo=contentWrite.split(' ',Qt::SkipEmptyParts);
                swcInfo[5]=QString::number(r);
                contentWrite=swcInfo.join(' ');
                streamWrite<<contentWrite<<"\n";
            }
        }
        qf.close();
    }
    logOut << "setSWCRadius end\n";
}

QStringList CollDetection::getSWCSpecNInfo(QString filePath, int val){
    logOut << "begin getSWCSpecNCount...\n";
    QStringList result;
    if (filePath.endsWith(".swc") || filePath.endsWith(".SWC") || filePath.endsWith(".eswc") || filePath.endsWith(".ESWC"))
    {
        QFile qf(filePath);
        QString arryRead;
        if(!qf.open(QIODevice::ReadOnly|QIODevice::Text)){
            return result;
        }
        arryRead=qf.readAll();

        qf.close();

        QStringList arryListWrite= arryRead.split("\n");
        //        for(int i=0;i<arryListWrite.size();i++){
        //            qDebug()<<arryListWrite.at(i);
        //        }

        for(int i=0;i<arryListWrite.size()-1;i++){      //这里到arryListWrite.size()-1是因为arryListWrite数组按照\n分 段时，最后一行尾部有个\n，所以数组最后一个值为空，需要将它去掉
            if(arryListWrite.at(i).contains("#")){
                continue;
            }else{
                QString contentWrite= arryListWrite.at(i);
                QStringList swcInfo=contentWrite.split(' ',Qt::SkipEmptyParts);
                if(swcInfo[6] == QString::number(val)){
                    result.push_back(contentWrite);
                }
            }
        }
    }
    return result;
    logOut << "getSWCSpecNCount end\n";
}

void CollDetection::getImageRES(){
    logOut<<"begin getImageRES...\n";

    QJsonObject json;
    QJsonObject userMap;
    userMap.insert("name", "zackzhy");
    userMap.insert("passwd", "123456");
    QJsonObject imageMap;
    imageMap.insert("name", image);
    imageMap.insert("detail", "");
    json.insert("user", userMap);
    json.insert("Image", imageMap);

    QJsonDocument document;
    document.setObject(json);
    QString str=QString(document.toJson());
    QByteArray byteArray=str.toUtf8();

    // 创建一个QNetworkRequest对象，设置URL和请求方法
    QNetworkRequest request(QUrl(BrainTellHostAddress+"/image/getimagelist"));
    request.setHeader(QNetworkRequest::ContentTypeHeader,"application/json");

    // 发送HTTP POST请求
    QNetworkReply* reply = accessManager->post(request, byteArray);

    QEventLoop eventLoop;
    connect(reply, &QNetworkReply::finished, &eventLoop, &QEventLoop::quit);
    eventLoop.exec(QEventLoop::ExcludeUserInputEvents);

    if(reply->error())
    {
        logOut << "ERROR!\n";
        logOut << reply->errorString() << "\n";
    }
    int code=reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
    logOut << "getImageRes "<< code << "\n";
    QByteArray responseData = reply->readAll();

    QString maxResString;
    QString subMaxResString;
    if(code==200){
        //解析json
        QJsonParseError json_error;
        QJsonDocument doucment = QJsonDocument::fromJson(responseData, &json_error);
        if (json_error.error == QJsonParseError::NoError) {
            if (doucment.isArray()) {
                const QJsonArray array = doucment.array();
//                qDebug() << array;
                int nSize = array.size();
                for (int i = 0; i < nSize; ++i) {
                    QJsonValue mapValue = array.at(i);
                    if (mapValue.isObject()) {
                        QJsonObject info = mapValue.toObject();
                        QString name="";
                        QString details;
                        if(info.contains("name")){
                            QJsonValue nameVal=info.value("name");
                            if(nameVal.isString()){
                                name=nameVal.toString();
                            }
                        }

                        if(name!=image)
                            continue;
                        if(info.contains("detail")){
                            QJsonValue detailsVal=info.value("detail");
                            if(detailsVal.isString()){
                                details=detailsVal.toString();
                                QJsonDocument detailsDoc=QJsonDocument::fromJson(details.toUtf8());
                                if(detailsDoc.isArray()){
                                    QJsonArray detailsArray=detailsDoc.array();
                                    maxResString=detailsArray.at(0).toString();
                                    subMaxResString=detailsArray.at(1).toString();
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    logOut << maxResString << " " << subMaxResString << "\n";

    QRegularExpression regex("^RES\\((\\d+)x(\\d+)x(\\d+)\\)");
    QRegularExpressionMatchIterator matches = regex.globalMatch(maxResString);

    while (matches.hasNext()) {
        QRegularExpressionMatch match = matches.next();
        QString matchedText = match.captured(0);
        QString matchedY = match.captured(1);
        QString matchedX = match.captured(2);
        QString matchedZ = match.captured(3);
        maxRes.x=matchedX.toFloat();
        maxRes.y=matchedY.toFloat();
        maxRes.z=matchedZ.toFloat();
    }

    matches=regex.globalMatch(subMaxResString);
    while (matches.hasNext()) {
        QRegularExpressionMatch match = matches.next();
        QString matchedText = match.captured(0);
        QString matchedY = match.captured(1);
        QString matchedX = match.captured(2);
        QString matchedZ = match.captured(3);
        subMaxRes.x=matchedX.toFloat();
        subMaxRes.y=matchedY.toFloat();
        subMaxRes.z=matchedZ.toFloat();
    }

    //清理资源
    reply->deleteLater();
    logOut<<"getImageRes end\n";
}

//void CollDetection::getImageMaxRES(){
//    // 定义正则表达式
//    QRegularExpression regex("RES\\((\\d+)x(\\d+)x(\\d+)\\)");

//    // 进行匹配
//    QRegularExpressionMatch match = regex.match(myServer->RES);

//    // 检查匹配是否成功
//    if (match.hasMatch()) {
//        // 提取宽度、高度和深度
//        QString yStr = match.captured(1);
//        QString xStr = match.captured(2);
//        QString zStr = match.captured(3);

//        maxRes.x=xStr.toInt();
//        maxRes.y=yStr.toInt();
//        maxRes.z=zStr.toInt();
//    } else {
//        qDebug() << "No match found";
//    }
//}

void CollDetection::removeErrorSegs(V_NeuronSWC_list& segments){
    logOut << "begin removeErrorSegs...\n";

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC seg = segments.seg[i];
        if(seg.row.size()==1){
            segments.seg[i].to_be_deleted = true;
            errorSegsNum++;
            for (V3DLONG p=0;p<seg.row.size();p++)
            {
                V_NeuronSWC_unit v = seg.row.at(p);
                logOut <<"row ["<< p <<"] : "<<v.data[0]<<" "<<v.data[1]<<" "<<v.data[2]<<" "<<v.data[3]<<" "<<v.data[4]<<" "<<v.data[5]<<" "<<v.data[6] << "\n";
            }
            continue;
        }
        set<string> coors;
        for(size_t j=0; j<seg.row.size(); j++){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            coors.insert(gridKey);
        }

        if(coors.size() < seg.row.size())
        {
            segments.seg[i].to_be_deleted = true;
            for (V3DLONG p=0;p<seg.row.size();p++)
            {
                V_NeuronSWC_unit v = seg.row.at(p);
                logOut <<"row ["<< p <<"] : "<<v.data[0]<<" "<<v.data[1]<<" "<<v.data[2]<<" "<<v.data[3]<<" "<<v.data[4]<<" "<<v.data[5]<<" "<<v.data[6] << "\n";
            }
            errorSegsNum++;
            continue;
        }
    }

    auto iter = segments.seg.begin();
    while (iter != segments.seg.end())
        if (iter->to_be_deleted){;
            iter = segments.seg.erase(iter);
        }
        else
            ++iter;

    logOut << "removeErrorSegs end\n";
}

void CollDetection::tuneErrorSegs(V_NeuronSWC_list& segments){
    logOut << "begin tuneErrorSegs...\n";

    for(size_t i=0; i<segments.seg.size(); ++i){
        V_NeuronSWC& seg = segments.seg[i];
        if(seg.row.size()==1){
            segments.seg[i].to_be_deleted = true;
            errorSegsNum++;
            for (V3DLONG p=0;p<seg.row.size();p++)
            {
                V_NeuronSWC_unit v = seg.row.at(p);
                logOut <<"row ["<< p <<"] : "<<v.data[0]<<" "<<v.data[1]<<" "<<v.data[2]<<" "<<v.data[3]<<" "<<v.data[4]<<" "<<v.data[5]<<" "<<v.data[6] << "\n";
            }
            continue;
        }
        set<string> coors;
        for(size_t j=0; j<seg.row.size(); j++){
            float xLabel = seg.row[j].x;
            float yLabel = seg.row[j].y;
            float zLabel = seg.row[j].z;
            QString gridKeyQ = QString::number(xLabel) + "_" + QString::number(yLabel) + "_" + QString::number(zLabel);
            string gridKey = gridKeyQ.toStdString();
            coors.insert(gridKey);
        }

        if(coors.size() < seg.row.size())
        {
            for (V3DLONG p=0;p<seg.row.size();p++)
            {
                V_NeuronSWC_unit v = seg.row.at(p);
                logOut <<"row ["<< p <<"] : "<<v.data[0]<<" "<<v.data[1]<<" "<<v.data[2]<<" "<<v.data[3]<<" "<<v.data[4]<<" "<<v.data[5]<<" "<<v.data[6] << "\n";
            }
            errorSegsNum++;
            for (auto it = seg.row.begin(); it!=seg.row.end() - 1 && it!=seg.row.end();)
            {
                V_NeuronSWC_unit v1 = *it;
                V_NeuronSWC_unit v2 = *(it+1);
                if(!(fabs(v1.x - v2.x) < 1e-5) || !(fabs(v1.y - v2.y) < 1e-5) || !(fabs(v1.z - v2.z) < 1e-5)){
                    it++;
                }
                else{
                    it = seg.row.erase(it);
                }
            }

            if(seg.row.size()==1){
                segments.seg[i].to_be_deleted = true;
                continue;
            }

            int count = 1;
            for (V3DLONG p=0;p<seg.row.size();p++)
            {
                V_NeuronSWC_unit& v = seg.row.at(p);
                v.n = count++;
                v.parent = count;
            }
            seg.row[seg.row.size() - 1].parent = -1;
        }
    }

    auto iter = segments.seg.begin();
    while (iter != segments.seg.end())
        if (iter->to_be_deleted){;
            iter = segments.seg.erase(iter);
        }
        else
            ++iter;

    logOut << "tuneErrorSegs end\n";
}

void CollDetection::detectOverlapSegs(V_NeuronSWC_list inputSegList, double dist_thres){
    logOut << "begin detectOverlapSegs...\n";
    set<size_t> overlapSegIds;
    //检测overlap线段
    for(size_t i=0; i<inputSegList.seg.size(); ++i){
        float xLabel1 = inputSegList.seg[i].row[0].x;
        float yLabel1 = inputSegList.seg[i].row[0].y;
        float zLabel1 = inputSegList.seg[i].row[0].z;
        float xLabel2 = inputSegList.seg[i].row[inputSegList.seg[i].row.size()-1].x;
        float yLabel2 = inputSegList.seg[i].row[inputSegList.seg[i].row.size()-1].y;
        float zLabel2 = inputSegList.seg[i].row[inputSegList.seg[i].row.size()-1].z;

        if(isSomaExists){
            if(distance(xLabel1, somaCoordinate.x, yLabel1, somaCoordinate.y, zLabel1, somaCoordinate.z) < dist_thres)
                continue;
            if(distance(xLabel2, somaCoordinate.x, yLabel2, somaCoordinate.y, zLabel2, somaCoordinate.z) < dist_thres)
                continue;
        }

        for(size_t j=i+1; j<inputSegList.seg.size(); j++){
            float xLabel3 = inputSegList.seg[j].row[0].x;
            float yLabel3 = inputSegList.seg[j].row[0].y;
            float zLabel3 = inputSegList.seg[j].row[0].z;
            float xLabel4 = inputSegList.seg[j].row[inputSegList.seg[j].row.size()-1].x;
            float yLabel4 = inputSegList.seg[j].row[inputSegList.seg[j].row.size()-1].y;
            float zLabel4 = inputSegList.seg[j].row[inputSegList.seg[j].row.size()-1].z;

            if(isSomaExists){
                if(distance(xLabel3, somaCoordinate.x, yLabel3, somaCoordinate.y, zLabel3, somaCoordinate.z) < dist_thres)
                    continue;
                if(distance(xLabel4, somaCoordinate.x, yLabel4, somaCoordinate.y, zLabel4, somaCoordinate.z) < dist_thres)
                    continue;
            }
            int result = isOverlapOfTwoSegs(logOut, inputSegList.seg[i], inputSegList.seg[j]);
            if(result == 1)
                overlapSegIds.insert(i);
            if(result == 2)
                overlapSegIds.insert(j);
        }
    }

    set<string> outputSpecialPoints;
    for(auto it=overlapSegIds.begin(); it!=overlapSegIds.end(); it++){
        V_NeuronSWC seg = inputSegList.seg[*it];
        float xLabel1 = seg.row[0].x;
        float yLabel1 = seg.row[0].y;
        float zLabel1 = seg.row[0].z;
        QString gridKeyQ1 = QString::number(xLabel1) + "_" + QString::number(yLabel1) + "_" + QString::number(zLabel1);
        string gridKey1 = gridKeyQ1.toStdString();
        float xLabel2 = seg.row[seg.row.size()-1].x;
        float yLabel2 = seg.row[seg.row.size()-1].y;
        float zLabel2 = seg.row[seg.row.size()-1].z;
        QString gridKeyQ2 = QString::number(xLabel2) + "_" + QString::number(yLabel2) + "_" + QString::number(zLabel2);
        string gridKey2 = gridKeyQ2.toStdString();
        outputSpecialPoints.insert(gridKey1);
        outputSpecialPoints.insert(gridKey2);
    }

    handleOverlapSegs(outputSpecialPoints);
    logOut << "detectOverlapSegs end\n";
}

void CollDetection::addMarkers(vector<CellAPO> cellApos){
    logOut<<"begin addMarkers...\n";

    for(auto it=cellApos.begin(); it!=cellApos.end(); it++){
        bool isSucess = true;
//        for(auto markerIt=markers.begin();markerIt!=markers.end(); ++markerIt)
//        {
//            if(abs(it->x-markerIt->x)<1&&abs(it->y-markerIt->y)<1&&abs(it->z-markerIt->z)<1)
//            {
//                logOut<<"the marker has already existed\n";
//                isSucess = false;
//                break;
//            }
//        }
        if(isSucess)
            markers.append(*it);
    }

    logOut<<"addMarkers end\n";
}

