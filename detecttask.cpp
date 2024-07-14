#include "detecttask.h"

DetectTask::DetectTask(){}

DetectTask::DetectTask(QString infilePath, QString fileName, QString logPath, QString apoPath, QString anoPath, QString swcOutputPath, QString tmpSwcPath, QString resultPath,
                       QString croppedApoMulfurPath, QString croppedApoBifurPath, QString croppedApoLoopPath, QString croppedApoMissingPath,
                       QString croppedApoCrossingPath, QString croppedApoOverlapPath, QString croppedApoDissoPath, QString croppedSwcMulfurPath,
                       QString croppedSwcBifurPath, QString croppedSwcLoopPath, QString croppedSwcMissingPath,
                       QString croppedSwcCrossingPath, QString croppedSwcOverlapPath, QString croppedSwcDissoPath)
{
    d_inFile=infilePath;
    d_fileName=fileName;
    d_logPath=logPath;
    d_apoPath=apoPath;
    d_anoPath=anoPath;
    d_resultPath=resultPath;
    d_swcOutputPath = swcOutputPath;
    d_tmpswcPath = tmpSwcPath;
    d_croppedApoMulfurPath = croppedApoMulfurPath;
    d_croppedApoBifurPath = croppedApoBifurPath;
    d_croppedApoLoopPath = croppedApoLoopPath;
    d_croppedApoMissingPath = croppedApoMissingPath;
    d_croppedApoCrossingPath = croppedApoCrossingPath;
    d_croppedApoOverlapPath = croppedApoOverlapPath;
    d_croppedApoDissoPath = croppedApoDissoPath;
    d_croppedSwcMulfurPath = croppedSwcMulfurPath;
    d_croppedSwcBifurPath = croppedSwcBifurPath;
    d_croppedSwcLoopPath = croppedSwcLoopPath;
    d_croppedSwcMissingPath = croppedSwcMissingPath;
    d_croppedSwcCrossingPath = croppedSwcCrossingPath;
    d_croppedSwcOverlapPath = croppedSwcOverlapPath;
    d_croppedSwcDissoPath = croppedSwcDissoPath;
}

QString DetectTask::getImageName(QString fileName){
    QString subStr1 = ".tif";
    QString subStr2 = ".v3d";
    int pos1 = fileName.indexOf(subStr1);
    int pos2 = fileName.indexOf(subStr2);
    QString imageName;
    if(pos1 != -1){
        imageName = fileName.left(pos1);
    }
    else if(pos2 != -1){
        imageName = fileName.left(pos2);
    }
    return imageName;
}

void DetectTask::run()
{
    // 在这里编写需要在后台线程中执行的代码
    qDebug() << "Task executed in thread: " << QThread::currentThread();
    CollDetection collDetection(d_inFile, d_fileName, d_logPath, d_apoPath, d_anoPath, d_swcOutputPath, d_tmpswcPath, d_resultPath,
                                d_croppedApoMulfurPath,
                                d_croppedApoBifurPath,
                                d_croppedApoLoopPath,
                                d_croppedApoMissingPath,
                                d_croppedApoCrossingPath,
                                d_croppedApoOverlapPath,
                                d_croppedApoDissoPath,
                                d_croppedSwcMulfurPath,
                                d_croppedSwcBifurPath,
                                d_croppedSwcLoopPath,
                                d_croppedSwcMissingPath,
                                d_croppedSwcCrossingPath,
                                d_croppedSwcOverlapPath,
                                d_croppedSwcDissoPath);
//    QStringList parts1 = d_fileName.split("_");
//    QStringList parts2 = d_fileName.split(".");
//    if(d_fileName.startsWith("brainID")){
//        collDetection.image=parts1[1];
//    }
//    else{
//        collDetection.image=parts2[0];
//    }

//    collDetection.image = getImageName(d_fileName);
//    colldetection->image=parts[1];

    QStringList parts = d_fileName.split("_");
    if(parts[1].size() == 1){
        collDetection.image = parts[0] + "_" + parts[1];
    }
    else{
        collDetection.image = parts[0];
    }
    qDebug()<<collDetection.image;

    //路径、不带后缀的文件名
    collDetection.detectAll();
    collDetection.generateResult();
//    collDetection.getApoAndCroppedSwc();
}


