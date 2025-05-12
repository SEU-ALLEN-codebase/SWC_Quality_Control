#include <QCoreApplication>
#include <stdio.h>
#include "utils.h"
#include <signal.h>
#include <unistd.h>
#include "detecttask.h"
#include <filesystem>

//QString swcInputDirPath = R"(/home/seu/Desktop/quality_control_swc_input)";

int main(int argc, char *argv[])
{
    //控制台应用程序
    QCoreApplication a(argc, argv);

    QString swcInputDirPath = argv[1];
//    bool needSort = argv[2];
    bool needSort = false;

    QString outputDirPath = swcInputDirPath + "_QC_Result";
    QString logDirPath = outputDirPath + "/" + "log";
    QString resultDirPath = outputDirPath + "/" + "result";
    QString resultSumCSVDirPath = outputDirPath + "/" + "result_sum";
    QString resultSumCSVAllPath = outputDirPath + "/" + "all.csv";
    QString errorResultSumCSVAllPath = outputDirPath + "/" + "error_all.csv";
    QString swcOutputDirPath = outputDirPath + "/" + "output";
    QString swcOutputAllDirPath = swcOutputDirPath + "/" + "all";
    QString swcOutputDefinedDirPath = swcOutputDirPath + "/" + "soma_defined_swc";

    QString apoDirPath = swcOutputDirPath;
    QString anoDirPath = swcOutputDirPath;

    CollDetection::swcInputDirPath = swcInputDirPath;
    CollDetection::needSort = needSort;
    QThreadPool threadPool;
    threadPool.setMaxThreadCount(QThread::idealThreadCount()); // 设置线程池的最大线程数

    QStringList swcFilters;
    swcFilters << "*.swc" << "*.eswc";
    QStringList allSwcPaths = getAllTargetPaths(swcInputDirPath, swcFilters);
    QDir swcInputDir(swcInputDirPath);
    for(auto swcPath:allSwcPaths){
        //获取相对路径
        QString relativePath = swcInputDir.relativeFilePath(swcPath);
        QString middlePath = "";
        int lastSlashIndex = relativePath.lastIndexOf("/");
        if(lastSlashIndex != -1){
            middlePath = relativePath.left(lastSlashIndex);
        }
        QString subLogDirPath;
        QString subResultDirPath;
        QString subSwcOutputDirPath;
        QString subSwcOutputAllDirPath;
        QString subSwcOutputDefinedDirPath;

        QString subResultSumCSVDirPath;

        QString croppedApoDirPath;
        QString croppedSwcDirPath;
        QString croppedApoDirMulfurPath;
        QString croppedApoDirBifurPath;
        QString croppedApoDirLoopPath;
        QString croppedApoDirMissingPath;
        QString croppedApoDirCrossingPath;
        QString croppedApoDirOverlapPath;
        QString croppedApoDirDissoPath;
        QString croppedApoDirAnglePath;
        QString croppedSwcDirMulfurPath;
        QString croppedSwcDirBifurPath;
        QString croppedSwcDirLoopPath;
        QString croppedSwcDirMissingPath;
        QString croppedSwcDirCrossingPath;
        QString croppedSwcDirOverlapPath;
        QString croppedSwcDirDissoPath;
        QString croppedSwcDirAnglePath;

        if(middlePath == ""){
            subLogDirPath = logDirPath;
            subResultDirPath = resultDirPath;
            subSwcOutputDirPath = swcOutputDirPath;
            subSwcOutputAllDirPath = swcOutputAllDirPath;
            subSwcOutputDefinedDirPath = swcOutputDefinedDirPath;
            subResultSumCSVDirPath = resultSumCSVDirPath;
        }
        else{
            subLogDirPath = logDirPath + "/" + middlePath;
            subResultDirPath = resultDirPath + "/" + middlePath;
            subSwcOutputDirPath = swcOutputDirPath + "/" + middlePath;
            subSwcOutputAllDirPath = swcOutputAllDirPath + "/" + middlePath;
            subSwcOutputDefinedDirPath = swcOutputDefinedDirPath + "/" + middlePath;
            subResultSumCSVDirPath = resultSumCSVDirPath + "/" + middlePath;
        }

        croppedApoDirPath = subSwcOutputDirPath + "/croppedApo";
        croppedSwcDirPath = subSwcOutputDirPath + "/croppedSwc";

        croppedApoDirMulfurPath = croppedApoDirPath + "/Multifurcation";
        croppedApoDirBifurPath = croppedApoDirPath + "/Approaching_bifurcation";
        croppedApoDirLoopPath = croppedApoDirPath + "/Loop";
        croppedApoDirMissingPath = croppedApoDirPath + "/Missing";
        croppedApoDirCrossingPath = croppedApoDirPath + "/Crossing";
        croppedApoDirOverlapPath = croppedApoDirPath + "/Overlap_seg";
        croppedApoDirDissoPath = croppedApoDirPath + "/Dissociative_seg";
        croppedApoDirAnglePath = croppedApoDirPath + "/Angle";

        croppedSwcDirMulfurPath = croppedSwcDirPath + "/Multifurcation";
        croppedSwcDirBifurPath = croppedSwcDirPath + "/Approaching_bifurcation";
        croppedSwcDirLoopPath = croppedSwcDirPath + "/Loop";
        croppedSwcDirMissingPath = croppedSwcDirPath + "/Missing";
        croppedSwcDirCrossingPath = croppedSwcDirPath + "/Crossing";
        croppedSwcDirOverlapPath = croppedSwcDirPath + "/Overlap_seg";
        croppedSwcDirDissoPath = croppedSwcDirPath + "/Dissociative_seg";
        croppedSwcDirAnglePath = croppedSwcDirPath + "/Angle";

        filesystem::create_directories(subLogDirPath.toStdString());
        filesystem::create_directories(subResultDirPath.toStdString());
        filesystem::create_directories(subSwcOutputDirPath.toStdString());
        if(needSort){
            filesystem::create_directories(subSwcOutputAllDirPath.toStdString());
            filesystem::create_directories(subSwcOutputDefinedDirPath.toStdString());
        }
        filesystem::create_directories(subResultSumCSVDirPath.toStdString());

//        filesystem::create_directories(croppedApoDirMulfurPath.toStdString());
//        filesystem::create_directories(croppedApoDirBifurPath.toStdString());
//        filesystem::create_directories(croppedApoDirLoopPath.toStdString());
//        filesystem::create_directories(croppedApoDirMissingPath.toStdString());
//        filesystem::create_directories(croppedApoDirCrossingPath.toStdString());
//        filesystem::create_directories(croppedApoDirOverlapPath.toStdString());
//        filesystem::create_directories(croppedApoDirDissoPath.toStdString());
//        filesystem::create_directories(croppedApoDirAnglePath.toStdString());

        QFileInfo swcFileInfo(swcPath);
        QString swcFileBaseName;
        if(swcPath.endsWith(".ano.eswc")){
            swcFileBaseName = swcFileInfo.fileName().left(swcFileInfo.fileName().size() - 9);
        }
        else{
            swcFileBaseName= swcFileInfo.completeBaseName();
        }
        QString inputPath = swcPath;

        QString logName = swcFileBaseName + ".log";
        QString logPath = subLogDirPath + "/" + logName;
        QString resultName = swcFileBaseName + ".txt";
        QString resultCSVName = swcFileBaseName + ".csv";
        QString resultPath = subResultDirPath + "/" + resultName;

        QString apoPath, anoPath, swcOutputPath;
        QString apoName = swcFileBaseName + ".ano.apo";
        QString anoName = swcFileBaseName + ".ano";
        QString swcName = swcFileBaseName + ".ano.eswc";

        if(needSort){
            apoPath = subSwcOutputAllDirPath + "/" + apoName;
            anoPath = subSwcOutputAllDirPath + "/" + anoName;
            swcOutputPath = subSwcOutputAllDirPath + "/" + swcName;
        }
        else{
            apoPath = subSwcOutputDirPath + "/" + apoName;
            anoPath = subSwcOutputDirPath + "/" + anoName;
            swcOutputPath = subSwcOutputDirPath + "/" + swcName;
        }

        QString croppedApoMulfurPath = croppedApoDirMulfurPath + "/" + swcFileBaseName + ".apo";
        QString croppedApoBifurPath = croppedApoDirBifurPath + "/" + swcFileBaseName + ".apo";
        QString croppedApoLoopPath = croppedApoDirLoopPath + "/" + swcFileBaseName + ".apo";
        QString croppedApoMissingPath = croppedApoDirMissingPath + "/" + swcFileBaseName + ".apo";
        QString croppedApoCrossingPath = croppedApoDirCrossingPath + "/" + swcFileBaseName + ".apo";
        QString croppedApoOverlapPath = croppedApoDirOverlapPath + "/" + swcFileBaseName + ".apo";
        QString croppedApoDissoPath = croppedApoDirDissoPath + "/" + swcFileBaseName + ".apo";
        QString croppedApoAnglePath = croppedApoDirAnglePath + "/" + swcFileBaseName + ".apo";

        QString croppedSwcMulfurPath = croppedSwcDirMulfurPath + "/" + swcFileBaseName;
        QString croppedSwcBifurPath = croppedSwcDirBifurPath + "/" + swcFileBaseName;
        QString croppedSwcLoopPath = croppedSwcDirLoopPath + "/" + swcFileBaseName;
        QString croppedSwcMissingPath = croppedSwcDirMissingPath + "/" + swcFileBaseName;
        QString croppedSwcCrossingPath = croppedSwcDirCrossingPath + "/" + swcFileBaseName;
        QString croppedSwcOverlapPath = croppedSwcDirOverlapPath + "/" + swcFileBaseName;
        QString croppedSwcDissoPath = croppedSwcDirDissoPath + "/" + swcFileBaseName;
        QString croppedSwcAnglePath = croppedSwcDirAnglePath + "/" + swcFileBaseName;

//        filesystem::create_directories(croppedSwcMulfurPath.toStdString());
//        filesystem::create_directories(croppedSwcBifurPath.toStdString());
//        filesystem::create_directories(croppedSwcLoopPath.toStdString());
//        filesystem::create_directories(croppedSwcMissingPath.toStdString());
//        filesystem::create_directories(croppedSwcCrossingPath.toStdString());
//        filesystem::create_directories(croppedSwcOverlapPath.toStdString());
//        filesystem::create_directories(croppedSwcDissoPath.toStdString());
//        filesystem::create_directories(croppedSwcAnglePath.toStdString());

        QString tmpSwcPath = subSwcOutputDirPath + "/" + swcFileBaseName + "_tmp.eswc";
        QString somaDefinedSwcPath = subSwcOutputDefinedDirPath + "/" + swcFileBaseName + ".ano.eswc.attachment.eswc";

//        DetectTask* task = new DetectTask();

        DetectTask* task = new DetectTask(swcPath, swcFileBaseName, logPath, apoPath, anoPath, swcOutputPath, somaDefinedSwcPath, tmpSwcPath, resultPath,
                                          croppedApoMulfurPath, croppedApoBifurPath, croppedApoLoopPath, croppedApoMissingPath,
                                          croppedApoCrossingPath, croppedApoOverlapPath, croppedApoDissoPath, croppedApoAnglePath,
                                          croppedSwcMulfurPath, croppedSwcBifurPath, croppedSwcLoopPath, croppedSwcMissingPath,
                                          croppedSwcCrossingPath, croppedSwcOverlapPath, croppedSwcDissoPath, croppedSwcAnglePath);
        threadPool.start(task);
    }

    // 等待线程池中的任务完成
    threadPool.waitForDone();

    // 汇总各个文件夹的结果
    map<QString, QString> map;
    QStringList csvFilters;
    csvFilters << "*.csv";
    QStringList allCsvPaths = getAllTargetPaths(resultDirPath, csvFilters);
    QDir resultDir(resultDirPath);
    for(auto csvPath:allCsvPaths){
        //获取相对路径
        QString relativePath = resultDir.relativeFilePath(csvPath);
        QString middlePath = "";
        int lastSlashIndex = relativePath.lastIndexOf("/");
        if(lastSlashIndex != -1){
            middlePath = relativePath.left(lastSlashIndex);
        }
        QString subResultDirPath;
        QString subResultSumDirPath;
        if(middlePath == ""){
            subResultDirPath = resultDirPath;
            subResultSumDirPath = resultSumCSVDirPath;
        }
        else{
            subResultDirPath = resultDirPath + "/" + middlePath;
            subResultSumDirPath = resultSumCSVDirPath + "/" + middlePath;
        }

        QString subResultSumCsvPath = subResultSumDirPath + "/" + "all.csv";
        map[subResultDirPath] = subResultSumCsvPath;
    }

    for(auto it = map.begin(); it != map.end(); it++){
        mergeResultCSVFiles(it->first.toStdString(), it->second.toStdString());
    }

    // 将结果汇总到一个文件中
    mergeResultCSVFilesAll(allCsvPaths, resultSumCSVAllPath.toStdString());
    // 将存在错误的结果汇总到一个文件中
//    mergeErrorResultCSVFilesAll(allCsvPaths, errorResultSumCSVAllPath.toStdString());
//    while(i<list.size()){
//        QFileInfo fileinfo = list.at(i);
//        QString infile = fileinfo.filePath();
//        qDebug()<<"infile:"<<infile;
//        if(!infile.endsWith(".swc")&&!infile.endsWith(".eswc"))
//        {
//            i++;
//            continue;
//        }
//        QString fileName = fileinfo.completeBaseName();
//        qDebug()<<"filename"<<fileName;
//        QString inputPath = infile;
//        QString logName = fileName + ".log";
//        QString logPath = logDirPath + "/" + logName;
//        QString resultName = fileName + ".txt";
//        QString resultPath = resultDirPath + "/" + resultName;
//        QString apoName = fileName + ".apo";
//        QString apoPath = apoDirPath + "/" + apoName;
//        QString anoName = fileName + ".ano";
//        QString anoPath = anoDirPath+ "/" + anoName;
//        QString swcOutputPath = swcOutputDirPath + "/" + fileinfo.fileName();
//        QString tmpSwcPath = swcOutputDirPath + "/" + fileName + "_tmp.eswc";
//        qDebug()<<tmpSwcPath;

//        CollDetection* task = new CollDetection(infile, fileName, logPath, apoPath, anoPath, swcOutputPath, tmpSwcPath, resultPath);
//        QStringList parts = fileName.split("_");
//        task->image = parts[0];
//        task->detectAll();
//        task->generateResult();
//        delete task;
//        i++;
//    }

    return 0;

}
