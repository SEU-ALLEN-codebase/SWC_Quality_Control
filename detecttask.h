#ifndef DETECTTASK_H
#define DETECTTASK_H
#include <QDebug>
#include <QThreadPool>
#include <QRunnable>
#include "colldetection.h"


class DetectTask : public QRunnable
{
public:
    QString d_inFile;
    QString d_fileName;
    QString d_logPath;
    QString d_apoPath;
    QString d_anoPath;
    QString d_resultPath;
    QString d_swcOutputPath;
    QString d_tmpswcPath;
    QString d_croppedApoMulfurPath;
    QString d_croppedApoBifurPath;
    QString d_croppedApoLoopPath;
    QString d_croppedApoMissingPath;
    QString d_croppedApoCrossingPath;
    QString d_croppedApoOverlapPath;
    QString d_croppedApoDissoPath;
    QString d_croppedSwcMulfurPath;
    QString d_croppedSwcBifurPath;
    QString d_croppedSwcLoopPath;
    QString d_croppedSwcMissingPath;
    QString d_croppedSwcCrossingPath;
    QString d_croppedSwcOverlapPath;
    QString d_croppedSwcDissoPath;

public:
    DetectTask();
    DetectTask(QString infilePath, QString fileName, QString logPath, QString apoPath, QString anoPath, QString swcOutputPath, QString tmpSwcPath, QString resultPath,
               QString croppedApoMulfurPath, QString croppedApoBifurPath, QString croppedApoLoopPath, QString croppedMissingPath,
               QString croppedApoCrossingPath, QString croppedApoOverlapPath, QString croppedApoDissoPath, QString croppedSwcMulfurPath,
               QString croppedSwcBifurPath, QString croppedSwcLoopPath, QString croppedSwcMissingPath,
               QString croppedSwcCrossingPath, QString croppedSwcOverlapPath, QString croppedSwcDissoPath);
    QString getImageName(QString fileName);
    void run() override;
};

#endif // DETECTTASK_H
