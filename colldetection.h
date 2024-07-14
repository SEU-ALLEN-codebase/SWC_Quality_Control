#ifndef COLLDETECTION_H
#define COLLDETECTION_H
#include <vector>
#include <set>
#include <unordered_set>
#include "utils.h"
#include "basic_c_fun/basic_surf_objs.h"
#include "neuron_editing/neuron_format_converter.h"
#include <QNetworkRequest>
#include <QEventLoop>
#include <QNetworkReply>
#include <QHttpMultiPart>
#include <QFile>
#include "json.hpp"
#include <QJsonArray>
#include <QJsonObject>
#include <QJsonDocument>
#include <algorithm>
#include <cmath>

class CollServer;
class CollDetection : public QObject
{
    Q_OBJECT
private:
    QNetworkAccessManager* accessManager;
    QString SuperUserHostAddress;
    QString BrainTellHostAddress;
    vector<NeuronSWC> tipPoints;
    QString inFile;
    QString inFilename;
    QString logPath;
    QString apoPath;
    QString anoPath;
    QString swcOutputPath;

    QString croppedApoMulfurPath;
    QString croppedApoBifurPath;
    QString croppedApoLoopPath;
    QString croppedApoMissingPath;
    QString croppedApoCrossingPath;
    QString croppedApoOverlapPath;
    QString croppedApoDissoPath;
    QString croppedSwcMulfurPath;
    QString croppedSwcBifurPath;
    QString croppedSwcLoopPath;
    QString croppedSwcMissingPath;
    QString croppedSwcCrossingPath;
    QString croppedSwcOverlapPath;
    QString croppedSwcDissoPath;

    QString resultPath;
    QString tmpInFile;
    QFile* logFile;
    QTextStream logOut;

    V_NeuronSWC_list segments;
    QList<CellAPO> markers;

    bool isSomaExists;
    bool isSortedSwc;

    XYZ somaCoordinate;

    int errorSegsNum = 0;
    int loopNum = 0;
    vector<CellAPO> overlapSegsMarkers;
    vector<CellAPO> mulFurcationMarkers;
    vector<CellAPO> nearBifurcationMarkers;
    vector<CellAPO> loopMarkers;
    vector<CellAPO> dissociativeSegsMarkers;
    vector<CellAPO> anglesMarkers;
    vector<CellAPO> tipUndoneMarkers;
    vector<CellAPO> crossingMarkers;
    vector<CellAPO> tipAllMarkers;
    vector<CellAPO> colorMutationMarkers;

public:
    vector<CellAPO> crossingAllMarkers;
    static QString swcInputDirPath;
    XYZ maxRes;
    XYZ subMaxRes;
    QString image;

    explicit CollDetection(QString infilepath, QString infilename, QString logpath, QString apopath, QString anopath, QString swcoutputpath, QString tmpswcPath, QString resultpath,
                           QString croppedApomulfurPath, QString croppedApobifurPath, QString croppedApoloopPath, QString croppedApomissingPath,
                           QString croppedApocrossingPath, QString croppedApooverlapPath, QString croppedApodissoPath, QString croppedSwcmulfurPath,
                           QString croppedSwcbifurPath, QString croppedSwcloopPath, QString croppedSwcmissingPath,
                           QString croppedSwccrossingPath, QString croppedSwcoverlapPath, QString croppedSwcdissoPath, QObject* parent = nullptr);
    ~CollDetection(){logFile->flush(); logFile->close(); delete logFile;}
    XYZ getSomaCoordinate(QString apoPath);
    vector<NeuronSWC> specStructsDetection(V_NeuronSWC_list& inputSegList, double adjacent_dist_thre=1.5, double soma_dist_thre=8);
    vector<NeuronSWC> loopDetection(V_NeuronSWC_list& inputSegList);
    vector<NeuronSWC> tipDetection(V_NeuronSWC_list& inputSegList, bool removeFlag, map<string, set<size_t>> allPoint2SegIdMap, double dist_thresh=20);
    QJsonArray crossingDetection();
    void handleMulFurcation(vector<NeuronSWC>& outputSpecialPoints, double dist_thre=8);
    void handleLoop(vector<NeuronSWC>& outputSpecialPoints);
    void handleNearBifurcation(vector<NeuronSWC>& bifurPoints);
    void handleTip(vector<NeuronSWC>& tipPoints);
    void handleCrossing(QJsonArray& json);
    void handleOverlapSegs(set<string>& outputSpecialPoints);

    void sortSWC(QString fileOpenName, QString fileSaveName, double thres=1000000000, V3DLONG rootid=1000000000);
    bool sortSWCAndDetectLoop(QString fileOpenName, QString fileSaveName, V3DLONG rootid=1000000000);
    void setSWCRadius(QString filePath, int r);
    QStringList getSWCSpecNInfo(QString filePath, int val);
    void getImageRES();
//    void getImageMaxRES();

    void detectAll();
    void detectLoops();
    void detectTips();
    void detectCrossings();
    void detectOthers();
    void detectOverlapSegs(V_NeuronSWC_list inputSegList, double dist_thres = 1);
    void removeErrorSegs(V_NeuronSWC_list& segments);
    void tuneErrorSegs(V_NeuronSWC_list& segments);
    void generateResult();
    void getApoAndCroppedSwc();
    void getApoForCrop(QString fileSaveName, vector<CellAPO> points);
    void getCropedSwc(QString fileInputPath, QString fileSavePath, XYZ coor1, XYZ coor2);
    void convertCoordInCropedSwc(QString filePath, XYZ coor1);

    void addMarkers(vector<CellAPO> cellApos);//加marker

};

#endif // COLLDETECTION_H
