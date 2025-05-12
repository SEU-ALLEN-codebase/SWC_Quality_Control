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
    QString somaDefinedSwcPath;

    QString croppedApoMulfurPath;
    QString croppedApoBifurPath;
    QString croppedApoLoopPath;
    QString croppedApoMissingPath;
    QString croppedApoCrossingPath;
    QString croppedApoOverlapPath;
    QString croppedApoDissoPath;
    QString croppedApoAnglePath;
    QString croppedSwcMulfurPath;
    QString croppedSwcBifurPath;
    QString croppedSwcLoopPath;
    QString croppedSwcMissingPath;
    QString croppedSwcCrossingPath;
    QString croppedSwcOverlapPath;
    QString croppedSwcDissoPath;
    QString croppedSwcAnglePath;

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
    int removedShortSegNum = 0;
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
    double consumedTime;
    double othersConsumedTime;
    double loopConsumedTime;
    double tipConsumedTime;
    double crossingConsumedTime;
    double overlapConsumedTime;
    double dissociativeConsumedTime;
    double angleConsumedTime;
    double colorConsumedTime;

public:
    vector<CellAPO> crossingAllMarkers;
    static QString swcInputDirPath;
    static bool needSort;

    struct TipCoorPredictedResult{
        QString storeDirName;
        XYZ maxResCoor;
        int y_pred;
        TipCoorPredictedResult(){storeDirName=""; maxResCoor=XYZ(); y_pred=-1;}
    };

    struct MissingForSegData{
        XYZ maxResCoor;
        QString storeDirName;
        XYZ centerCoor;
        XYZ edgeCoor;
        int type;
        V_NeuronSWC addedSeg;
        MissingForSegData(){maxResCoor=XYZ(); storeDirName=""; centerCoor=XYZ(); edgeCoor=XYZ();}
    };

    struct CrossingInfo{
        pair<pair<QString, int>, pair<QString, int>> coor2SegIndexPair;
        pair<pair<QString, pair<int, int>>, pair<QString, pair<int, int>>> coor2RowIndexRangePair;
        pair<pair<QString, vector<XYZ>>, pair<QString, vector<XYZ>>> fiberCoorInfoPair;
        bool isAbleCorrect = true;
    };

    struct FiberCoorData{
        int segID;
        pair<int, int> rangePair;
        vector<XYZ> coorVec;
        V_NeuronSWC addedSeg;
    };

    map<QString, MissingForSegData> tipInfoMap;
    map<QString, CrossingInfo> crossingInfoMap;
    bool isAutoCorrect = false;

    XYZ maxRes;
    XYZ subMaxRes;
    QString image;

    explicit CollDetection(QString infilepath, QString infilename, QString logpath, QString apopath, QString anopath, QString swcoutputpath, QString somadefinedswcpath, QString tmpswcPath, QString resultpath,
                           QString croppedApomulfurPath, QString croppedApobifurPath, QString croppedApoloopPath, QString croppedApomissingPath,
                           QString croppedApocrossingPath, QString croppedApooverlapPath, QString croppedApodissoPath, QString croppedApoAnglePath, QString croppedSwcmulfurPath,
                           QString croppedSwcbifurPath, QString croppedSwcloopPath, QString croppedSwcmissingPath,
                           QString croppedSwccrossingPath, QString croppedSwcoverlapPath, QString croppedSwcdissoPath, QString croppedSwcAnglePath, QObject* parent = nullptr);
    ~CollDetection(){logFile->flush(); logFile->close(); delete logFile;}
    XYZ getSomaCoordinate(QString apoPath);
    vector<NeuronSWC> specStructsDetection(V_NeuronSWC_list& inputSegList, double adjacent_dist_thre=0.2, double soma_dist_thre=8);
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
    void setSomaCondition();
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
    void removeShortSegs(V_NeuronSWC_list& segments, double dist_thres);
    void generateResult();
    void generateSortedSwc(QString somaDefinedSwcPath);
    void getApoAndCroppedSwc();
    void getApoForCrop(QString fileSaveName, vector<CellAPO> points);
    void getCropedSwc(QString fileInputPath, QString fileSavePath, XYZ coor1, XYZ coor2);
    void convertCoordInCropedSwc(QString filePath, XYZ coor1);

    void addMarkers(vector<CellAPO> cellApos);//加marker

};

#endif // COLLDETECTION_H
