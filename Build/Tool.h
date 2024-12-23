#ifndef TOOL_H
#define TOOL_H

#include <QFile>
#include <QString>
#include <QVector>

#include <limits>
#include <stdio.h>
#include <vector>

class KVParams;
struct GeomMap;
struct ShankMap;

class QFileInfo;

#define I64MAX  std::numeric_limits<qint64>::max()

/* ---------------------------------------------------------------- */
/* Types ---------------------------------------------------------- */
/* ---------------------------------------------------------------- */

struct TblRow {
    quint32  nspike, pkchan;
};

struct CTable {
    TblRow  *data;
    int     ndata;
};

struct Wrkspc {
    std::vector<int>    wsums;  // rows=nsamp, cols=nN
    std::vector<qint16> medpk;  // rows=nsamp, cols=maxwaves
    double              snr;
    qint64              s2;
    int                 lbl,
                        step,
                        cntdwn,
                        nsum;
    Wrkspc( int lbl, int nspike, int nNeurChan );
    bool wantIt();
    void addSpike(
        const std::vector<int>  &disk,
        const qint16            *src,
        int                     nC,
        int                     nN,
        int                     pkchan );
    void stats( const std::vector<int> &disk, int nN );
    void sampMedian( int it );
    void getMeans( std::vector<double> &D, double uV, int nN );
    void getMedian( std::vector<double> &D, double uV );
};

struct MyNPY {
    FILE                *fp;
    std::vector<char>   block;
    std::vector<size_t> shape;
    size_t              word_size;
    size_t              remVals;
    int                 nRead;
    char                type;
    MyNPY() : fp(0)     {}
    virtual ~MyNPY()    {if( fp ) fclose( fp );}
    bool openAndParseHdr( const char *path, int typeBytes );
    void parseHdr();
    bool readBlock();
};

class Tool
{
private:
    QString                 outpath;
    CTable                  clustbl;
    std::vector<Wrkspc>     vW;
    std::vector<int>        L2W;
    std::vector<std::vector<int> >  TSM;
    QFile                   *fbin;
    qint16                  *mbin,
                            *wbin;
    MyNPY                   pytim,
                            pylbl;
    double                  fCount,
                            uV;
    quint64                 file_ntpts;
    qint64                  startsmps,
                            endsmps;
    int                     smpBytes,
                            spikeBytes,
                            nC,
                            nN;

public:
    Tool()
        :   fbin(0), mbin(0), fCount(1.0),
            startsmps(-I64MAX), endsmps(I64MAX) {}
    virtual ~Tool();

    void entrypoint();

private:
    bool parseMeta();
    bool TSM_fromGeomMap( const QFileInfo &fim, const KVParams &kvp );
    bool TSM_fromShankMap( const QFileInfo &fim, const KVParams &kvp );
    bool getSavedChans(
        QVector<uint>   &snsFileChans,
        const QFileInfo &fim,
        const KVParams  &kvp );
    void snrTable_fromGeomMap( const GeomMap &GM );
    void snrTable_fromShankMap( const ShankMap &SM );
    void createWorkspaces();
    bool openFiles();
    bool getSpike( qint64 T );
    void sumWaves();
    void fixPath();
    void writeMeans();
    void writeMedians();
    void writeSNRs();
};

#endif  // TOOL_H




