#ifndef TOOL_H
#define TOOL_H

#include <QFile>
#include <QString>
#include <QVector>

#include <stdio.h>
#include <vector>

struct ShankMap;

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
        int                     nN );
    void stats( const std::vector<int> &disk, int nN );
    void getMeans( std::vector<double> &D, double uV, int nN );
};

struct MyNPY {
    FILE                *fp;
    std::vector<char>   block;
    std::vector<size_t> shape;
    size_t              word_size;
    size_t              remVals;
    int                 nRead;
    MyNPY() : fp(0)     {}
    virtual ~MyNPY()    {if( fp ) fclose( fp );}
    bool openAndParseHdr( const char *path );
    void parseHdr();
    bool readBlock();
};

class Tool
{
private:
    ShankMap                *shankMap;
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
    double                  uV;
    quint64                 file_ntpts;
    int                     smpBytes,
                            spikeBytes,
                            nC,
                            nN;

public:
    Tool();
    virtual ~Tool();

    void entrypoint();

private:
    bool parseMeta();
    void createWorkspaces();
    void snrTable( int nAP );
    bool openFiles();
    bool getSpike( quint64 T );
    void sumWaves();
    void fixPath();
    void writeMeans();
    void writeSNRs();
};

#endif  // TOOL_H




