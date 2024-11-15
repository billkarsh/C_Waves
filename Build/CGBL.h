#ifndef CGBL_H
#define CGBL_H

#include <QMap>
#include <QVector>

/* ---------------------------------------------------------------- */
/* Types ---------------------------------------------------------- */
/* ---------------------------------------------------------------- */

class CGBL
{
public:
    double          startsecs,
                    endsecs;
    QVector<uint>   vexc;
    const char      *sglbin,
                    *tblnpy,
                    *timenpy,
                    *lblnpy,
                    *dest,
                    *prefix;
    int             nsamp,
                    lhsamp,
                    maxwaves,
                    snrrad,
                    snrradum;
    bool            debug_npy;

public:
    CGBL()
        :   startsecs(-1.0), endsecs(-1.0), sglbin(0), tblnpy(0),
            timenpy(0), lblnpy(0), dest(0), prefix(0), nsamp(0),
            lhsamp(0), maxwaves(0), snrrad(8), snrradum(-1),
            debug_npy(false)    {}

    bool SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern CGBL GBL;

#endif  // CGBL_H


