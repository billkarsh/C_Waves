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
                    snrrad;
    bool            debug_npy;

public:
    CGBL()
        :   sglbin(0), tblnpy(0), timenpy(0), lblnpy(0),
            dest(0), prefix(0), nsamp(0), lhsamp(0),
            maxwaves(0), snrrad(8), debug_npy(false)    {}

    bool SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern CGBL GBL;

#endif  // CGBL_H


