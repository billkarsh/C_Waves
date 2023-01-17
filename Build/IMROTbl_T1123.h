#ifndef IMROTBL_T1123_H
#define IMROTBL_T1123_H

#include "IMROTbl_T1100.h"

/* ---------------------------------------------------------------- */
/* Types ---------------------------------------------------------- */
/* ---------------------------------------------------------------- */

// UHD phase 3 (layout 4) 12x32 (4.5um pitch)
//
struct IMROTbl_T1123 : public IMROTbl_T1100
{
    enum imLims_T1123 {
        imType1123Type      = 1123,
        imType1123Col       = 12
    };

    IMROTbl_T1123() {type=imType1123Type;}

    virtual int typeConst() const       {return imType1123Type;}
    virtual int nCol() const            {return imType1123Col;}
    virtual int nRow() const            {return imType1100Elec/imType1123Col;}

    virtual void locFltRadii( int &rin, int &rout, int iflt ) const;    // iflt = {1,2}
};

#endif  // IMROTBL_T1123_H


