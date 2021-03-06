
#include "IMROTbl_T0base.h"
#include "Util.h"

#include <QFileInfo>
#include <QStringList>
#include <QRegExp>
#include <QTextStream>

/* ---------------------------------------------------------------- */
/* struct IMRODesc_T0base ----------------------------------------- */
/* ---------------------------------------------------------------- */

int IMRODesc_T0base::chToEl( int ch ) const
{
    return (ch >= 0 ? ch + bank * 384 : 0);
}


// Pattern: "(chn bank refid apgn lfgn apflt)"
//
QString IMRODesc_T0base::toString( int chn ) const
{
    return QString("(%1 %2 %3 %4 %5 %6)")
            .arg( chn )
            .arg( bank ).arg( refid )
            .arg( apgn ).arg( lfgn )
            .arg( apflt );
}


// Pattern: "chn bank refid apgn lfgn apflt"
//
// Note: The chn field is discarded.
//
IMRODesc_T0base IMRODesc_T0base::fromString( const QString &s )
{
    const QStringList   sl = s.split(
                                QRegExp("\\s+"),
                                QString::SkipEmptyParts );

    return IMRODesc_T0base(
            sl.at( 1 ).toInt(), sl.at( 2 ).toInt(),
            sl.at( 3 ).toInt(), sl.at( 4 ).toInt(),
            sl.at( 5 ).toInt() );
}

/* ---------------------------------------------------------------- */
/* struct IMROTbl_T0base ------------------------------------------ */
/* ---------------------------------------------------------------- */

void IMROTbl_T0base::fillDefault()
{
    type = typeConst();
    e.clear();
    e.resize( nAP() );
}


void IMROTbl_T0base::fillShankAndBank( int shank, int bank )
{
    Q_UNUSED( shank )

    for( int i = 0, n = e.size(); i < n; ++i )
        e[i].bank = qMin( bank, maxBank( i ) );
}


// Return true if two tables are same w.r.t connectivity.
//
bool IMROTbl_T0base::isConnectedSame( const IMROTbl *rhs ) const
{
    const IMROTbl_T0base    *RHS    = (const IMROTbl_T0base*)rhs;
    int                     n       = nChan();

    for( int i = 0; i < n; ++i ) {

        if( e[i].bank != RHS->e[i].bank )
            return false;
    }

    return true;
}


// Pattern: (type,nchan)(chn bank refid apgn lfgn apflt)()()...
//
QString IMROTbl_T0base::toString() const
{
    QString     s;
    QTextStream ts( &s, QIODevice::WriteOnly );
    int         n = nChan();

    ts << "(" << type << "," << n << ")";

    for( int i = 0; i < n; ++i )
        ts << e[i].toString( i );

    return s;
}


// Pattern: (type,nchan)(chn bank refid apgn lfgn apflt)()()...
//
// Return true if file type compatible.
//
bool IMROTbl_T0base::fromString( QString *msg, const QString &s )
{
    QStringList sl = s.split(
                        QRegExp("^\\s*\\(|\\)\\s*\\(|\\)\\s*$"),
                        QString::SkipEmptyParts );
    int         n  = sl.size();

// Header

    QStringList hl = sl[0].split(
                        QRegExp("^\\s+|\\s*,\\s*"),
                        QString::SkipEmptyParts );

    if( hl.size() != 2 ) {
        type = -3;      // 3A type
        if( msg )
            *msg = "Wrong imro header size (should be 2)";
        return false;
    }

    type = hl[0].toInt();

    if( type != typeConst() ) {
        if( msg ) {
            *msg = QString("Wrong imro type[%1] for probe type[%2]")
                    .arg( type ).arg( typeConst() );
        }
        return false;
    }

// Entries

    e.clear();
    e.reserve( n - 1 );

    for( int i = 1; i < n; ++i )
        e.push_back( IMRODesc_T0base::fromString( sl[i] ) );

    if( e.size() != nAP() ) {
        if( msg ) {
            *msg = QString("Wrong imro entry count [%1] (should be %2)")
                    .arg( e.size() ).arg( nAP() );
        }
        return false;
    }

    return true;
}


bool IMROTbl_T0base::loadFile( QString &msg, const QString &path )
{
    QFile       f( path );
    QFileInfo   fi( path );

    if( !fi.exists() ) {

        msg = QString("Can't find '%1'").arg( fi.fileName() );
        return false;
    }
    else if( f.open( QIODevice::ReadOnly | QIODevice::Text ) ) {

        QString reason;

        if( fromString( &reason, f.readAll() ) ) {

            msg = QString("Loaded (type=%1) file '%2'")
                    .arg( type ).arg( fi.fileName() );
            return true;
        }
        else {
            msg = QString("Error: %1 in file '%2'")
                    .arg( reason ).arg( fi.fileName() );
            return false;
        }
    }
    else {
        msg = QString("Error opening '%1'").arg( fi.fileName() );
        return false;
    }
}


bool IMROTbl_T0base::saveFile( QString &msg, const QString &path ) const
{
    QFile       f( path );
    QFileInfo   fi( path );

    if( f.open( QIODevice::WriteOnly | QIODevice::Text ) ) {

        int n = f.write( STR2CHR( toString() ) );

        if( n > 0 ) {

            msg = QString("Saved (type=%1) file '%2'")
                    .arg( type )
                    .arg( fi.fileName() );
            return true;
        }
        else {
            msg = QString("Error writing '%1'").arg( fi.fileName() );
            return false;
        }
    }
    else {
        msg = QString("Error opening '%1'").arg( fi.fileName() );
        return false;
    }
}


int IMROTbl_T0base::elShankAndBank( int &bank, int ch ) const
{
    bank = e[ch].bank;
    return 0;
}


int IMROTbl_T0base::elShankColRow( int &col, int &row, int ch ) const
{
    int el = e[ch].chToEl( ch ),
        nc = nCol();

    row = el / nc;
    col = el - nc * row;

    return 0;
}


void IMROTbl_T0base::eaChansOrder( QVector<int> &v ) const
{
    QMap<int,int>   el2Ch;
    int             _nAP    = nAP(),
                    order   = 0;

    v.resize( 2 * _nAP + 1 );

// Order the AP set

    for( int ic = 0; ic < _nAP; ++ic )
        el2Ch[e[ic].chToEl( ic )] = ic;

    QMap<int,int>::iterator it;

    for( it = el2Ch.begin(); it != el2Ch.end(); ++it )
        v[it.value()] = order++;

// The LF set have same order but offset by nAP

    for( it = el2Ch.begin(); it != el2Ch.end(); ++it )
        v[it.value() + _nAP] = order++;

// SY is last

    v[order] = order;
}


// refid [0]    ext, shank=0, bank=0.
// refid [1]    tip, shank=0, bank=0.
// refid [2..4] int, shank=0, bank=id-2.
//
int IMROTbl_T0base::refTypeAndFields( int &shank, int &bank, int ch ) const
{
    int rid = e[ch].refid;

    shank = 0;

    if( rid == 0 ) {
        bank = 0;
        return 0;
    }
    else if( rid == 1 ) {
        bank = 0;
        return 1;
    }

    bank = rid - 2;
    return 2;
}


static int i2gn[IMROTbl_T0base::imType0baseGains]
            = {50,125,250,500,1000,1500,2000,3000};


bool IMROTbl_T0base::chIsRef( int ch ) const
{
    return ch == 191;
}


int IMROTbl_T0base::idxToGain( int idx ) const
{
    return (idx >= 0 && idx < 8 ? i2gn[idx] : i2gn[3]);
}


int IMROTbl_T0base::gainToIdx( int gain ) const
{
    switch( gain ) {
        case 50:    return 0;
        case 125:   return 1;
        case 250:   return 2;
        case 500:   return 3;
        case 1000:  return 4;
        case 1500:  return 5;
        case 2000:  return 6;
        case 3000:  return 7;
        default:    return 3;
    }
}


void IMROTbl_T0base::locFltRadii( int &rin, int &rout, int iflt ) const
{
    switch( iflt ) {
        case 2:     rin = 2, rout = 8; break;
        default:    rin = 0, rout = 2; break;
    }
}


void IMROTbl_T0base::muxTable( int &nADC, int &nGrp, std::vector<int> &T ) const
{
    nADC = 32;
    nGrp = 12;

    T.resize( imType0baseChan );

// Generate by pairs of columns

    int ch = 0;

    for( int icol = 0; icol < nADC; icol += 2 ) {

        for( int irow = 0; irow < nGrp; ++irow ) {
            T[nADC*irow + icol]     = ch++;
            T[nADC*irow + icol + 1] = ch++;
        }
    }
}


