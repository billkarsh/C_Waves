
#include "Tool.h"
#include "CGBL.h"
#include "cnpy.h"
#include "IMROTbl.h"
#include "KVParams.h"
#include "ShankMap.h"
#include "Subset.h"
#include "Util.h"

#include <QFileInfo>

#include <math.h>


/* ---------------------------------------------------------------- */
/* Wrkspc --------------------------------------------------------- */
/* ---------------------------------------------------------------- */

// SNR(clus) = (Vmax-Vmin)/(2*sd)
//
// S = num spikes  = Wrkspc::nsum,
// C = num chans   = nNeurChan
// T = num samples = GBL.nsamp = "timepoints"
// W = SUMs(X)
// M = W/S
//
// sd  = sqrt(var)
// var = chisq(unit sigmas) for 2D model over C-T surface
// var = (1/(N-degfreedom))*SUM(residuals)^2
// var = (1/(S*C*T - C*T))*SUMsct(Xsct - Mct)^2
//
// SUMsct(Xsct - Mct)^2 =
// SUMsct(Xsct)^2 - S*SUMct(Mct)^2 =
// SUMsct(Xsct)^2 - (1/S)*SUMct(Wct)^2 =
// s2 - (1/S)*SUMct(Wct)^2
//

Wrkspc::Wrkspc( int lbl, int nspike, int nNeurChan )
    :   snr(-1), s2(0), lbl(lbl), cntdwn(1), nsum(0)
{
    wsums.resize( GBL.nsamp * nNeurChan, 0 );

    step = qMax( 1, int(ceil( float(nspike) / GBL.maxwaves )) );
}


bool Wrkspc::wantIt()
{
    if( --cntdwn > 0 )
        return false;

    cntdwn = step;
    ++nsum;

    return true;
}


void Wrkspc::addSpike(
    const std::vector<int>  &disk,
    const qint16            *src,
    int                     nC,
    int                     nN )
{
    const int   *dsk    = &disk[0];
    int         *sum    = &wsums[0];
    int         nd      = disk.size();

    for( int it = 0, nt = GBL.nsamp; it < nt; ++it, src += nC ) {

//--------------------------------------------------------------------
// Calculate noise only over first 15 points in tail
if( it < 15 ) {
//--------------------------------------------------------------------
        // SNR sum of squares over disk
        for( int id = 0; id < nd; ++id ) {
            int v = src[dsk[id]];
            s2 += v * v;
        }
//--------------------------------------------------------------------
// Calculate noise only over first 15 points in tail
}
//--------------------------------------------------------------------

        // mean wave sums over all nN
        for( int ic = 0; ic < nN; ++ic, ++sum )
            *sum += src[ic];
    }
}


void Wrkspc::stats( const std::vector<int> &disk, int nN )
{
    const int   *dsk    = &disk[0];
    int         *sum    = &wsums[0];
    qint64      sumsq   = 0;
    int         nd      = disk.size(),
                vmax    = -65536,
                vmin    = 65536;

// SNR

    for( int it = 0, nt = GBL.nsamp; it < nt; ++it, sum += nN ) {

        for( int id = 0; id < nd; ++id ) {

            int v = sum[dsk[id]];

//--------------------------------------------------------------------
// Calculate noise only over first 15 points in tail
if( it < 15 ) {
//--------------------------------------------------------------------
            sumsq += v * v;
//--------------------------------------------------------------------
// Calculate noise only over first 15 points in tail
}
//--------------------------------------------------------------------

            if( v > vmax )
                vmax = v;
            else if( v < vmin )
                vmin = v;
        }
    }

//--------------------------------------------------------------------
// Calculate noise only over first 15 points in tail
int temp = GBL.nsamp;
GBL.nsamp = 15;
//--------------------------------------------------------------------
    if( nsum > 1 ) {

        double  var = double((s2 - sumsq/nsum))
                        / ((nsum - 1) * nd * GBL.nsamp);

        if( var > 0.0 )
            snr = (vmax - vmin) / (2.0 * nsum * sqrt( var ));
    }
//--------------------------------------------------------------------
// Calculate noise only over first 15 points in tail
GBL.nsamp = temp;
//--------------------------------------------------------------------
}


void Wrkspc::getMeans( std::vector<double> &D, double uV, int nN )
{
    double  *dst    = &D[0];
    int     *src    = &wsums[0];
    int     ns      = GBL.nsamp;

    uV /= nsum; // convert sums to means here

    for( int ic = 0; ic < nN; ++ic ) {

        for( int is = 0; is < ns; ++is, ++dst )
            *dst = uV * src[nN*is + ic];
    }
}

/* ---------------------------------------------------------------- */
/* MyNPY ---------------------------------------------------------- */
/* ---------------------------------------------------------------- */

bool MyNPY::openAndParseHdr( const char *path )
{
    fp = fopen( path, "rb" );

    if( !fp ) {
        Log() << QString("Unable to open npy '%1'").arg( path );
        return false;
    }

    parseHdr();

    return true;
}


void MyNPY::parseHdr()
{
#define BLOCKVALS   256

    bool fortran_order;

    cnpy::parse_npy_header( fp, word_size, shape, fortran_order );

    remVals = 1;
    for( int is = 0, ns = shape.size(); is < ns; ++is )
        remVals *= shape[is];

    block.resize( BLOCKVALS * word_size );
}


bool MyNPY::readBlock()
{
    nRead = 0;

    int nBytes = word_size * qMin( remVals, size_t(BLOCKVALS) );

    if( !nBytes )
        return false;

    size_t nb = fread( &block[0], 1, nBytes, fp );

    nRead    = nb / word_size;
    remVals -= nRead;

    return nRead > 0;
}

/* ---------------------------------------------------------------- */
/* Tool ----------------------------------------------------------- */
/* ---------------------------------------------------------------- */

Tool::Tool() : shankMap(0), fbin(0), mbin(0)
{
}


Tool::~Tool()
{
    if( shankMap )
        delete shankMap;

    if( fbin ) {

        if( mbin )
            fbin->unmap( (uchar*)mbin );

        fbin->close();
    }
}


void Tool::entrypoint()
{
// --------------
// Organize input
// --------------

    fixPath();

    QFileInfo   fi( GBL.tblnpy );

    if( !fi.exists() ) {
        Log() << QString("Cluster table not found '%1'.").arg( fi.filePath() );
        return;
    }

    cnpy::NpyArray  tbl = cnpy::npy_load( GBL.tblnpy );
    clustbl.data  = tbl.data<TblRow>();
    clustbl.ndata = tbl.shape[0];

    if( GBL.debug_npy ) {
        Log() << "Cluster Table: " << clustbl.ndata << " entries";
        for( int i = 0; i < 10; ++i )
            Log() << QString("%1 %2").arg( clustbl.data[i].nspike ).arg( clustbl.data[i].pkchan );
    }

    if( !parseMeta() )
        return;

// -------------------
// Initialize workflow
// -------------------

    createWorkspaces();

    if( !openFiles() )
        return;

    if( GBL.debug_npy ) {

        Log() << "Cluster Times: " << pytim.remVals << " entries";
        pytim.readBlock();
        quint64 *T = (quint64*)&pytim.block[0];
        for( int i = 0; i < 10; ++i )
            Log() << QString("%1").arg( T[i] );

        Log() << "Cluster Labels: " << pylbl.remVals << " entries";
        pylbl.readBlock();
        quint32 *L = (quint32*)&pylbl.block[0];
        for( int i = 0; i < 10; ++i )
            Log() << QString("%1").arg( L[i] );

        return;
    }

// -------
// Process
// -------

    sumWaves();

    writeMeans();
    writeSNRs();
}


bool Tool::parseMeta()
{
    QString inMeta;
    QRegExp re("bin$");
    re.setCaseSensitivity( Qt::CaseInsensitive );
    inMeta = QString(GBL.sglbin).replace( re, "meta" );

    QFileInfo   fim( inMeta );

    if( !fim.exists() ) {
        Log() << QString("Meta file not found '%1'.").arg( fim.filePath() );
        return false;
    }

    KVParams    kvp;

    if( !kvp.fromMetaFile( inMeta ) ) {
        Log() << QString("Meta file is corrupt '%1'.").arg( fim.fileName() );
        return false;
    }

// --------------------------
// V = i * Vmax / Imax / gain
// --------------------------

    KVParams::const_iterator    it_kvp = kvp.find( "imDatPrb_type" );
    int                         prbType = -999;

    if( it_kvp != kvp.end() )
        prbType = it_kvp.value().toInt();
    else if( kvp.contains( "imProbeOpt" ) )
        prbType = -3;

    if( prbType != -999 ) {
        IMROTbl *R = IMROTbl::alloc( prbType );
        R->fromString( 0, kvp["~imroTbl"].toString() );
        uV = 1E+6 * kvp["imAiRangeMax"].toDouble() / R->maxInt() / R->apGain( 0 );
        delete R;
    }

// -------------------
// Total channel count
// -------------------

    nC = kvp["nSavedChans"].toInt();

    smpBytes    = nC * sizeof(qint16);
    spikeBytes = GBL.nsamp * smpBytes;

// -----------------------------
// Parse neural channel count nN
// -----------------------------

    QStringList sl = kvp["snsApLfSy"].toString().split(
                        QRegExp("^\\s+|\\s*,\\s*"),
                        QString::SkipEmptyParts );

    nN = sl[0].toInt();

// ---------------------
// Saved channel ID list
// ---------------------

    QVector<uint>   chanIds;
    QString         chnstr = kvp["snsSaveChanSubset"].toString();

    if( Subset::isAllChansStr( chnstr ) )
        Subset::defaultVec( chanIds, nC );
    else if( !Subset::rngStr2Vec( chanIds, chnstr ) ) {
        Log() << QString("Bad snsSaveChanSubset tag '%1'.")
                    .arg( fim.fileName() );
        return false;
    }

// ---------------------------------------------
// Graph (saved channel) to acq channel mappings
// ---------------------------------------------

    sl = kvp["acqApLfSy"].toString().split(
            QRegExp("^\\s+|\\s*,\\s*"),
            QString::SkipEmptyParts );

    int nAcqChan = sl[0].toInt() + sl[1].toInt() + sl[2].toInt();

    ig2ic.resize( nC );
    ic2ig.fill( -1, nAcqChan );

    for( int ig = 0; ig < nC; ++ig ) {

        int &C = ig2ic[ig];

        C           = chanIds[ig];
        ic2ig[C]    = ig;
    }

// --------
// ShankMap
// --------

    if( shankMap )
        delete shankMap;

    shankMap = new ShankMap;

    it_kvp = kvp.find( "~snsShankMap" );

    if( it_kvp != kvp.end() )
        shankMap->fromString( it_kvp.value().toString() );
    else {
        Log() << QString("Missing ~snsShankMap tag '%1'.")
                    .arg( fim.fileName() );
        return false;
    }

// -----------------
// Excluded channels
// -----------------

    for( int ie = 0, ne = GBL.vexc.size(); ie < ne; ++ie ) {

        // IMPORTANT:
        // User chan label not necessarily within span of
        // true channels, hence, span of ic2ig[], so don't
        // use ic2ig for this lookup.

        int ig = ig2ic.indexOf( GBL.vexc[ie] );

        if( ig >= 0 )
            shankMap->e[ig].u = 0;
    }

// ---
// TSM
// ---

    snrTable( nN );

    return true;
}


void Tool::createWorkspaces()
{
    for( int i = 0; i < clustbl.ndata; ++i ) {

        int nspike = clustbl.data[i].nspike;

        if( nspike ) {
            L2W.push_back( vW.size() );
            vW.push_back( Wrkspc( i, nspike, nN ) );
        }
        else
            L2W.push_back( -1 );
    }
}


// For each channel [0,nAP), calculate an 8-way
// neighborhood of indices into a timepoint's channels.
// - Disk with radius {GBL.snr_radius}.
// - The list is sorted for cache friendliness.
//
void Tool::snrTable( int nAP )
{
    TSM.clear();
    TSM.resize( nAP );

    QMap<ShankMapDesc,uint> ISM;
    shankMap->inverseMap( ISM );

    int R = GBL.snrrad;

    for( int ig = 0; ig < nAP; ++ig ) {

        const ShankMapDesc  &E = shankMap->e[ig];

        if( !E.u )
            continue;

        // ----------------------
        // Fill with disk members
        // ----------------------

        std::vector<int>    &V = TSM[ig];

        int xL  = qMax( int(E.c)  - R, 0 ),
            xH  = qMin( uint(E.c) + R + 1, shankMap->nc ),
            yL  = qMax( int(E.r)  - R, 0 ),
            yH  = qMin( uint(E.r) + R + 1, shankMap->nr );

        for( int ix = xL; ix < xH; ++ix ) {

            for( int iy = yL; iy < yH; ++iy ) {

                QMap<ShankMapDesc,uint>::iterator   it;

                it = ISM.find( ShankMapDesc( E.s, ix, iy, 1 ) );

                if( it != ISM.end() )
                    V.push_back( it.value() );
            }
        }

        qSort( V );
    }
}


bool Tool::openFiles()
{
// SGL

    fbin = new QFile( GBL.sglbin );

    if( !(fbin->open( QIODevice::ReadOnly )) ) {
        Log() << QString("Unable to open SGL binary %1").arg( GBL.sglbin );
        return false;
    }

    file_ntpts = fbin->size();

    if( !(mbin = (qint16*)fbin->map( 0, file_ntpts )) ) {
        Log() << "File mapping failed";
        return false;
    }

    file_ntpts /= smpBytes;

// NPY

    if( !pytim.openAndParseHdr( GBL.timenpy ) )
        return false;

    if( !pylbl.openAndParseHdr( GBL.lblnpy ) )
        return false;

    return true;
}


bool Tool::getSpike( quint64 T )
{
    if( T < GBL.lhsamp )
        return false;

    if( T - GBL.lhsamp + GBL.nsamp > file_ntpts )
        return false;

    wbin = &mbin[(T - GBL.lhsamp) * nC];

    return true;
}


void Tool::sumWaves()
{
// ----------
// Accumulate
// ----------

    while( pytim.readBlock() ) {

        pylbl.readBlock();

        quint64 *T = (quint64*)&pytim.block[0];
        quint32 *L = (quint32*)&pylbl.block[0];

        for( int it = 0; it < pytim.nRead; ++it, ++T, ++L ) {

            Wrkspc  &W = vW[L2W[*L]];

            if( !W.wantIt() )
                continue;

            if( getSpike( *T ) )
                W.addSpike( TSM[clustbl.data[*L].pkchan], wbin, nC, nN );
        }
    }

// -----
// Stats
// -----

    for( int iw = 0, nw = vW.size(); iw < nw; ++iw ) {
        Wrkspc  &W = vW[iw];
        W.stats( TSM[clustbl.data[W.lbl].pkchan], nN );
    }
}


void Tool::fixPath()
{
// ---------------------
// Trim trailing slashes
// ---------------------

    outpath = GBL.dest;

    QRegExp re("[/\\\\]+$");
    int     i;

    while( (i = outpath.indexOf( re )) > 0 )
        outpath.truncate( i );

// --------------------
// Only forward slashes
// --------------------

    outpath.replace( "\\", "/" );

// ------------------
// Add terminal slash
// ------------------

    outpath += "/";
}


void Tool::writeMeans()
{
// ----
// Open
// ----

    QString name;

    if( GBL.prefix )
        name = QString("%1%2_mean_waveforms.npy").arg( outpath ).arg( GBL.prefix );
    else
        name = QString("%1mean_waveforms.npy").arg( outpath );

    FILE *fp = fopen( STR2CHR( name ), "wb" );

    if( !fp ) {
        Log() << QString("Unable to open mean_waveforms.npy here '%1'").arg( outpath );
        return;
    }

// -------
// Buffers
// -------

    int n = nC * GBL.nsamp;

    std::vector<double> D( n, 0 ), Z( n, 0 );

// ------
// Header
// ------

    std::vector<size_t> shape;

    shape.push_back( clustbl.ndata );
    shape.push_back( nC );
    shape.push_back( GBL.nsamp );

    std::vector<char>   H = cnpy::create_npy_header<double>( shape );

    fwrite( &H[0], 1, H.size(), fp );

// ----
// Data
// ----

    for( int ik = 0; ik < clustbl.ndata; ++ik ) {

        int i = L2W[ik];

        if( i >= 0 ) {

            vW[i].getMeans( D, uV, nN );
            fwrite( &D[0], sizeof(double), n, fp );
        }
        else
            fwrite( &Z[0], sizeof(double), n, fp );
    }

// -----
// Close
// -----

    fclose( fp );
}


void Tool::writeSNRs()
{
// ----
// Open
// ----

    QString name;

    if( GBL.prefix )
        name = QString("%1%2_cluster_snr.npy").arg( outpath ).arg( GBL.prefix );
    else
        name = QString("%1cluster_snr.npy").arg( outpath );

    FILE *fp = fopen( STR2CHR( name ), "wb" );

    if( !fp ) {
        Log() << QString("Unable to open cluster_snr.npy here '%1'").arg( outpath );
        return;
    }

// -------
// Buffers
// -------

    std::vector<double> D( 2 * clustbl.ndata, 0 );

    double  *d = &D[0];

    for( int ic = 0; ic < clustbl.ndata; ++ic )
        d[2*ic] = -1.0;

// ------
// Header
// ------

    std::vector<size_t> shape;

    shape.push_back( clustbl.ndata );
    shape.push_back( 2 );

    std::vector<char>   H = cnpy::create_npy_header<double>( shape );

    fwrite( &H[0], 1, H.size(), fp );

// ----
// Data
// ----

    for( int iw = 0, nw = vW.size(); iw < nw; ++iw ) {

        const Wrkspc  &W = vW[iw];

        d[2*W.lbl]      = W.snr;
        d[2*W.lbl + 1]  = W.nsum;
    }

    fwrite( &D[0], sizeof(double), 2 * clustbl.ndata, fp );

// -----
// Close
// -----

    fclose( fp );
}


