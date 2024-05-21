
#include "Tool.h"
#include "CGBL.h"
#include "cnpy.h"
#include "GeomMap.h"
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

    if( !nd )
        return;

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

bool MyNPY::openAndParseHdr( const char *path , int typeBytes )
{
    fp = fopen( path, "rb" );

    if( !fp ) {
        Log() << QString("Unable to open npy '%1'").arg( path );
        return false;
    }

    parseHdr();

    if( shape.size() == 1 )
        ;
    else if( shape.size() != 2 || shape[1] > 1 ) {
        Log() << QString("NPY file must have dimensions: (N spikes)X(1 col) '%1'.")
                    .arg( path );
        return false;
    }

    if( typeBytes == 8 ) {
        // times
        if( !(type == 'i' || type == 'u') || word_size != 8 ) {
            Log() <<
            QString("NPY file must be int64 or uint64 instead of <%1%2> '%3'.")
            .arg( type ).arg( 8*word_size ).arg( path );
            return false;
        }
    }
    else {
        // labels
        if( type != 'u' || word_size != 4 ) {
            Log() <<
            QString("NPY file must be data type uint32 instead of <%1%2> '%3'.")
            .arg( type ).arg( 8*word_size ).arg( path );
            return false;
        }
    }

    return true;
}


void MyNPY::parseHdr()
{
#define BLOCKVALS   256

    bool fortran_order;

    type = cnpy::parse_npy_header( fp, word_size, shape, fortran_order );

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

Tool::~Tool()
{
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

    if( tbl.shape.size() != 2 || tbl.shape[1] != 2 ) {
        Log() << QString("Cluster table must have dimensions: (N clusters)X(2 cols) '%1'.")
                    .arg( fi.filePath() );
        return;
    }

    if( tbl.type != 'u' || tbl.word_size != 4 ) {
        Log() << QString("Cluster table must be data type uint32 instead of <%1%2> '%3'.")
                    .arg( tbl.type ).arg( 8*tbl.word_size ).arg( fi.filePath() );
        return;
    }

    if( GBL.debug_npy ) {
        Log() << "Cluster Table: " << clustbl.ndata << " entries";
        for( int ic = 0, nc = qMin( size_t(10), tbl.shape[0] ); ic < nc; ++ic )
            Log() << QString("%1 %2").arg( clustbl.data[ic].nspike ).arg( clustbl.data[ic].pkchan );
    }

    if( !parseMeta() )
        return;

// -------------------
// Initialize workflow
// -------------------

    createWorkspaces();

    if( !openFiles() )
        return;

    // Spike-count agreement check
    {
        quint64  tblSpikes = 0;
        for( int ic = 0, nc = tbl.shape[0]; ic < nc; ++ic )
            tblSpikes += clustbl.data[ic].nspike;

        if( tblSpikes != pytim.shape[0] || tblSpikes != pylbl.shape[0] ) {
            Log() << QString("Spike-count disagreement between cluster-table(%1), times(%2), labels(%3).")
                        .arg( tblSpikes ).arg( pytim.shape[0] ).arg( pylbl.shape[0] );
            return;
        }
    }

    if( GBL.debug_npy ) {

        Log() << "Cluster Times: " << pytim.remVals << " entries";
        pytim.readBlock();
        qint64 *T = (qint64*)&pytim.block[0];
        for( int ic = 0, nc = qMin( 10, pytim.nRead ); ic < nc; ++ic )
            Log() << QString("%1").arg( T[ic] );

        Log() << "Cluster Labels: " << pylbl.remVals << " entries";
        pylbl.readBlock();
        quint32 *L = (quint32*)&pylbl.block[0];
        for( int ic = 0, nc = qMin( 10, pylbl.nRead ); ic < nc; ++ic )
            Log() << QString("%1").arg( L[ic] );

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

    KVParams::const_iterator    it_kvp = kvp.find( "imDatPrb_pn" );
    QString                     pn;

    if( it_kvp != kvp.end() )
        pn = it_kvp.value().toString();
    else if( kvp.contains( "imProbeOpt" ) )
        pn = "Probe3A";

    if( !pn.isEmpty() ) {
        IMROTbl *R = IMROTbl::alloc( pn );
        if( !R ) {
            Log() << QString("Unknown probe type '%1'.").arg( pn );
            return false;
        }
        R->fromString( 0, kvp["~imroTbl"].toString() );
        uV = 1E+6 * R->maxVolts() / R->maxInt() / R->apGain( 0 );
        delete R;
    }

// -------------------
// Total channel count
// -------------------

    nC = kvp["nSavedChans"].toInt();

    smpBytes    = nC * sizeof(qint16);
    spikeBytes  = GBL.nsamp * smpBytes;

// -----------------------------
// Parse neural channel count nN
// -----------------------------

    QStringList sl = kvp["snsApLfSy"].toString().split(
                        QRegExp("^\\s+|\\s*,\\s*"),
                        QString::SkipEmptyParts );
    nN = sl[0].toInt();

// --------
// SNR disk
// --------

// Default disks are radius=0

    TSM.clear();
    TSM.resize( nN );

// Select geom- or shank- map method

    if( GBL.snrradum < 0 ) {
old_way:
        if( !TSM_fromShankMap( fim, kvp ) )
            return false;
    }
    else if( GBL.snrradum == 0 )
        ;
    else if( kvp.find( "~snsGeomMap" ) != kvp.end() ) {
        if( !TSM_fromGeomMap( fim, kvp ) )
            return false;
    }
    else {
        Log() << QString("Missing ~snsGeomMap tag '%1'...Will try ShankMap.")
                    .arg( fim.fileName() );
        goto old_way;
    }

    return true;
}


bool Tool::TSM_fromGeomMap(  const QFileInfo &fim, const KVParams &kvp  )
{
    QVector<uint>   snsFileChans;

    if( !getSavedChans( snsFileChans, fim, kvp ) )
        return false;

// -------
// GeomMap
// -------

    GeomMap GM;

    KVParams::const_iterator    it_kvp = kvp.find( "~snsGeomMap" );

    GM.fromString( it_kvp.value().toString() );

// -----------------
// Excluded channels
// -----------------

    for( int ie = 0, ne = GBL.vexc.size(); ie < ne; ++ie ) {

        int ig = snsFileChans.indexOf( GBL.vexc[ie] );

        if( ig >= 0 )
            GM.e[ig].u = 0;
    }

// ---
// TSM
// ---

    snrTable_fromGeomMap( GM );

    return true;
}


bool Tool::TSM_fromShankMap(  const QFileInfo &fim, const KVParams &kvp  )
{
    if( GBL.snrrad == 0 )
        return true;

    QVector<uint>   snsFileChans;

    if( !getSavedChans( snsFileChans, fim, kvp ) )
        return false;

// --------
// ShankMap
// --------

    ShankMap    SM;

    KVParams::const_iterator    it_kvp = kvp.find( "~snsShankMap" );

    if( it_kvp != kvp.end() )
        SM.fromString( it_kvp.value().toString() );
    else {
        Log() << QString("Missing ~snsShankMap tag '%1'.")
                    .arg( fim.fileName() );
        return false;
    }

// -----------------
// Excluded channels
// -----------------

    for( int ie = 0, ne = GBL.vexc.size(); ie < ne; ++ie ) {

        int ig = snsFileChans.indexOf( GBL.vexc[ie] );

        if( ig >= 0 )
            SM.e[ig].u = 0;
    }

// ---
// TSM
// ---

    snrTable_fromShankMap( SM );

    return true;
}


bool Tool::getSavedChans(
    QVector<uint>   &snsFileChans,
    const QFileInfo &fim,
    const KVParams  &kvp )
{
    QString chnstr = kvp["snsSaveChanSubset"].toString();

    if( Subset::isAllChansStr( chnstr ) )
        Subset::defaultVec( snsFileChans, nC );
    else if( !Subset::rngStr2Vec( snsFileChans, chnstr ) ) {
        Log() << QString("Bad snsSaveChanSubset tag '%1'.")
                    .arg( fim.fileName() );
        return false;
    }

    return true;
}


// For each channel [0,nN), calculate neighborhood
// of indices into a timepoint's channels.
// - Disk with radius {GBL.snrradiusum}.
// - The list is sorted for cache friendliness.
//
void Tool::snrTable_fromGeomMap( const GeomMap &GM )
{
    float   R2 = GBL.snrradum * GBL.snrradum;

    for( int ig = 0; ig < nN; ++ig ) {

        const GeomMapDesc   &E = GM.e[ig];

        if( !E.u )
            continue;

        std::vector<int>    &V = TSM[ig];

        for( int ie = 0; ie < nN; ++ie ) {

            const GeomMapDesc   &e = GM.e[ie];

            if( e.u && e.s == E.s ) {

                float   dx = e.x - E.x,
                        dz = e.z - E.z;

                if( dx*dx + dz*dz <= R2 )
                    V.push_back( ie );
            }
        }

        std::sort( V.begin(), V.end() );
    }
}


// For each channel [0,nN), calculate an 8-way
// neighborhood of indices into a timepoint's channels.
// - Disk with radius {GBL.snrradius}.
// - The list is sorted for cache friendliness.
//
void Tool::snrTable_fromShankMap( const ShankMap &SM )
{
    QMap<ShankMapDesc,uint> ISM;
    SM.inverseMap( ISM );

    int R = GBL.snrrad;

    for( int ig = 0; ig < nN; ++ig ) {

        const ShankMapDesc  &E = SM.e[ig];

        if( !E.u )
            continue;

        // ----------------------
        // Fill with disk members
        // ----------------------

        std::vector<int>    &V = TSM[ig];

        int xL  = qMax( int(E.c)  - R, 0 ),
            xH  = qMin( uint(E.c) + R + 1, SM.nc ),
            yL  = qMax( int(E.r)  - R, 0 ),
            yH  = qMin( uint(E.r) + R + 1, SM.nr );

        for( int ix = xL; ix < xH; ++ix ) {

            for( int iy = yL; iy < yH; ++iy ) {

                QMap<ShankMapDesc,uint>::iterator   it;

                it = ISM.find( ShankMapDesc( E.s, ix, iy, 1 ) );

                if( it != ISM.end() )
                    V.push_back( it.value() );
            }
        }

        std::sort( V.begin(), V.end() );
    }
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

    if( !pytim.openAndParseHdr( GBL.timenpy, 8 ) )
        return false;

    if( !pylbl.openAndParseHdr( GBL.lblnpy, 4 ) )
        return false;

    return true;
}


bool Tool::getSpike( qint64 T )
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

        qint64  *T = (qint64*)&pytim.block[0];
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


