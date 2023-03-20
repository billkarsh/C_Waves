
#include "CGBL.h"
#include "Cmdline.h"
#include "Subset.h"
#include "Util.h"

#include <QFileInfo>


/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

CGBL    GBL;

/* --------------------------------------------------------------- */
/* PrintUsage  --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void PrintUsage()
{
    Log();
    Log() << "*** ERROR: MISSING CRITICAL PARAMETERS ***\n";
    Log() << "------------------------";
    Log() << "Purpose:";
    Log() << "+ Mean waveforms calculator.";
    Log() << "Run messages are appended to C_Waves.log in the current working directory.\n";
    Log() << "Output:";
    Log() << "+ 'mean_waveforms.npy': Array of mean waveforms (uV) with major-to-minor";
    Log() << "+    dims {nClusters, nChannels, nSamples}. nChannels matches the channel";
    Log() << "+    count of the spikeglx bin file. Non-neural 'waveforms' are zeroed.";
    Log() << "+ 'cluster_snr.npy': 3-col table of snr for the peak channels with cols";
    Log() << "+    {clus_lbl, snr, nSpikes_in_snr}.\n";
    Log() << "Usage:";
    Log() << ">C_Waves < required_parameters > [ options ]\n";
    Log() << "Required parameters:";
    Log() << "-spikeglx_bin=<filepath>    ;spikeglx metafile required in same dir";
    Log() << "-clus_table_npy=<filepath>  ;2-col table (uint32): {num_spike, pk-chan}";
    Log() << "-clus_time_npy=<filepath>   ;1-col table (uint64): {spike_time} (ASCENDING ORDER)";
    Log() << "-clus_lbl_npy=<filepath>    ;1-col table (uint32): {clus_lbl per spike_time}";
    Log() << "-dest=<path>                ;output dir (must exist)";
    Log() << "-samples_per_spike=82       ;waveform timepoints";
    Log() << "-pre_samples=30             ;subset of samples_per_spike to left of peak";
    Log() << "-num_spikes=1000            ;max waveforms included in averages";
    Log() << "-snr_radius=8               ;disk radius (rows/columns) about pk-chan, or zero\n";
    Log() << "Options:";
    Log() << "-snr_radius_um=140          ;disk radius (microns) about pk-chan, or zero";
    Log() << "-chnexcl=0,3:5              ;exclude these acq chans from snr";
    Log() << "-prefix=string              ;output files are named:";
    Log() << "                            ;  'prefix_mean_waveforms.npy'";
    Log() << "                            ;  'prefix_cluster_snr.npy'";
    Log() << "-debug_npy                  ;log 1st 10 {table, time, lbl} entries then quit";
    Log() << "------------------------\n";
}

/* ---------------------------------------------------------------- */
/* CGBL ----------------------------------------------------------- */
/* ---------------------------------------------------------------- */


bool CGBL::SetCmdLine( int argc, char* argv[] )
{
// Parse args

    const char  *sarg = 0;

    for( int i = 1; i < argc; ++i ) {

        if( GetArgStr( sglbin, "-spikeglx_bin=", argv[i] ) )
            ;
        else if( GetArgStr( tblnpy, "-clus_table_npy=", argv[i] ) )
            ;
        else if( GetArgStr( timenpy, "-clus_time_npy=", argv[i] ) )
            ;
        else if( GetArgStr( lblnpy, "-clus_lbl_npy=", argv[i] ) )
            ;
        else if( GetArgStr( dest, "-dest=", argv[i] ) )
            ;
        else if( GetArg( &nsamp, "-samples_per_spike=%d", argv[i] ) )
            ;
        else if( GetArg( &lhsamp, "-pre_samples=%d", argv[i] ) )
            ;
        else if( GetArg( &maxwaves, "-num_spikes=%d", argv[i] ) )
            ;
        else if( GetArg( &snrrad, "-snr_radius=%d", argv[i] ) )
            ;
        else if( GetArg( &snrradum, "-snr_radius_um=%d", argv[i] ) )
            ;
        else if( GetArgStr( sarg, "-chnexcl=", argv[i] ) ) {

            if( !Subset::rngStr2Vec( vexc, sarg ) )
                goto bad_param;
        }
        else if( GetArgStr( prefix, "-prefix=", argv[i] ) )
            ;
        else if( IsArg( "-debug_npy", argv[i] ) )
            debug_npy = true;
        else {
bad_param:
            Log() <<
            QString("Did not understand option '%1'.").arg( argv[i] );
            return false;
        }
    }

// Check args

    if( !sglbin || !tblnpy || !timenpy || !lblnpy || !dest
        || !nsamp || !lhsamp || !maxwaves ) {

        PrintUsage();
        return false;
    }

// Echo

    QString ssnrradum   = "",
            schnexc     = "",
            sprefix     = "",
            sdebug      = "";

    if( snrradum >= 0 )
        ssnrradum = QString(" -snr_radius_um=%1").arg( snrradum );

    if( vexc.size() )
        schnexc = " -chnexcl=" + Subset::vec2RngStr( vexc );

    if( prefix )
        sprefix = QString(" -prefix=%1").arg( prefix );

    if( debug_npy )
        sdebug = "-debug_npy";

    Log() <<
        QString(
        "Cmdline: C_Waves -spikeglx_bin=%1 -clus_table_npy=%2"
        " -clus_time_npy=%3 -clus_lbl_npy=%4 -dest=%5"
        " -samples_per_spike=%6 -pre_samples=%7 -num_spikes=%8"
        " -snr_radius=%9%10%11%12%13")
        .arg( sglbin )
        .arg( tblnpy )
        .arg( timenpy )
        .arg( lblnpy )
        .arg( dest )
        .arg( nsamp )
        .arg( lhsamp )
        .arg( maxwaves )
        .arg( snrrad )
        .arg( ssnrradum )
        .arg( schnexc )
        .arg( sprefix )
        .arg( sdebug );

    return true;
}


