
#include "CGBL.h"
#include "Util.h"
#include "Tool.h"


int main( int argc, char *argv[] )
{
    setLogFileName( "C_Waves.log" );

//double  T = getTime();

    if( !GBL.SetCmdLine( argc, argv ) ) {
        Log();
        return 42;
    }

    Tool    tool;
    tool.entrypoint();

//Log() << QString("Seconds: %1").arg( getTime() - T );

    Log();
    return 0;
}


