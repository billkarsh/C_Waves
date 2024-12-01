
TEMPLATE = app
TARGET   = C_Waves

win32 {
    CONFIG(debug, debug|release) {
        DESTDIR = C:/Users/labadmin/Desktop/SGLBUILD/FIXU/C_Waves/Debug
    }
    else {
        DESTDIR = C:/Users/labadmin/Desktop/SGLBUILD/FIXU/C_Waves/C_Waves-win
    }
}

unix {
    DESTDIR = /home/billkarsh/Code/C_Waves/C_Waves-linux
}

QT += widgets

HEADERS += \
    CGBL.h \
    Cmdline.h \
    cnpy.h \
    GeomMap.h \
    IMROTbl.h \
    IMROTbl_T0.h \
    IMROTbl_T0base.h \
    IMROTbl_T1020.h \
    IMROTbl_T1030.h \
    IMROTbl_T1100.h \
    IMROTbl_T1110.h \
    IMROTbl_T1120.h \
    IMROTbl_T1121.h \
    IMROTbl_T1122.h \
    IMROTbl_T1123.h \
    IMROTbl_T1200.h \
    IMROTbl_T1221.h \
    IMROTbl_T1300.h \
    IMROTbl_T2003.h \
    IMROTbl_T2013.h \
    IMROTbl_T2020.h \
    IMROTbl_T21.h \
    IMROTbl_T21base.h \
    IMROTbl_T24.h \
    IMROTbl_T24base.h \
    IMROTbl_T3010.h \
    IMROTbl_T3010base.h \
    IMROTbl_T3020.h \
    IMROTbl_T3020base.h \
    IMROTbl_T3A.h \
    KVParams.h \
    SGLTypes.h \
    ShankMap.h \
    Subset.h \
    Tool.h \
    Util.h

SOURCES += \
    CGBL.cpp \
    Cmdline.cpp \
    cnpy.cpp \
    GeomMap.cpp \
    IMROTbl.cpp \
    IMROTbl_T0base.cpp \
    IMROTbl_T1100.cpp \
    IMROTbl_T1110.cpp \
    IMROTbl_T1120.cpp \
    IMROTbl_T1121.cpp \
    IMROTbl_T1122.cpp \
    IMROTbl_T1123.cpp \
    IMROTbl_T1200.cpp \
    IMROTbl_T1221.cpp \
    IMROTbl_T2003.cpp \
    IMROTbl_T2013.cpp \
    IMROTbl_T2020.cpp \
    IMROTbl_T21.cpp \
    IMROTbl_T21base.cpp \
    IMROTbl_T24.cpp \
    IMROTbl_T24base.cpp \
    IMROTbl_T3010base.cpp \
    IMROTbl_T3020base.cpp \
    IMROTbl_T3A.cpp \
    KVParams.cpp \
    main.cpp \
    ShankMap.cpp \
    Subset.cpp \
    Tool.cpp \
    Util.cpp \
    Util_osdep.cpp

win32 {
    LIBS    += -lWS2_32 -lUser32 -lwinmm
    DEFINES += _CRT_SECURE_NO_WARNINGS WIN32
}

QMAKE_TARGET_COMPANY = Bill Karsh
QMAKE_TARGET_PRODUCT = C_Waves
QMAKE_TARGET_DESCRIPTION = Mean waveform calculator
QMAKE_TARGET_COPYRIGHT = (c) 2024, Bill Karsh, All rights reserved
VERSION = 2.9


