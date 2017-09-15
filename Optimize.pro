#### openCV ####
CONFIG += link_pkgconfig
PKGCONFIG += opencv

#### Eigen3 ####
INCLUDEPATH += /usr/local/include/eigen3

#### thirdparty Sophus ####
INCLUDEPATH += ../thirdparty/Sophus/

##### GLEW #####
LIBS += -lGLEW

##### OPENGL #####
QT += opengl

####### pangolin ########
INCLUDEPATH += /usr/local/include/pangolin
LIBS += -L"/usr/local/lib" -lpangolin

##### libzip #####
INCLUDEPATH += /usr/local/include/
LIBS += -L"/usr/local/lib" -lzip

HEADERS += \
    util/include_3000world.h \
    util/sse2_3000world.h \
    util/typedef_3000world.h \
    util/multithread_3000world.h \
    Frontend/Part.h \
    Backend/CalculateHessian/AccumulatedSCHessian.h \
    FullSystem/Settings.h \
    util/settings_3000world.h \
    util/all_util_include.h \
    Frontend/Mapping.h \
    Frontend/Tracking.h \
    Backend/BackPart.h \
    Backend/BackEnd.h \
    FullSystem/FullSystem.h \
    util/projection.h \
    util/interpolated_3000world.h \
    Backend/CalculateHessian/AccumulatedTopHessian.h \
    util/enum_3000world.h

SOURCES += \
    Frontend/Part.cpp \
    FullSystem/Settings.cpp \
    util/settings_3000world.cpp \
    Test/testmultithread.cpp \
    main/main.cpp \
    Test/cond_var.cpp \
    Test/concurrency.cpp \
    Test/parallel_accumulate.cpp \
    Backend/BackEnd.cpp \
    Frontend/Tracking.cpp \
    FullSystem/FullSystem.cpp \
    Frontend/Mapping.cpp \
    Backend/BackPart.cpp \
    Test/threadsafeDatastructure.cpp \
    Backend/BackOptimize.cpp \
    Backend/CalculateHessian/AccumulatedTopHessian.cpp \
    Backend/CalculateHessian/AccumulatedSCHessian.cpp

