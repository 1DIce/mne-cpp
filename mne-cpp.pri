############################################## GLOBAL FUNCTIONS ###############################################
#Define minQtVersion Test
defineTest(minQtVersion) {
    maj = $$1
    min = $$2
    patch = $$3
    isEqual(QT_MAJOR_VERSION, $$maj) {
        isEqual(QT_MINOR_VERSION, $$min) {
            isEqual(QT_PATCH_VERSION, $$patch) {
                return(true)
            }
            greaterThan(QT_PATCH_VERSION, $$patch) {
                return(true)
            }
        }
        greaterThan(QT_MINOR_VERSION, $$min) {
            return(true)
        }
    }
    greaterThan(QT_MAJOR_VERSION, $$maj) {
        return(true)
    }
    return(false)
}

defineReplace(MacDeployArgs) {
    target = $$1
    target_ext = $$2
    mne_binary_dir = $$3
    mne_library_dir = $$4
    extra_args = $$5

    isEmpty($${target_ext}) {
        target_custom_ext = .app
    } else {
        target_custom_ext = $${target_ext}
    }

    deploy_cmd = macdeployqt

    deploy_target = $$shell_quote($$shell_path($${mne_binary_dir}/$${target}$${target_custom_ext}))

    deploy_libs_to_copy = -libpath=$${mne_library_dir}
    !isEmpty(extra_args) {
      deploy_libs_to_copy += $${extra_args}
    }
    return($$deploy_cmd $$deploy_target $$deploy_libs_to_copy)
}

defineReplace(WinDeployArgs) {
    target = $$1
    target_ext = $$2
    mne_binary_dir = $$3
    libs_to_deploy = $$4
    extra_args = $$5

    isEmpty($${target_ext}) {
        target_custom_ext = .exe
    } else {
        target_custom_ext = $${target_ext}
    }

    deploy_cmd = windeployqt

    deploy_target = $$shell_quote($$shell_path($${mne_binary_dir}/$${target}$${target_custom_ext}))

    deploy_extra_args

    !isEmpty(extra_args) {
      deploy_extra_args += $${extra_args}
    }

    final_deploy_command += $$deploy_cmd $$deploy_target $$deploy_extra_args $$escape_expand(\\n\\t)

    # Parse libs from libs_to_deploy copy them to the bin folder and deploy qt dependecies for each of them
    for(FILE, libs_to_deploy) {
        FILE ~= s,-lMNE,MNE,g
        FILE = $${FILE}.dll
        TRGTDIR = $$shell_quote($$shell_path($${MNE_BINARY_DIR}))
        FILEPATH= $$shell_quote($$shell_path($${MNE_LIBRARY_DIR}/$${FILE}))

        exists($${FILEPATH}) {
            # Copy library
            final_deploy_command += $${QMAKE_COPY} $$quote($${FILEPATH}) $$quote($${TRGTDIR}) $$escape_expand(\\n\\t)

            # Deploy Qt dependencies for the library
            deploy_target = $$shell_quote($$shell_path($${TRGTDIR}/$${FILE}))
            final_deploy_command += windeployqt $${deploy_target} $$deploy_extra_args $$escape_expand(\\n\\t)

            #warning(Deploying $${FILEPATH} to $${TRGTDIR})
            #warning($${deploy_target})
        }
    }

    return( $${final_deploy_command} )
}


############################################### GLOBAL DEFINES ################################################

MNE_CPP_VERSION = 1.0.0
MNE_LIB_VERSION = 1

QMAKE_TARGET_PRODUCT = mne-cpp
QMAKE_TARGET_DESCRIPTION = MNE Qt 5 based C++ library.
QMAKE_TARGET_COPYRIGHT = Copyright (C) 2018 Authors of mne-cpp. All rights reserved.


########################################### PROJECT CONFIGURATION #############################################

## To build only the minimal version, i.e, for mne_rt_server run: qmake MNECPP_CONFIG+=minimalVersion
## To set CodeCov coverage compiler flag run: qmake MNECPP_CONFIG+=withCodeCov
## To disable tests run: qmake MNECPP_CONFIG+=noTests
## To disable examples run: qmake MNECPP_CONFIG+=noExamples
## To disable applications run: qmake MNECPP_CONFIG+=noApplications
## To build basic MNE Scan version run: qmake MNECPP_CONFIG+=buildBasicMneScanVersion
## To build MNE-CPP libraries as static libs: qmake MNECPP_CONFIG+=buildStaticLibraries

## Build MNE-CPP Deep library
MNECPP_CONFIG += buildDeep

#Build minimalVersion for qt versions < 5.10.0
!minQtVersion(5, 10, 0) {
    message("Building minimal version due to Qt version $${QT_VERSION}.")
    MNECPP_CONFIG += minimalVersion
}


########################################### DIRECTORY DEFINITIONS #############################################

# Eigen
EIGEN_INCLUDE_DIR = $$EIGEN_INCLUDE_DIR
isEmpty(EIGEN_INCLUDE_DIR) {
    EIGEN_INCLUDE_DIR = $${PWD}/include/3rdParty/eigen3
}

#CNTK
CNTK_INCLUDE_DIR = $$CNTK_INCLUDE_DIR
isEmpty( CNTK_INCLUDE_DIR ) {
    # Check CNTK Path options
    exists($$(CNTKPATH)/cntk/Include/Eval.h) {
        CNTK_TEST_DIR = $$(CNTKPATH)/cntk
    }
    exists($$(CNTKPATH)/Include/Eval.h) {
        CNTK_TEST_DIR = $$(CNTKPATH)
    }
    exists($$(MYCNTKPATH)/cntk/Include/Eval.h) {
        CNTK_TEST_DIR = $$(MYCNTKPATH)/cntk
    }
    exists($$(MYCNTKPATH)/Include/Eval.h) {
        CNTK_TEST_DIR = $$(MYCNTKPATH)
    }
    # Set CNTK path variables
    !isEmpty( CNTK_TEST_DIR ) {
        CNTK_INCLUDE_DIR = $${CNTK_TEST_DIR}/Include
        CNTK_LIBRARY_DIR = $${CNTK_TEST_DIR}/cntk
    }
}

# include
MNE_INCLUDE_DIR = $$MNE_INCLUDE_DIR
isEmpty( MNE_INCLUDE_DIR ) {
    MNE_INCLUDE_DIR = $${PWD}/libraries
}
MNE_SCAN_INCLUDE_DIR = $$MNE_SCAN_INCLUDE_DIR
isEmpty( MNE_SCAN_INCLUDE_DIR ) {
    MNE_SCAN_INCLUDE_DIR = $${PWD}/applications/mne_scan/libs
}
MNE_ANALYZE_INCLUDE_DIR = $$MNE_ANALYZE_INCLUDE_DIR
isEmpty( MNE_ANALYZE_INCLUDE_DIR ) {
    MNE_ANALYZE_INCLUDE_DIR = $${PWD}/applications/mne_analyze/libs
}
MNE_ANALYZE_EXTENSIONS_DIR = $$MNE_ANALYZE_EXTENSIONS_DIR
isEmpty( MNE_ANALYZE_EXTENSIONS_DIR ) {
    MNE_ANALYZE_EXTENSIONS_DIR = $${PWD}/applications/mne_analyze/extensions
}

# lib
MNE_LIBRARY_DIR = $$MNE_LIBRARY_DIR
isEmpty( MNE_LIBRARY_DIR ) {
    MNE_LIBRARY_DIR = $${PWD}/lib
}
contains(MNECPP_CONFIG, buildDeep) {
    CNTK_LIBRARY_DIR = $$CNTK_LIBRARY_DIR
    isEmpty( CNTK_LIBRARY_DIR ) {
        CNTK_LIBRARY_DIR = C:/local/cntk/cntk
    }
}

# bin
MNE_BINARY_DIR = $$MNE_BINARY_DIR
isEmpty( MNE_BINARY_DIR ) {
    MNE_BINARY_DIR = $${PWD}/bin
}

# repository dir
ROOT_DIR = $${PWD}

# install
MNE_INSTALL_INCLUDE_DIR = $$MNE_INSTALL_INCLUDE_DIR
isEmpty( MNE_INSTALL_INCLUDE_DIR ) {
    MNE_INSTALL_INCLUDE_DIR = $${PWD}/include
}

