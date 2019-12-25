#!/usr/bin/env bash

# Build "GrokConnect" server
#
# Parameters:
#   - mode: test/run
#

DART_SDK=/usr/lib/dart
DART_PUB=${DART_SDK}/bin/pub

# Base project directory and executable
GROK_CONNECT_DIR=grok_connect
GROK_CONNECT=grok_connect-1.0.3.jar
TARGET_DIR=${GROK_CONNECT_DIR}/target

# Remove target
rm -rf ${TARGET_DIR}

if [ "$1" == "test" ]; then
    # Get dart test dependencies
    cd ../../ddt
    ${DART_PUB} get

    cd ../grok_shared
    ${DART_PUB} get

    cd ../ddt/bin/serialization_test
    ${DART_PUB} get

    cd ../../../public/connectors

    # Build project
    mvn package
else
    # Build project
    mvn -Dmaven.test.skip=true package
fi

# Pack into zip
ZIP_TMP_DIR=${TARGET_DIR}/grok_connect
mkdir ${ZIP_TMP_DIR}
echo "java -Xmx4g -classpath lib/*:${GROK_CONNECT} grok_connect.GrokConnect" > ${ZIP_TMP_DIR}/run_grok_connect.sh
cp ${TARGET_DIR}/${GROK_CONNECT} ${ZIP_TMP_DIR}/${GROK_CONNECT}
cp -R ${GROK_CONNECT_DIR}/lib ${ZIP_TMP_DIR}/lib
CUR_DIR=$(pwd)
cd ${ZIP_TMP_DIR}
zip -r ../grok_connect.zip *
cd ${CUR_DIR}
rm -rf ${ZIP_TMP_DIR}

# Run connector server with shared libraries
if [ "$1" == "run" ]; then
    java -Xmx4g -classpath ${GROK_CONNECT_DIR}/lib/*:${TARGET_DIR}/${GROK_CONNECT} grok_connect.GrokConnect
fi
