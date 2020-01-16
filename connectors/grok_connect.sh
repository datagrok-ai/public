#!/usr/bin/env bash

# Builds "GrokConnect" server or runs as shell application
#
# Parameters:
#   - mode: test/run/shell
#
#
# "shell" mode:
#
#    Template:
#        bash grok_connect_shell.sh <mode> <args>
#
#    where:
#      * mode - run/build
#      * args:
#          -o, --output <arg>   Output CSV file path
#          -q, --query  <arg>   Query JSON file path
#
#    see query json example in "query_example.json"
#
#    Examples:
#        bash grok_connect.sh shell build -q examples/query.json -o output.csv
#        bash grok_connect.sh shell build -q examples/query.json
#        bash grok_connect.sh shell -q examples/query.json
#
DART_SDK=/usr/lib/dart
DART_PUB=${DART_SDK}/bin/pub

# Base project directory and executable
GROK_CONNECT_DIR=grok_connect
GROK_CONNECT=grok_connect-1.0.3.jar
TARGET_DIR=${GROK_CONNECT_DIR}/target

if [ "$1" == "shell" ]; then
    shift

    if [ "$1" == "build" ]; then
        shift
        rm -rf ${TARGET_DIR}
        mvn -Dmaven.test.skip=true package
    fi

    java -Xmx4g -classpath ${GROK_CONNECT_DIR}/lib/*:${TARGET_DIR}/${GROK_CONNECT} grok_connect.GrokConnectShell "$@"
else
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
fi
