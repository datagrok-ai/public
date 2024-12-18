#!/bin/bash

set -xe

command=$1
expected_result=$2
exact_match=${3:-false}

result=$(eval "$command")

if [[ "$exact_match" == "true" ]]; then
    if [[ "$result" == "$expected_result" ]]; then
        exit 0
    else
        exit 1
    fi
else
    echo -e "${result}"
    echo -e "${result}" | grep -i "${expected_result}"
fi
