#!/bin/bash

mvn clean compile assembly:single
cd ..
rm -rf bin
mkdir bin
cp src/target/tika-extractor-0.1.0-jar-with-dependencies.jar bin/tika-extractor.jar
