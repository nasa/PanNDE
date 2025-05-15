#!/bin/bash

source build_and_run_resources/dome1_env/env.sh

export dt=$(date '+%Y_%m_%d_%H_%M_%S')
export buildId="build_$dt"

echo "build Id: $buildId"

export CC="sourceanalyzer -b $buildId gcc"
export CXX="sourceanalyzer -b $buildId g++"

echo $CC

cmake3 -B build -S . -DGO_TEST:bool=False
make -C build -j

echo "analyze"
sourceanalyzer -b $buildId -scan -f result.fpr
echo "make report"
ReportGenerator -format pdf -f report.pdf -source result.fpr -template ReportAll.xml