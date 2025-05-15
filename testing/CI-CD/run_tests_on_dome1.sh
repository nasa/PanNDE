#!/bin/bash

source build_and_run_resources/dome1_env/env.sh

./bin/HostDataTests > HostDataTestResults.log
./bin/StubsTests > StubsTestResults.log
