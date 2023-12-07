#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export KB_AUTH_TOKEN=`cat /kb/module/work/token`
export PYTHONPATH=$script_dir/../lib:$PATH:$PYTHONPATH

# Set TEST_PATH to run a specific test. Eg: TEST_PATH=test.core.update_taxon_assignments_test
export TEST_PATH=.

cd $script_dir/../test
python -m nose --with-coverage --cover-package=GenomeFileUtil --cover-html --cover-html-dir=/kb/module/work/test_coverage --cover-xml --cover-xml-file=/kb/module/work/test_coverage/coverage.xml --nocapture --nologcapture $TEST_PATH
cp /kb/module/.coveragerc .
cp .coverage /kb/module/work/
mkdir -p /kb/module/work/kb/module/lib/
cp -R /kb/module/lib/GenomeFileUtil/ /kb/module/work/kb/module/lib/P
