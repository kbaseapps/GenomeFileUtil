#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export KB_AUTH_TOKEN=`cat /kb/module/work/token`
export PYTHONPATH=$script_dir/../lib:$PYTHONPATH
export PYTHONPATH=$script_dir/../test:$PYTHONPATH

# Set TEST_PATH to run a specific test. Eg: TEST_PATH=test.core.update_taxon_assignments_test
export TEST_PATH=test.problematic_tests.impl_test

cd $script_dir/../test
python -m nose --with-coverage --cover-package=GenomeFileUtil --cover-html --cover-html-dir=/kb/module/work/test_coverage --cover-xml --cover-xml-file=/kb/module/work/test_coverage/coverage.xml --nocapture --nologcapture $TEST_PATH
returncode=$?
if [ $returncode != 0 ]; then exit $returncode; fi
cp .coverage /kb/module/work/
mkdir -p /kb/module/work/kb/module/lib/
cp -R /kb/module/lib/GenomeFileUtil/ /kb/module/work/kb/module/lib/P
