#!/bin/bash
version=`python -c "from DevTools.Utilities.utilities import getCMSSWVersion; print getCMSSWVersion()"`
if [ "$version" == "76X" ]; then
    runPeriod="Collisions15"
elif [ "$version" == "80X" ]; then
    runPeriod="Collisions16"
else
    runPeriod="Collisions17"
fi
lumimask=`python -c "from DevTools.Utilities.utilities import getJson; print getJson('$runPeriod')"`
normtag=`python -c "from DevTools.Utilities.utilities import getNormtag; print getNormtag('$runPeriod')"`
pileupjson="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/$runPeriod/13TeV/PileUp/pileup_latest.txt"
mkdir -p pileup

brilcalc lumi -i $lumimask --normtag $normtag -b "STABLE BEAMS" -u /pb
