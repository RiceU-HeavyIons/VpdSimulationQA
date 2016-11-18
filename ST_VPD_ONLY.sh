#!/bin/bash

rm vpdlog

# Quick, vpdSim Only chain to produce histograms
root4star -b -q -l 'bfc.C ( 20, "sdt20150303.000000 Simu McEvent event MuDst Tree vpdSim fzin","/star/data03/pwg/jdb/AuAu14p5/StarSim/999.fz", "Vpd_Only" )' > vpdlog
