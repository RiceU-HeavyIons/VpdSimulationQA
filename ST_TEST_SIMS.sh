#!/bin/bash

rm vpd_btof_log

# Quick, VpdSimMaker and BTofSimMaker chain to produce histograms
root4star -b -q -l 'bfc.C ( 20, "sdt20150303.000000 Simu McEvent event MuDst Tree vpdSim btofSim fzin","/star/data03/pwg/jdb/AuAu14p5/StarSim/999.fz", "Both" )' > vpd_btof_log
