#!/bin/bash

rm FastChain_log

# Full chain with tpc, btof, and vpd (also uses the btofMixer)
root4star -b -q -l 'bfc.C ( 20, "sdt20140303.000000 TpcFastSim Simu y2014a sfs ssdfast McEvOut IdTruth McAna tpc_T globT tls db tpcDB ssdIT ITTF VFMinuit Idst event analysis EventQA tags EvOut StarMagField FieldOn IAna CMuDst vpdSim btofSim btofMixer btofMatch fzin","/star/data03/pwg/jdb/AuAu14p5/StarSim/999.fz", "FastChain" )' > FastChain_log
