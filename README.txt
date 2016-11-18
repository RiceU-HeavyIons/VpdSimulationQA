README - StVpdSimMaker and StBTofMixerMaker

Below are included step-by-step commands needed to run the new simulation code. Note that ALL of these directories are needed for
the code to perform properly.

1) `starver SL16k`
2) `cvs co StRoot/StBFChain`
		`cvs co StRoot/StVpdSimMaker`
		`cvs co StRoot/StVpdCalibMaker`
		`cvs co StRoot/StBTofSimMaker`
		`cvs co StRoot/StBTofCalibMaker`
		`cvs co StRoot/StBTofMixerMaker`
		`cvs co StRoot/StTofUtil`
3) `cons`
4) `./ST_VPD_ONLY.sh`   or another shell script
5) Run root on output files

Testing Scripts:
	ST_VPD_ONLY.sh		Shell script to run the StVpdSimMaker only and produce output histograms (Vpd_Only.VpdSim.root)
	ST_TEST_SIMS.sh	  Shell script to run both the StVpdSimMaker and StBTofSimMaker and produce output histograms (Both.VpdSim.root & Both.BTofSim.root)
	ST_FASTCHAIN.sh		Shell script to run a chain with TpcFastSim, both vpd and btof makers, and the StBTofMixerMaker. Produces a MuDst.root file as well as the histogram files (FastChain.MuDst.root, FastChain.VpdSim.root, FastChain.BTofSim.root)

These scripts currently reference a particular .fz file (Au+Au 14.5 GeV UrQmd) that only has 20 events contained within it. There are lists of such files in:
lists/AuAu14p5.lis
lists/AuAu62p4.lis

Input and Output:
	The input .fz files use a uniform distribution between +/- 30 cm for the Primary Vertex. This should be apparent in the VpdSim.root histograms.
	Descriptions of histograms are present in header files.
	In general the MuDst produced should have branches populated for PrimaryVertex, PrimaryEvent, and BTofHit that are of interest.
