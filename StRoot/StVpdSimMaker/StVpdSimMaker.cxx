/***************************************************************************
 *
 * $Id$
 *
 * Author: Nickolas Luttrell (Rice University), November 2016
 ***************************************************************************
 *
 * Description: StVpdSimMaker.cpp   - This simulation maker was adapted from
 * existing code in the StBTofSimMaker circa 2016. It takes Vpd hits from
 * geant, determines the earliest hit in each tube, then passes these hits to
 * StEvent. For QA purposes it also calculates tStart and the Vpd Vertex.
 *
 ***************************************************************************
 *
 * $Log$
 *
 ***************************************************************************/

#include <Stiostream.h>
#include <string>         // std::string
#include <cstddef>         // std::size_t
#include <math.h>
#include "StVpdSimMaker.h"
#include "StVpdSimConfig.h"

#include <TRandom3.h>
#include <math.h>
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include <vector>
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TNtuple.h"

// Tables
#include "tables/St_g2t_vpd_hit_Table.h"

#include "Random.h"
#include "RanluxEngine.h"
#include "RandGauss.h"

#include "StMcEvent/StMcVertex.hh"
#include "StEvent/StPrimaryVertex.h"
#include "StMcEvent/StMcBTofHit.hh"
#include "StEventTypes.h"
#include "StEvent/StBTofCollection.h"
#include "StChain/StChainOpt.h"

static RanluxEngine engine;
static RandGauss ranGauss(engine);

TRandom3 randengine(0);

ClassImp(StVpdSimMaker)

//_____________________________________________________________________________
StVpdSimMaker::StVpdSimMaker(const char *name):StMaker(name)
{
    //set default values
    mMcBTofHitCollection = 0;
    mBookHisto=kTRUE;// histograms
    mUseFileParameters=kFALSE;
    Reset();
}


//_____________________________________________________________________________
StVpdSimMaker::~StVpdSimMaker()
{
    //Destructor, delete any DaqMaps or similar that are no longer needed.
    delete mSimConfig;
}

//_____________________________________________________________________________
Int_t StVpdSimMaker::Init()
{
    Reset();
    if (mBookHisto) {
        bookHistograms();
    }
    return StMaker::Init();
}

//_____________________________________________________________________________
void StVpdSimMaker::Reset()
{
    mGeantData = 0;
    mEvent = 0;
    mMcEvent = 0;
    mSimConfig = 0;
    delete mMcBTofHitCollection;
}
//_____________________________________________________________________________
Int_t StVpdSimMaker::InitRun(Int_t runnumber)
{
    LOG_DEBUG << "StVpdSimMaker::InitRun  -- reading configuration file --" << endm;

    mSimConfig = new StVpdSimConfig;
    if (mUseFileParameters) {
        setParamsFile();
        mSimConfig->loadVpdSimParams(getParamsFileName());
        LOG_DEBUG << "Vpd Simulation Parameters loaded from file." << endm;
    }
    else {
        mSimConfig->loadVpdSimParams();
        LOG_DEBUG << "Vpd Simulation Parameters loaded from database." << endm;
    }
    mSimParams = mSimConfig->getParams();
    return kStOK;
}

//_____________________________________________________________________________
Int_t StVpdSimMaker::FinishRun(Int_t runnumber)
{

    return kStOk;
}


//_____________________________________________________________________________
Int_t StVpdSimMaker::Finish()
{
    if (mBookHisto) {
        LOG_DEBUG << "StVpdSimMaker::Finish  writing VpdSim.root ..." << endm;
        mHistoFileName = pullHistFileName();
        if (mHistoFileName == "") {
            LOG_DEBUG << "Nothing stored in mHistoFileName!" << endm;
            mHistoFileName = "VPDSim.root";
        }
        else {
            LOG_DEBUG << "The Filename is: " << mHistoFileName.c_str() << endm;
        }
        TFile aFile(mHistoFileName.c_str(),"RECREATE","vpdsim");
        aFile.cd();
        writeHistograms();
        aFile.Close();
    }
    return kStOK;
}


//_____________________________________________________________________________
Int_t StVpdSimMaker::Make()
{
    mMcBTofHitCollection = new StMcBTofHitCollection();

    VpdSingleHit singleHit;
    Int_t nWest = 0;
    Int_t nEast = 0;
    Double_t tubeTimeWest = 0;
    Double_t tubeTimeEast = 0;
    std::vector<VpdSingleHit> hitsWest;
    std::vector<VpdSingleHit> hitsEast;

    std::vector<Int_t> tubeCountsWest(19);
    std::vector<Int_t> tubeCountsEast(19);

    // Check to see that there are GEANT hits
    mGeantData = GetInputDS("geant"); // in bfc chain
    if(!mGeantData) { // when reading the geant.root file
        mGeantData = GetInputDS("geantBranch");
    }
    if(!mGeantData) {
        LOG_WARN << " No GEANT data loaded. Exit! " << endm;
        return kStWarn;
    }
    LOG_INFO << " Found GEANT data -- loading Vpd hits... " << endm;

    mMcEvent = (StMcEvent*)GetInputDS("StMcEvent");
    if (!mMcEvent) {
        LOG_ERROR << "No StMcEvent! Bailing out ..." << endm;
    }
    else if (mMcEvent->primaryVertex() != NULL){
        StThreeVectorF Vtx = mMcEvent->primaryVertex()->position();
        mVz = Vtx.z();       // VertexZ in cm
        LOG_DEBUG << "The simulated vZ is: " << mVz << endm;
    }

    // Look for Vpd hits
    St_g2t_vpd_hit* g2t_vpd_hits = 0;
    g2t_vpd_hits = dynamic_cast<St_g2t_vpd_hit*>(mGeantData->Find("g2t_vpd_hit"));
    if(!g2t_vpd_hits){
        LOG_WARN << " No Vpd hits in GEANT" << endm;
    }
    else {
        Int_t nVpdHits = g2t_vpd_hits->GetNRows();
        LOG_DEBUG << " Found Vpd hits: " << nVpdHits << endm;
        g2t_vpd_hit_st* vpd_hit = g2t_vpd_hits->begin();

        // Check for vpd_hit (only applies to the first one)
        if(!vpd_hit) {
            LOG_WARN << " No Vpd hit!" << endm;
            return kStWarn;
        }
        else {
            for(Int_t i=0; i<nVpdHits; i++, vpd_hit++) {
                Int_t vId = vpd_hit->volume_id;
                Int_t iTray = vId/1000 + 120;    //! 121: west, 122: east

                if ((iTray == 121) && (mSimParams[(vId % 100)-1].tubeStatusFlag == 1)) {   // West Vpd, check that tube is active
                    singleHit.tof = mSimConfig->getVpdDistance()*1e9/C_C_LIGHT - mVz*1e9/C_C_LIGHT;
                    vpdResponse(singleHit, vpd_hit, vId);
                    tubeCountsWest[singleHit.tubeId - 1] += 1;
                    hitsWest.push_back(singleHit);
                    if (mBookHisto) { mLeTimeWest->Fill(singleHit.tof); }
				}
                else if ((iTray == 122) && (mSimParams[(vId % 100)-1+19].tubeStatusFlag) == 1) { // East Vpd, check that tube is active
                    singleHit.tof = mSimConfig->getVpdDistance()*1e9/C_C_LIGHT + mVz*1e9/C_C_LIGHT;
					vpdResponse(singleHit, vpd_hit, vId);
					tubeCountsEast[singleHit.tubeId - 1] += 1;
					hitsEast.push_back(singleHit);
                    if (mBookHisto) { mLeTimeEast->Fill(singleHit.tof); }
				}
            }
        }
    }

    if (mBookHisto) {
        for(unsigned int i = 0; i < tubeCountsWest.size(); i++ ){
            mNRawHitsWest->Fill( i + 1,  tubeCountsWest[i]  );
            mLeTubeHitsWest->Fill(tubeCountsWest[i], singleHit.tof);
        }
        for(unsigned int i = 0; i < tubeCountsEast.size(); i++ ){
            mNRawHitsEast->Fill( i + 1,  tubeCountsEast[i] );
            mLeTubeHitsEast->Fill(tubeCountsEast[i], singleHit.tof);
        }
    }

    if(hitsWest.size() > 0) {       // Make sure there is at least one hit on the west.
        mTubeTAvgWest = thresholdCut(hitsWest, tubeCountsWest, mTubeHitsWest, mNHitsWest);
        tubeTimeWest = mSumTubeTime;
        nWest = mNHits;
    }
    else mTubeTAvgWest = -9999;

    if(hitsEast.size() > 0) {       // Make sure there is at least one hit on the east.
        mTubeTAvgEast = thresholdCut(hitsEast, tubeCountsEast, mTubeHitsEast, mNHitsEast);
        tubeTimeEast = mSumTubeTime;
        nEast = mNHits;
    }
    else mTubeTAvgEast = -9999;

    if ((mTubeTAvgEast == -9999) or (mTubeTAvgWest == -9999)) {
        mVpdVertex = -9999;
    }
    else {
        mVpdVertex = C_C_LIGHT*1.e-9*(mTubeTAvgEast - mTubeTAvgWest)/2;
        mTStart = ((tubeTimeEast+tubeTimeWest)-(nEast-nWest)*(mVz)/(C_C_LIGHT*1.e9))/(nEast+nWest) - mSimConfig->getVpdDistance()/(C_C_LIGHT*1.e-9);
        LOG_DEBUG << "The vpdVz is: " << mVpdVertex << " cm" << endm;
        LOG_DEBUG << "The tStart is " << mTStart << " ns" << endm;
    }

    fillEvent();

    if (mBookHisto) {
        mTStartHist->Fill(mTStart);
        mVpdVertexHist->Fill(mVpdVertex);
        mVpdResHist->Fill(mVz-mVpdVertex);

        mResVsNumHits->Fill(nWest, nEast, mVz-mVpdVertex);
    }

    return kStOK;
}



//_____________________________________________________________________________
// vpdResponse takes a single Vpd hit in a PMT and stores the relevant information
// for that hit, including tray, tubeId, and tof
//
Int_t StVpdSimMaker::vpdResponse(VpdSingleHit &singleHit, g2t_vpd_hit_st* vpd_hit, Int_t vId)
{
	Double_t randNum;

    singleHit.tray = vId/1000 + 120;    //! 121: west, 122: east
    singleHit.tubeId = vId%100;           //! 1-19

    if (singleHit.tray == 122) {
        randNum = randengine.Gaus(0, mSimParams[singleHit.tubeId-1+19].singleTubeRes);
    }
    else {
        randNum = randengine.Gaus(0, mSimParams[singleHit.tubeId-1].singleTubeRes);
    }

    singleHit.tof += randNum/1000;    // Applies single-tube resolution smearing
    singleHit.t0 = vpd_hit->tof*1./nanosecond;     // Stores value in nanoseconds.

    singleHit.de = vpd_hit->de;
    singleHit.pathL = -9999;
    singleHit.q = 0.;

    return kStOK;
}



//_____________________________________________________________________________
// thresholdCut applies cuts on the raw Hits(West/East) based on an assigned threshold value, then finds the earliest time in each tube and returns average time information from a given Vpd.

Double_t StVpdSimMaker::thresholdCut(std::vector<VpdSingleHit> singleHitsVec, std::vector<Int_t> tubeCounts, TH1F* TubeHits, TH1F* NHits) {

    std::vector<Double_t> timesVec;
    mCounter = 0;
    mSumTubeTime = 0;
    mTubeTAvg = 0;

    for(int i = 0; i < (int)tubeCounts.size(); i++) {   // Iterate through each of the 19 tubes
        Int_t dex = -1;
        Double_t lowTime = 1.e+9;

        if (tubeCounts[i] >= mSimConfig->getThreshold()) {     // Check if the tube count is over the minimum threshold (default is 1 hit)
            if (mBookHisto) { TubeHits->Fill( i + 1,  1  ); }
            for(int j = 0; j < (int)singleHitsVec.size(); j++) {  // Iterate through all the hits
                if ((i == singleHitsVec[j].tubeId-1) and (singleHitsVec[j].tof > 1.e-4)) {    // Find the hits corresponding to the current tube (that's over threshold)
                    if (singleHitsVec[j].tof < lowTime) {
                        lowTime = singleHitsVec[j].tof;  // Stores the lowest ToF amongst the hits in a tube
                        dex = j;
                    }
                }
            }
        }
        else if (mBookHisto) { NHits->Fill(0); } // Populate zero bin

        if (dex >= 0) {
            StMcBTofHit *mcHit = new StMcBTofHit(singleHitsVec[dex].tray,0,singleHitsVec[dex].tubeId,singleHitsVec[dex].de,singleHitsVec[dex].pathL,singleHitsVec[dex].t0,singleHitsVec[dex].tof,singleHitsVec[dex].q);
            storeMcVpdHit(mcHit); // Pass the hit with the lowest time value to store in the hit collection.

            mCounter += 1;
            mSumTubeTime += singleHitsVec[dex].tof;
            timesVec.push_back(singleHitsVec[dex].tof);
        }
    }

    if (timesVec.size() != 0) {        // This segment adapted from VpdCalibMaker::vzVpdFinder(), for QA purposes ONLY
        mNHits = timesVec.size();
//        Double_t TotalTime = mSumTubeTime;
//        for (Double_t timeVal : timesVec) {
//            Double_t vpdtime = (timeVal*timesVec.size() - TotalTime)/(timesVec.size() - 1);
//            if (fabs(vpdtime) > mSimConfig->getTDiffCut) {
//                mSumTubeTime -= timeVal;
//                mNHits -= 1;
//                mCounter -= 1;
//            }
//        }
        mTubeTAvg = mSumTubeTime/mNHits;
    }
    else { mTubeTAvg = -9999; }

    if (mBookHisto) { NHits->Fill(mCounter); }

    return mTubeTAvg;
}



//_____________________________________________________________________________
// storeMcVpdHit replaces duplicate hits (hits that match the same cell location),
// or it stores the new hit (the last part below)
Int_t StVpdSimMaker::storeMcVpdHit(StMcBTofHit* mcVpdHit)
{
    Bool_t hitFound = kFALSE;

    for(size_t j=0;j<mMcBTofHitCollection->hits().size();j++) {
        StMcBTofHit *tempHit = mMcBTofHitCollection->hits()[j];
        if(!tempHit) {
            LOG_DEBUG << " No hit stored in mMcBTofHitCollection." << endm;
        }
        if(mcVpdHit->sameCell(*tempHit)) {
            LOG_WARN << " Multiple hits passed to same cell. Exit! " << endm;
            return kStWarn;
        }
    }

    if(!hitFound) {
        LOG_DEBUG << " Storing mcVpdHit to Collection." << endm;
        mMcBTofHitCollection->addHit(mcVpdHit);
    } else {
        delete mcVpdHit;
    }
    return kStOk;
}


//___________________________________________________________________________
// fillEvent stores the btofCollection from McEvent into StEvent
//
Int_t StVpdSimMaker::fillEvent()
{
    LOG_DEBUG << "Filling McEvent and Event"<<endm;

    // send off to StMcEvent
    mMcEvent = (StMcEvent*)GetInputDS("StMcEvent");
    if (!mMcEvent) {
        LOG_ERROR << "No StMcEvent! Bailing out ..." << endm;
    }
    else {
        mMcEvent->setBTofHitCollection(mMcBTofHitCollection);       // Replaces existing collection with the passed argument
        LOG_INFO << " ... StMcVpdHitCollection stored in StMcEvent" << endm;
    }

    // send off to StEvent
    mEvent = (StEvent*)GetInputDS("StEvent");
    if (!mEvent) {
        LOG_ERROR << "No StEvent! Bailing out ..." << endm;
    }
    else { // mEvent non-zero
        LOG_DEBUG << "mEvent = " << mEvent << endm;

        // Check for the simulated vertex
        if (mMcEvent->primaryVertex() != NULL){
            StThreeVectorF Vtx = mMcEvent->primaryVertex()->position();
            mVx = Vtx.x();
            mVy = Vtx.y();
            mVz = Vtx.z();
            LOG_DEBUG << " mVx, mVy, vZ: "
            << mVx << " "
            << mVy << " "
            << mVz << " " << endm;
            if (mBookHisto) { mZVertexHist->Fill(mVz); }
        }

        //Store Collections
        mBTofCollection = mEvent->btofCollection();
        if(!mBTofCollection) {
            LOG_INFO << "Creating new StBTofCollection" << endm;
            mBTofCollection = new StBTofCollection();
            mEvent->setBTofCollection(mBTofCollection);
        }

        for(Int_t jj = 0; jj < (Int_t)mMcBTofHitCollection->hits().size(); jj++) {
            StMcBTofHit *aMcVpdHit = mMcBTofHitCollection->hits()[jj];

            if(!aMcVpdHit) continue;

            //Fill the StBTofHit
            StBTofHit aVpdHit;
            aVpdHit.Clear();

            Float_t mcTof = aMcVpdHit->tof();      //
            mcTof -= mSimConfig->getVpdDistance()*1e9/C_C_LIGHT;

            aVpdHit.setHardwarePosition(kBTofId);
            aVpdHit.setTray((Int_t)aMcVpdHit->tray());
            aVpdHit.setModule((unsigned char)aMcVpdHit->module());
            aVpdHit.setCell((Int_t)aMcVpdHit->cell());
            aVpdHit.setLeadingEdgeTime((Double_t)mcTof);
            aVpdHit.setTrailingEdgeTime((Double_t)mcTof);
            aVpdHit.setAssociatedTrack(NULL);          //done in StBTofMatchMaker
            aVpdHit.setIdTruth(aMcVpdHit->parentTrackId(), 1);  // Set qaTruth to 1 so simulated hits can be tracked in the Calib Makers
            mBTofCollection->addHit(new StBTofHit(aVpdHit));

            //Fill the StBTofRawHit
            StBTofRawHit aVpdRawHit;
            aVpdRawHit.Clear();
            aVpdRawHit.setTray((Int_t)aMcVpdHit->tray());
            aVpdRawHit.setChannel(6*(aMcVpdHit->module() - 1) + (Int_t)aMcVpdHit->cell());
            aVpdRawHit.setFlag(1);
            mBTofCollection->addRawHit(new StBTofRawHit(aVpdRawHit));
        }

        //Fill StBTofHeader --
        StBTofHeader aHead;
        mBTofCollection->setHeader(new StBTofHeader(aHead));

        LOG_INFO << "... StBTofCollection Stored in StEvent! " << endm;

    }
    return kStOK;
}


//_____________________________________________________________________________
// pullHistFileName uses the string argument from the chain being run to set
// the name of the output histogram file.
//
string StVpdSimMaker::pullHistFileName() {

    string extension = ".VpdSim.root";

    if (GetChainOpt()->GetFileOut() != NULL) {
        TString outFile = GetChainOpt()->GetFileOut();
        mHistoFileName = (string)outFile;
        size_t lastindex = mHistoFileName.find_last_of(".");
        mHistoFileName = mHistoFileName.substr(0, lastindex);
        lastindex = mHistoFileName.find_last_of("/");
        mHistoFileName = mHistoFileName.substr(lastindex+1, mHistoFileName.length());
        mHistoFileName = mHistoFileName + extension;
    }

    return mHistoFileName;
}



//_____________________________________________________________________________
Int_t StVpdSimMaker::bookHistograms()
{
    // Create new histograms

    mNRawHitsWest = new TH1F("mNRawHitsWest_tubeId","mNRawHitsWest vs. tubeId; Tube Number; # of Hits", 21, 0, 21);
    mNRawHitsEast  = new TH1F("mNRawHitsEast_tubeId","mNRawHitsEast vs. tubeId; Tube Number; # of Hits",21, 0, 21);

    mTubeHitsWest = new TH1F("mTubeHitsWest_tubeId","mTubeHitsWest vs. tubeId (Threshold Cut); Tube Number; # of Hits", 21, 0, 21);
    mTubeHitsEast = new TH1F("mTubeHitsEast_tubeId","mTubeHitsEast vs. tubeId (Threshold Cut); Tube Number; # of Hits",21, 0, 21);

    mNHitsWest = new TH1F("mNHitsWest","# Tubes Hit per Event West; # Tubes; Counts", 21, 0, 21);
    mNHitsEast = new TH1F("mNHitsEast","# Tubes Hit per Event East; # Tubes; Counts",21, 0, 21);

    mLeTimeWest = new TH1F("mLeTimeWest","Leading Edge Times West (No Cuts); LE (ns); Counts", 300, 0, 30);
    mLeTimeEast = new TH1F("mLeTimeEast","Leading Edge Times East (No Cuts); LE (ns); Counts", 300, 0, 30);
    mTStartHist = new TH1F("mTStart","mTStart; Time (ns); Counts", 300, -30, 30);

    mLeTubeHitsWest = new TH2F("mLeTubeHitsWest"," mTubeHitsWest & LE Times (No Cuts); LE (ns); # of Hits; Counts", 300, 0, 30, 25, 0, 25);
    mLeTubeHitsEast = new TH2F("mLeTubeHitsEast","mTubeHitsEast & LE Times (No Cuts); LE (ns); # of Hits; Counts", 300, 0, 30, 25, 0, 25);

    mZVertexHist = new TH1F("mZVertexHist","True Vertices; Position (cm); Counts", 300, -50, 50);
    mVpdVertexHist = new TH1F("mVpdVertexHist","Calculated Vpd Vertices; Position (cm); Counts", 300, -50, 50);
    mVpdResHist = new TH1F("mVpdResHist","True - Vpd Vertex Position; TruePos - CalcPos (cm); Counts", 300, -15, 15);

    mResVsNumHits = new TH3F("mResVsNumHits","mResVsNumHits; numTubesWest; numTubesEast; vZ-vpdVz", 21, 0, 21, 21, 0, 21, 300, -10, 10);

    return kStOK;
}


//_____________________________________________________________________________
Int_t StVpdSimMaker::writeHistograms()
{
    //Write Histograms

    mNRawHitsWest->Write();
    mNRawHitsEast->Write();

    mTubeHitsEast->Write();
    mTubeHitsWest->Write();

    mNHitsWest->Write();
    mNHitsEast->Write();

    mLeTimeWest->Write();
    mLeTimeEast->Write();
    mTStartHist->Write();
    mLeTubeHitsWest->Write();
    mLeTubeHitsEast->Write();

    mZVertexHist->Write();
    mVpdVertexHist->Write();
    mVpdResHist->Write();
    mResVsNumHits->Write();


return kStOk;
}


// End StVpdSimMaker.cpp
