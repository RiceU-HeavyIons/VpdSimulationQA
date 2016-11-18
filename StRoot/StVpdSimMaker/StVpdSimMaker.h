//
//  StVpdSimMaker.h
//
//
//  Created by Nickolas Luttrell on 6/3/16.
//
//

#ifndef StVpdSimMaker_HH
#define StVpdSimMaker_HH

#include "StMaker.h"
#include "St_DataSet.h"
class TH1F;
class TH2F;
class TH3F;
class TNtuple;
class StEvent;
class StBTofCollection;
class StVpdSimConfig;

#include "tables/St_g2t_vpd_hit_Table.h"
#include "StMcEvent/StMcEvent.hh"
#include "StMcEvent/StMcBTofHitCollection.hh"
#include "StVpdSimConfig.h"

#include <vector>
#ifndef ST_NO_NAMESPACES

#endif

class StVpdSimMaker : public StMaker {

public:
	StVpdSimMaker(const char *name = "VpdSim");
	virtual ~StVpdSimMaker();

	void           Reset();
	virtual Int_t  Init();
	Int_t          InitRun(Int_t);
	Int_t          FinishRun(Int_t);
	virtual Int_t  Make();
	virtual Int_t  Finish();

	StBTofCollection*  GetBTofCollection()  const { return mBTofCollection; }
	StMcBTofHitCollection* GetMcBTofHitCollection() const { return mMcBTofHitCollection; }
    

	string   pullHistFileName();
    string   getParamsFileName() { return mParamsFileName; }

	virtual const char *GetCVS() const{
		static const char cvs[] = "Tag $Name:  $ $Id: StVpdSimMaker.h,v 1.0 2016/06/03 00:00:00 nluttrel Exp $ built " __DATE__ " " __TIME__; return cvs;
	}

protected:


	StMcBTofHitCollection * mMcBTofHitCollection = nullptr;
	St_DataSet            * mGeantData           = nullptr;
	StEvent               * mEvent               = nullptr;
	StMcEvent             * mMcEvent             = nullptr;
	StBTofCollection      * mBTofCollection      = nullptr;
    std::map<Int_t, StVpdSimConfig::SingleTubeParams>	mSimParams;

	struct VpdSingleHit {
		Int_t tray;
		Int_t tubeId;
		Double_t tof;
		Double_t t0;
		Double_t de;
		Double_t pathL;
		Double_t q;
	};

	//  Various general use variables

	StVpdSimConfig* mSimConfig;
	Bool_t mBookHisto;          // Default is kFALSE.
    Bool_t mUseFileParameters;        // Default is kFALSE
    string mParamsFileName;

	Int_t       mCounter;
	Int_t       mNHits;
	Double_t    mVx;
	Double_t    mVy;
	Double_t    mVz;
	Double_t    mSumTubeTime;           // Tracks the time measured by each Vpd tube and sums them
	Double_t    mTubeTAvg;           // Average time lapse seen by the east or west Vpd

	Double_t    mTStart;
	Double_t    mTubeTAvgWest;      // Corrected event time for the west Vpd
	Double_t    mTubeTAvgEast;      // Corrected event time for the east Vpd
	Float_t     mVpdVertex;         // The calculated vertex as seen by the Vpd

	// Histogram variables
	string   mHistoFileName;    //! histogram file name
	TH1F* mNRawHitsWest;     // Number of hits on each west Vpd tube before threshold cuts
	TH1F* mNRawHitsEast;
	TH1F* mTubeHitsWest;     // Number of hits on each west Vpd tube after threshold cuts
	TH1F* mTubeHitsEast;
	TH1F* mNHitsWest;        // Number of tubes hit across events
	TH1F* mNHitsEast;
	TH1F* mLeTimeWest;      // Leading edge times (currently equal to Time of Flight)
	TH1F* mLeTimeEast;
	TH1F* mTStartHist;       // mTStart times (in ns)
	TH2F* mLeTubeHitsWest;  // Leading edge times and Number of tubes hit across events for West Vpd
	TH2F* mLeTubeHitsEast;
	TH1F* mZVertexHist;     // Histogram of the provided Mc Vertices for all events
	TH1F* mVpdVertexHist;   // Histogram of the calculated Vpd vertices for all events
	TH1F* mVpdResHist;       // Histogram of the difference between Mc and calculated vertices
	TH3F* mResVsNumHits;

							// Classes

	Int_t        vpdResponse(VpdSingleHit &Hit, g2t_vpd_hit_st* vpd_hit, Int_t vId);
	Double_t     thresholdCut(std::vector<VpdSingleHit> Hits, std::vector<Int_t> Tube_Counts, TH1F* TubeHits, TH1F* NHits);
	Int_t        storeMcVpdHit(StMcBTofHit* mcVpdHit);
	Int_t        fillEvent();
	Int_t        bookHistograms();
	Int_t        writeHistograms();
    
    void    setParamsFile(string fileName = "db/vpdSimParams/vpdSimParams.dat") { mParamsFileName = fileName; }

	ClassDef(StVpdSimMaker, 2)
};

#endif /* StVpdSimMaker_h */



// end of StVpdSimMaker.h
