/***************************************************************************
 *
 * $Id: StBTofSimMaker.h,v 1.6 2015/07/28 22:49:55 smirnovd Exp $
 *
 * Author:  Frank Geurts
 ***************************************************************************
 *
 * Description: StBTofSimMaker virtual base class for Barrel TOF Simulations
 *
 ***************************************************************************
 *
 * $Log: StBTofSimMaker.h,v $
 * Revision 1.6  2015/07/28 22:49:55  smirnovd
 *  Initialize static constants outside of class definition
 *
 *  C++ forbids initialization of non-integral static const members within the class
 *  definition. The syntax is allowed only for integral type variables.
 *
 * Revision 1.5  2014/08/06 11:42:54  jeromel
 * Suffix on literals need to be space (later gcc compiler makes it an error) - first wave of fixes
 *
 * Revision 1.4  2011/02/03 19:01:01  geurts
 * Introduce option to switch writing simulated hits to StEvent. Default behavior is OFF.
 *
 * Revision 1.3  2010/08/10 19:18:32  geurts
 * Look for geant data in bfc ("geant") or geant.root ("geantBranch"); Protect storing BTofMcHitCollection in case McEvent is NULL.  [Xin]
 *
 * Revision 1.2  2010/07/14 20:32:58  geurts
 * remove geometry initialization (not used)
 *
 * Revision 1.1  2009/12/09 21:56:41  dthein
 * First version of StBTofSimMaker
 * 
 *
 **************************************************************************/
#ifndef STBTOFSIMMAKER_HH
#define STBTOFSIMMAKER_HH
#include "StMaker.h"

#include "St_DataSet.h"
class TH1F;
class TH2F;
class TNtuple;

class StEvent;
class StBTofCollection;
class StTofSimParam;
class StBTofDaqMap;
struct g2t_ctf_hit_st;

// g2t tables
#include "tables/St_g2t_ctf_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_tpc_hit_Table.h"

#include "StMcEvent/StMcEvent.hh"
#include "StMcEvent/StMcBTofHitCollection.hh"
#include "StMcEvent/StMcBTofHit.hh"
#include "StThreeVectorF.hh"
#include <vector>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class StBTofSimMaker : public StMaker{
	protected:


		StTofSimParam*      mSimDb;          //!
		StBTofDaqMap*       mDaqMap;         //! Tof Daq map
		StMcBTofHitCollection *mMcBTofHitCollection; //! barrel tof hit

		St_DataSet        *mGeantData;        //! geant table
		StEvent           *mEvent;            //!
		StMcEvent         *mMcEvent;
		StBTofCollection   *mBTofCollection = nullptr;

		//define some constants
		enum {
			mNTOF = 192,    //! TOF channels per tray
			mNTray = 120,   //! 120 TOF trays
			mNModule = 32,  //! 32 modules per tray
			mNCell = 6,     //! 6 cells per module
			mAMP = 50000,     
			mADCBINWIDTH = 25,
			mTDCBINWIDTH = 50
		};
		static const float mVHRBIN2PS;    //! Very High resolution mode, ps/bin
		static const float mHRBIN2PS;     //! High resolution mode, ps/bin
		static const float mBTofPadWidth; //! Pad Width


		Bool_t mCellXtalk;     //! switch for cell xtalk
		Bool_t mSlow;           //If True, runs the slow Tof Simulation, including CellResponse and CellTimePassTh
		Bool_t mBookHisto;
		Bool_t mWriteStEvent;  //! switch to enable Maker to write out simulated hits to StEvent

		Int_t     mTofHitFlag[mNTray][mNTOF];   //! hit flag for tof geant hits

		struct TrackHit{
			Int_t          tray;
			Int_t          module;
			Int_t          cell;
			Int_t          trkId;
			Double_t       dE;
			Double_t       dQ;
			Double_t       dQdt[600];//this 600 (nTimebins) comes from the /TofUtil/StTofParam file
			Double_t       tof;
			Double_t       s_track;
			Double_t          t0;      //! t0 (in ps) as the start of this tof hit -- was ns - changed for consistency
			StThreeVectorF position;
		};


		typedef vector<TrackHit, allocator<TrackHit> > TrackVec;
		typedef vector<Int_t> IntVec;


		string mHistoFile;//for QA histograms
        string mHistoFileName;
    
        TNtuple* ntuple;
    
        TH1F* mRawBetaHist; //
		TH1F* mBetaHist;    //! speed of particles hitting tof
        TH2F* mRawBetaVsMom;
        TH2F* mCalcBetaVsMom;
        TH2F* mBetaVsMom;       // 1/beta vs momentum
    
        TH2F* Electron_BetaVsMom;
        TH2F* Muon_BetaVsMom;
        TH2F* Pion_BetaVsMom;
            TH2F* Kaon_BetaVsMom;
        TH2F* Proton_BetaVsMom;
        TH2F* maybeLambda_BetaVsMom;
    
        TH1F* mPathLHist;    //! speed of particles hitting tof
        TH1F* mRawTofHist;      //! total time of flight of particle before resolution smearing
        TH1F* mTofHist;    //! total time of flight of partilce
        TH1F* mRecMass;    //! reconstructed mass of particle
    
        TH1F* massHist;
        TH2F* m2VsP;        // Mass Squared versus momentum
        TH1F* mTofCalculated;
        TH2F* tof_RealVsCalc;
    
        TH1F* momBinRaw1;
        TH1F* momBinRaw2;
        TH1F* momBinRaw3;
        TH1F* momBinRaw4;
        TH1F* momBinRaw5;
        TH1F* momBinRaw6;
        TH1F* momBinRaw7;
        TH1F* momBinRaw8;
        
        TH1F* momBin1;
        TH1F* momBin2;
        TH1F* momBin3;
        TH1F* momBin4;
        TH1F* momBin5;
        TH1F* momBin6;
        TH1F* momBin7;
        TH1F* momBin8;

		TH2F* mCellGeant;    //! cellId of geant hit
		TH2F* mNCellGeant;   //! # of cells of geant hit
		TH1F* mDeGeant;      //! deposited-energy in geant hit
		TH1F* mTofGeant;     //! tof in geant hit

		TH2F* mCellSeen;     //! cellId after DetectorResponse
		TH2F* mNCellSeen;    //! # of cells after DetectorResponse
		TH1F* mDeSeen;       //! deposited-energy after DetectorResponse
		TH1F* mT0Seen;      //! 
		TH1F* mTofSeen;      //! smeared-tof after DetectorResponse
		TH1F* mTofResSeen;   //! time resolution after Detector Response

		TH2F* mCellReco;     //! cellId after recon
		TH2F* mNCellReco;    //! # of cells after recon
		TH1F* mTDCReco;      //! TDC recon
		TH1F* mADCReco;      //! ADC recon -- empty
		TH1F* mT0Reco;   //! 
		TH1F* mTofResReco;   //! time resolution after recon
		TH2F* mTACorr;       //! T-A Slewing Correlation
		TH1F* mModHist;       //! T-A Slewing Correlation

		/// TOFp histograms
		TH1F* mdE;           //!
		TH1F* mdS;           //!
		TH1F* mT;            //!
		TH1F* mTime;         //!
		TH1F* mTime1;        //!
		TH1F* mPMlength;     //!
		TH1F* mAdc;          //!
		TH1F* mTdc;          //!

		TVolume *starHall;

		Int_t CellResponse(g2t_ctf_hit_st* tof_hit,
				TrackVec& trackVec);   //! Slow simulation step one
		Int_t CellTimePassTh(TrackVec& trackVec);        //! Slow simulation step two

		Int_t FastCellResponse(g2t_ctf_hit_st* tof_hit);

		IntVec    CalcCellId(Int_t volume_id, Float_t ylocal);
		Int_t CellXtalk(Int_t icell, Float_t ylocal, Float_t& wt, Int_t& icellx);
		Int_t      storeMcBTofHit(StMcBTofHit *mcCellHit);

		Int_t        fillRaw(void);
		Int_t        electronicNoise(void);
		Float_t      slatResponseExp(Float_t&);
		Double_t GammaRandom();


		Int_t        fillEvent();
		Int_t        bookHistograms();
		Int_t        writeHistograms();
		Int_t        ResetFlags();


	public:
		StBTofSimMaker(const char *name="TofSim");


		virtual ~StBTofSimMaker();

		void           Reset();
		virtual Int_t  Init();
		Int_t          InitRun(Int_t);
		Int_t          FinishRun(Int_t);
		virtual Int_t  Make();
		virtual Int_t  Finish();

		StTofSimParam*    GetSimParam()       const { return mSimDb; }
		StBTofCollection*  GetBTofCollection()  const { return mBTofCollection; }
		StMcBTofHitCollection* GetMcBTofHitCollection() const { return mMcBTofHitCollection; }

		void   setCellXtalk(Bool_t val) { mCellXtalk = val; }
		string   setHistFileName();
		void   setBookHist(Bool_t val) { mBookHisto = val; }
		void   writeStEvent(Bool_t val = kTRUE) {mWriteStEvent = val;}

		virtual const char *GetCVS() const
		{static const char cvs[]="Tag $Name:  $ $Id: StBTofSimMaker.h,v 1.6 2015/07/28 22:49:55 smirnovd Exp $ built " __DATE__ " " __TIME__ ; return cvs;}

		ClassDef(StBTofSimMaker,1)
};
#endif
