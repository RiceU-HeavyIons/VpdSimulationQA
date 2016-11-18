//
//  StVpdSimConfig.h
//
//
//  Created by Nickolas Luttrell on 6/21/16.
//
//

#ifndef StVpdSimConfig_h
#define StVpdSimConfig_h

#include <iostream>
#include <fstream>
#include <vector>
#include "St_db_Maker/St_db_Maker.h"
#include "tables/St_vpdSimParams_Table.h"

using std::string;
class vpdSimParams_st;


class StVpdSimConfig {
public:

		StVpdSimConfig() {}
		~StVpdSimConfig() {}

	struct SingleTubeParams{
		Float_t singleTubeRes;
		Int_t tubeId;
		Int_t tubeStatusFlag;
		Int_t tubeTriggerFlag;
	};


    // Loads Vpd Sim Params from database

	void loadVpdSimParams(const int date = 20160913, const int time = 175725, const char* Default_time = "2016-09-13 17:57:25")
	{

		St_db_Maker *dbMk=new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
		dbMk->SetDebug();
		dbMk->SetDateTime(date,time); // event or run start time, set to your liking
		dbMk->SetFlavor("ofl");

		dbMk->Init();
		dbMk->Make();

		TDataSet *DB = 0;

		DB = dbMk->GetInputDB("Calibrations/tof/vpdSimParams");
		if (!DB) {
			LOG_WARN << "Failed to connect to Database!" << endm;
			}

		St_vpdSimParams *dataset = 0;
		dataset = (St_vpdSimParams*) DB->Find("vpdSimParams");
		Int_t rows = dataset->GetNRows();
		if (rows > 1) {
			LOG_INFO << "INFO: found table with " << endm;
		}

		if (dataset) {
			TDatime val[3];
			dbMk->GetValidity((TTable*)dataset,val);
			vpdSimParams_st* table = static_cast<vpdSimParams_st*>(dataset->GetTable());

			const Int_t MAX_ARRAY_INDEX = 38;

			SingleTubeParams params;
			for (Int_t i = 0; i < MAX_ARRAY_INDEX; i++) {
				params.tubeId = table->tubeID[i];
				params.singleTubeRes = table->tubeRes[i];
				params.tubeStatusFlag = table->tubeStatusFlag[i];
				params.tubeTriggerFlag = table->tubeTriggerFlag[i];
				mSimParams[table->tubeID[i]] = params;
			}

			return;

		}
		else {
			LOG_WARN << "ERROR: dataset does not contain requested table" << endm;
			return;
		}
	}

	/* Reads VPD Sim Params from a file for DEBUG purposes
	 * TODO: add some safety for badly formed files
	 */
	void loadVpdSimParams(string params_filename ) {

		Int_t MAX_DB_INDEX = 38;

		Int_t vpdTubeRes;
		Int_t vpdTubeId;
		Int_t vpdTubeStatusFlag;
		Int_t vpdTubeTriggerFlag;

		SingleTubeParams params;

		std::ifstream inData;
		inData.open( params_filename.c_str() );

		for (int i = 0; i < MAX_DB_INDEX; i++) {
			inData >> vpdTubeId >> vpdTubeRes >> vpdTubeStatusFlag >> vpdTubeTriggerFlag;
			params.tubeId = vpdTubeId;
			params.singleTubeRes = vpdTubeRes;
			params.tubeStatusFlag = vpdTubeStatusFlag;
			params.tubeTriggerFlag = vpdTubeTriggerFlag;
			mSimParams[params.tubeId] = params;
		}

		inData.close();
		return;
	}

	std::map<Int_t, SingleTubeParams> getParams(){
		return mSimParams;
	}
    
    Float_t getThreshold()  const { return mThreshold; }
    Float_t getVpdDistance()    const { return VPDDISTANCE; }
    Float_t getTDiffCut()   const { return TDIFFCUT; }

protected:

	// stores a map of the single tube params indexed on tubeId
	std::map<Int_t, SingleTubeParams> mSimParams;
	Int_t mThreshold = 1;         // Threshold value for a tube to recognize it as a true hit.

	const Float_t VPDDISTANCE = 570;       // Distance (in cm) of each Vpd from the zero point
	const Float_t TDIFFCUT = 0.8;       // Cut value for eliminating times with a significant deviation from avg.

};

#endif /* Config_h */
