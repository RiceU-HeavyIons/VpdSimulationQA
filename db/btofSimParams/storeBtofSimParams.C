// Macro to store BTof Simulation Parameters table to Calibrations/tof/tofSimPars database
//
// based on
//  http://www.star.bnl.gov/STAR/comp/db/StoreDbTable.cc.html
//
// Nickolas Luttrell 11/14/2016
//

#include <iostream>
#include <string>
#include <fstream>
using namespace std;

void storeBtofSimParams( bool test = true ) {
    //-- load dBase and Table definition libraries
    gSystem->Load("St_base");
    gROOT->Macro("LoadLogger.C");
    gSystem->Load("StUtilities");
    gSystem->Load("St_Tables.so");
    
    gSystem->Load("StDbLib.so");
    gSystem->Load("libStDb_Tables.so");
    
    //-- get the singleton manager
    StDbManager* dbManager = StDbManager::Instance();
    
    //-- connect to the db & get an empty container
    StDbConfigNode* configNode = dbManager->initConfig("Calibrations_tof");
    
    // get the db table
    StDbTable* simParams = configNode->addDbTable( "tofSimPars" );
    
    const Int_t MAX_DB_INDEX = 4;
    
    ifstream inData;
    inData.open("tofSimPars.dat");
    
    tofSimPars_st btofMap;
    
    // set the request time
    string writeTime = time;
    string time;
    inData >> writeTime;
    inData >> time;
    writeTime = writeTime + " " + time;
    
    cout << "The time is: " << writeTime << endl;
    
    float par_nClusters[4];
    float er;
    float d_inner;
    float d_outer;
    float d_gap;
    float nGaps;
    float alpha;
    float nClustersMax;
    float vDriftMean;
    float vDriftErr;
    float nEperCluster;
    float nEMax;
    float vLight;
    float qElectron;
    float Adcthreshold;
    float timeRes;
    float nTbins;
    float TbinWidth;
    float T0;
    float amp;
    float AdcBinWidth;
    float TdcBinWidth;
    
    ///////////////////////////////////////////////////////
    
    inData >> par_nClusters[0] >> par_nClusters[1] >> par_nClusters[2] >> par_nClusters[3];
    for ( int i=0; i<MAX_DB_INDEX; i++) {
        btofMap.par_nClusters[ i ] = par_nClusters[i];
    }
    
    inData >> er;
    btofMap.er = er;
    inData >> d_inner;
    btofMap.d_inner = d_inner;
    inData >> d_outer;
    btofMap.d_outer = d_outer;
    inData >> d_gap;
    btofMap.d_gap = d_gap;
    inData >> nGaps;
    btofMap.nGaps = nGaps;
    inData >> alpha;
    btofMap.alpha = alpha;
    inData >> nClustersMax;
    btofMap.nClustersMax = nClustersMax;
    inData >> vDriftMean;
    btofMap.vDriftMean = vDriftMean;
    inData >> vDriftErr;
    btofMap.vDriftErr = vDriftErr;
    inData >> nEperCluster;
    btofMap.nEperCluster = nEperCluster;
    inData >> nEMax;
    btofMap.nEMax = nEMax;
    inData >> vLight;
    btofMap.vLight = vLight;
    inData >> qElectron;
    btofMap.qElectron = qElectron;
    inData >> Adcthreshold;
    btofMap.Adcthreshold = Adcthreshold;
    inData >> timeRes;
    btofMap.timeRes = timeRes;
    inData >> nTbins;
    btofMap.nTbins = nTbins;
    inData >> TbinWidth;
    btofMap.TbinWidth = TbinWidth;
    inData >> T0;
    btofMap.T0 = T0;
    inData >> amp;
    btofMap.amp = amp;
    inData >> AdcBinWidth;
    btofMap.AdcBinWidth = AdcBinWidth;
    inData >> TdcBinWidth;
    btofMap.TdcBinWidth = TdcBinWidth;
    
    ///////////////////////////////////////////////////////
    
    inData.close();
    
    cout << "Setting store time to : " << writeTime << endl;
    dbManager->setStoreTime( writeTime.c_str() );
    
    cout << "Uploading BTof Simulation Parameters map" << endl;
    simParams->SetTable( (char*) &btofMap, 1 );
    
    if ( false == test ){
        dbManager->storeDbTable( simParams );
        cout << "vpdSimParams Uploaded to DB" << endl;
    } else {
        cout << "TEST mode : Not storing to DB" << endl;
    }
    
}