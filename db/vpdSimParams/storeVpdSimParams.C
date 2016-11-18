// Macro to store Vpd Simulation Parameters table to Calibrations/tof/vpdSimParams database
//
// based on
//  http://www.star.bnl.gov/STAR/comp/db/StoreDbTable.cc.html
//
// Update by Nickolas Luttrell 09/19/2016
//
// Daniel Brandenburg   06/19/2014 - original version for use with tof
//

#include <iostream>
#include <string>
#include <fstream>
using namespace std;

void storeVpdSimParams( bool test = false ) {
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
    StDbTable* simParams = configNode->addDbTable( "vpdSimParams" );
    
    const Int_t MAX_DB_INDEX = 38;
    
    ifstream inData;
    inData.open("vpdSimParams.dat");
    
    vpdSimParams_st vpdMap;
    
    // set the request time
    string writeTime = time;
    string time;
    inData >> writeTime;
    inData >> time;
    writeTime = writeTime + " " + time;
    
    cout << "The time is: " << writeTime << endl;
    
    int vpdTubeRes_int;
    int vpdTubeId_int;
    int vpdTubeStatusFlag_int;
    int vpdTubeTriggerFlag_int;
    
    for ( int i = 0; i < MAX_DB_INDEX; i++ ){
        inData >> vpdTubeId_int >> vpdTubeRes_int >> vpdTubeStatusFlag_int >> vpdTubeTriggerFlag_int;
        
        vpdMap.tubeRes[ i ] = (unsigned char)vpdTubeId_int;
        vpdMap.tubeID[ i ] = (short)vpdTubeRes_int;
        vpdMap.tubeStatusFlag[ i ] = (unsigned char)vpdTubeStatusFlag_int;
        vpdMap.tubeTriggerFlag[ i ] = (unsigned char)vpdTubeTriggerFlag_int;
        
        cout << "tubeRes is: " << (int)vpdMap.tubeRes[ i ] << endl;
        cout << "tubeId is: " << (int)vpdMap.tubeID[ i ] << endl;
        cout << "tubeStatusFlag is: " << (int)vpdMap.tubeStatusFlag[ i ] << endl;
        cout << "tubeTriggerFlag is: " << (int)vpdMap.tubeTriggerFlag[ i ] << endl;
        cout << "_________________________" << endl;
    }
    
    inData.close();
    
    cout << "Setting store time to : " << writeTime << endl;
    dbManager->setStoreTime( writeTime.c_str() );
    
    cout << "Uploading VPD Simulation Parameters map" << endl;
    simParams->SetTable( (char*) &vpdMap, MAX_DB_INDEX );
    
    if ( false == test ){
        dbManager->storeDbTable( simParams );
        cout << "vpdSimParams Uploaded to DB" << endl;
    } else {
        cout << "TEST mode : Not storing to DB" << endl;
    }
    
}