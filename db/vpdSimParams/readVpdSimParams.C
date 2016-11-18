// Macro to read Vpd Simulation Parameters from Calibrations/tof/vpdSimParams database and write to file.
//
// based on
//  http://www.star.bnl.gov/STAR/comp/db/StoreDbTable.cc.html
//
// Nickolas Luttrell   09/13/2016
//

#include <iostream>
#include <fstream>
//#include <algorithm>
//#include "readVpdSimParams.h"
using namespace std;

void readVpdSimParams(const char* time = "2016-02-03 00:00:02")
{
  
  //-- load dBase and Table definition libraries
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("StUtilities");
  gSystem->Load("St_Tables.so");

  gSystem->Load("StDbLib.so");
  gSystem->Load("libStDb_Tables.so");

  //-- get the singleton manager
  StDbManager* dbManager = StDbManager::Instance();

  //-- connect to the db & get an empty container
  StDbConfigNode* configNode = dbManager->initConfig("Calibrations_tof");

  // set the request time
  string readTime = time;
  dbManager->setRequestTime( readTime.c_str() );

  // get the db table
  StDbTable* simParams = configNode->addDbTable( "vpdSimParams" );
  dbManager->fetchDbTable( simParams );

  // output some table details
  cout<<simParams->getVersion()<<endl;
  cout<<simParams->getBeginDateTime()<<endl;
  cout<<simParams->getEndDateTime()<<endl;


  vpdSimParams_st* table = static_cast<vpdSimParams_st*>(simParams->GetTable());
  
  if( !table ) {
    cout << " Table is invalid, exiting! " << endl;
    return;
  }

      Int_t nRows = simParams->GetNRows();
  
      cout << " Number of rows = " << nRows << endl;

  cout<<"<---------------- Read out from DataBase -------------->"<<endl;

  ofstream outData;
  outData.open("vpdSimParams_read.dat");
    
    const Int_t MAX_DB_INDEX = 38;

    cout << "Preparing to read out table data..." << endl;
    
  for ( int i = 0; i < MAX_DB_INDEX; i++ ){
    outData << table.tubeID[i] << "\t" << table.tubeRes[i] << "\t" << table.tubeStatusFlag[i] << "\t" << table.tubeTriggerFlag[i] << endl;
  }
  
  outData.close();

  cout << " Vpd Simulation Parameters written to vpdSimParams_read.dat" << endl;

}