// Macro to read BTof Simulation Parameters from Calibrations/tof/tofSimPars database and write to file.
//
// based on
//  http://www.star.bnl.gov/STAR/comp/db/StoreDbTable.cc.html
//
// Nickolas Luttrell   11/14/2016
//

#include <iostream>
#include <fstream>

using namespace std;

void readBtofSimParams(const char* time = "2016-02-03 00:00:02")
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
  StDbTable* simParams = configNode->addDbTable( "tofSimPars" );
  dbManager->fetchDbTable( simParams );

  // output some table details
  cout<<simParams->getVersion()<<endl;
  cout<<simParams->getBeginDateTime()<<endl;
  cout<<simParams->getEndDateTime()<<endl;


  tofSimPars_st* table = static_cast<tofSimPars_st*>(simParams->GetTable());
  
  if( !table ) {
    cout << " Table is invalid, exiting! " << endl;
    return;
  }

      Int_t nRows = simParams->GetNRows();
  
      cout << " Number of rows = " << nRows << endl;

  cout<<"<---------------- Read out from DataBase -------------->"<<endl;

  ofstream outData;
  outData.open("tofSimPars_read.dat");
    
    const Int_t MAX_DB_INDEX = 4;

    cout << "Preparing to read out table data..." << endl;
    
  for ( int i = 0; i < MAX_DB_INDEX; i++ ){
    outData << table.par_nClusters[i] << endl;
  }
    
    outData << table.er << endl;
    outData << table.d_inner << endl;
    outData << table.d_outer << endl;
    outData << table.d_gap << endl;
    outData << table.nGaps << endl;
    outData << table.alpha << endl;
    outData << table.nClustersMax << endl;
    outData << table.vDriftMean << endl;
    outData << table.vDriftErr << endl;
    outData << table.nEperCluster << endl;
    outData << table.nEMax << endl;
    outData << table.vLight << endl;
    outData << table.qElectron << endl;
    outData << table.Adcthreshold << endl;
    outData << table.timeRes << endl;
    outData << table.nTbins << endl;
    outData << table.TbinWidth << endl;
    outData << table.T0 << endl;
    outData << table.amp << endl;
    outData << table.AdcBinWidth << endl;
    outData << table.TdcBinWidth << endl;
  
  outData.close();

  cout << " BTof Simulation Parameters written to tofSimPars_read.dat" << endl;

}