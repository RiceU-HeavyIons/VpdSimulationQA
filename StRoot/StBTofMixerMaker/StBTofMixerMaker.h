//
//  StBTofMixerMaker.h
//  
//
//  Created by Nickolas Luttrell on 6/3/16.
//
//

#ifndef StBTofMixerMaker_HH
#define StBTofMixerMaker_HH

#include "StMaker.h"
#include "St_DataSet.h"

class StEvent;
class StBTofCollection;

#include <vector>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class StBTofMixerMaker : public StMaker{
    protected:
    
        // All the internal definitions should be placed here.
    
        StEvent           *mEvent;
        StBTofCollection   *mBTofCollection = nullptr;
    
    public:
    StBTofMixerMaker(const char *name="Simulation");
    
    virtual ~StBTofMixerMaker();
    
    virtual Int_t  Init();
    Int_t          InitRun(Int_t);
    Int_t          FinishRun(Int_t);
    virtual Int_t  Make();
    virtual Int_t  Finish();

    virtual const char *GetCVS() const
    {static const char cvs[]="Tag $Name:  $ $Id: StBTofMixerMaker.h,v 1.0 2016/06/03 00:00:00 nluttrel Exp $ built " __DATE__ " " __TIME__ ; return cvs;}
    
    ClassDef(StBTofMixerMaker,2)
    
};

#endif /* StBTofMixerMaker_h */



// end of StBTofMixerMaker.h