/***************************************************************************
 *
 * $Id$
 *
 * Author: Nickolas Luttrell (Rice University), November 2016
 ***************************************************************************
 *
 * Description: StBTofMixerMaker.cpp - This maker accepts the BTof Collection
 * and checks for duplicate hits (hits in the same cell of the BTof). Upon
 * finding a duplicate, it takes the hit with the lowest time of flight. The
 * original BTof collection is replaced with the new, clean collection. This
 * can be used for embedding, as both simulation and data hits are parsed.
 *
 ***************************************************************************
 *
 * $Log$
 *
 ***************************************************************************/

#include <Stiostream.h>
#include "StBTofMixerMaker.h"
#include "StEventTypes.h"
#include "StEvent/StBTofCollection.h"
#include <map>

ClassImp(StBTofMixerMaker)

//_____________________________________________________________________________
StBTofMixerMaker::StBTofMixerMaker(const char *name):StMaker(name)
{
    //set default values
}

StBTofMixerMaker::~StBTofMixerMaker()
{
    //Destructor
}

//_____________________________________________________________________________
Int_t StBTofMixerMaker::Init()
{

    return StMaker::Init();
}


//_____________________________________________________________________________
Int_t StBTofMixerMaker::InitRun(Int_t runnumber)
{
    return kStOK;
}

//_____________________________________________________________________________
Int_t StBTofMixerMaker::FinishRun(Int_t runnumber)
{
    return kStOk;
}


//_____________________________________________________________________________
Int_t StBTofMixerMaker::Finish()
{
    return kStOK;
}


Int_t StBTofMixerMaker::Make()
{
    Int_t counter = 0;  // Number of duplicate entries discovered.
    Int_t keyId = 0;    // Unique identifier for a given cell. Scaling factors insure that no digits can be overwritten.

    mEvent = (StEvent*)GetInputDS("StEvent");
    if (!mEvent) {
        LOG_ERROR << "No StEvent! Bailing out ..." << endm;
    }

    mBTofCollection = mEvent->btofCollection();
    if(!mBTofCollection) {
        LOG_ERROR << "No BTofCollection! Bailing out ..." << endm;
    }

    LOG_DEBUG << "The original size of the collection was " << mBTofCollection->tofHits().size() << endm;

    // Iterate through the Hits to find duplicates, take hit with earliest time
    std::map<Int_t,StBTofHit*> Hits;
    std::pair<std::map<Int_t,StBTofHit*>::iterator,Bool_t> duplicateFlag;

    StBTofHit *tempHit;

    // Iterate through the BTofCollection, inserting hits into map.
    for (Int_t j=0; j < (Int_t)mBTofCollection->tofHits().size(); j++) {

        tempHit = mBTofCollection->tofHits()[j];
        keyId = (tempHit->tray() * 10000) + (tempHit->module() * 100) + (tempHit->cell());
        duplicateFlag = Hits.insert(std::pair<Int_t,StBTofHit*>(keyId,tempHit));    // Inserts into map, or if key already exists, returns False.
        if (!duplicateFlag.second) {
            counter += 1;
            StBTofHit *hit1 = Hits[keyId];  // The existing hit
            StBTofHit *hit2 = tempHit;                // The duplicate hit

            LOG_DEBUG << "The keyId and LE of the existing hit is: " << keyId << " : " << hit1->leadingEdgeTime() << endm;
            LOG_DEBUG << "The keyId and LE of the duplicate is: " << keyId << " : " << hit2->leadingEdgeTime() << endm;

            // Check which hit has the earlier LE time and replace if necessary.
            if ( hit1->leadingEdgeTime() > hit2->leadingEdgeTime() ) {
                Hits[keyId] = hit2;
            }
            LOG_DEBUG << "The value in the map is now: " << Hits[keyId]->leadingEdgeTime() << endm;
        }

    }

    LOG_DEBUG << "Duplicates found: " << counter << endm;

    mBTofCollection = new StBTofCollection();

    std::map<Int_t,StBTofHit*>::iterator item;

    // Build the new BTofCollection.
    for (item=Hits.begin(); item!=Hits.end(); ++item) {
        StBTofHit &aBTofHit = *item->second;
        mBTofCollection->addHit(new StBTofHit(aBTofHit));
    }

    //Fill StBTofHeader --

    StBTofHeader *tofHeader = 0;
    StBTofHeader aHead;

    if(!mBTofCollection) {
        LOG_DEBUG << " No StEvent/btofCollection, creating new... " << endm;
        mBTofCollection->setHeader(new StBTofHeader(aHead));
    }
    else {
        tofHeader = (StBTofHeader *) mBTofCollection->tofHeader();
    }

    LOG_DEBUG << "... Modified StBTofCollection Stored in StEvent! " << endm;

    LOG_DEBUG << "The size of the collection is now " << mBTofCollection->tofHits().size() << endm;

    return kStOK;
}


// End StBTofMixerMaker.cpp
