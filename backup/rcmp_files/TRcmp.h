#ifndef TRCMP_H
#define TRCMP_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <map>
#include <utility>
#include <set>

#include "TDetector.h"
#include "TChannel.h"
#include "TRcmpHit.h"
#include "TMath.h"
#include "TVector3.h"
#include "Globals.h"
#include "TFragment.h"
#include "TGRSIMnemonic.h"

class TRcmp : public TDetector {
    public:
        TRcmp();               
        TRcmp(const TRcmp&);   
        TRcmp(TRcmp&&) noexcept = default;
        TRcmp& operator=(const TRcmp&);
        TRcmp& operator=(TRcmp&&) noexcept = default;
        ~TRcmp() override;   

#ifndef __CINT__

        void AddFragment(const std::shared_ptr<const TFragment>&, TChannel*) override;   
#endif
        void BuildHits() override {}     

        TRcmpHit* GetRcmpHit(const int& i);         

        void Copy(TObject&) const override;              
        void Clear(Option_t* opt = "all") override;      
        void Print(Option_t* opt = "") const override;   
        void Print(std::ostream& out) const override;    

        ClassDefOverride(TRcmp, 2) 
};
#endif
