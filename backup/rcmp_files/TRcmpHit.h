#ifndef RCMPHIT_H
#define RCMPHIT_H

#include <cstdio>
#include <cmath>
#include <vector>

#include "TVector3.h"

#include "TFragment.h"
#include "TChannel.h"

#include "TDetectorHit.h"

class TRcmpHit : public TDetectorHit {
public:
   TRcmpHit() = default;        
   TRcmpHit(const TRcmpHit&);   
   TRcmpHit(TRcmpHit&&) noexcept            = default;
   TRcmpHit& operator=(const TRcmpHit&)     = default;
   TRcmpHit& operator=(TRcmpHit&&) noexcept = default;

   explicit TRcmpHit(const TFragment& frag);                          
   ~TRcmpHit() override = default; 
   void Clear(Option_t* opt = "") override;         
   void Print(Option_t* opt = "") const override;   
   void Print(std::ostream& out) const override;    
   void Copy(TObject&) const override;              
   void Copy(TObject&, bool) const override;        
   
   ClassDefOverride(TRcmpHit, 2);};

#endif
