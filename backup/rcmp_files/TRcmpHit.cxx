#include <iostream>
#include "TRcmpHit.h"

#include <cstdlib>
#include <cmath>

TRcmpHit::TRcmpHit(const TRcmpHit& rhs) : TDetectorHit(rhs)
{
    rhs.Copy(*this);
}

TRcmpHit::TRcmpHit(const TFragment& frag)
{
    frag.Copy(*this);
}

void TRcmpHit::Copy(TObject& rhs) const
{
    TDetectorHit::Copy(rhs);
}

void TRcmpHit::Copy(TObject& rhs, bool waveform) const
{
    Copy(rhs);

    if(waveform) 
    {
        CopyWave(rhs);
    }
}

void TRcmpHit::Clear(Option_t* opt)
{
    TDetectorHit::Clear(opt);   
}

void TRcmpHit::Print(Option_t*) const
{
    Print(std::cout);
}

void TRcmpHit::Print(std::ostream& out) const
{
    std::ostringstream str;
    out << str.str();
}
