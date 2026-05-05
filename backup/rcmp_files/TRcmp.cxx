#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "TRcmp.h"

TRcmp::TRcmp()
{
    Clear();
}

TRcmp::~TRcmp()
{
}

TRcmp& TRcmp::operator=(const TRcmp& rhs)
{
    rhs.Copy(*this);
    return *this;
}

TRcmp::TRcmp(const TRcmp& rhs) : TDetector(rhs)
{
    rhs.Copy(*this);
}

void TRcmp::Copy(TObject& rhs) const
{
    TDetector::Copy(rhs);   
}

void TRcmp::Clear(Option_t* opt)
{
    TDetector::Clear(opt);
}

void TRcmp::AddFragment(const std::shared_ptr<const TFragment>& frag, TChannel* chan)
{
    if(frag == nullptr || chan == nullptr) return; 
         
    auto* hit = new TRcmpHit(*frag);

    Hits().push_back(hit);
}

TRcmpHit* TRcmp::GetRcmpHit(const int& i)
{
    try { return static_cast<TRcmpHit*>(Hits().at(i)); }
    
    catch(const std::out_of_range& oor) 
    {
        std::cerr << ClassName() << " is out of range: " << oor.what() << std::endl;

        throw grsi::exit_exception(1);
    }

    return nullptr;
}

void TRcmp::Print(Option_t*) const
{
    Print(std::cout);
}

void TRcmp::Print(std::ostream& out) const
{
    std::ostringstream str;
    out << str.str();
}
