#ifndef RCMPHELPER_HH
#define RCMPHELPER_HH

#include "/home/emile/GRSISort/include/TGRSIHelper.h"

#include "/home/emile/GRSISort/GRSIData/include/TRcmp.h"

class RCMPHelper : public TGRSIHelper, public ROOT::Detail::RDF::RActionImpl<RCMPHelper> {
public:

   RCMPHelper(TList* list) : TGRSIHelper(list)
   {
      Prefix("RCMP");
      Setup();
   }

   ROOT::RDF::RResultPtr<std::map<std::string, TList>> Book(ROOT::RDataFrame* d) override
   {
      return d->Book<TRcmp>(std::move(*this), {"TRcmp"});
   }

   void CreateHistograms(unsigned int slot) override;

   void Exec(unsigned int slot, TRcmp& rcmp);

   void EndOfSort(std::shared_ptr<std::map<std::string, TList>>& list) override;

private:
   Long64_t fCycleLength{0};
};

extern "C" RCMPHelper* CreateHelper(TList* list) { return new RCMPHelper(list); }

extern "C" void DestroyHelper(TGRSIHelper* helper) { delete helper; }

#endif
