#ifndef RCMPHELPER_HH
#define RCMPHELPER_HH

#include "/home/emile/GRSISort/include/TGRSIHelper.h"

#include "/home/emile/GRSISort/GRSIData/include/TRcmp.h"
#include "/home/emile/GRSISort/GRSIData/include/TGriffin.h"
#include "/home/emile/GRSISort/GRSIData/include/TGriffinBgo.h"

class RCMPHelper : public TGRSIHelper, public ROOT::Detail::RDF::RActionImpl<RCMPHelper> {
public:

   RCMPHelper(TList* list) : TGRSIHelper(list)
   {
      Prefix("RCMP");
      Setup();
   }

   ROOT::RDF::RResultPtr<std::map<std::string, TList>> Book(ROOT::RDataFrame* d) override
   {
      return d->Book<TRcmp, TGriffin, TGriffinBgo>(std::move(*this), {"TRcmp", "TGriffin", "TGriffinBgo"});
   }

   void CreateHistograms(unsigned int slot) override;

   void Exec(unsigned int slot, TRcmp& rcmp, TGriffin& griffin, TGriffinBgo& griffinbgo);

   void EndOfSort(std::shared_ptr<std::map<std::string, TList>>& list) override;

private:
   std::vector<int> frontstrip;
   std::vector<double> frontenergy;
   std::vector<int> backstrip;
   std::vector<double> backenergy;
   std::vector<double> frontdet;
   std::vector<double> backdet;

   std::vector<double> backaddbackvector;
   std::vector<double> frontaddbackvector;

   int firststripfront;
   int firststripback;
   double firstenergyfront;
   double firstenergyback;
   double addbackenergyfront;
   double addbackenergyback;
};

extern "C" RCMPHelper* CreateHelper(TList* list) { return new RCMPHelper(list); }

extern "C" void DestroyHelper(TGRSIHelper* helper) { delete helper; }

#endif
