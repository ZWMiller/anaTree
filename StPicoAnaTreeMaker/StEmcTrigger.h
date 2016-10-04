#ifndef StEmcTrigger_hh
#define StEmcTrigger_hh

#include "TObject.h"
#include "stdio.h"


class StEmcTrigger : public TObject {
 public:
  StEmcTrigger();
  ~StEmcTrigger();
  StEmcTrigger(int, int, int, int, int);
  void    Clear(const Option_t *opt="");
  virtual void Print(const Char_t *option = "") const;  ///< Print trigger info
 
  Int_t   flag() const           { return (Int_t)mFlag; }
  Int_t   id() const             { return (Int_t)mId; }
  Int_t   adc() const            { return (Int_t)mAdc; }
  Int_t   eId() const            { return (Int_t)mEId; }
  Int_t   adc0() const            { return (Int_t)mAdc0; }
  
  void    setEId(int eid)        { mEId = eid; }
  void    setAdc0(int adc0)        { mAdc0 = adc0; }

 protected:

  UChar_t mFlag;   // 0x1: ht0, 0x2: ht1, 0x4: ht2; 0x8: ht3
  UShort_t mId;    // soft id.  bjp: 1-18, ht: 1-4800
  UShort_t mAdc;   // adc
  UShort_t mEId;   // id to matched electron
  UShort_t mAdc0;   // adc0, matched adc0
  

  friend class StPicoDst;

  ClassDef(StEmcTrigger, 1)
};

#endif
