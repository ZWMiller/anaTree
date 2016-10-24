/* 
 * File:   prescales.h
 * Author: mustafa
 *
 * Created on November 7, 2011, 10:10 PM
 */

#ifndef PRESCALES_H
#define	PRESCALES_H

#include "TNamed.h"

#include <vector>
#include <map>

class TH1F;

class prescales : public TNamed
{
private:
    typedef std::vector<float> vecPrescales;
    std::map<int,vecPrescales> mTable;

     std::map<int,vecPrescales>::iterator mLastQuery;
    
    static prescales* mInstance;

    void ReadList(int trg);

protected:
    prescales();
public:
    float GetPrescale(int run,int trg);
    int GetNumberOfRuns();
    int* GetListOfRuns();
    int RunIndex(int run);
    bool RunExists(int run);
    void FillPrescalesHist(TH1F*,int trg);

    static prescales* Instance();
    virtual ~prescales();

    ClassDef(prescales,1)
};

#endif	/* PRESCALES_H */

