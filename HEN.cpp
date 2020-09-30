#include <bits/stdc++.h>

using namespace std;

struct stream
{
 string name; int id;
 double Ti, Tf, CP;
};

struct HeatExchanger
{
 stream H, C;
 double q;
};

bool comp(struct stream a, struct stream b)
{
 if(a.CP>b.CP) return true;
 else return false;
}

bool TempFinalSmaller(struct stream a, struct stream b)
{
 if(a.Tf<b.Tf) return true;
 else if(a.Tf>b.Tf) return false;
 else if(a.CP>b.CP) return true;
 else return false;
}

bool TempInitialSmaller(struct stream a, struct stream b)
{
 if(a.Ti<b.Ti) return true;
 else if(a.Ti>b.Ti) return false;
 else if(a.CP>b.CP) return true;
 else return false;
}

bool TempFinalGreater(struct stream a, struct stream b)
{
 return (not TempFinalSmaller(a, b));
}

bool TempInitialGreater(struct stream a, struct stream b)
{
 return (not TempInitialSmaller(a, b));
}

//is there intersection in set (x1, x2) & (x3, x4)?
bool intersection(double x1, double x2, double x3, double x4)
{
 if(max(x1, x3)<min(x2, x4)) return true;
 else return false;
}

void reverseVector(vector<double> &v)
{
 double temp;
 for(int i=0; i<v.size()/2; i++)
 {
  temp = v[i];
  v[i] = v[v.size()-1-i];
  v[v.size()-1-i] = temp;
 }
}


void removeDuplicates(vector<double> &A)
{
 if(A.size()<2) return;
 int count=0;
 for(int i=1; i<A.size(); i++)
    if(A[i]!=A[count]) A[++count] = A[i];
 A.resize(count+1);
 return;
}

void printVector(vector<double> &A, string name)
{
 cout << name+":  ";
 for(int i=0; i<A.size(); i++)
    cout << A[i] << "  ";
 cout << endl << endl;
}

void printStreams(vector<stream> &A, string name)
{
 cout << name+":" << endl << endl;
 cout << "Ti          Tf          CP          name" << endl;
 for(int i=0; i<A.size(); i++)
    cout << A[i].Ti << "          " << A[i].Tf << "          " << A[i].CP << "          " << A[i].name << endl;
 cout << endl;
}

void printHeatExchangers(vector<HeatExchanger> A, string name)
{
 cout << name+":" << endl << endl;
 cout << "TiH          TfH           CPh          TiC          TfC          CPc          HSname          CSname           HeatDuty" << endl;
 for(int i=0; i<A.size(); i++)
  {
    cout << A[i].H.Ti << "          " << A[i].H.Tf << "          " << A[i].H.CP << "          ";
    cout << A[i].C.Ti << "          " << A[i].C.Tf << "          " <<  A[i].C.CP << "          ";
    cout << A[i].H.name << "          " << A[i].C.name << "          " <<  A[i].q << "          " << endl;
  }
 cout << endl;
}

void giveStreamsAbovePinch(vector<stream> &HS, vector<stream> &CS, vector<stream> &HStemp, vector<stream> &CStemp, 
                                                 vector<stream> &HSnp, vector<stream> &CSnp, double Tci, double Tcf, double Thi, double Thf)
{
  for(int i=0; i<HS.size(); i++)
 {
   if(intersection(Thf, Thi, HS[i].Tf, HS[i].Ti))
   {
    stream temp = HS[i];
    temp.Tf = max(temp.Tf, Thf);
    temp.Ti = min(temp.Ti, Thi);
    if(abs(temp.Tf-Thf)<0.01) HStemp.push_back(temp);
    else HSnp.push_back(temp);
   }
 }
 for(int i=0; i<CS.size(); i++)
 {
   if(intersection(Tci, Tcf, CS[i].Ti, CS[i].Tf))
   {
    stream temp = CS[i];
    temp.Ti = max(temp.Ti, Tci);
    temp.Tf = min(temp.Tf, Tcf); 
    if(abs(temp.Ti-Tci)<0.01) CStemp.push_back(temp);
    else if(CS[i].Ti>Tci) CSnp.push_back(temp);
   }
 }
 
 sort(HSnp.begin(), HSnp.end(), TempFinalSmaller);
 sort(CSnp.begin(), CSnp.end(), TempInitialSmaller);
}

void giveStreamsBelowPinch(vector<stream> &HS, vector<stream> &CS, vector<stream> &HStemp, vector<stream> &CStemp, 
                                                 vector<stream> &HSnp, vector<stream> &CSnp, double Tci, double Tcf, double Thi, double Thf)
{
  for(int i=0; i<HS.size(); i++)
 {
   if(intersection(Thf, Thi, HS[i].Tf, HS[i].Ti))
   {
    stream temp = HS[i];
    temp.Tf = max(temp.Tf, Thf);
    temp.Ti = min(temp.Ti, Thi);
    if(abs(temp.Ti-Thi)<0.01) HStemp.push_back(temp);
    else HSnp.push_back(temp);
   }
 }
 for(int i=0; i<CS.size(); i++)
 {
   if(intersection(Tci, Tcf, CS[i].Ti, CS[i].Tf))
   {
    stream temp = CS[i];
    temp.Ti = max(temp.Ti, Tci);
    temp.Tf = min(temp.Tf, Tcf); 
    if(abs(temp.Tf-Tcf)<0.01) CStemp.push_back(temp);
    else if(CS[i].Ti>Tci) CSnp.push_back(temp);
   }
 }
 
 sort(HSnp.begin(), HSnp.end(), TempInitialGreater);
 sort(CSnp.begin(), CSnp.end(), TempFinalGreater);
}

void findTempDiff(vector<double> &CPh, vector<double> &CPc, vector<double> &T, int p, 
                   double dTmin, vector<double> &nearPinchesC, vector<double> &nearPinchesH)
{
 vector<double> Tc, Th;
 double tempC = T[p], tempH = T[p];
 int i=p, j=p, n = T.size(); 
 Tc.push_back(T[p]-(dTmin/2.0));
 Th.push_back(T[p]+(dTmin/2.0));
 while(i<n-1 && j<n-1)
 {
  if(CPh[i]==0.0)
  { i++; tempH = T[i]; continue; }
  if(CPc[j]==0.0)
  { j++;tempC = T[j]; continue; }

  if(CPh[i]*(tempH-T[i+1]) < CPc[j]*(tempC-T[j+1]))
  {
    tempC = tempC - (CPh[i]*(tempH-T[i+1]))/CPc[j];
    tempH = T[i+1];
    i++;
  }
  else if(CPh[i]*(tempH-T[i+1]) > CPc[j]*(tempC-T[j+1]))
  {
    tempH = tempH - (CPc[j]*(tempC-T[j+1]))/CPh[i];
    tempC = T[j+1];
    j++;
  }
  else 
  {
   tempH = T[i+1]; 
   tempC = T[j+1];
   i++; j++;
  }
  Tc.push_back(tempC-(dTmin/2.0));
  Th.push_back(tempH+(dTmin/2.0));
 }
 
 reverseVector(Tc);
 reverseVector(Th);
 
 i=p-1; j=p-1;
 tempC = T[p]; tempH = T[p];
 while(i>=0 && j>=0)
 {
  if(CPh[i]==0.0)
  { i--; tempH = T[i+1]; continue; }
  if(CPc[j]==0.0)
  { j--; tempC = T[j+1]; continue; }

  if(CPh[i]*(T[i]-tempH) < CPc[j]*(T[j]-tempC))
  {
    tempC = tempC + (CPh[i]*(T[i]-tempH))/CPc[j];
    tempH = T[i];
    i--;
  }
  else if(CPh[i]*(T[i]-tempH) > CPc[j]*(T[j]-tempC))
  {
    tempH = tempH + (CPc[j]*(T[j]-tempC))/CPh[i];
    tempC = T[j];
    j--;
  }
  else 
  {
   tempH = T[i]; 
   tempC = T[j];
   i--; j--;
  }
  Tc.push_back(tempC-(dTmin/2.0));
  Th.push_back(tempH+(dTmin/2.0));
 }

 //printVector(Th, "Hot Temp");
 //printVector(Tc, "cold Temp");
 
 for(int i=0; i<Th.size(); i++) 
 {
  if(Th.size()<=1) break;
  if(i==0 && Th[0]-Tc[0]<Th[1]-Tc[1] && (Th[0]-Tc[0])/dTmin <1.5) {nearPinchesC.push_back(Tc[i]); nearPinchesH.push_back(Th[i]);}
  else if(i==Th.size()-1 && Th[i]-Tc[i]<Th[i-1]-Tc[i-1] && (Th[i]-Tc[i])/dTmin <1.5) 
                                                                   {nearPinchesC.push_back(Tc[i]); nearPinchesH.push_back(Th[i]);}
  else if(Th[i]-Tc[i]<Th[i+1]-Tc[i+1] && Th[i]-Tc[i]<Th[i-1]-Tc[i-1] && (Th[i]-Tc[i])/dTmin < 1.5) 
                                                                   {nearPinchesC.push_back(Tc[i]); nearPinchesH.push_back(Th[i]);}
 }
}

bool fitMatch(vector<stream> &A, vector<stream> &B)
{
 for(int i=0; i<A.size(); i++) 
  if(A[i].CP>B[i].CP) return false;
 return true;
}

HeatExchanger matchColdstream(vector<stream> &HS, int i, vector<stream> &CS, int j, double dTmin, double stopT)
{
 if(HS[i].Tf < CS[j].Ti + dTmin)
 {
  HeatExchanger result;
  result.q = 0.0;
  return result;
 }

 if(HS[i].CP <= CS[j].CP)
 {
  HeatExchanger result;
  result.q = min(HS[i].CP*(HS[i].Ti-HS[i].Tf), CS[j].CP*(stopT-CS[j].Ti)); 
  result.H = HS[i];
  result.C = CS[j];
  result.H.Ti = result.H.Tf + (result.q/result.H.CP);
  result.C.Tf = result.C.Ti + (result.q/result.C.CP);
  HS[i].Tf = result.H.Ti;
  CS[j].Ti = result.C.Tf;
  return result;
 }
 
 else
 {
  double temp = (HS[i].CP*(HS[i].Tf-dTmin)-CS[j].CP*CS[j].Ti)/(HS[i].CP-CS[j].CP);
  HeatExchanger result;
  result.q = min(min(HS[i].CP*(HS[i].Ti-HS[i].Tf), CS[j].CP*(stopT-CS[j].Ti)), CS[j].CP*(temp-CS[j].Ti));
  result.H = HS[i];
  result.C = CS[j];
  result.H.Ti = result.H.Tf + (result.q/result.H.CP);
  result.C.Tf = result.C.Ti + (result.q/result.C.CP);
  HS[i].Tf = result.H.Ti;
  CS[j].Ti = result.C.Tf;
  return result;
 }
}

HeatExchanger matchHotstream(vector<stream> &HS, int i, vector<stream> &CS, int j, double dTmin, double stopT)
{
 if(HS[i].Ti < CS[j].Tf + dTmin)
 {
  HeatExchanger result;
  result.q = 0.0;
  return result;
 }

 if(HS[i].CP >= CS[j].CP)
 {
  HeatExchanger result;
  result.q = min(HS[i].CP*(HS[i].Ti-stopT), CS[j].CP*(CS[j].Tf-CS[j].Ti)); 
  result.H = HS[i];
  result.C = CS[j];
  result.H.Tf = result.H.Ti - (result.q/result.H.CP);
  result.C.Ti = result.C.Tf - (result.q/result.C.CP);
  HS[i].Ti = result.H.Tf;
  CS[j].Tf = result.C.Ti;
  return result;
 }
 
 else
 {
  double temp = (CS[j].CP*CS[j].Tf-HS[i].CP*(HS[i].Ti-dTmin))/(CS[j].CP-HS[i].CP);
  HeatExchanger result;
  result.q = min(min(CS[j].CP*(CS[j].Tf-CS[j].Ti), HS[i].CP*(HS[i].Ti-stopT)), CS[j].CP*(CS[j].Tf-temp));
  result.H = HS[i];
  result.C = CS[j];
  result.H.Tf = result.H.Ti - (result.q/result.H.CP);
  result.C.Ti = result.C.Tf - (result.q/result.C.CP);
  HS[i].Ti = result.H.Tf;
  CS[j].Tf = result.C.Ti;
  return result;
 }
}

void matchAllAP(vector<stream> &HS, vector<stream> &CS, vector<stream> &HSnp, 
                                 vector<stream> &CSnp, double dTmin, vector<HeatExchanger> &result)
{
 if(HS.size()+HSnp.size() == 0|| CS.size()+CSnp.size() == 0) return;

 //SPLIT COLDSTREAM
 if(HS.size()>CS.size())
 {
  CS.insert(CS.begin(), CS[0]);
  CS[0].name += "-s1";
  CS[1].name += "-s2";
  CS[0].CP *= HS[0].CP/(HS[0].CP+HS[1].CP);
  CS[1].CP *= HS[1].CP/(HS[0].CP+HS[1].CP);
  sort(CS.begin(), CS.end(), comp);
  matchAllAP(HS, CS, HSnp, CSnp, dTmin, result);
 }
 
 //SPLIT HOTSTREAM
 if(fitMatch(HS, CS) == false)
 {
  HS.insert(HS.begin(), HS[0]);
  HS[0].name += "-s1";
  HS[1].name += "-s2";
  HS[0].CP /= 2.0;
  HS[1].CP /= 2.0;
  sort(HS.begin(), HS.end(), comp);
  matchAllAP(HS, CS, HSnp, CSnp, dTmin, result);
 }

 double stopT, sumHcp=0.0, sumCcp=0.0;
 for(int i=0; i<CS.size(); i++)  sumCcp+=CS[i].CP;
 for(int i=0; i<HS.size(); i++)  sumHcp+=HS[i].CP;

 while(HS.size()+HSnp.size() > 0)
 {
  if(HS.size() == 0){vector<stream> k(0); matchAllAP(HSnp, CS , k, CSnp, dTmin, result);}
  stopT = 100000.0;
  for(int i=0; i<HSnp.size(); i++) stopT = min(stopT, HSnp[i].Tf-dTmin);
  for(int i=0; i<HS.size(); i++) stopT = min(stopT, CS[i].Tf);

  int i=0, j=0;
  while(i<HS.size())
  {
   result.push_back(matchColdstream(HS, i, CS, j, dTmin, stopT));
   cout << "Match" + result[result.size()-1].H.name + "(" << result[result.size()-1].H.CP << ") and " + result[result.size()-1].C.name+"(" <<
	result[result.size()-1].C.CP << ") with heat duty = " << result[result.size()-1].q;
  cout << "Hotstream " << result[result.size()-1].H.Ti << " -> " << result[result.size()-1].H.Tf << " and " << 
          "Coldstream" << result[result.size()-1].C.Ti << " -> " << result[result.size()-1].C.Tf << endl;
   if(abs(HS[i].Ti-HS[i].Tf)<0.00001)
   {  sumHcp -= HS[i].CP; HS.erase(HS.begin()+i); i--;}
   if(abs(CS[j].Ti-CS[j].Tf)<0.00001)
   {  sumCcp -= CS[j].CP; CS.erase(CS.begin()+j); j--;}
   i++; j++;
  }

  while(HSnp.size()>0)
  {
   if(HSnp[0].Tf <=stopT+dTmin)
    {
     HS.push_back(HSnp[0]);
     sumHcp += HSnp[0].CP;
     HSnp.erase(HSnp.begin());
    }
   else break;
  }
  
  while(CSnp.size()>0)
  {
   if(CSnp[0].Ti <= stopT)
    {
     CS.push_back(CSnp[0]);
     sumCcp += CSnp[0].CP;
     CSnp.erase(CSnp.begin());
    }
   else break;
  }
  
  if(sumCcp<sumHcp)
  {
   for(i=0; i<HSnp.size(); i++) {HS.push_back(HSnp[i]); sumHcp+=HSnp[i].CP;}
   HSnp.clear();
   sort(HS.begin(), HS.end(), TempInitialSmaller);
   break;
  }
 }
 
 while(HS.size()!=0)
 {
  if(HS.size() == 1) stopT = 100000.0;
  else stopT = HS[1].Tf-dTmin;
  int val = CS.size()-1, q = 0.0;
  for(int j=CS.size()-1; j>=0; j--)
  {
    stream h = HS[0], c = CS[j];
    HeatExchanger hc = matchColdstream(HS, 0, CS, j, dTmin, stopT);
    if(hc.q>q) {q = hc.q; val = j;}
    HS[0] = h; CS[j]=c;
  }
  result.push_back(matchColdstream(HS, 0, CS, val, dTmin, stopT));
  cout << "Match" + result[result.size()-1].H.name + "(" << result[result.size()-1].H.CP << ") and " + result[result.size()-1].C.name+"(" <<
	result[result.size()-1].C.CP << ") with heat duty = " << result[result.size()-1].q;
  cout << "Hotstream " << result[result.size()-1].H.Ti << " -> " << result[result.size()-1].H.Tf << " and " << 
          "Coldstream" << result[result.size()-1].C.Ti << " -> " << result[result.size()-1].C.Tf << endl;
  if(abs(HS[val].Ti-HS[val].Tf)<0.00001) HS.erase(HS.begin());
  if(abs(CS[val].Ti-CS[val].Tf)<0.00001) CS.erase(CS.begin());

  while(CSnp.size()!=0)
  {
   if(CSnp[0].Ti<=stopT)
     {
      CS.push_back(CSnp[0]);
      CSnp.erase(CSnp.begin());
     }
   else break;
  }
  
  sort(CS.begin(), CS.end(), comp);
  sort(HS.begin(), HS.end(), TempInitialSmaller);
 }
 return;
}

void matchAllBP(vector<stream> &HS, vector<stream> &CS, vector<stream> &HSnp, 
                                 vector<stream> &CSnp, double dTmin, vector<HeatExchanger> &result)
{
 if(HS.size()+HSnp.size() == 0|| CS.size()+CSnp.size() == 0) return;

 //SPLIT HOTSTREAM
 if(HS.size()<CS.size())
 {
  HS.insert(HS.begin(), HS[0]);
  HS[0].name += "-s1";
  HS[1].name += "-s2";
  HS[0].CP *= CS[0].CP/(CS[0].CP+CS[1].CP);
  HS[1].CP *= CS[1].CP/(CS[0].CP+CS[1].CP);
  sort(HS.begin(), HS.end(), comp);
  matchAllBP(HS, CS, HSnp, CSnp, dTmin, result);
  return;
 }
 
 //SPLIT COLDSTREAM
 if(fitMatch(CS, HS) == false)
 {
  CS.insert(CS.begin(), CS[0]);
  CS[0].name += "-s1";
  CS[1].name += "-s2";
  CS[0].CP /= 2.0;
  CS[1].CP /= 2.0;
  sort(CS.begin(), CS.end(), comp);
  matchAllBP(HS, CS, HSnp, CSnp, dTmin, result);
  return;
 }

 double stopT, sumHcp=0.0, sumCcp=0.0;
 for(int i=0; i<CS.size(); i++)  sumCcp+=CS[i].CP;
 for(int i=0; i<HS.size(); i++)  sumHcp+=HS[i].CP;

 while(CS.size()+CSnp.size() > 0)
 {
  if(CS.size() == 0){vector<stream> k(0); matchAllBP(HS, CSnp , HSnp, k, dTmin, result);}
  stopT = -100000.0;
  for(int i=0; i<CSnp.size(); i++) stopT = max(stopT, CSnp[i].Tf+dTmin);
  for(int i=0; i<CS.size(); i++) stopT = max(stopT, HS[i].Tf);

  int i=0, j=0;
  while(j<CS.size())
  {
   result.push_back(matchHotstream(HS, i, CS, j, dTmin, stopT));
   //cout << result[result.size()-1].H.name << " " << result[result.size()-1].H.Ti << " " << result[result.size()-1].H.Tf << " " << 
           //result[result.size()-1].C.name << " " << result[result.size()-1].C.Ti << " " << result[result.size()-1].C.Tf << " " << endl;

   if(abs(HS[i].Ti-HS[i].Tf)<0.00001)
   {  sumHcp -= HS[i].CP; HS.erase(HS.begin()+i); i--;}
   if(abs(CS[j].Ti-CS[j].Tf)<0.00001)
   {  sumCcp -= CS[j].CP; CS.erase(CS.begin()+j); j--;}
   i++; j++;
  }

  while(HSnp.size()>0)
  {
   if(HSnp[0].Tf >=stopT)
    {
     HS.push_back(HSnp[0]);
     sumHcp += HSnp[0].CP;
     HSnp.erase(HSnp.begin());
    }
   else break;
  }
  
  while(CSnp.size()>0)
  {
   if(CSnp[0].Tf >= stopT-dTmin)
    {
     CS.push_back(CSnp[0]);
     sumCcp += CSnp[0].CP;
     CSnp.erase(CSnp.begin());
    }
   else break;
  }
  
  if(sumCcp>sumHcp)
  {
   for(i=0; i<CSnp.size(); i++) {CS.push_back(CSnp[i]); sumCcp+=CSnp[i].CP;}
   CSnp.clear();
   sort(CS.begin(), CS.end(), TempFinalGreater);
   break;
  }
  
 }
 //printHeatExchangers(result, "Below pinch");
 while(CS.size()!=0)
 {
  if(CS.size() == 1) stopT = -100000.0;
  else stopT = CS[1].Tf+dTmin;
  int val = HS.size()-1, q = 0.0;
  for(int j=HS.size()-1; j>=0; j--)
  {
    stream c = CS[0], h = HS[j];
    HeatExchanger hc = matchHotstream(HS, j, CS, 0, dTmin, stopT);
    if(hc.q>q) {q = hc.q; val = j;}
    CS[0] = c; HS[j]=h;
  }
  result.push_back(matchHotstream(HS, val, CS, 0, dTmin, stopT));
  //cout << result[result.size()-1].H.name << " " << result[result.size()-1].H.Ti << " " << result[result.size()-1].H.Tf << " " << 
          // result[result.size()-1].C.name << " " << result[result.size()-1].C.Ti << " " << result[result.size()-1].C.Tf << " " << endl;
  if(abs(HS[val].Ti-HS[val].Tf)<0.01) HS.erase(HS.begin());
  if(abs(CS[val].Ti-CS[val].Tf)<0.01) CS.erase(CS.begin());

  while(HSnp.size()!=0)
  {
   if(HSnp[0].Ti>=stopT)
     {
      HS.push_back(HSnp[0]);
      HSnp.erase(HSnp.begin());
     }
   else break;
  }
  
  sort(HS.begin(), HS.end(), comp);
  sort(CS.begin(), CS.end(), TempInitialGreater);
 }
 return;
}

void matchBP(vector<stream> &HS, vector<stream> &CS, vector<stream> &HSnp, vector<stream> &CSnp, 
              double dTmin, double Tci, double Tcf, double Thi, double Thf, vector<HeatExchanger> &result)
{
 if(HS.size()+HSnp.size() == 0|| CS.size()+CSnp.size() == 0) return;

 //SPLIT HOTSTREAM
 if(HS.size()<CS.size())
 {
  HS.insert(HS.begin(), HS[0]);
  HS[0].name += "-s1";
  HS[1].name += "-s2";
  HS[0].CP *= CS[0].CP/(CS[0].CP+CS[1].CP);
  HS[1].CP *= CS[1].CP/(CS[0].CP+CS[1].CP);
  sort(HS.begin(), HS.end(), comp);
  matchBP(HS, CS, HSnp, CSnp, dTmin, Tci, Tcf, Thi, Thf, result);
  return;
 }
 
 //SPLIT COLDSTREAM
 if(fitMatch(CS, HS) == false)
 {
  CS.insert(CS.begin(), CS[0]);
  CS[0].name += "-s1";
  CS[1].name += "-s2";
  CS[0].CP /= 2.0;
  CS[1].CP /= 2.0;
  sort(CS.begin(), CS.end(), comp);
  matchBP(HS, CS, HSnp, CSnp, dTmin, Tci, Tcf, Thi, Thf, result);
  return;
 }

 double stopT, sumHcp=0.0, sumCcp=0.0;
 for(int i=0; i<CS.size(); i++)  sumCcp+=CS[i].CP;
 for(int i=0; i<HS.size(); i++)  sumHcp+=HS[i].CP;

  if(CS.size() == 0){vector<stream> k(0); matchAllBP(HS, CSnp , HSnp, k, dTmin, result);}
  stopT = -100000.0;
  for(int i=0; i<CSnp.size(); i++) stopT = max(stopT, CSnp[i].Tf+dTmin);
  for(int i=0; i<CS.size(); i++) stopT = max(stopT, HS[i].Tf);

  int i=0, j=0;
  while(i<CS.size())
  {
   result.push_back(matchHotstream(HS, i, CS, j, dTmin, stopT));
   //cout << result[result.size()-1].H.name << " " << result[result.size()-1].H.Ti << " " << result[result.size()-1].H.Tf << " " << 
           //result[result.size()-1].C.name << " " << result[result.size()-1].C.Ti << " " << result[result.size()-1].C.Tf << " " << endl;

   if(abs(HS[i].Ti-HS[i].Tf)<0.01)
   {  sumHcp -= HS[i].CP; HS.erase(HS.begin()+i); i--;}
   if(abs(CS[j].Ti-CS[j].Tf)<0.01)
   {  sumCcp -= CS[j].CP; CS.erase(CS.begin()+j); j--;}
   i++; j++;
  }

  while(HSnp.size()>0)
  {
   if(HSnp[0].Tf >=stopT)
    {
     HS.push_back(HSnp[0]);
     sumHcp += HSnp[0].CP;
     HSnp.erase(HSnp.begin());
    }
   else break;
  }
  
  while(CSnp.size()>0)
  {
   if(CSnp[0].Tf >= stopT-dTmin)
    {
     CS.push_back(CSnp[0]);
     sumCcp += CSnp[0].CP;
     CSnp.erase(CSnp.begin());
    }
   else break;
  }
  
  if(sumCcp>sumHcp)
  {
   for(i=0; i<HSnp.size(); i++) HS.push_back(HSnp[i]);
   for(i=0; i<CSnp.size(); i++) CS.push_back(CSnp[i]);
   HSnp.clear(); CSnp.clear();
   vector<stream> HStemp, CStemp;
   giveStreamsAbovePinch(HS, CS, HStemp, CStemp, HSnp, CSnp, Tci, Tcf, Thi, Thf);
   matchAllAP(HStemp, CStemp, HSnp, CSnp, dTmin, result);
   return;
  }
  
  else
 {
  matchBP(HS, CS, HSnp, CSnp, dTmin, Tci, Tcf, Thi, Thf, result);
  return;
 }

 return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matchAP(vector<stream> &HS, vector<stream> &CS, vector<stream> &HSnp, vector<stream> &CSnp, 
              double dTmin, double Tci, double Tcf, double Thi, double Thf, vector<HeatExchanger> &result)
{
 if(HS.size()+HSnp.size() == 0|| CS.size()+CSnp.size() == 0) return;

 //SPLIT COLDSTREAM
 if(HS.size()>CS.size())
 {
  CS.insert(CS.begin(), CS[0]);
  CS[0].name += "-s1";
  CS[1].name += "-s2";
  CS[0].CP *= HS[0].CP/(HS[0].CP+HS[1].CP);
  CS[1].CP *= HS[1].CP/(HS[0].CP+HS[1].CP);
  sort(CS.begin(), CS.end(), comp);
  matchAP(HS, CS, HSnp, CSnp, dTmin, Tci, Tcf, Thi, Thf, result);
 }
 
 //SPLIT HOTSTREAM
 if(fitMatch(HS, CS) == false)
 {
  HS.insert(HS.begin(), HS[0]);
  HS[0].name += "-s1";
  HS[1].name += "-s2";
  HS[0].CP /= 2.0;
  HS[1].CP /= 2.0;
  sort(HS.begin(), HS.end(), comp);
  matchAP(HS, CS, HSnp, CSnp, dTmin, Tci, Tcf, Thi, Thf, result);
 }

 double stopT, sumHcp=0.0, sumCcp=0.0;
 for(int i=0; i<CS.size(); i++)  sumCcp+=CS[i].CP;
 for(int i=0; i<HS.size(); i++)  sumHcp+=HS[i].CP;

  if(HS.size() == 0){vector<stream> k(0); matchAP(HSnp, CS, k, CSnp, dTmin, Tci, Tcf, Thi, Thf, result);}
  stopT = 100000.0;
  for(int i=0; i<HSnp.size(); i++) stopT = min(stopT, HSnp[i].Tf-dTmin);
  for(int i=0; i<HS.size(); i++) stopT = min(stopT, CS[i].Tf);

  int i=0, j=0;
  while(i<HS.size())
  {
   result.push_back(matchColdstream(HS, i, CS, j, dTmin, stopT));
   //cout << result[result.size()-1].H.name << " " << result[result.size()-1].H.Ti << " " << result[result.size()-1].H.Tf << " " << 
           //result[result.size()-1].C.name << " " << result[result.size()-1].C.Ti << " " << result[result.size()-1].C.Tf << " " << endl;

   if(abs(HS[i].Ti-HS[i].Tf)<0.00001)
   {  sumHcp -= HS[i].CP; HS.erase(HS.begin()+i); i--;}
   if(abs(CS[j].Ti-CS[j].Tf)<0.00001)
   {  sumCcp -= CS[j].CP; CS.erase(CS.begin()+j); j--;}
   i++; j++;
  }

  while(HSnp.size()>0)
  {
   if(HSnp[0].Tf <=stopT+dTmin)
    {
     HS.push_back(HSnp[0]);
     sumHcp += HSnp[0].CP;
     HSnp.erase(HSnp.begin());
    }
   else break;
  }
  
  while(CSnp.size()>0)
  {
   if(CSnp[0].Ti <= stopT)
    {
     CS.push_back(CSnp[0]);
     sumCcp += CSnp[0].CP;
     CSnp.erase(CSnp.begin());
    }
   else break;
  }
  
  if(sumCcp<sumHcp)
  {
   for(i=0; i<HSnp.size(); i++) HS.push_back(HSnp[i]);
   for(i=0; i<CSnp.size(); i++) CS.push_back(CSnp[i]);
   HSnp.clear(); CSnp.clear();
   vector<stream> HStemp, CStemp;
   giveStreamsBelowPinch(HS, CS, HStemp, CStemp, HSnp, CSnp, Tci, Tcf, Thi, Thf);
   matchAllBP(HStemp, CStemp, HSnp, CSnp, dTmin, result);
   return;
  }
  else
  {
   matchAP(HS, CS, HSnp, CSnp, dTmin, Tci, Tcf, Thi, Thf, result);
   return;
  }
}


int main()
{
 int n, index=-1; double dTmin, minH = 0.0, p;
 vector<stream> HS, CS, HStemp, HSnp, CStemp, CSnp;
 vector<double> T, dH, nearPinchesC, nearPinchesH, CPh, CPc;
 cout << "Enter number of streams you want: ";
 cin >> n;
 cout << "Enter value of Delta T: ";
 cin >> dTmin;

 for(int i=0; i<n; i++)
 {
   stream temp;
   cout << "Enter stream name, Ti, Tf and CP (in order) of stream " << i+1 << ": ";
   cin>>temp.name; cin>>temp.Ti; cin>>temp.Tf; cin>>temp.CP; temp.id = i+1;
   if(temp.Ti<temp.Tf){ CS.push_back(temp); T.push_back(temp.Ti+(dTmin/2.0)); T.push_back(temp.Tf+(dTmin/2.0));}
   else{ HS.push_back(temp); T.push_back(temp.Ti-(dTmin/2.0)); T.push_back(temp.Tf-(dTmin/2.0));}
 }
 
 //printStreams(HS, "Hotstreams");
 //printStreams(CS,  "Coldstreams");
 
 sort(T.begin(), T.end(), greater<double>());
 removeDuplicates(T);

 dH.resize(T.size(), 0.0);
 CPh.resize(T.size(), 0.0);
 CPc.resize(T.size(), 0.0);
 for(int i=0; i<HS.size(); i++)
   for(int j=0; j<T.size()-1; j++)
    {
      if(HS[i].Ti-(dTmin/2)>=T[j] && HS[i].Tf-(dTmin/2)<=T[j+1])
       {
        dH[j+1] = dH[j+1] + HS[i].CP*(T[j]-T[j+1]);
        CPh[j] += HS[i].CP;
       }
      
      else continue;
    }
 for(int i=0; i<CS.size(); i++)
   for(int j=0; j<T.size()-1; j++)
    {
      if(CS[i].Tf+(dTmin/2)>=T[j] && CS[i].Ti+(dTmin/2)<=T[j+1])
       {
        dH[j+1] = dH[j+1] - CS[i].CP*(T[j]-T[j+1]);
        CPc[j] += CS[i].CP;
       }
      
      else continue;
    }

 for(int i=1; i<dH.size(); i++) 
 {
  dH[i] += dH[i-1];
  if(dH[i] < minH)
   { minH=dH[i]; index=i;}
 }
 

 if(index==-1) cout << "No Pinch Temperature." << endl;
 else 
 {
   findTempDiff(CPh, CPc, T, index, dTmin, nearPinchesC, nearPinchesH);
   printVector(nearPinchesC, "COLD near pinches found");
   printVector(nearPinchesH, "HOT near pinches found");
   cout << "Pinch Temperature = " << T[index] << endl << endl;
   cout << "Hot Utility = " << -minH << endl << endl;
   cout << "Cold Utility = " << dH[dH.size()-1]-minH << endl << endl;
 }

 for(int i=0; i<nearPinchesC.size(); i++) 
   if(nearPinchesC[i]+(dTmin/2.0) == T[index])
     {
       p=i; break;
     }

 sort(HS.begin(), HS.end(), comp);
 sort(CS.begin(), CS.end(), comp);

 stringstream ss1, ss2;
 
 vector<HeatExchanger> result;
 for(int i=p; i<nearPinchesC.size(); i++)
 {
  if(i == nearPinchesC.size()-1)
  {
   giveStreamsAbovePinch(HS, CS, HStemp, CStemp, HSnp, CSnp, nearPinchesC[i], 100000.0, 100000.0, nearPinchesH[i]);
   matchAllAP(HStemp, CStemp, HSnp, CSnp, dTmin, result);
   ss1.str("");
   ss1<<i+1;
   string t = "Heat Exchanger Network above " + ss1.str() + "th pinch";
   printHeatExchangers(result, t);
  }
  else
  {
   giveStreamsAbovePinch(HS, CS, HStemp, CStemp, HSnp, CSnp, nearPinchesC[i], nearPinchesC[i+1], nearPinchesH[i+1], nearPinchesH[i]);
   matchAP(HStemp, CStemp, HSnp, CSnp, dTmin, nearPinchesC[i], nearPinchesC[i+1], nearPinchesH[i+1], nearPinchesH[i], result);
   ss1.str(""); ss2.str("");
   ss1<<i+1;    ss2<<i+2;
   string t = "Heat Exchanger Network between " + ss1.str() + "th and" + ss2.str() + "th pinch"; 
   printHeatExchangers(result, t);
  }
  HStemp.clear(); CStemp.clear();
  HSnp.clear(); CSnp.clear();
  result.clear();
 }
 
 for(int i=p; i>=0; i--)
 {
  if(i == 0)
  {
   giveStreamsBelowPinch(HS, CS, HStemp, CStemp, HSnp, CSnp, -100000.0, nearPinchesC[i], nearPinchesH[i], -100000.0);
   matchAllBP(HStemp, CStemp, HSnp, CSnp, dTmin, result);
   ss1.str("");
   ss1<<i+1;
   string t = "Heat Exchanger Network below " + ss1.str() + "th pinch";
   printHeatExchangers(result, t);
  }
  else
  {
   giveStreamsBelowPinch(HS, CS, HStemp, CStemp, HSnp, CSnp, nearPinchesC[i-1], nearPinchesC[i], nearPinchesH[i], nearPinchesH[i-1]);
   matchBP(HStemp, CStemp, HSnp, CSnp, dTmin, nearPinchesC[i-1], nearPinchesC[i], nearPinchesH[i], nearPinchesH[i-1], result);
   ss1.str(""); ss2.str("");
   ss1<<i;    ss2<<i+1;
   string t = "Heat Exchanger Network between " + ss1.str() + "th and" + ss2.str() + "th pinch"; 
   printHeatExchangers(result, t);
  }
  HStemp.clear(); CStemp.clear();
  HSnp.clear(); CSnp.clear();
  result.clear();
 }

 return 0;
}

