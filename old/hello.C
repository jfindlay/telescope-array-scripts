{
//TBrowser b;

//gSystem->Load("libAstro");
//gSystem->Load("libMisc");
//gSystem->Load("libEvent");
//gSystem->Load("libDst");
//gSystem->Load("libMd_Calib.so");
//TFile f("data/scott/ledtree.root");
//TTree*t = f.Get("T");
//TBranch*b = t->GetBranch("led");
//TLED*led = 0;
//b->SetAddress(&led);
//b->GetEntry(0);
//led->GetFlashes(1);

//TLED*led2 = 0;
//TBranch*B2 = t->GetBranch("rxf");
//B2->SetAddress(&led2);
//t->Print();
//led->GetDate();
//led->GetGain(2,-1,1);
//TH1F*h = led->GetHist(1,1);
//h->Draw();

#include <stdlib.h>
#include <vector>
typedef vector<TArrayD> vector_TArrayD;

gSystem->Load("libAstro");
gSystem->Load("libMisc");
gSystem->Load("libEvent");
gSystem->Load("libDst");
gSystem->Load("libMd_Calib.so");
file = TFile("tmp/middle_drum/20081008/y2008m10d08p01.led355.root");
TTree*tree = file.Get("T");
TBranch*branch = tree->GetBranch("event");
THPKT1_DST_EVENT*event = 0;
branch->SetAddress(&event);
cout << branch->GetEntries() << endl;

enum eConstants {MIN_CRATE = 0, MAX_CRATE = 14, NCRATE = 15, MIN_TUBE = 0, MAX_TUBE = 255, NTUBE = 256};
vector<TArrayD>  fQdcb[NCRATE];
TArrayD flash(NTUBE) = {0};
Int_t br_entries = branch->GetEntries();
for(Int_t nentry = 0; nentry < br_entries; nentry++)
{
  branch->GetEntry(nentry);
  if(event->ntubes > 240)
  {
    if(event->pktHdr_type == THPKT1_DST_EVENT::eEVENT)
    {
      for(Int_t ntube = 0; ntube < event->ntubes; ntube++)
      {
        // QDC truncates value to integer.  This spreads out the data
        // continuously throughout the interval like it originally was
        flash[event->tube_num[ntube]] = (double)event->qdcB[ntube] + (double)rand()/(double)RAND_MAX;
      }
      fQdcb[event->pktHdr_crate].push_back(flash);
    }
  }
  else cout << "bad flash: " << event->pktHdr_crate << endl;
}
}

void FillQdcb(TBranch *br_event)
{
  srand(time(NULL));
  THPKT1_DST_EVENT *event = 0;
  br_event->SetAddress(&event);
  TArrayD flash(NTUBE);
  Int_t br_entries = br_event->GetEntries();

  for(Int_t nentry = 0; nentry < br_entries; nentry++) {
    br_event->GetEntry(nentry);
    if(event->ntubes > 240) {
      if(event->pktHdr_type == THPKT1_DST_EVENT::eEVENT) {
        for(Int_t ntube = 0; ntube < NTUBE; ntube++)
          flash[ntube] = -1.0;
        for(Int_t ntube = 0; ntube < event->ntubes; ntube++) {
          if(event->qdcB[ntube] > 3900)
            fErrFlag[event->pktHdr_crate][event->tube_num[ntube]] = -1;
          flash[event->tube_num[ntube]] = (double)event->qdcB[ntube]+
            (double)rand()/(double)RAND_MAX;
        }
        fQdcb[event->pktHdr_crate].push_back(flash);
      }
    }
    else fBadFlash[event->pktHdr_crate]++;
  }
}
