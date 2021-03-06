{
  gROOT->ProcessLine(".L macros/treefill.C");
  TFile file("ppetree.root","recreate");
  TTree tree("T","rxf data through June 9th, 2008");
  TRXF *rxf = 0;
  TBranch *br_rxf = tree.Branch("rxf","TRXF",&rxf,64000,0);
  Char_t date[15] = "y2007m12d29";
  Int_t rxfpart[15];
  for(int i = 0; i <15; i++) rxfpart[i] = 1;
  rxfpart[1] = rxfpart[2] = 2;
  fillrxf(&rxf,br_rxf,2007,12,29,date,rxfpart,2);
  strcpy(date,"y2008m01d28");
  rxfpart[1] = 3;
  rxfpart[2] = 1;
  fillrxf(&rxf,br_rxf,2008,01,28,date,rxfpart,1);
  strcpy(date,"y2008m02d01");
  rxfpart[1] = rxfpart[4] = rxfpart[6] = 3;
  fillrxf(&rxf,br_rxf,2008,02,13,date,rxfpart,2);
  strcpy(date,"y2008m02d13");
  rxfpart[1] = rxfpart[4] = rxfpart[6] = 1;
  fillrxf(&rxf,br_rxf,2008,02,13,date,rxfpart,2);
  strcpy(date,"y2008m03d02");
  fillrxf(&rxf,br_rxf,2008,03,02,date,rxfpart,1);
  strcpy(date,"y2008m04d26");
  rxfpart[1] = 2;
  fillrxf(&rxf,br_rxf,2008,04,26,date,rxfpart,2);
  strcpy(date,"y2008m04d30");
  fillrxf(&rxf,br_rxf,2008,04,30,date,rxfpart,1);
  strcpy(date,"y2008m06d05");
  fillrxf(&rxf,br_rxf,2008,06,05,date,rxfpart,2);
  strcpy(date,"y2008m06d09");
  fillrxf(&rxf,br_rxf,2008,06,09,date,rxfpart,2);
  br_rxf = 0;
  file.cd();
  tree->Write();
  file.Close();
}
