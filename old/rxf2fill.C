#include <string>

rxf2fill(TLED **led, TBranch *branch, string date, Int_t parts,
	Int_t year, Int_t month, Int_t day) {
  Int_t part = -1;
  string noise[parts], rxffile;
  ostringstream oss;
  oss.str("");
  oss << date << "/" << date << ".rxf.all.root";
  rxffile = oss.str();
  for(Int_t i = 0; i < parts; i++) {
    oss.str("");
    oss << date << "/" << date << "p0" << i+1 << ".noise-closed.root";
    noise[i] = oss.str();
  }
  cout << "rxffile: " << rxffile << endl;
  for(Int_t i = 0; i < parts; i++) {
    cout << "noise["<< i << "]: " << noise[i] << endl;
  }
  *led = new TLED(year,month,day,part,rxffile,noise,parts);
  if(led == 0) {
    cerr << "Error: TLED for " << date << " not created\n";
  }
  branch->Fill();
  cout << "branch filled for " << date << endl << endl;
  delete *led;
  *led = 0;
}
