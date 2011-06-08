#include <string>

ledfill(TLED **led, TBranch *branch, string date, Int_t parts,
      Int_t year, Int_t month, Int_t day, Int_t part) {
  string folder, noise[parts], ledfile;
  ostringstream oss;
  for(Int_t i = 0; i < 11; i++) {
    folder[i] = date[i];
  }
  oss.str("");
  oss << folder << "/" << date << ".led355.root";
  ledfile = oss.str();
  for(Int_t i = 0; i < parts; i++) {
    oss.str("");
    oss << folder << "/" << folder << "p0" << i+1 << ".noise-closed.root";
    noise[i] = oss.str();
  }
  cout << "ledfile: " << ledfile << endl;
  for(Int_t i = 0; i < parts; i++) {
    cout << "noise["<< i << "]: " << noise[i] << endl;
  }
  cout << "*led = new TLED(" << year << ',' << month << ',' << day << ',' << part << ',' << *ledfile << ',' << *noise << ',' << parts << ");" << endl;
//  *led = new TLED(year,month,day,part,ledfile,noise,parts);
//  if(led == 0) {
//    cerr << "Error: TLED for " << date << " not created\n";
//  }
//  branch->Fill();
//  cout << "branch filled for " << date << endl << endl;
//  delete *led;
//  *led = 0;
}
