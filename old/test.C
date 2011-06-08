{
gSystem->Load("libMd_Calib");
TFile file("ledtree.root","recreate");
TTree tree("T","led data through Sept. 4th, 2008");
TLED *led = 0;
TBranch *branch = tree.Branch("led","TLED",&led,65536,0);
string date = "y2008m07d29p02";
Int_t year = 2008;
Int_t month = 7;
Int_t day = 29;
Int_t part = 2;
Int_t parts = 2;
folder = string("tmp/middle_drum/20080729/");
string ledfile, noise[parts];
oss = ostringstream("");
oss << folder << date << ".00.led355.root";
ledfile = oss.str();
for(Int_t i = 0; i < parts; i++)
{
  oss.str("");
  oss << folder << "p0" << i+1 << ".noise-closed.root";
  noise[i] = oss.str();
}
cout << "*led = new TLED(" << year << ',' << month << ',' << day << ',' << part << ",\"" << ledfile << "\",\"" << *noise << "\"," << parts << ");" << endl;
*led = new TLED(year,month,day,part,ledfile,noise,parts);
if(led == 0) {
  cerr << "Error: TLED for " << date << " not created\n";
}
//branch->Fill();
//cout << "branch filled for " << date << endl << endl;
//delete *led;
//*led = 0;
}
