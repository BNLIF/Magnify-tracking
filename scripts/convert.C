void convert(){
  ifstream infile("organized.dat");
  TFile *file = new TFile("temp.root","RECREATE");
  TTree *T = new TTree("T","T");
  Int_t run_no, subrun_no, event_no;
  Double_t ks1, ks2;
  Double_t norm1, norm2;
  Int_t flag;
  T->Branch("run_no",&run_no,"run_no/I");
  T->Branch("subrun_no",&subrun_no,"subrun_no/I");
  T->Branch("event_no",&event_no,"event_no/I");
  T->Branch("ks1",&ks1,"ks1/D");
  T->Branch("ks2",&ks2,"ks2/D");
  T->Branch("norm1",&norm1,"norm1/D");
  T->Branch("norm2",&norm2,"norm2/D");
  T->Branch("flag",&flag,"flag/I");
  
  while (infile >> run_no >> subrun_no >> event_no >> ks1 >> ks2 >> norm1 >> norm2 >> flag) {
    T->Fill();
  }
  file->Write();
  file->Close();
}
