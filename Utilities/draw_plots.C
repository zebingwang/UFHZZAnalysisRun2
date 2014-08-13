{

  char* filename = 
    "GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT_testAna_eidTight_FSRrecovery.root"
  ;

  char plotfile1[300];
  char plotfile2[300];
  sprintf(plotfile1, "%s.ps", filename);
  sprintf(plotfile2, "%s.pdf", filename);

  char name[1000];

  TFile* file = TFile::Open( filename );

  TTree* tree = (TTree*)file->Get("AnaAfterHlt/passedEvents");
 
  TH1D* h1;
  TH1D* h2;
  TH1D* h3;
  TH1D* h4;
  TH1D* h5;
  TH1D* h6;
   
  TH2D* h2d1;
  TH2D* h2d2;


  TLegend* lg;

  //TH1D* hmass4l = new TH1D("hmass4l", "hmass4l", 100, 100, 150);
  //TH1D* hmass4lNoFSR = new TH1D("hmass4lNoFSR", "hmass4lNoFSR", 100, 100, 150);
  //TH1D* hptl1 = new TH1D("hptl1", "hptl1", 100, 0, 100); 
  //TH1D* hptl2 = new TH1D("hptl2", "hptl2", 100, 0, 100); 
  //TH1D* hptl3 = new TH1D("hptl3", "hptl3", 100, 0, 100); 
  //TH1D* hptl4 = new TH1D("hptl4", "hptl4", 100, 0, 100); 
  //TH2D* hmassz1z2 = new TH2D("hmassz1z2", "hmassz1z2", 100, 100, 150);

  TCanvas* plots = new TCanvas("plots", "plots");
  sprintf(name, "%s[", plotfile1);
  plots->Print(name);

  // mass
  h1 = new TH1D("h1", "hmass4l", 100, 100, 150);
  h2 = new TH1D("h2", "hmass4lNoFSR", 100, 100, 150);
  tree->Draw("mass4l>>h1");
  tree->Draw("mass4lNoFSR>>h2");
  h1->SetName("hmass4l");
  h2->SetName("hmass4lNoFSR");
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  
  lg = new TLegend(0.56,0.7, 0.9, 0.9);
  lg->SetName("lgmass4l");
  lg->AddEntry(h1, "Mass 4l with FSR recovery", "pl");
  lg->AddEntry(h2, "Mass 4l without FSR recovery", "pl");

  plots->Clear();
  h1->Draw();
  h2->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();

  
  // fsr pt
  h1 = new TH1D("h1", "hFSRPhot1_Pt", 100, 0, 50);
  h2 = new TH1D("h2", "hFSRPhot2_Pt", 100, 0, 50);
  tree->Draw("FSRPhot1_Pt>>h1", "FSRPhot1_Pt>0");
  tree->Draw("FSRPhot2_Pt>>h2", "FSRPhot2_Pt>0");
  h1->SetName("hFSRPhot1_Pt");
  h2->SetName("hFSRPhot2_Pt");
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h1->GetXaxis()->SetTitle("FSR Pt (GeV)");

  lg = new TLegend(0.56,0.7, 0.9, 0.9);
  lg->SetName("lgFSRPhot_Pt");
  lg->AddEntry(h1, "FSR added to Z1", "pl");
  lg->AddEntry(h2, "FSR added to Z2", "pl");

  plots->Clear();
  h1->Draw();
  h2->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();

  // pt l 2e2mu
  h1 = new TH1D("h1", "hptl1_2e2mu", 100, 0, 100);
  h2 = new TH1D("h2", "hptl2_2e2mu", 100, 0, 100);
  h3 = new TH1D("h3", "hptl3_2e2mu", 100, 0, 100);
  h4 = new TH1D("h4", "hptl4_2e2mu", 100, 0, 100);
  tree->Draw("pTL1>>h1", "mass2e2mu>0");
  tree->Draw("pTL2>>h2", "mass2e2mu>0");
  tree->Draw("pTL3>>h3", "mass2e2mu>0");
  tree->Draw("pTL4>>h4", "mass2e2mu>0");
  h1->SetName("hptl1_2e2mu");
  h2->SetName("hptl2_2e2mu");
  h3->SetName("hptl3_2e2mu");
  h4->SetName("hptl4_2e2mu");
  h1->GetXaxis()->SetTitle("pT(l) (GeV)");
  h2->GetXaxis()->SetTitle("pT(l) (GeV)");
  h3->GetXaxis()->SetTitle("pT(l) (GeV)");
  h4->GetXaxis()->SetTitle("pT(l) (GeV)");
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(6);
  h4->SetLineColor(8);
  
  lg = new TLegend(0.56,0.7, 0.9, 0.9);
  lg->SetName("lgptl_2e2mu");
  lg->AddEntry(h1, "Lepton 1", "pl");
  lg->AddEntry(h2, "Lepton 2", "pl");
  lg->AddEntry(h3, "Lepton 3", "pl");
  lg->AddEntry(h4, "Lepton 4", "pl");

  plots->Clear();
  h4->Draw();
  h2->Draw("same");
  h3->Draw("same");
  h1->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();


  // pt l 4e
  h1 = new TH1D("h1", "hptl1_4e", 100, 0, 100);
  h2 = new TH1D("h2", "hptl2_4e", 100, 0, 100);
  h3 = new TH1D("h3", "hptl3_4e", 100, 0, 100);
  h4 = new TH1D("h4", "hptl4_4e", 100, 0, 100);
  tree->Draw("pTL1>>h1", "mass4e>0");
  tree->Draw("pTL2>>h2", "mass4e>0");
  tree->Draw("pTL3>>h3", "mass4e>0");
  tree->Draw("pTL4>>h4", "mass4e>0");
  h1->SetName("hptl1_4e");
  h2->SetName("hptl2_4e");
  h3->SetName("hptl3_4e");
  h4->SetName("hptl4_4e");
  h1->GetXaxis()->SetTitle("pT(l) (GeV)");
  h2->GetXaxis()->SetTitle("pT(l) (GeV)");
  h3->GetXaxis()->SetTitle("pT(l) (GeV)");
  h4->GetXaxis()->SetTitle("pT(l) (GeV)");
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(6);
  h4->SetLineColor(8);

  lg = new TLegend(0.56,0.7, 0.9, 0.9);
  lg->SetName("lgptl_4e");
  lg->AddEntry(h1, "Lepton 1", "pl");
  lg->AddEntry(h2, "Lepton 2", "pl");
  lg->AddEntry(h3, "Lepton 3", "pl");
  lg->AddEntry(h4, "Lepton 4", "pl");

  plots->Clear();
  h4->Draw();
  h2->Draw("same");
  h3->Draw("same");
  h1->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();


  // pt l 4mu
  h1 = new TH1D("h1", "hptl1_4mu", 100, 0, 100);
  h2 = new TH1D("h2", "hptl2_4mu", 100, 0, 100);
  h3 = new TH1D("h3", "hptl3_4mu", 100, 0, 100);
  h4 = new TH1D("h4", "hptl4_4mu", 100, 0, 100);
  tree->Draw("pTL1>>h1", "mass4mu>0");
  tree->Draw("pTL2>>h2", "mass4mu>0");
  tree->Draw("pTL3>>h3", "mass4mu>0");
  tree->Draw("pTL4>>h4", "mass4mu>0");
  h1->SetName("hptl1_4mu");
  h2->SetName("hptl2_4mu");
  h3->SetName("hptl3_4mu");
  h4->SetName("hptl4_4mu");
  h1->GetXaxis()->SetTitle("pT(l) (GeV)");
  h2->GetXaxis()->SetTitle("pT(l) (GeV)");
  h3->GetXaxis()->SetTitle("pT(l) (GeV)");
  h4->GetXaxis()->SetTitle("pT(l) (GeV)");
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(6);
  h4->SetLineColor(8);

  lg = new TLegend(0.56,0.7, 0.9, 0.9);
  lg->SetName("lgptl_4mu");
  lg->AddEntry(h1, "Lepton 1", "pl");
  lg->AddEntry(h2, "Lepton 2", "pl");
  lg->AddEntry(h3, "Lepton 3", "pl");
  lg->AddEntry(h4, "Lepton 4", "pl");

  plots->Clear();
  h4->Draw();
  h3->Draw("same");
  h2->Draw("same");
  h1->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();


  // mass Z1 vs Z2
  h2d1 = new TH2D("h2d1", "hmassZ1Z2", 50, 10, 110, 50, 0, 100);
  tree->Draw("massZ2:massZ1>>h2d1", "", "colz");
  h2d1->SetName("hmassZ1Z2");

  h2d1->GetXaxis()->SetTitle("mass Z1 (GeV)");  
  h2d1->GetYaxis()->SetTitle("mass Z2 (GeV)");  

  plots->Clear();
  h2d1->Draw("colz");    
  plots->Print(plotfile1);
  plots->Clear();


  // end
  sprintf(name, "%s]", plotfile1);
  plots->Print(name);
  sprintf(name, ".! ps2pdf %s %s", plotfile1, plotfile2);
  gROOT->ProcessLine(name);

}
