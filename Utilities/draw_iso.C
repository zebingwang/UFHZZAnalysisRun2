{

  TFile* file1 = TFile::Open("GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_8TeV_1.root");
  //TFile* file2 = TFile::Open("GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT_testAna_eidTight_FSRrecovery_dumpTrees.root");
  TFile* file2 = TFile::Open("GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT_testAna.root");
  TFile* file3 = TFile::Open("GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU_S14_PAT_big_testAna.root");

  TTree* tree1 = (TTree*)file1->Get("AnaAfterHlt/muons/muonDumpTree");
  TTree* tree2 = (TTree*)file2->Get("AnaAfterHlt/muons/muonDumpTree");
  TTree* tree3 = (TTree*)file3->Get("AnaAfterHlt/muons/muonDumpTree");

  TTree* tree4 = (TTree*)file1->Get("AnaAfterHlt/electrons/electronDumpTree");
  TTree* tree5 = (TTree*)file2->Get("AnaAfterHlt/electrons/electronDumpTree");
  TTree* tree6 = (TTree*)file3->Get("AnaAfterHlt/electrons/electronDumpTree");


  char* filename = "draw_iso";

  char plotfile1[300];
  char plotfile2[300];
  sprintf(plotfile1, "%s.ps", filename);
  sprintf(plotfile2, "%s.pdf", filename);

  char name[1000];
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TH1D* h1;
  TH1D* h2;
  TH1D* h3;
  TH1D* h4;
  TH1D* h5;
  TH1D* h6;

  TH2D* h2d1;
  TH2D* h2d2;

  TLegend* lg;

  TCanvas* plots = new TCanvas("plots", "plots");
  sprintf(name, "%s[", plotfile1);
  plots->Print(name);

  h1 = new TH1D("h1", "h1", 10, 0, 4);
  h2 = new TH1D("h2", "h2", 10, 0, 4);
  h3 = new TH1D("h3", "h3", 10, 0, 4);
  
  tree1->Draw("(isoCH+max(isoPhot+isoNH-0.5*dB,0))/pT>>h1");
  h1->Scale(1.0/h1->Integral());
  tree2->Draw("(isoCH+max(isoPhot+isoNH-0.5*dB,0))/pT>>h2");
  h2->Scale(1.0/h2->Integral());
  tree3->Draw("(isoCH+max(isoPhot+isoNH-0.5*dB,0))/pT>>h3");
  h3->Scale(1.0/h3->Integral());

  h1->SetName("hIsoMu8TeV");
  h2->SetName("hIsoMu13TeV25ns");
  h3->SetName("hIsoMu13TeV50ns");
  
  h1->SetTitle("hIsoMu");
  h2->SetTitle("hIsoMu");
  h3->SetTitle("hIsoMu");

  h1->GetXaxis()->SetTitle("ISO(muon)");
  h2->GetXaxis()->SetTitle("ISO(muon)");
  h3->GetXaxis()->SetTitle("ISO(muon)");

  h1->GetYaxis()->SetTitle("norm.");
  h2->GetYaxis()->SetTitle("norm.");
  h3->GetYaxis()->SetTitle("norm.");

  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(8);

  lg = new TLegend(0.5,0.6,0.9,0.9);
  lg->SetName("lg_IsoMuon");
  lg->AddEntry(h1, "8 TeV  50ns bx", "l");
  lg->AddEntry(h2, "13 TeV PU20bx25", "l");
  lg->AddEntry(h3, "13 TeV PU 14S", "l");

  plots->Clear();
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();


  plots->Clear();
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  lg->Draw("same");
  plots->SetLogy(1);
  plots->Print(plotfile1);
  plots->SetLogy(0);
  plots->Clear();

  // el ISO
  h1 = new TH1D("h1", "h1", 10, 0, 4);
  h2 = new TH1D("h2", "h2", 10, 0, 4);
  h3 = new TH1D("h3", "h3", 10, 0, 4);

  tree1->Draw("relIso>>h1");
  h1->Scale(1.0/h1->Integral());
  tree2->Draw("relIso>>h2");
  h2->Scale(1.0/h2->Integral());
  tree3->Draw("relIso>>h3");
  h3->Scale(1.0/h3->Integral());

  h1->SetName("hIsoEl8TeV");
  h2->SetName("hIsoEl13TeV25ns");
  h3->SetName("hIsoEl13TeV50ns");

  h1->SetTitle("hIsoEl");
  h2->SetTitle("hIsoEl");
  h3->SetTitle("hIsoEl");

  h1->GetXaxis()->SetTitle("ISO(e)");
  h2->GetXaxis()->SetTitle("ISO(e)");
  h3->GetXaxis()->SetTitle("ISO(e)");

  h1->GetYaxis()->SetTitle("norm.");
  h2->GetYaxis()->SetTitle("norm.");
  h3->GetYaxis()->SetTitle("norm.");

  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(8);


  lg = new TLegend(0.5,0.6,0.9,0.9);
  lg->SetName("lg_IsoEl");
  lg->AddEntry(h1, "8 TeV  50ns bx", "l");
  lg->AddEntry(h2, "13 TeV PU20bx25", "l");
  lg->AddEntry(h3, "13 TeV PU 14S", "l");

  plots->Clear();
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();


  plots->Clear();
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  lg->Draw("same");
  plots->SetLogy(1);
  plots->Print(plotfile1);
  plots->SetLogy(0);
  plots->Clear();


  // nvtx
  h1 = (TH1D*)file1->Get("AnaAfterHlt/nVtx_ReWeighted");
  h2 = (TH1D*)file2->Get("AnaAfterHlt/nVtx_ReWeighted");
  h3 = (TH1D*)file3->Get("AnaAfterHlt/nVtx_ReWeighted");
  h1->SetName("hnvtxrw8tev");
  h2->SetName("hnvtxrw13tev25ns");
  h3->SetName("hnvtxrw13tev50ns");

  h1->Scale(1.0/h1->Integral());
  h2->Scale(1.0/h2->Integral());
  h3->Scale(1.0/h3->Integral());

  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(8);

  lg = new TLegend(0.6,0.7,0.9,0.9);
  lg->SetName("lg_nvtxrw");
  lg->AddEntry(h1, "8 TeV  50ns bx", "l");
  lg->AddEntry(h2, "13 TeV PU20bx25", "l");
  lg->AddEntry(h3, "13 TeV PU 14S", "l");

  plots->Clear();
  h2->Draw();
  h1->Draw("same");
  h3->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();  


  // nvtx
  h1 = (TH1D*)file1->Get("AnaAfterHlt/nVtx");
  h2 = (TH1D*)file2->Get("AnaAfterHlt/nVtx");
  h3 = (TH1D*)file3->Get("AnaAfterHlt/nVtx");
  h1->SetName("hnvtx8tev");
  h2->SetName("hnvtx13tev25ns");
  h3->SetName("hnvtx13tev50ns");

  h1->Scale(1.0/h1->Integral());
  h2->Scale(1.0/h2->Integral());
  h3->Scale(1.0/h3->Integral());

  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(8);

  lg = new TLegend(0.6,0.7,0.9,0.9);
  lg->SetName("lg_nvtx");
  lg->AddEntry(h1, "8 TeV  50ns bx", "l");
  lg->AddEntry(h2, "13 TeV PU20bx25", "l");
  lg->AddEntry(h3, "13 TeV PU 14S", "l");

  plots->Clear();
  h2->Draw();
  h1->Draw("same");
  h3->Draw("same");
  lg->Draw("same");
  plots->Print(plotfile1);
  plots->Clear();


  // end
  sprintf(name, "%s]", plotfile1);
  plots->Print(name);
  sprintf(name, ".! ps2pdf %s %s", plotfile1, plotfile2);
  gROOT->ProcessLine(name);


}
