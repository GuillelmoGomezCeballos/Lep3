void HMuMu() {

  TCut cut3("Lep1pdgId+Lep2pdgId==0&&abs(Lep1pdgId-Lep2pdgId)==26&&abs(mMiss-95)<14&&acol>100&&ptMiss>11&&abs(pMiss*ctMiss)<54&&nJets>1&&max(Jet1fEle+Jet1fGamma,Jet2fEle+Jet2fGamma)<0.92");
  TFile *h = TFile::Open("HMuMu.root");
  TH1F* higgsbb = new TH1F("higgsbb","higgsbb",12,114.,138.);
  TH1F* higgsvis = (TH1F*)(higgsbb->Clone("higgsvis"));
  HMuMuTreeProducer_HMuMuAnalyzer->Project("higgsvis","ZMass",cut3);
  higgsvis->Scale(1.);

  TFile *z = TFile::Open("ZZMuMu.root");
  TH1F* zzvis = (TH1F*)(higgsbb->Clone("zzvis"));
  HMuMuTreeProducer_HMuMuAnalyzer->Project("zzvis","ZMass",cut3);
  zzvis->Scale(0.87);

  TFile *m = TFile::Open("llXMUMU.root");
  TH1F* mmvis = (TH1F*)(higgsbb->Clone("mmvis"));
  LLXTreeProducer_LLXAnalyzer->Project("mmvis","ZMass",cut3);
  mmvis->Scale(2.2*4);

  TFile *e = TFile::Open("llXELEELE.root");
  TH1F* eevis = (TH1F*)(higgsbb->Clone("eevis"));
  LLXTreeProducer_LLXAnalyzer->Project("eevis","ZMass",cut3);
  eevis->Scale(2.2*4);

  TFile *t = TFile::Open("llXTAUTAU.root");
  TH1F* ttvis = (TH1F*)(higgsbb->Clone("ttvis"));
  LLXTreeProducer_LLXAnalyzer->Project("ttvis","ZMass",cut3);
  ttvis->Scale(2.2*4);

  TFile *q = TFile::Open("llXQQ.root");
  TH1F* qqvis = (TH1F*)(higgsbb->Clone("qqvis"));
  LLXTreeProducer_LLXAnalyzer->Project("qqvis","ZMass",cut3);
  qqvis->Scale(4.0*4);

  TFile *w = TFile::Open("llXWW.root");
  TH1F* wwvis = (TH1F*)(higgsbb->Clone("wwvis"));
  LLXTreeProducer_LLXAnalyzer->Project("wwvis","ZMass",cut3);
  wwvis->Scale(4.0*4);

  TFile *wenu = TFile::Open("llXWENU.root");
  TH1F* wenuvis = (TH1F*)(higgsbb->Clone("wenuvis"));
  LLXTreeProducer_LLXAnalyzer->Project("wenuvis","ZMass",cut3);
  wenuvis->Scale(0.70*4);

  TFile *zee = TFile::Open("llXZEE.root");
  TH1F* zeevis = (TH1F*)(higgsbb->Clone("zeevis"));
  LLXTreeProducer_LLXAnalyzer->Project("zeevis","ZMass",cut3);
  zeevis->Scale(1.*4);

  TFile *znnb = TFile::Open("llXZNNB.root");
  TH1F* znnbvis = (TH1F*)(higgsbb->Clone("znnbvis"));
  LLXTreeProducer_LLXAnalyzer->Project("znnbvis","ZMass",cut3);
  znnbvis->Scale(0.83*4);


  TCanvas* c1 = new TCanvas();
  mmvis->Add(ttvis);
  eevis->Add(mmvis);
  qqvis->Add(eevis);
  zeevis->Add(qqvis);
  wenuvis->Add(zeevis);
  znnbvis->Add(wenuvis);
  wwvis->Add(znnbvis);
  zzvis->Add(wwvis);
  higgsvis->Add(zzvis);

  TF1 *fzvis = new TF1("fzvis","[0]+[1]*x+[2]*x*x+[3]*x*x*x",114.,138.);
  fzvis->SetParameters(0.,0.,0.,0.);
  zzvis->Fit("fzvis","EL","",114,138);  
  double p0z = fzvis->GetParameter(0);
  double p1z = fzvis->GetParameter(1);
  double p2z = fzvis->GetParameter(2);
  double p3z = fzvis->GetParameter(3);
  std::cout << p0z << " " << p1z << " " << p2z << " " << p3z << std::endl;

  TF1 *fhvis = new TF1("fhvis","39624.7-923.32*x+7.23078*x*x-0.0189505*x*x*x+[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",114.,138.);
  //TF1 *fhvis = new TF1("fhvis","60269.9-1394.01*x+10.8067*x*x-0.0279997*x*x*x+[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",114.,138.);
  fhvis->SetParameters(50.,125.,1.);
  higgsvis->Fit("fhvis","EL","",114,138);  
  //higgsvis->Fit("fhvis","EL","",110,150);  
  //higgsvis->Fit("fhvis","E","",110,150);  

  fzvis->SetLineStyle(13);
  fhvis->SetLineColor(6);
  higgsvis->SetMarkerStyle(21);
  higgsvis->SetMarkerColor(2);
  higgsvis->SetMarkerSize(1.2);
  higgsvis->SetLineColor(2);
  higgsvis->SetLineWidth(4);

  zzvis->SetLineColor(4);
  zzvis->SetLineWidth(4);
  zzvis->SetFillColor(4);
  zzvis->SetFillStyle(3013);

  znnbvis->SetLineColor(1);
  znnbvis->SetLineWidth(4);

  wwvis->SetLineColor(6);
  wwvis->SetLineWidth(4);

  eevis->SetLineColor(3);
  eevis->SetLineWidth(4);

  higgsvis->SetTitle( "Higgs -> mu+ mu-" );
  higgsvis->SetXTitle( "Higgs mass (GeV)" );
  higgsvis->SetYTitle( "Events / 2 GeV" );
  higgsvis->GetYaxis()->SetTitleOffset(1.4);
  higgsvis->GetYaxis()->SetLabelSize(0.045);
  higgsvis->GetXaxis()->SetLabelSize(0.045);
  higgsvis->SetStats(0);
  higgsvis->SetMaximum(320);
  higgsvis->SetMinimum(100);

  higgsvis->Draw("erro");
  fhvis->Draw("same");
  zzvis->Draw("same");
  fzvis->Draw("same");
  //wwvis->Draw("same");
  //znnbvis->Draw("same");
  //wenuvis->Draw("same");
  //zeevis->Draw("same");
  //qqvis->Draw("same");
  //eevis->Draw("same");

  TLegend *leg0=new TLegend(0.25,0.75,0.46,0.85);
  leg0->AddEntry( higgsvis, "Signal", "l0");
  leg0->AddEntry( zzvis, "All backgrounds", "lf");
  //leg0->AddEntry( zzvis, "ZZ", "l");
  //leg0->AddEntry( wwvis, "WW", "l");
  //leg0->AddEntry( znnbvis, "Zvv,Zee,Wev", "l");
  //leg0->AddEntry( eevis, "l+l-", "l");
  leg0->SetTextSize(0.03);
  leg0->Draw();

  TText *text = new TText(60,850,"L = 500 fb-1");
  //text->Draw("same");
  TPaveText* cmslumi = new TPaveText(0.6, 0.78, 0.85, 0.88, "NDC");
  cmslumi->AddText("LEP3, 500 fb^{-1},  #sqrt{s}=240 GeV");
  cmslumi->SetBorderSize(   0 );
  cmslumi->SetFillStyle(    0 );
  cmslumi->SetTextAlign(   12 );
  cmslumi->SetTextSize ( 0.03 );
  cmslumi->SetTextColor(    1 );
  cmslumi->SetTextFont (   62 );
  cmslumi->Draw("same");

  TPaveText* cmsdet = new TPaveText(0.65, 0.73, 0.85, 0.83, "NDC");
  cmsdet->AddText("Four detectors");
  cmsdet->SetBorderSize(   0 );
  cmsdet->SetFillStyle(    0 );
  cmsdet->SetTextAlign(   12 );
  cmsdet->SetTextSize ( 0.03 );
  cmsdet->SetTextColor(    1 );
  cmsdet->SetTextFont (   62 );
  cmsdet->Draw("same");

  TText *cms = new TText(130,325,"CMS Simulation");
  cms->Draw("same");

  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SaveAs("HMuMuvis.png");
  gPad->SaveAs("HMuMuvis.pdf");

}

Double_t CrystalBall(Double_t *x,Double_t *par) {

//Crystal ball function for signal, parameters are 0:alpha,1:n,2:mean,3:sigma,4:normalization;

  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);

  if (t >= -absAlpha) {
    return par[4]*exp(-0.5*t*t) + par[5] + t*par[6]+t*t*par[7] + t*t*t*par[8]; 
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha;

    return par[4]*(a/TMath::Power(b - t, par[1])) + par[5] + t*par[6] + t*t*par[7] + t*t*t*par[8];
  }
}
