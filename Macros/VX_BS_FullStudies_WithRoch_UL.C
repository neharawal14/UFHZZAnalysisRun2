#include "/afs/cern.ch/work/f/ferrico/private/HZZ_mass/CMSSW_10_2_18/src/UFHZZAnalysisRun2/Macros/Higgs_Mass_setup.h"

void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TString x_name, bool LogY, int modification = 0);
TH1F* RecoursiveFit(TString histo_name, std::vector<double> binning, TH2F* h2, bool Mean = 0);
TH1F* RecoursiveUnbinnedFit(TString histo_name, std::vector<double> binning, RooRealVar* rv_Res, RooDataSet* Data_Res, TString save, bool Mean = 0);
TH2F* RecoursiveFit3D(TString histo_name, std::vector<double> binning_x, std::vector<double> binning_y, TH3F* h2, bool Mean = 1);
TH2F* MassFit(TString histo_name, std::vector<double> binning_x, std::vector<double> binning_y, TH3F* h2, bool GEN);
std::vector<float> SingleMassFit(TString histo_name, TH1F* h2, TString save_name);
TH2F* MassScale(TString histo_name, TH2F* h1,  TH2F* h2);
void Draw_TH2F(TH2F* h1, TString nome_canvas, TString save, bool colz = 1);
void Draw_TH2F_FinalScale(TH2F* h1, TString nome_canvas, TString save, bool colz = 1);
float dataMC_SF(float pt, float eta, TH2F* h2);
TH2F* Draw_pTRes_Jake(TH2F* A, TH2F* B, TString nome_canvas, TString save, bool LogY, bool impro);

TLegend *legend_reco_gen = new TLegend(0.75,0.75,0.9,0.9);
TLegend *legend_year = new TLegend(0.1,0.7,0.2,0.9);
TLegend *legend_yearMass = new TLegend(0.8,0.7,0.9,0.9);

bool muon;

using namespace std;

void VX_BS_FullStudies_WithRoch_UL(TString year="2018", TString resonance="DY", bool deltaR=0, bool ScaleDone=0){


	
//	std::string DATAPATH;
//	RoccoR  *calibrator;



// 	if(year == "2016")
// 		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/HZZ_mass/CMSSW_10_2_18/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/final_HZZ_muon_SF_2016RunB2H_legacy_newLoose_newIso_paper.root";
// 	else if(year == "2017")
// 		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/HZZ_mass/CMSSW_10_2_18/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/final_HZZ_muon_SF_2017_newLooseIso_mupogSysts_paper.root";
// 	else
// 		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/HZZ_mass/CMSSW_10_2_18/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/final_HZZ_muon_SF_2018RunA2D_ER_newLoose_newIso_paper.root";

//     edm::FileInPath mu_scalefacFileInPath(DATAPATH);
//     TFile *fMuScalFac = TFile::Open(mu_scalefacFileInPath.fullPath().c_str());
//     hMuScaleFac = (TH2F*)fMuScalFac->Get("FINAL");
//     hMuScaleFacUnc = (TH2F*)fMuScalFac->Get("ERROR");
	
	gStyle->SetOptStat(0);
	gROOT->Reset();
	gROOT->SetBatch();
	
	std::vector<TString> GenMcData;
	GenMcData.clear();
	GenMcData.push_back("GEN_");
	GenMcData.push_back("RECO_");

	std::vector<TString> Samples;
	Samples.clear();
	Samples.push_back("UL_");
	Samples.push_back("Run3_");

	std::vector<TString> reco_type;
	reco_type.clear();
	reco_type.push_back("Roch_");
	reco_type.push_back("Roch_VX+BS_");

	std::vector<TString> leptons;
	leptons.clear();
	leptons.push_back("Inclusive_");
//	leptons.push_back("Negative_");
//	leptons.push_back("Positive_");

	pT_bins.clear();
	pT_bins.push_back(5);
	pT_bins.push_back(10);
	pT_bins.push_back(15);
	pT_bins.push_back(20);
	pT_bins.push_back(30);
	pT_bins.push_back(40);
	pT_bins.push_back(50);
	pT_bins.push_back(60);
	pT_bins.push_back(100);
// 	pT_bins.push_back(200);

	eta_bins.clear();
	eta_bins.push_back(0);
// 	eta_bins.push_back(0.4);
// 	eta_bins.push_back(0.8);
// 	eta_bins.push_back(1.2);
// 	eta_bins.push_back(1.6);
// 	eta_bins.push_back(2.0);
// 	eta_bins.push_back(2.4);
	eta_bins.push_back(0.9);
	eta_bins.push_back(1.4);
	eta_bins.push_back(2.4);

	phi_bins.clear();
	for(int i = 0; i < 20; i++)
		phi_bins.push_back(-3.14 + 0.314*i);
	phi_bins.push_back(3.14);

	std::vector<Double_t> d0_bins;
	d0_bins.clear();
	d0_bins.push_back(-0.01);
	d0_bins.push_back(-0.007);
	d0_bins.push_back(-0.005);
	d0_bins.push_back(-0.0035);
	d0_bins.push_back(-0.002);
	d0_bins.push_back(-0.001);
	d0_bins.push_back(0.0);
	d0_bins.push_back(0.001);
	d0_bins.push_back(0.002);
	d0_bins.push_back(0.0035);
	d0_bins.push_back(0.005);
	d0_bins.push_back(0.007);
	d0_bins.push_back(0.01);
	
	std::vector<Double_t> pTRes_bin;
	pTRes_bin.clear();
// 	for(int i = 0; i < 2000; i++)
// 		pTRes_bin.push_back(-0.25 + 0.00025*i);
// 	pTRes_bin.push_back(0.25);
	for(int i = 0; i < 1000; i++)
		pTRes_bin.push_back(-0.2 + 0.0004*i);
	pTRes_bin.push_back(0.2);

	std::vector<Double_t> scale_bins;
	scale_bins.clear();
	for(int i = 0; i < 200; i++)
		scale_bins.push_back(-0.5 + 0.005*i);
	scale_bins.push_back(0.5);

	std::vector<Double_t> mass_bins;
	mass_bins.clear();
// 	for(int i = 0; i < 2000; i++)
// 		mass_bins.push_back(81 + 0.01*i);
// 	mass_bins.push_back(101);
	for(int i = 0; i < 480; i++)
		mass_bins.push_back(60 + 0.125*i);
	mass_bins.push_back(120);

	std::vector<Double_t> deltapT_bins;
	deltapT_bins.clear();
	for(int i = 0; i < 400; i++)
		deltapT_bins.push_back(-0.05 + i*0.00025);
	deltapT_bins.push_back(0.05);
	
	std::vector<Double_t> vtx_bins;
	vtx_bins.clear();
	vtx_bins.push_back(0);
	if(year != "2018") vtx_bins.push_back(11);
	if(year != "2018") vtx_bins.push_back(14);
	if(year != "2018") 	vtx_bins.push_back(16);
	vtx_bins.push_back(17);
	if(year != "2018") 	vtx_bins.push_back(18);
	if(year != "2018") 	vtx_bins.push_back(19);
	vtx_bins.push_back(20);
	if(year != "2018") 	vtx_bins.push_back(21);
	vtx_bins.push_back(22);
	if(year != "2018") 	vtx_bins.push_back(23);
	vtx_bins.push_back(24);
	vtx_bins.push_back(25);
	vtx_bins.push_back(26);
	vtx_bins.push_back(27);
	vtx_bins.push_back(28);
	vtx_bins.push_back(30);
	vtx_bins.push_back(32);
	vtx_bins.push_back(34);
	if(year != "2018")
	vtx_bins.push_back(37);
	else{
		vtx_bins.push_back(36);
		vtx_bins.push_back(38);
	}
	vtx_bins.push_back(40);
	if(year != "2018")		
		vtx_bins.push_back(45);
	else{
		vtx_bins.push_back(42);
		vtx_bins.push_back(44);
		vtx_bins.push_back(46);
		vtx_bins.push_back(48);
		vtx_bins.push_back(50);
	}
	if(year == "2017") vtx_bins.push_back(50);
	if(year == "2017") vtx_bins.push_back(55);
	vtx_bins.push_back(60);
	if(year == "2017") vtx_bins.push_back(70);


	TH1F* PtFromTrackerDistribution[2][2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* EtaFromTrackerDistribution[2][2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* PhiFromTrackerDistribution[2][2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) x (incl x neg x pos)

	TH2F* MuonPV_xVSy[2][3]; //(Run3 vs UL) x (incl x neg x pos)
	TH1F* MuonPV_z[2][3]; //(Run3 vs UL) x (incl x neg x pos)





	TH1F* PtDistribution[2][2][3][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* PtErrDistribution[2][2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* EtaDistribution[2][2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* PhiDistribution[2][2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* DeltaRDistribution[2][2][3]; //(Run3 vs UL) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* d0Distribution[2][2][3]; //(Run3 vs UL) x (Baseline - VX+BS) x (incl x neg x pos)
	TH1F* MassDistribution[2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) 
	TH1F* MassErrDistribution[2][2][2]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) 
	
	TH2F* PtResolution_vsPt[2][3][3][3];
	
	RooRealVar* rv_PtResolution_vsPt[pT_bins.size()-1];
	RooArgSet*  rastmp_PtResolution_vsPt[pT_bins.size()-1];
	RooDataSet* Data_PtResolution_vsPt[pT_bins.size()-1];
	TH1F* h_PtResolution_vsPt_Unbinned[pT_bins.size()-1];
	
	for(int i = 0; i < pT_bins.size()-1; i++){
		TString nome = Form("DeltaPt_%.0f_%.0f", pT_bins.at(i), pT_bins.at(i+1));
		rv_PtResolution_vsPt[i] = new RooRealVar(nome, nome, -0.2, 0.2);
		rastmp_PtResolution_vsPt[i] = new RooArgSet(*rv_PtResolution_vsPt[i]);
		Data_PtResolution_vsPt[i] = new RooDataSet("Data" + nome, "Data" + nome, *rastmp_PtResolution_vsPt[i]);			
		h_PtResolution_vsPt_Unbinned[i] = new TH1F("PtResVsPt", "PtResVsPt", pT_bins.size()-1, &pT_bins[0]);
	}

	
	TH2F* PtResolution_vsEta[2][3][3][3];
	TH2F* PtResolution_vsPhi[2][3][3][3];
	TH2F* PtResolution_vsVertex[2][3][3][3];

	TH1F* h_PtResolution_vsPt[2][3][3][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) 
	TH1F* h_PtResolution_vsEta[2][3][3][3]; 
	TH1F* h_PtResolution_vsPhi[2][3][3][3]; 
	TH1F* h_PtResolution_vsVertex[2][3][3][3];

	TH1F* h_PtScale_vsPt[2][2][3]; //(Run3 vs UL) x (GEN vs RECO) x (Baseline - VX+BS) 
	TH1F* h_PtScale_vsEta[2][2][3]; 
	TH1F* h_PtScale_vsPhi[2][2][3]; 
	TH1F* h_PtScale_vsVertex[2][2][3];
	
	TH2F* qd0vsDeltaPtoverPt[2][2][3];//(Run3 vs UL) x (Baseline - VX+BS) x (incl x neg x pos)
	
	std::vector<TString> whered0;
	whered0.push_back("d0Inclusive");		

	
	for(int s = 0; s < Samples.size(); s++){
		for(int gmd = 0; gmd < GenMcData.size(); gmd++){
			for(int r = 0; r < reco_type.size(); r++){
				for(int l = 0; l < leptons.size(); l++){
				
					histo_name = Samples.at(s) + GenMcData.at(gmd) + reco_type.at(r) + leptons.at(l);

// 					if(gmd == 0 && r == 0){
// 						MuonPV_xVSy[s][l] = new TH2F(histo_name + "MuonPV_xVSy", histo_name + "MuonPV_xVSy", 200, -0.2, 0.2, 200, -0.2, 0.2);
// 						MuonPV_z[s][l] = new TH1F(histo_name + "MuonPV_z", histo_name + "MuonPV_z", 200, -10, 10);
// 					}
// 
					PtDistribution[s][gmd][r][l] = new TH1F(histo_name + "Pt", histo_name + "Pt", 200, 0, 100);
// 					PtErrDistribution[s][gmd][r][l] = new TH1F(histo_name + "PtErr", histo_name + "PtErr", 200, 0, 5);
// 					EtaDistribution[s][gmd][r][l] = new TH1F(histo_name + "Eta", histo_name + "Eta", 48, 0, 2.4);
// 					PhiDistribution[s][gmd][r][l] = new TH1F(histo_name + "Phi", histo_name + "Phi", 64, -3.2, 3.2);
// 
// 					PtFromTrackerDistribution[s][gmd][r][l] = new TH1F(histo_name + "PtFromTracker", histo_name + "PtFromTracker", 200, 0, 100);
// 					EtaFromTrackerDistribution[s][gmd][r][l] = new TH1F(histo_name + "EtaFromTracker", histo_name + "EtaFromTracker", 48, 0, 2.4);
// 					PhiFromTrackerDistribution[s][gmd][r][l] = new TH1F(histo_name + "PhiFromTracker", histo_name + "PhiFromTracker", 64, -3.2, 3.2);
// 
					if(l == 0){
						MassDistribution[s][gmd][r] = new TH1F(histo_name + "Mass", histo_name + "Mass", 100, 60, 120);
// 						MassErrDistribution[s][gmd][r] = new TH1F(histo_name + "MassErr", histo_name + "MassErr", 100, 0, 5);
					}
					if(gmd == 1){
// 						DeltaRDistribution[s][r][l] = new TH1F(histo_name + "DeltaR", histo_name + "DeltaR", 100, 0, 0.05);
// 						d0Distribution[s][r][l] = new TH1F(histo_name + "DeltaR", histo_name + "DeltaR", 200, -0.01, 0.01);
						
						for(int d0 = 0; d0 < whered0.size(); d0++){						
							histo_name = "h2_Resolution_" + Samples.at(s) + reco_type.at(r) + leptons.at(l) + whered0.at(d0);
							PtResolution_vsPt[s][r][l][d0] = new TH2F(histo_name + "vsPt", histo_name + "vsPt", pT_bins.size()-1, &pT_bins[0], pTRes_bin.size()-1, &pTRes_bin[0]);
							PtResolution_vsEta[s][r][l][d0] = new TH2F(histo_name + "vsEta", histo_name + "vsEta", eta_bins.size()-1, &eta_bins[0], pTRes_bin.size()-1, &pTRes_bin[0]);
							PtResolution_vsPhi[s][r][l][d0] = new TH2F(histo_name + "vsPhi", histo_name + "vsPhi", phi_bins.size()-1, &phi_bins[0], pTRes_bin.size()-1, &pTRes_bin[0]);
							PtResolution_vsVertex[s][r][l][d0] = new TH2F(histo_name + "vsVtx", histo_name + "vsVtx", vtx_bins.size()-1, &vtx_bins[0], pTRes_bin.size()-1, &pTRes_bin[0]);
							histo_name = "Resolution_" + Samples.at(s) + reco_type.at(r) + leptons.at(l) + whered0.at(d0);
							h_PtResolution_vsPt[s][r][l][d0] = new TH1F(histo_name + "vsPt", histo_name + "vsPt", pT_bins.size()-1, &pT_bins[0]);
							h_PtResolution_vsEta[s][r][l][d0] = new TH1F(histo_name + "vsEta", histo_name + "vsEta", eta_bins.size()-1, &eta_bins[0]);
							h_PtResolution_vsPhi[s][r][l][d0] = new TH1F(histo_name + "vsPhi", histo_name + "vsPhi", phi_bins.size()-1, &phi_bins[0]);
							h_PtResolution_vsVertex[s][r][l][d0] = new TH1F(histo_name + "vsVtx", histo_name + "vsVtx", vtx_bins.size()-1, &vtx_bins[0]);
						}

// 						histo_name = "Scale_" + Samples.at(s) + reco_type.at(r) + leptons.at(l);
// 						h_PtScale_vsPt[s][r][l] = new TH1F(histo_name + "vsPt", histo_name + "vsPt", pT_bins.size()-1, &pT_bins[0]);
// 						h_PtScale_vsEta[s][r][l] = new TH1F(histo_name + "vsEta", histo_name + "vsEta", eta_bins.size()-1, &eta_bins[0]);
// 						h_PtScale_vsPhi[s][r][l] = new TH1F(histo_name + "vsPhi", histo_name + "vsPhi", phi_bins.size()-1, &phi_bins[0]);
// 						h_PtScale_vsVertex[s][r][l] = new TH1F(histo_name + "vsVtx", histo_name + "vsVtx", vtx_bins.size()-1, &vtx_bins[0]);
// 
// 						histo_name = "qd0vsDeltaPtoverPt" + Samples.at(s) + reco_type.at(r) + leptons.at(l);						
// 						qd0vsDeltaPtoverPt[s][r][l] = new TH2F(histo_name, histo_name, 200, -.05, .05, 200, -0.01, 0.01);
					}
				}
			}
		}	
	}

// 	TH2F* h_MassVSmass[3];// (GEN - MC - DATA) 
// 	TH1F* h_DeltaMass[3];//  (GEN - MC - DATA) 
// 	TH1F* h_PtCheck[3][2];// (GEN - MC - DATA) x (Baseline - VX+BS)
// 	TH2F* h_ptVSpt[3][2][3];// (GEN - MC - DATA) 
// 	TH2F* h_GenPtVSpt[2][2][3];	
// 	TH1F* h_EtaCheck[3][2];// (GEN - MC - DATA) x (Baseline - VX+BS)
// 	TH2F* h_etaVSeta[3][2];// (GEN - MC - DATA) 
// 	TH1F* h_PhiCheck[3][2];// (GEN - MC - DATA) x (Baseline - VX+BS)
// 	TH2F* h_phiVSphi[3][2];// (GEN - MC - DATA) 
// 	TH1F* h_d0Check[3][2]; // (GEN - MC - DATA) x (Baseline - VX+BS)
// 	TH2F* h_d0VSd0[3][2][3];// (GEN - MC - DATA) 
	
	float massMin, massMax, maxDR;
	
	if(resonance == "DY"){
		massMin = 60;
		massMax = 120;
		maxDR = 0.002;
		maxDR = 0.05;
	}
	if(resonance == "JPsi"){
		massMin = 2.9;
		massMax = 3.3;
		maxDR = 0.005;
	}
	if(resonance == "Upsilon"){
		massMin = 8.5;
		massMax = 11;
		maxDR = 0.005;
	}

	if(resonance == "JpsiUpsilon"){
		massMin = 2.9;
		massMax = 11;
		maxDR = 0.005;
	}
	
	TFile* input_scale;
	TFile* new_scale;	
		
// 	for(int s = 0; s < 3; s++){
	for(int s = 1; s < 2; s++){
// 	for(int s = 0; s < 2; s++){
	
// 		if(deltaR) directory = "./Production_10_2_18_RecoGen_pt_Width_check_DeltaR";
// 		else directory = "./Production_10_2_18_RecoGen_pt_Width_check";		
// 		directory = "./Production_10_6_12_UL_Summer19_VX+BS_FullStudies";
// 		directory = "./Production_10_6_12_UL_onlyMC_notScaled";
		directory = "./Run3_VXBS_studies";
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		directory += "/" + year;
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		directory += "/" + resonance;
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		if(s == 0) directory += "/NoScaled";
		else if(s == 1)  directory += "/Scaled";
		else directory += "/SecondScaled";
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		directory += "/";

			
		int JPU = 2;
		if(resonance == "JpsiUpsilon") JPU = 4;
			
		for(int h = 0; h < JPU; h++){
/*		
if(h == 0){		
	if(year == "20160")
		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/UL_HZZ/CMSSW_10_6_12/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2016aUL.txt";
	else if(year == "20165")
		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/UL_HZZ/CMSSW_10_6_12/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2016bUL.txt";
	else if(year == "2017")
		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/UL_HZZ/CMSSW_10_6_12/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2017UL.txt";
	else
		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/UL_HZZ/CMSSW_10_6_12/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2018UL.txt";

	calibrator = new RoccoR(DATAPATH); 	
}
else{		
	if(year == "20160" || year == "20165")
		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/HZZ_mass/CMSSW_10_2_18/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2016.txt";
	else if(year == "2017")
		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/HZZ_mass/CMSSW_10_2_18/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2017.txt";
	else
		DATAPATH = "/afs/cern.ch/work/f/ferrico/private/HZZ_mass/CMSSW_10_2_18/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2018.txt";

	calibrator = new RoccoR(DATAPATH); 	
}
*/		
		
			int count_genFromReco = 0;
			TString filename;
			if(h == 0){
				filename = "/eos/user/f/ferrico/www/Chenguang/DY_muons_UL_20_s5m0/DY_madgraph_" + year + "_skimmed.root";
			}
			else
				filename = "DY_madgraph_skimmed_Run3_cancella_2018.root"; 
	
			std::cout<<h<<"\t"<<filename<<"\t"<<h<<std::endl;
//			std::cout<<h<<"\t"<<DATAPATH<<std::endl;

			 _file0 = new TFile(filename);
	
			if(_file0){
					tree = (TTree*)_file0->Get("passedEvents");		
			}
			else std::cout<<"ERROR could not find the file"<<std::endl;
		
			tree->SetBranchStatus("*",0);
			tree->SetBranchStatus("NInt", 1);
			tree->SetBranchStatus("NVtx", 1);

			tree->SetBranchStatus("massZ", 1);
			tree->SetBranchStatus("massZ_FSR", 1);
			tree->SetBranchStatus("massZ_single_BS", 1);
// 			tree->SetBranchStatus("massZ_single_BS_FSR", 1);
			tree->SetBranchStatus("massZ_vtx_BS", 1);
			tree->SetBranchStatus("massZ_vtx_BS_FSR", 1);

			tree->SetBranchStatus("massZErr", 1);
			tree->SetBranchStatus("massZErr_FSR", 1);
// 			tree->SetBranchStatus("massZErr_single_BS", 1);
// 			tree->SetBranchStatus("massZErr_single_BS_FSR", 1);
			tree->SetBranchStatus("massZErr_vtx_BS", 1);
			tree->SetBranchStatus("massZErr_vtx_BS_FSR", 1);		       	

			tree->SetBranchStatus("pT1_FromMuonBestTrack", 1);
			tree->SetBranchStatus("pT2_FromMuonBestTrack", 1);
			tree->SetBranchStatus("eta1_FromMuonBestTrack", 1);
			tree->SetBranchStatus("eta2_FromMuonBestTrack", 1);
			tree->SetBranchStatus("phi1_FromMuonBestTrack", 1);
			tree->SetBranchStatus("phi2_FromMuonBestTrack", 1);

			tree->SetBranchStatus("muonPV_x1", 1);
			tree->SetBranchStatus("muonPV_x2", 1);
			tree->SetBranchStatus("muonPV_y1", 1);
			tree->SetBranchStatus("muonPV_y2", 1);
			tree->SetBranchStatus("muonPV_z1", 1);
			tree->SetBranchStatus("muonPV_z2", 1);

			tree->SetBranchStatus("pT1_genFromReco", 1);
			tree->SetBranchStatus("pT2_genFromReco", 1);
			tree->SetBranchStatus("Tracker1", 1);
			tree->SetBranchStatus("Tracker2", 1);
			tree->SetBranchStatus("pT1", 1);
			tree->SetBranchStatus("pterr1", 1);
			tree->SetBranchStatus("eta1", 1);
			tree->SetBranchStatus("phi1", 1);
			tree->SetBranchStatus("m1", 1);
			tree->SetBranchStatus("d0BS1", 1);
			tree->SetBranchStatus("Id1", 1);
			tree->SetBranchStatus("pT2", 1);
			tree->SetBranchStatus("pterr2", 1);
			tree->SetBranchStatus("eta2", 1);
			tree->SetBranchStatus("phi2", 1);
			tree->SetBranchStatus("m2", 1);
			tree->SetBranchStatus("Id2", 1);
			tree->SetBranchStatus("d0BS2", 1);
			tree->SetBranchStatus("pT_FSR1", 1);
			tree->SetBranchStatus("eta_FSR1", 1);
			tree->SetBranchStatus("phi_FSR1", 1);
			tree->SetBranchStatus("m_FSR1", 1);
			tree->SetBranchStatus("pT_FSR2", 1);
			tree->SetBranchStatus("eta_FSR2", 1);
			tree->SetBranchStatus("phi_FSR2", 1);
			tree->SetBranchStatus("m_FSR2", 1);

			tree->SetBranchStatus("single_BS_pT1", 1);
// 			tree->SetBranchStatus("pterr1_single", 1);
			tree->SetBranchStatus("single_BS_eta1", 1);
			tree->SetBranchStatus("single_BS_phi1", 1);
// 			tree->SetBranchStatus("single_BS_m1", 1);
			tree->SetBranchStatus("single_BS_pT2", 1);
// 			tree->SetBranchStatus("pterr2_single", 1);
			tree->SetBranchStatus("single_BS_eta2", 1);
			tree->SetBranchStatus("single_BS_phi2", 1);
// 			tree->SetBranchStatus("single_BS_m2", 1);
// 			tree->SetBranchStatus("single_BS_pT_FSR1", 1);
// 			tree->SetBranchStatus("single_BS_eta_FSR1", 1);
// 			tree->SetBranchStatus("single_BS_phi_FSR1", 1);
// 			tree->SetBranchStatus("single_BS_m_FSR1", 1);
// 			tree->SetBranchStatus("single_BS_pT_FSR2", 1);
// 			tree->SetBranchStatus("single_BS_eta_FSR2", 1);
// 			tree->SetBranchStatus("single_BS_phi_FSR2", 1);
// 			tree->SetBranchStatus("single_BS_m_FSR2", 1);

			tree->SetBranchStatus("vtx_BS_pT1", 1);
			tree->SetBranchStatus("pterr1_VX_BS", 1);
			tree->SetBranchStatus("vtx_BS_eta1", 1);
			tree->SetBranchStatus("vtx_BS_phi1", 1);
			tree->SetBranchStatus("vtx_BS_m1", 1);
			tree->SetBranchStatus("d0BS_vtx_BS1", 1);
			tree->SetBranchStatus("vtx_BS_pT2", 1);
			tree->SetBranchStatus("pterr2_VX_BS", 1);
			tree->SetBranchStatus("vtx_BS_eta2", 1);
			tree->SetBranchStatus("vtx_BS_phi2", 1);
			tree->SetBranchStatus("vtx_BS_m2", 1);
			tree->SetBranchStatus("d0BS_vtx_BS2", 1);
			tree->SetBranchStatus("vtx_BS_pT_FSR1", 1);
			tree->SetBranchStatus("vtx_BS_eta_FSR1", 1);
			tree->SetBranchStatus("vtx_BS_phi_FSR1", 1);
			tree->SetBranchStatus("vtx_BS_m_FSR1", 1);
			tree->SetBranchStatus("vtx_BS_pT_FSR2", 1);
			tree->SetBranchStatus("vtx_BS_eta_FSR2", 1);
			tree->SetBranchStatus("vtx_BS_phi_FSR2", 1);
			tree->SetBranchStatus("vtx_BS_m_FSR2", 1);

// 			if(h == 0 || h == 2){		
			tree->SetBranchStatus("weight", 1);
			tree->SetBranchStatus("GENmass2l", 1);
			tree->SetBranchStatus("genLep_pt1", 1);
			tree->SetBranchStatus("genLep_eta1", 1);
			tree->SetBranchStatus("genLep_phi1", 1);
			tree->SetBranchStatus("genLep_pt2", 1);
			tree->SetBranchStatus("genLep_eta2", 1);
			tree->SetBranchStatus("genLep_phi2", 1);
// 			}


			tree->SetBranchAddress("NInt", &nInt);
			tree->SetBranchAddress("NVtx", &nVtx);

			tree->SetBranchAddress("massZ", &massZ);
			tree->SetBranchAddress("massZ_FSR", &massZ_FSR);
						// 			tree->SetBranchAddress("massZ_single_BS", &massZ);
// 			tree->SetBranchAddress("massZ_single_BS_FSR", &massZ_single_BS_FSR);
			tree->SetBranchAddress("massZ_vtx_BS", &massZ_vtx_BS);
			tree->SetBranchAddress("massZ_vtx_BS_FSR", &massZ_vtx_BS_FSR);

			tree->SetBranchAddress("massZErr", &massZErr);
			tree->SetBranchAddress("massZErr_FSR", &massZErr_FSR);
// 			tree->SetBranchAddress("massZErr_single_BS", &massZErr_single_BS);
// 			tree->SetBranchAddress("massZErr_single_BS_FSR", &massZErr_single_BS_FSR);
			tree->SetBranchAddress("massZErr_vtx_BS", &massZErr_vtx_BS);
			tree->SetBranchAddress("massZErr_vtx_BS_FSR", &massZErr_vtx_BS_FSR);
				
			tree->SetBranchAddress("pT1_FromMuonBestTrack", &pT1_FromMuonBestTrack); 
			tree->SetBranchAddress("pT2_FromMuonBestTrack", &pT2_FromMuonBestTrack);
			tree->SetBranchAddress("eta1_FromMuonBestTrack", &eta1_FromMuonBestTrack); 
			tree->SetBranchAddress("eta2_FromMuonBestTrack", &eta2_FromMuonBestTrack);
			tree->SetBranchAddress("phi1_FromMuonBestTrack", &phi1_FromMuonBestTrack); 
			tree->SetBranchAddress("phi2_FromMuonBestTrack", &phi2_FromMuonBestTrack);

			tree->SetBranchAddress("muonPV_x1", &muonPV_x1);
			tree->SetBranchAddress("muonPV_x2", &muonPV_x2);
			tree->SetBranchAddress("muonPV_y1", &muonPV_y1);
			tree->SetBranchAddress("muonPV_y2", &muonPV_y2);
			tree->SetBranchAddress("muonPV_z1", &muonPV_z1);
			tree->SetBranchAddress("muonPV_z2", &muonPV_z2);

			tree->SetBranchAddress("pT1_genFromReco", &pT1_genFromReco); 
			tree->SetBranchAddress("pT2_genFromReco", &pT2_genFromReco);
			tree->SetBranchAddress("Tracker1", &Tracker1); 
			tree->SetBranchAddress("Tracker2", &Tracker2);
 
			tree->SetBranchAddress("pT1", &pT1); 
			tree->SetBranchAddress("pterr1", &pterr1); 
			tree->SetBranchAddress("eta1", &eta1);   
			tree->SetBranchAddress("phi1", &phi1); 
			tree->SetBranchAddress("m1", &m2); 
			tree->SetBranchAddress("Id1", &Id1); 
			tree->SetBranchAddress("d0BS1", &d0BS1); 		  
			tree->SetBranchAddress("pT2", &pT2);   
			tree->SetBranchAddress("pterr2", &pterr2); 
			tree->SetBranchAddress("eta2", &eta2);   
			tree->SetBranchAddress("phi2", &phi2);   
			tree->SetBranchAddress("m1", &m2);   
			tree->SetBranchAddress("Id2", &Id2);   
			tree->SetBranchAddress("d0BS2", &d0BS2); 		  
			tree->SetBranchAddress("pT_FSR1", &pT_FSR1); 
			tree->SetBranchAddress("eta_FSR1", &eta_FSR1);   
			tree->SetBranchAddress("phi_FSR1", &phi_FSR1);   
			tree->SetBranchAddress("pT_FSR2", &pT_FSR2);   
			tree->SetBranchAddress("eta_FSR2", &eta_FSR2);   
			tree->SetBranchAddress("phi_FSR2", &phi_FSR2);   				
// 			tree->SetBranchAddress("vtx_pT1", &vtx_pT1); 
// 			tree->SetBranchAddress("pterr1_VX", &pterr1_VX); 
// 			tree->SetBranchAddress("vtx_eta1", &vtx_eta1);   
// 			tree->SetBranchAddress("vtx_phi1", &vtx_phi1);   
// 			tree->SetBranchAddress("vtx_pT2", &vtx_pT2);   
// 			tree->SetBranchAddress("pterr2_VX", &pterr2_VX); 
// 			tree->SetBranchAddress("vtx_eta2", &vtx_eta2);   
// 			tree->SetBranchAddress("vtx_phi2", &vtx_phi2);   
						// 			tree->SetBranchAddress("single_BS_pT1", &pT1); 
						// // 			tree->SetBranchAddress("pterr1_single", &pterr1_single); 
						// 			tree->SetBranchAddress("single_BS_eta1", &eta1);   
						// 			tree->SetBranchAddress("single_BS_phi1", &phi1);   
						// 			tree->SetBranchAddress("single_BS_pT2", &pT2);   
						// // 			tree->SetBranchAddress("pterr2_single", &pterr2_single); 
						// 			tree->SetBranchAddress("single_BS_eta2", &eta2);   
						// 			tree->SetBranchAddress("single_BS_phi2", &phi2);   
		
			tree->SetBranchAddress("vtx_BS_pT1", &vtx_BS_pT1); 
			tree->SetBranchAddress("pterr1_VX_BS", &pterr1_VX_BS); 
			tree->SetBranchAddress("vtx_BS_eta1", &vtx_BS_eta1);   
			tree->SetBranchAddress("vtx_BS_phi1", &vtx_BS_phi1);   
			tree->SetBranchAddress("vtx_BS_m1", &vtx_BS_m1);   
			tree->SetBranchAddress("d0BS_vtx_BS1", &d0BS_vtx_BS1);   
			tree->SetBranchAddress("vtx_BS_pT2", &vtx_BS_pT2);   
			tree->SetBranchAddress("pterr2_VX_BS", &pterr2_VX_BS); 
			tree->SetBranchAddress("vtx_BS_eta2", &vtx_BS_eta2);   
			tree->SetBranchAddress("vtx_BS_phi2", &vtx_BS_phi2);   
			tree->SetBranchAddress("vtx_BS_m2", &vtx_BS_m2);   
			tree->SetBranchAddress("d0BS_vtx_BS2", &d0BS_vtx_BS2);   

// 			tree->SetBranchAddress("vtx_pT_FSR1", &vtx_pT_FSR1); 
// 			tree->SetBranchAddress("vtx_eta_FSR1", &vtx_eta_FSR1);   
// 			tree->SetBranchAddress("vtx_phi_FSR1", &vtx_phi_FSR1);   
// 			tree->SetBranchAddress("vtx_pT_FSR2", &vtx_pT_FSR2);   
// 			tree->SetBranchAddress("vtx_eta_FSR2", &vtx_eta_FSR2);   
// 			tree->SetBranchAddress("vtx_phi_FSR2", &vtx_phi_FSR2);   
// 			tree->SetBranchAddress("single_BS_pT_FSR1", &single_BS_pT_FSR1); 
// 			tree->SetBranchAddress("single_BS_eta_FSR1", &single_BS_eta_FSR1);   
// 			tree->SetBranchAddress("single_BS_phi_FSR1", &single_BS_phi_FSR1);   
// 			tree->SetBranchAddress("single_BS_pT_FSR2", &single_BS_pT_FSR2);   
// 			tree->SetBranchAddress("single_BS_eta_FSR2", &single_BS_eta_FSR2);   
// 			tree->SetBranchAddress("single_BS_phi_FSR2", &single_BS_phi_FSR2);   

			tree->SetBranchAddress("vtx_BS_pT_FSR1", &vtx_BS_pT_FSR1); 
			tree->SetBranchAddress("vtx_BS_eta_FSR1", &vtx_BS_eta_FSR1);   
			tree->SetBranchAddress("vtx_BS_phi_FSR1", &vtx_BS_phi_FSR1);   
			tree->SetBranchAddress("vtx_BS_m_FSR1", &vtx_BS_m_FSR1);   
			tree->SetBranchAddress("vtx_BS_pT_FSR2", &vtx_BS_pT_FSR2);   
			tree->SetBranchAddress("vtx_BS_eta_FSR2", &vtx_BS_eta_FSR2);   
			tree->SetBranchAddress("vtx_BS_phi_FSR2", &vtx_BS_phi_FSR2);   
			tree->SetBranchAddress("vtx_BS_m_FSR2", &vtx_BS_m_FSR2);   

// 			if(h == 0 || h == 2){
			tree->SetBranchAddress("weight", &weight);   
			tree->SetBranchAddress("GENmass2l", &GENmass2l);   
			tree->SetBranchAddress("genLep_pt1", &genLep_pt1);   
			tree->SetBranchAddress("genLep_eta1", &genLep_eta1);       	
			tree->SetBranchAddress("genLep_phi1", &genLep_phi1);       	
			tree->SetBranchAddress("genLep_pt2", &genLep_pt2);   
			tree->SetBranchAddress("genLep_eta2", &genLep_eta2);  
			tree->SetBranchAddress("genLep_phi2", &genLep_phi2);  
// 			}

			Long64_t nentries = tree->GetEntries();

			std::cout<<nentries<<std::endl;	
		
			float deltaEta, deltaPhi, DR, DR1, DR2;
// 			if(deltaR) maxDR = 0.1;
// 			else 
			maxDR = 1000;
			int kpt, keta;

			TRandom3 rand;
			double u1; 
			bool FirstSecond; 
		
			int max_event = nentries;
// 			max_event = (int)nentries/2;
// 			max_event = (int)nentries/5;
// 			max_event = (int)nentries/10;
// 			max_event = (int)nentries/20;
// 			max_event = (int)nentries/50;
// 			max_event = (int)nentries/100;
// 			max_event = (int)nentries/1000;
// 			max_event = 100;
// 			max_event = 20;

			if(h == 0) max_event = (int)nentries / (int) 100;
			
			TH2F* tmp[2][3];
			if(ScaleDone){
				std::cout<<GenMcData.at(h+1) + "muN"<<std::endl;
				tmp[0][0] = (TH2F*)input_scale->Get(GenMcData.at(h+1) + "muN_lowPU");
				tmp[1][0] = (TH2F*)input_scale->Get(GenMcData.at(h+1) + "muP_lowPU");

				tmp[0][1] = (TH2F*)input_scale->Get(GenMcData.at(h+1) + "muN_mediumPU");
				tmp[1][1] = (TH2F*)input_scale->Get(GenMcData.at(h+1) + "muP_mediumPU");

				tmp[0][2] = (TH2F*)input_scale->Get(GenMcData.at(h+1) + "muN_highPU");
				tmp[1][2] = (TH2F*)input_scale->Get(GenMcData.at(h+1) + "muP_highPU");
			}		
				
// 			if(h == 0) max_event = nentries;
			
		
			for(int entry = 0; entry < max_event; entry++){

				tree->GetEntry(entry);  
			

				if(entry % 1000000 == 0)       
					std::cout<<entry<<" --- Dentro il TREE --- "<<year<<std::endl;  

// 				std::cout<<weight<<"\t";
				
// 				if(h == 0){
// 					if(year == "20160") weight *= 19000 * 7181 * 1. / (187485108);
// 					if(year == "20165") weight *= 16500 * 7181 * 1. / (95063845);
// 					if(year == "2017") weight *= 41855.2 * 7181 * 1. / (202549488);
// 					if(year == "2018") weight *= 58756.41 * 7181 * 1. / (203972117);
// 				}
// 				if(h == 1){
// 					if(year == "20160" || year == "20165") weight *= 35505.55 * 6225.4 * 1.3 / (121613296);
// 					if(year == "2017") weight *= 41855.2 * 6225.4 * 1.3 / (182359906);
// 					if(year == "2018") weight *= 58756.41 * 6225.4 * 1.4 / (193215674);
// 				}
				weight = 1;
// 				std::cout<<weight<<std::endl;


				float SF_n_First = 1;
				float SF_p_First = 1;				
				float SF_n_Second = 1;
				float SF_p_Second = 1;				

// 				std::cout<<vtx_BS_pT1<<"\t"<<fabs(eta1)<<"\t"<<phi1<<"\t"<<pT1<<std::endl;
// 				std::cout<<vtx_BS_pT2<<"\t"<<fabs(eta2)<<"\t"<<phi2<<"\t"<<pT2<<std::endl;
				
// 				TRandom3 rand;                                                                                                                                                                                                             
// 		       rand.SetSeed(abs(static_cast<int>(sin(mu.phi())*100000))); 
				if(s > 0){
					/*
					if(pT1_genFromReco < 0 || pT2_genFromReco < 0){
						count_genFromReco++;
						TRandom3 rand, rand2;                                                                                                                                                                                                             
						rand.SetSeed(abs(static_cast<int>(sin(vtx_BS_phi1)*100000)));                                                                                          
						rand2.SetSeed(abs(static_cast<int>(sin(vtx_BS_phi2)*100000)));                                                                                          
						double u1, u2;
						u1 = rand.Uniform(1.);
						u2 = rand2.Uniform(1.);
						if(h == 0){
							if(Id1 == 13){
								SF_n_First = calibrator->kSpreadMC(-1, vtx_BS_pT1, vtx_BS_eta1, vtx_BS_phi1, Tracker1, u1);
								SF_p_First = calibrator->kSpreadMC(1, vtx_BS_pT2, vtx_BS_eta2, vtx_BS_phi2, Tracker2, u2);
							}
							else{
								SF_n_First = calibrator->kSpreadMC(-1, vtx_BS_pT2, vtx_BS_eta2, vtx_BS_phi2, Tracker2, u2);
								SF_p_First = calibrator->kSpreadMC(1, vtx_BS_pT1, vtx_BS_eta1, vtx_BS_phi1, Tracker1, u1);
							}
						}
					}							
					else{						
                                                if(h == 0){

							if(Id1 == 13){
								SF_n_First = calibrator->kSpreadMC(-1, vtx_BS_pT1, vtx_BS_eta1, vtx_BS_phi1, pT1_genFromReco);
								SF_p_First = calibrator->kSpreadMC(1, vtx_BS_pT2, vtx_BS_eta2, vtx_BS_phi2, pT2_genFromReco);
							}
							else{
								SF_n_First = calibrator->kSpreadMC(-1, vtx_BS_pT2, vtx_BS_eta2, vtx_BS_phi2, pT2_genFromReco);
								SF_p_First = calibrator->kSpreadMC(1, vtx_BS_pT1, vtx_BS_eta1, vtx_BS_phi1, pT1_genFromReco);
							}
						}
					}*/

// 					std::cout<<SF_n_First<<"\t"<<SF_p_First<<std::endl;
			
					if(Id1 == 13){
						vtx_BS_pT1 = vtx_BS_pT1 * SF_n_First;
						vtx_BS_pT2 = vtx_BS_pT2 * SF_p_First;					
					}
					else{
						vtx_BS_pT1 = vtx_BS_pT1 * SF_p_First;
						vtx_BS_pT2 = vtx_BS_pT2 * SF_n_First;					
					}

// 					std::cout<<vtx_BS_pT1<<"\t"<<SF_n<<"\t"<<pT1<<std::endl;
// 					std::cout<<vtx_BS_pT2<<"\t"<<SF_p<<"\t"<<pT2<<std::endl;
					TLorentzVector lep1, lep2;
					lep1.SetPtEtaPhiM(vtx_BS_pT1, vtx_BS_eta1, vtx_BS_phi1, vtx_BS_m1);
					lep2.SetPtEtaPhiM(vtx_BS_pT2, vtx_BS_eta2, vtx_BS_phi2, vtx_BS_m2);
					massZ_vtx_BS = (lep1+lep2).M();					
				}	


/*
				// scale for pT from bestMuon
				if(s > 0){
					if(h == 0){
						if(pT1_genFromReco < 0 || pT2_genFromReco < 0){
							TRandom3 rand, rand2;                                                                                                                                                                                                             
							rand.SetSeed(abs(static_cast<int>(sin(phi1)*100000)));                                                                                          
							rand2.SetSeed(abs(static_cast<int>(sin(phi2)*100000)));                                                                                          
							double u1, u2;
							u1 = rand.Uniform(1.);
							u2 = rand2.Uniform(1.);
							if(Id1 == 13){
								SF_n_First = calibrator->kSpreadMC(-1, pT1_FromMuonBestTrack, eta1, phi1, Tracker1, u1);
								SF_p_First = calibrator->kSpreadMC(1, pT2_FromMuonBestTrack, eta2, phi2, Tracker2, u2);
							}
							else{
								SF_n_First = calibrator->kSpreadMC(-1, pT2_FromMuonBestTrack, eta2, phi2, Tracker2, u2);
								SF_p_First = calibrator->kSpreadMC(1, pT1_FromMuonBestTrack, eta1, phi1, Tracker1, u1);
							}
						}							
						else{
							if(Id1 == 13){
								SF_n_First = calibrator->kSpreadMC(-1, pT1_FromMuonBestTrack, eta1, phi1, pT1_genFromReco);
								SF_p_First = calibrator->kSpreadMC(1, pT2_FromMuonBestTrack, eta2, phi2, pT2_genFromReco);
							}
							else{
								SF_n_First = calibrator->kSpreadMC(-1, pT2_FromMuonBestTrack, eta2, phi2, pT2_genFromReco);
								SF_p_First = calibrator->kSpreadMC(1, pT1_FromMuonBestTrack, eta1, phi1, pT1_genFromReco);
							}
						}
					}
					else{
						if(Id1 == 13){
							SF_n_First = calibrator->kScaleDT(-1, pT1_FromMuonBestTrack, eta1, phi1);
							SF_p_First = calibrator->kScaleDT(1, pT2_FromMuonBestTrack, eta2, phi2);
						}
						else{
							SF_n_First = calibrator->kScaleDT(-1, pT1_FromMuonBestTrack, eta2, phi2);
							SF_p_First = calibrator->kScaleDT(1, vtx_BS_pT1, eta1, phi1);
						}
					}

// 					std::cout<<SF_n_First<<"\t"<<SF_p_First<<std::endl;
			
					if(Id1 == 13){
						pT1_FromMuonBestTrack = pT1_FromMuonBestTrack * SF_n_First;
						pT2_FromMuonBestTrack = pT2_FromMuonBestTrack * SF_p_First;					
					}
					else{
						pT1_FromMuonBestTrack = pT1_FromMuonBestTrack * SF_p_First;
						pT2_FromMuonBestTrack = pT2_FromMuonBestTrack * SF_n_First;					
					}
				}
				// scale for pT from bestMuon
*/
// 				std::cout<<h<<"\t"<<massZ<<"\t"<<massZ_vtx_BS<<"\t"<<weight<<std::endl;


				float d0_BS_charge1, d0_BS_charge2;
				if(Id1 == 13){
					d0_BS_charge1 = -d0BS1;
					d0_BS_charge2 = d0BS1;
				}
				else{
					d0_BS_charge1 = d0BS1;
					d0_BS_charge2 = -d0BS1;
				}

 				std::vector<float> mass;
 				mass.push_back(massZ);
 				mass.push_back(massZ_vtx_BS);
 				
 				MassDistribution[h][0][0]->Fill(GENmass2l);
 				MassDistribution[h][1][0]->Fill(massZ, weight);
 				MassDistribution[h][1][1]->Fill(massZ_vtx_BS, weight);
// 
//  				MassErrDistribution[h][1][0]->Fill(massZErr);
//  				MassErrDistribution[h][1][1]->Fill(massZErr_vtx_BS);
 				
/*
				// mass distribution;
 				if(h == 0 || h == 2){
	 				MassDistribution[0][0]->Fill(BW_mean_PDG);
	 				MassDistribution[1][0]->Fill(massZ, weight);
	 				MassDistribution[1][1]->Fill(massZ_vtx_BS, weight);	 				
//  					if(pT1 > 20 && pT2 > 20){
 					if(pT1 > 0 && pT2 > 0){
						h_MassVSmass[1]->Fill(massZ, massZ_vtx_BS, weight);
						h_etaVSeta[1][0]->Fill(eta1, vtx_BS_eta1, weight);
						h_etaVSeta[1][0]->Fill(eta2, vtx_BS_eta2, weight);
						h_phiVSphi[1][0]->Fill(phi1, vtx_BS_phi1, weight);
						h_phiVSphi[1][0]->Fill(phi2, vtx_BS_phi2, weight);
						float deltaM = (massZ_vtx_BS - massZ)/massZ;
						h_DeltaMass[1]->Fill(deltaM, weight);					
						if(fabs(massZ_vtx_BS - massZ) > 5 && massZ > 85 && massZ < 95){

							h_ptVSpt[1][1]->Fill(pT1, vtx_BS_pT1, weight);
							h_ptVSpt[1][1]->Fill(pT2, vtx_BS_pT2, weight);
							h_etaVSeta[1][1]->Fill(eta1, vtx_BS_eta1, weight);
							h_etaVSeta[1][1]->Fill(eta2, vtx_BS_eta2, weight);
							h_phiVSphi[1][1]->Fill(phi1, vtx_BS_phi1, weight);
							h_phiVSphi[1][1]->Fill(phi2, vtx_BS_phi2, weight);
							h_d0VSd0[1][1]->Fill(d0BS1, d0BS_vtx_BS1, weight);
							h_d0VSd0[1][1]->Fill(d0BS2, d0BS_vtx_BS2, weight);

							h_PtCheck[1][0]->Fill(pT1, weight);
							h_PtCheck[1][0]->Fill(pT2, weight);
							h_EtaCheck[1][0]->Fill(eta1, weight);
							h_EtaCheck[1][0]->Fill(eta2, weight);
							h_PhiCheck[1][0]->Fill(phi1, weight);
							h_PhiCheck[1][0]->Fill(phi2, weight);
							h_d0Check[1][0]->Fill(d0BS1, weight);
							h_d0Check[1][0]->Fill(d0BS2, weight);
							h_PtCheck[1][1]->Fill(vtx_BS_pT1, weight);
							h_PtCheck[1][1]->Fill(vtx_BS_pT2, weight);
							h_EtaCheck[1][1]->Fill(vtx_BS_eta1, weight);
							h_EtaCheck[1][1]->Fill(vtx_BS_eta2, weight);
							h_PhiCheck[1][1]->Fill(vtx_BS_phi1, weight);
							h_PhiCheck[1][1]->Fill(vtx_BS_phi2, weight);
							h_d0Check[1][1]->Fill(d0BS_vtx_BS1, weight);
							h_d0Check[1][1]->Fill(d0BS_vtx_BS2, weight);
						}
					}
		 		}
		 		else{
	 				MassDistribution[2][0]->Fill(massZ);
	 				MassDistribution[2][1]->Fill(massZ_vtx_BS);
//  					if(pT1 > 20 && pT2 > 20){
 					if(pT1 > 0 && pT2 > 0){
						h_MassVSmass[2]->Fill(massZ, massZ_vtx_BS);
						h_etaVSeta[2][0]->Fill(eta1, vtx_BS_eta1);
						h_etaVSeta[2][0]->Fill(eta2, vtx_BS_eta2);
						h_phiVSphi[2][0]->Fill(phi1, vtx_BS_phi1);
						h_phiVSphi[2][0]->Fill(phi2, vtx_BS_phi2);
						float deltaM = (massZ_vtx_BS - massZ)/massZ;
						if(massZ > 60 && massZ < 120) h_DeltaMass[2]->Fill(deltaM);
	// 					if(fabs(deltaM) < 0.01){
						if(fabs(massZ_vtx_BS - massZ) > 5 && massZ > 85 && massZ < 95){

							h_ptVSpt[2][1]->Fill(pT1, vtx_BS_pT1);
							h_ptVSpt[2][1]->Fill(pT2, vtx_BS_pT2);
							h_etaVSeta[2][1]->Fill(eta1, vtx_BS_eta1);
							h_etaVSeta[2][1]->Fill(eta2, vtx_BS_eta2);
							h_phiVSphi[2][1]->Fill(phi1, vtx_BS_phi1);
							h_phiVSphi[2][1]->Fill(phi2, vtx_BS_phi2);
							h_d0VSd0[2][1]->Fill(d0BS1, d0BS_vtx_BS1);
							h_d0VSd0[2][1]->Fill(d0BS2, d0BS_vtx_BS2);

							h_PtCheck[2][0]->Fill(pT1);
							h_PtCheck[2][0]->Fill(pT2);
							h_EtaCheck[2][0]->Fill(eta1);
							h_EtaCheck[2][0]->Fill(eta2);
							h_PhiCheck[2][0]->Fill(phi1);
							h_PhiCheck[2][0]->Fill(phi2);
							h_d0Check[2][0]->Fill(d0BS1);
							h_d0Check[2][0]->Fill(d0BS2);
							h_PtCheck[2][1]->Fill(vtx_BS_pT1);
							h_PtCheck[2][1]->Fill(vtx_BS_pT2);
							h_EtaCheck[2][1]->Fill(vtx_BS_eta1);
							h_EtaCheck[2][1]->Fill(vtx_BS_eta2);
							h_PhiCheck[2][1]->Fill(vtx_BS_phi1);
							h_PhiCheck[2][1]->Fill(vtx_BS_phi2);
							h_d0Check[2][1]->Fill(d0BS_vtx_BS1);
							h_d0Check[2][1]->Fill(d0BS_vtx_BS2);
						}
					}
		 		}
*/


 				std::vector<float> PT1;
 				PT1.push_back(pT1);
 				PT1.push_back(vtx_BS_pT1);

 				std::vector<float> PT2;
 				PT2.push_back(pT2);
 				PT2.push_back(vtx_BS_pT2);

 				std::vector<float> PT1err;
 				PT1err.push_back(pterr1);
 				PT1err.push_back(pterr1_VX_BS);

 				std::vector<float> PT2err;
 				PT2err.push_back(pterr2);
 				PT2err.push_back(pterr2_VX_BS);

 				std::vector<float> ETA1;
 				ETA1.push_back(fabs(eta1));
 				ETA1.push_back(fabs(vtx_BS_eta1));

 				std::vector<float> ETA2;
 				ETA2.push_back(fabs(eta2));
 				ETA2.push_back(fabs(vtx_BS_eta2));

 				std::vector<float> PHI1;
 				PHI1.push_back(phi1);
 				PHI1.push_back(vtx_BS_phi1);

 				std::vector<float> PHI2;
 				PHI2.push_back(phi2);
 				PHI2.push_back(vtx_BS_phi2);

 				std::vector<float> d01;
 				d01.push_back(d0BS1);
 				d01.push_back(d0BS_vtx_BS1);

 				std::vector<float> d02;
 				d02.push_back(d0BS2);
 				d02.push_back(d0BS_vtx_BS2);
				
 				for(int i = 0; i < reco_type.size(); i++){
 					 					
					if(mass.at(i) < massMin || mass.at(i) > massMax) continue;
					
// 					float max_DeltaPt = 0.05;
/*					
				 	if(h == 0 || h == 2){
					 	deltaEta = ETA1.at(i) - genLep_eta1;
						deltaPhi = PHI1.at(i) - genLep_phi1;
						DR1 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));

					 	deltaEta = ETA2.at(i) - genLep_eta2;
						deltaPhi = PHI2.at(i) - genLep_phi2;
						DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
					
						if(DR1 > maxDR || DR2 > maxDR) continue;
						
						DeltaRDistribution[1][0][i]->Fill(DR1, weight);
						DeltaRDistribution[1][0][i]->Fill(DR2, weight);						
					}
					else{
					 	deltaEta = ETA1.at(i) - genLep_eta1;
						deltaPhi = PHI1.at(i) - genLep_phi1;
						DR1 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));

					 	deltaEta = ETA2.at(i) - genLep_eta2;
						deltaPhi = PHI2.at(i) - genLep_phi2;
						DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
					
						if(DR1 > maxDR || DR2 > maxDR) continue;
						
						DeltaRDistribution[2][0][i]->Fill(DR1, weight);
						DeltaRDistribution[2][0][i]->Fill(DR2, weight);

					}
*/					
// 					std::cout<<genLep_pt1<<"\t"<<genLep_pt2<<std::endl;
// 					std::cout<<PT1.at(0)<<"\t"<<PT2.at(0)<<std::endl;
// 					std::cout<<PT1.at(1)<<"\t"<<PT2.at(1)<<std::endl;
// 					std::cout<<"\t\t"<<std::endl;
// 					int FirstPositive = -1;
// 					if(Id1 == 13) FirstPositive = 0;
// 					else FirstPositive = 1;
					
// 					deltaEta = ETA1.at(i) - genLep_eta1;
// 					deltaPhi = PHI1.at(i) - genLep_phi1;
// 					DR1 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
// 
// 					deltaEta = ETA2.at(i) - genLep_eta2;
// 					deltaPhi = PHI2.at(i) - genLep_phi2;
// 					DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
					
// 					if(DR1 > maxDR || DR2 > maxDR) continue; 

					
// 					PtDistribution[h][0][0][0]->Fill(genLep_pt1);
// 					PtDistribution[h][0][0][0]->Fill(genLep_pt2);
// 					
					PtDistribution[h][1][i][0]->Fill(PT1.at(i));
					PtDistribution[h][1][i][0]->Fill(PT2.at(i));
// 					PtErrDistribution[h][1][i][0]->Fill(PT1err.at(i));
// 					PtErrDistribution[h][1][i][0]->Fill(PT2err.at(i));
// 					EtaDistribution[h][1][i][0]->Fill(ETA1.at(i));
// 					EtaDistribution[h][1][i][0]->Fill(ETA2.at(i));
// 					PhiDistribution[h][1][i][0]->Fill(PHI2.at(i));
// 					PhiDistribution[h][1][i][0]->Fill(PHI1.at(i));

// 					if(i == 0){
// 						PtFromTrackerDistribution[h][1][i][0]->Fill(pT1_FromMuonBestTrack);
// 						PtFromTrackerDistribution[h][1][i][0]->Fill(pT2_FromMuonBestTrack);
// 						if(year == "2018"){
// 							EtaFromTrackerDistribution[h][1][i][0]->Fill(eta1_FromMuonBestTrack);
// 							EtaFromTrackerDistribution[h][1][i][0]->Fill(eta2_FromMuonBestTrack);
// 							PhiFromTrackerDistribution[h][1][i][0]->Fill(phi1_FromMuonBestTrack);
// 							PhiFromTrackerDistribution[h][1][i][0]->Fill(phi2_FromMuonBestTrack);
// 							MuonPV_xVSy[h][0]->Fill(muonPV_x1, muonPV_y1);
// 							MuonPV_xVSy[h][0]->Fill(muonPV_x2, muonPV_y2);
// 							MuonPV_z[h][0]->Fill(muonPV_z1);
// 							MuonPV_z[h][0]->Fill(muonPV_z2);
// 						}
// 					}


// 					DeltaRDistribution[h][i][0]->Fill(DR1);
// 					DeltaRDistribution[h][i][0]->Fill(DR2);
// 					d0Distribution[h][i][0]->Fill(d01.at(i));
// 					d0Distribution[h][i][0]->Fill(d02.at(i));
// 					
// 					qd0vsDeltaPtoverPt[h][i][0]->Fill((PT1.at(i) - genLep_pt1)/genLep_pt1, Id1*d01.at(i)/(13));
// 					qd0vsDeltaPtoverPt[h][i][0]->Fill((PT2.at(i) - genLep_pt2)/genLep_pt2, Id2*d02.at(i)/(13));
					
					PtResolution_vsPt[h][i][0][0]->Fill(PT1.at(i), -(genLep_pt1 - PT1.at(i))/genLep_pt1);
					PtResolution_vsPt[h][i][0][0]->Fill(PT2.at(i), -(genLep_pt2 - PT2.at(i))/genLep_pt2);
					PtResolution_vsEta[h][i][0][0]->Fill(ETA1.at(i), -(genLep_pt1 - PT1.at(i))/genLep_pt1);
					PtResolution_vsEta[h][i][0][0]->Fill(ETA2.at(i), -(genLep_pt2 - PT2.at(i))/genLep_pt2);
					PtResolution_vsPhi[h][i][0][0]->Fill(PHI1.at(i), -(genLep_pt1 - PT1.at(i))/genLep_pt1);
					PtResolution_vsPhi[h][i][0][0]->Fill(PHI2.at(i), -(genLep_pt2 - PT2.at(i))/genLep_pt2);
					PtResolution_vsVertex[h][i][0][0]->Fill(nVtx, -(genLep_pt1 - PT1.at(i))/genLep_pt1);
					PtResolution_vsVertex[h][i][0][0]->Fill(nVtx, -(genLep_pt2 - PT2.at(i))/genLep_pt2);
					
					if(i == 0){
						if(fabs(genLep_pt1 - PT1.at(i))/genLep_pt1 < 0.25){
							if(PT1.at(i) > 5 && PT1.at(i) < 20){
								rv_PtResolution_vsPt[0]->setVal(-(genLep_pt1 - PT1.at(i))/genLep_pt1);
								Data_PtResolution_vsPt[0]->add(*rastmp_PtResolution_vsPt[0]);
							}
							if(PT1.at(i) > 20 && PT1.at(i) < 30){
								rv_PtResolution_vsPt[1]->setVal(-(genLep_pt1 - PT1.at(i))/genLep_pt1);
								Data_PtResolution_vsPt[1]->add(*rastmp_PtResolution_vsPt[1]);
							}

							if(PT1.at(i) > 30 && PT1.at(i) < 40){
								rv_PtResolution_vsPt[2]->setVal(-(genLep_pt1 - PT1.at(i))/genLep_pt1);
								Data_PtResolution_vsPt[2]->add(*rastmp_PtResolution_vsPt[2]);
							}

							if(PT1.at(i) > 40 && PT1.at(i) < 50){
								rv_PtResolution_vsPt[3]->setVal(-(genLep_pt1 - PT1.at(i))/genLep_pt1);
								Data_PtResolution_vsPt[3]->add(*rastmp_PtResolution_vsPt[3]);
							}

							if(PT1.at(i) > 50 && PT1.at(i) < 60){
								rv_PtResolution_vsPt[4]->setVal(-(genLep_pt1 - PT1.at(i))/genLep_pt1);
								Data_PtResolution_vsPt[4]->add(*rastmp_PtResolution_vsPt[4]);
							}
							if(PT1.at(i) > 60 && PT1.at(i) < 100){
								rv_PtResolution_vsPt[5]->setVal(-(genLep_pt1 - PT1.at(i))/genLep_pt1);
								Data_PtResolution_vsPt[5]->add(*rastmp_PtResolution_vsPt[5]);
							}

							
						}

						if(fabs(genLep_pt2 - PT2.at(i))/genLep_pt2 < 0.25){
							if(PT2.at(i) > 5 && PT2.at(i) < 20){
								rv_PtResolution_vsPt[0]->setVal(-(genLep_pt2 - PT2.at(i))/genLep_pt2);
								Data_PtResolution_vsPt[0]->add(*rastmp_PtResolution_vsPt[0]);
							}
							if(PT2.at(i) > 20 && PT2.at(i) < 30){
								rv_PtResolution_vsPt[1]->setVal(-(genLep_pt2 - PT2.at(i))/genLep_pt2);
								Data_PtResolution_vsPt[1]->add(*rastmp_PtResolution_vsPt[1]);
							}

							if(PT2.at(i) > 30 && PT2.at(i) < 40){
								rv_PtResolution_vsPt[2]->setVal(-(genLep_pt2 - PT2.at(i))/genLep_pt2);
								Data_PtResolution_vsPt[2]->add(*rastmp_PtResolution_vsPt[2]);
							}

							if(PT2.at(i) > 40 && PT2.at(i) < 50){
								rv_PtResolution_vsPt[3]->setVal(-(genLep_pt2 - PT2.at(i))/genLep_pt2);
								Data_PtResolution_vsPt[3]->add(*rastmp_PtResolution_vsPt[3]);
							}

							if(PT2.at(i) > 50 && PT2.at(i) < 60){
								rv_PtResolution_vsPt[4]->setVal(-(genLep_pt2 - PT2.at(i))/genLep_pt2);
								Data_PtResolution_vsPt[4]->add(*rastmp_PtResolution_vsPt[4]);
							}
							if(PT2.at(i) > 60 && PT2.at(i) < 100){
								rv_PtResolution_vsPt[5]->setVal(-(genLep_pt2 - PT2.at(i))/genLep_pt2);
								Data_PtResolution_vsPt[5]->add(*rastmp_PtResolution_vsPt[5]);
							}
						}


					}


/*					
					if(!FirstPositive){
						PtDistribution[h][0][0][1]->Fill(genLep_pt1);
						PtDistribution[h][0][0][2]->Fill(genLep_pt2);

						PtDistribution[h][1][i][1]->Fill(PT1.at(i));
						PtDistribution[h][1][i][2]->Fill(PT2.at(i));
						PtErrDistribution[h][1][i][1]->Fill(PT1err.at(i));
						PtErrDistribution[h][1][i][2]->Fill(PT2err.at(i));
						EtaDistribution[h][1][i][1]->Fill(ETA1.at(i));
						EtaDistribution[h][1][i][2]->Fill(ETA2.at(i));
						PhiDistribution[h][1][i][1]->Fill(PHI1.at(i));
						PhiDistribution[h][1][i][2]->Fill(PHI2.at(i));
						
						if(i == 0){						
							PtFromTrackerDistribution[h][1][i][1]->Fill(pT1_FromMuonBestTrack);
							PtFromTrackerDistribution[h][1][i][2]->Fill(pT2_FromMuonBestTrack);
							if(year == "2018"){
								EtaFromTrackerDistribution[h][1][i][1]->Fill(eta1_FromMuonBestTrack);
								EtaFromTrackerDistribution[h][1][i][2]->Fill(eta2_FromMuonBestTrack);
								PhiFromTrackerDistribution[h][1][i][1]->Fill(phi1_FromMuonBestTrack);
								PhiFromTrackerDistribution[h][1][i][2]->Fill(phi2_FromMuonBestTrack);
								MuonPV_xVSy[h][1]->Fill(muonPV_x1, muonPV_y1);
								MuonPV_xVSy[h][2]->Fill(muonPV_x2, muonPV_y2);
								MuonPV_z[h][1]->Fill(muonPV_z1);
								MuonPV_z[h][2]->Fill(muonPV_z2);
							}
						}

						DeltaRDistribution[h][i][1]->Fill(DR1);
						DeltaRDistribution[h][i][2]->Fill(DR2);
						d0Distribution[h][i][1]->Fill(d01.at(i));
						d0Distribution[h][i][2]->Fill(d02.at(i));
						PtResolution_vsPt[h][i][1]->Fill(PT1.at(i), (genLep_pt1 - PT1.at(i))/genLep_pt1);
						PtResolution_vsPt[h][i][2]->Fill(PT2.at(i), (genLep_pt2 - PT2.at(i))/genLep_pt2);
						PtResolution_vsEta[h][i][1]->Fill(ETA1.at(i), (genLep_pt1 - PT1.at(i))/genLep_pt1);
						PtResolution_vsEta[h][i][2]->Fill(ETA2.at(i), (genLep_pt2 - PT2.at(i))/genLep_pt2);
						PtResolution_vsPhi[h][i][1]->Fill(PHI1.at(i), (genLep_pt1 - PT1.at(i))/genLep_pt1);
						PtResolution_vsPhi[h][i][2]->Fill(PHI2.at(i), (genLep_pt2 - PT2.at(i))/genLep_pt2);
						PtResolution_vsVertex[h][i][1]->Fill(nVtx, (genLep_pt1 - PT1.at(i))/genLep_pt1);
						PtResolution_vsVertex[h][i][2]->Fill(nVtx, (genLep_pt2 - PT2.at(i))/genLep_pt2);

						qd0vsDeltaPtoverPt[h][i][1]->Fill((PT1.at(i) - genLep_pt1)/genLep_pt1, Id1*d01.at(i)/(13));
						qd0vsDeltaPtoverPt[h][i][2]->Fill((PT2.at(i) - genLep_pt2)/genLep_pt2, Id2*d02.at(i)/(13));

					}
					else{
						PtDistribution[h][0][0][1]->Fill(genLep_pt2);
						PtDistribution[h][0][0][2]->Fill(genLep_pt1);

						PtDistribution[h][1][i][1]->Fill(PT2.at(i));
						PtDistribution[h][1][i][2]->Fill(PT1.at(i));
						PtErrDistribution[h][1][i][1]->Fill(PT2err.at(i));
						PtErrDistribution[h][1][i][2]->Fill(PT1err.at(i));
						EtaDistribution[h][1][i][1]->Fill(ETA2.at(i));
						EtaDistribution[h][1][i][2]->Fill(ETA1.at(i));
						PhiDistribution[h][1][i][1]->Fill(PHI2.at(i));
						PhiDistribution[h][1][i][2]->Fill(PHI1.at(i));

						if(i == 0){												
							PtFromTrackerDistribution[h][1][i][1]->Fill(pT2_FromMuonBestTrack);
							PtFromTrackerDistribution[h][1][i][2]->Fill(pT1_FromMuonBestTrack);
							if(year == "2018"){
								EtaFromTrackerDistribution[h][1][i][1]->Fill(eta2_FromMuonBestTrack);
								EtaFromTrackerDistribution[h][1][i][2]->Fill(eta1_FromMuonBestTrack);
								PhiFromTrackerDistribution[h][1][i][1]->Fill(phi2_FromMuonBestTrack);
								PhiFromTrackerDistribution[h][1][i][2]->Fill(phi1_FromMuonBestTrack);
								MuonPV_xVSy[h][1]->Fill(muonPV_x2, muonPV_y2);
								MuonPV_xVSy[h][2]->Fill(muonPV_x1, muonPV_y1);
								MuonPV_z[h][1]->Fill(muonPV_z2);
								MuonPV_z[h][2]->Fill(muonPV_z1);
							}
						}

						DeltaRDistribution[h][i][1]->Fill(DR2);
						DeltaRDistribution[h][i][2]->Fill(DR1);
						d0Distribution[h][i][1]->Fill(d02.at(i));
						d0Distribution[h][i][2]->Fill(d01.at(i));
						PtResolution_vsPt[h][i][1]->Fill(PT2.at(i), (genLep_pt2 - PT2.at(i))/genLep_pt2);
						PtResolution_vsPt[h][i][2]->Fill(PT1.at(i), (genLep_pt1 - PT1.at(i))/genLep_pt1);
						PtResolution_vsEta[h][i][1]->Fill(ETA2.at(i), (genLep_pt2 - PT2.at(i))/genLep_pt2);
						PtResolution_vsEta[h][i][2]->Fill(ETA1.at(i), (genLep_pt1 - PT1.at(i))/genLep_pt1);
						PtResolution_vsPhi[h][i][1]->Fill(PHI2.at(i), (genLep_pt2 - PT2.at(i))/genLep_pt2);
						PtResolution_vsPhi[h][i][2]->Fill(PHI1.at(i), (genLep_pt1 - PT1.at(i))/genLep_pt1);
						PtResolution_vsVertex[h][i][1]->Fill(nVtx, (genLep_pt2 - PT2.at(i))/genLep_pt2);
						PtResolution_vsVertex[h][i][2]->Fill(nVtx, (genLep_pt1 - PT1.at(i))/genLep_pt1);

						qd0vsDeltaPtoverPt[h][i][1]->Fill((PT2.at(i) - genLep_pt2)/genLep_pt2, Id2*d02.at(i)/(13));
						qd0vsDeltaPtoverPt[h][i][2]->Fill((PT1.at(i) - genLep_pt1)/genLep_pt1, Id1*d01.at(i)/(13));
					}										
*/
				} // for on Roch - VX+BS

				PT1.clear(); PT2.clear();
				ETA1.clear(); ETA2.clear();
				PHI1.clear(); PHI2.clear();
				d01.clear(); d02.clear();

			} // for on entries
		
			std::cout<<"CAZZO di genFromReco = "<<count_genFromReco<<std::endl;

		} // for on MC - DATA
		
	
		TString nome_canvas;

		save_nome = "Distribution";
		if(!h_blank) std::cout<<"h_blank proble."<<std::endl;
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf[", "", 100);

// 		nome_canvas = "GEN Pt inclusive";
// 		Draw(PtDistribution[1][0][0][0], PtDistribution[0][0][0][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
// 
// 		nome_canvas = "Pt from tracker";
// 		Draw(PtFromTrackerDistribution[1][1][0][0], PtFromTrackerDistribution[0][1][0][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
// 		nome_canvas = "Eta from tracker";
// 		Draw(EtaFromTrackerDistribution[1][1][0][0], EtaFromTrackerDistribution[0][1][0][0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0);
// 		nome_canvas = "Phi from tracker";
// 		Draw(PhiFromTrackerDistribution[1][1][0][0], PhiFromTrackerDistribution[0][1][0][0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0);
//  		nome_canvas = "MuonPV x vs y: " + Samples.at(0) + leptons.at(0);
// 		Draw_TH2F(MuonPV_xVSy[0][0], nome_canvas, directory + save_nome + ".pdf", 0);			
//  		nome_canvas = "MuonPV x vs y: " + Samples.at(1) + leptons.at(0);
// 		Draw_TH2F(MuonPV_xVSy[1][0], nome_canvas, directory + save_nome + ".pdf", 0);			
//  		nome_canvas = "MuonPV z: " + leptons.at(0);
// 		Draw(MuonPV_z[1][0], MuonPV_z[0][0],  nome_canvas, directory + save_nome + ".pdf", "cm", 0);

		nome_canvas = "Pt inclusive: Roch";
		Draw(PtDistribution[1][1][0][0], PtDistribution[0][1][0][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
		nome_canvas = "Pt inclusive: VX+BS";
		Draw(PtDistribution[1][1][1][0], PtDistribution[0][1][1][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
// 		nome_canvas = "PtErr inclusive: Roch";
// 		Draw(PtErrDistribution[1][1][0][0], PtErrDistribution[0][1][0][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
// 		nome_canvas = "PtErr inclusive: VX+BS";
// 		Draw(PtErrDistribution[1][1][1][0], PtErrDistribution[0][1][1][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
// 		nome_canvas = "Eta inclusive: Roch";
// 		Draw(EtaDistribution[1][1][0][0], EtaDistribution[0][1][0][0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0);
// 		nome_canvas = "Eta inclusive: VX+BS";
// 		Draw(EtaDistribution[1][1][1][0], EtaDistribution[0][1][1][0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0);
// 		nome_canvas = "Phi inclusive: Roch";
// 		Draw(PhiDistribution[1][1][0][0], PhiDistribution[0][1][0][0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0);
// 		nome_canvas = "Phi inclusive: VX+BS";
// 		Draw(PhiDistribution[1][1][1][0], PhiDistribution[0][1][1][0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0);
// 		nome_canvas = "DeltaR inclusive: Roch";
// 		Draw(DeltaRDistribution[1][0][0], DeltaRDistribution[0][0][0],  nome_canvas, directory + save_nome + ".pdf", "#DeltaR", 1);
// 		nome_canvas = "DeltaR inclusive: VX+BS";
// 		Draw(DeltaRDistribution[1][1][0], DeltaRDistribution[0][1][0],  nome_canvas, directory + save_nome + ".pdf", "#DeltaR", 1);
// 		nome_canvas = "d0 inclusive: Roch";
// 		Draw(d0Distribution[1][0][0], d0Distribution[0][0][0],  nome_canvas, directory + save_nome + ".pdf", "d0(BS)", 0);
// 		nome_canvas = "d0 inclusive: VX+BS";
// 		Draw(d0Distribution[1][1][0], d0Distribution[0][1][0],  nome_canvas, directory + save_nome + ".pdf", "d0(BS)", 0);
		nome_canvas = "Mass inclusive: Roch";
		Draw(MassDistribution[1][1][0], MassDistribution[0][1][0],  nome_canvas, directory + save_nome + ".pdf", "m_{2l} [GeV]", 0, -1);
		nome_canvas = "Mass inclusive: VX+BS";
		Draw(MassDistribution[1][1][1], MassDistribution[0][1][1],  nome_canvas, directory + save_nome + ".pdf", "m_{2l} [GeV]", 0, -1);
// 		nome_canvas = "MassErr inclusive: Roch";
// 		Draw(MassErrDistribution[1][1][0], MassErrDistribution[0][1][0],  nome_canvas, directory + save_nome + ".pdf", "m_{2l} [GeV]", 1);
// 		nome_canvas = "MassErr inclusive: VX+BS";
// 		Draw(MassErrDistribution[1][1][1], MassErrDistribution[0][1][1],  nome_canvas, directory + save_nome + ".pdf", "m_{2l} [GeV]", 1);
/*
		for(int l = 1; l < leptons.size(); l++){
			
			nome_canvas = leptons.at(l) + "GEN Pt inclusive";
			Draw(PtDistribution[1][0][0][l], PtDistribution[0][0][0][l],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);

			nome_canvas = leptons.at(l) + "Pt from tracker";
			Draw(PtFromTrackerDistribution[1][1][0][l], PtFromTrackerDistribution[0][1][0][l],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
			nome_canvas = leptons.at(l) + "Eta from tracker";
			Draw(EtaFromTrackerDistribution[1][1][0][l], EtaFromTrackerDistribution[0][1][0][l],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0);
			nome_canvas = leptons.at(l) + "Phi from tracker";
			Draw(PhiFromTrackerDistribution[1][1][0][l], PhiFromTrackerDistribution[0][1][0][l],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0);
	 		nome_canvas = leptons.at(l) + "MuonPV x vs y " + Samples.at(0);
			Draw_TH2F(MuonPV_xVSy[0][l], nome_canvas, directory + save_nome + ".pdf", 0);			
			nome_canvas = leptons.at(l) + "MuonPV x vs y " + Samples.at(1);
			Draw_TH2F(MuonPV_xVSy[1][l], nome_canvas, directory + save_nome + ".pdf", 0);			
			nome_canvas = leptons.at(l) + "MuonPV z";
			Draw(MuonPV_z[1][l], MuonPV_z[0][l],  nome_canvas, directory + save_nome + ".pdf", "cm", 0);

			nome_canvas = leptons.at(l) + "Pt inclusive: Roch";
			Draw(PtDistribution[1][1][0][l], PtDistribution[0][1][0][l],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
			nome_canvas = leptons.at(l) + "Pt inclusive: VX+BS";
			Draw(PtDistribution[1][1][1][l], PtDistribution[0][1][1][l],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
			nome_canvas = leptons.at(l) + "PtErr inclusive: Roch";
			Draw(PtErrDistribution[1][1][0][l], PtErrDistribution[0][1][0][l],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
			nome_canvas = leptons.at(l) + "PtErr inclusive: VX+BS";
			Draw(PtErrDistribution[1][1][1][l], PtErrDistribution[0][1][1][l],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 1);
			nome_canvas = leptons.at(l) + "Eta inclusive: Roch";
			Draw(EtaDistribution[1][1][0][l], EtaDistribution[0][1][0][l],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0);
			nome_canvas = leptons.at(l) + "Eta inclusive: VX+BS";
			Draw(EtaDistribution[1][1][1][l], EtaDistribution[0][1][1][l],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0);
			nome_canvas = leptons.at(l) + "Phi inclusive: Roch";
			Draw(PhiDistribution[1][1][0][l], PhiDistribution[0][1][0][l],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0);
			nome_canvas = leptons.at(l) + "Phi inclusive: VX+BS";
			Draw(PhiDistribution[1][1][1][l], PhiDistribution[0][1][1][l],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0);
			nome_canvas = leptons.at(l) + "DeltaR inclusive: Roch";
			Draw(DeltaRDistribution[1][0][l], DeltaRDistribution[0][0][l],  nome_canvas, directory + save_nome + ".pdf", "#DeltaR", 1);
			nome_canvas = leptons.at(l) + "DeltaR inclusive: VX+BS";
			Draw(DeltaRDistribution[1][1][l], DeltaRDistribution[0][1][l],  nome_canvas, directory + save_nome + ".pdf", "#DeltaR", 1);
			nome_canvas = leptons.at(l) + "d0 inclusive: Roch";
			Draw(d0Distribution[1][0][l], d0Distribution[0][0][l],  nome_canvas, directory + save_nome + ".pdf", "d0(BS)", 0);
			nome_canvas = leptons.at(l) + "d0 inclusive: VX+BS";
			Draw(d0Distribution[1][1][l], d0Distribution[0][1][l],  nome_canvas, directory + save_nome + ".pdf", "d0(BS)", 0);
		}	
*/
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf]", "", 100);
				
		
		save_nome = "UnbinnedResolution_" + Samples.at(s);
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf[", "", 100);		
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf", "", 100);		
		for(int i = 0; i < pT_bins.size()-1; i++){
			TString name = Form("PtRes_vsPt_Unbinned_%.0f", pT_bins.at(i));
			std::cout<<name<<std::endl;
			h_PtResolution_vsPt_Unbinned[i] = RecoursiveUnbinnedFit(name, pT_bins, rv_PtResolution_vsPt[i], Data_PtResolution_vsPt[i], directory + save_nome + ".pdf");		
		}
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf]", "", 100);


		std::vector<TString> list_histo;
		list_histo.push_back("Roch");
		list_histo.push_back("VX+BS");
                list_histo.push_back("VX+BS");

		// resolution		
 		for(int s = 0; s < Samples.size(); s++){
//		for(int s = 1; s < Samples.size(); s++){
			for(int r = 0; r < reco_type.size(); r++){
				for(int l = 0; l < leptons.size(); l++){
					for(int d0 = 0; d0 < whered0.size(); d0++){
						histo_name = "h2_Resolution_" + Samples.at(s) + reco_type.at(r) + leptons.at(l) + whered0.at(d0);
						TString name = histo_name + "PtRes_vsPt";
						h_PtResolution_vsPt[s][r][l][d0] = RecoursiveFit(name, pT_bins, PtResolution_vsPt[s][r][l][d0]);
						name = histo_name + "PtRes_vsEta";
						h_PtResolution_vsEta[s][r][l][d0] = RecoursiveFit(name, eta_bins, PtResolution_vsEta[s][r][l][d0]);
						name = histo_name + "PtRes_vsPhi";
						h_PtResolution_vsPhi[s][r][l][d0] = RecoursiveFit(name, phi_bins, PtResolution_vsPhi[s][r][l][d0]);
						name = histo_name + "PtRes_vsVtx";
						h_PtResolution_vsVertex[s][r][l][d0] = RecoursiveFit(name, vtx_bins, PtResolution_vsVertex[s][r][l][d0]);
					}
				}
			}		
			save_nome = "Resolution_" + Samples.at(s);
			Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf[", "", 100);
			for(int d0 = 0; d0 < whered0.size(); d0++){
				nome_canvas = "pT inclusive" + whered0.at(d0);
				DrawResolution(h_PtResolution_vsPt[s][0][0][d0], h_PtResolution_vsPt[s][1][0][d0], h_PtResolution_vsPt[s][1][0][d0], list_histo, nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]");

				nome_canvas = "#eta inclusive" + whered0.at(d0);
				DrawResolution(h_PtResolution_vsEta[s][0][0][d0], h_PtResolution_vsEta[s][1][0][d0], h_PtResolution_vsEta[s][1][0][d0], list_histo, nome_canvas, directory + save_nome + ".pdf", "#eta");

				nome_canvas = "#phi inclusive" + whered0.at(d0);
				DrawResolution(h_PtResolution_vsPhi[s][0][0][d0], h_PtResolution_vsPhi[s][1][0][d0], h_PtResolution_vsPhi[s][1][0][d0], list_histo, nome_canvas, directory + save_nome + ".pdf", "#phi");

				nome_canvas = "reco Vertex inclusive" + whered0.at(d0);
				DrawResolution(h_PtResolution_vsVertex[s][0][0][d0], h_PtResolution_vsVertex[s][1][0][d0], h_PtResolution_vsVertex[s][1][0][d0], list_histo, nome_canvas, directory + save_nome + ".pdf", "# vtx");
			}
			for(int d0 = 0; d0 < whered0.size(); d0++){
				for(int i = 0; i < leptons.size(); i++){
					nome_canvas = leptons.at(i) + "pT" + whered0.at(d0);
					Draw(h_PtResolution_vsPt[s][0][i][d0], h_PtResolution_vsPt[s][1][i][d0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 0, 2);
					nome_canvas = leptons.at(i) + "#eta" + whered0.at(d0);
					Draw(h_PtResolution_vsEta[s][0][i][d0], h_PtResolution_vsEta[s][1][i][d0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0, 2);
					nome_canvas = leptons.at(i) + "#phi" + whered0.at(d0);
					Draw(h_PtResolution_vsPhi[s][0][i][d0], h_PtResolution_vsPhi[s][1][i][d0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0, 2);
					nome_canvas = leptons.at(i) + "reco Vertex" + whered0.at(d0);
					Draw(h_PtResolution_vsVertex[s][0][i][d0], h_PtResolution_vsVertex[s][1][i][d0],  nome_canvas, directory + save_nome + ".pdf", "# vtx", 0, 2);
				}
			}
			Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf]", "", 100);
		}
		save_nome = "Resolution";
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf[", "", 100);
		for(int r = 0; r < reco_type.size(); r++){
			nome_canvas = "pT inclusive" + reco_type.at(r);
			Draw(h_PtResolution_vsPt[0][r][0][0], h_PtResolution_vsPt[1][r][0][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 0, -2);
			nome_canvas = "#eta inclusive" + reco_type.at(r);
			Draw(h_PtResolution_vsEta[0][r][0][0], h_PtResolution_vsEta[1][r][0][0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0, -2);
			nome_canvas = "#phi inclusive" + reco_type.at(r);
			Draw(h_PtResolution_vsPhi[0][r][0][0], h_PtResolution_vsPhi[1][r][0][0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0, -2);
			nome_canvas = "reco Vertex inclusive" + reco_type.at(r);
			Draw(h_PtResolution_vsVertex[0][r][0][0], h_PtResolution_vsVertex[1][r][0][0],  nome_canvas, directory + save_nome + ".pdf", "# vtx", 0, -2);
			for(int i = 1; i < leptons.size(); i++){
				nome_canvas = leptons.at(i) + "pT" + reco_type.at(r);
				Draw(h_PtResolution_vsPt[0][r][i][0], h_PtResolution_vsPt[1][r][i][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 0, -2);
				nome_canvas = leptons.at(i) + "#eta" + reco_type.at(r);
				Draw(h_PtResolution_vsEta[0][r][i][0], h_PtResolution_vsEta[1][r][i][0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0, -2);
				nome_canvas = leptons.at(i) + "#phi" + reco_type.at(r);
				Draw(h_PtResolution_vsPhi[0][r][i][0], h_PtResolution_vsPhi[1][r][i][0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0, -2);
				nome_canvas = leptons.at(i) + "reco Vertex" + reco_type.at(r);
				Draw(h_PtResolution_vsVertex[0][r][i][0], h_PtResolution_vsVertex[1][r][i][0],  nome_canvas, directory + save_nome + ".pdf", "# vtx", 0, -2);
			}
		}
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf]", "", 100);

		/*	
		save_nome = "qd0vsDeltaPtOverPt";
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf[", "", 100);
		for(int i = 0; i < leptons.size(); i++){
			for(int r = 0; r < reco_type.size(); r++){
				for(int s = 0; s < Samples.size(); s++){
					nome_canvas = Samples.at(s) + reco_type.at(r) + leptons.at(i);
					Draw_TH2F(qd0vsDeltaPtoverPt[s][r][i], nome_canvas, directory + save_nome + ".pdf", 0);			
				}
			}
		}	
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf]", "", 100);
			
		// scale
		for(int s = 0; s < Samples.size(); s++){
			for(int r = 0; r < reco_type.size(); r++){
				for(int l = 0; l < leptons.size(); l++){
					histo_name = "h2_Scale_" + Samples.at(s) + reco_type.at(r) + leptons.at(l);
					TString name = histo_name + "PtRes_vsPt";
					h_PtScale_vsPt[s][r][l] = RecoursiveFit(name, pT_bins, PtResolution_vsPt[s][r][l], 1);
					name = histo_name + "PtRes_vsEta";
					h_PtScale_vsEta[s][r][l] = RecoursiveFit(name, eta_bins, PtResolution_vsEta[s][r][l], 1);
					name = histo_name + "PtRes_vsPhi";
					h_PtScale_vsPhi[s][r][l] = RecoursiveFit(name, phi_bins, PtResolution_vsPhi[s][r][l], 1);
					name = histo_name + "PtRes_vsVtx";
					h_PtScale_vsVertex[s][r][l] = RecoursiveFit(name, vtx_bins, PtResolution_vsVertex[s][r][l], 1);
				}
			}		

			save_nome = "Scale_" + Samples.at(s);
			Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf[", "", 100);
			nome_canvas = "pT inclusive";
			Draw(h_PtScale_vsPt[s][0][0], h_PtScale_vsPt[s][1][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 0, 3);
			nome_canvas = "#eta inclusive";
			Draw(h_PtScale_vsEta[s][0][0], h_PtScale_vsEta[s][1][0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0, 3);
			nome_canvas = "#phi inclusive";
			Draw(h_PtScale_vsPhi[s][0][0], h_PtScale_vsPhi[s][1][0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0, 3);
			nome_canvas = "reco Vertex inclusive";
			Draw(h_PtScale_vsVertex[s][0][0], h_PtScale_vsVertex[s][1][0],  nome_canvas, directory + save_nome + ".pdf", "# vtx", 0, 3);
			for(int i = 1; i < leptons.size(); i++){
				nome_canvas = leptons.at(i) + "pT";
				Draw(h_PtScale_vsPt[s][0][i], h_PtScale_vsPt[s][1][i],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 0, 3);
				nome_canvas = leptons.at(i) + "#eta";
				Draw(h_PtScale_vsEta[s][0][i], h_PtScale_vsEta[s][1][i],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0, 3);
				nome_canvas = leptons.at(i) + "#phi";
				Draw(h_PtScale_vsPhi[s][0][i], h_PtScale_vsPhi[s][1][i],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0, 3);
				nome_canvas = leptons.at(i) + "reco Vertex";
				Draw(h_PtScale_vsVertex[s][0][i], h_PtScale_vsVertex[s][1][i],  nome_canvas, directory + save_nome + ".pdf", "# vtx", 0, 3);
			}
			Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf]", "", 100);
		}

		save_nome = "Scale";
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf[", "", 100);
		for(int r = 0; r < reco_type.size(); r++){
			nome_canvas = "pT inclusive" + reco_type.at(r);
			Draw(h_PtScale_vsPt[0][r][0], h_PtScale_vsPt[1][r][0],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 0, -3);
			nome_canvas = "#eta inclusive" + reco_type.at(r);
			Draw(h_PtScale_vsEta[0][r][0], h_PtScale_vsEta[1][r][0],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0, -3);
			nome_canvas = "#phi inclusive" + reco_type.at(r);
			Draw(h_PtScale_vsPhi[0][r][0], h_PtScale_vsPhi[1][r][0],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0, -3);
			nome_canvas = "reco Vertex inclusive" + reco_type.at(r);
			Draw(h_PtScale_vsVertex[0][r][0], h_PtScale_vsVertex[1][r][0],  nome_canvas, directory + save_nome + ".pdf", "# vtx", 0, -3);
			for(int i = 1; i < leptons.size(); i++){
				nome_canvas = leptons.at(i) + "pT" + reco_type.at(r);
				Draw(h_PtScale_vsPt[0][r][i], h_PtScale_vsPt[1][r][i],  nome_canvas, directory + save_nome + ".pdf", "p_{T} [GeV]", 0, -3);
				nome_canvas = leptons.at(i) + "#eta" + reco_type.at(r);
				Draw(h_PtScale_vsEta[0][r][i], h_PtScale_vsEta[1][r][i],  nome_canvas, directory + save_nome + ".pdf", "#eta", 0, -3);
				nome_canvas = leptons.at(i) + "#phi" + reco_type.at(r);
				Draw(h_PtScale_vsPhi[0][r][i], h_PtScale_vsPhi[1][r][i],  nome_canvas, directory + save_nome + ".pdf", "#phi", 0, -3);
				nome_canvas = leptons.at(i) + "reco Vertex" + reco_type.at(r);
				Draw(h_PtScale_vsVertex[0][r][i], h_PtScale_vsVertex[1][r][i],  nome_canvas, directory + save_nome + ".pdf", "# vtx", 0, -3);
			}
		}
		Draw(h_blank, h_blank, "blank", directory + save_nome + ".pdf]", "", 100);
		*/				
	} //scaled - no scaled

}



void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TString x_name, bool LogY, int modification = 0){

	std::vector<float> tmp_1;
	std::vector<float> tmp_2;
	
	if(modification == 1 || modification == 0){
		h1->Scale(1/h1->Integral());
		h2->Scale(1/h2->Integral());
	}

	if(modification == -1){
// 		tmp_1 = SingleMassFit("Data", h1, save);
		tmp_1 = SingleMassFit("Run3", h1, save);
		tmp_2 = SingleMassFit("UL", h2, save);		
	}		

	
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11;
   	pad11 = new TPad("pad1", "pad1", 0, 0.2, 1, 1);
	if(modification == 0)
	   	pad11 = new TPad("pad1", "pad1", 0, 0.1, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
//    	pad11->SetTicks();
   	
	Double_t min, max;

   	if(LogY){
   		pad11->SetLogy();
		min = Min(h1->GetMinimum(), h2->GetMinimum());		
		if(min == 0)
			min = 0.00001;
   	}
	else
		min = 0;
		
	if(modification == 3){
		min = Min(h1->GetMinimum(), h2->GetMinimum());		
		min = 1.1 * min;
	}
	
	max = Max(h1->GetMaximum(), h2->GetMaximum());		
	if(max > 0)
		max = max * 1.1;
	else
		max = max / 1.1;
	
	if(max == 0) max = 1;	
	
	h1->SetLineColor(kBlack);
	h2->SetLineColor(kRed);
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(22);
	h1->SetMarkerSize(0.75);
	h2->SetMarkerSize(0.75);
	h1->SetMarkerColor(kBlack);
	h2->SetMarkerColor(kRed);

	TH1F *ratio_2 = (TH1F*) h1->Clone();
	
	canvas->Update();

		
	if(modification < 1 && modification != -2 && modification != -3){

		THStack* MC = new THStack("UL", "UL");		
		h2->SetFillColor(kRed);
		MC->Add(h2, "HIST");
		h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
		h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);	
		MC->Draw();
// 		h2->Draw("Same B");
		h1->Draw("same E");

		h1->SetTitle(nome_canvas);
		MC->SetTitle(nome_canvas);
		h1->GetXaxis()->SetTitle(x_name);
		MC->GetXaxis()->SetTitle(x_name);

		if(modification == -1){
// 		if(modification == 11){
// 			TLatex* tex1 = new TLatex(0.125,0.65, "Data:");
			TLatex* tex1 = new TLatex(0.125,0.65, "Run3:");
			tex1->SetTextSize(0.04);
			tex1->SetNDC();
			tex1->Draw();
	
			TString ciao = Form("Mean = %.3f +/- %.3f", tmp_1.at(0), tmp_1.at(2));
			tex1 = new TLatex(0.125, 0.6, ciao);
			tex1->SetTextSize(0.04);
			tex1->SetNDC();
			tex1->Draw();

			ciao = Form("Sigma = %.3f +/- %.3f", tmp_1.at(1), tmp_1.at(3));
			tex1 = new TLatex(0.125, 0.55, ciao);
			tex1->SetTextSize(0.04);
			tex1->SetNDC();
			tex1->Draw();

// 			tex1 = new TLatex(0.6, 0.65, "MC:");
// 			tex1->SetTextSize(0.04);
// 			tex1->SetTextColor(kRed);
// 			tex1->SetNDC();
// 			tex1->Draw();
// 	
// 			ciao = Form("Mean = %.3f +/- %.3f", tmp_2.at(0), tmp_2.at(2));
// 			tex1 = new TLatex(0.6, 0.6, ciao);
// 			tex1->SetTextSize(0.04);
// 			tex1->SetTextColor(kRed);
// 			tex1->SetNDC();
// 			tex1->Draw();
// 
// 			ciao = Form("Sigma = %.3f +/- %.3f", tmp_2.at(1), tmp_2.at(3));
// 			tex1 = new TLatex(0.6, 0.55, ciao);
// 			tex1->SetTextSize(0.04);
// 			tex1->SetTextColor(kRed);
// 			tex1->SetNDC();
// 			tex1->Draw();
		
		}		
	}
	else{
		h1->Draw("E");
		h2->Draw("Same E");

		h1->SetTitle(nome_canvas);
		h2->SetTitle(nome_canvas);
		h1->GetXaxis()->SetTitle(x_name);
		h2->GetXaxis()->SetTitle(x_name);
		
		if(fabs(modification) == 2){
			h1->GetYaxis()->SetRangeUser(0, 0.03);//SetMaximum(max);
			h2->GetYaxis()->SetRangeUser(0, 0.03);//SetMaximum(max);	
		}
		else{
			h1->GetYaxis()->SetRangeUser(-0.0025, 0.0025);//SetMaximum(max);
			h2->GetYaxis()->SetRangeUser(-0.0025, 0.0025);//SetMaximum(max);
		}
		
	}

	canvas->Update();		

	TLegend *legend;
	if(modification < 1 || modification == 3){
		legend = new TLegend(0.65,0.75,0.9,0.9);
		legend->AddEntry(h1, "Run3");
		legend->AddEntry(h2, "UL");
	}
	if(modification == 2 || modification == 3){
		legend = new TLegend(0.65,0.15,0.9,0.35);
		legend->AddEntry(h1, "Roch");
		legend->AddEntry(h2, "Roch + VX+BS");
	}
	if(modification == -2 || modification == -3){
		legend = new TLegend(0.65,0.15,0.9,0.35);
		legend->AddEntry(h1, "UL");
		legend->AddEntry(h2, "Run3");
	}
	legend->Draw();	
		
	
	canvas->Update();
   	canvas->cd();
   	
	TPad* pad22;
	pad22 = new TPad("pad2", "pad2", 0, 0.001, 1, 0.2);
	if(modification == 0)
	   	pad22 = new TPad("pad2", "pad2", 0, 0.001, 1, 0.1);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
// 	pad22->SetTicks();
	
	if(fabs(modification) == 3)
		ratio_2->Add(h2,-1);
	else		
		ratio_2->Divide(h2);
	ratio_2->GetYaxis()->SetRangeUser(0.7, 1.2);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("Run3 / UL");
	if(modification == 2 || modification == 3)
		ratio_2->GetYaxis()->SetTitle("Roch / Roch+VX+BS");
	ratio_2->GetYaxis()->SetTitleOffset(0.50);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	ratio_2->SetLineColor(kBlack);
	if(fabs(modification) == 2) ratio_2->GetYaxis()->SetRangeUser(0.7, 1.2);
	if(fabs(modification) == 3) ratio_2->GetYaxis()->SetRangeUser(-0.001, 0.001);

	ratio_2->Draw();
	
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
// 	TLegend *legend_comparison = new TLegend(0.75,0.75,0.9,0.9);
// 	legend_comparison->AddEntry(ratio_2, "Roch+VX+BS");
// 	legend_comparison->Draw();

	canvas->Print(save);
	
}

TH1F* RecoursiveFit(TString histo_name, std::vector<double> binning, TH2F* h2, bool Mean = 0){

	TH1F* tmp = new TH1F("tmp", "tmp", binning.size()-1, &binning[0]);
	
	int bin = 1;

	for(int b = 1; b < binning.size(); b++){
		RooRealVar var("var", "var", -0.2, 0.2);
		RooRealVar MeanA("MeanA", "MeanA", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussA("Sigma_GaussA", "Sigma_GaussA", 0.01, 0.0001, 0.05);
		RooRealVar MeanB("MeanB", "MeanB", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussB("Sigma_GaussB", "Sigma_GaussB", 0.01, 0.0001, 0.05);
		RooRealVar MeanC("MeanC", "MeanC", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussC("Sigma_GaussC", "Sigma_GaussC", 0.01, 0.0001, 0.05);
		RooRealVar MeanD("MeanD", "MeanD", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussD("Sigma_GaussD", "Sigma_GaussD", 0.01, 0.0001, 0.05);
		RooRealVar MeanE("MeanE", "MeanE", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussE("Sigma_GaussE", "Sigma_GaussE", 0.01, 0.0001, 0.05);
		RooRealVar MeanF("MeanF", "MeanF", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussF("Sigma_GaussF", "Sigma_GaussF", 0.01, 0.0001, 0.05);

		RooGaussian GaussA("GaussA", "GaussA", var, MeanA, Sigma_GaussA);
		RooGaussian GaussB("GaussB", "GaussB", var, MeanB, Sigma_GaussB);
		RooGaussian GaussC("GaussC", "GaussC", var, MeanC, Sigma_GaussC);
		RooGaussian GaussD("GaussD", "GaussD", var, MeanD, Sigma_GaussD);
		RooGaussian GaussE("GaussE", "GaussE", var, MeanE, Sigma_GaussE);
		RooGaussian GaussF("GaussF", "GaussF", var, MeanF, Sigma_GaussF);


		RooRealVar CB_mean("CB_mean", "CB_mean", 0, -0.05, 0.05);
		RooRealVar Sigma("Sigma", "Sigma", 0.001, 0, 0.1);
		RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0, 30);
		RooRealVar ExpL("ExpL", "ExpL", 10, 0, 30);
		RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0, 30);
		RooRealVar ExpR("ExpR", "ExpR", 10, 1, 30);
//		RooMyPDF_DSCB Final_DY("Final_DY", "Final_DY", var, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
	        RooCBShape Final_DY("Final_DY", "Final_DY", var, CB_mean, Sigma, AlphaL, ExpL);

	
		TString nome = Form("%s_%.0f_%.0f", histo_name.Data(), binning.at(b-1), binning.at(b));
		TH1F* new_1 = (TH1F*) h2->ProjectionY(nome, b, b);
		new_1->SetTitle(nome);
		new_1->GetXaxis()->SetTitle("(p_{T}^{reco} - p_{T}^{GEN})/p_{T}^{GEN}");
		RooDataHist histo(nome, nome, var, new_1);
		TCanvas* c1 = new TCanvas(nome, nome, 700, 500);
// 		TPad* pad11 = new TPad("pad1", "pad1", 0, 0.1, 1, 1);
// 	   	pad11->SetGrid();
//    		pad11->SetBottomMargin(0.1);
// 	   	pad11->Draw();
//    		pad11->cd();
//    		pad11->SetLogy();
//    		new_1->GetYaxis()->SetRangeUser(1, new_1->GetMaximum()*1.1);
//    		new_1->Draw();
//    		c1->Update();
   		
   		RooPlot* xframe = var.frame(Title(nome));
		histo.plotOn(xframe);
	
		if(	new_1->Integral() > 0 && !Mean){
				
			GaussA.fitTo(histo, Range(-0.2, 0.2));
			GaussB.fitTo(histo, Range(MeanA.getVal()-2*Sigma_GaussA.getVal(), MeanA.getVal()+2*Sigma_GaussA.getVal()));
			GaussC.fitTo(histo, Range(MeanB.getVal()-2*Sigma_GaussB.getVal(), MeanA.getVal()+2*Sigma_GaussB.getVal()));
			GaussD.fitTo(histo, Range(MeanC.getVal()-2*Sigma_GaussC.getVal(), MeanA.getVal()+2*Sigma_GaussC.getVal()));
			GaussE.fitTo(histo, Range(MeanD.getVal()-2*Sigma_GaussD.getVal(), MeanA.getVal()+2*Sigma_GaussD.getVal()));
// 			GaussF.fitTo(histo, Range(MeanE.getVal()-2*Sigma_GaussE.getVal(), MeanA.getVal()+2*Sigma_GaussE.getVal()));
			GaussE.plotOn(xframe,RooFit::LineColor(kRed+2));
			GaussE.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
// 			GaussA.plotOn(xframe,RooFit::LineColor(kRed+2));
// 			GaussA.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kRed+2);	
			GaussD.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.85));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kGreen+2);	
			if(!Mean){					
// 				tmp->SetBinContent(bin, Sigma_GaussA.getVal());
// 				tmp->SetBinError(bin, Sigma_GaussA.getError());
// 				tmp->SetBinContent(bin, Sigma_GaussF.getVal());
// 				tmp->SetBinError(bin, Sigma_GaussF.getError());
				tmp->SetBinContent(bin, Sigma_GaussE.getVal());
				tmp->SetBinError(bin, Sigma_GaussE.getError());
			}
			
			Double_t integral = new_1->Integral();
			TString integ = Form("Integral = %.0f", integral);
			TLatex* tex1 = new TLatex(0.6,0.3, integ);
			tex1->SetNDC();
		
			xframe->Draw();

			tex1->Draw();
// 			else{
// 				tmp->SetBinContent(bin, MeanA.getVal());
// 				tmp->SetBinError(bin, MeanA.getError());
// 			}
		}
		else{
			Final_DY.fitTo(histo);
			Final_DY.plotOn(xframe,RooFit::LineColor(kRed+2));
			Final_DY.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kRed+2);	
			tmp->SetBinContent(bin, CB_mean.getVal());
			tmp->SetBinError(bin, CB_mean.getError());
			
			Double_t integral = new_1->Integral();
			TString integ = Form("Integral = %.0f", integral);
			TLatex* tex1 = new TLatex(0.6,0.3, integ);
			tex1->SetNDC();
		
			xframe->Draw();

			tex1->Draw();													
		}

		bin++;				

		if(b == 1)
			c1->Print(directory + histo_name + ".pdf[");
		c1->Print(directory + histo_name + ".pdf");
		if(b == binning.size()-1)
			c1->Print(directory + histo_name + ".pdf]");
	}
	
	return tmp;
}


TH1F* RecoursiveUnbinnedFit(TString histo_name, std::vector<double> binning, RooRealVar* rv_Res, RooDataSet* Data_Res, TString save, bool Mean = 0){

	TH1F* tmp = new TH1F("tmp", "tmp", binning.size()-1, &binning[0]);

	RooDataSet histo = RooDataSet(Data_Res->GetName(), Data_Res->GetTitle(), Data_Res, *Data_Res->get());

	int bin = 1;
	int b = 1;

// 	for(int b = 1; b < binning.size(); b++){
		RooRealVar var("var", "var", -0.2, 0.2);
		RooRealVar MeanA("MeanA", "MeanA", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussA("Sigma_GaussA", "Sigma_GaussA", 0.01, 0.0001, 0.05);
		RooRealVar MeanB("MeanB", "MeanB", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussB("Sigma_GaussB", "Sigma_GaussB", 0.01, 0.0001, 0.05);
		RooRealVar MeanC("MeanC", "MeanC", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussC("Sigma_GaussC", "Sigma_GaussC", 0.01, 0.0001, 0.05);
		RooRealVar MeanD("MeanD", "MeanD", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussD("Sigma_GaussD", "Sigma_GaussD", 0.01, 0.0001, 0.05);
		RooRealVar MeanE("MeanE", "MeanE", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussE("Sigma_GaussE", "Sigma_GaussE", 0.01, 0.0001, 0.05);
		RooRealVar MeanF("MeanF", "MeanF", 0, -0.05, 0.05);
		RooRealVar Sigma_GaussF("Sigma_GaussF", "Sigma_GaussF", 0.01, 0.0001, 0.05);

		RooGaussian GaussA("GaussA", "GaussA", *rv_Res, MeanA, Sigma_GaussA);
		RooGaussian GaussB("GaussB", "GaussB", *rv_Res, MeanB, Sigma_GaussB);
		RooGaussian GaussC("GaussC", "GaussC", *rv_Res, MeanC, Sigma_GaussC);
		RooGaussian GaussD("GaussD", "GaussD", *rv_Res, MeanD, Sigma_GaussD);
		RooGaussian GaussE("GaussE", "GaussE", *rv_Res, MeanE, Sigma_GaussE);
		RooGaussian GaussF("GaussF", "GaussF", *rv_Res, MeanF, Sigma_GaussF);

		RooRealVar CB_mean("CB_mean", "CB_mean", 0, -0.05, 0.05);
		RooRealVar Sigma("Sigma", "Sigma", 0.001, 0, 0.1);
		RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0, 30);
		RooRealVar ExpL("ExpL", "ExpL", 10, 0, 30);
		RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0, 30);
		RooRealVar ExpR("ExpR", "ExpR", 10, 1, 30);
//		RooMyPDF_DSCB Final_DY("Final_DY", "Final_DY", *rv_Res, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
                RooCBShape Final_DY("Final_DY", "Final_DY", *rv_Res, CB_mean, Sigma, AlphaL, ExpL);
	
		TString nome = Form("%s_%.0f_%.0f", histo_name.Data(), binning.at(b-1), binning.at(b));
		TCanvas* c1 = new TCanvas(nome, nome, 700, 500);
   		
   		RooPlot* xframe = rv_Res->frame(Title(nome));
		histo.plotOn(xframe);
		
		Double_t integral = histo.numEntries();
		TString integ = Form("Integral = %.0f", integral);
		TLatex* tex1 = new TLatex(0.6,0.3, integ);
		tex1->SetNDC();
	
		if(!Mean){
				
			GaussA.fitTo(histo, Range(-0.2, 0.2));
			GaussB.fitTo(histo, Range(MeanA.getVal()-2*Sigma_GaussA.getVal(), MeanA.getVal()+2*Sigma_GaussA.getVal()));
			GaussC.fitTo(histo, Range(MeanB.getVal()-2*Sigma_GaussB.getVal(), MeanA.getVal()+2*Sigma_GaussB.getVal()));
			GaussD.fitTo(histo, Range(MeanC.getVal()-2*Sigma_GaussC.getVal(), MeanA.getVal()+2*Sigma_GaussC.getVal()));
			GaussE.fitTo(histo, Range(MeanD.getVal()-2*Sigma_GaussD.getVal(), MeanA.getVal()+2*Sigma_GaussD.getVal()));
// 			GaussF.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
			GaussE.plotOn(xframe,RooFit::LineColor(kRed+2));
			GaussE.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
// 			GaussA.plotOn(xframe,RooFit::LineColor(kRed+2));
// 			GaussA.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kRed+2);	
			GaussD.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.85));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kGreen+2);	
			if(!Mean){					
				tmp->SetBinContent(bin, Sigma_GaussE.getVal());
				tmp->SetBinError(bin, Sigma_GaussE.getError());
			}
			
			xframe->Draw();

// 			tex1->Draw();										
		}
		else{
			Final_DY.fitTo(histo);
			Final_DY.plotOn(xframe,RooFit::LineColor(kRed+2));
			Final_DY.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kRed+2);	
			tmp->SetBinContent(bin, CB_mean.getVal());
			tmp->SetBinError(bin, CB_mean.getError());
					
			xframe->Draw();

// 			tex1->Draw();										

		}

		bin++;	
		
		std::cout<<"CIAO = "<<save<<std::endl;

		c1->Print(save);


// 		if(b == 1)
// 			c1->Print(directory + histo_name + "_Unbinned.pdf[");
// 		c1->Print(directory + histo_name + "_Unbinned.pdf");
// 		if(b == binning.size()-1)
// 			c1->Print(directory + histo_name + "_Unbinned.pdf]");


// 	}
	
	return tmp;
}


TH2F* RecoursiveFit3D(TString histo_name, std::vector<double> binning_x, std::vector<double> binning_y, TH3F* h2, bool Mean = 1){

	TH2F* tmp = new TH2F(histo_name, histo_name, binning_x.size()-1, &binning_x[0], binning_y.size()-1, &binning_y[0]);
		
	for(int x = 1; x < binning_x.size(); x++){
		for(int y = 1; y < binning_y.size(); y++){
			RooRealVar var("var", "var", -0.1, 0.1);
			RooRealVar MeanA("MeanA", "MeanA", 0, -0.01, 0.01);
			RooRealVar Sigma_GaussA("Sigma_GaussA", "Sigma_GaussA", 0.001, 0.0001, 0.05);

			RooGaussian GaussA("GaussA", "GaussA", var, MeanA, Sigma_GaussA);
			RooGaussian GaussB("GaussB", "GaussB", var, MeanA, Sigma_GaussA);
			RooGaussian GaussC("GaussC", "GaussC", var, MeanA, Sigma_GaussA);
			RooGaussian GaussD("GaussD", "GaussD", var, MeanA, Sigma_GaussA);
			RooGaussian GaussE("GaussE", "GaussE", var, MeanA, Sigma_GaussA);
			RooGaussian GaussF("GaussF", "GaussF", var, MeanA, Sigma_GaussA);
			TString nome = Form("%s_%.3f_%.3f_%.3f_%.3f", histo_name.Data(), binning_x.at(x-1), binning_x.at(x), binning_y.at(y-1), binning_y.at(y));
			TH1F* new_1 = (TH1F*) h2->ProjectionZ(nome, x, x, y, y);
			new_1->SetTitle(nome);
			new_1->GetXaxis()->SetTitle("(p_{T}^{VX+BS} - p_{T}^{Roch})/p_{T}^{Roch}");
			RooDataHist histo(nome, nome, var, new_1);
			TCanvas* c1 = new TCanvas(nome, nome, 700, 500);
			RooPlot* xframe = var.frame(Title(nome));
			histo.plotOn(xframe);
			GaussA.fitTo(histo, Range(-0.1, 0.1));
			GaussB.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
			GaussC.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
			GaussD.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
			GaussE.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
			GaussF.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
			GaussF.plotOn(xframe,RooFit::LineColor(kRed+2));
			GaussF.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kRed+2);	
			if(Mean){
				tmp->SetBinContent(x, y, MeanA.getVal());
				tmp->SetBinError(x, y, MeanA.getError());
			}
			else{
				tmp->SetBinContent(x, y, Sigma_GaussA.getVal());
				tmp->SetBinError(x, y, Sigma_GaussA.getError());
			}

			Double_t integral = new_1->Integral();
			TString integ = Form("Integral = %.0f", integral);
			TLatex* tex1 = new TLatex(0.6,0.3, integ);
			tex1->SetNDC();
		
			xframe->Draw();

			tex1->Draw();	
			
			if(x == 1 && y == 1)
				c1->Print(directory + histo_name + ".pdf[");
			c1->Print(directory + histo_name + ".pdf");
			if(x == binning_x.size()-1 && y == binning_y.size()-1)
				c1->Print(directory + histo_name + ".pdf]");
		}
	}
	
	return tmp;

}

TH2F* MassFit(TString histo_name, std::vector<double> binning_x, std::vector<double> binning_y, TH3F* h2, bool GEN){

	TH2F* tmp = new TH2F(histo_name, histo_name, binning_x.size()-1, &binning_x[0], binning_y.size()-1, &binning_y[0]);
		
	for(int x = 1; x < binning_x.size(); x++){
		for(int y = 1; y < binning_y.size(); y++){

// 			RooRealVar var("M_{ll} [GeV]", "M_{ll} [GeV]", 81, 101);
			RooRealVar var("M_{ll} [GeV]", "M_{ll} [GeV]", 60, 120);
			//// BW
			RooRealVar BW_mean_DY("BW_mean_DY", "BW_mean_DY", BW_mean_PDG);
			RooRealVar BW_sigma_DY("BW_sigma_DY", "BW_sigma_DY", BW_sigma_PDG);
			RooBreitWigner BW_DY("BW_DY", "BW_DY", var, BW_mean_DY, BW_sigma_DY);
			//// BW

/*	
			//// CB			
			RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);
			RooRealVar CB_sigma("CB_sigma", "CB_sigma", 1, 0., 10);
			RooRealVar CB_alpha("CB_alpha", "CB_alpha", 1, 0., 10);
			RooRealVar CB_exp("CB_exp", "CB_exp", 5, 0., 30);
		 	RooCBShape CB("CB", "CB", var, CB_mean, CB_sigma, CB_alpha, CB_exp);
			//// CB				

			RooAddPdf Final_DY("Final_DY","Final_DY", RooArgList(BWxCB, bkg), fsig);
*/
			// expo	
			RooRealVar tau("tau", "tau", 0, -5., 5);
			RooExponential bkg("bkg","bkg", var, tau);
			// expo	

			RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);//BW_mean_PDG, 81, 101);
			RooRealVar Sigma("Sigma", "Sigma", 1, 0, 30);//sigma[decay]);
			RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0, 30);//alphaL[decay]);
			RooRealVar ExpL("ExpL", "ExpL", 1, 0, 30);//expL[decay]);
			RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0, 30);//alphaR[decay]);
			RooRealVar ExpR("ExpR", "ExpR", 1, 1, 50);//expR[decay]);
// 			RooMyPDF_DSCB Final_DY("Final_DY", "Final_DY", var, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);

			RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
//			RooMyPDF_DSCB DSCB("DSCB", "DSCB", var, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
		        RooCBShape DSCB("DSCB", "DSCB", var, CB_mean, Sigma, AlphaL, ExpL);

			RooFFTConvPdf BWxDSCB("BWxDSCB","BWxDSCB", var, BW_DY, DSCB);		

// 			RooAddPdf Final_DY("Final_DY","Final_DY", RooArgList(DSCB, bkg), fsig);
			RooAddPdf Final_DY("Final_DY","Final_DY", RooArgList(BWxDSCB, bkg), fsig);

			TString nome = Form("%s_%.3f_%.3f_%.3f_%.3f", histo_name.Data(), binning_x.at(x-1), binning_x.at(x), binning_y.at(y-1), binning_y.at(y));
			TH1F* new_1 = (TH1F*) h2->ProjectionZ(nome, x, x, y, y);
			std::cout<<nome<<std::endl;
			new_1->SetTitle(nome);
			new_1->GetXaxis()->SetTitle("#Delta m / m");
			RooDataHist histo(nome, nome, var, new_1);
			TCanvas* c1 = new TCanvas(nome, nome, 700, 500);
			RooPlot* xframe = var.frame(Title(nome));
			histo.plotOn(xframe);
			TLatex* tex1 = new TLatex(0.6,0.3, "");
			if(new_1->Integral() > 0){	
				if(!GEN){			
					Final_DY.fitTo(histo);
					Final_DY.plotOn(xframe,RooFit::LineColor(kRed+2));
					Final_DY.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
					xframe->getAttText()->SetTextSize(0.025);
					xframe->getAttText()->SetTextColor(kRed+2);						
					tmp->SetBinContent(x, y, CB_mean.getVal());
					tmp->SetBinError(x, y, CB_mean.getError());
				
// 					Double_t integral = new_1->Integral();
// 					TString integ = Form("Int = %.1f", integral);
// 					tex1 = new TLatex(0.6,0.3, integ);
// 					tex1->SetNDC();
// 					std::cout<<nome<<"\t"<<integral<<std::endl;
				}
				else{
					tmp->SetBinContent(x, y, 0.001);
					tmp->SetBinError(x, y, 0);
				}

			}

			xframe->Draw();	
// 			tex1->Draw();
						
			if(x == 1 && y == 1)
				c1->Print(directory + histo_name + ".pdf[");
			c1->Print(directory + histo_name + ".pdf");
			if(x == binning_x.size()-1 && y == binning_y.size()-1)
				c1->Print(directory + histo_name + ".pdf]");
		}
	}
	
	return tmp;
}

std::vector<float> SingleMassFit(TString histo_name, TH1F* h2, TString save_name){
		
	std::cout<<"Inside the fit"<<std::endl;
	
	h2->Scale(1/h2->Integral());

	std::vector<float> tmp;

	RooRealVar var("M_{ll} [GeV]", "M_{ll} [GeV]", 60, 120);
	// BW
	RooRealVar BW_mean_DY("BW_mean_DY", "BW_mean_DY", BW_mean_PDG);
	RooRealVar BW_sigma_DY("BW_sigma_DY", "BW_sigma_DY", BW_sigma_PDG);
	RooBreitWigner BW_DY("BW_DY", "BW_DY", var, BW_mean_DY, BW_sigma_DY);
	// BW

	RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);
	RooRealVar Sigma("Sigma", "Sigma", 1, 0, 30);//sigma[decay]);
	RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0, 30);//alphaL[decay]);
	RooRealVar ExpL("ExpL", "ExpL", 1, 0, 30);//expL[decay]);
	RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0, 30);//alphaR[decay]);
	RooRealVar ExpR("ExpR", "ExpR", 1, 1, 50);//expR[decay]);
//	RooMyPDF_DSCB DSCB("DSCB", "DSCB", var, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
        RooCBShape DSCB("DSCB", "DSCB", var, CB_mean, Sigma, AlphaL, ExpL);
	
	RooFFTConvPdf BWxDSCB("BWxDSCB","BWxDSCB", var, BW_DY, DSCB);		

/*	

	//// CB			
	RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);
	RooRealVar CB_sigma("CB_sigma", "CB_sigma", 1, 0., 10);
	RooRealVar CB_alpha("CB_alpha", "CB_alpha", 1, 0., 10);
	RooRealVar CB_exp("CB_exp", "CB_exp", 5, 0., 30);
	RooCBShape CB("CB", "CB", var, CB_mean, CB_sigma, CB_alpha, CB_exp);
	//// CB				

	RooAddPdf Final_DY("Final_DY","Final_DY", RooArgList(BWxCB, bkg), fsig);
*/
	// expo	
	RooRealVar tau("tau", "tau", 0, -5., 5);
	RooExponential bkg("bkg","bkg", var, tau);
	// expo	


	RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
// 	RooMyPDF_DSCB DSCB("DSCB", "DSCB", var, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
	RooAddPdf Final_DY("Final_DY","Final_DY", RooArgList(BWxDSCB, bkg), fsig);

	RooDataHist histo(histo_name, histo_name, var, h2);

	TCanvas* c1 = new TCanvas(histo_name, histo_name, 700, 500);
	RooPlot* xframe = var.frame(Title(histo_name));
	histo.plotOn(xframe);
	if(h2->Integral() > 0){	
		Final_DY.fitTo(histo);
		Final_DY.plotOn(xframe,RooFit::LineColor(kRed+2));
		Final_DY.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
		xframe->getAttText()->SetTextSize(0.025);
		xframe->getAttText()->SetTextColor(kRed+2);						
		tmp.push_back(CB_mean.getVal());
		tmp.push_back(Sigma.getVal());
		tmp.push_back(CB_mean.getError());
		tmp.push_back(Sigma.getError());
	}
	else{
		tmp.push_back(0);
		tmp.push_back(0);
		tmp.push_back(0);
		tmp.push_back(0);
	}

	xframe->Draw();	

	c1->Print(save_name);// + ".pdf");
							
	return tmp;
}

TH2F* MassScale(TString histo_name, TH2F* h1,  TH2F* h2){

	TH2F* tmp = (TH2F*)h1->Clone();
	tmp->Reset("ICESM");
	float value = -999;
		
	for(int x_bin = 1; x_bin < h1->GetNbinsX()+1; x_bin++){
		for(int y_bin = 1; y_bin < h1->GetNbinsY()+1; y_bin++){

			if(h2->GetBinContent(x_bin,y_bin) != 0 && h1->GetBinContent(x_bin,y_bin) != 0){
// 				value = 100 * (h1->GetBinContent(x_bin,y_bin) - h2->GetBinContent(x_bin,y_bin))/h2->GetBinContent(x_bin,y_bin);
// 				value = 100 * (h1->GetBinContent(x_bin,y_bin) - h2->GetBinContent(x_bin,y_bin))/BW_mean_PDG;
				value = h1->GetBinContent(x_bin,y_bin) - h2->GetBinContent(x_bin,y_bin);				
				tmp->SetBinContent(x_bin, y_bin, value);
// 				std::cout<<x_bin<<"\t"<<y_bin<<"\t"<<h1->GetBinContent(x_bin,y_bin)<<"\t"<<h2->GetBinContent(x_bin,y_bin)<<"\t"<<value<<std::endl;
			}
			else{
// 				std::cout<<x_bin<<"\t"<<y_bin<<"\t"<<h1->GetBinContent(x_bin,y_bin)<<"\t"<<h2->GetBinContent(x_bin,y_bin)<<"\t"<<value<<std::endl;
				continue;
			}
		}
	
	}
	
	tmp->SetTitle(histo_name);
	
	return tmp;
}

void Draw_TH2F(TH2F* h1, TString nome_canvas, TString save, bool colz = 1){

	if(colz) gStyle->SetPalette(70);
// 	gStyle->SetOptStat(1111111);
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	TPad* pad11 = new TPad("pad1", "pad1", 0.01, 0.01, 1, 1);
//    	pad11->SetGrid();
   	pad11->SetBottomMargin(0.05);
  	pad11->SetRightMargin(0.2);
   	pad11->Draw();
   	pad11->cd();
	h1->SetTitle(nome_canvas);
	if(colz) h1->Draw("COLZ1 TEXT45");
	else h1->Draw("COLZ1");
	if(colz){
		h1->GetZaxis()->SetRangeUser(-1, 1);
// 		gStyle->SetPaintTextFormat("3.2f %%");
		gStyle->SetPaintTextFormat("4.3f");
	}
	canvas->Update();
   	canvas->cd();
	canvas->Print(save);
}

void Draw_TH2F_FinalScale(TH2F* h1, TString nome_canvas, TString save, bool colz = 1){

	if(colz) gStyle->SetPalette(70);
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	TPad* pad11 = new TPad("pad1", "pad1", 0.01, 0.01, 1, 1);
//    	pad11->SetGrid();
   	pad11->SetBottomMargin(0.05);
  	pad11->SetRightMargin(0.2);
   	pad11->Draw();
   	pad11->cd();
	h1->SetTitle(nome_canvas);
	if(colz) h1->Draw("COLZ1 TEXT45");
	else h1->Draw("COLZ1");
	canvas->Update();
   	canvas->cd();
	canvas->Print(save);
}

TH2F* Draw_pTRes_Jake(TH2F* A, TH2F* B, TString nome_canvas, TString save, bool LogY, bool impro){
	
	TH2F* tmp;
	if(impro){		
		tmp = (TH2F*)A->Clone();
		tmp->Add(B, -1);
		tmp->Divide(A);
	
		TCanvas* c1 = new TCanvas(nome_canvas, nome_canvas, 700, 500);
		TPad* pad11 = new TPad("pad1", "pad1", 0.01, 0.01, 1, 1);
		pad11->SetBottomMargin(0.05);
		pad11->SetRightMargin(0.2);
		pad11->Draw();
		pad11->cd();
		pad11->SetLogx();
		tmp->Draw("COLZ1 TEXT0");
		gStyle->SetPaintTextFormat("0.2f");
		tmp->GetZaxis()->SetRangeUser(0, 0.22);
		c1->Update();
		c1->cd();
		c1->Print(save);	
	}
	else{
		TCanvas* c1 = new TCanvas(nome_canvas, nome_canvas, 700, 500);
		TPad* pad11 = new TPad("pad1", "pad1", 0.01, 0.01, 1, 1);
		pad11->SetBottomMargin(0.05);
		pad11->SetRightMargin(0.2);
		pad11->Draw();
		pad11->cd();
		pad11->SetLogx();
		A->Draw("COLZ1 TEXT0");
		gStyle->SetPaintTextFormat("0.3f");
		A->GetZaxis()->SetRangeUser(0, 0.1);
		c1->Update();
		c1->cd();
		c1->Print(save);		
	}	
	
	return tmp;

}








