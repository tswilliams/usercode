{
	tsw::Samples2012 modIsoAnaTuples("modIsoZCandTree");
	tsw::Samples2012 noIsoAnaTuples("noIsoZCandTree");
	tsw::Samples2012 emuAnaTuples("eleMuTree");

	//-> 1) Data vs MC plots
	tsw::DistPlotter mcComparison;
	mcComparison.setTree("modIsoZBosonTree");
	mcComparison.setSelection("(Zmass<105 && Zmass>75) && (abs(dPhi)>=0.3 || abs(dEta)>=0.07) && (eleA_modHeepStdThr==0 && eleB_modHeepStdThr==0)");
	mcComparison.setWeight("genWeight*puWeight");
	mcComparison.rescaleMC();
//	mcComparison.descriptiveText("#sqrt{s} = 8 TeV, #int L dt = 12 fb^{-1} approx, Runs 12A-D; HEEPModIsoStdThr, #phi-road veto, 75-105, PU reweighted");
//	mcComparison.outFilePrefix("results/20121121/DataVsMc_ee/ZCand_ModIso75to105_phiRdVeto_PUweight_8TeV_12fb");

	mcComparison.add( modIsoAnaTuples.wJets() );
	mcComparison.add( modIsoAnaTuples.dyTauTau_powheg() );
	mcComparison.add( modIsoAnaTuples.topBkgds() );
	mcComparison.add( modIsoAnaTuples.vzBkgds()  );
	mcComparison.add( modIsoAnaTuples.dyEE_mg_merged() );

	mcComparison.add_signal( modIsoAnaTuples.qStarGI_M500() );
	mcComparison.add_signal( modIsoAnaTuples.qStarGI_M1500() );

	mcComparison.data( modIsoAnaTuples.data2012() );

	mcComparison.add( tsw::AxisDefn("ZpT", 50, 0.0, 1000.0, "p_{T,ee} [GeV]"), true, true );
//	mcComparison.add( tsw::AxisDefn("Zp4.P()", 70, 0.0, 1400.0, "Z momentum [GeV]"), true, true);
//	mcComparison.add( tsw::AxisDefn("Zmass", 55, 65.0, 120.0, "M_{ee} [GeV]"));
//	mcComparison.add( tsw::AxisDefn("Zp4.Eta()", 40, -4.0, +4.0, "Z boson #eta"));
//	mcComparison.add( tsw::AxisDefn("Zp4.Phi()", 20, -3.14, +3.14, "Z boson #phi"));
//	mcComparison.add( tsw::AxisDefn("eleA_p4.Pt()", 45, 0.0, 900.0, "Leading ele p_{T} [GeV]"));
//	mcComparison.add( tsw::AxisDefn("eleB_p4.Pt()", 45, 0.0, 450.0, "Sub-leading ele p_{T} [GeV]"), true, true);
//	mcComparison.add( tsw::AxisDefn("nVtx", 40, -0.5, 39.5, "N_{vtx}"), false, true);
//
//	mcComparison.add( tsw::AxisDefn("dR", 30, 0.0, 6.0, "#DeltaR_{ee}"));
//	mcComparison.add( tsw::AxisDefn("dEta", 32, -3.2, +3.2, "#Delta#eta_{ee}"));
//	mcComparison.add( tsw::AxisDefn("dPhi", 20, -3.14, +3.14, "#Delta#phi_{ee}"));

	mcComparison.run( );



//	//-> 2) Acceptance plots
//   tsw::CompositeMC dyJetsMC4Acc("DY[ee]+Jets, MG", tsw::Blue, "Diff cuts", "");
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-MG_*_zEffiTree.root",      "MADGRAPH", tsw::Blue, /*TWiki NNLO*/(3503.7e3)/3.0, "ZpT<160" ) , 0.0);
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-Pt100-MG_zEffiTree.root", "pT>160", tsw::Blue, /*TWiki NNLO*/(3503.7e3)*(34.1/2950.0)/3.0,  "ZpT>160" ) , 1.0);
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-Pt100-MG_zEffiTree.root", "pT>600", tsw::Blue, /*TWiki NNLO*/(3503.7e3)*(34.1/2950.0)/3.0,  "ZpT>600" ) , 2.0);
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-Pt100-MG_zEffiTree.root", "pT>600 && dR<0.3", tsw::Blue, /*TWiki NNLO*/(3503.7e3)*(34.1/2950.0)/3.0,  "ZpT>600 && ZdR<0.3" ) , 3.0);
//
//	tsw::AcceptancePlotter accPlotter("zBosonEffiTree", "1");
//	std::string effiSeln_acceptWithMPhiVeto = "( mcAccept_ebeb && ( abs(ZdEta)>0.07 || abs(ZdPhi)>0.3) )";
//	std::string effiSeln_reco               = "bothRecod";
//	std::string effiSeln_recoAndAccept      = tsw::AndOfCuts(effiSeln_acceptWithMPhiVeto, effiSeln_reco);
//	std::string effiSeln_heepId             = "cut_both_heepId";
//	std::string effiSeln_modCaloIso_normThr = tsw::AndOfCuts("((eleA_inrMVetoModEmH1-eleA_EmH1RhoCorrn)<(2.0+0.03*eleA_p4.Et()))",
//                                                            "((eleB_inrMVetoModEmH1-eleB_EmH1RhoCorrn)<(2.0+0.03*eleB_p4.Et()))");
//	std::string effiSeln_modCaloIso_loose   = tsw::AndOfCuts("((eleA_inrMVetoModEmH1-eleA_EmH1RhoCorrn)<(10.0+0.04*eleA_p4.Et()))",
//			                                                   "((eleB_inrMVetoModEmH1-eleB_EmH1RhoCorrn)<(10.0+0.04*eleB_p4.Et()))");
//	std::string effiSeln_modTrkIso_normThr  = tsw::AndOfCuts("( eleA_inrXSVetoModTrk<5.0 )",
//	                                                         "( eleB_inrXSVetoModTrk<5.0 )");
//	std::string effiSeln_modTrkIso_loose    = tsw::AndOfCuts("( eleA_inrXSVetoModTrk<(8.0+0.06*eleA_p4.Et()) )",
//	                                                         "( eleB_inrXSVetoModTrk<(8.0+0.06*eleB_p4.Et()) )");
//	accPlotter.add( dyJetsMC4Acc );
//	accPlotter.add( tsw::EffiDefn(effiSeln_acceptWithMPhiVeto, effiSeln_reco, "reco", tsw::Black) );
//	accPlotter.add( tsw::EffiDefn(effiSeln_recoAndAccept, effiSeln_heepId, "IdNoIso", tsw::Blue) );
//	accPlotter.add( tsw::EffiDefn(tsw::AndOfCuts(effiSeln_recoAndAccept, effiSeln_heepId), tsw::AndOfCuts(effiSeln_modCaloIso_normThr, effiSeln_modTrkIso_normThr), "ModIso, norm", tsw::Red) );
//	accPlotter.add( tsw::EffiDefn(tsw::AndOfCuts(effiSeln_recoAndAccept, effiSeln_heepId), tsw::AndOfCuts(effiSeln_modCaloIso_loose  , effiSeln_modTrkIso_loose  ), "ModIso, loose", tsw::Green) );
//	accPlotter.run();




}
