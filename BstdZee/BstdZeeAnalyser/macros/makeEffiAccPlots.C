{
	tsw::Samples2012 anaTuples("zEffiTree");

	//-> 0) Standard cuts applied to EffiCalcTree ...
	std::string seln_accept                 = "( (eleA_p4.Pt()>35.0 && abs(eleA_p4.Eta())<1.44 ) && (eleB_p4.Pt()>35.0 && abs(eleB_p4.Eta())<1.44) && (Zp4.M()>60.0 && Zp4.M()<120.0) )";
	std::string seln_reco                   = "bothRecod";
	std::string seln_MPhiVeto               = "(abs(ZdEta)>0.07 || abs(ZdPhi)>0.3)";
	std::string seln_acceptRecoAndPhiVeto  = tsw::AndOfCuts( tsw::AndOfCuts(seln_accept, seln_reco), seln_MPhiVeto);
	std::string seln_heepId                 = "((eleA_stdHeep & 0x08ff)==0 && (eleB_stdHeep & 0x08ff)==0)";
	std::string seln_acceptRecoPhiVetoAndId = tsw::AndOfCuts(seln_acceptRecoAndPhiVeto, seln_heepId);
	std::string seln_stdHeep = "(eleA_stdHeep==0 && eleB_stdHeep==0)";
	std::string seln_stdHeepTrk = "((eleA_stdHeep & 0x0cff)==0 && (eleB_stdHeep & 0x0cff)==0)";
	std::string seln_stdHeepCalo = "((eleA_stdHeep & 0x09ff)==0 && (eleB_stdHeep & 0x09ff)==0)";
	std::string seln_modHeepStdThr = "(eleA_modHeepStdThr==0 && eleB_modHeepStdThr==0)";
	std::string seln_modHeepTrkStdThr = "((eleA_modHeepStdThr & 0x0cff)==0 && (eleB_modHeepStdThr & 0x0cff)==0)";
	std::string seln_modHeepCaloStdThr = "((eleA_modHeepStdThr & 0x09ff)==0 && (eleB_modHeepStdThr & 0x09ff)==0)";
	std::string seln_modHeepCaloStdThr_SC = "((eleA_modHeepStdThr & 0x09ff)==0 && (eleB_modHeepStdThr & 0x09ff)==0)";
	std::string seln_modHeepCaloStdThr_5xtal = "((eleA_modHeepStdThr & 0x09ff)==0 && (eleB_modHeepStdThr & 0x09ff)==0)";
	std::string seln_modHeepCaloStdThr_4xtal = "((eleA_modHeepStdThr & 0x09ff)==0 && (eleB_modHeepStdThr & 0x09ff)==0)";
	std::string seln_modHeepCaloStdThr_3xtal = "((eleA_modHeepStdThr & 0x09ff)==0 && (eleB_modHeepStdThr & 0x09ff)==0)";
	std::string seln_modHeepColThr = "(eleA_modHeepColThr==0 && eleB_modHeepColThr==0)";
	std::string seln_genHadPtSumLeq5 = "(eleA_ptSumGenHadronsDr04<=5.0 && eleB_ptSumGenHadronsDr04<=5.0)";

	// relative effis
	tsw::EffiDefn relEffi_accept("1", seln_accept, "Accept", tsw::LightGreen);
	tsw::EffiDefn relEffi_reco(seln_accept, seln_reco, "reco", tsw::Black);
	tsw::EffiDefn relEffi_MPhiVeto(relEffi_reco.denomAndNumerCuts(), seln_MPhiVeto, "#phi road", tsw::Pink);
	tsw::EffiDefn relEffi_heepId(seln_acceptRecoAndPhiVeto, seln_heepId, "HEEP ID", tsw::Blue);
	tsw::EffiDefn relEffi_stdHeep(seln_acceptRecoPhiVetoAndId, seln_stdHeep, "Both standard isol", tsw::Red);
	tsw::EffiDefn relEffi_stdHeepTrk(seln_acceptRecoPhiVetoAndId, seln_stdHeepTrk, "Standard I_{trk}", tsw::Black);
	tsw::EffiDefn relEffi_stdHeepCalo(seln_acceptRecoPhiVetoAndId, seln_stdHeepCalo, "Standard I_{calo}", tsw::Black);
	//
	tsw::EffiDefn relEffi_modHeepStdThr(seln_acceptRecoPhiVetoAndId, seln_modHeepStdThr, "Both modified isol", tsw::Red);
	tsw::EffiDefn relEffi_modHeepTrkStdThr(seln_acceptRecoPhiVetoAndId, seln_modHeepTrkStdThr, "Modified I_{trk}", tsw::Green);
	tsw::EffiDefn relEffi_modHeepCaloStdThr(      seln_acceptRecoPhiVetoAndId, seln_modHeepCaloStdThr, "Modified I_{calo}", tsw::Pink);
	tsw::EffiDefn relEffi_modHeepCaloStdThr_SC(   seln_acceptRecoPhiVetoAndId, tsw::AndOfCuts("(eleA_scModEmH1Iso-eleA_EmH1RhoCorrn)<(2.0+0.03*eleA_p4.Et())", "(eleB_scModEmH1Iso-eleB_EmH1RhoCorrn)<(2.0+0.03*eleB_p4.Et())"), "SC-based", tsw::Black);
	tsw::EffiDefn relEffi_modHeepCaloStdThr_5xtal(seln_acceptRecoPhiVetoAndId, tsw::AndOfCuts("(eleA_inrLVetoModEmH1-eleA_EmH1RhoCorrn)<(2.0+0.03*eleA_p4.Et())", "(eleB_inrLVetoModEmH1-eleB_EmH1RhoCorrn)<(2.0+0.03*eleB_p4.Et())"), "5 xtal veto", tsw::Green);
	tsw::EffiDefn relEffi_modHeepCaloStdThr_4xtal(seln_acceptRecoPhiVetoAndId, tsw::AndOfCuts("(eleA_inrMVetoModEmH1-eleA_EmH1RhoCorrn)<(2.0+0.03*eleA_p4.Et())", "(eleB_inrMVetoModEmH1-eleB_EmH1RhoCorrn)<(2.0+0.03*eleB_p4.Et())"), "4 xtal veto", tsw::Red);
	tsw::EffiDefn relEffi_modHeepCaloStdThr_3xtal(seln_acceptRecoPhiVetoAndId, tsw::AndOfCuts("(eleA_inrXSVetoModEmH1-eleA_EmH1RhoCorrn)<(2.0+0.03*eleA_p4.Et())", "(eleB_inrXSVetoModEmH1-eleB_EmH1RhoCorrn)<(2.0+0.03*eleB_p4.Et())"), "3 xtal veto", tsw::Blue);
	//
	tsw::EffiDefn relEffi_modHeepColThr(seln_acceptRecoPhiVetoAndId, seln_modHeepColThr, "ModIso (col thr)", tsw::LightBlue);
	tsw::EffiDefn relEffi_genHadPtSumLeq5(seln_acceptRecoPhiVetoAndId, seln_genHadPtSumLeq5, "genHad ptSum < 5", tsw::Blue);

	tsw::EffiDefn absEffi_accept("1", seln_accept, "Accept", tsw::LightGreen);
	tsw::EffiDefn absEffi_reco("1", relEffi_reco.denomAndNumerCuts(), "reco", tsw::Black);
	tsw::EffiDefn absEffi_MPhiVeto("1", relEffi_MPhiVeto.denomAndNumerCuts(), "#phi road", tsw::Pink);
	tsw::EffiDefn absEffi_heepId("1", relEffi_heepId.denomAndNumerCuts(), "heepId", tsw::Blue);
	tsw::EffiDefn absEffi_stdHeep("1", relEffi_stdHeep.denomAndNumerCuts(), "StdIso", tsw::Red);
	tsw::EffiDefn absEffi_modHeepStdThr("1", relEffi_modHeepStdThr.denomAndNumerCuts(), "ModIso (std thr)", tsw::Green);
	tsw::EffiDefn absEffi_modHeepColThr("1", relEffi_modHeepColThr.denomAndNumerCuts(), "ModIso (col thr)", tsw::LightBlue);

	std::string seln_modCaloIso_normThr = tsw::AndOfCuts("((eleA_inrMVetoModEmH1-eleA_EmH1RhoCorrn)<(2.0+0.03*eleA_p4.Et()))",
                                                            "((eleB_inrMVetoModEmH1-eleB_EmH1RhoCorrn)<(2.0+0.03*eleB_p4.Et()))");
	std::string seln_modCaloIso_loose   = tsw::AndOfCuts("((eleA_inrMVetoModEmH1-eleA_EmH1RhoCorrn)<(10.0+0.04*eleA_p4.Et()))",
			                                                   "((eleB_inrMVetoModEmH1-eleB_EmH1RhoCorrn)<(10.0+0.04*eleB_p4.Et()))");
	std::string seln_modTrkIso_normThr  = tsw::AndOfCuts("( eleA_inrXSVetoModTrk<5.0 )",
	                                                         "( eleB_inrXSVetoModTrk<5.0 )");
	std::string seln_modTrkIso_loose    = tsw::AndOfCuts("( eleA_inrXSVetoModTrk<(8.0+0.06*eleA_p4.Et()) )",
	                                                         "( eleB_inrXSVetoModTrk<(8.0+0.06*eleB_p4.Et()) )");


	//------------------------------------//
	// ---  'Acceptance'-style plots  --- //
	//------------------------------------//

//   tsw::CompositeMC dyJetsMC4Acc("DY[ee]+Jets, MG", tsw::Blue, "Diff cuts", "");
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-MG_*_zEffiTree.root",      "MADGRAPH", tsw::Blue, /*TWiki NNLO*/(3503.7e3)/3.0, "ZpT<160" ) , 0.0);
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-Pt100-MG_zEffiTree.root", "pT>160", tsw::Blue, /*TWiki NNLO*/(3503.7e3)*(34.1/2950.0)/3.0,  "ZpT>160" ) , 1.0);
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-Pt100-MG_zEffiTree.root", "pT>600", tsw::Blue, /*TWiki NNLO*/(3503.7e3)*(34.1/2950.0)/3.0,  "ZpT>600" ) , 2.0);
//   dyJetsMC4Acc.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120723/DYJetsToEE-Pt100-MG_zEffiTree.root", "pT>600 && dR<0.3", tsw::Blue, /*TWiki NNLO*/(3503.7e3)*(34.1/2950.0)/3.0,  "ZpT>600 && ZdR<0.3" ) , 3.0);
//
//	tsw::AcceptancePlotter accPlotter("zBosonEffiTree", "1", 0.46);
//	accPlotter.add( qStar_GI ).add( qStar_GI_42X );//.add( qStar_CI )
//			//.add( dyEe_calchep ).add( dyEe_mg_merged );
////	accPlotter.add( absEffi_accept )
////			.add( absEffi_reco ).add( absEffi_MPhiVeto )
////			.add( absEffi_heepId )
////			.add( absEffi_stdHeep )
////			.add( absEffi_modHeepColThr )
////			.add( absEffi_modHeepStdThr )
//	accPlotter//.add( relEffi_modHeepColThr ).add( relEffi_modHeepStdThr )
//			.add( relEffi_genHadPtSumLeq5 )
//			.outFilePrefix("results/20121203/signalMc_isolEffiVsMass").run();


	//----------------------------//
	// --- MOD ISO EFFI PLOTS --- //
	//----------------------------//

	tsw::AxisDefn zPt_calchep("Zp4.Pt()", 30, 280.0, 1480.0, "p_{T,ee} [GeV]");
	tsw::AxisDefn zE_calchep("Zp4.E()", 28, 280.0, 2520.0, "Z boson energy, E_{ee} [GeV]");
	tsw::AxisDefn zMom_calchep("Zp4.P()", 30, 280.0, 1480.0, "p_{ee} [GeV]");
	tsw::AxisDefn dRee_calchep("ZdR", 15, 0.1, 0.7, "#Delta R_{ee}");
	double zPtLims_mg[] = {280., 320., 360., 400., 440., 480., 520., 600., 700., 1000.};
	double dRLims_mg[] = {.1, .3, .34, .38, .42, .44, .46, .48, .5, .52, .54, .56, .58, .6, .62, .64, .66, .68, .7};
	tsw::AxisDefn zPt_mg("Zp4.E()", 9, zPtLims_mg, "Z boson energy, E_{ee} [GeV]");
	tsw::AxisDefn zMom_mg("Zp4.P()", 9, zPtLims_mg, "p_{ee} [GeV]");
	tsw::AxisDefn dRee_mg("ZdR", 18, dRLims_mg, "#Delta R_{ee}");


	// Effis
	tsw::EffiPlotter effiPlotter_recoAndId_b4PhiVeto("zBosonEffiTree", 0.0, true);
	effiPlotter_recoAndId_b4PhiVeto//.add( dyEe_calchep ).add( zPt_calchep ).add( dRee_calchep )
		.add( tsw::EffiDefn(seln_accept, seln_reco, "RECO", tsw::Black) )
		.add( tsw::EffiDefn(seln_accept, tsw::AndOfCuts(seln_heepId, seln_reco), "RECO + HEEP ID", tsw::Blue) )
		.gridLinesY().descriptiveText("Within acceptance").outFilePrefix("results/20121203/effi/effiAfterAcc_RecoId");
	effiPlotter_recoAndId_b4PhiVeto
	.drawPlot( anaTuples.dyEE_calchep(), zPt_calchep ).drawPlot( anaTuples.dyEE_calchep(), zE_calchep ).drawPlot( anaTuples.dyEE_calchep(), dRee_calchep );

	tsw::EffiPlotter effiPlotter_recoAndId_b4PhiVetoSmallDr("zBosonEffiTree", 0.0);
	effiPlotter_recoAndId_b4PhiVetoSmallDr
			.add( tsw::EffiDefn(tsw::AndOfCuts(seln_accept, "ZdR<0.3"), seln_reco, "RECO", tsw::Black) )
			.add( tsw::EffiDefn(tsw::AndOfCuts(seln_accept, "ZdR<0.3"), tsw::AndOfCuts(seln_heepId, seln_reco), "RECO + HEEP ID", tsw::Blue) )
//			.add( tsw::EffiDefn(tsw::AndOfCuts(tsw::AndOfCuts(seln_accept, "ZdR<0.3"),seln_reco), seln_heepId, "HEEP ID wrt RECO", tsw::Pink) )
			.gridLinesY().descriptiveText("Within acceptance, #Delta R < 0.3").outFilePrefix("results/20121203/effi/effiAfterAccDrLeThan03_RecoId");
	effiPlotter_recoAndId_b4PhiVetoSmallDr
			.drawPlot( anaTuples.dyEE_calchep(), zPt_calchep ).drawPlot( anaTuples.dyEE_calchep(), zE_calchep ).drawPlot( anaTuples.dyEE_calchep(), dRee_calchep )
			.drawPlot( anaTuples.dyEE_calchep(), tsw::AxisDefn("abs(ZdEta)", 14, 0.0, 0.35, "#Delta#eta_{ee}") );

	tsw::EffiPlotter effiPlotter_recoAndPhiRdVeto("zBosonEffiTree", 0.0);
	effiPlotter_recoAndPhiRdVeto
			.add( tsw::EffiDefn(seln_accept, seln_reco, "RECO", tsw::Black) )
			.add( tsw::EffiDefn(seln_accept, seln_acceptRecoAndPhiVeto, "RECO + #phi rd veto", tsw::Pink) )
			.gridLinesY().descriptiveText("Within acceptance").outFilePrefix("results/20121203/effi/effiAfterAcc_RecoPhiRoads");
	effiPlotter_recoAndPhiRdVeto.drawPlot( anaTuples.dyEE_calchep(), zPt_calchep ).drawPlot( anaTuples.dyEE_calchep(), zE_calchep )
			.drawPlot( anaTuples.dyEE_calchep(), dRee_calchep );

	tsw::EffiPlotter effiPlotter_IdAndStdIsosRelToRecoWPhiVeto("zBosonEffiTree", 0.0);
	effiPlotter_IdAndStdIsosRelToRecoWPhiVeto.gridLinesY().legendPosition(tsw::LOWER_LEFT)
			.add( relEffi_heepId ).add( tsw::EffiDefn(seln_acceptRecoAndPhiVeto, tsw::AndOfCuts(seln_heepId,seln_stdHeepTrk), "ID + track iso", tsw::Red) )
			.add( tsw::EffiDefn(seln_acceptRecoAndPhiVeto, tsw::AndOfCuts(seln_heepId,seln_stdHeepCalo), "ID + cal iso", tsw::Green) )
			.descriptiveText("Recontructed ee pairs after #phi road veto").outFilePrefix("results/20121203/effi/effiAfterAccRecoPhiRd_IdStdIso");
	effiPlotter_IdAndStdIsosRelToRecoWPhiVeto.drawPlot( anaTuples.dyEE_calchep(), zPt_calchep ).drawPlot( anaTuples.dyEE_calchep(), zE_calchep )
			.drawPlot( anaTuples.dyEE_calchep(), dRee_calchep );


	tsw::EffiPlotter effiPlotter_modTrkIso("zBosonEffiTree", 0.7);
	tsw::EffiPlotter effiPlotter_modCalIso(effiPlotter_modTrkIso);
	effiPlotter_modTrkIso.add( relEffi_stdHeepTrk ).add( relEffi_modHeepTrkStdThr )
			.gridLinesY().outFilePrefix("results/20121203/effi/modTrkIsoEffiAfterId");
	effiPlotter_modTrkIso
			.drawPlot( anaTuples.dyEE_calchep(), zPt_calchep ).drawPlot( anaTuples.dyEE_calchep(), zE_calchep ).drawPlot( anaTuples.dyEE_calchep(), dRee_calchep )
			;/*.drawPlot( anaTuples.dyEE_mg_merged(), zPt_mg )   .drawPlot( anaTuples.dyEE_mg_merged(), dRee_mg );*/

	effiPlotter_modCalIso.add( relEffi_modHeepCaloStdThr_SC ).add( relEffi_modHeepCaloStdThr_5xtal )
			.add( relEffi_modHeepCaloStdThr_4xtal ).add( relEffi_modHeepCaloStdThr_3xtal )
			.gridLinesY().legendPosition(tsw::LOWER_LEFT).outFilePrefix("results/20121203/effi/modCalIsoEffiAfterId");
	effiPlotter_modCalIso
			.drawPlot( anaTuples.dyEE_calchep(), zPt_calchep ).drawPlot( anaTuples.dyEE_calchep(), zE_calchep ).drawPlot( anaTuples.dyEE_calchep(), dRee_calchep )
			;/*.drawPlot( anaTuples.dyEE_mg_merged(), zPt_mg )   .drawPlot( anaTuples.dyEE_mg_merged(), dRee_mg );*/

//	tsw::EffiPlotter effiPlotter_chep("zBosonEffiTree", 0.7);
//	tsw::EffiPlotter effiPlotter_mg(effiPlotter_chep);
//	effiPlotter_chep.add( anaTuples.dyEE_calchep() ).add( zPt_calchep ).add( dRee_calchep );
//	effiPlotter_mg.add( anaTuples.dyEE_mg_pt100() ).add( zPt_calchep ).add( dRee_calchep );
//
//	tsw::EffiPlotter effiPlotter_chep2( effiPlotter_chep );
//	tsw::EffiPlotter effiPlotter_mg2( effiPlotter_mg );
//
//	effiPlotter_chep.add( relEffi_stdHeepTrk  ).add( relEffi_modHeepTrkStdThr )
//			.outFilePrefix("results/20121203/accEffi/modTrkIsoEffiDyCalcHEP_").run();
//	effiPlotter_mg.add( relEffi_stdHeepTrk  ).add( relEffi_modHeepTrkStdThr )
//			.outFilePrefix("results/20121203/accEffi/modTrkIsoEffiDyMG_").run();
//	effiPlotter_chep2
//			.add( relEffi_modHeepCaloStdThr_SC ).add( relEffi_modHeepCaloStdThr_5xtal )
//			.add( relEffi_modHeepCaloStdThr_4xtal ).add( relEffi_modHeepCaloStdThr_3xtal )
//			.outFilePrefix("results/20121203/accEffi/modCaloIsoEffiDyCalcHEP_").run();
//	effiPlotter_mg2
//			.add( relEffi_modHeepCaloStdThr_SC ).add( relEffi_modHeepCaloStdThr_5xtal )
//			.add( relEffi_modHeepCaloStdThr_4xtal ).add( relEffi_modHeepCaloStdThr_3xtal )
//			.outFilePrefix("results/20121203/accEffi/modCaloIsoEffiDyMG_").run();

}
