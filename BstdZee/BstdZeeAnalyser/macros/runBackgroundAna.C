{
	tsw::Samples2012 modIsoAnaTuples("modIsoZCandTree");
	tsw::Samples2012 emuAnaTuples("eleMuTree");
	tsw::Samples2012 abcdAnaTuples("abcdDiGsfTree");

	// EMU CONTROL REGION //

	tsw::EMuControlPlotMaker emuPlotter("eleMuTree");
	emuPlotter.selection("(eleMu_p4.M()<105 && eleMu_p4.M()>75)");
	emuPlotter.mcWeight("genWeight*puWeight");
	emuPlotter.chargeBranches("ele_charge", "muon_charge");
//	emuPlotter.outFilePrefix("results/20121121/DataVsMc_emu/emu_heepV41_interimMuId_approx13fb");
//	emuPlotter.descriptiveText("#sqrt{s} = 8 TeV, #int L dt = 12 fb^{-1} approx, Runs 12A-D; HEEPModIso & MuHighPt, 75-105, PU reweighted");

//	tsw::MCSample wJets_emu( emuAnaTuples.wJets() );
//	wJets_emu.mBaseSeln = tsw::AndOfCuts( wJets_emu.mBaseSeln, "eleMu_p4.Pt()<220.");
//	emuPlotter.add( wJets_emu );
//	emuPlotter.add( emuAnaTuples.dyTauTau_powheg() );
	emuPlotter.add( emuAnaTuples.vzBkgds()  );
	emuPlotter.add( emuAnaTuples.topBkgds() );

	emuPlotter.data( emuAnaTuples.muEG2012() );
	emuPlotter.includeRatioPlots(true);
	emuPlotter.drawPlots( tsw::AxisDefn("eleMu_p4.Pt()", 33, 40.0, 700.0, "p_{T,e#mu} [GeV]") );
	emuPlotter.includeRatioPlots(false);
	emuPlotter.drawPlots( tsw::AxisDefn("eleMu_p4.M()", 25, 71.0, 121.0, "Mass, M_{e#mu} [GeV]") );
//	emuPlotter.add( tsw::AxisDefn("ele_p4.Pt()", 24, 20.0, 500.0, "Electron p_{T} [GeV]") );
//	emuPlotter.add( tsw::AxisDefn("muon_p4.Pt()", 24, 20.0, 500.0, "Muon p_{T} [GeV]") );
//	emuPlotter.run();


	// QCD -- SIDEBAND-BASED //

	/*tsw::SidebandJetsEstimator sidebandEstimator("modIsoZBosonTree", "1", "1");
	sidebandEstimator.ptBranch("ZpT");
	sidebandEstimator.massBranch("Zmass");
//	sidebandEstimator.smearedMassBranch4mc("(tsw::GausianMultSmear(Zmass, 0.9935, 0.0082))");
	sidebandEstimator.setBaselineSelection( "(abs(dPhi)>=0.3 || abs(dEta)>=0.07)" );

	sidebandEstimator.add(  modIsoAnaTuples.dyTauTau_powheg() );
	sidebandEstimator.add(  modIsoAnaTuples.topBkgds() );
	sidebandEstimator.add(  modIsoAnaTuples.vzBkgds() );
	sidebandEstimator.add(  modIsoAnaTuples.dyEE_mg_merged() );
	sidebandEstimator.data( modIsoAnaTuples.data2012() );

	sidebandEstimator.outputFileDir( "" );
	sidebandEstimator.outputFileTag( "results/20121121/Jets_sideband/ZCandSideband_8TeV_12fb" );
	sidebandEstimator.run();*/


	// QCD -- ABCD METHOD //

	/*tsw::AbcdJetsEstimator abcdEstimator("abcdDiGsfFrPreTree");
	abcdEstimator
			.mcWeights("genWeight*puWeight")
			.baselineSeln("(Zmass<120.0 && Zmass>60.0) && trgDecision && (eleA_frPre==0 && eleB_frPre==0)")
			.data( abcdAnaTuples.data2012() )
			.add(  abcdAnaTuples.topBkgds() )
			.add(  abcdAnaTuples.vzBkgds() )
			.add(  abcdAnaTuples.dyTauTau_powheg() )
			.add(  abcdAnaTuples.dyEE_mg_merged() )
			.add(  abcdAnaTuples.wJets()  )
			.outFilePrefix("results/20121121/Jets_ABCD/ABCD_12fb_");
	abcdEstimator.run();*/

	// QCD -- COMPARISON //

}
