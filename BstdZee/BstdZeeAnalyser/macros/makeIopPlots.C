{

   tsw::Samples2012 anaTuples("modIsoZCandTree");

   tsw::DistPlotter comp1;
   comp1.setTree("modIsoZBosonTree");
   comp1.setSelection("(Zmass<105 && Zmass>75) && (abs(dPhi)>=0.3 || abs(dEta)>=0.07)");
   comp1.setWeight("genWeight");
//	mcComparison.rescaleMC();
//	mcComparison.descriptiveText("#sqrt{s} = 8 TeV, #int L dt = 8.4 fb^{-1}, Runs 12A,B&C; HEEPModIso, #phi-road veto, 75-105, PU reweighted");

   comp1.add( tsw::AxisDefn("ZpT", 43, 0.0, 860.0, "p_{T,ee} [GeV]") );
   comp1.data( anaTuples.data2012() );
   comp1.add_signal( anaTuples.qStarGI_M1000() );
   comp1.add_signal( anaTuples.qStarGI_M1500() );

   tsw::DistPlotter comp2( comp1 );

   comp2.add( anaTuples.wJets() );
   comp2.add( anaTuples.dyTauTau_mg() );
   comp2.add( anaTuples.topBkgds() );
   comp2.add( anaTuples.vzBkgds() );
   comp1.add( anaTuples.dyEE_mg_merged() );
   comp2.add( anaTuples.dyEE_mg_merged() );

   comp1.outFilePrefix("results/20130401/IopMcOnly/ZCandMCOnlyDY_20fb8TeV_ModIso_phiRdVeto");
   comp2.outFilePrefix("results/20130401/IopMcOnly/ZCandMCAll_20fb8TeV_ModIso_phiRdVeto");

   comp1.run();
   comp2.run();

	
   /*	mcComparison.add( tsw::AxisDefn("Zp4.P()", 70, 0.0, 1400.0, "Z momentum [GeV]") );
	mcComparison.add( tsw::AxisDefn("Zmass", 55, 65.0, 120.0, "M_{ee} [GeV]") );
	mcComparison.add( tsw::AxisDefn("Zp4.Eta()", 40, -4.0, +4.0, "Z boson #eta") );
	mcComparison.add( tsw::AxisDefn("Zp4.Phi()", 20, -3.14, +3.14, "Z boson #phi") );
	mcComparison.add( tsw::AxisDefn("eleA_p4.Pt()", 45, 0.0, 900.0, "Leading ele p_{T} [GeV]") );
	mcComparison.add( tsw::AxisDefn("eleB_p4.Pt()", 45, 0.0, 450.0, "Sub-leading ele p_{T} [GeV]") );
//	mcComparison.add( tsw::AxisDefn("nVtx", 40, -0.5, 39.5, "N_{vtx}"));

	mcComparison.add( tsw::AxisDefn("dR", 30, 0.0, 6.0, "#DeltaR_{ee}"));
	mcComparison.add( tsw::AxisDefn("dEta", 32, -3.2, +3.2, "#Delta#eta_{ee}"));
	mcComparison.add( tsw::AxisDefn("dPhi", 20, -3.14, +3.14, "#Delta#phi_{ee}"));

	mcComparison.run( );

	tsw::DistPlotter sigMcComp(true);
	sigMcComp.rescaleMC();
	sigMcComp.setTree("modIsoZBosonTree");
	sigMcComp.setSelection("(Zmass<105 && Zmass>75) && (abs(dPhi)>=0.3 || abs(dEta)>=0.07)");
	sigMcComp.setWeight("genWeight");
	sigMcComp.add(qStar_GI_M500);
	sigMcComp.add(qStar_GI_M1500);
	sigMcComp.add(qStar_GI_M2000);
	sigMcComp.add(qStar_GI_M2500);

	sigMcComp.add( tsw::AxisDefn("ZpT", 50, 0.0, 1000.0, "p_{T,ee} [GeV]") );
	sigMcComp.add( tsw::AxisDefn("Zp4.P()", 50, 0.0, 2000.0, "Z momentum [GeV]") );*/
//	sigMcComp.run();
}
