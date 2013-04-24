{
   tsw::Samples2012 modIsoAnaTuples("modIsoZCandTree");
   tsw::Samples2012 emuAnaTuples("eleMuTree");
   tsw::Samples2012 abcdAnaTuples("abcdDiGsfTree");

   // EMU CONTROL REGION //

   tsw::EMuControlPlotMaker emuPlotter("eleMuTree");
   emuPlotter.selection("(eleMu_p4.M()<105 && eleMu_p4.M()>75) && (muon_p4.Pt()>35.0 && ele_p4.Pt()>35.0) && abs(muon_p4.Eta())<1.442 && abs(ele_p4.Eta())<1.442");
   emuPlotter.mcWeight("genWeight * puWeight * 0.995 * tsw::effi_hltMu22( muon_p4.Eta() ) * 0.995");
   emuPlotter.chargeBranches("ele_charge", "muon_charge");
   emuPlotter.outFilePrefix("results/dummy/DataVsMc_emu/emu_heepV41_hiPtMuId_8TeV_all2012");
// emuPlotter.descriptiveText("#sqrt{s} = 8 TeV, #int L dt = 12 fb^{-1} approx, Runs 12A-D; HEEPModIso & MuHighPt, 75-105, PU reweighted");

   tsw::MCSample wJets_emu( emuAnaTuples.wJets() );
   wJets_emu.mBaseSeln = tsw::AndOfCuts( wJets_emu.mBaseSeln, "eleMu_p4.Pt()<220.");
   emuPlotter.add( wJets_emu );
   emuPlotter.add( emuAnaTuples.dyEE_mg_merged() );
   emuPlotter.add( emuAnaTuples.dyMuMu_mg_merged() );
   emuPlotter.add( emuAnaTuples.dyTauTau_mg() );
     //emuPlotter.add( emuAnaTuples.dyLL_mg() );
   emuPlotter.add( emuAnaTuples.vzBkgds()  );
   emuPlotter.add( emuAnaTuples.topBkgds() );

   emuPlotter.data( emuAnaTuples.muEG2012() );
   emuPlotter.includeRatioPlots(true);
   emuPlotter.drawPlots( tsw::AxisDefn("eleMu_p4.Pt()", "0,10,20,30,40,50,60,70,80,90, 100,110,120,130,140,150,160,170,180,190, 200,210,220,240,260,320,700", "p_{T,e#mu} [GeV]", 10.0) );
   //emuPlotter.includeRatioPlots(false);
   //emuPlotter.drawPlots( tsw::AxisDefn("eleMu_p4.M()", 25, 71.0, 121.0, "Mass, M_{e#mu} [GeV]") );
   emuPlotter.add( tsw::AxisDefn("ele_p4.Pt()", 48, 20.0, 500.0, "Electron p_{T} [GeV]") );
   emuPlotter.add( tsw::AxisDefn("muon_p4.Pt()", 48, 20.0, 500.0, "Muon p_{T} [GeV]") );
   emuPlotter.run();

   // QCD -- SIDEBAND-BASED //

   /*tsw::SidebandJetsEstimator sidebandEstimator("modIsoZBosonTree", "1", "1");
   sidebandEstimator.ptBranch("ZpT");
   sidebandEstimator.massBranch("Zmass");
// sidebandEstimator.smearedMassBranch4mc("(tsw::GausianMultSmear(Zmass, 0.9935, 0.0082))");
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
      .add(  abcdAnaTuples.dyTauTau_pythia() )
      .add(  abcdAnaTuples.dyEE_mg_merged() )
      .add(  abcdAnaTuples.wJets()  )
      .outFilePrefix("results/20130313/qcd/ABCD_full2012_");
      abcdEstimator.run();*/

   // QCD -- COMPARISON //

}
