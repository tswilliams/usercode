{
	tsw::Samples2012 anaTuples("noIsoZCandTree");

	std::string out_file_base = "results/20130213/dataVsMc_otherEleInrVetoIso/isoOtherEleInrVeto_ZpeakLeadEle";

	tsw::DistPlotter p;
	p.setTree("noIsoZBosonTree");
	p.setWeight("genWeight*puWeight");
	p.rescaleMC();
	p.add( anaTuples.wJets() );
	p.add( anaTuples.dyTauTau_powheg() );
	p.add( anaTuples.topBkgds() );
	p.add( anaTuples.vzBkgds()  );
	p.add( anaTuples.dyEE_mg_merged() );
	p.data( anaTuples.data2012() );

	p.add( tsw::AxisDefn("eleA_modIsoTk_otherEleVetoForSelf", 70, 0.0, 14.0, "Track isolation, I_{trk} [GeV]") );
	p.add( tsw::AxisDefn("eleA_modIsoEcal_otherEleVetoForSelf+eleA_modIsoHcalD1_otherEleVetoForSelf", 70, 0.0, 14.0, "Calo isolation, I_{calo} [GeV]") );
	p.add( tsw::AxisDefn("eleA_modIsoEcal_otherEleVetoForSelf+eleA_modIsoHcalD1_otherEleVetoForSelf", 70, 0.0, 14.0, "Calo isolation (corrected in MC), I_{calo} [GeV]").with_option("var_mc" , "eleA_modIsoEcal_otherEleVetoForSelf+eleA_modIsoHcalD1_otherEleVetoForSelf + tsw::GausianRandom(0.6,0.5)") );

	tsw::DistPlotter p2(p), p3(p), p4(p), p5(p), p6(p);
	p.descriptiveText("Leading ele, 35#leqp_{T}<50GeV");
	p.setSelection("abs(eleA_p4.Eta())<1.45 && abs(eleB_p4.Eta())<1.45 && trgDecision && eleA_p4.Pt()>=35 && eleA_p4.Pt()<50");
	p.outFilePrefix(out_file_base+"35To50GeV");
	p.run();

	p2.descriptiveText("Leading ele, 50#leqp_{T}<75GeV");
	p2.setSelection("abs(eleA_p4.Eta())<1.45 && abs(eleB_p4.Eta())<1.45 && trgDecision && eleA_p4.Pt()>=50 && eleA_p4.Pt()<75");
	p2.outFilePrefix(out_file_base+"50To75GeV");
	p2.run();

	p3.descriptiveText("Leading ele, 75#leqp_{T}<100GeV");
	p3.setSelection("abs(eleA_p4.Eta())<1.45 && abs(eleB_p4.Eta())<1.45 && trgDecision && eleA_p4.Pt()>=75 && eleA_p4.Pt()<100");
	p3.outFilePrefix(out_file_base+"75To100GeV");
	p3.run();

	p4.descriptiveText("Leading ele, 100#leqp_{T}<200GeV");
	p4.setSelection("abs(eleA_p4.Eta())<1.45 && abs(eleB_p4.Eta())<1.45 && trgDecision && eleA_p4.Pt()>=100 && eleA_p4.Pt()<200");
	p4.outFilePrefix(out_file_base+"100To200GeV");
	p4.run();

	p5.descriptiveText("Leading ele, 200#leqp_{T}<300GeV");
	p5.setSelection("abs(eleA_p4.Eta())<1.45 && abs(eleB_p4.Eta())<1.45 && trgDecision && eleA_p4.Pt()>=200 && eleA_p4.Pt()<300");
	p5.outFilePrefix(out_file_base+"200To300GeV");
	p5.run();

	p6.descriptiveText("Leading ele, p_{T}#geq300GeV");
	p6.setSelection("abs(eleA_p4.Eta())<1.45 && abs(eleB_p4.Eta())<1.45 && trgDecision && eleA_p4.Pt()>=300");
	p6.outFilePrefix(out_file_base+"Above300GeV");
	p6.run();

}
