{
	// SIGNAL DATASETS

	tsw::CompositeMC qStar_GI("q*, GI", tsw::Lilac, "M_{q*}", "TeV");
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M500_prodReco_heepTpTree.root", "GI, 0.5TeV", tsw::Green), 0.5 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M600_prodReco_heepTpTree.root", "GI, 0.6TeV", tsw::Green), 0.6 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M1500_prodReco_heepTpTree.root", "GI, 1.5TeV", tsw::Green), 1.5 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M1600_prodReco_heepTpTree.root", "GI, 1.6TeV", tsw::Green), 1.6 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M1900_prodReco_heepTpTree.root", "GI, 1.9TeV", tsw::Green), 1.9 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M2000_prodReco_heepTpTree.root", "GI, 2.0TeV", tsw::Green), 2.0 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M2500_prodReco_heepTpTree.root", "GI, 2.5TeV", tsw::Green), 2.5 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M2700_prodReco_heepTpTree.root", "GI, 2.7TeV", tsw::Green), 2.7 );
	qStar_GI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarGI_M2900_prodReco_heepTpTree.root", "GI, 2.9TeV", tsw::Green), 2.9 );

	tsw::CompositeMC qStar_CI("q*, CI", tsw::LightGreen, "M_{q*}", "TeV");
	qStar_CI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarCI_M1300_prodReco_heepTpTree.root", "CI, 1.3TeV", tsw::Blue), 1.3 );
	qStar_CI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarCI_M1500_prodReco_heepTpTree.root", "CI, 1.5TeV", tsw::Blue), 1.5 );
	qStar_CI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarCI_M1700_prodReco_heepTpTree.root", "CI, 1.7TeV", tsw::Blue), 1.7 );
	qStar_CI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarCI_M2600_prodReco_heepTpTree.root", "CI, 2.6TeV", tsw::Blue), 2.6 );
	qStar_CI.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/QstarCI_M2900_prodReco_heepTpTree.root", "CI, 2.9TeV", tsw::Blue), 2.9 );

	// DY DATASETS

	tsw::MCSample dyEe_mg_single("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DYJetsToEE-MG_?_heepTpTree.root", "DY(ee), MG", tsw::LightBlue, /*TWiki NNLO*/(3503.7e3)/3.0 );

	tsw::CompositeMC dyEe_mg_merged("DY(ee), MG", tsw::LightBlue);
	dyEe_mg_merged.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DYJetsToEE-MG_?_heepTpTree.root", "MADGRAPH", tsw::Blue, /*TWiki NNLO*/(3503.7e3)/3.0, "pair_p4.Pt()<240" ) , 0.0);
	dyEe_mg_merged.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DYJetsToEE-Pt100-MG_heepTpTree.root", "MG pT>100", tsw::Blue, /*TWiki NNLO*/(3503.7e3)*(34.1/2950.0)/3.0,  "pair_p4.Pt()>240" ) , 1.0);

	tsw::MCSample dyEe_powheg("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DYToEE-powheg_heepTpTree.root",  "DY(ee), POWHEG", tsw::Blue, /*TWiki NNLO*/(5745.3e3)/3.0 );
	tsw::MCSample dyEe_pythia("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DYToEE-pythia6_?_heepTpTree.root", "DY(ee), PYTHIA", tsw::Blue, /*TWiki NNLO*/(5745.3e3)/3.0 );

	tsw::CompositeMC dyEe_calchep("DY(ee), CalcHEP", tsw::LightBlue, "p_{T,ee}", "GeV");
	dyEe_calchep.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/DYJetsToEE-calchep-PtZ100To300_prodReco_heepTpTree.root", "DY(ee), CalcHEP, 100 #leq p_{T} < 300", tsw::LightBlue, 18582.6/2.), 200. );
	dyEe_calchep.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/DYJetsToEE-calchep-PtZ300To500_prodReco_heepTpTree.root", "DY(ee), CalcHEP, 300 #leq p_{T} < 500", tsw::LightBlue, 294.645/2.), 400. );
	dyEe_calchep.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/DYJetsToEE-calchep-PtZ500To700_prodReco_heepTpTree.root", "DY(ee), CalcHEP, 500 #leq p_{T} < 700", tsw::LightBlue, 22.744/2.), 600. );
	dyEe_calchep.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/DYJetsToEE-calchep-PtZ700To900_prodReco_heepTpTree.root", "DY(ee), CalcHEP, 700 #leq p_{T} < 900", tsw::LightBlue, 3.055/2.), 800. );
	dyEe_calchep.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/int/DYJetsToEE-calchep-PtZ900_prodReco_heepTpTree.root", "DY(ee), CalcHEP, p_{T} #geq 900", tsw::LightBlue, 0.668/2.), 1000. );

	// OTHER BACKGROUNDS

	tsw::CompositeMC topBkgds("t#bar{t}, tW, #bar{t}W & WW", tsw::Red);
	topBkgds.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/TTbarJets_?_heepTpTree.root", "ttbar", tsw::Red, /*TWiki NLO*/225.2e3), 0.0);
	topBkgds.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/TW_heepTpTree.root", "tW", tsw::Red, 11.1e3/*TWiki approx NNLO*/), 0.1);
	topBkgds.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/TbarW_heepTpTree.root", "#bar{t}W", tsw::Red, 11.1e3/*TWiki approx NNLO*/), 0.2);
	topBkgds.add( tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/WW_?_heepTpTree.root", "WW", tsw::Red, 54.8e3/*TWiki NLO*/), 1.0 );

	tsw::CompositeMC vzBkgds("ZZ & ZW", tsw::Green);
	vzBkgds.add(  tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/ZZ_?_heepTpTree.root", "ZZ", tsw::Black, 8.3e3/*TWiki NLO*/), 0.0  );
	vzBkgds.add(  tsw::MCSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/WZ_?_heepTpTree.root", "WZ", tsw::Black, 22.0e3/*Z'(ee) ICHEP */), 1.0  );

	tsw::MCSample wJets_mg("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/WJetsToLNu_?_heepTpTree.root", "W(l#nu)+Jets", tsw::Blue, 36257.2e3/*TWiki NNLO*/);
	tsw::MCSample dyTauTau_mg("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DYJetsToTauTau-MG_?_heepTpTree.root", "DY(#tau#tau)", tsw::Pink, (3503.7e3)/3.0/*TWiki NNLO*/);

	// DATA

	tsw::CompositeData data2012("Data", tsw::Black);
	data2012.add( tsw::DataSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/Photon-12A-PromptV1_*_heepTpTree.root", "Run12A", tsw::Black, 8.4/3., "trgDecision"), 1.0 );
	data2012.add( tsw::DataSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DblPhotonHighPt-12B-PromptV1_*_heepTpTree.root", "Run12B", tsw::Black, 8.4/3., "trgDecision"), 2.0 );
	data2012.add( tsw::DataSample("/opt/ppd/newscratch/williams/Datafiles/AnaTuples/bstdZeeAna_8TeV_20120913/batch/DblPhotonHighPt-12C-PromptV*_heepTpTree.root", "Run12C", tsw::Black, 8.4/3., "trgDecision"), 3.0 );


	//-> 0) Setup objects
	/*  ---  AXES  ---  */
	double ptBinLims[] = {35., 36., 37., 38., 39., 40., 41., 43., 45., 47.5, 50.,
									55., 60., 70., 80., 90., 100., 150., 250., 500.};
	tsw::AxisDefn axis_probePt("probe_p4.Pt()", 19, ptBinLims, "p_{T,probe} [GeV]", -1.0, true);
	tsw::AxisDefn axis_nVtx("nVtx", 18, 0.5, 18.5, "N_{vtx}");
	double zPtBinLims[] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.,
									120., 140., 160., 180., 200., 300., 400., 600., 1000.};
	tsw::AxisDefn axis_zPt("pair_p4.Pt()", 23, zPtBinLims, "p_{T,ee} [GeV]", -1.0, true);
	double dRBinLims[] = {0.2, 0.35, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.25};
	tsw::AxisDefn axis_dRee("pair_dR", 13, dRBinLims, "#Delta R", -1.0, true);

	/*  ---  TAGS  ---  */
	tsw::TagProbeEffiCalc::TagDefn tag_heepStdIso("EB HeepStdIso tag, Z peak");
	tag_heepStdIso.mSelectionCuts_tag   = "(tag_stdHeep & 0x0dff)==0 && abs(tag_scEta)<1.44";
	tag_heepStdIso.mSelectionCuts_event = "pair_p4.M()<120 && pair_p4.M()>60";
	tag_heepStdIso.mFakeRateFunc = "tsw::FakeRate_heepStdIso";
	//	tpEffiCalc.add_tagDefn("EB HEEPStdIso tag, Z peak", , "(tag_stdHeep & 0x0dff)==0 && abs(tag_scEta)<1.44", "tsw::FakeRate_heepStdIso");

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakZeroFake("Z peak");
	tag_zPeakZeroFake.mSelectionCuts_tag = "";
	tag_zPeakZeroFake.mSelectionCuts_event = "pair_p4.M()<120 && pair_p4.M()>60";
	tag_zPeakZeroFake.mFakeRateFunc = "zero";
	//	tpEffiCalc.add_tagDefn("EB HEEPNoIso tag, Z peak", "pair_p4.M()<120 && pair_p4.M()>60", "(tag_stdHeep & 0x08ff)==0 && abs(tag_scEta)<1.44", "tsw::FakeRate_heepStdIso");

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakZeroFakePhiRdVeto = tag_zPeakZeroFake;
	tag_zPeakZeroFakePhiRdVeto.mName = "Z peak, #phi road veto";
	tag_zPeakZeroFakePhiRdVeto.mSelectionCuts_event = "(pair_p4.M()<120 && pair_p4.M()>60) && (abs(pair_dEta)>=0.07 || abs(pair_dPhi)>=0.3)";

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakNarrowZeroFakePhiRdVeto = tag_zPeakZeroFake;
	tag_zPeakNarrowZeroFakePhiRdVeto.mName = "Narrow Z peak, #phi road veto";
	tag_zPeakNarrowZeroFakePhiRdVeto.mSelectionCuts_event = "(pair_p4.M()<105 && pair_p4.M()>75) && (abs(pair_dEta)>=0.07 || abs(pair_dPhi)>=0.3)";

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakNarrowZeroFakePhiRdVetoNonDupl = tag_zPeakZeroFake;
	tag_zPeakNarrowZeroFakePhiRdVetoNonDupl.mName = "Narrow Z peak, #phi road veto, non-dupl";
	tag_zPeakNarrowZeroFakePhiRdVetoNonDupl.mSelectionCuts_event = "(pair_p4.M()<105 && pair_p4.M()>75) && (abs(pair_dEta)>=0.07 || abs(pair_dPhi)>=0.3) && tag_scEta<probe_scEta";

	/*  ---  PROBES  ---  */
	tsw::TagProbeEffiCalc::ProbeCutDefn prb_heepStdIso;
	prb_heepStdIso.mEffiDefn = tsw::EffiDefn("((probe_stdHeep & 0x0007)==0) && (abs(probe_scEta)<1.44)", "(probe_stdHeep & 0x0dff)==0", "AllStdHEEP, EB", tsw::Black);
	prb_heepStdIso.mFakeRateFunc_allProbe = "NO";
	prb_heepStdIso.mFakeRateFunc_passProbe = "tsw::FakeRate_heepStdIso";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_heepNoIso;
	prb_heepNoIso.mEffiDefn = tsw::EffiDefn("((probe_stdHeep & 0x0007)==0) && (abs(probe_scEta)<1.44)", "(probe_stdHeep & 0x08ff)==0", "HEEP ID, EB", tsw::Black);
	prb_heepNoIso.mFakeRateFunc_allProbe = "NO";
	prb_heepNoIso.mFakeRateFunc_passProbe = "tsw::FakeRate_heepStdIso";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_bothHeepModIsoStdThr;
	prb_bothHeepModIsoStdThr.mEffiDefn = tsw::EffiDefn("((tag_stdHeep & 0x08ff)==0) && (abs(tag_scEta)<1.44) && ((probe_stdHeep & 0x08ff)==0) && (abs(probe_scEta)<1.44)", "(tag_modHeepStdThr & 0x0dff)==0 && (probe_modHeepStdThr & 0x0dff)==0", "ModIsoStdThr, 2-leg, EB", tsw::Black);
	prb_bothHeepModIsoStdThr.mFakeRateFunc_allProbe = "zero";
	prb_bothHeepModIsoStdThr.mFakeRateFunc_passProbe = "zero";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_bothHeepModIsoColThr;
	prb_bothHeepModIsoColThr.mEffiDefn = tsw::EffiDefn("((tag_stdHeep & 0x08ff)==0) && (abs(tag_scEta)<1.44) && ((probe_stdHeep & 0x08ff)==0) && (abs(probe_scEta)<1.44)", "(tag_modHeepColThr & 0x0dff)==0 && (probe_modHeepColThr & 0x0dff)==0", "ModIsoColThr, 2-leg, EB", tsw::Black);
	prb_bothHeepModIsoColThr.mFakeRateFunc_allProbe = "zero";
	prb_bothHeepModIsoColThr.mFakeRateFunc_passProbe = "zero";


	//-> 1) Run TagProbeEffiCalc

	tsw::TagProbeEffiCalc tpEffiCalc("tagProbeTree", "qcdGsfGsfTree");
	tpEffiCalc.baselineSelection("1");
	tpEffiCalc.eventWeight("genWeight*puWeight");
	tpEffiCalc.stdFakeRate("tsw::FakeRate_heepStdIso", "(tag_stdHeep & 0x08ff)==0");
	tpEffiCalc.fakeRatePreSeln("tag_fakePreCutCode==0");

	tpEffiCalc.drellYan( dyEe_mg_merged );
	tpEffiCalc.add_background( wJets_mg );
	tpEffiCalc.add_background( dyTauTau_mg );
	tpEffiCalc.add_background( topBkgds );
	tpEffiCalc.add_background( vzBkgds  );

	tpEffiCalc.data( data2012 );

	tpEffiCalc.outFilePrefix("results/20120926/tagProbe/tagProbe_dyMG");
	tsw::TagProbeEffiCalc tpEffiCalc_powheg(tpEffiCalc);
	tpEffiCalc_powheg.drellYan( dyEe_powheg );
	tpEffiCalc_powheg.outFilePrefix("results/20120926/tagProbe/tagProbe_dyPowheg");
	tsw::TagProbeEffiCalc tpEffiCalc_pythia(tpEffiCalc);
	tpEffiCalc_pythia.drellYan( dyEe_pythia );
	tpEffiCalc_pythia.outFilePrefix("results/20120926/tagProbe/tagProbe_dyPythia");

	tpEffiCalc.run(tag_heepStdIso, prb_heepStdIso, axis_probePt )
			.run(tag_heepStdIso, prb_heepStdIso, axis_nVtx);
	tpEffiCalc_powheg.run(tag_heepStdIso, prb_heepStdIso, axis_probePt )
			.run(tag_heepStdIso, prb_heepStdIso, axis_nVtx);
	tpEffiCalc_pythia.run(tag_heepStdIso, prb_heepStdIso, axis_probePt )
			.run(tag_heepStdIso, prb_heepStdIso, axis_nVtx);

	tpEffiCalc.run(tag_heepStdIso, prb_heepNoIso, axis_probePt )
			.run(tag_heepStdIso, prb_heepNoIso, axis_nVtx )
			.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoStdThr, axis_zPt)
			.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoStdThr, axis_dRee);
	/*	tpEffiCalc.run(tag_zPeakNarrowZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_zPt);
	tpEffiCalc.run(tag_zPeakNarrowZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_dRee);
	tpEffiCalc.run(tag_zPeakZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_zPt);
	tpEffiCalc.run(tag_zPeakZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_dRee);
	tpEffiCalc.run(tag_zPeakZeroFake, prb_bothHeepModIsoStdThr, axis_zPt);
	tpEffiCalc.run(tag_zPeakZeroFake, prb_bothHeepModIsoStdThr, axis_dRee);*/

	tpEffiCalc.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoColThr, axis_zPt)
			.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoColThr, axis_dRee);

}
