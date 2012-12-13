{
	tsw::Samples2012 anaTuples("heepTpTree");

	//-> 0) Setup objects
	/*  ---  AXES  ---  */
	tsw::AxisDefn axis_probePt("probe_p4.Pt()", "[35,36,37,38,39,40,41,43,45,47.5,50, 55,60,70,80,90,100, 150,250,500,1000]", "p_{T,probe} [GeV]", -1.0, true);
	tsw::AxisDefn axis_mass("pair_p4.M()", 30, 60., 120., "M_{tag-probe} [GeV]");
	tsw::AxisDefn axis_nVtx("nVtx", 30, 0.5, 30.5, "N_{vtx}");
	tsw::AxisDefn axis_zPt("pair_p4.Pt()", "[0,5,10,15,20,25,30,35,40,50, 60,70,80,90,100, 120,140,160,180,200, 300,400,600,1000]", "p_{T,ee} [GeV]", -1.0, true);
	tsw::AxisDefn axis_dRee("pair_dR", "[0.2,0.35,0.5,0.6,0.7,0.8,1, 1.25,1.5,1.75,2, 2.5,3,3.25]", "#Delta R_{ee}", -1.0, true);

	/*  ---  TAGS  ---  */
	tsw::TagProbeEffiCalc::TagDefn tag_heepStdIso("EB StdHeep, Z peak, probe EB");
	tag_heepStdIso.mSelectionCuts_tag   = pass_elefull("tag_stdHeep");
	tag_heepStdIso.mSelectionCuts_event = "pair_p4.M()<120 && pair_p4.M()>60 && abs(tag_scEta)<1.44 && abs(probe_scEta)<1.44";
	tag_heepStdIso.mFakeRateFunc = "tsw::fr_heepV41_full";
	//	tpEffiCalc.add_tagDefn("EB HEEPStdIso tag, Z peak", , "(tag_stdHeep & 0x0dff)==0 && abs(tag_scEta)<1.44", "tsw::FakeRate_heepStdIso");
	tsw::TagProbeEffiCalc::TagDefn tag_heepStdIso_dRAbove05 = tag_heepStdIso;
	tag_heepStdIso_dRAbove05.mName = "EB StdHeep, Z peak, probe EB, #DeltaR_{ee}>0.4"
	tag_heepStdIso_dRAbove05.mSelectionCuts_event = "pair_p4.M()<120 && pair_p4.M()>60 && abs(tag_scEta)<1.44 && abs(probe_scEta)<1.44 && pair_dR>0.4";

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakZeroFake("Z peak");
	tag_zPeakZeroFake.mSelectionCuts_tag = "";
	tag_zPeakZeroFake.mSelectionCuts_event = "pair_p4.M()<120 && pair_p4.M()>60";
	tag_zPeakZeroFake.mFakeRateFunc = "zero";
	//	tpEffiCalc.add_tagDefn("EB HEEPNoIso tag, Z peak", "pair_p4.M()<120 && pair_p4.M()>60", "(tag_stdHeep & 0x08ff)==0 && abs(tag_scEta)<1.44", "tsw::FakeRate_heepStdIso");

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakZeroFakePhiRdVeto = tag_zPeakZeroFake;
	tag_zPeakZeroFakePhiRdVeto.mName = "Z peak, #phi road veto";
	tag_zPeakZeroFakePhiRdVeto.mSelectionCuts_event = "(pair_p4.M()<120 && pair_p4.M()>60) && (abs(pair_dEta)>=0.07 || abs(pair_dPhi)>=0.3)";

	/*tsw::TagProbeEffiCalc::TagDefn tag_zPeakNarrowZeroFakePhiRdVeto = tag_zPeakZeroFake;
	tag_zPeakNarrowZeroFakePhiRdVeto.mName = "Nrrw Z peak, #phi road veto";
	tag_zPeakNarrowZeroFakePhiRdVeto.mSelectionCuts_event = "(pair_p4.M()<105 && pair_p4.M()>75) && (abs(pair_dEta)>=0.07 || abs(pair_dPhi)>=0.3)";*/

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakNarrowZeroFakePhiRdVetoNonDupl = tag_zPeakZeroFake;
	tag_zPeakNarrowZeroFakePhiRdVetoNonDupl.mName = "Nrrw Z peak, #phi road veto, non-dupl";
	tag_zPeakNarrowZeroFakePhiRdVetoNonDupl.mSelectionCuts_event = "(pair_p4.M()<105 && pair_p4.M()>75) && (abs(pair_dEta)>=0.07 || abs(pair_dPhi)>=0.3) && tag_scEta<probe_scEta";


	/*  ---  PROBES  ---  */

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_stdHeep;
	prb_stdHeep.mEffiDefn = tsw::EffiDefn("((probe_stdHeep & 0x0007)==0) && (abs(probe_scEta)<1.44)", pass_elefull("probe_stdHeep"), "AllStdHEEP, EB", tsw::Black);
	prb_stdHeep.mFakeRateFunc_allProbe = "NO";
	prb_stdHeep.mFakeRateFunc_passProbe = "tsw::fr_heepV41_full";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_heepNoIso;
	prb_heepNoIso.mEffiDefn = tsw::EffiDefn("((probe_stdHeep & 0x0007)==0) && (abs(probe_scEta)<1.44)", pass_eleid("probe_stdHeep"), "HEEP ID, EB", tsw::Black);
	prb_heepNoIso.mFakeRateFunc_allProbe = "NO";
	prb_heepNoIso.mFakeRateFunc_passProbe = "tsw::fr_heepV41_noIso";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_bothHeepModIsoStdThr;
	prb_bothHeepModIsoStdThr.mEffiDefn = tsw::EffiDefn("("+pass_eleid("tag_stdHeep")+" && (abs(tag_scEta)<1.44) && "+pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44) )", pass_elefull("tag_modHeepStdThr")+" && "+pass_elefull("probe_modHeepStdThr"), "ModIsoStdThr, 2-leg, EB", tsw::Black);
	prb_bothHeepModIsoStdThr.mFakeRateFunc_allProbe  = "zero";
	prb_bothHeepModIsoStdThr.mFakeRateFunc_passProbe = "zero";

//	tsw::TagProbeEffiCalc::ProbeCutDefn prb_bothHeepModIsoColThr;
//	prb_bothHeepModIsoColThr.mEffiDefn = tsw::EffiDefn( "("+pass_eleid("tag_stdHeep")+" && (abs(tag_scEta)<1.44) && ("+pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44) )", pass_elefull("tag_modHeepColThr")+" && "+pass_elefull("probe_modHeepColThr"), "ModIsoColThr, 2-leg, EB", tsw::Black);
//	prb_bothHeepModIsoColThr.mFakeRateFunc_allProbe  = "zero";
//	prb_bothHeepModIsoColThr.mFakeRateFunc_passProbe = "zero";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_stdHeepIsoWrtId;
	prb_stdHeepIsoWrtId.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+ " && (abs(probe_scEta)<1.44)", pass_elefull("probe_stdHeep"), "HEEP iso wrt ID, EB", tsw::Black);
	prb_stdHeepIsoWrtId.mFakeRateFunc_allProbe = "tsw::fr_heepV41_noIso";
	prb_stdHeepIsoWrtId.mFakeRateFunc_passProbe = "tsw::fr_heepV41_full";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_phantomIso_dR005To010 = prb_stdHeepIsoWrtId;
	prb_phantomIso_dR005To010.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_005To010, probe_modEmHad1IsoPhantom_005To010, probe_calEaCorr)", "HEEP Iso, phantom, 0.05<#DeltaR<0.1", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_phantomIso_dR010To015 = prb_phantomIso_dR005To010;
	prb_phantomIso_dR010To015.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_010To015, probe_modEmHad1IsoPhantom_010To015, probe_calEaCorr)", "HEEP Iso, phantom, 0.10<#DeltaR<0.15", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_phantomIso_dR015To020 = prb_phantomIso_dR005To010;
	prb_phantomIso_dR015To020.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_015To020, probe_modEmHad1IsoPhantom_015To020, probe_calEaCorr)", "HEEP Iso, phantom, 0.15<#DeltaR<0.2", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_phantomIso_dR020To025 = prb_phantomIso_dR005To010;
	prb_phantomIso_dR020To025.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_020To025, probe_modEmHad1IsoPhantom_020To025, probe_calEaCorr)", "HEEP Iso, phantom, 0.20<#DeltaR<0.25", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_phantomIso_dR025To030 = prb_phantomIso_dR005To010;
	prb_phantomIso_dR025To030.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_025To030, probe_modEmHad1IsoPhantom_025To030, probe_calEaCorr)", "HEEP Iso, phantom, 0.25<#DeltaR<0.3", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_phantomIso_dR030To035 = prb_phantomIso_dR005To010;
	prb_phantomIso_dR030To035.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_030To035, probe_modEmHad1IsoPhantom_030To035, probe_calEaCorr)", "HEEP Iso, phantom, 0.30<#DeltaR<0.35", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_phantomIso_dR035To040 = prb_phantomIso_dR005To010;
	prb_phantomIso_dR035To040.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_035To040, probe_modEmHad1IsoPhantom_035To040, probe_calEaCorr)", "HEEP Iso, phantom, 0.35<#DeltaR<0.4", tsw::Black);


	/* --- TAG-PROBE CALCULATORS --- */

	tsw::TagProbeEffiCalc tpCalc("tagProbeTree", "qcdGsfGsfTree", false);
	tpCalc.baselineSelection("trgDecision");
	tpCalc.eventWeight("genWeight*puWeight");
	tpCalc.stdFakeRate("tsw::fr_heepV41_full", pass_elefull("tag_stdHeep") );
	tpCalc.fakeRatePreSeln("tag_fakePreCutCode==0");
	tpCalc.qcdSysUncertScaleRange( tsw::Range(0.6, 1.4) );

	tpCalc.drellYan( anaTuples.dyEE_mg_single() );
	tpCalc.add_background( anaTuples.wJets() );
	tpCalc.add_background( anaTuples.dyTauTau_powheg() );
	tpCalc.add_background( anaTuples.topBkgds() );
	tpCalc.add_background( anaTuples.vzBkgds()  );

	tpCalc.data( anaTuples.data2012() );

	tpCalc.outFilePrefix("results/20121209/tagProbe/A/tagPrb_dyMG");
	tsw::TagProbeEffiCalc tpCalc_powheg(tpCalc);
	tpCalc_powheg.drellYan( anaTuples.dyEE_powheg() );
	tpCalc_powheg.outFilePrefix("results/20121209/tagProbe/A/tagPrb_dyPowheg");
//	tsw::TagProbeEffiCalc tpCalc_pythia(tpCalc);
//	tpCalc_pythia.drellYan( anaTuples.dyEe_pythia() );
//	tpCalc_pythia.outFilePrefix("results/20121209/tagProbe/A/tagPrb_dyPythia");

	// A) Plots for full standard HEEP [cross-check with Z'(ee) note]
	/*tpCalc.run(tag_heepStdIso, prb_stdHeep, axis_probePt )
					.run(tag_heepStdIso, prb_stdHeep, axis_nVtx)
					.run(tag_heepStdIso, prb_stdHeep, axis_mass);
	tpCalc_powheg.run(tag_heepStdIso, prb_stdHeep, axis_probePt )
					.run(tag_heepStdIso, prb_stdHeep, axis_nVtx)
					.run(tag_heepStdIso, prb_stdHeep, axis_mass);*/
//	tpCalc_pythia.run(tag_heepStdIso, prb_stdHeep, axis_probePt )
//					.run(tag_heepStdIso, prb_stdHeep, axis_nVtx)
//					.run(tag_heepStdIso, prb_stdHeep, axis_mass);

	// B) HEEP ID plots
	/*tpCalc.run(tag_heepStdIso, prb_heepNoIso, axis_probePt )
					.run(tag_heepStdIso, prb_heepNoIso, axis_nVtx )
					.run(tag_heepStdIso, prb_heepNoIso, axis_mass );

	// C) Mod iso plots (probe-probe)
	tpCalc.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoStdThr, axis_zPt)
			.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoStdThr, axis_dRee);*/

	// D) Mod iso plots from phantom ele iso values
	tpCalc.run(tag_heepStdIso_dRAbove05, prb_phantomIso_dR005To010, axis_probePt )
//					.run(tag_heepStdIso_dRAbove05, prb_phantomIso_dR010To015, axis_probePt)
////					.run(tag_heepStdIso_dRAbove05, prb_phantomIso_dR015To020, axis_probePt)
//					.run(tag_heepStdIso_dRAbove05, prb_phantomIso_dR020To025, axis_probePt)
////					.run(tag_heepStdIso_dRAbove05, prb_phantomIso_dR025To030, axis_probePt)
//					.run(tag_heepStdIso_dRAbove05, prb_phantomIso_dR030To035, axis_probePt)
////					.run(tag_heepStdIso_dRAbove05, prb_phantomIso_dR035To040, axis_probePt)
					.run(tag_heepStdIso_dRAbove05, prb_stdHeepIsoWrtId, axis_probePt);



	// --- OLD --- //
//	/*	tpEffiCalc.run(tag_zPeakNarrowZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_zPt);
//	tpEffiCalc.run(tag_zPeakNarrowZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_dRee);
//	tpEffiCalc.run(tag_zPeakZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_zPt);
//	tpEffiCalc.run(tag_zPeakZeroFakePhiRdVeto, prb_bothHeepModIsoStdThr, axis_dRee);
//	tpEffiCalc.run(tag_zPeakZeroFake, prb_bothHeepModIsoStdThr, axis_zPt);
//	tpEffiCalc.run(tag_zPeakZeroFake, prb_bothHeepModIsoStdThr, axis_dRee);*/
//
//	tpEffiCalc.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoColThr, axis_zPt)
//			.run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoColThr, axis_dRee);

}
