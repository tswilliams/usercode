{
	tsw::Samples2012 anaTuples("heepTpTree");

	std::string OUT_DIR = "results/20130415/TagProbe/";

	//-> 0) Setup objects
	/*  ---  AXES  ---  */
	tsw::AxisDefn axis_probePt("probe_p4.Pt()", "[35,36,37,38,39,40,41,43,45,47.5,50,  55,60,65,70,80,90,100,  150,200,300,500,1000]", "p_{T,probe} [GeV]", 1.0, true);
	tsw::AxisDefn axis_mass("pair_p4.M()", 50, 65., 115., "M_{tag-probe} [GeV]");
	tsw::AxisDefn axis_nVtx("nVtx", 30, 0.5, 30.5, "N_{vtx}");
	tsw::AxisDefn axis_zPt("pair_p4.Pt()", "[0,5,10,15,20,25,30,35,40,50, 60,70,80,90,100, 120,140,160,180,200, 300,400,600,1000]", "p_{T,ee} [GeV]", 1.0, true);
	tsw::AxisDefn axis_dRee("pair_dR", "[0.2,0.3,0.4,0.5,0.6,0.7,0.8,1, 1.25,1.5,1.75,2, 2.5,3,3.25]", "#Delta R_{ee}", 0.1, true);

	/*  ---  TAGS  ---  */
	tsw::TagProbeEffiCalc::TagDefn tag_heepStdIso("EB StdHeep, Z peak, probe EB");
	tag_heepStdIso.mSelectionCuts_tag   = pass_elefull("tag_stdHeep");
	tag_heepStdIso.mSelectionCuts_event = "pair_p4.M()<105 && pair_p4.M()>75 && abs(tag_scEta)<1.44 && abs(probe_scEta)<1.44";
	tag_heepStdIso.mFakeRateFunc = "tsw::fr_heepV41_full";
	//	tpEffiCalc.add_tagDefn("EB HEEPStdIso tag, Z peak", , "(tag_stdHeep & 0x0dff)==0 && abs(tag_scEta)<1.44", "tsw::FakeRate_heepStdIso");
	tsw::TagProbeEffiCalc::TagDefn tag_heepStdIso_dRAbove05 = tag_heepStdIso;
	tag_heepStdIso_dRAbove05.mName = "EB StdHeep, Z peak, probe EB, dR_{ee}>0.4";
	tag_heepStdIso_dRAbove05.mSelectionCuts_event = "pair_p4.M()<105 && pair_p4.M()>75 && abs(tag_scEta)<1.44 && abs(probe_scEta)<1.44 && pair_dR>0.4";

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakZeroFake("Z peak");
	tag_zPeakZeroFake.mSelectionCuts_tag = "";
	tag_zPeakZeroFake.mSelectionCuts_event = "pair_p4.M()<105 && pair_p4.M()>75";
	tag_zPeakZeroFake.mFakeRateFunc = "zero";
	//	tpEffiCalc.add_tagDefn("EB HEEPNoIso tag, Z peak", "pair_p4.M()<105 && pair_p4.M()>75", "(tag_stdHeep & 0x08ff)==0 && abs(tag_scEta)<1.44", "tsw::FakeRate_heepStdIso");

	tsw::TagProbeEffiCalc::TagDefn tag_zPeakZeroFakePhiRdVeto = tag_zPeakZeroFake;
	tag_zPeakZeroFakePhiRdVeto.mName = "Z peak, #phi road veto";
	tag_zPeakZeroFakePhiRdVeto.mSelectionCuts_event = "(pair_p4.M()<105 && pair_p4.M()>75) && (abs(pair_dEta)>=0.07 || abs(pair_dPhi)>=0.3)";

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
//	prb_stdHeepIsoWrtId.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+ " && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_stdTrkIso, probe_stdEmHad1Iso + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)", "HEEP iso wrt ID from func, EB", tsw::Black);
	prb_stdHeepIsoWrtId.mFakeRateFunc_allProbe = "tsw::fr_heepV41_noIso";
	prb_stdHeepIsoWrtId.mFakeRateFunc_passProbe = "tsw::fr_heepV41_full";
	/*tsw::TagProbeEffiCalc::ProbeCutDefn prb_stdHeepIsoWrtIdWFunc;
	prb_stdHeepIsoWrtIdWFunc.mEffiDefn = tsw::EffiDefn(pass_eleid("probe_stdHeep")+ " && (abs(probe_scEta)<1.44)", "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_stdTrkIso, probe_stdEmHad1Iso, probe_calEaCorr)", "HEEP iso wrt ID from func, EB", tsw::Black);
	prb_stdHeepIsoWrtIdWFunc.mFakeRateFunc_allProbe = "zero";//"tsw::fr_heepV41_noIso";
	prb_stdHeepIsoWrtIdWFunc.mFakeRateFunc_passProbe = "zero";//"tsw::fr_heepV41_full";*/

   // PHANTOM ISO plots //
	std::string cut_mimickIso_dR005To010 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_005To010, probe_modEmHad1IsoPhantom_005To010, probe_calEaCorr)";
	std::string cut_mimickIso_dR010To015 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_010To015, probe_modEmHad1IsoPhantom_010To015, probe_calEaCorr)";
	std::string cut_mimickIso_dR015To020 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_015To020, probe_modEmHad1IsoPhantom_015To020, probe_calEaCorr)";
	std::string cut_mimickIso_dR020To025 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_020To025, probe_modEmHad1IsoPhantom_020To025, probe_calEaCorr)";
	std::string cut_mimickIso_dR025To030 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_025To030, probe_modEmHad1IsoPhantom_025To030, probe_calEaCorr)";
	std::string cut_mimickIso_dR030To035 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_030To035, probe_modEmHad1IsoPhantom_030To035, probe_calEaCorr)";
	std::string cut_mimickIso_dR035To040 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_035To040, probe_modEmHad1IsoPhantom_035To040, probe_calEaCorr)";

	std::string cut_mimickIsoSmeared_dR005To010 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_005To010, probe_modEmHad1IsoPhantom_005To010 + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)";
	std::string cut_mimickIsoSmeared_dR010To015 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_010To015, probe_modEmHad1IsoPhantom_010To015 + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)";
	std::string cut_mimickIsoSmeared_dR015To020 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_015To020, probe_modEmHad1IsoPhantom_015To020 + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)";
	std::string cut_mimickIsoSmeared_dR020To025 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_020To025, probe_modEmHad1IsoPhantom_020To025 + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)";
	std::string cut_mimickIsoSmeared_dR025To030 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_025To030, probe_modEmHad1IsoPhantom_025To030 + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)";
	std::string cut_mimickIsoSmeared_dR030To035 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_030To035, probe_modEmHad1IsoPhantom_030To035 + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)";
	std::string cut_mimickIsoSmeared_dR035To040 = "tsw::heep_iso_cut(probe_p4.Pt(), probe_scEta, probe_modTrkIsoPhantom_035To040, probe_modEmHad1IsoPhantom_035To040 + tsw::GausianRandom(0.6,0.5), probe_calEaCorr)";

	std::string cut_prbHeepId = pass_eleid("probe_stdHeep");


	tsw::TagProbeEffiCalc::ProbeCutDefn prb_mimickIso_dR005To010 = prb_stdHeepIsoWrtId;
	prb_mimickIso_dR005To010.mEffiDefn       =     tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIso_dR005To010,        "HEEP Iso, phantom, 0.05<#DeltaR<0.1", tsw::Black);
	prb_mimickIso_dR005To010.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIsoSmeared_dR005To010, "HEEP Iso, phantom smeared, 0.05<#DeltaR<0.1", tsw::Black);
	prb_mimickIso_dR005To010.mFakeRateFunc_allProbe  = "zero";
	prb_mimickIso_dR005To010.mFakeRateFunc_passProbe = "zero";

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_mimickIso_dR010To015 = prb_mimickIso_dR005To010;
	prb_mimickIso_dR010To015.mEffiDefn       =     tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIso_dR010To015, "HEEP Iso, phantom, 0.10<#DeltaR<0.15", tsw::Black);
	prb_mimickIso_dR010To015.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIsoSmeared_dR010To015, "HEEP Iso, phantom smeared, 0.10<#DeltaR<0.15", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_mimickIso_dR015To020 = prb_mimickIso_dR005To010;
	prb_mimickIso_dR015To020.mEffiDefn       =     tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIso_dR015To020, "HEEP Iso, phantom, 0.15<#DeltaR<0.2", tsw::Black);
	prb_mimickIso_dR015To020.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIsoSmeared_dR015To020, "HEEP Iso, phantom smeared, 0.15<#DeltaR<0.2", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_mimickIso_dR020To025 = prb_mimickIso_dR005To010;
	prb_mimickIso_dR020To025.mEffiDefn       =     tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIso_dR020To025, "HEEP Iso, phantom, 0.20<#DeltaR<0.25", tsw::Black);
	prb_mimickIso_dR020To025.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIsoSmeared_dR020To025, "HEEP Iso, phantom smeared, 0.20<#DeltaR<0.25", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_mimickIso_dR025To030 = prb_mimickIso_dR005To010;
	prb_mimickIso_dR025To030.mEffiDefn       =     tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIso_dR025To030, "HEEP Iso, phantom, 0.25<#DeltaR<0.3", tsw::Black);
	prb_mimickIso_dR025To030.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIsoSmeared_dR025To030, "HEEP Iso, phantom smeared, 0.25<#DeltaR<0.3", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_mimickIso_dR030To035 = prb_mimickIso_dR005To010;
	prb_mimickIso_dR030To035.mEffiDefn       =     tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIso_dR030To035, "HEEP Iso, phantom, 0.30<#DeltaR<0.35", tsw::Black);
	prb_mimickIso_dR030To035.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIsoSmeared_dR030To035, "HEEP Iso, phantom smeared, 0.30<#DeltaR<0.35", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_mimickIso_dR035To040 = prb_mimickIso_dR005To010;
	prb_mimickIso_dR035To040.mEffiDefn       =     tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIso_dR035To040, "HEEP Iso, phantom, 0.35<#DeltaR<0.4", tsw::Black);
	prb_mimickIso_dR035To040.mMcEffiDefn4Sys = new tsw::EffiDefn(pass_eleid("probe_stdHeep")+" && (abs(probe_scEta)<1.44)", cut_mimickIsoSmeared_dR035To040, "HEEP Iso, phantom smeared, 0.35<#DeltaR<0.4", tsw::Black);



	tsw::TagProbeEffiCalc::ProbeCutDefn prb_idMimickIso_dR005To010 = prb_mimickIso_dR005To010;
	prb_idMimickIso_dR005To010.mEffiDefn       =     tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIso_dR005To010),        "ID and Iso, phantom, 0.05<#DeltaR<0.10", tsw::Black);
	prb_idMimickIso_dR005To010.mMcEffiDefn4Sys = new tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIsoSmeared_dR005To010), "ID and Iso, phantom smeared, 0.05<#DeltaR<0.10", tsw::Black);
	prb_idMimickIso_dR005To010.mFakeRateFunc_allProbe  = "NO";
	prb_idMimickIso_dR005To010.mFakeRateFunc_passProbe = "tsw::fr_heepV41_full";


	tsw::TagProbeEffiCalc::ProbeCutDefn prb_idMimickIso_dR010To015 = prb_idMimickIso_dR005To010;
	prb_idMimickIso_dR010To015.mEffiDefn       =     tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIso_dR010To015),        "ID and Iso, phantom, 0.10<#DeltaR<0.15", tsw::Black);
	prb_idMimickIso_dR010To015.mMcEffiDefn4Sys = new tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIsoSmeared_dR010To015), "ID and Iso, phantom smeared, 0.10<#DeltaR<0.15", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_idMimickIso_dR015To020 = prb_idMimickIso_dR005To010;
	prb_idMimickIso_dR015To020.mEffiDefn       =     tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIso_dR015To020),        "ID and Iso, phantom, 0.15<#DeltaR<0.20", tsw::Black);
	prb_idMimickIso_dR015To020.mMcEffiDefn4Sys = new tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIsoSmeared_dR015To020), "ID and Iso, phantom smeared, 0.15<#DeltaR<0.20", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_idMimickIso_dR020To025 = prb_idMimickIso_dR005To010;
	prb_idMimickIso_dR020To025.mEffiDefn       =     tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIso_dR020To025),        "ID and Iso, phantom, 0.20<#DeltaR<0.25", tsw::Black);
	prb_idMimickIso_dR020To025.mMcEffiDefn4Sys = new tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIsoSmeared_dR020To025), "ID and Iso, phantom smeared, 0.20<#DeltaR<0.25", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_idMimickIso_dR025To030 = prb_idMimickIso_dR005To010;
	prb_idMimickIso_dR025To030.mEffiDefn       =     tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIso_dR025To030),        "ID and Iso, phantom, 0.25<#DeltaR<0.30", tsw::Black);
	prb_idMimickIso_dR025To030.mMcEffiDefn4Sys = new tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIsoSmeared_dR025To030), "ID and Iso, phantom smeared, 0.25<#DeltaR<0.30", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_idMimickIso_dR030To035 = prb_idMimickIso_dR005To010;
	prb_idMimickIso_dR030To035.mEffiDefn       =     tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIso_dR030To035),        "ID and Iso, phantom, 0.30<#DeltaR<0.35", tsw::Black);
	prb_idMimickIso_dR030To035.mMcEffiDefn4Sys = new tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIsoSmeared_dR030To035), "ID and Iso, phantom smeared, 0.30<#DeltaR<0.35", tsw::Black);

	tsw::TagProbeEffiCalc::ProbeCutDefn prb_idMimickIso_dR035To040 = prb_idMimickIso_dR005To010;
	prb_idMimickIso_dR035To040.mEffiDefn       =     tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIso_dR035To040),        "ID and Iso, phantom, 0.35<#DeltaR<0.40", tsw::Black);
	prb_idMimickIso_dR035To040.mMcEffiDefn4Sys = new tsw::EffiDefn("(probe_stdHeep & 0x007)==0 && (abs(probe_scEta)<1.44)", tsw::AndOfCuts(cut_prbHeepId, cut_mimickIsoSmeared_dR035To040), "ID and Iso, phantom smeared, 0.35<#DeltaR<0.40", tsw::Black);
	
	

	/* --- TAG-PROBE CALCULATORS --- */

	tsw::TagProbeEffiCalc tpCalc("tagProbeTree", "qcdGsfGsfTree", true);
	tpCalc.baselineSelection("trgDecision");
	tpCalc.eventWeight("genWeight*puWeight");
	tpCalc.stdFakeRate("tsw::fr_heepV41_full", pass_elefull("tag_stdHeep") );
	tpCalc.fakeRatePreSeln("tag_fakePreCutCode==0");
  	tpCalc.qcdSysUncertScaleRange( tsw::Range(0.6, 1.4) );

	tpCalc.drellYan( anaTuples.dyEE_mg_merged() );
	tpCalc.add_background( anaTuples.wJets() );
	tpCalc.add_background( anaTuples.dyTauTau_powheg() );
	tpCalc.add_background( anaTuples.topBkgds() );
	tpCalc.add_background( anaTuples.vzBkgds()  );

	tpCalc.data( anaTuples.data2012() );

	tpCalc.outFilePrefix(OUT_DIR+"/tagPrb_dyMG");
	tsw::TagProbeEffiCalc tpCalc_powheg(tpCalc);
	tpCalc_powheg.drellYan( anaTuples.dyEE_powheg() );
	tpCalc_powheg.outFilePrefix(OUT_DIR+"/tagPrb_dyPowheg");
	tsw::TagProbeEffiCalc tpCalc_pythia(tpCalc);
	tpCalc_pythia.drellYan( anaTuples.dyEE_pythia() );
	tpCalc_pythia.outFilePrefix(OUT_DIR+"tagPrb_dyPythia");

	// A) Plots for full standard HEEP [cross-check with Z'(ee) note]
	tpCalc.run(tag_heepStdIso, prb_stdHeep, axis_probePt )
	  .run(tag_heepStdIso, prb_stdHeep, axis_nVtx)
	  .run(tag_heepStdIso, prb_stdHeep, axis_mass);
	tpCalc_powheg.run(tag_heepStdIso, prb_stdHeep, axis_probePt )
	   			.run(tag_heepStdIso, prb_stdHeep, axis_nVtx)
	   			.run(tag_heepStdIso, prb_stdHeep, axis_mass);
	tpCalc_pythia.run(tag_heepStdIso, prb_stdHeep, axis_probePt )
	  			.run(tag_heepStdIso, prb_stdHeep, axis_nVtx)
	      			.run(tag_heepStdIso, prb_stdHeep, axis_mass);

//	// B) HEEP ID plots
	tpCalc.run(tag_heepStdIso, prb_heepNoIso, axis_probePt )
	  			.run(tag_heepStdIso, prb_heepNoIso, axis_nVtx )
	      			.run(tag_heepStdIso, prb_heepNoIso, axis_mass );

	// C) Mod iso plots (probe-probe)
	tpCalc.effiLegendPos(tsw::LOWER_LEFT)
         .run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoStdThr, axis_zPt)
         .effiLegendPos(tsw::LOWER_RIGHT)
         .histLegendPos(tsw::UPPER_LEFT)
         .run(tag_zPeakNarrowZeroFakePhiRdVetoNonDupl, prb_bothHeepModIsoStdThr, axis_dRee)
         .histLegendPos(tsw::UPPER_RIGHT);

	// D) Mod iso plots from phantom ele iso values
	tpCalc.effiLegendPos(tsw::LOWER_LEFT);
	tpCalc.descriptiveText("0.05<#DeltaR<0.1"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_mimickIso_dR005To010, axis_probePt );
	tpCalc.descriptiveText("0.1<#DeltaR<0.15"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_mimickIso_dR010To015, axis_probePt);
	tpCalc.descriptiveText("0.15<#DeltaR<0.2"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_mimickIso_dR015To020, axis_probePt);
	tpCalc.descriptiveText("0.2<#DeltaR<0.25"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_mimickIso_dR020To025, axis_probePt);
	tpCalc.descriptiveText("0.25<#DeltaR<0.3"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_mimickIso_dR025To030, axis_probePt);
	tpCalc.descriptiveText("0.3<#DeltaR<0.35"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_mimickIso_dR030To035, axis_probePt);
	tpCalc.descriptiveText("0.35<#DeltaR<0.4"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_mimickIso_dR035To040, axis_probePt);
	tpCalc.descriptiveText("Standard isolation values"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_stdHeepIsoWrtId, axis_probePt);

	tpCalc.descriptiveText("ID + iso, 0.05<#DeltaR<0.1"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_idMimickIso_dR005To010, axis_probePt);
	tpCalc.descriptiveText("ID + iso, 0.1<#DeltaR<0.15"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_idMimickIso_dR010To015, axis_probePt);
	tpCalc.descriptiveText("ID + iso, 0.15<#DeltaR<0.2"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_idMimickIso_dR015To020, axis_probePt);
	tpCalc.descriptiveText("ID + iso, 0.2<#DeltaR<0.25"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_idMimickIso_dR020To025, axis_probePt);
	tpCalc.descriptiveText("ID + iso, 0.25<#DeltaR<0.3"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_idMimickIso_dR025To030, axis_probePt);
	tpCalc.descriptiveText("ID + iso, 0.3<#DeltaR<0.35"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_idMimickIso_dR030To035, axis_probePt);
	tpCalc.descriptiveText("ID + iso, 0.35<#DeltaR<0.4"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_idMimickIso_dR035To040, axis_probePt);
	tpCalc.descriptiveText("Standard isolation values");  tpCalc.run(tag_heepStdIso_dRAbove05, prb_stdHeep, axis_probePt);

	
	tpCalc.descriptiveText("Standard isolation vals w func"); tpCalc.run(tag_heepStdIso_dRAbove05, prb_stdHeepIsoWrtIdWFunc, axis_probePt);


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
