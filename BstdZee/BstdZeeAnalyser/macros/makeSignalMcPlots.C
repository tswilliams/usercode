{

   tsw::Samples2012 samples("zEffiTree");
   std::string selnString = "(abs(mcZ_ele1_p4.Eta())<1.44 && abs(mcZ_ele2_p4.Eta())<1.44) && (mcZ_ele1_p4.Et()>35.0 && mcZ_ele2_p4.Et()>35.0)";
   tsw::DistPlotter plotter(true);
   plotter.setTree("zBosonEffiTree");
   plotter.setSelection(selnString);
   plotter.descriptiveText("MC truth after acceptance cuts; gauge interaction.");
   plotter.outFilePrefix("results/20121121/qStarSignalMcAfterAcc_");
   plotter.rescaleMC();

   tsw::MCSample qStarM1500( samples.qStarGI_M1500() );
   qStarM1500.mName = "1.5TeV q*";
   tsw::MCSample qStarM2500( samples.qStarGI_M2500() );
   qStarM2500.mName = "2.5TeV q*";

   plotter.add( qStarM1500 );
   //plotter.add( samples.qStarGI_M2000() );
   plotter.add( qStarM2500 );
   plotter.add( tsw::AxisDefn("ZpT", 40, 400.0, 1400, "Z boson p_{T} [GeV]") );
   plotter.add( tsw::AxisDefn("ZdR", 96, 0.0, 1.2, "#DeltaR_{ee}") );
   plotter.run();

}
