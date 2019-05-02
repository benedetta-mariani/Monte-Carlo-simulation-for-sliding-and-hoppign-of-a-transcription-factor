void compilemyclass(TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  gSystem->CompileMacro("RandomWalker.cxx",opt.Data());
  gSystem->CompileMacro("Data.cxx",opt.Data());
  gSystem->CompileMacro("TFactor.cxx",opt.Data());
  gSystem->CompileMacro("Simulation.cpp",opt.Data());
}
