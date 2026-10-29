// Shower/decay an LHE sample with Pythia8 and write HepMC3 for Rivet.
// Compile with the helper script run_gluon_d8_collider_v0_29.sh.
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include "HepMC3/WriterAscii.h"
#include <iostream>

int main(int argc,char** argv){
  if(argc<3){std::cerr<<"usage: shower_gluon_d8_pythia8_v0_29 input.lhe output.hepmc\n"; return 2;}
  Pythia8::Pythia p;
  p.readString("Beams:frameType = 4");
  p.readString(std::string("Beams:LHEF = ")+argv[1]);
  p.readString("HadronLevel:Decay = on");
  p.readString("PartonLevel:ISR = on");
  p.readString("PartonLevel:FSR = on");
  p.readString("PartonLevel:MPI = on");
  p.readString("Random:setSeed = on");
  p.readString("Random:seed = 290529");
  if(!p.init()) return 3;
  HepMC3::Pythia8ToHepMC3 conv;
  HepMC3::WriterAscii writer(argv[2]);
  long n=0;
  while(p.next()){
    auto evt=std::make_shared<HepMC3::GenEvent>();
    conv.fill_next_event(p,evt.get());
    writer.write_event(*evt); ++n;
  }
  writer.close(); p.stat();
  std::cout<<"WROTE "<<n<<" events to "<<argv[2]<<"\n";
  return n>0?0:4;
}
