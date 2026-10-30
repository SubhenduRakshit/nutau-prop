#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: pythia_lhe_to_hepmc_v0_51 input.lhe output.hepmc [seed]\n";
    return 2;
  }
  const std::string input = argv[1];
  const std::string output = argv[2];
  const int seed = argc > 3 ? std::atoi(argv[3]) : 5101;

  Pythia8::Pythia pythia;
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + input);
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = on");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Main:timesAllowErrors = 20");
  if (!pythia.init()) {
    std::cerr << "Pythia init failed\n";
    return 3;
  }

  HepMC::IO_GenEvent ascii_io(output, std::ios::out);
  HepMC::Pythia8ToHepMC to_hepmc;
  long accepted = 0;
  long failed = 0;
  while (true) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break;
      ++failed;
      if (failed > 20) {
        std::cerr << "Too many Pythia failures\n";
        return 4;
      }
      continue;
    }
    auto evt = std::make_unique<HepMC::GenEvent>();
    to_hepmc.fill_next_event(pythia, evt.get());
    ascii_io << evt.get();
    ++accepted;
  }
  pythia.stat();
  std::cout << "PYTHIA_ACCEPTED_EVENTS=" << accepted << "\n";
  std::cout << "PYTHIA_FAILED_EVENTS=" << failed << "\n";
  return accepted > 0 ? 0 : 5;
}
