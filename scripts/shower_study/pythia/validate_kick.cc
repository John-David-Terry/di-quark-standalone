// Check of the SRC_patch_v2 kick options used for Figure 3 (kTdist = 2, kickMode = 1).
// In the gamma* p rest frame, compares the transverse momentum of the whole scattered
// system (status 62 copies, which carry the initiator kT) with that of the remnant d
// (status 63), for valence-u events with a separate u + d remnant. Also counts how often
// hadronization re-merges the remnant u and d into a diquark (status 74).
// Usage: validate_kick <primordial kT widths 0=zero,1=default> <showers 0/1> <src 0/1> <nEvents>
#include "Pythia8/Pythia.h"
#include <cmath>
using namespace Pythia8;

int main(int argc, char* argv[]) {
  int widths = (argc > 1) ? atoi(argv[1]) : 0;
  int shower = (argc > 2) ? atoi(argv[2]) : 0;
  int srcOn  = (argc > 3) ? atoi(argv[3]) : 1;
  int nEvent = (argc > 4) ? atoi(argv[4]) : 20000;

  Pythia pythia(XMLDIR, false);
  pythia.settings.addFlag("SRC:on", false);
  pythia.settings.addFlag("SRC:splitDiquarkDIS", false);
  pythia.settings.addMode("SRC:partnerId", 1, true, true, 1, 2);
  pythia.settings.addMode("SRC:kTdist", 0, true, true, 0, 2);
  pythia.settings.addMode("SRC:kickMode", 0, true, true, 0, 1);
  pythia.settings.addParm("SRC:kT", 1.0, true, false, 0., 0.);
  pythia.settings.addParm("SRC:kbarWidth", 0.2, true, false, 0., 0.);
  for (string c : {"Beams:frameType = 2", "Beams:idA = 11", "Beams:idB = 2212",
    "Beams:eA = 5.", "Beams:eB = 41.", "WeakBosonExchange:ff2ff(t:gmZ) = on",
    "PhaseSpace:Q2Min = 16.", "PDF:lepton = off", "PartonLevel:MPI = off",
    "SpaceShower:dipoleRecoil = on", "TimeShower:QEDshowerByL = off",
    "BeamRemnants:remnantMode = 1", "ColourReconnection:mode = 1",
    "ColourReconnection:reconnect = off", "BeamRemnants:beamJunction = on",
    "SRC:partnerId = 1", "SRC:kTdist = 2", "SRC:kickMode = 1",
    "SRC:kT = 0.4", "SRC:kbarWidth = 0.", "Next:numberCount = 0",
    "Next:numberShowEvent = 0", "Random:setSeed = on", "Random:seed = 777"})
    pythia.readString(c);
  pythia.readString(srcOn ? "SRC:on = on" : "SRC:on = off");
  pythia.readString(shower ? "PartonLevel:ISR = on" : "PartonLevel:ISR = off");
  pythia.readString(shower ? "PartonLevel:FSR = on" : "PartonLevel:FSR = off");
  if (widths == 0) {
    pythia.readString("BeamRemnants:primordialKTsoft = 0.");
    pythia.readString("BeamRemnants:primordialKThard = 0.");
    pythia.readString("BeamRemnants:primordialKTremnant = 0.");
  }
  if (!pythia.init()) return 1;

  long n = 0, nFail = 0, nMerged = 0, nOff = 0, nShow = 0;
  double sumK = 0., maxK = 0., sumCos = 0., maxSum = 0., maxDp = 0.;
  for (int iEv = 0; iEv < nEvent; ++iEv) {
    if (!pythia.next()) { ++nFail; continue; }
    Event& ev = pythia.event;
    Vec4 pIn = ev[1].p() + ev[2].p(), pOut;
    for (int i = 0; i < ev.size(); ++i) if (ev[i].isFinal()) pOut += ev[i].p();
    maxDp = max(maxDp, max(abs((pIn - pOut).e()), (pIn - pOut).pAbs()));
    if (!(pythia.beamB[0].isValence() && pythia.beamB[0].id() == 2)) continue;
    int iLep = -1, iD = -1;
    Vec4 pScat;
    for (int i = 0; i < ev.size(); ++i) {
      // Scattered electron: final copy of the hard-process outgoing lepton.
      if (abs(ev[i].status()) == 23 && ev[i].id() == 11 && iLep < 0) iLep = ev[i].iBotCopyId();
      if (abs(ev[i].status()) == 62) pScat += ev[i].p();
      if (abs(ev[i].status()) == 63 && ev[i].id() == 1 && iD < 0) iD = i;
    }
    if (iD < 0) continue;
    // Remnant d absorbed into a diquark by hadronization (status 74 daughter of the d).
    int iDau = ev[iD].daughter1();
    if (iDau > 0 && abs(ev[iDau].status()) == 74 && ev[iDau].idAbs() > 1000) ++nMerged;
    // Pythia's own hadronic rest frame: remnant vs scattered system before the kT boost
    // (the status 62 copies' mothers), as in BeamRemnants::setOneRemnKinematics.
    Vec4 pScatPre;
    for (int i = 0; i < ev.size(); ++i)
      if (abs(ev[i].status()) == 62) pScatPre += ev[ev[i].mother1()].p();
    Vec4 pRemPre = ev[1].p() + ev[2].p() - ev[iLep].p() - pScatPre;
    RotBstMatrix M;
    M.toCMframe(pScatPre, pRemPre);
    Vec4 a = pScat, b = ev[iD].p();
    a.rotbst(M);
    b.rotbst(M);
    double sumAB = sqrt(pow2(a.px() + b.px()) + pow2(a.py() + b.py()));
    if (sumAB > 0.01) {
      ++nOff;
      if (getenv("SHOW_BAD") && sumAB > atof(getenv("SHOW_BAD")) && nShow++ < 1) {
        printf("BAD: |k_scat|=%.3f |k_d|=%.3f |sum|=%.3f x=%.3f Q2=%.2f\n", a.pT(), b.pT(),
          sumAB, pythia.info.x2(), pythia.info.Q2Fac());
        ev.list();
      }
    }
    ++n;
    sumK += a.pT();
    maxK = max(maxK, a.pT());
    sumCos += (a.px() * b.px() + a.py() * b.py()) / max(1e-12, a.pT() * b.pT());
    maxSum = max(maxSum, sqrt(pow2(a.px() + b.px()) + pow2(a.py() + b.py())));
  }
  printf("widths=%-7s showers=%d src=%d: n=%ld fail=%ld  <|k_scat|>=%.4f  max|k_scat|=%.4f  "
    "<cos(scat,d)>=%+.4f  max|k_scat+k_d|=%.1e (>0.01 GeV in %.2f%%)  d merged into diquark=%.1f%%"
    "  max|dP|=%.1e GeV\n",
    widths ? "default" : "zero", shower, srcOn, n, nFail, sumK / n, maxK, sumCos / n, maxSum,
    100. * nOff / n, 100. * nMerged / n, maxDp);
  return 0;
}
