// Figure 3 observable R_{pi-/pi+}(phi_hJ) generated with the SRC-patched Pythia 8.311.
//
// Differences from the original Track A pipeline (generate -> surgery v2 -> reinject):
//  - the SRC kick is applied inside Pythia's beam-remnant step (SRC_patch_v2.diff), so the
//    order is kick -> shower -> hadronize in one pythia.next();
//  - the remnant ud diquark is broken up with the stock junction remnant model
//    (BeamRemnants:remnantMode = 1, beamJunction = on) instead of (ud) -> d + (us) + sbar;
//  - the jet is a Breit-frame e+e- anti-kT jet (fjcore ee_genkt, p = -1), not the struck parton.
//    The struck-parton "jet" of the original analysis is also written (phi_hJ_parton).
// Everything else follows Figure 3: 5 x 41 GeV, x > 0.05, 16 < Q^2 < 25 GeV^2, struck u,
// x <= 0.495 (surgery-v2 rejection), leading charged pion in the Breit target hemisphere
// with proton-forward lab 2 < eta < 3.5, z_LC >= 0.1.
//
// Output rows use the column names of jet_hadron_pi_pm_observables.csv, so the existing
// plot_jet_hadron_pi_pm_observables.py can read them unchanged.
//
// Usage: fig3_patched key=value ...
//   out=<csv>  sample=unchanged|altered  src=0|1  showers=0|1  remnant=junction|diquark
//   nacc=<SRC-eligible events>  seed=<int>  eidoffset=<int>  jetR=1.0  kt=0.4  xmldoc=<dir>
//   ktdist=0|1|2: SRC:kTdist, the law of the back-to-back pair kick k (+k to the struck u, -k to
//     the remnant d, on top of their primordial kT): 2 (default, Figure 3) |k| = kt sin(theta),
//     theta ~ U[0, pi]; 0 a 2D Gaussian with per-component width kt; 1 fixed |k| = kt.
//   ptmaxmatch=<0|1|2>: optional SpaceShower:pTmaxMatch (empty keeps the Pythia default, 0; Pythia's
//     DIS example uses 2, which lets the shower fill the whole phase space).
//   primkt=<sigma>: optional per-component Gaussian width (GeV) used for all three of
//     BeamRemnants:primordialKTsoft/Hard/remnant; empty keeps the Pythia defaults.
//   mjoinjunction=<GeV>: optional FragmentationSystems:mJoinJunction, the excess-mass threshold
//     below which hadronization fuses two junction quarks into a diquark (status 74). Set with
//     forceParm, so values below the documented 0.5 GeV minimum are allowed; empty = 1.0 default.
//   struck=u|any: u (default) keeps the Figure 3 struck-u, SRC-eligible selection; any keeps
//     every struck quark/antiquark (no eligibility cut), for flavour-tagging studies.
//   struckid=<pdg>: with struck=any, keep only this struck-parton id (e.g. 2 = u, 1 = d),
//     valence and sea alike, vetoed at process level.
//   hadrons=<csv>: also write one row per final-state charged hadron of every event with a
//     current jet: event_id, charge, 3D opening angle to the jet axis in the Breit frame,
//     whether it is in the current hemisphere (Breit pz < 0), and whether it is a jet constituent.
//   qqc=<csv>: also write one row per final-state charged hadron of every selected event (jet
//     or not): event_id, charge, Breit azimuth and transverse momentum, current hemisphere flag
//     (Breit pz < 0), proton-forward lab eta, z_LC = q.P_h / q.P, SIDIS z_h = P.P_h / P.q, and
//     whether it is the leading (largest z_h) charged hadron of the current hemisphere.
//   utagonly=1: write event and qqc rows only for events passing the experimental struck-u tag
//     (leading current-hemisphere charged hadron is a pi+ or K+ with z_h > 0.5); nacc still
//     counts all selected DIS events.
// Besides the Figure 3 columns, each row has the struck-parton id and, whenever a current jet
// exists, its charge Q_kappa = sum_i (E_i/E_J)^kappa Q_i (Breit frame) and its leading
// (highest-energy) constituent and leading charged constituent. Independently of the jet, it
// also has the leading charged hadron of the Breit current hemisphere (pz < 0), i.e. the one
// with the largest SIDIS z_h = P.P_h / P.q, and that z_h.

#include "Pythia8/Pythia.h"
#include "Pythia8/FJcore.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>

using namespace Pythia8;

// Figure 3 DIS window, applied at process level so rejected events skip showers and
// hadronization. Same quantities as the original: Q^2 and x from the electron, struck u
// (or any struck parton with struck=any, for flavour-tagging studies).
class Fig3ProcessCuts : public UserHooks {
public:
  Fig3ProcessCuts(bool anyStruckIn, int struckIdIn) : anyStruck(anyStruckIn), struckId(struckIdIn) {}
  bool anyStruck;
  int struckId;
  bool canVetoProcessLevel() override {return true;}
  bool doVetoProcessLevel(Event& process) override {
    int iEin = -1, iEout = -1, iQin = -1;
    for (int i = 1; i < process.size(); ++i) {
      const Particle& p = process[i];
      if (p.status() == -21 && p.id() == 11) iEin = i;
      else if (p.status() == -21 && p.idAbs() <= 6) iQin = i;
      else if (p.status() == 23 && p.id() == 11) iEout = i;
    }
    if (iEin < 0 || iEout < 0 || iQin < 0) return true;
    if (!anyStruck && process[iQin].id() != 2) return true;
    if (struckId != 0 && process[iQin].id() != struckId) return true;
    Vec4 q = process[iEin].p() - process[iEout].p();
    double Q2 = -q.m2Calc();
    double x  = Q2 / (2. * (process[2].p() * q));
    return !(Q2 > 16. && Q2 < 25. && x > 0.05 && x <= 0.495);
  }
};

// Lab -> Breit: proton along +z, q = (0, 0, 0, -Q). Target hemisphere is pz > 0.
static RotBstMatrix breitFrame(const Vec4& pP, const Vec4& q) {
  RotBstMatrix M;
  M.toCMframe(pP, q);
  Vec4 qc = q;
  qc.rotbst(M);
  M.bst(0., 0., -qc.e() / qc.pz());
  return M;
}

static double angle2D(double ax, double ay, double bx, double by) {
  double ma = sqrt(ax * ax + ay * ay), mb = sqrt(bx * bx + by * by);
  if (ma < 1e-6 || mb < 1e-6) return NAN;
  double c = (ax * bx + ay * by) / (ma * mb);
  return acos(max(-1., min(1., c)));
}

int main(int argc, char* argv[]) {
  map<string, string> opt = {
    {"out", "fig3.csv"}, {"sample", "unchanged"}, {"src", "0"}, {"showers", "0"},
    {"remnant", "junction"}, {"nacc", "1000"}, {"seed", "12345"}, {"eidoffset", "0"},
    {"jetR", "1.0"}, {"kt", "0.4"}, {"ktdist", "2"}, {"ptmaxmatch", ""}, {"primkt", ""}, {"mjoinjunction", ""}, {"struck", "u"},
    {"struckid", "0"}, {"hadrons", ""}, {"qqc", ""}, {"utagonly", "0"}, {"xmldoc", XMLDIR}};
  for (int i = 1; i < argc; ++i) {
    string a = argv[i];
    size_t eq = a.find('=');
    if (eq == string::npos || !opt.count(a.substr(0, eq))) {
      fprintf(stderr, "bad argument: %s\n", a.c_str());
      return 1;
    }
    opt[a.substr(0, eq)] = a.substr(eq + 1);
  }
  bool srcOn    = stoi(opt["src"]) == 1;
  bool showers  = stoi(opt["showers"]) == 1;
  bool junction = opt["remnant"] == "junction";
  long nAcc     = stol(opt["nacc"]);
  long eidOff   = stol(opt["eidoffset"]);
  double jetR   = stod(opt["jetR"]);
  bool anyStruck = opt["struck"] == "any";
  if (!anyStruck && opt["struck"] != "u") {
    fprintf(stderr, "struck must be u or any\n");
    return 1;
  }
  int struckIdCut = stoi(opt["struckid"]);
  if (struckIdCut != 0 && !anyStruck) {
    fprintf(stderr, "struckid needs struck=any\n");
    return 1;
  }
  if (srcOn && !junction) {
    fprintf(stderr, "src=1 needs remnant=junction (the partner d must be a separate parton)\n");
    return 1;
  }

  Pythia pythia(opt["xmldoc"], false);

  // SRC settings read by the patched BeamRemnants (registered by the user program).
  pythia.settings.addFlag("SRC:on", false);
  pythia.settings.addFlag("SRC:splitDiquarkDIS", false);
  pythia.settings.addMode("SRC:partnerId", 1, true, true, 1, 2);
  pythia.settings.addMode("SRC:kTdist", 0, true, true, 0, 2);
  pythia.settings.addMode("SRC:kickMode", 0, true, true, 0, 1);
  pythia.settings.addParm("SRC:kT", 1.0, true, false, 0., 0.);
  pythia.settings.addParm("SRC:kbarWidth", 0.2, true, false, 0., 0.);

  vector<string> cards = {
    // Figure 3 beams and process.
    "Beams:frameType = 2", "Beams:idA = 11", "Beams:idB = 2212",
    "Beams:eA = 5.", "Beams:eB = 41.",
    "WeakBosonExchange:ff2ff(t:gmZ) = on", "PhaseSpace:Q2Min = 4.", "PDF:lepton = off",
    "PartonLevel:MPI = off",
    // The only switch that differs between the shower-off and shower-on samples.
    showers ? "PartonLevel:ISR = on" : "PartonLevel:ISR = off",
    showers ? "PartonLevel:FSR = on" : "PartonLevel:FSR = off",
    // DIS-appropriate ISR recoil (keeps the scattered lepton fixed); QCD radiation only.
    "SpaceShower:dipoleRecoil = on",
    "TimeShower:QEDshowerByL = off", "TimeShower:QEDshowerByQ = off",
    "TimeShower:QEDshowerByGamma = off",
    "SpaceShower:QEDshowerByL = off", "SpaceShower:QEDshowerByQ = off",
    // Primordial kT at the Pythia defaults, as in the original generation, unless primkt is set.
    "BeamRemnants:primordialKT = on",
    // SRC kick: back-to-back, added on top of the primordial kT of the struck u and the remnant d,
    // no pair kbar. Default ktdist = 2 is Figure 3 surgery v2, |k| = kt * sin(theta), theta ~ U[0, pi].
    srcOn ? "SRC:on = on" : "SRC:on = off",
    "SRC:splitDiquarkDIS = off", "SRC:partnerId = 1", "SRC:kTdist = " + opt["ktdist"], "SRC:kickMode = 1",
    "SRC:kT = " + opt["kt"], "SRC:kbarWidth = 0.",
    "Next:numberCount = 0", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0",
    "Next:numberShowProcess = 0",
    "Random:setSeed = on", "Random:seed = " + opt["seed"]};
  if (!opt["ptmaxmatch"].empty()) cards.push_back("SpaceShower:pTmaxMatch = " + opt["ptmaxmatch"]);
  if (!opt["primkt"].empty()) {
    // Same per-component width for initiators (soft and hard limits) and remnant partons.
    cards.push_back("BeamRemnants:primordialKTsoft = " + opt["primkt"]);
    cards.push_back("BeamRemnants:primordialKThard = " + opt["primkt"]);
    cards.push_back("BeamRemnants:primordialKTremnant = " + opt["primkt"]);
  }
  if (junction) {
    // Stock way to keep the remnant u and d separate (junction instead of ud diquark).
    // Colour reconnection is never applied to DIS remnants, so reconnect = off is cosmetic.
    cards.push_back("BeamRemnants:remnantMode = 1");
    cards.push_back("ColourReconnection:mode = 1");
    cards.push_back("ColourReconnection:reconnect = off");
    cards.push_back("BeamRemnants:beamJunction = on");
  }
  for (const string& c : cards)
    if (!pythia.readString(c)) { fprintf(stderr, "bad card: %s\n", c.c_str()); return 1; }
  if (!opt["mjoinjunction"].empty())
    pythia.settings.forceParm("FragmentationSystems:mJoinJunction", stod(opt["mjoinjunction"]));
  pythia.setUserHooksPtr(make_shared<Fig3ProcessCuts>(anyStruck, struckIdCut));
  if (!pythia.init()) return 1;

  FILE* out = fopen(opt["out"].c_str(), "w");
  if (!out) { fprintf(stderr, "cannot open %s\n", opt["out"].c_str()); return 1; }
  fprintf(out, "event_id,sample,pion,pdg_id,Q2,xB,ok,reason,zLC,phi_hJ,phi_hJ_parton,"
    "jet_pjt_GeV,parton_pjt_GeV,pion_pht_GeV,pion_lab_eta,jet_E_breit,jet_y_breit,"
    "pion_y_breit,jet_nconst,remnant_code,kt_scat_GeV,kt_d_GeV,remnant_merged,"
    "struck_id,tagjet_E,tagjet_nconst,jetQ_k0,jetQ_k0p3,jetQ_k0p5,jetQ_k1,"
    "lead_id,lead_z,leadch_id,leadch_z,leadcur_id,leadcur_z\n");

  FILE* hadOut = nullptr;
  if (!opt["hadrons"].empty()) {
    hadOut = fopen(opt["hadrons"].c_str(), "w");
    if (!hadOut) { fprintf(stderr, "cannot open %s\n", opt["hadrons"].c_str()); return 1; }
    fprintf(hadOut, "event_id,charge,theta_hJ,current,in_jet\n");
  }
  FILE* qqcOut = nullptr;
  if (!opt["qqc"].empty()) {
    qqcOut = fopen(opt["qqc"].c_str(), "w");
    if (!qqcOut) { fprintf(stderr, "cannot open %s\n", opt["qqc"].c_str()); return 1; }
    fprintf(qqcOut, "event_id,charge,phi,pt,current,eta_eic,zlc,zh,is_lead\n");
  }
  bool uTagOnly = stoi(opt["utagonly"]) == 1;

  auto num = [](double v) {
    if (std::isnan(v)) return string();
    char buf[32];
    snprintf(buf, sizeof(buf), "%.8g", v);
    return string(buf);
  };
  auto writeRow = [&](const vector<string>& fields) {
    for (size_t k = 0; k < fields.size(); ++k) {
      if (k) fputc(',', out);
      fputs(fields[k].c_str(), out);
    }
    fputc('\n', out);
  };

  fjcore::JetDefinition jetDef(fjcore::ee_genkt_algorithm, jetR, -1.);
  map<string, long> count;
  long nRows = 0, nNextFail = 0;
  auto t0 = std::chrono::steady_clock::now();
  const Event& ev = pythia.event;

  while (nRows < nAcc) {
    if (!pythia.next()) { ++nNextFail; continue; }
    ++count["accepted_process_window"];

    // Event-level DIS kinematics from the electron (re-check of the process-level window).
    // Scattered electron = final copy of the hard outgoing lepton (status -23 -> 44 recoil copy
    // when showers are on); "highest-energy e-" is not safe at 5 x 41 GeV, where Dalitz
    // electrons from hard target-side pi0/eta can win.
    int iLep = -1;
    for (int i = 0; i < ev.size(); ++i)
      if (abs(ev[i].status()) == 23 && ev[i].id() == 11) { iLep = ev[i].iBotCopyId(); break; }
    if (iLep < 0 || !ev[iLep].isFinal()) { ++count["no_scattered_lepton"]; continue; }
    Vec4 pP = ev[2].p(), q = ev[1].p() - ev[iLep].p();
    double Q2 = -q.m2Calc(), xB = Q2 / (2. * (pP * q));
    if (!(Q2 > 16. && Q2 < 25. && xB > 0.05 && xB <= 0.495)) {
      ++count["event_window_mismatch"];
      continue;
    }

    // Remnant composition code: diquark id, or 10*n_d + n_u for separate quarks.
    int nRem = 0, nRemD = 0, nRemU = 0, idDiq = 0;
    for (int i = 0; i < ev.size(); ++i) if (abs(ev[i].status()) == 63) {
      ++nRem;
      if (ev[i].id() == 1) ++nRemD;
      else if (ev[i].id() == 2) ++nRemU;
      else if (ev[i].idAbs() > 1000 && ev[i].idAbs() < 10000) idDiq = ev[i].id();
    }
    int remnantCode = idDiq ? idDiq : 10 * nRemD + nRemU;

    // SRC eligibility, identical in both arms: valence-u initiator (the patch's condition).
    // The default DIS remnant marks a sea initiator as valence after merging its companion
    // into a hadron, so there the remnant itself is used: exactly one ud diquark (the
    // original surgery required ud0 = 2101; remnant_code keeps 2101 vs 2103 for filtering).
    bool eligible = anyStruck || (pythia.info.id2() == 2 && (junction
      ? pythia.beamB[0].isValence() && pythia.beamB[0].id() == 2
      : nRem == 1 && (idDiq == 2101 || idDiq == 2103)));
    if (!eligible) { ++count["not_valence_u_initiator"]; continue; }
    int struckId = pythia.info.id2();
    ++count["struck_" + to_string(struckId)];

    // Diagnostics only (no effect on the selection). Realized transverse momenta in Pythia's
    // own hadronic rest frame (remnant vs scattered system before the kT boost, as in
    // BeamRemnants::setOneRemnKinematics): the scattered system (status 62 copies, i.e. the
    // primordial kT plus SRC kick given to the initiator) and the remnant d. Also whether
    // hadronization recombined the separate remnant d into a diquark (status 74).
    Vec4 pScat, pScatPre;
    int iRemD = -1;
    for (int i = 0; i < ev.size(); ++i) {
      if (abs(ev[i].status()) == 62) { pScat += ev[i].p(); pScatPre += ev[ev[i].mother1()].p(); }
      if (abs(ev[i].status()) == 63 && ev[i].id() == 1 && iRemD < 0) iRemD = i;
    }
    double ktScat = NAN, ktD = NAN;
    if (pScatPre.e() > 0.) {
      RotBstMatrix toHad;
      toHad.toCMframe(pScatPre, ev[1].p() + ev[2].p() - ev[iLep].p() - pScatPre);
      Vec4 a = pScat;
      a.rotbst(toHad);
      ktScat = a.pT();
      if (iRemD >= 0) {
        Vec4 b = ev[iRemD].p();
        b.rotbst(toHad);
        ktD = b.pT();
      }
    }
    string merged;
    if (iRemD >= 0) {
      int iDau = ev[ev[iRemD].iBotCopyId()].daughter1();
      bool isMerged = iDau > 0 && abs(ev[iDau].status()) == 74 && ev[iDau].idAbs() > 1000;
      merged = isMerged ? "1" : "0";
      ++count[isMerged ? "remnant_d_merged" : "remnant_d_separate"];
    }
    long eventId = eidOff + nRows;
    ++nRows;

    RotBstMatrix toBreit = breitFrame(pP, q);
    double qDotP = q * pP;

    // Leading charged pion: proton-forward lab eta in (2.0, 3.5), Breit target hemisphere,
    // highest Breit energy. Pythia's lab has the electron along +z, so eta_EIC = -eta.
    int iPi = -1, leadCurId = 0, iLeadCur = -1;
    double eBest = -1., leadCurZ = NAN;
    vector<fjcore::PseudoJet> inputs;
    for (int i = 0; i < ev.size(); ++i) {
      const Particle& p = ev[i];
      if (!p.isFinal() || i == iLep || !p.isVisible()) continue;
      Vec4 pb = p.p();
      pb.rotbst(toBreit);
      inputs.emplace_back(pb.px(), pb.py(), pb.pz(), pb.e());
      inputs.back().set_user_index(i);
      // Leading charged hadron of the current hemisphere by SIDIS z_h = P.P_h / P.q.
      if (p.isHadron() && p.isCharged() && pb.pz() < 0.) {
        double zh = (pP * p.p()) / qDotP;
        if (std::isnan(leadCurZ) || zh > leadCurZ) { leadCurZ = zh; leadCurId = p.id(); iLeadCur = i; }
      }
      if (p.idAbs() != 211) continue;
      double etaEIC = -p.eta();
      if (!(etaEIC > 2.0 && etaEIC < 3.5)) continue;
      if (pb.pz() <= 0.) continue;
      if (pb.e() > eBest) { eBest = pb.e(); iPi = i; }
    }

    // Experimental struck-u tag: leading current charged hadron is a pi+ or K+ with z_h > 0.5.
    bool uTagged = (leadCurId == 211 || leadCurId == 321) && leadCurZ > 0.5;
    if (uTagged) ++count["utagged"];
    bool writeEvent = !uTagOnly || uTagged;

    // Charged final-state hadrons for the azimuthal charge-charge correlation.
    if (qqcOut && writeEvent) {
      for (int i = 0; i < ev.size(); ++i) {
        const Particle& p = ev[i];
        if (!p.isFinal() || !p.isHadron() || !p.isCharged()) continue;
        Vec4 pb = p.p();
        pb.rotbst(toBreit);
        fprintf(qqcOut, "%ld,%d,%.5f,%.5f,%d,%.4f,%.5f,%.5f,%d\n", eventId, int(round(p.charge())),
          pb.phi(), pb.pT(), pb.pz() < 0. ? 1 : 0, -p.eta(), (q * p.p()) / qDotP,
          (pP * p.p()) / qDotP, i == iLeadCur ? 1 : 0);
      }
    }

    // Struck-parton "jet" of the original analysis: final copy of the hard outgoing quark.
    Vec4 kOut;
    for (int i = 0; i < ev.size(); ++i)
      if (abs(ev[i].status()) == 23 && ev[i].idAbs() <= 6) {
        kOut = ev[ev[i].iBotCopyId()].p();
        break;
      }
    kOut.rotbst(toBreit);

    // Current-hemisphere jet: highest-energy e+e- anti-kT jet with Breit pz < 0.
    fjcore::ClusterSequence cs(inputs, jetDef);
    vector<fjcore::PseudoJet> jets = cs.inclusive_jets(0.);
    int iJet = -1;
    for (int j = 0; j < int(jets.size()); ++j)
      if (jets[j].pz() < 0. && (iJet < 0 || jets[j].E() > jets[iJet].E())) iJet = j;

    // Flavour-tagging quantities of the current jet (filled whenever the jet exists).
    const double kappas[4] = {0., 0.3, 0.5, 1.};
    double jetQ[4] = {NAN, NAN, NAN, NAN}, tagE = NAN, leadZ = NAN, leadChZ = NAN;
    int tagN = -1, leadId = 0, leadChId = 0;
    if (iJet >= 0) {
      tagE = jets[iJet].E();
      vector<fjcore::PseudoJet> cons = jets[iJet].constituents();
      tagN = int(cons.size());
      for (int k = 0; k < 4; ++k) jetQ[k] = 0.;
      double eLead = -1., eLeadCh = -1.;
      for (const fjcore::PseudoJet& c : cons) {
        const Particle& p = ev[c.user_index()];
        double z = c.E() / tagE;
        for (int k = 0; k < 4; ++k) jetQ[k] += pow(z, kappas[k]) * p.charge();
        if (c.E() > eLead) { eLead = c.E(); leadId = p.id(); leadZ = z; }
        if (p.isCharged() && c.E() > eLeadCh) { eLeadCh = c.E(); leadChId = p.id(); leadChZ = z; }
      }
      // Charged final-state hadrons: opening angle to the jet axis (Breit frame).
      if (hadOut) {
        vector<bool> inJet(ev.size(), false);
        for (const fjcore::PseudoJet& c : cons) inJet[c.user_index()] = true;
        const fjcore::PseudoJet& jet = jets[iJet];
        double jAbs = sqrt(jet.px() * jet.px() + jet.py() * jet.py() + jet.pz() * jet.pz());
        for (int i = 0; i < ev.size(); ++i) {
          const Particle& p = ev[i];
          if (!p.isFinal() || !p.isHadron() || !p.isCharged()) continue;
          Vec4 pb = p.p();
          pb.rotbst(toBreit);
          double cosT = (pb.px() * jet.px() + pb.py() * jet.py() + pb.pz() * jet.pz())
            / max(1e-12, pb.pAbs() * jAbs);
          fprintf(hadOut, "%ld,%d,%.6g,%d,%d\n", eventId, int(round(p.charge())),
            acos(max(-1., min(1., cosT))), pb.pz() < 0. ? 1 : 0, inJet[i] ? 1 : 0);
        }
      }
    }

    // Same order of checks as the original: pion candidate, then z_LC; then the jet.
    string pion, reason = "ok";
    int pdg = 0, nConst = -1;
    double zLC = NAN, phiJet = NAN, phiParton = NAN, jetPt = NAN, pionPt = NAN,
      etaEIC = NAN, jetE = NAN, jetY = NAN, pionY = NAN;
    if (iPi < 0) {
      reason = "no_pion_candidate";
    } else {
      Vec4 pib = ev[iPi].p();
      pib.rotbst(toBreit);
      pdg = ev[iPi].id();
      pion = pdg > 0 ? "piplus" : "piminus";
      zLC = (q * ev[iPi].p()) / qDotP;
      etaEIC = -ev[iPi].eta();
      pionPt = pib.pT();
      pionY = pib.rap();
      phiParton = angle2D(kOut.px(), kOut.py(), pib.px(), pib.py());
      if (zLC < 0.1 || zLC > 1. + 1e-6) reason = "zLC_range";
      else if (iJet < 0) reason = "no_jet";
      else {
        const fjcore::PseudoJet& jet = jets[iJet];
        phiJet = angle2D(jet.px(), jet.py(), pib.px(), pib.py());
        jetPt = jet.pt();
        jetE = jet.E();
        jetY = jet.rap();
        nConst = int(jet.constituents().size());
      }
    }
    bool ok = reason == "ok";
    ++count[ok ? "ok_" + pion : reason];
    if (writeEvent) writeRow({to_string(eventId), opt["sample"], pion, to_string(pdg), num(Q2), num(xB),
      ok ? "True" : "False", reason, num(zLC), num(phiJet), num(phiParton), num(jetPt),
      num(kOut.pT()), num(pionPt), num(etaEIC), num(jetE), num(jetY), num(pionY),
      nConst < 0 ? "" : to_string(nConst), to_string(remnantCode), num(ktScat), num(ktD), merged,
      to_string(struckId), num(tagE), tagN < 0 ? "" : to_string(tagN), num(jetQ[0]), num(jetQ[1]),
      num(jetQ[2]), num(jetQ[3]), to_string(leadId), num(leadZ), to_string(leadChId), num(leadChZ),
      to_string(leadCurId), num(leadCurZ)});
  }
  fclose(out);
  if (hadOut) fclose(hadOut);
  if (qqcOut) fclose(qqcOut);

  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  FILE* js = fopen((opt["out"] + ".summary.json").c_str(), "w");
  fprintf(js, "{\n  \"pythia_version\": %.3f,\n", pythia.parm("Pythia:versionNumber"));
  for (const auto& kv : opt)
    if (kv.first != "xmldoc") fprintf(js, "  \"%s\": \"%s\",\n", kv.first.c_str(), kv.second.c_str());
  for (const auto& kv : count) fprintf(js, "  \"n_%s\": %ld,\n", kv.first.c_str(), kv.second);
  fprintf(js, "  \"n_rows\": %ld,\n  \"n_next_failed\": %ld,\n", nRows, nNextFail);
  fprintf(js, "  \"n_tried\": %ld,\n", pythia.info.nTried());
  fprintf(js, "  \"seconds\": %.1f\n}\n", secs);
  fclose(js);
  return 0;
}
