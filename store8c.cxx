#include <iostream>
#include <TRandom3.h>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TCanvas.h> 
#include "TLorentzVector.h"
#include "TImage.h"

//constants
const double M_RES = 1.11568; // Resonance mass in (GeV)
const double WIDTH_RES = 0.0249; // Resonance width in (Gev)
const double M_PROTON = 0.939; // Proton mass in (GeV)
const double M_PION = 0.140; // Pion mass (GeV)
const double SQRT_SNN = 200.0; // Center-of-mass energy in GeV
const int N_SIGNAL_EVENTS = 10000; //number of signal events
const int N_BACKGROUND_EVENTS = 7000; // number of background events

//function to generate random number according to Breit-Wigner distribution
double breit_wigner(double mean, double gamma) {
    TRandom3 rng;
    rng.SetSeed(0);
    return rng.BreitWigner(mean, gamma);
}
class Event {
private:
    TLorentzVector resonance;
    TLorentzVector proton;
    TLorentzVector pion;

public:
    Event() {}

    //function to generate a resonance particle
    //void generateResonance() {
      //  double m_res = breit_wigner(M_RES, WIDTH_RES);
        //resonance.SetPxPyPzE(0, 0, SQRT_SNN / 2.0, sqrt(SQRT_SNN * SQRT_SNN / 4.0 + m_res * m_res));
    //}

void generateResonance() {
        double m_res = breit_wigner(M_RES, WIDTH_RES);
        //resonance.SetPxPyPzE(0, 0, SQRT_SNN / 2.0, sqrt(SQRT_SNN * SQRT_SNN / 4.0 + m_res * m_res));
        //double m_res = M_RES; // fixed Res mass
        //double m_res = breit_wigner(M_RES, WIDTH_RES); // BW noise Res mass
        double fixed_pT=1.0;
        double eta1 = 1.0 * (gRandom->Rndm() - 0.5);
        //calculate theta from eta
        //double theta = 2 * atan(exp(-eta));
        //generate random azimuthal angle phi between 0 and 2*pi (isotropic)
        double phi1 = 2.0 * TMath::Pi() * gRandom->Rndm();
        
        //resonance.SetPtEtaPhiE(pt, eta1, phi1, sqrt(SQRT_SNN * SQRT_SNN / 4.0 + m_res * m_res));

        double px_resonance = fixed_pT * cos(phi1);
        double py_resonance = fixed_pT * sin(phi1);
        double pz_resonance = fixed_pT * sinh(eta1);
        double energy_resonance = sqrt(px_resonance * px_resonance + py_resonance * py_resonance + pz_resonance * pz_resonance + m_res*m_res);
       resonance.SetPxPyPzE(px_resonance, py_resonance, pz_resonance, energy_resonance);
    }

    //to perform resonance decay
    void decayResonance() {

double eta = 1.0 * (gRandom->Rndm() - 0.5);
 // Calculate theta from eta
  double theta = 2.0 * atan(exp(-eta));
//generate random azimuthal angle phi between 0 and 2*pi (isotropic)
double phi = 2.0 * TMath::Pi() * gRandom->Rndm();

//calculate momenta of daughter particles in resonance rest frame
double p = sqrt((M_RES * M_RES - M_PROTON * M_PROTON - M_PION * M_PION) *
                (M_RES * M_RES - M_PROTON * M_PROTON - M_PION * M_PION) -
                4.0 * M_PROTON * M_PROTON * M_PION * M_PION) /
           (2.0 * M_RES);
double E_proton = sqrt(p * p + M_PROTON * M_PROTON);
double E_pion = sqrt(p * p + M_PION * M_PION);

//create TLorentzVectors for daughter particles
proton.SetPxPyPzE(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), E_proton);
pion.SetPxPyPzE(-p * sin(theta) * cos(phi), -p * sin(theta) * sin(phi), -p * cos(theta), E_pion);


    }

    //to boost daughter particles to lab frame
    void boostToLabFrame() {
        proton.Boost(resonance.BoostVector());
        pion.Boost(resonance.BoostVector());
    }

    // Setter function to set momentum and energy components for the proton
    void setProtonMomentum(double px, double py, double pz, double energy) {
        proton.SetPxPyPzE(px, py, pz, energy);
    }

    // Setter function to set momentum and energy components for the pion
    void setPionMomentum(double px, double py, double pz, double energy) {
        pion.SetPxPyPzE(px, py, pz, energy);
    }

    //to get  the invariant mass of daughter particles
    double getInvariantMass() const {
        TLorentzVector total = proton + pion;
        return total.M();
    }

    //to get eta values for the decay products
    //double getEta() const {
      //  double eta_proton = proton.Eta();
        //double eta_pion = pion.Eta();
        //return (eta_proton + eta_pion) / 2.0;
    //}

    //to get momentum components and energy of  the proton
    double getProtonPx() const { return proton.Px(); }
    double getProtonPy() const { return proton.Py(); }
    double getProtonPz() const { return proton.Pz(); }
    double getProtonE() const { return proton.E(); }
    double getProtonPhi() const { return proton.Phi(); }
    double getProtonEta() const { return proton.Eta(); }
    //to get momentum components and energy of the pion
    double getPionPx() const { return pion.Px(); }
    double getPionPy() const { return pion.Py(); }
    double getPionPz() const { return pion.Pz(); }
    double getPionE() const { return pion.E(); }
    double getPionPhi() const { return pion.Phi(); }
    double getPionEta() const { return pion.Eta(); }
};

int main() {
    //create 2 histograms for resonance events + background events
    TH1F* h_signal_mass = new TH1F("h_signal_mass", "Invariant mass distribution", 200, 1.05, 1.2); 
    TH1F* h_background_mass = new TH1F("h_background_mass", "Background", 200, 1.05, 1.2);

    //generate background events
    for (int i = 0; i < N_BACKGROUND_EVENTS; ++i) {
      // generate resonance particle -method1
        //double bg_mass = breit_wigner(1.04, 0.1);
       //double inv_mass = m_res;
       //double smeared_background = gRandom->Gaus(inv_mass, inv_mass * 0.004);
        // Add background invariant mass to histogram
//method2
        double bg_mass;
        //double bg_mass2
        //bg_mass2 = gRandom->Uniform(1.08, 1.15);
//method3
        bg_mass = gRandom->Gaus(1.09, 0.009);
        double smeared_bg_mass = gRandom->Gaus(bg_mass, bg_mass * 0.005);
        h_background_mass->Fill(smeared_bg_mass);

    }


    // create a TTree to store signal event data
    TFile* file = new TFile("store82.root", "RECREATE");
    TTree* tree = new TTree("Events", "Events");

    //define variables to store the event data
    double inv_mass;
    double proton_px, proton_py, proton_pz, proton_e, proton_phi, proton_eta;
    double pion_px, pion_py, pion_pz, pion_e, pion_phi, pion_eta;

    //set branch addresses for the TTree
    tree->Branch("InvariantMass", &inv_mass);
    tree->Branch("ProtonPx", &proton_px);
    tree->Branch("ProtonPy", &proton_py);
    tree->Branch("ProtonPz", &proton_pz);
    tree->Branch("ProtonE", &proton_e);
    tree->Branch("ProtonPhi", &proton_phi);
    tree->Branch("ProtonEta", &proton_eta);
    tree->Branch("PionPx", &pion_px);
    tree->Branch("PionPy", &pion_py);
    tree->Branch("PionPz", &pion_pz);
    tree->Branch("PionE", &pion_e);
    tree->Branch("PionPhi", &pion_phi);
    tree->Branch("PionEta", &pion_eta);

    // Loop for resonance events generation
    for (int i = 0; i < N_SIGNAL_EVENTS; ++i) {
        // Generate resonance particle
        double m_res = breit_wigner(M_RES, WIDTH_RES);
        TLorentzVector resonance(0, 0, SQRT_SNN / 2.0, sqrt(SQRT_SNN * SQRT_SNN / 4.0 + m_res * m_res)); // Declare resonance here

        //perform resonance decay
        double eta = 1.0 * (gRandom->Rndm() - 0.5);
        double theta = 2 * atan(exp(-eta));
        double phi = 2.0 * TMath::Pi() * gRandom->Rndm(); // Use TMath::Pi() for Pi constant
        double p = sqrt((M_RES * M_RES - M_PROTON * M_PROTON - M_PION * M_PION) *
                        (M_RES * M_RES - M_PROTON * M_PROTON - M_PION * M_PION) -
                        4.0 * M_PROTON * M_PROTON * M_PION * M_PION) /
                   (2.0 * M_RES);
        double E_proton = sqrt(p * p + M_PROTON * M_PROTON);
        double E_pion = sqrt(p * p + M_PION * M_PION);
        TLorentzVector proton(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), E_proton); 
        TLorentzVector pion(-p * sin(theta) * cos(phi), -p * sin(theta) * sin(phi), -p * cos(theta), E_pion); 

        //boost daughter particles to lab frame
        proton.Boost(resonance.BoostVector());
        pion.Boost(resonance.BoostVector());

        //calculate invariant mass of daughter particles
        inv_mass = (proton + pion).M();

        //store kinematic variables of daughters
        proton_px = proton.Px();
        proton_py = proton.Py();
        proton_pz = proton.Pz();
        proton_e = proton.E();
        proton_phi = proton.Phi();
        proton_eta = proton.Eta();
        pion_px = pion.Px();
        pion_py = pion.Py();
        pion_pz = pion.Pz();
        pion_e = pion.E();
        pion_phi = pion.Phi();
        pion_eta = pion.Eta();

        //fill event data into the TTree
        tree->Fill();
        double smeared_inv_mass = gRandom->Gaus(inv_mass, inv_mass * 0.005);
        // add resonance invariant mass to histogram
        h_signal_mass->Fill(smeared_inv_mass);
    }

    // create canvas for plotting
    TCanvas* canvas = new TCanvas("canvas", "Invariant Mass", 800, 600);

    // histograms with both signal and background
    h_signal_mass->Draw();
    h_background_mass->Draw("same"); 

    h_signal_mass->GetXaxis()->SetTitle("Invariant Mass (GeV)");
    h_signal_mass->GetYaxis()->SetTitle("Events");

    h_signal_mass->SetLineColor(kRed);
    h_signal_mass->SetFillColorAlpha(kRed, 0.5); 
    h_background_mass->SetLineColor(kBlue); 
    h_background_mass->SetFillColorAlpha(kBlue, 0.5); 

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // legend coordinates
    legend->AddEntry(h_signal_mass, "Lambda resonance", "l"); 
    legend->AddEntry(h_background_mass, "Background", "l"); 
    legend->Draw(); 

    // save canvas plot to png
    canvas->SaveAs("InvariantMassPlot10.png");

    // write histograms and TTree to root output file
    h_signal_mass->Write();
    h_background_mass->Write();
    tree->Write();

    file->Close();

    return 0;
}
