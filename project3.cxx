// 6860 Final Project
// Relative Resonance of Lambda-0
// Dana Ali, Alec Ferensic, Jamar Philip, Akash Raj

// C++ Headers
#include <iostream>
#include <fstream>
#include <random>
#include <vector>

// ROOT Headers
#include <TRandom.h>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;

//resonance properties
const double M_RES = 1.115683; // Resonance mass in GeV/c^2
const double WIDTH_RES = 0.0249; // Resonance width in GeV/c^2
const double LIFETIME_RES = 2.632e-10; // Resonance lifetime in seconds

//decay products and their properties
const double M_PROTON = 0.938272; // Proton mass in GeV/c^2
const double M_PION = 0.139570; // Pion mass in GeV/c^2

const double SQRT_SNN = 200.0; // Center-of-mass energy in GeV

const int N_EVENTS = 10000; // Number of events
const int N_BG = 10; // Number of background particles per event

class FourMomentum {
   private:
      TLorentzVector vec;
      
   public:
      FourMomentum(double px, double py, double pz, double e) : vec(px, py, pz, e) {}
      FourMomentum(TLorentzVector v) : vec(v) {}

      double Px() const { return vec.Px(); }
      double Py() const { return vec.Py(); }
      double Pz() const { return vec.Pz(); }
      double E() const { return vec.E(); }

      FourMomentum Boost(const FourMomentum& resonance) {
         TLorentzVector v1 (resonance.Px(), resonance.Py(), resonance.Pz(), resonance.E());
         TVector3 v2 = v1.BoostVector();
         TLorentzVector myparticle_newsys=vec;
	     myparticle_newsys.Boost(v2);
	    //cout << vec.M() << endl;
         return FourMomentum(myparticle_newsys);
      }
     
      
      // Setters
    void SetPx(double px) { vec.SetPx(px); }
    void SetPy(double py) { vec.SetPy(py); }
    void SetPz(double pz) { vec.SetPz(pz); }
    void SetE(double e) { vec.SetE(e); }
    
    //functions
    double Pt() const {return sqrt(vec.Px()*(vec.Px()) + vec.Py()*(vec.Py())); }
    double P() const { return sqrt(vec.Px()*(vec.Px()) + vec.Py()*(vec.Py()) + vec.Pz()*(vec.Pz())); }
    double M() const { return sqrt(vec.E()*(vec.E()) - vec.P()*(vec.P())); }
	
	//operators
	FourMomentum operator + (FourMomentum &);

      void Print() const {
         cout << "(" << Px() << ", " << Py() << ", " << Pz() << ", " << E() << ")" << endl;
      }
};

FourMomentum FourMomentum::operator+(FourMomentum &m)
{
	return FourMomentum(vec.Px()+m.Px(),vec.Py()+m.Py(),vec.Pz()+m.Pz(),vec.E()+m.E());
}


// Function to generate a random number according to a Breit-Wigner distribution

double breit_wigner(double mean, double gamma) {
    gRandom->SetSeed (0);
    double bw = (float)gRandom->BreitWigner(mean, gamma);
    return bw;
}

// Function to perform the resonance decay
    vector<FourMomentum> ResonanceDecay(const FourMomentum& resonance) {
    vector<FourMomentum> products;

    // Get the momentum of the resonance in its rest frame
    
    //FourMomentum resonance_lab = resonance.Boost(beta_res, gamma_res);
    double P_res = resonance.P();
    double E_res = resonance.E();
    double M = resonance.M();
    double beta_res = P_res / E_res;
    double gamma_res = E_res / M;
    
    //to boost to lab frame using TLorentzVector    

    // Generate random decay angles in the rest frame
    double theta = acos(2.0*rand()/RAND_MAX - 1.0);
    double phi = 2.0*3.14159*rand()/RAND_MAX;

    // Calculate the momentum of the daughter particles in the rest frame
    double P_d = sqrt((M*M - M_PROTON*M_PROTON - M_PION*M_PION)*(M*M - M_PROTON*M_PROTON - M_PION*M_PION)
            - 4.0*M_PROTON*M_PROTON*M_PION*M_PION)/(2.0*M_RES);
    double E_d_PROTON = sqrt(P_d*P_d + M_PROTON*M_PROTON);
    double E_d_PION = sqrt(P_d*P_d + M_PION*M_PION);
    //double beta_d = P_d / E_d_;
    //double gamma_d = E_d / M_PROTON;

    //cout<<P_d<<" "<<E_d_PROTON<<endl;
    // Calculate the 4-momenta of the daughter particles in the rest frame
    FourMomentum proton(P_d*sin(theta)*cos(phi), P_d*sin(theta)*sin(phi), P_d*cos(theta), E_d_PROTON);
    FourMomentum pion(-P_d*sin(theta)*cos(phi), -P_d*sin(theta)*sin(phi), -P_d*cos(theta), E_d_PION);

    //cout<< pion.Pz()<<" " <<endl;
    // Boost the daughter particles to the lab frame
    FourMomentum proton_lab = proton.Boost(resonance);
    FourMomentum pion_lab = pion.Boost(resonance);
   //cout<<proton_lab.Pz()<<endl;
    
    products.push_back(proton_lab);
    products.push_back(pion_lab);

    return products;
}

int project() {
    // Initialize histograms for the invariant mass and background
    TFile *f=new TFile("lorenz.root","RECREATE");
    TH1F* h_mass = new TH1F("h_mass", "Invariant Mass", 100, 1.08, 1.15);
    TH1F* h_bg = new TH1F("h_bg", "Background", 100, 1.08, 1.15);

    // Initialize the random number generator
    TRandom3 rng;

    // Loop over events
    for (int i_event = 0; i_event < N_EVENTS; i_event++) {
        // Generate a resonance particle
        double m_res = breit_wigner(M_RES, WIDTH_RES);
        FourMomentum resonance(0,0, SQRT_SNN/(2),sqrt((SQRT_SNN*SQRT_SNN/4)+m_res*m_res));

        // Decay the resonance particle
       // std::vector
    FourMomentum parent = resonance;
    vector<FourMomentum> daughters = ResonanceDecay(parent);

	
	
    // Add Gaussian smearing to the daughter particle momenta (data is too sensitive to handle, so didn't do it)
    for (auto& daughter : daughters) {
        //daughter.SetPx(rng.Gaus(daughter.Px(), daughter.Px()*0.00005));
        //daughter.SetPy(rng.Gaus(daughter.Py(), daughter.Py()*0.00005));
        //daughter.SetPz(rng.Gaus(daughter.Pz(), daughter.Pz()*0.00005));
        //daughter.SetE(abs(rng.Gaus(daughter.E(), daughter.E()*0.00005)));
    }
   
  

    // Calculate the invariant mass of the daughter particles
    double inv_mass = (daughters[0] + daughters[1]).M();
	//cout<<"result "<<inv_mass<<endl;
    // Fill the invariant mass histogram
    //h_mass->Fill(inv_mass);
    h_mass->Fill(rng.Gaus(inv_mass, inv_mass*0.005));

    // Generate uncorrelated background events
    int j=0;
    while(j<N_BG){
    	double bg_inv_mass = rng.Uniform(h_mass->GetXaxis()->GetXmin(), h_mass->GetXaxis()->GetXmax());
    	h_bg->Fill(bg_inv_mass);
    	j++;
    //h_mass->Fill(bg_inv_mass);
	}
	}
	// histograms

	//TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
	h_mass->SetLineColor(kRed);
	h_mass->GetXaxis()->SetTitle("Invariant Mass [GeV/c^{2}]");
	h_mass->GetYaxis()->SetTitle("Counts");
	h_mass->Add(h_bg);
	h_mass->Draw();
	
	h_bg->SetLineColor(kBlue);
	h_bg->Draw("SAME");
	//f->Write();
	//f->Close();
	// Print the number of signal and background events in the signal region
	double signal = h_mass->Integral(h_mass->FindBin(M_RES - 3.0*WIDTH_RES), h_mass->FindBin(M_RES + 3.0*WIDTH_RES));
	double bg = h_bg->Integral(h_bg->FindBin(M_RES - 3.0*WIDTH_RES), h_bg->FindBin(M_RES + 3.0*WIDTH_RES));
	cout << "Number of signal events: " << signal << endl;
	cout << "Number of background events: " << bg << endl;
	
	return 0;
}

/*
This code should generate events for the decay of a Lambda particle resonance, apply Gaussian smearing to the daughter particle momenta, and calculate the invariant mass of the daughter particles. It also includes a Breit-Wigner distribution for generating the resonance mass, uncorrelated background events, and histograms for the invariant mass and background. Finally, it calculates the number of signal and background events in the signal region and prints them to the console.
*/

