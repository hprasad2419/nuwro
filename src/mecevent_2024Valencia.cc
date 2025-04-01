#include "mecevent_2024Valencia.h"

void mecevent_2024Valencia (params & p, event & e, nucleus & t, bool cc)
{
  e.par = p;            // save params in the event
  e.flag.cc = cc;       // set flags for the event
  e.flag.nc = !cc;
  e.flag.dis = false;
  e.flag.qel = false;
  e.flag.coh = false;
  e.flag.mec = true;
  int ile_pb = p.mec_pb_trials;
  double mec_central = p.mec_central_motion;
  double mec_smearing = p.mec_back_to_back_smearing;
  double binding = p.kaskada_w;
  int mecskalowanie = p.mec_scaling;

  double mc_sampling[3] = {0,0,0};          //sampling direction for pp, np, and pn respectively
  int mc_strength[3] = {1,1,1};      //strength for sampling direction pp, np, and pn respectively
  mc_sampling[0] = p.MEC_cm_direction_pp;
  mc_strength[0] = p.MEC_cm_strength_pp;

  mc_sampling[1] = p.MEC_cm_direction_np;
  mc_strength[1] = p.MEC_cm_strength_np;
  
  mc_sampling[2] = p.MEC_cm_direction_pn;
  mc_strength[2] = p.MEC_cm_strength_pn;
 

  // sadly, only CC events available so far...
  if(e.flag.nc)
  {
    cerr<<" MEC error: Wrong Settings!\n";
    e.weight = 0;
    return;
  }

  particle meclepton;
  particle meclepton_3p3h;
  ap=(e.in[0].pdg<0);


  meclepton_3p3h.pdg = meclepton.pdg = e.in[0].pdg-1+2*ap;
  meclepton.set_mass (PDG::mass (meclepton.pdg)); //set mass coresponding to pdg
  meclepton_3p3h.set_mass (PDG::mass(meclepton_3p3h.pdg));
  ml=meclepton.mass();
  ml2=ml*ml;

 
  PB=p.MEC_pauli_blocking;

  // Binding energy / Correlation within the medium 

  if (p.nucleus_p > 12) {
    Bmec = ap ? E_corr[5] : E_corr[4];
  } else if ( (p.nucleus_p <=12) and (p.nucleus_p > 6) ) {
      Bmec = ap ? E_corr[3] : E_corr[2];
    } else {
        Bmec = ap ? E_corr[1] : E_corr[0];
      } 


  double q0max = e.in[0].energy() - ml - Bmec;
  width_q0 = q0max;

  if( q0max > qmax_Nieves)
  {
    q0max = qmax_Nieves;
  }

  
  particle mecnucleon_in[2];
  particle mecnucleon_out[2];

  particle mecnucleon_3p3h_in[3];
  particle mecnucleon_3p3h_out[3];

  double individual_weight[4] = {0,0,0,0};
  double weight_2p2h=0;
  double weight_3p3h=0;

  if(width_q0>0)
  {
    weight_2p2h=Valencia2020_kin_and_weight_2p2h (e.in[0].E(), individual_weight, meclepton, mecnucleon_in, mecnucleon_out, flag_2p2h_pn, t, mec_central, mec_smearing, binding, ile_pb, mc_sampling, mc_strength);
    
    //weight_3p3h = Valencia2020_kin_and_weight_3p3h (e.in[0].E(), individual_weight, meclepton_3p3h, mecnucleon_3p3h_in, mecnucleon_3p3h_out, t, binding, ile_pb);
  }

  double weight = weight_2p2h + weight_3p3h;
  e.weight = weight;
  if(weight > 0)
  {
    double ratio_2p2h  = weight_2p2h/weight;
    flag_2p2h = (bool)(frandom() < ratio_2p2h);
    flag_3p3h = (!flag_2p2h);
  
    if(flag_2p2h)
    {
        e.in.push_back (mecnucleon_in[0]);
        e.in.push_back (mecnucleon_in[1]);
        e.out.push_back (meclepton);
        e.out.push_back (mecnucleon_out[0]);
        e.out.push_back (mecnucleon_out[1]);
    }
    else if (flag_3p3h) {
        e.in.push_back (mecnucleon_3p3h_in[0]);
        e.in.push_back (mecnucleon_3p3h_in[1]);
        e.in.push_back (mecnucleon_3p3h_in[2]);
        e.out.push_back (meclepton_3p3h);
        e.out.push_back (mecnucleon_3p3h_out[0]);
        e.out.push_back (mecnucleon_3p3h_out[1]);
        e.out.push_back (mecnucleon_3p3h_out[2]);
    }
    else {
      cout << "Impossible event topology !! \n Stopping .. \n";
      std::exit(1);
    }
  }
  // Reset the event topology flags
  flag_2p2h = false;
  flag_2p2h_pn = false;
  flag_3p3h = false;
  
}
