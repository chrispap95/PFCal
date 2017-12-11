#include "HGCSSMipHit.hh"

void HGCSSMipHit::set_xref_neighlay(const unsigned idx,const double & val){
  if (idx>4) {
    std::cout << " -- Problem with neighbour index, too large!" << std::endl;
    exit(1);
  }
  xref_neighlay_[idx] = val;
}

void HGCSSMipHit::set_yref_neighlay(const unsigned idx,const double & val){
  if (idx>4) {
    std::cout << " -- Problem with neighbour index, too large!" << std::endl;
    exit(1);
  }
  yref_neighlay_[idx] = val;
}

void HGCSSMipHit::set_neigh_e_samelayer(const unsigned idx,const double & e){
    if (idx>8) {
      std::cout << " -- Problem with neighbour index, too large!" << std::endl;
      exit(1);
    }
    neigh_e_samelayer_[idx] = e;
}

void HGCSSMipHit::set_neigh_e_prevlayer(const unsigned idx,const double & e){
    if (idx>8) {
      std::cout << " -- Problem with neighbour index, too large!" << std::endl;
      exit(1);
    }
    neigh_e_prevlayer_[idx] = e;
}

void HGCSSMipHit::set_neigh_e_nextlayer(const unsigned idx,const double & e){
    if (idx>8) {
      std::cout << " -- Problem with neighbour index, too large!" << std::endl;
      exit(1);
    }
    neigh_e_nextlayer_[idx] = e;
}

void HGCSSMipHit::set_neigh_e_prev2layer(const unsigned idx,const double & e){
    if (idx>8) {
      std::cout << " -- Problem with neighbour index, too large!" << std::endl;
      exit(1);
    }
    neigh_e_prev2layer_[idx] = e;
}

void HGCSSMipHit::set_neigh_e_next2layer(const unsigned idx,const double & e){
    if (idx>8) {
      std::cout << " -- Problem with neighbour index, too large!" << std::endl;
      exit(1);
    }
    neigh_e_next2layer_[idx] = e;
}

double HGCSSMipHit::getMaxEnergy(int idx){
  if (abs(idx)>2) {
    std::cout << " Problem, max allowed is 2... exiting..." << std::endl;
    exit(1);
  }
  double maxE = 0;
  double *e = idx==-1?neigh_e_prevlayer_ : idx==1?neigh_e_nextlayer_ :
    idx==-2 ? neigh_e_prev2layer_ : idx==2? neigh_e_next2layer_ : neigh_e_samelayer_;
  if (!e) return 0;
  for (unsigned i(0); i<9;++i){
    if (idx==0 && i==0) continue;
    //if (idx==-1 && e[i]!=0) std::cout << i << " " << e[i] << " " << neigh_e_prevlayer_[i] << std::endl;
    if (e[i]>maxE) maxE=e[i];
  }
  //if (idx==-1 && maxE!=0) std::cout << " max : " << maxE << std::endl;
  return maxE;
}

double HGCSSMipHit::getSumEnergy(int idx){
  if (abs(idx)>2) {
    std::cout << " Problem, max allowed is 2... exiting..." << std::endl;
    exit(1);
  }
  double sumE = 0;
  double *e = idx==-1?neigh_e_prevlayer_ : idx==1?neigh_e_nextlayer_ :
    idx==-2 ? neigh_e_prev2layer_ : idx==2? neigh_e_next2layer_ : neigh_e_samelayer_;
  if (!e) return 0;
  for (unsigned i(0); i<9;++i){
    if (idx==0 && i==0) continue;
    //if (idx==-1 && e[i]!=0) std::cout << i << " " << e[i] << " " << neigh_e_prevlayer_[i] << std::endl;
    sumE+=e[i];
  }
  //if (idx==-1 && maxE!=0) std::cout << " max : " << maxE << std::endl;
  return sumE;
}
