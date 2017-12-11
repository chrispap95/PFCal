#ifndef _hgcssmiphit_hh_
#define _hgcssmiphit_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <iostream>


class HGCSSMipHit{
  
public:
  
  HGCSSMipHit():
    e_(0),
    l_(0),
    x_(0),
    y_(0),
    z_(0),
    nfrac_(0){
    for (unsigned i(0); i<9; ++i){
      neigh_e_samelayer_[i] = 0;
      neigh_e_prevlayer_[i] = 0;
      neigh_e_nextlayer_[i] = 0;
      neigh_e_prev2layer_[i] = 0;
      neigh_e_next2layer_[i] = 0;
    }
    for (unsigned i(0); i<5; ++i){
      xref_neighlay_[i] = 0;
      yref_neighlay_[i] = 0;
    }
  };
  ~HGCSSMipHit(){
    
  };

  inline double e(){
    return e_;
  };
  inline double x(){
    return x_;
  };
  inline double y(){
    return y_;
  };
  inline double z(){
    return z_;
  };
  inline double noiseFrac(){
    return nfrac_;
  };
  inline unsigned l(){
    return l_;
  };

  inline void setE(const double & val){
    e_ = val;
  };
  inline void setx(const double & val){
    x_ = val;
  };
  inline void sety(const double & val){
    y_ = val;
  };
  inline void setz(const double & val){
    z_ = val;
  };
  inline void setnoiseFrac(const double & val){
    nfrac_ = val;
  };
  inline void setLayer(const unsigned val){
    l_ = val;
  };

  void set_xref_neighlay(const unsigned idx,const double & val);
  void set_yref_neighlay(const unsigned idx,const double & val);

  void set_neigh_e_samelayer(const unsigned idx,const double & e);
  void set_neigh_e_prevlayer(const unsigned idx,const double & e);
  void set_neigh_e_nextlayer(const unsigned idx,const double & e);
  void set_neigh_e_prev2layer(const unsigned idx,const double & e);
  void set_neigh_e_next2layer(const unsigned idx,const double & e);

  inline double xref_neighlay(const unsigned idx){
    if (idx>4) return 0;
    return xref_neighlay_[idx];
  };

  inline double yref_neighlay(const unsigned idx){
    if (idx>4) return 0;
    return yref_neighlay_[idx];
  };

  inline double neigh_e_samelayer(const unsigned idx){
    if (idx>8) return 0;
    return neigh_e_samelayer_[idx];
  };
  inline double neigh_e_prevlayer(const unsigned idx){
    if (idx>8) return 0;
    return neigh_e_prevlayer_[idx];
  };
  inline double neigh_e_nextlayer(const unsigned idx){
    if (idx>8) return 0;
    return neigh_e_nextlayer_[idx];
  };
  inline double neigh_e_prev2layer(const unsigned idx){
    if (idx>8) return 0;
    return neigh_e_prev2layer_[idx];
  };
  inline double neigh_e_next2layer(const unsigned idx){
    if (idx>8) return 0;
    return neigh_e_next2layer_[idx];
  };
  
  double getMaxEnergy(int idx);
  double getSumEnergy(int idx);

 private:
  double e_;
  unsigned l_;
  double x_;
  double y_;
  double z_;
  double nfrac_;

  //layer-energy of closest neighbours
  double neigh_e_samelayer_[9];
  double neigh_e_prevlayer_[9];
  double neigh_e_nextlayer_[9];
  double neigh_e_prev2layer_[9];
  double neigh_e_next2layer_[9];

  double xref_neighlay_[5];
  double yref_neighlay_[5];

  ClassDef(HGCSSMipHit,4);


};//class

typedef std::vector<HGCSSMipHit> HGCSSMipHitVec;

#endif
