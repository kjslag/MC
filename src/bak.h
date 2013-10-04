

// Value

struct Value
{
  Value() : _count(0) {}
  
  void operator<<(double x)
  {
    uint k = 0;
    next_bin:
    
    if ( k < _bins.size() ) {
      Bin &bin = _bins[k];
      
      bin.x2 += x*x;
      if ( !full(k) )
        bin.x = x;
      else {
        x += bin.x;
        ++k;
        goto next_bin;
      }
    }
    else
      _bins.push_back(x);
    
    ++_count;
  }
  
  double sum() const
  {
    double sum = 0;
    const uint n = _bins.size();
    FOR(k,n)
      if ( full(k) )
        sum += _bins[k].x;
    return sum;
  }
  
  uint64_t count() const { return _count; }
  
private:
  struct Bin {
    Bin(double y) : x(y), x2(y*y) {}
    
    double x, x2;
  };
  
  bool full(uint k) const
  { return _count & (uint64_t(1) << k); }
  
  vector<Bin> _bins;
  uint64_t    _count;
};



template<uint n>
struct ExternalFieldPotential_ : public SpinFunc_<n>
{
  typedef Spin_<n> Spin;
  
  ExternalFieldPotential_(const Spin &h_) : h(h_) {}
  
  virtual Float operator()(Spin s) const final { return -(s|h); }
  
  Spin h;
};


template<uint n>
using SpinMatrix_ = array<Spin_<n>,n>;

template<uint n>
Spin_<n> operator|(const SpinMatrix_<n> &M, const Spin_<n> s0) {
  Spin_<n> s = s0;
  FOR(k,n) s[k] = M[k] | s0;
  return s;
}

template<uint n>
struct RotateSpin_ : public SpinFlipper_<n> {
  typedef Spin_<n>       Spin;
  typedef SpinMatrix_<n> SpinMatrix;
  
  RotateSpin_(const vector<SpinMatrix> &matrices) : Rs(matrices), R(Rs.front()), _dist(0, Rs.size()-1) { }
  
  virtual void  reset() final { R = Rs[_dist(random_engine)]; }
  virtual void  flip(Spin &s) const final { s = R|s; }
  virtual Float delta(Spin s1, Spin s2) const final { return s1|(R1|s2); }
  
  const vector<SpinMatrix> Rs;
  
protected:
  void set_R1() {
    R1 = R;
    FOR(k,n) R[k][k] = R[k][k] - 1;
  }
  
private:
  SpinMatrix &R;
  SpinMatrix  R1;
  std::uniform_int_distribution<uint> _dist;
};
