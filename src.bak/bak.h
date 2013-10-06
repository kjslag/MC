



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

  virtual void fast_update() final __attribute__((hot))
  {
    FOR(d, dim)
      Assert( L[d]%2 == 0, L[d] );
    
    if ( _update_sublattice == 2 )
      _update_sublattice = bool_dist(random_engine);
    else
      _update_sublattice = !_update_sublattice;
    
    for (Index i0=0; i0<N; i0+=N/L[dim-1]) {
      const Pos p0 = pos(i0);
      Index i = _update_sublattice;
      FOR(d, dim-1)
        i += p0[d]%2; // check this
      
      for (i = i%2; i<i0+L[dim-1]; i+=2) {
        //const Spin s1 = _spins[i];
        Spin sn(Spin::zero);
        
        for(Index nn : nearestNeighbors(i))
          sn += _spins[nn];
        
        //TODO
      }
    }
    
    ++_nSweeps;
  }
  
/*
run: MC
        @echo > jobs
        @all=""; \
           for method in smart; \
        do for L in 20; \
        do for J in `seq  0 .25 +3`; \
        do for u in `seq -3 .25 +3`; \
        do      f="results/vison_hexagon_sVBS_$${L}_$${J}_$${u}_$${method}"; \
                echo -e "$$f:\n\t./MC --L $$L $$L $$L --J $$J --J $$u --potential 'vison hexagon s-VBS' --sweeps 1000 --update-method $$method --file $$f\n" >> jobs; \
                all="$$all $$f"; \
        done; done; done; done; \
        echo "all:$$all" >> jobs
        @mkdir -p results
        make -j4 -f jobs all
*/

