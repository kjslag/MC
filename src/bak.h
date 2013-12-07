virtual Spin ideal_spin() const final
  {
    const Float c=cos(pi/8), s=sin(pi/8), r=sqrt(.5);
    Spin candidates[] = {
      {{1,0,0,0,0,0}}, // c-VBS
      {{c,s,0,0,0,0}}, // plaquette VBS
      {{r,r,0,0,1,0}}, // s-VBS
      {{r,r,0,0,0,1}}, // caterpillar VBS
      {{c,s,c,s,c,s}}, // star VBS
      {{s,c,s,c,s,c}}  // triangle VBS
    };
    Spin *min_s = nullptr;
    Float min_V = infinity<Float>();
    for (Spin &s0: candidates) {
      s0.normalize();
      const Float V0 = _this(s0);
      if ( V0 < min_V ) {
        min_s = &s0;
        min_V = V0;
      }
    }
    Assert (min_s);
    return *min_s;
  }


    if ( dynamic_cast<const cos_Potential_<3>*>(_V) )
      _sum_cos3.add( cos_potential<3>(sum) );
    if ( dynamic_cast<const OnZ2Z2_Potential_<n/2>*>(_V) ) {
      pair<LongFloat,LongFloat> V_uv = OnZ2Z2_potential(sum);
      _sum_OnZ2Z2_u     .add( V_uv.first  - sum_OnZ2Z2_u_sub()*sum_4 );
      _sum_OnZ2Z2_v     .add( V_uv.second - sum_OnZ2Z2_v_sub()*sum_4 );
      _sum_OnZ2Z2_sigma2.add( V_uv.first );
      _sum_OnZ2Z2_sigma4 << sq(V_uv.first);
    }
    if ( dynamic_cast<const VisonTrianglePotential*>(_V) ) {
      pair<LongFloat,LongFloat> V_uv = visonTriangle_potential(sum);
      _sum_tri_u.add( V_uv.first - sum_tri_u_sub()*sum_4 );
      _sum_tri_v.add( V_uv.second );
    }

template<uint n>
struct ExternalFieldPotential_ : public SpinFunc_<n>
{
  typedef Spin_<n> Spin;
  
  ExternalFieldPotential_(const Spin &h_) : h(h_) {}
  
  virtual Float operator()(Spin s) const final { return -(s|h); }
  
  Spin h;
};


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


#include <execinfo.h> // backtrace_symbols_fd
static void print_backtrace()
{
  //g_on_error_query("MC");
  constexpr int size = 20;
  void *buffer[size];
  backtrace_symbols_fd(buffer, backtrace(buffer, size), STDOUT_FILENO);
}

static void single_handler(int sig)
{
  printf("Error: signal %d:\n", sig);
  
  print_backtrace();
  
  exit(1);
}

  // http://www.gnu.org/s/hello/manual/libc/Signal-Handling.html
  signal( SIGSEGV, single_handler );
  //signal( SIGABRT, single_handler );
  