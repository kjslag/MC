
struct Accumulator
{
  Accumulator() : _n_samples(0) {}
  
  typedef function<void(LongFloat,uint)> Func;
  
  void set_func(const Func &f) { _f = f; }
  
  void operator<<(LongFloat x)
  {
    uint k = 0;
    next_bin: {
      if ( k < _sums.size() ) {
        if (_f) _f(x, k);
        Sum &sum = _sums[k];
        const bool full = _n_samples & exp2i(k);
        if ( !full )
          sum.x = x;
        else {
          x += sum.x;
          sum.x = 0;
          ++k;
          goto next_bin;
        }
      } else
        _sums.push_back(x);
    }
  }
  
  LongFloat operator[](uint k) const { return _sums[k].x; }
  LongFloat      first(uint k) const { return _sums[k].first; }
  
  uint64_t n_samples() const { return _n_samples; }
  uint size() const { return _sums.size(); }
  
  LongFloat total(uint k=0) const {
    LongFloat tot=0;
    const uint n = size();
    for(; k<n; ++k)
      tot += _sums[k].x;
    return tot;
  }
  
protected:
  struct Sum {
    Sum(LongFloat x_) : x(x_), first(x_) {}
    
    LongFloat x;
    const LongFloat first;
  };
  typedef vector<Sum> Sums;
  
  friend ostream& operator<<(ostream &os, const Sum &sum);
  friend ostream& operator<<(ostream &os, const Accumulator &acc);
  
private:
  uint64_t _n_samples;
  Sums     _sums;
  Func     _f;
};

ostream& operator<<(ostream &os, const Accumulator::Sum &sum)
{ return os << "sum[" << sum.x << ", " << sum.first << "]"; }

ostream& operator<<(ostream &os, const Accumulator &acc)
{ return os << "accumulator[" << acc._sums << "]"; } // todo fix [{}]

struct Value
{
  Value() : _means_bin_num(0) {
    _x_acc.set_func( [this](LongFloat x, uint k){this->x_acc_func(x, k);} );
    _means.reserve(max_means);
  }
  
  void operator<<(const LongFloat x)
  {
    Assert( _means_bin_num == calc_means_bin_num(),
            max_means, n_samples(), _means_bin_num, calc_means_bin_num() );
    
    _x_acc << x;
  }
  
  uint64_t n_samples() const { return _x_acc.n_samples(); }
  
  uint default_thermalization_exp() const
  { return max(int(_x_acc.size()) - 6, 0); } // TODO
  
  LongFloat mean() const
  {
    // TODO: calc from _binned_vec
    const uint n_means = _means.size();
    LongFloat mean_ = 0;
    if ( n_means >= 2 ) {
      for (uint i=1; i<n_means; ++i)
        mean_ += _means[i];
      mean_ /= n_means - 1;
    } else
      mean_ = NAN;
    
    return mean_;
  }
  
  const vector<LongFloat>& means() const { return _means; }
  
  LongFloat error() const { return n_samples()>1 ? error(default_thermalization_exp()) : INFINITY; }
  
  LongFloat error(uint bin_num) const
  {
    /*const LongFloat mean_ = mean(bin_num);
    const LongFloat x2 = _binned_vec[bin_num].x2 - _binned_vec[bin_num].first_x2;
    const uint64_t n_samples = error_bin_count(bin_num);
    return sqrt(max(LongFloat(0), x2 - n_samples*bin_size(bin_num)*mean_*mean_)) / n_samples;*/
  }
  
  vector<LongFloat> binned_errors() const
  {
    vector<LongFloat> errors_;
    const int n_bins = _x2_binned_acc.size();
    FOR(bin_num, n_bins-1)
      errors_.push_back( error(bin_num) );
    return errors_;
  }
  
protected:
  typedef vector<Accumulator> BinnedAccumulators;
  
  friend ostream& operator<<(ostream &os, const BinnedAccumulators &binned_accs);
  friend ostream& operator<<(ostream &os, const Value &value);
  
  void x_acc_func(LongFloat x, uint k)
  {
    if ( k == _means_bin_num ) {
      if ( _means.size() < max_means )
        _means.push_back(x/exp2i(_means_bin_num));
      else {
        FOR(i, max_means/2)
          _means[i] = (_means[2*i] + _means[2*i+1])/2;
        _means.resize(max_means/2);
        //Assert( !full(k), k, n_samples() ); // else there would be double counting
        ++_means_bin_num;
         // x will be added later
      }
    }
    
    if ( k == _x2_binned_acc.size() )
      _x2_binned_acc.emplace_back();
    
    _x2_binned_acc[k] << x*x;
  }
  
  // number of samples per bin
  static uint64_t bin_size(uint bin_num)
  { return exp2i(bin_num); }
  
  // number of samples that went into _binned_vec[bin_num].x2
  uint64_t error_bin_count(uint bin_num) const
  { return (n_samples() & ~(bin_size(bin_num)-1)) - bin_size(bin_num); }
  
  uint calc_means_bin_num() const {
    uint n = 0;
    while ( exp2i(n)*max_means < n_samples() )
      ++n;
    return n;
  }
  
private:
  Accumulator        _x_acc;
  BinnedAccumulators _x2_binned_acc;
  vector<LongFloat>  _means;
  uint               _means_bin_num;
  
  static const uint max_means;
};
const uint Value::max_means = 1u<<4; // should be a power of 2 // TODO change to 1<<8

ostream& operator<<(ostream &os, const Value::BinnedAccumulators &binned_accs)
{ return os << "binnedAccumulators[" << binned_accs << "]"; } // TODO: fix [{}]

ostream& operator<<(ostream &os, const Value &value)
{
  return os << "value["                  << value.mean()
            << ", "                      << value.error()
            << ", \"means\" -> "         << value.means()
            << ", \"binned errors\" -> " << value.binned_errors()
            << ", \"# samples\" -> "     << value.n_samples()
            << "]"; // TODO Accumulators
}
