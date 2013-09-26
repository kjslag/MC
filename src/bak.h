
struct Value
{
  Value() : _count(0) {}
  
  void operator<<(double x)
  {
    uint k = 0;
    next_bin:
    
    if ( k < _bins.size() )
    {
      Bin &bin = _bins[k];
      
      bin.x2 += x*x;
      if ( !full(k) )
        bin.x = x;
      else
      {
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
  struct Bin
  {
    Bin(double y) : x(y), x2(y*y) {}
    
    double x, x2;
  };
  
  bool full(uint k) const
  { return _count & (uint64_t(1) << k); }
  
  vector<Bin> _bins;
  uint64_t    _count;
};
