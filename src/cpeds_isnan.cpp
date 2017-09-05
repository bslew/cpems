//****************************************************************************************************************

struct LongDouble {
  unsigned char mant[8];
  unsigned int exp;
};

union Real {
  struct LongDouble sld;
  long double ld;
};

int CPEDS_isnan (double x) {
  union Real real;
  real.ld = x;
  if((real.sld.exp &0x7FFF) == 0x7FFF)
    return 1;

  return 0;

}

//****************************************************************************************************************
