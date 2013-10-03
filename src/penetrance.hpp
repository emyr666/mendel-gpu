class Penetrance:public MendelGPU{
public:
  Penetrance();
  ~Penetrance();
  virtual void compute_penetrance()=0;
};

class PenetranceGlf:public Penetrance{
public:
  PenetranceGlf();
  ~PenetranceGlf();
  void compute_penetrance();
private:
  void precompute_penetrance();
  void impute_penetrance_matrix();
  void precompute_penetrance_fast(int * center_snp);
};

class PenetranceReads:public Penetrance{
public:
  PenetranceReads();
  ~PenetranceReads();
  void compute_penetrance();
};
