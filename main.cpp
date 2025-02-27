#include <iostream>
#include <complex>
#include <cmath>
#include <fftw3.h>

class Transformata
{
  private:
    Transformata(const Transformata& s);//zabezpieczenie konstruktora (konstr. kopiujacym) i nizej operatora przypisania
    Transformata& operator=(const Transformata& s);
    std::complex<double> *fi, *FFTfi, *C, *IDFTfi, *IFFTfi;
    std::complex<double> **fi2d;
    int licznikDFT, licznikFFT, wymiar, N, M;
    double *DFTmodules, *FFTmodules;
  public:
    Transformata();
    ~Transformata();
    void Set();
    void Zero(std::complex<double> *tab, int r);
    void Build1D();
    void Build2D();
    void DFT();
    void IDFT();
    void FFT(std::complex<double>* f, int r);
    void IFFT(std::complex<double>* f, int r);
    void DFT2D();
    double CountModule(double re, double im){return sqrt(pow(re, 2) + pow(im, 2));}
    void Get1D();
    void Start();
};

Transformata::Transformata()
{
  licznikDFT = 0;
  licznikFFT = 0;
  wymiar = 0;
  N = 0;
  M = 0;
}

void Transformata::Set()
{
  std::cin >> wymiar;
  if(wymiar == 1) std::cin >> N;
  else std::cin >> N >> M;
}

void Transformata::Zero(std::complex<double> *tab, int r)
{
  for(int i=0;i<r;i++) tab[i] = 0.0;
}

void Transformata::Build1D()
{
  fi = new std::complex<double>[N]; FFTfi = new std::complex<double>[N]; C = new std::complex<double>[N];
  IDFTfi = new std::complex<double>[N]; IFFTfi = new std::complex<double>[N];
  DFTmodules = new double[N]; FFTmodules = new double[N];
  for(int i=0;i<N;i++)
  {
    std::cin >> fi[i];
    FFTfi[i] = fi[i];
  }
}

void Transformata::Build2D()
{
  fi2d = new std::complex<double>*[N];
  for(int i=0;i<N;i++){fi2d[i] = new std::complex<double>[M];}
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<M;j++) std::cin >> fi2d[i][j];
  }
}

void Transformata::DFT()
{
  Zero(C, N);
  std::complex<double> j(0.0, 1.0);
  for(int n=0;n<N;n++)
  {
    for(int i=0;i<N;i++)
    {
      C[n] += fi[i] * exp((-j * std::complex<double>((2*M_PI*i*n), 0.0))/static_cast<double>(N));
      licznikDFT++;
    }
    if(fabs(C[n].real()) < 5e-4) C[n].real(0.0);
    if(fabs(C[n].imag()) < 5e-4) C[n].imag(0.0);
    DFTmodules[n] = CountModule(C[n].real(), C[n].imag());
  }
}

void Transformata::IDFT()
{
  std::complex<double> j(0.0, 1.0);
  for(int i=0;i<N;i++)
  {
    IDFTfi[i] = std::complex<double>(0.0, 0.0);
    for(int n=0;n<N;n++)
    {
      IDFTfi[i] += C[n] * exp((j * std::complex<double>((2*M_PI*i*n), 0.0))/static_cast<double>(N));
    }
    IDFTfi[i] *= (1/static_cast<double>(N));
  }
}

void Transformata::FFT(std::complex<double>* f, int r)
{
  if(r <= 1) return;

  std::complex<double>* even = new std::complex<double>[r / 2];
  std::complex<double>* odd = new std::complex<double>[r / 2];

  for (int i=0;i<r/2;++i)
  {
    even[i] = f[i*2];
    odd[i] = f[i*2+1];
  }

  FFT(even, r/2);
  FFT(odd, r/2);

  for(int i=0;i<r/2;++i)
  {
    std::complex<double> t = std::polar(1.0, -2.0 * M_PI * i / r) * odd[i];
    f[i] = even[i] + t;
    f[i+r/2] = even[i] - t;
    licznikFFT+=2;
  }

  delete[]even;
  delete[]odd;
}

void Transformata::IFFT(std::complex<double>* f, int r)
{
  if(r <= 1) return;

  std::complex<double>* even = new std::complex<double>[r / 2];
  std::complex<double>* odd = new std::complex<double>[r / 2];

  for (int i=0; i<r/2; ++i)
  {
    even[i] = f[i*2];
    odd[i] = f[i*2 + 1];
  }

  IFFT(even, r/2);
  IFFT(odd, r/2);

  for(int i=0; i<r/2; ++i)
  {
    std::complex<double> t = std::polar(1.0, 2.0 * M_PI * i / r) * odd[i];
    f[i] = even[i] + t;
    f[i + r/2] = even[i] - t;
  }

  delete[]even;
  delete[]odd;
}

void Transformata::DFT2D()
{
  bool REAL = 0, IMAG = 1;
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * M);
  fftw_plan plan = fftw_plan_dft_2d(N, M, (fftw_complex*)fi2d, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<M;j++)
    {
      //std::cout << out[i*M+j][REAL] << " " << out[i*M+j][IMAG] << " ";
      std::cout << CountModule(out[i*M+j][REAL], out[i*M+j][IMAG]) << " ";
    }
    std::cout << std::endl;
  }

  fftw_destroy_plan(plan);
  fftw_free(out);
}

void Transformata::Get1D()
{
  std::cout << "Zlozonosc DFT: " << licznikDFT << std::endl;
  std::cout << "Zlozonosc FFT: " << licznikFFT << std::endl;
  for(int i=0;i<N;i++)
  {
    if(fabs(FFTfi[i].real()) < 0.0005) FFTfi[i].real(0.0);
    if(fabs(FFTfi[i].imag()) < 0.0005) FFTfi[i].imag(0.0);
    FFTmodules[i] = CountModule(FFTfi[i].real(), FFTfi[i].imag());
  }
  for(int i=0;i<N;i++) std::cout << i << " " << fi[i].real() << " " << DFTmodules[i] << " " << IDFTfi[i].real() << " " << FFTmodules[i] << " " << IFFTfi[i].real()/static_cast<double>(N) << std::endl;
}

void Transformata::Start()
{
  Set();
  if(wymiar == 1)
  {
    Build1D();
    DFT();
    IDFT();
    FFT(FFTfi, N);
    for(int i=0;i<N;i++) IFFTfi[i] = FFTfi[i];
    IFFT(IFFTfi, N);
    Get1D();
  }
  else
  {
    Build2D();
    DFT2D();
  }
}

Transformata::~Transformata()
{
  if(fi) delete[]fi;
  if(FFTfi) delete[]FFTfi;
  if(IDFTfi) delete[]IDFTfi;
  if(IFFTfi) delete[]IFFTfi;
  if(C) delete[]C;
  if(DFTmodules) delete[]DFTmodules;
  if(FFTmodules) delete[]FFTmodules;
  if(fi2d)
  {
    for(int i=0;i<N;i++) delete[]fi2d[i];
    delete[]fi2d;
  }
}

int main()
{
  Transformata T;
  T.Start();
  return 0;
}