#include <iostream>
#include <complex>
#include <cmath>
#include <fftw3.h>

void Zero(std::complex<double> *C, int N_in)
{
  for(int i=0;i<N_in;i++) C[i] = std::complex<double>(0.0, 0.0);
}

double CountModule(double Re, double Im)
{
  return sqrt(pow(Re, 2) + pow(Im, 2));
}

void DFT(std::complex<double> *C, double *f, int N_in, double *Mods)
{
  int licznik = 0;
  std::complex<double> j(0.0, 1.0);
  for(int n=0;n<N_in;n++)
  {
    for(int i=0;i<N_in;i++)
    {
      C[n] += std::complex<double>(f[i], 0.0) * exp((-j * std::complex<double>((2*M_PI*i*n), 0.0))/static_cast<double>(N_in));
      licznik++;
    }
    if(fabs(C[n].real()) < 5e-4) C[n].real(0.0);
    if(fabs(C[n].imag()) < 5e-4) C[n].imag(0.0);
    Mods[n] = CountModule(C[n].real(), C[n].imag());
  }
  std::cout << "Licznik DFT: " << licznik << std::endl;
}

void IDFT(std::complex<double> *f, std::complex<double> *C, int N_in)
{
  std::complex<double> j(0.0, 1.0);
  for(int i=0;i<N_in;i++)
  {
    for(int n=0;n<N_in;n++)
    {
      f[i] += C[n] * exp((j * std::complex<double>((2*M_PI*i*n), 0.0))/static_cast<double>(N_in));
    }
    f[i] *= (1/static_cast<double>(N_in));
  }
}

void FFT(std::complex<double> *C, double *f, int N_in, double *Mods)
{
  int licznik = 0;
  std::complex<double> j(0.0, 1.0);
  for(int n=0;n<(N_in/2);n++)
  {
    for(int i=0;i<(N_in/2);i++)
    {
      C[n] += std::complex<double>(f[2*i], 0.0) * exp(-j * std::complex<double>((2*M_PI*i*n), 0.0)/(static_cast<double>(N_in)/2)) + exp(-j * std::complex<double>((2*M_PI*n), 0.0)/static_cast<double>(N_in)) * std::complex<double>(f[2*i+1], 0.0) * exp(-j * std::complex<double>((2*M_PI*i*n),0.0)/(static_cast<double>(N_in)/2));
      C[n + (N_in/2)] += std::complex<double>(f[2*i], 0.0) * exp(-j * std::complex<double>((2*M_PI*i*n), 0.0)/(static_cast<double>(N_in)/2)) - exp(-j * std::complex<double>((2*M_PI*n), 0.0)/static_cast<double>(N_in)) * std::complex<double>(f[2*i+1], 0.0) * exp(-j * std::complex<double>((2*M_PI*i*n),0.0)/(static_cast<double>(N_in)/2));
      licznik++;
    }
  }
  for(int n=0;n<N_in;n++)
  {
    if(fabs(C[n].real()) < 5e-4) C[n].real(0.0);
    if(fabs(C[n].imag()) < 5e-4) C[n].imag(0.0);
    Mods[n] = CountModule(C[n].real(), C[n].imag());
  }
  std::cout << "Licznik FFT: " << licznik << std::endl;
}

void IFFT(std::complex<double> *f, std::complex<double> *C, int N_in)
{
  std::complex<double> j(0.0, 1.0);
  for(int i=0;i<N_in/2;i++)
  {
    for(int n=0;n<N_in/2;n++)
    {
      f[i] += (1/static_cast<double>(N_in)) * ( C[2*n] * exp(j * std::complex<double>((2*M_PI*i*n), 0.0)/(static_cast<double>(N_in)/2)) + exp(j * std::complex<double>((2*M_PI*n), 0.0)/static_cast<double>(N_in)) * C[2*n+1] * exp(j * std::complex<double>((2*M_PI*i*n), 0.0)/(static_cast<double>(N_in)/2)) );
      f[i+(N_in/2)] += (1/static_cast<double>(N_in)) * ( C[2*n]* exp(j * std::complex<double>((2*M_PI*i*n), 0.0)/(static_cast<double>(N_in)/2)) - exp(j * std::complex<double>((2 * M_PI * n),0.0)/static_cast<double>(N_in)) * C[2*n+1] * exp(j * std::complex<double>((2*M_PI*i*n),0.0)/(static_cast<double>(N_in)/2)) );
    }
  }
}

int main()
{
  int wymiar = 0;
  int N = 0, M = 0;
  double *fi;
  double *DFTmodules;
  double *FFTmodules;
  std::complex<double> *Cn, *IDFTfi, *IFFTfi;
  double **fi2d;
  
  std::cin >> wymiar;
  if(wymiar == 1)
  {
    std::cin >> N;
    fi = new double[N];
    Cn = new std::complex<double>[N], IDFTfi = new std::complex<double>[N], IFFTfi = new std::complex<double>[N];
    DFTmodules = new double[N];
    FFTmodules = new double[N];
    for(int i=0;i<N;i++) std::cin >> fi[i];
    DFT(Cn, fi, N, DFTmodules);
    IDFT(IDFTfi, Cn, N);
    Zero(Cn, N);
    FFT(Cn, fi, N, FFTmodules);
    IFFT(IFFTfi, Cn, N);
    for(int i=0;i<N;i++) std::cout << i << " " << DFTmodules[i] << " " << IDFTfi[i].real() << " " << FFTmodules[i] << " " << IFFTfi[i].real() << std::endl;
  }
  else
  {
    std::cin >> N >> M;
    fi2d = new double*[N];
    for(int i=0;i<N;i++) fi2d[i] = new double[M];
    for(int i=0;i<N;i++)
    {
      for(int j=0;j<M;j++) std::cin >> fi2d[i][j];
    }

    // Allocate memory for FFTW plans and output data
    fftw_complex* in = reinterpret_cast<fftw_complex*>(fi2d);
    fftw_complex* out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N * M));

    // Create FFTW plans
    fftw_plan plan = fftw_plan_dft_2d(N, M, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute the FFT
    fftw_execute(plan);

    // Output the result as complex numbers and their modulus
    /*for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << std::fixed << "(" << out[i * N + j][0] << " + " << out[i * N + j][1] << "i)   ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;*/
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double modulus = sqrt(pow(out[i * N + j][0], 2) + pow(out[i * N + j][1], 2));
            std::cout << std::fixed << modulus << "   ";
        }
        std::cout << std::endl;
    }

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(out);
  }
  
  delete[]fi;
  delete[]IDFTfi;
  delete[]IFFTfi;
  delete[]Cn;
  delete[]DFTmodules;
  delete[]FFTmodules;
  for(int i=0;i<M;i++) delete[]fi2d[i];
  delete[]fi2d;
  return 0;
}