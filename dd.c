/* Time-stamp: <Fri Apr 27 07:59:37 JST 2018> 

   1D DRIFT-DIFFUSION SIMULATOR / Nobuya Mori (Osaka Univ)

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define P2(x) ((x) * (x))
#define P3(x) ((x) * (x) * (x))
#define P4(x) ((x) * (x) * (x) * (x))
#define P5(x) ((x) * (x) * (x) * (x) * (x))

/* physical constants */
// 物理定数

#define PV 8.854187817620e-12 // electric constant (F/m)電気定数
#define EC 1.6021766208e-19   // elementary charge (C)電荷
#define KB 1.38064852e-23     // Boltzmann constant (J/K)ボルツマン定数

/* material parameters and conditions */
// 材料パラメータ、条件

#define TL 300.0          // lattice temperature (K)格子温度
#define VT (KB * TL / EC) // thermal voltage (V)熱電圧
#define NI 1.08e16        // intrinsic density (/m3)固有密度
#define DC (11.9 * PV)    // dielectric constant (F/m)誘電率
#define MUP 0.050         // hole mobility (m2/Vs)ホール移動度
#define MUN 0.147         // electron mobility (m2/Vs)電子移動度
#define TAUP 10.0e-9      // hole life time (s)ホール寿命
#define TAUN 10.0e-9      // electron life time (s)電子寿命

/* Bernoulli function */
// ベルヌーイ関数

double bernoulli(double x)
{
  if (fabs(x) < 1e-2) //xの絶対値が0.01より小さい時
    return 1.0 - x * (0.5 - x * (1.0 / 12 - P2(x) / 720));
  else
    return x / (exp(x) - 1.0);
}

/* generation and recombination rate */
// 生成および再結合率

double gr_rate(double p, double n)
{
  return -NI * (p * n - 1.0) / (TAUP * (p + 1.0) + TAUN * (n + 1.0));
}
/*

  solve tridiagonal linear systems A x = b
  tridg(a, n)
    a[4][n]       double precision matrix
    a[0][1...n]   [in]  lower sub-diagonal part of A
    a[1][0...n]   [in]  diagonal part of A
    a[2][0...n-1] [in]  upper sub-diagonal of A
    a[3][0...n]   [in]  right hand side, b
                  [out] solution, x
  see https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
*/

void tridg(double **a, int n)
{
  int i;

  for (i = 1; i < n; ++i)
  {
    a[0][i] /= a[1][i - 1];
    a[1][i] -= a[0][i] * a[2][i - 1];
    a[3][i] -= a[0][i] * a[3][i - 1];
  }

  a[3][n - 1] /= a[1][n - 1];
  for (i = n - 2; i >= 0; --i)
    a[3][i] = (a[3][i] - a[2][i] * a[3][i + 1]) / a[1][i];
}

/* solve Poisson equation */

double *zeros(int n)
{
  double *a;

  if ((a = (double *)calloc((size_t)n, sizeof(double))) == NULL)
    exit(-1);

  return a;
}

void poisson(double eta, double *p, double *n, double *c,
             double *psi, double *dpsi, int m)
{

  double *a[4];
  int i;

  for (i = 0; i < 4; ++i)
    a[i] = zeros(m - 1);

  for (i = 1; i < m - 1; ++i)
    a[0][i] = -1.0;

  for (i = 0; i < m - 1; ++i)
    a[1][i] = 2.0 + eta * (p[i + 1] + n[i + 1]);

  for (i = 0; i < m - 2; ++i)
    a[2][i] = -1.0;

  for (i = 0; i < m - 1; ++i)
    a[3][i] = eta * (p[i + 1] - n[i + 1] + c[i + 1]) + psi[i] - 2.0 * psi[i + 1] + psi[i + 2];

  tridg(a, m - 1);
  for (i = 0; i < m - 1; ++i)
    dpsi[i + 1] = a[3][i];

  for (i = 0; i < 4; ++i)
    free(a[i]);
}

/* solve continuity equation */

void continuity(int sign, double *psi, double *gr, double *pn, int m)
{
  double *a[4];
  int i;

  for (i = 0; i < 4; ++i)
    a[i] = zeros(m - 1);

  for (i = 1; i < m - 1; ++i)
    a[0][i] = -bernoulli(sign * (psi[i + 1] - psi[i]));

  for (i = 0; i < m - 1; ++i)
    a[1][i] = bernoulli(sign * (psi[i + 2] - psi[i + 1])) + bernoulli(sign * (psi[i] - psi[i + 1]));

  for (i = 0; i < m - 2; ++i)
    a[2][i] = -bernoulli(sign * (psi[i + 1] - psi[i + 2]));

  for (i = 0; i < m - 1; ++i)
    a[3][i] = gr[i + 1];

  // boundary condition

  a[3][0] += bernoulli(sign * (psi[1] - psi[0])) * pn[0];
  a[3][m - 2] += bernoulli(sign * (psi[m - 1] - psi[m])) * pn[m];

  tridg(a, m - 1);
  for (i = 0; i < m - 1; ++i)
    pn[i + 1] = a[3][i];

  for (i = 0; i < 4; ++i)
  {
    free(a[i]);
  }
}

/* evaluate current */

void current(int sign, double *psi, double *pn, double *Jpn, int m)
{
  int i;

  for (i = 0; i < m; ++i)
    Jpn[i] = bernoulli(sign * (psi[i] - psi[i + 1])) * pn[i + 1] - bernoulli(sign * (psi[i + 1] - psi[i])) * pn[i];
}

/* read condition */

#define STRBUF 512

double read_double()
{
  char line[STRBUF];
  double d;

  if (fgets(line, STRBUF, stdin) == NULL)
    exit(-1);

  sscanf(line, "%lf", &d);

  return d;
}

double read_int()
{
  char line[STRBUF];
  int i;

  if (fgets(line, STRBUF, stdin) == NULL)
    exit(-1);

  sscanf(line, "%d", &i);

  return i;
}

/* main routine */

main()
{
  double tolerance; // potential tolerance (V)

  double Vstp; // applied bias step (V)
  int Vnum;    // number of bias points
  double L;    // device length (m)
  int N;       // number of mesh points
  double Nd;   // doping density in n-region (/m3)
  double Na;   // doping density in p-region (/m3)

  double *p;    // hole density (/m3)
  double *n;    // electron density (/m3)
  double *c;    // net impurity density Nd - Na (/m3)
  double *psi;  // potential (V)
  double *dpsi; // potential update (V)
  double *Jp;   // hole current density (A/m2)
  double *Jn;   // electron current density (A/m2)
  double *gr;   // net GR rate (/s)
  double *gp;   // normalized GR rate, Gamma_p
  double *gn;   // normalized GR rate, Gamma_n

  double h, Vapp, psi0, eta, nni, ujp, ujn, pfp, pfn, residue;
  int loop, i, iv;
  FILE *fp;

  // read conditions

  L = read_double();
  N = read_int();
  Nd = read_double();
  Na = read_double();
  Vstp = read_double();
  Vnum = read_int();
  tolerance = read_double();

  // set up the system

  h = L / N; // mesh width (m)

  p = zeros(N + 1);
  n = zeros(N + 1);
  c = zeros(N + 1);

  psi = zeros(N + 1);
  dpsi = zeros(N + 1);

  Jp = zeros(N);
  Jn = zeros(N);

  gr = zeros(N + 1);
  gp = zeros(N + 1);
  gn = zeros(N + 1);

  for (i = 0; i < (N + 1) / 2; ++i)
    c[i] = -Na; // p-region
  for (i = (N + 1) / 2; i < N + 1; ++i)
    c[i] = Nd; // n-region

  // initial guess for p & n

  for (i = 0; i < N + 1; ++i)
  {
    p[i] = 0.5 * (sqrt(P2(c[i]) + 4 * P2(NI)) - c[i]);
    n[i] = 0.5 * (sqrt(P2(c[i]) + 4 * P2(NI)) + c[i]);
  }

  // constants for normalization

  eta = EC / h / DC / VT; // for Poisson

  ujp = EC * MUP * VT / P4(h); // unit of hole current
  ujn = EC * MUN * VT / P4(h); // unit of electron current

  pfp = P5(h) / (MUP * VT); // pre-factor for Gamma_p
  pfn = P5(h) / (MUN * VT); // pre-factor for Gamma_n

  nni = NI * P3(h); // normalized intrinsic density

  // normalization

  for (i = 0; i < N + 1; ++i)
  {
    p[i] *= P3(h);
    n[i] *= P3(h);
    c[i] *= P3(h);
  }

  // initial guess for normalized potential

  for (i = 0; i < N + 1; ++i)
    psi[i] = log(n[i] / nni);
  psi0 = psi[0];

  // main loop

  for (iv = 0; iv <= Vnum; ++iv)
  {

    Vapp = iv * Vstp;
    psi[0] = psi0 + Vapp / VT; // normalized applied voltage

    // self-consistent loop

    for (loop = 1;; ++loop)
    {

      // generation-recombination rate

      for (i = 0; i < N + 1; ++i)
      {
        gr[i] = gr_rate(p[i] / nni, n[i] / nni);
        gp[i] = pfp * gr[i];
        gn[i] = pfn * gr[i];
      }

      // solve continuity equations

      continuity(1, psi, gp, p, N);
      continuity(-1, psi, gn, n, N);

      // solve Poisson equation

      poisson(eta, p, n, c, psi, dpsi, N);
      for (i = 1; i < N; ++i)
        psi[i] = psi[i] + dpsi[i];

      // converge?

      residue = 0.0;
      for (i = 1; i < N; ++i)
        residue += P2(dpsi[i]);
      residue = sqrt(residue / (N - 1)) * VT; // (V)
      if (residue < tolerance)
        break;
    }

    // evaluate current

    current(1, psi, p, Jp, N);
    for (i = 0; i < N; ++i)
      Jp[i] *= -ujp;
    current(-1, psi, n, Jn, N);
    for (i = 0; i < N; ++i)
      Jn[i] *= ujn;

    printf("%15.8e %15.8e %15.8e %15.8e %3d\n",
           Vapp, Jp[0], Jn[0], Jp[0] + Jn[0], loop);
  }

  // restore the units

  for (i = 0; i < N + 1; ++i)
  {
    p[i] /= P3(h);
    n[i] /= P3(h);
    c[i] /= P3(h);
    psi[i] *= VT;
  }

  // save results

  if ((fp = fopen("dd_out_m.txt", "w")) == NULL)
    exit(-1);
  for (i = 0; i < N + 1; ++i)
    fprintf(fp, "%15.8e %15.8e %15.8e %15.8e %15.8e\n",
            i * h, psi[i], p[i], n[i], p[i] - n[i] + c[i]);
  fclose(fp);

  if ((fp = fopen("dd_out_a.txt", "w")) == NULL)
    exit(-1);
  for (i = 0; i < N; ++i)
    fprintf(fp, "%15.8e %15.8e %15.8e %15.8e %15.8e\n",
            (i + 0.5) * h, -(psi[i + 1] - psi[i]) / h,
            Jp[i], Jn[i], Jp[i] + Jn[i]);
  fclose(fp);

  return 0;
}
