/*
###################################################################################
#
# RIAM-COMPACT
#
# Copyright (C) 2015-2017 Research Institute for Applied Mechanics(RIAM)
#                       / Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
# Copyright (C) 2015-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
*/

/**
 * @file   riamc_stretch.cpp
 * @brief  Vinokur's stretching function
 * @author RIIT
 */

#include "stretch.h"

#define EPS 1.0e-5
#define ITER_MAX 1000


// @brief 格子分布を反復的に求める
// @param [out] coord 座標値
// @note coordのヌルポインタチェックは事前に
bool Stretch::distribution(REAL_TYPE* coord)
{
  FILE* fp;
  double n1 = (double)(numNodes-1);

  // [0,1]に正規化
  for (int j=0; j<numNodes; j++)
  {
    x [j] = (double)j / n1;
    wk[j] = (double)j / n1;
  }

  // [0,1]で分布を計算するためlengthで正規化
  double sp1 = sp_s / length;
  double sp2 = sp_e / length;

  if ( sp1>1.0 || sp2>1.0) return false;

  double s11 = 1.0 / n1 / sp1;
  double s21 = 1.0 / n1 / sp2;

  stretching(s11, s21);
  int n = numNodes;
  double dx1 = wk[1]   - wk[0];
  double dx2 = wk[n-1] - wk[n-2];

#ifdef _DEBUG
  if ( !(fp=fopen("coord.txt", "w")) )
  {
    printf("\tSorry, can't open 'coord.txt' file.\n");
    exit(1);
  }
#endif

  if ( fabs(sp1-dx1)/sp1<EPS && fabs(sp2-dx2)/sp2<EPS ) goto L10;

  int itr;
  for (itr=0; itr<ITER_MAX; itr++) {
    double ds1 = s11 * 0.1 ;
    double ds2 = s21 * 0.1 ;
    double s12 = s11 + ds1 ;
    double s22 = s21 + ds2 ;

    stretching(s12, s21);
    double dx1_ds1 = (wk[1]  -wk[0]  -dx1) / ds1;
    double dx2_ds1 = (wk[n-1]-wk[n-2]-dx2) / ds1;

    stretching(s11, s22);
    double dx1_ds2 = (wk[1]  -wk[0]  -dx1) / ds2;
    double dx2_ds2 = (wk[n-1]-wk[n-2]-dx2) / ds2;

    double df1 = sp1 - dx1;
    double df2 = sp2 - dx2;

    double d = dx1_ds1 * dx2_ds2 - dx1_ds2 * dx2_ds1;
    ds1 = (df1 * dx2_ds2 - df2 * dx1_ds2) / d;
    ds2 = (df2 * dx1_ds1 - df1 * dx2_ds1) / d;

    s12 = s11 + ds1 ;
    s22 = s21 + ds2 ;

    stretching(s12, s22);
    dx1 = wk[1]   - wk[0];
    dx2 = wk[n-1] - wk[n-2];

#ifdef _DEBUG
    fprintf(fp,"%4d: ", itr);
    for (int i=0; i<numNodes; i++)
    {
      fprintf(fp,"%8.5lf ", wk[i]);
    }
    fprintf(fp,"\n");
#endif

    if ( fabs(sp1-dx1)/sp1<EPS && fabs(sp2-dx2)/sp2<EPS ) break;

    s11 = s12;
    s21 = s22;
  }

  L10:
    wk[0]   = 0.0;
    wk[n-1] = 1.0;


  for (int j=0; j<numNodes; j++)
  {
    coord[j] = (REAL_TYPE) (x_s + wk[j] * length);
  }

#ifdef _DEBUG
  fprintf(fp, "\n\nIteration = %d\n", itr);
  for (int i=0; i<numNodes; i++)
  {
    fprintf(fp,"%3d %8.5lf %8.5f ", i, wk[i], coord[i]);
    if (i==0) {
      fprintf(fp, "\n");
    }
    else {
      fprintf(fp, "%8.5f\n", coord[i]-coord[i-1]);
    }
  }
  fclose(fp);
#endif

  if (itr == ITER_MAX) return false;
  return true;
}


/**
  * @brief 両端がs0, s1の間隔の分布を近似する
  * @param [in]  s0   width of grid at start position [0, 1]
  * @param [in]  s1   width of grid at end position [0, 1]
  */
void Stretch::stretching(const double s0, const double s1)
{
  const double epsilon = 0.001;
  double dz, tanx, tanhx, u;

  double b = sqrt(s0*s1);
  double a = b / s1;

  if ( b < 1.0-epsilon )
  {
    dz = asin(b);
    for (int j=0; j<numNodes; j++)
    {
      tanx = tan( dz * x[j] );
      wk[j] = tanx / (a*sin(dz) + (1.0-a*cos(dz))*tanx);
    }
  }
  else
  {
    if ( b > 1.0+epsilon )
    {
      dz = asinh(b);
      for (int j=0; j<numNodes; j++)
      {
        tanhx = tanh( dz * x[j] ) ;
        wk[j] = tanhx / (a*sinh(dz) + (1.0-a*cosh(dz)) *tanhx);
      }
    }
    else
    {
      for (int j=0; j<numNodes; j++)
      {
        u = x[j] * (1.0+2.0*(b-1.0)*(x[j]-5.0)*(1.0-x[j])); // 5.0 > 0.5?
        wk[j] = u / ( a+(1.0-a)*u);
      }
    }
  }
}

/* deprecated
// asin
double Stretch::f_asin(const double y)
{
  const double a3 = -2.6449341;
  const double a4 = 6.794732;
  const double a5 = -13.205501;
  const double a6 = 11.726095;
  const double b1 = 0.15;
  const double b2 = 0.057321429;
  const double b3 = 0.048774283;
  const double b4 = -0.053337753;
  const double b5 = 0.075845134;
  const double u1 = 0.2693897165164;
  const double pi = 3.14159265358981;
  double ret, yb;

  if ( y < u1 )
  {
    ret = pi*((((((a6*y+a5)*y + a4) *y +a3)*y+1.0)*y-1.0)*y+1.0);
  }
  else
  {
    yb = 1.0 - y;
   ret = sqrt(6.0*yb) * (((((b5*yb+b4)*yb+b3)*yb+b2)*yb+b1)*yb+1.0);
}

  return ret;
}


// asinh
double Stretch::f_asinh(const double y)
{
  const double a1 = 0.15;
  const double a2 = 0.0573214285714;
  const double a3 = -0.024907294878;
  const double a4 = 0.0077424460899;
  const double a5 = -0.0010794122691;
  const double b0 = -0.0204176930892;
  const double b1 = 0.2490272170591;
  const double b2 = 1.9496443322775;
  const double b3 = -2.629454725241;
  const double b4 = 8.5679591096315;
  const double u1 = 2.7829681178603;
  const double u2 = 0.028527431;
  double ret, yb, v, w;

  if ( y < u1 )
  {
    yb = y - 1.0000000;
    ret = sqrt(6.0*yb)*(((((a5*yb+a4)*yb+a3)*yb+a2)*yb+a1)*yb+1.0);
  }
  else
  {
    v = log(y);
    w = 1.0 / y - u2;
    ret = v + log(2.0 * v) *(1.0 + 1.0/v)
        + (((b4 * w +b3)*w+b2)*w+b1)*w+b0;
  }

  return ret;
}
*/
