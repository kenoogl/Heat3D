#ifndef _RIAMC_STRETCH_H_
#define _RIAMC_STRETCH_H_
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
 * @file   riamc_stretch.h
 * @brief  Vinokur's stretching function
 * @author RIIT
 * @ref https://github.com/jsitaraman/tioga/
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define REAL_TYPE float

class Stretch {

private:
  double x_s;    ///< 開始点座標
  double x_e;    ///< 終了点座標
  double sp_s;   ///< 開始点側の格子点間隔
  double sp_e;   ///< 終了点側の格子点間隔
  double length; ///< 線分の長さ
  int numNodes;  ///< 分割格子点数
  double* x;     ///< 座標軸配列
  double* wk;    ///< 作業用


public:

  Stretch(const int m_node,
          const REAL_TYPE xs,
          const REAL_TYPE xe,
          const REAL_TYPE sps,
          const REAL_TYPE spe)
  {
    numNodes = m_node;
    x_s  = (double)xs;
    x_e  = (double)xe;
    sp_s = (double)sps;
    sp_e = (double)spe;
    length = fabs(x_e - x_s);

    x  = new double[numNodes];
    wk = new double[numNodes];
  };

  ~Stretch()
  {
    delete [] x;
    delete [] wk;
  };


  // @brief 格子分布を反復的に求める
  // @param [out] coord 座標値
  bool distribution(REAL_TYPE* coord);


private:

  void stretching(const double s0, const double s1);
  //double f_asin(const double y);
  //double f_asinh(const double y);
};

#endif // _RIAMC_STRETCH_H_
