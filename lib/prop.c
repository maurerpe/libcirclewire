/*
  Copyright (C) 2024 Paul Maurer

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  
  1. Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
  
  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
  
  3. Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <math.h>
#include <string.h>

#include "util.h"
#include "wire.h"

#define X  LCW_X
#define Y  LCW_Y
#define XY LCW_XY

static void Add_Ext(struct lcw_properties *prop, float val, int axis) {
  if (val < prop->min[axis])
    prop->min[axis] = val;
  if (val > prop->max[axis])
    prop->max[axis] = val;
}

static void Ext_Straight(struct lcw_properties *prop, const float *p1, const float *p2) {
  Add_Ext(prop, p1[X], X);
  Add_Ext(prop, p1[Y], Y);
  Add_Ext(prop, p2[X], X);
  Add_Ext(prop, p2[Y], Y);
}

static void Com_Straight(struct lcw_properties *prop, const float *ref, const float *p1, const float *p2) {
  float x1, y1, x2, y2, aa;

  x1 = p1[X] - ref[X];
  x2 = p2[X] - ref[X];
  y1 = p1[Y] - ref[Y];
  y2 = p2[Y] - ref[Y];
  aa = x1 * y2 - x2 * y1;
  prop->area += 0.5f * aa;
  prop->center_of_mass[X] += 1.0/6 * (x2 + x1) * aa;
  prop->center_of_mass[Y] += 1.0/6 * (y2 + y1) * aa;
}

static void Mom2_Straight(struct lcw_properties *prop, const float *p1, const float *p2) {
  float x1, x2, y1, y2, aa, xx, yy;
  
  x1 = p1[X] - prop->center_of_mass[X];
  x2 = p2[X] - prop->center_of_mass[X];
  y1 = p1[Y] - prop->center_of_mass[Y];
  y2 = p2[Y] - prop->center_of_mass[Y];
  
  aa = x1 * y2 - x2 * y1;
  xx = x2 * x2 + x1 * x2 + x1 * x1;
  yy = y2 * y2 + y1 * y2 + y1 * y1;
  prop->second_moment[X] += 1.0f/12 * aa * yy;
  prop->second_moment[Y] += 1.0f/12 * aa * xx;
  prop->second_moment[XY] += 1.0f/12 * aa * (x2 * y2 + x1 * y1 + 0.5f * (x1 * y2 + x2 * y1));
}

static void Mom3x_Straight(struct lcw_properties *prop, const float *p1, const float *p2) {
  struct lcw_properties pp;
  float cx, cy, x1, x2, y1, y2, aa, x1_2, x2_2, y1_2, y2_2;
  
  cx = prop->center_of_mass[X];
  cy = prop->center_of_mass[Y];
  
  memset(&pp, 0, sizeof(pp));
  pp.center_of_mass[X] = cx;
  pp.center_of_mass[Y] = cy;
  Mom2_Straight(&pp, p1, p2);
  
  x1 = p1[X] - cx;
  x2 = p2[X] - cx;
  y1 = p1[Y] - cy;
  y2 = p2[Y] - cy;
  
  aa = x1 * y2 - x2 * y1;
  x1_2 = x1 * x1;
  x2_2 = x2 * x2;
  y1_2 = y1 * y1;
  y2_2 = y2 * y2;
  
  prop->third_moment_x[X]  += 1.0f/20 * aa * (y2 + y1) * (y2_2 + y1_2) + cy * pp.second_moment[X];
  prop->third_moment_x[Y]  += 1.0f/60 * aa * (3 * (x2_2 * y2 + x1_2 * y1) + 2 * x1 * x2 * (y2 + y1) + x1_2 * y2 + x2_2 * y1) + cy * pp.second_moment[Y];
  prop->third_moment_x[XY] += 1.0f/60 * aa * (3 * (x2 * y2_2 + x1 * y1_2) + 2 * y1 * y2 * (x1 + x2) + x1 * y2_2 + x2 * y1_2) + cy * pp.second_moment[XY];
}

static float EvalPoly(const float *coefs, float x, int num) {
  float ans = *coefs++;
  
  while (--num > 0)
    ans = ans * x + *coefs++;
  
  return ans;
}

#define EvalPolyM(c, x) EvalPoly((c),(x),sizeof(c)/sizeof(*(c)))

static int Between(float val, float lim1, float lim2) {
  float min, max;
  if (lim1 < lim2) {
    min = lim1;
    max = lim2;
  } else {
    min = lim2;
    max = lim1;
  }
  
  return val > min && val < max;
}

static void Ext_Circle(struct lcw_properties *prop, const float *p1, const float *p2, float alpha) {
  float dist, sq, center[2], rad;
  
  dist = Dist(p1, p2);
  sq = sqrtf(1 - alpha * alpha) / alpha;
  center[X] = 0.5f * (p1[X] + p2[X] - sq * (p2[Y] - p1[Y]));
  center[Y] = 0.5f * (p1[Y] + p2[Y] + sq * (p2[X] - p1[X]));
  rad = dist / (2 * alpha);
  if (Between(center[X], p1[X], p2[X]))
    Add_Ext(prop, center[Y] + (p2[X] > p1[X] ? -rad :  rad), Y);
  if (Between(center[Y], p1[Y], p2[Y]))
    Add_Ext(prop, center[X] + (p2[Y] > p1[Y] ?  rad : -rad), X);
}

static const float area_coefs[] =
  {0.006160960477941176f,
   0.00751953125f,
   0.00946514423076923f,
   0.01242897727272727f,
   0.01736111111111111f,
   0.02678571428571428f,
   0.05f,
   0.1666666666666667f};

static const float momx_coefs[] =
  {0.001970647171885562f,
   0.002340143516614105f,
   0.002841602841602842f,
   0.003552003552003552f,
   0.004617604617604618f,
   0.006349206349206349f,
   0.009523809523809525f,
   0.01666666666666667f};
   
static void Com_Circle(struct lcw_properties *prop, const float *ref, const float *p1, const float *p2, float alpha) {
  float dx, dy, d2, a2, asn, sqa, d, d3, a3, area, momx, mx, my;
  
  dx = p2[X] - p1[X];
  dy = p2[Y] - p1[Y];
  d2 = dx * dx + dy * dy;
  d  = sqrtf(d2);
  d3 = d2 * d;
  a2 = alpha * alpha;
  
  if (fabsf(alpha) < 0.3f) {
    area = d2 * EvalPolyM(area_coefs, a2) * alpha;
    momx = d3 * EvalPolyM(momx_coefs, a2) * a2;
  } else {
    a3 = alpha * a2;
    asn = asinf(alpha);
    sqa = sqrtf(1 - a2);
    area = 0.25f * d2 * (asn - alpha * sqa) / a2;
    momx = 1.0f/24 * d3 * (3 * alpha - a3 - 3 * sqa * asn) / a3;
  }

  mx = 0.5f * (p1[X] + p2[X]) - ref[X];
  my = 0.5f * (p1[Y] + p2[Y]) - ref[Y];
  prop->area += area;
  prop->center_of_mass[X] += mx * area + dy / d * momx;
  prop->center_of_mass[Y] += my * area - dx / d * momx;
}

static float momxx_coefs[] =
  {2.595495294641565e-4f,
   2.936617533365886e-4f,
   3.359974294468976e-4f,
   3.896440778459821e-4f,
   4.593698601973684e-4f,
   5.529067095588236e-4f,
   6.8359375e-4f,
   8.764022435897436e-4f,
   0.001183712121212121f,
   0.001736111111111111f,
   0.002976190476190476f,
   0.008333333333333333};

static float momyy_coefs[] =
  {1.029471438286337e-4f,
   1.136350482997966e-4f,
   1.262968699246397e-4f,
   1.414784860192964e-4f,
   1.599334701454077e-4f,
   1.827261008126603e-4f,
   2.113978878783035e-4f,
   2.482430146559123e-4f,
   2.967683132721398e-4f,
   3.624307835022121e-4f,
   4.535095856524428e-4f,
   5.793650793650793e-4f,
   7.142857142857143e-4f};

static void Mom2_Circle(struct lcw_properties *prop, const float *p1, const float *p2, float alpha) {
  float dx, dy, d2, a2, asn, sqa, d, d3, d4, a3, a4;
  float area, momx, comx, momxx, momyy;
  float mx, my, rx, ry, rx2, ry2, rxy;
  
  dx = p2[X] - p1[X];
  dy = p2[Y] - p1[Y];
  d2 = dx * dx + dy * dy;
  d  = sqrtf(d2);
  d3 = d2 * d;
  d4 = d2 * d2;
  a2 = alpha * alpha;
  a3 = alpha * a2;
  a4 = a2 * a2;

  if (fabsf(alpha) < 0.3f) {
    area = d2 * EvalPolyM(area_coefs, a2) * alpha;
    momx = d3 * EvalPolyM(momx_coefs, a2) * a2;
  } else {
    a3 = alpha * a2;
    asn = asinf(alpha);
    sqa = sqrtf(1 - a2);
    area = 0.25f * d2 * (asn - alpha * sqa) / a2;
    momx = 1.0f/24 * d3 * (3 * alpha - a3 - 3 * sqa * asn) / a3;
  }
  comx = momx / area;
  if (fabsf(alpha) < 0.5f) {
    momxx = d4 * EvalPolyM(momxx_coefs, a2) * alpha;
    momyy = d4 * EvalPolyM(momyy_coefs, a2) * a3;
  } else {
    asn = asinf(alpha);
    sqa = sqrtf(1 - a2);
    momxx = 1.0f/192 * d4 * (3 * asn - sqa * alpha * (2 * a2 + 3)) / a4;
    momyy = 1.0f/192 * d4 * (-12 * a2 * asn + 15 * asn + sqa * alpha * (2 * a2 - 15)) / a4 - momx * comx;
  }
  
  rx = dx / d;
  ry = dy / d;
  rx2 = rx * rx;
  ry2 = ry * ry;
  rxy = rx * ry;
  mx = 0.5f * (p1[X] + p2[X]) - prop->center_of_mass[X] + ry * comx;
  my = 0.5f * (p1[Y] + p2[Y]) - prop->center_of_mass[Y] - rx * comx;
  
  prop->second_moment[X]  += ry2 * momxx + rx2 * momyy + area * my * my;
  prop->second_moment[Y]  += rx2 * momxx + ry2 * momyy + area * mx * mx;
  prop->second_moment[XY] += rxy * (momxx - momyy) + area * mx * my;
}

static float momyyraw_coefs[] =
  {3.735174680136291e-4f,
   4.116992536328001e-4f,
   4.568071718569156e-4f,
   5.107160927592844e-4f,
   5.759955933375389e-4f,
   6.56242657424812e-4f,
   7.566091815015479e-4f,
   8.846507352941176e-4f,
   0.001051682692307692f,
   0.0012747668997669f,
   0.001578282828282828f,
   0.001984126984126984f,
   0.002380952380952381f};

static float momxxy_coefs[] =
  {5.814951529750177e-5f,
   6.43798205079484e-5f,
   7.180826133578859e-5f,
   8.078429400276216e-5f,
   9.180033409404791e-5f,
   1.055703842081551e-4f,
   1.231654482428476e-4f,
   1.462589697883816e-4f,
   1.776001776001776e-4f,
   2.22000222000222e-4f,
   2.886002886002886e-4f,
   3.968253968253968e-4f,
   5.952380952380953e-4f};

static void Mom3x_Circle(struct lcw_properties *prop, const float *p1, const float *p2, float alpha) {
  struct lcw_properties pp;
  float cx, cy, v, qx, qy;
  float dx, dy, d2, d, d3, d4, d5, a2, a3, a4, a5;
  float A, Ix, Ixx, Iyyraw, Ixxy, asn, sqa;
  float cxx, my, cyy, st, ct, st2, ct2, cst;
  float prod1, prod2, sub;
  
  cx  = prop->center_of_mass[X];
  cy  = prop->center_of_mass[Y];

  memset(&pp, 0, sizeof(pp));
  pp.center_of_mass[X] = cx;
  pp.center_of_mass[Y] = cy;
  Mom2_Circle(&pp, p1, p2, alpha);
  
  v  = pp.area * cy;
  qx = pp.second_moment[X]  / v + cx;
  qy = pp.second_moment[XY] / v + cy;
  
  dx = p2[X] - p1[X];
  dy = p2[Y] - p1[Y];
  d2 = dx * dx + dy * dy;
  d  = sqrtf(d2);
  d3 = d2 * d;
  d4 = d2 * d2;
  d5 = d4 * d;
  a2 = alpha * alpha;
  a3 = a2 * alpha;
  a4 = a2 * a2;
  a5 = a4 * alpha;
  
  if (fabsf(alpha) < 0.3f) {
    A = d2 * EvalPolyM(area_coefs, a2) * alpha;
    Ix = d3 * EvalPolyM(momx_coefs, a2) * a2;
  } else {
    asn = asinf(alpha);
    sqa = sqrtf(1 - a2);
    A = 0.25f * d2 * (asn - alpha * sqa) / a2;
    Ix = 1.0f/24 * d3 * (3 * alpha - a3 - 3 * sqa * asn) / a3;
  }
  if (fabsf(alpha) < 0.5f) {
    Ixx    = d4 * EvalPolyM(momxx_coefs, a2) * alpha;
    Iyyraw = d4 * EvalPolyM(momyyraw_coefs, a2) * a3;
    Ixxy   = d5 * EvalPolyM(momxxy_coefs, a2) * a2;
  } else {
    asn = asinf(alpha);
    sqa = sqrtf(1 - a2);
    Ixx    = 1.0f/192 * d4 * (3 * asn - sqa * alpha * (2 * a2 + 3)) / a4;
    Iyyraw = 1.0f/192 * d4 * (-12 * a2 * asn + 15 * asn + 2 * a3 * sqa - 15 * alpha * sqa) / a4;
    Ixxy   = 1.0f/1920 * d5 * (15 * alpha + 15 * sqa * asn - 2 * a5 - 5 * a3) / a5;
  }
  
  cxx = 0.5f * (p1[X] + p2[X]) - qx;
  my  = 0.5f * (p1[Y] + p2[Y]);
  cyy = my - qy;
  st  = dx / d;
  ct  = dy / d;
  st2 = st * st;
  ct2 = ct * ct;
  cst = ct * st;
  
  prod1 = (2 * cyy - d * st) * (2 * cyy + d * st);
  
  prop->third_moment_x[X] +=
    0.5f * A * cyy * d2 * st2 +
    0.25f * A * my * prod1 +
    0.25f * Ix * st * prod1 +
    -2 * Ix * my * cyy * st +
    2 * Ixx * cyy * (4 * ct2 - 5) +
    Ixx * my * (4 * st2 + 1) +
    Ixxy * st * (4 * st2 + 3);
  
  prod2 = (2 * cxx - d * ct) * (2 * cxx + d * ct);
  
  prop->third_moment_x[Y] +=
    0.5f * A * cxx * d2 * cst +
    0.25f * A * my * prod2 +
    0.25f * Ix * st * prod2 +
    -2 * Ix * my * cxx * ct +
    -8 * Ixx * cxx * cst +
    Ixx * my * (4 * ct2 + 1) +
    Ixxy * st * (4 * ct2 + 1);
  
  sub = 4 * cxx * cyy - d2 * cst;
  
  prop->third_moment_x[XY] +=
    0.25f * A * cyy * d2 * cst +
    0.25f * A * my * sub +
    0.25f * Ix * st * sub +
    -Ix * my * (cxx * st + cyy * ct) +
    -Ixx * ct * (4 * cyy * st + cxx * ct) +
    4 * Ixx * my * cst +
    -Iyyraw * cxx * st2 +
    Ixxy * (4 * st2 + 1);
}

static void Prop(struct lcw_properties *prop, const struct lcw_wire *wire) {
  struct lcw_seg *cur, *next;
  size_t count;
  float ref[2];
  
  memset(prop, 0, sizeof(*prop));
  prop->min[X] = INFINITY;
  prop->min[Y] = INFINITY;
  prop->max[X] = -INFINITY;
  prop->max[Y] = -INFINITY;
  for (count = 0; count < wire->num_seg[0] - 1; count++) {
    cur  = &wire->seg[0][count];
    next = &wire->seg[0][count + 1];
    Ext_Straight(prop, cur->pt, next->pt);
    if (next->alpha != 0)
      Ext_Circle(prop, cur->pt, next->pt, next->alpha);
  }
  
  ref[X] = 0.5f * (prop->max[X] + prop->min[X]);
  ref[Y] = 0.5f * (prop->max[Y] + prop->min[Y]);
  
  for (count = 0; count < wire->num_seg[0] - 1; count++) {
    cur  = &wire->seg[0][count];
    next = &wire->seg[0][count + 1];
    Com_Straight(prop, ref, cur->pt, next->pt);
    if (next->alpha != 0)
      Com_Circle(prop, ref, cur->pt, next->pt, next->alpha);
  }
  
  prop->center_of_mass[X] = ref[X] + prop->center_of_mass[X] / prop->area;
  prop->center_of_mass[Y] = ref[Y] + prop->center_of_mass[Y] / prop->area;
  
  for (count = 0; count < wire->num_seg[0] - 1; count++) {
    cur  = &wire->seg[0][count];
    next = &wire->seg[0][count + 1];
    Mom2_Straight(prop, cur->pt, next->pt);
    if (next->alpha != 0)
      Mom2_Circle(prop, cur->pt, next->pt, next->alpha);
  }
  
  for (count = 0; count < wire->num_seg[0] - 1; count++) {
    cur  = &wire->seg[0][count];
    next = &wire->seg[0][count + 1];
    Mom3x_Straight(prop, cur->pt, next->pt);
    if (next->alpha != 0)
      Mom3x_Circle(prop, cur->pt, next->pt, next->alpha);
  }
}

int LCW_Properties(struct lcw_properties *prop, const struct lcw_wire *wire) {
  struct lcw_wire *cp;
  int ret;
  
  if (!LCW_IsClosed(wire))
    return LCW_WIRE_NOT_CLOSED;
  
  if (!IS_VATTI(wire)) {
    Prop(prop, wire);
  } else {
    if ((cp = LCW_Copy(wire)) == NULL) {
      ret = LCW_OUT_OF_MEMORY;
      goto err;
    }
    
    if ((ret = LCW_WireLoop(cp)) != LCW_NO_ERROR)
      goto err2;
  
    Prop(prop, cp);
    
    LCW_Free(cp);
  }
  
  return LCW_NO_ERROR;

 err2:
  LCW_Free(cp);
 err:
  return ret;
}