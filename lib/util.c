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

#include "util.h"

float Dot(const float *a, const float *b) {
  return a[0] * b[0] + a[1] * b[1];
}

float Dist2(const float *a, const float *b) {
  float dx, dy;

  dx = b[0] - a[0];
  dy = b[1] - a[1];

  return dx * dx + dy * dy;
}

float Dist(const float *a, const float *b) {
  return sqrtf(Dist2(a, b));
}

float Norm2(const float *a) {
  float x, y;

  x = a[0];
  y = a[1];

  return x * x + y * y;
}

float Norm(const float *a) {
  return sqrtf(Norm2(a));
}

void Normalize(float *a) {
  float norm;

  norm = Norm(a);
  if (norm == 0)
    return;
  a[0] /= norm;
  a[1] /= norm;
}

float ArcCenter(float *cent, const float *a, const float *b, float alpha) {
  float delta[2], dd, rad, leg, perp[2];
  
  delta[0] = b[0] - a[0];
  delta[1] = b[1] - a[1];
  dd  = Norm(delta);
  rad = dd / (2 * alpha);
  leg = rad * sqrtf(1 - alpha * alpha);
  perp[0] = -delta[1] / dd;
  perp[1] =  delta[0] / dd;
  cent[0] = (a[0] + b[0]) / 2 + leg * perp[0];
  cent[1] = (a[1] + b[1]) / 2 + leg * perp[1];
  
  return rad;
}

size_t RevRequiredSegs(float ang, float rad, float tol) {
  float ratio, denom, num;
  
  if (ang == 0 || rad == 0)
    return 1;

  ratio = fabsf(tol / rad);
  if (ratio >= 2)
    return 1;
  if (ratio < 0.01)
    denom = sqrtf(2 * ratio) * ((3.0f/160 * ratio + 1.0f/12) * ratio + 1);
  else
    denom = acosf(1 - ratio);
  num = ceilf(0.5 * fabsf(ang) / denom);
  
  if (num < 1)
    return 1;
  if (num > SIZE_MAX)
    return SIZE_MAX;
  return num;
}

size_t RequiredSegs(float dist, float alpha, float tol) {
  return RevRequiredSegs(2*asinf(alpha), dist / (2 * alpha), tol);
}

static float SDA(float t, float a) {
  float a2, t2, c2, c1;
  
  a2 = a * a;
  t2 = t * t;
  c2 = (4 * t2 * t2 - 4) * t / 15;
  c1 = (-4 * t2 - 2) * t / 3;
  return (c2 * a2 + c1) * a2 + 2 * t;
}

static float CDA(float t, float a) {
  float a2, t2, c2, c1, c0;
  
  a2 = a * a;
  t2 = t * t;
  c2 = ((4 * t2 - 20) * t2 + 16) * t2 / 45;
  c1 = (2 * t2 - 2) * t2 / 3;
  c0 = 2 * t2;
  return -((c2 * a2 + c1) * a2 + c0) * a;
}

void EvalSeg(float *pt, const float *a, const float *b, float alpha, float t) {
  float sqa, tas, ssa, cca, sda, cda, in1, in2, xxx, ox, oy;

  sqa = sqrtf(1 - alpha * alpha);
  tas = 2 * t * asinf(alpha);
  ssa = sinf(tas);
  cca = cosf(tas);
  if (fabsf(alpha) >= 0.1) {
    sda = sqa * ssa / alpha;
    cda = sqa * (cca - 1) / alpha;
  } else {
    sda = SDA(t, alpha);
    cda = sqa * CDA(t, alpha);
  }
  in1 = 1 + cca - sda;
  in2 = 1 + sda - cca;
  xxx = ssa + cda;

  ox = in1 * a[0] + in2 * b[0] + xxx * (b[1] - a[1]);
  oy = in1 * a[1] + in2 * b[1] - xxx * (b[0] - a[0]);
  pt[0] = ox / 2;
  pt[1] = oy / 2;
}

void SplitSeg(float *pt, float *a1, float *a2, const float *a, const float *b, float alpha, float tt) {
  float hang;

  EvalSeg(pt, a, b, alpha, tt);

  hang = asinf(alpha);
  if (a1)
    *a1 = sinf(tt * hang);
  if (a2)
    *a2 = sinf((1 - tt) * hang);
}

float EvalAngle(const float *a, const float *b, float alpha, float t) {
  float dx, dy, base, as, ang;
  
  dx = b[0] - a[0];
  dy = b[1] - a[1];
  base = atan2f(dy, dx);
  as = 2 * asinf(alpha);
  ang = base + (t - 0.5) * as;
  
  if (ang > M_PI)
    ang -= 2 * M_PI;
  if (ang < -M_PI)
    ang += 2 * M_PI;
  
  return ang;
}

float ArcLen(const float *a, const float *b, float alpha) {
  float dd;

  dd = Dist(a, b);
  if (alpha == 0)
    return dd;
  return dd * asinf(alpha) / alpha;
}

static float Limit01(float val) {
  return fmaxf(0, fminf(1, val));
}

static float Wrap(float val) {
  float mm;
  
  mm = fmodf(val + M_PI, 2 * M_PI);
  if (mm < 0)
    mm += 2 * M_PI;
  return mm - M_PI;
}

float TatX(const float *a, const float *b, float alpha, float xval) {
  float asa, dd, idd, dx, dy, ra, nx, ny, ph, yy, ac, t1, t2;
  float tt, ddx, ddy, pt[2], dxdt;
  int cnt;
  
  if (alpha == 0) {
    if (b[0] == a[0])
      return 0.5f;
    return Limit01((xval - a[0]) / (b[0] - a[0]));
  }
  
  asa = asinf(alpha);
  dd = Dist(a, b);
  if (dd == 0)
    return 0.5f;
  
  if (fabsf(alpha) >= 0.1) {
    idd = 1 / dd;
    dx = b[0] - a[0];
    dy = b[1] - a[1];
    yy = (b[1] - a[1]) * sqrt(1 - alpha * alpha);
    ra = (yy + (2 * xval - b[0] - a[0]) * alpha) * idd;
    if (alpha < 0) {
      ra = -ra;
      nx = -dy;
      ny = dx;
    } else {
      nx = dy;
      ny = -dx;
    }
    ph = -atan2(ny, nx);
    if (ra > 1)
      ac = 0;
    else if (ra < -1)
      ac = M_PI;
    else
      ac = acosf(ra);
    t1 = Wrap(ph + ac) / (2 * asa);
    t2 = Wrap(ph - ac) / (2 * asa);
    return Limit01((fabsf(t1) < fabsf(t2) ? t1 : t2) + 0.5f);
  }

  if (b[0] == a[0])
    tt = 0.25;
  else
    tt = Limit01((xval - a[0]) / (b[0] - a[0]));
  ddx = b[0] - a[0];
  ddy = b[1] - a[1];
  for (cnt = 0; cnt < 10; cnt++) {
    EvalSeg(pt, a, b, alpha, tt);
    dx = pt[0] - xval;
    dxdt = ddx - (2 * tt - 1) * ddy * alpha;
    tt = Limit01(tt - dx / dxdt);
  }

  return tt;
}
