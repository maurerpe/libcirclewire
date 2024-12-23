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
#include "vatti.h"
#include "wire.h"

static float NormAngPos(float ang) {
  if (ang < 0)
    ang += 2 * M_PI;
  if (ang >= 2 * M_PI)
    ang -= 2 * M_PI;
  return ang;
}

static void TatTangent(const float *a, const float *m, const float *b, float alpha_a, float alpha_b, float *ta, float *tb) {
  float end_a, beg_b, diff, ha, hb, tta, ttb, pt_a[2], pt_b[2], ang;
  int cnt;
  
  end_a = EvalAngle(a, m, alpha_a, 1);
  beg_b = EvalAngle(m, b, alpha_b, 0);
  diff = beg_b - end_a;
  
  if (diff >= -1e-5) {
    *ta = 1;
    *tb = 0;
    return;
  }
  
  if (alpha_a == 0 && alpha_b == 0) {
    *ta = 0;
    *tb = 1;
    return;
  }
  
  ha = 2 * asin(alpha_a);
  hb = 2 * asin(alpha_b);
  tta = 0.5;
  ttb = 0.5;
  for (cnt = 0; cnt < 20; cnt++) {
    EvalSeg(pt_a, a, m, alpha_a, tta);
    EvalSeg(pt_b, m, b, alpha_b, ttb);
    ang = atan2(pt_b[1] - pt_a[1], pt_b[0] - pt_a[0]);
    if (alpha_a == 0)
      tta = 0;
    else
      tta = fmaxf(0, 1 - NormAngPos(end_a - ang) / ha);
    if (alpha_b == 0)
      ttb = 1;
    else
      ttb = fminf(1, NormAngPos(ang - beg_b) / hb);
  }
  
  *ta = tta;
  *tb = ttb;
}

static int HullAdd(struct lcw_wire *wire, const float *pt, float alpha, int slot) {
  struct lcw_seg *seg;
  float ta, tb, npt[2];
  int ret;
  
  if (alpha < 0)
    alpha = 0;
  
  while (wire->num_seg[slot] > 1) {
    seg = &wire->seg[slot][wire->num_seg[slot] - 1];
    TatTangent(seg[-1].pt, seg->pt, pt, seg->alpha, alpha, &ta, &tb);
    if (ta >= 1)
      // Already convex
      break;
    if (tb < 1) {
      SplitSeg(npt, NULL, &alpha, seg->pt, pt, alpha, tb);
    } else {
      alpha = 0;
    }
    if (ta > 0) {
      // Tangent found on existing curve
      SplitSeg(seg->pt, &seg->alpha, NULL, seg[-1].pt, seg->pt, seg->alpha, ta);
      if (tb < 1 && (ret = To(wire, npt, 0, slot)) != LCW_NO_ERROR)
	return ret;
      break;
    }
    
    // Need to roll back wire
    wire->num_seg[slot]--;
    if (tb < 1 && (ret = HullAdd(wire, npt, 0, slot)) != LCW_NO_ERROR)
      return ret;
  }
  
  return To(wire, pt, alpha, slot);
}

static void NegY(float *dst, const float *src) {
  dst[0] =  src[0];
  dst[1] = -src[1];
}

static int Hull(struct lcw_wire *hull, const struct lcw_wire *wire) {
  struct lcw_seg *seg;
  float pt[2];
  size_t idx;
  int ret;

  seg = &wire->seg[0][0];
  Reset(hull, seg->pt, 1);
  MirrorYTop(hull);
  
  for (idx = 1; idx < wire->num_seg[0]; idx++) {
    seg = &wire->seg[0][idx];
    if ((ret = HullAdd(hull, seg->pt, seg->alpha, 0)) < 0)
      return ret;
  }
  
  for (idx = 1; idx < wire->num_seg[1]; idx++) {
    seg = &wire->seg[1][idx];
    NegY(pt, seg->pt);
    if ((ret = HullAdd(hull, pt, -seg->alpha, 1)) < 0)
      return ret;
  }
  
  MirrorYTop(hull);
  
  return LCW_NO_ERROR;
}

struct lcw_wire *LCW_ConvexHull(const struct lcw_wire *wire) {
  struct lcw_wire *hull, *cp;
  
  if ((hull = LCW_New(wire->comb_tol)) == NULL)
    goto err;
  
  if (IS_VATTI(wire)) {
    if (Hull(hull, wire) < 0)
      goto err2;
  } else {
    if ((cp = LCW_Copy(wire)) == NULL)
      goto err2;
    
    if (LCW_VattiClip(cp) < 0)
      goto err3;
  
    if (Hull(hull, cp) < 0)
      goto err3;
    
    LCW_Free(cp);
  }
  
  return hull;

 err3:
  LCW_Free(cp);
 err2:
  LCW_Free(hull);
 err:
  return NULL;
}
