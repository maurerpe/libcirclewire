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
#include "vatti.h"

#define T_TOL 1e-6f

static int SegHasVert(const float *a, const float *b, float alpha, float *tt) {
  float hang, dx, dy, ctr;
  
  if (alpha == 0) {
    if (a[0] == b[0]) {
      if (tt)
	*tt = 0;
      return 1;
    }
    
    return 0;
  }
  
  hang = asinf(alpha);
  dx = b[0] - a[0];
  dy = b[1] - a[1];
  ctr = atan2f(dy, dx);
  if (fabsf(ctr - (float)M_PI_2) < hang) {
    if (tt)
      *tt = (M_PI_2 - (ctr - hang)) / (2 * hang);
    return 1;
  }
  
  if (fabsf(ctr + (float)M_PI_2) < hang) {
    if (tt)
      *tt = (-(float)M_PI_2 - (ctr - hang)) / (2 * hang);
    return 1;
  }

  return 0;
}

#define FARTHER_PT(nn, ff, cmp) ((nn)[0] cmp (ff)[0] || ((nn)[0] == (ff)[0] && (nn)[1] cmp (ff)[1]))
#define CHECK_EXTREME(nn, ix, t, ei, cmp)	\
  do {						\
    if (FARTHER_PT(nn, (ei).pt, cmp)) {		\
      (ei).idx = ix;				\
      (ei).tt = t;				\
      (ei).pt[0] = (nn)[0];			\
      (ei).pt[1] = (nn)[1];			\
    }						\
  } while (0)

struct extreme_info {
  size_t idx;
  float tt;
  float pt[2];
};

static void ExtremeX(const struct lcw_wire *wire, int slot, size_t start, size_t stop, struct extreme_info *left_ei, struct extreme_info *right_ei) {
  const struct lcw_seg *seg;
  struct extreme_info rei, lei;
  size_t idx;
  float tt, pt[2];
  
  memset(&rei, 0, sizeof(rei));
  memset(&lei, 0, sizeof(lei));
  rei.idx = lei.idx = start;
  rei.pt[0] = lei.pt[0] = wire->seg[slot][start].pt[0];
  rei.pt[1] = lei.pt[1] = wire->seg[slot][start].pt[1];
  
  for (idx = start; idx < stop; idx++) {
    seg = &wire->seg[slot][idx];
    if (SegHasVert(seg[0].pt, seg[1].pt, seg[1].alpha, &tt) && tt > 0 && tt < 1) {
      EvalSeg(pt, seg[0].pt, seg[1].pt, seg[1].alpha, tt);
      CHECK_EXTREME(pt, idx, tt, rei, >);
      CHECK_EXTREME(pt, idx, tt, lei, <);
    }
    CHECK_EXTREME(seg[0].pt, idx, 0, rei, >);
    CHECK_EXTREME(seg[0].pt, idx, 0, lei, <);
  }
  CHECK_EXTREME(wire->seg[slot][stop].pt, stop-1, 1, rei, >);
  CHECK_EXTREME(wire->seg[slot][stop].pt, stop-1, 1, lei, <);
  
  if (right_ei)
    memcpy(right_ei, &rei, sizeof(*right_ei));
  if (left_ei)
    memcpy(left_ei, &lei, sizeof(*left_ei));  
}

static void FarthestLeft(const struct lcw_wire *wire, int slot, struct extreme_info *ei) {
  ExtremeX(wire, slot, 0, wire->num_seg[0] - 1, ei, NULL);
}

static int ClipFwd(struct lcw_wire *nn, const struct lcw_wire *oo, size_t idx, int slot) {
  struct lcw_seg *seg;
  float tt, pt[2], a1, a2;
  size_t ns, cnt;
  int ret;
  
  ns = oo->num_seg[0] - 1;
  cnt = ns;
  while (1) {
    seg = &oo->seg[0][idx];
    if (SegHasVert(seg[0].pt, seg[1].pt, seg[1].alpha, &tt) && tt > T_TOL && tt < 1 - T_TOL) {
      SplitSeg(pt, &a1, &a2, seg[0].pt, seg[1].pt, seg[1].alpha, tt);
      if ((ret = To(nn, pt, a1, slot)) < 0)
	return ret;
      break;
    }
    if (seg[1].pt[0] < seg[0].pt[0])
      break;
    if ((ret = To(nn, seg[1].pt, seg[1].alpha, slot)) < 0)
      return ret;
    idx = (idx + 1) % ns;
    if (cnt-- == 0)
      return LCW_WIRE_NOT_VATTI_SINGLE;
  }

  return LCW_NO_ERROR;
}

static int ClipRev(struct lcw_wire *nn, const struct lcw_wire *oo, size_t idx, int slot) {
  struct lcw_seg *seg;
  float tt, pt[2], a1, a2;
  size_t ns, cnt;
  int ret;
  
  ns = oo->num_seg[0] - 1;
  cnt = ns;
  while (1) {
    seg = &oo->seg[0][idx];
    if (SegHasVert(seg[0].pt, seg[1].pt, seg[1].alpha, &tt) && tt > T_TOL && tt < 1 - T_TOL) {
      SplitSeg(pt, &a1, &a2, seg[0].pt, seg[1].pt, seg[1].alpha, tt);
      if ((ret = To(nn, pt, -a2, slot)) < 0)
	return ret;
      break;
    }
    if (seg[0].pt[0] < seg[1].pt[0])
      break;
    if ((ret = To(nn, seg[0].pt, -seg[1].alpha, slot)) < 0)
      return ret;
    idx = (idx + ns - 1) % ns;
    if (cnt-- == 0)
      return LCW_WIRE_NOT_VATTI_SINGLE;
  }

  return LCW_NO_ERROR;
}

struct merge_info {
  struct lcw_wire *wire;
  size_t idx;
  int slot;
};

static void ExtremeX_Merge(const struct merge_info *mi, float *min, float *max) {
  struct extreme_info rei, lei;
  ExtremeX(mi->wire, mi->slot, mi->idx, mi->idx + 1, &lei, &rei);
  if (min)
    *min = lei.pt[0];
  if (max)
    *max = rei.pt[0];
}

static int Compatible(const struct merge_info *tinfo, const struct merge_info *binfo) {
  /* To ensure vatti singular, check that tinfo is never below binfo */
  const struct lcw_seg *tseg, *bseg;
  float tol, tmin, tmax, bmin, bmax, val, ttt, btt, tpt[2], bpt[2];
  int cnt;

  tol = tinfo->wire->comb_tol;
  ExtremeX_Merge(tinfo, &tmin, &tmax);
  ExtremeX_Merge(binfo, &bmin, &bmax);
  if (bmin > tmin)
    tmin = bmin;
  if (bmax < tmax)
    tmax = bmax;

  tseg = &tinfo->wire->seg[1][tinfo->idx];
  bseg = &binfo->wire->seg[0][binfo->idx];
  for (cnt = 0; cnt <= 30; cnt++) {
    val = cnt * (tmax - tmin) / 30 + tmin;
    ttt = TatX(tseg[0].pt, tseg[1].pt, tseg[1].alpha, val);
    btt = TatX(bseg[0].pt, bseg[1].pt, bseg[1].alpha, val);
    EvalSeg(tpt, tseg[0].pt, tseg[1].pt, tseg[1].alpha, ttt);
    EvalSeg(bpt, bseg[0].pt, bseg[1].pt, bseg[1].alpha, btt);
    if (tpt[1] < bpt[1] - tol)
      return 0;
  }
  
  return 1;
}

int VattiVerify(const struct lcw_wire *wire) {
  struct merge_info tinfo, binfo;
  size_t tns, bns;

  if (wire->num_seg[0] < 2 || wire->num_seg[1] < 2)
    return LCW_WIRE_NOT_VATTI_SINGLE;
  
  tns = wire->num_seg[1] - 1;
  bns = wire->num_seg[0] - 1;
  if (fabsf(wire->seg[1][0  ].pt[0] - wire->seg[0][0  ].pt[0]) > wire->comb_tol ||
      fabsf(wire->seg[1][0  ].pt[1] - wire->seg[0][0  ].pt[1]) > wire->comb_tol ||
      fabsf(wire->seg[1][tns].pt[0] - wire->seg[0][bns].pt[0]) > wire->comb_tol ||
      fabsf(wire->seg[1][tns].pt[1] - wire->seg[0][bns].pt[1]) > wire->comb_tol) {
    return LCW_WIRE_NOT_VATTI_SINGLE;
  }
  
  wire->seg[1][0  ].pt[0] = wire->seg[0][0  ].pt[0];
  wire->seg[1][0  ].pt[1] = wire->seg[0][0  ].pt[1];
  wire->seg[1][tns].pt[0] = wire->seg[0][bns].pt[0];
  wire->seg[1][tns].pt[1] = wire->seg[0][bns].pt[1];
  
  memset(&tinfo, 0, sizeof(tinfo));
  memset(&binfo, 0, sizeof(binfo));
  tinfo.wire = (struct lcw_wire *) wire;
  binfo.wire = (struct lcw_wire *) wire;
  tinfo.slot = 1;
  binfo.slot = 0;
  while (tinfo.idx < tns - 1 || binfo.idx < bns - 1) {
    if (!Compatible(&tinfo, &binfo))
      return LCW_WIRE_NOT_VATTI_SINGLE;
    if (tinfo.idx >= tns - 1) {
      binfo.idx++;
    } else if (binfo.idx >= bns - 1) {
      tinfo.idx++;
    } else if (wire->seg[1][tinfo.idx+1].pt[0] < wire->seg[0][binfo.idx+1].pt[0] ||
	       (wire->seg[1][tinfo.idx+1].pt[0] == wire->seg[0][binfo.idx+1].pt[0] &&
		wire->seg[1][tinfo.idx].pt[0] <= wire->seg[0][binfo.idx].pt[0])) {
      tinfo.idx++;
    } else {
      binfo.idx++;
    }
  }
  
  return LCW_NO_ERROR;
}

int LCW_VattiClip(struct lcw_wire *wire) {
  struct lcw_wire *va;
  struct lcw_seg *seg, *pre;
  struct extreme_info ei;
  float pt[2], a1, a2, s1, s2;
  size_t ns, pidx;
  int ret = LCW_NO_ERROR, slot;
  
  if (IS_VATTI(wire))
    return LCW_NO_ERROR;
  
  if (!LCW_IsClosed(wire)) {
    ret = LCW_WIRE_NOT_CLOSED;
    goto err;
  }
  
  if (wire->num_seg[0] < 3) {
    ret = LCW_WIRE_TOO_FEW_SEGMENTS;
    goto err;
  }
  
  if ((va = LCW_New(wire->comb_tol)) == NULL) {
    ret = LCW_OUT_OF_MEMORY;
    goto err;
  }
  
  ns = wire->num_seg[0] - 1;
  FarthestLeft(wire, 0, &ei);
  seg = &wire->seg[0][ei.idx];
  if (ei.tt > T_TOL && ei.tt < 1 - T_TOL) {
    SplitSeg(pt, &a1, &a2, seg[0].pt, seg[1].pt, seg[1].alpha, ei.tt);
    Reset(va, pt, 1);
    slot = seg[0].pt[1] < seg[1].pt[1] ? 0 : 1;
    if ((ret = To(va, seg[0].pt, -a1, slot)) < 0)
      goto err2;
    if ((ret = To(va, seg[1].pt, a2, 1-slot)) < 0)
      goto err2;
    if ((ret = ClipRev(va, wire, (ei.idx + ns - 1) % ns, slot)) < 0)
      goto err2;
    if ((ret = ClipFwd(va, wire, (ei.idx + 1) % ns, 1-slot)) < 0)
      goto err2;
  } else {
    pidx = (ei.idx + ns - 1) % ns;
    pre = &wire->seg[0][pidx];
    s1 = EvalAngle(seg[0].pt, pre[0].pt, -seg[0].alpha, 0);
    s2 = EvalAngle(seg[0].pt, seg[1].pt,  seg[1].alpha, 0);
    Reset(va, seg[0].pt, 1);
    slot = s1 < s2 ? 0 : 1;
    if ((ret = ClipRev(va, wire, pidx, slot)) < 0)
      goto err2;
    if ((ret = ClipFwd(va, wire, ei.idx, 1-slot)) < 0)
      goto err2;
  }
  
  while (va->num_seg[1] > 1 && va->seg[1][va->num_seg[1] - 1].pt[0] == va->seg[1][va->num_seg[1] - 2].pt[0])
    va->num_seg[1]--;
  
  if ((ret = VattiVerify(va)) < 0)
    goto err2;
  
  LCW_Swap(wire, va);
  ret = LCW_NO_ERROR;
  /* Fall through */

 err2:
  LCW_Free(va);
 err:
  return ret;
}

int LCW_WireLoop(struct lcw_wire *wire) {
  struct lcw_wire *nn;
  struct lcw_seg *seg;
  float alpha;
  size_t idx, ns;
  int ret = LCW_NO_ERROR;
  
  if (!IS_VATTI(wire))
    return LCW_NO_ERROR;
  
  if ((nn = LCW_New(wire->comb_tol)) == NULL) {
    ret = LCW_OUT_OF_MEMORY;
    goto err;
  }
  
  ns = wire->num_seg[0] - 1;
  LCW_Reset(nn, wire->seg[0][0].pt);
  for (idx = 0; idx < ns; idx++) {
    seg = &wire->seg[0][idx + 1];
    if ((ret = LCW_To(nn, seg->pt, seg->alpha)) < 0)
      goto err2;
  }
  
  alpha = 0;
  ns = wire->num_seg[1] - 1;
  for (idx = ns + 1; idx > 0; idx--) {
    seg = &wire->seg[1][idx - 1];
    if ((ret = LCW_To(nn, seg->pt, -alpha)) < 0)
      goto err2;
    alpha = seg->alpha;
  }
  
  LCW_Swap(wire, nn);
  ret = LCW_NO_ERROR;
  
 err2:
  LCW_Free(nn);
 err:
  return ret;
}

int VattiInit(struct vatti_state *vs, const struct lcw_wire *wire) {
  int ret;
  
  memset(vs, 0, sizeof(*vs));
  vs->wire = (struct lcw_wire *) wire;
  
  if (!IS_VATTI(wire)) {
    if ((vs->wire = LCW_Copy(wire)) == NULL) {
      ret = LCW_OUT_OF_MEMORY;
      goto err;
    }
    vs->is_copy = 1;
    if ((ret = LCW_VattiClip(vs->wire)) < 0)
      goto err2;
  }
  
  vs->seg[0] = vs->wire->seg[0];
  vs->seg[1] = vs->wire->seg[1];

  return LCW_NO_ERROR;
  
 err2:
  if (vs->is_copy)
    LCW_Free(vs->wire);
 err:
  return ret;
}

int VattiAdv(struct vatti_state *vs, int slot) {
  struct lcw_seg *seg;
  
  if (vs->idx[slot] >= vs->wire->num_seg[slot] - 1)
    return 0;
  
  seg = &vs->wire->seg[slot][++vs->idx[slot]];
  vs->pt[slot] = vs->seg[slot]->pt;
  vs->seg[slot] = seg;
  
  return 1;
}

void VattiDestroy(struct vatti_state *vs) {
  if (vs->is_copy) {
    LCW_Free(vs->wire);
    vs->wire = NULL;
  }
}
