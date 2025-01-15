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

#include <limits.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "wire.h"

int EnsureAlloc(struct lcw_wire *wire, size_t num, int slot) {
  struct lcw_seg *new_seg;
  size_t num_alloc = 1;
  
  num += wire->num_seg[slot];
  if (num <= wire->num_alloc[slot])
    return LCW_NO_ERROR;
  
  while (num_alloc < num)
    num_alloc <<= 1;
  
  if ((new_seg = calloc(num_alloc, sizeof(*new_seg))) == NULL)
    return LCW_OUT_OF_MEMORY;
  
  if (wire->num_seg[slot] > 0)
    memcpy(new_seg, wire->seg[slot], wire->num_seg[slot] * sizeof(*new_seg));
  free(wire->seg[slot]);
  wire->seg[slot] = new_seg;
  wire->num_alloc[slot] = num_alloc;
  
  return LCW_NO_ERROR;
}

struct lcw_wire *LCW_New(float seg_combine_tol) {
  struct lcw_wire *wire;

  if ((wire = malloc(sizeof(*wire))) == NULL)
    return NULL;
  memset(wire, 0, sizeof(*wire));
  wire->comb_tol = seg_combine_tol;
  
  EnsureAlloc(wire, 16, 0);
  EnsureAlloc(wire, 16, 1);
  wire->num_seg[0] = 1;
  
  return wire;
}

void LCW_Free(struct lcw_wire *wire) {
  if (wire == NULL)
    return;

  free(wire->seg[0]);
  free(wire->seg[1]);
  free(wire);
}

struct lcw_wire *LCW_Copy(const struct lcw_wire *wire) {
  struct lcw_wire *copy;
  int slot;

  if ((copy = LCW_New(wire->comb_tol)) == NULL)
    goto err;

  for (slot = 0; slot < 2; slot++) {
    if (EnsureAlloc(copy, wire->num_seg[slot], slot) < 0)
      goto err2;
    
    memcpy(copy->seg[slot], wire->seg[slot], wire->num_seg[slot] * sizeof(*copy->seg[slot]));
    copy->num_seg[slot] = wire->num_seg[slot];
  }
  
  return copy;
  
 err2:
  LCW_Free(copy);
 err:
  return NULL;
}

void LCW_Swap(struct lcw_wire *a, struct lcw_wire *b) {
  struct lcw_wire tmp;

  memcpy(&tmp, a, sizeof(tmp));
  memcpy(a,    b, sizeof(*a));
  memcpy(b, &tmp, sizeof(*b));
}

void Reset(struct lcw_wire *wire, const float *start, int vatti) {
  struct lcw_seg *seg;
  
  seg = wire->seg[0];
  memset(seg, 0, sizeof(*seg));
  seg->pt[0] = start[0];
  seg->pt[1] = start[1];
  wire->num_seg[0] = 1;

  if (vatti) {
    seg = wire->seg[1];
    memset(seg, 0, sizeof(*seg));
    seg->pt[0] = start[0];
    seg->pt[1] = start[1];
    wire->num_seg[1] = 1;
  } else {
    wire->num_seg[1] = 0;
  }
}

void LCW_Reset(struct lcw_wire *wire, const float *start) {
  Reset(wire, start, 0);
}

static int CanCombineStraight(const float *a, const float *b, const float *c, float *new_alpha, float tol) {
  float ba[2], cb[2], mag;

  if (b[0] == a[0] && b[1] == a[1])
    return 0;
  
  ba[0] = b[0] - a[0];
  ba[1] = b[1] - a[1];
  
  cb[0] = c[0] - b[0];
  cb[1] = c[1] - b[1];
  
  Normalize(ba);
  mag = Dot(ba, cb);
  if (mag < 0)
    return 0;
  ba[0] *= mag;
  ba[1] *= mag;
  if (Dist2(ba, cb) > tol * tol)
    return 0;

  if (new_alpha)
    *new_alpha = 0;
  
  return 1;
}

static int CanCombine(const float *src, const float *mid, float malpha, const float *dst, float dalpha, float tol, float *new_alpha) {
  float a1, a2, ang, na, nm[2];
  
  if (tol == 0)
    return 0;
  
  if (malpha == 0 && dalpha == 0)
    return CanCombineStraight(src, mid, dst, new_alpha, tol);
  
  if (malpha * dalpha <= 0)
    return 0;
  
  a1 = 2 * asinf(malpha);
  a2 = 2 * asinf(dalpha);
  ang = a1 + a2;
  if (fabsf(ang) > M_PI)
    return 0;
  
  na = sinf(ang / 2);
  EvalSeg(nm, src, dst, na, a1 / ang);
  if (Dist2(mid, nm) > tol * tol)
    return 0;
  
  if (new_alpha)
    *new_alpha = na;
  
  return 1;
}

int To(struct lcw_wire *wire, const float *dst, float alpha, int slot) {
  struct lcw_seg *seg;
  float *src, new_alpha, dist;
  int ret;
  
  if (slot == 1 && !IS_VATTI(wire))
    return LCW_INTERNAL_ERROR;
  
  if ((ret = EnsureAlloc(wire, 1, slot)) < 0)
    return ret;
  
  seg = &wire->seg[slot][wire->num_seg[slot]];
  src = seg[-1].pt;
  if (src[0] == dst[0] && src[1] == dst[1])
    return LCW_NO_ERROR;
  if (fabsf(alpha) > 1) {
    dist = Dist(src, dst);
    if ((fabsf(alpha) - 1) * dist > wire->comb_tol)
      return LCW_RADIUS_TO_SMALL;
    alpha = 1;
  }
  if (wire->num_seg[slot] >= 2 && CanCombine(seg[-2].pt, src, seg[-1].alpha, dst, alpha, wire->comb_tol, &new_alpha)) {
    seg[-1].alpha = new_alpha;
    src[0] = dst[0];
    src[1] = dst[1];
    return LCW_NO_ERROR;
  }
  memset(seg, 0, sizeof(*seg));
  seg->pt[0] = dst[0];
  seg->pt[1] = dst[1];
  seg->alpha = alpha;
  wire->num_seg[slot]++;

  return LCW_NO_ERROR;
}

int LCW_To(struct lcw_wire *wire, const float *dst, float alpha) {
  int ret;
  
  if (IS_VATTI(wire) && (ret = LCW_WireLoop(wire)) < 0)
    return ret;
  
  return To(wire, dst, alpha, 0);
}

int LCW_LineTo(struct lcw_wire *wire, const float *dst) {
  return LCW_To(wire, dst, 0);
}

int LCW_CircleRadiusTo(struct lcw_wire *wire, const float *dst, float radius, int large, int cw) {
  float *src, delta[2], dist, hd, perp[2], scale, mid[2];
  int ret;
  
  if (IS_VATTI(wire) && (ret = LCW_WireLoop(wire)) < 0)
    return ret;
  
  src = wire->seg[0][wire->num_seg[0]-1].pt;
  delta[0] = dst[0] - src[0];
  delta[1] = dst[1] - src[1];
  dist = sqrtf(delta[0] * delta[0] + delta[1] * delta[1]);
  if (dist == 0)
    return LCW_START_AND_END_EQUAL;
  hd = dist / 2;
  if (hd > radius + wire->comb_tol)
    return LCW_RADIUS_TO_SMALL;
  if (cw)
    dist = -dist;
  
  if (!large)
    return LCW_To(wire, dst, 0.5f * dist / radius);
  
  perp[0] = -delta[1];
  perp[1] =  delta[0];
  
  scale = (sqrtf(radius * radius - hd * hd) + radius) / dist;
  mid[0] = 0.5 * (src[0] + dst[0]) - perp[0] * scale;
  mid[1] = 0.5 * (src[1] + dst[1]) - perp[1] * scale;
  
  if ((ret = LCW_CircleRadiusTo(wire, mid, radius, 0, cw)) < 0)
    return ret;
  
  return LCW_CircleRadiusTo(wire, dst, radius, 0, cw);
}

int LCW_CircleCenterTo(struct lcw_wire *wire, const float *dst, const float *center, int cw) {
  float *src, delta[2], dist, hd, perp[2], leg, radius;
  int ret;
  
  if (IS_VATTI(wire) && (ret = LCW_WireLoop(wire)) < 0)
    return ret;
  
  src = wire->seg[0][wire->num_seg[0]-1].pt;
  delta[0] = dst[0] - src[0];
  delta[1] = dst[1] - src[1];
  dist = sqrtf(delta[0] * delta[0] + delta[1] * delta[1]);
  if (dist == 0)
    return LCW_START_AND_END_EQUAL;
  hd = dist / 2;
  if (cw)
    dist = -dist;
  
  perp[0] = -delta[1];
  perp[1] =  delta[0];
  
  leg = (perp[0] * (center[0] - src[0]) + perp[1] * (center[1] - src[1])) / dist;
  radius = sqrtf(leg * leg + hd * hd);

  return LCW_CircleRadiusTo(wire, dst, radius, leg > 0, cw);
}

int LCW_Close(struct lcw_wire *wire) {
  if (LCW_IsClosed(wire))
    return LCW_NO_ERROR;
  
  return LCW_LineTo(wire, wire->seg[0][0].pt);
}

void LCW_Translate(struct lcw_wire *wire, float dx, float dy) {
  struct lcw_seg *seg;
  size_t idx, slot;

  for (slot = 0; slot < 2; slot++) {
    seg = wire->seg[slot];
    for (idx = 0; idx < wire->num_seg[slot]; idx++, seg++) {
      seg->pt[0] += dx;
      seg->pt[1] += dy;
    }
  }
}

void LCW_Scale(struct lcw_wire *wire, float cx, float cy, float scale) {
  struct lcw_seg *seg;
  size_t idx, slot;

  for (slot = 0; slot < 2; slot++) {
    seg = wire->seg[slot];
    for (idx = 0; idx < wire->num_seg[slot]; idx++, seg++) {
      seg->pt[0] = (seg->pt[0] - cx) * scale + cx;
      seg->pt[1] = (seg->pt[1] - cy) * scale + cy;
    }
  }
}

int LCW_Rotate(struct lcw_wire *wire, float cx, float cy, float ccw_rad) {
  struct lcw_seg *seg;
  size_t idx;
  float cosa, sina, x, y;
  int vatti, ret;

  vatti = IS_VATTI(wire);
  if (vatti && (ret = LCW_WireLoop(wire)) < 0)
    return ret;
  
  cosa = cosf(ccw_rad);
  sina = sinf(ccw_rad);
  seg = wire->seg[0];
  for (idx = 0; idx < wire->num_seg[0]; idx++, seg++) {
    x = seg->pt[0] - cx;
    y = seg->pt[1] - cy;
    seg->pt[0] = cosa * x - sina * y + cx;
    seg->pt[1] = cosa * y + sina * x + cy;
  }

  if (vatti && (ret = LCW_VattiClip(wire)) < 0)
    return ret;
  
  return LCW_NO_ERROR;
}

int LCW_MirrorX(struct lcw_wire *wire) {
  struct lcw_seg *seg;
  size_t idx;
  int vatti, ret;

  vatti = IS_VATTI(wire);
  if (vatti && (ret = LCW_WireLoop(wire)) < 0)
    return ret;
  
  seg = wire->seg[0];
  for (idx = 0; idx < wire->num_seg[0]; idx++, seg++) {
    seg->pt[0] = -seg->pt[0];
    seg->alpha = -seg->alpha;
  }

  if (vatti && (ret = LCW_VattiClip(wire)) < 0)
    return ret;
  
  return LCW_NO_ERROR;
}

void LCW_MirrorY(struct lcw_wire *wire) {
  struct lcw_seg *seg;
  size_t idx, slot, tmp;
  
  for (slot = 0; slot < 2; slot++) {
    seg = wire->seg[slot];
    for (idx = 0; idx < wire->num_seg[slot]; idx++, seg++) {
      seg->pt[1] = -seg->pt[1];
      seg->alpha = -seg->alpha;
    }
  }
  
  if (wire->num_seg[1] > 0) {
    seg = wire->seg[0];
    wire->seg[0] = wire->seg[1];
    wire->seg[1] = seg;
    tmp = wire->num_seg[0];
    wire->num_seg[0] = wire->num_seg[1];
    wire->num_seg[1] = tmp;
    tmp = wire->num_alloc[0];
    wire->num_alloc[0] = wire->num_alloc[1];
    wire->num_alloc[1] = tmp;
  }
}

void MirrorYTop(struct lcw_wire *wire) {
  struct lcw_seg *seg;
  size_t idx;
  
  seg = wire->seg[1];
  for (idx = 0; idx < wire->num_seg[1]; idx++, seg++) {
    seg->pt[1] = -seg->pt[1];
    seg->alpha = -seg->alpha;
  }
}

int LCW_IsClosed(const struct lcw_wire *wire) {
  struct lcw_seg *start, *stop;
  
  if (IS_VATTI(wire))
    return 1;
  
  start = &wire->seg[0][0];
  stop  = &wire->seg[0][wire->num_seg[0] - 1];
  return start->pt[0] == stop->pt[0] && start->pt[1] == stop->pt[1];
}

int LCW_IsVatti(const struct lcw_wire *wire) {
  return IS_VATTI(wire);
}

static int ToPoly(struct lcw_wire *dst, const struct lcw_wire *src, float tol) {
  struct lcw_seg *seg;
  int ret;
  size_t idx, slot;
  ssize_t num, cnt;
  float dist, alpha, t, pt[2];
  
  Reset(dst, src->seg[0][0].pt, IS_VATTI(src));
  for (slot = 0; slot < 2; slot++) {
    seg = src->seg[slot];
    for (idx = 0; idx + 1 < src->num_seg[slot]; idx++, seg++) {
      dist = Dist(seg[0].pt, seg[1].pt);
      alpha = seg[1].alpha;
      num = RequiredSegs(dist, alpha, tol);
      for (cnt = 1; cnt < num; cnt++) {
	t = ((float) cnt) / num;
	EvalSeg(pt, seg[0].pt, seg[1].pt, alpha, t);
	if ((ret = To(dst, pt, 0, slot)) < 0)
	  return ret;
      }
      
      if ((ret = To(dst, seg[1].pt, 0, slot)) < 0)
	return ret;
    }
  }
  
  return LCW_NO_ERROR;
}

struct lcw_wire *LCW_ToPolygon(const struct lcw_wire *wire, float tol) {
  struct lcw_wire *poly;
  
  if ((poly = LCW_New(wire->comb_tol)) == NULL)
    goto err;
  
  if (ToPoly(poly, wire, tol) != LCW_NO_ERROR)
    goto err2;
  
  return poly;

 err2:
  LCW_Free(poly);
 err:
  return NULL;
}

float LCW_TotalArcLen(const struct lcw_wire *wire) {
  size_t slot, idx;
  float tot = 0;

  for (slot = 0; slot < 2; slot++)
    for (idx = 1; idx < wire->num_seg[slot]; idx++)
      tot += ArcLen(wire->seg[slot][idx-1].pt, wire->seg[slot][idx].pt, wire->seg[slot][idx].alpha);
  
  return tot;
}

struct lp_vertex_list *LCW_Mesh(const struct lcw_wire *wire, float tol) {
  struct lcw_wire *poly;
  struct lp_vertex_list *vl, *vl2;
  size_t slot, idx;

  if ((poly = LCW_ToPolygon(wire, tol)) == NULL)
    goto err;
  
  if ((vl = LP_VertexList_New(2, lp_pt_line)) == NULL)
    goto err2;
  
  for (slot = 0; slot < 2; slot++) {
    for (idx = 1; idx < poly->num_seg[slot]; idx++) {
      if (LP_VertexList_Add(vl, poly->seg[slot][idx - 1].pt) == UINT_MAX)
	goto err3;
      if (LP_VertexList_Add(vl, poly->seg[slot][idx].pt) == UINT_MAX)
	goto err3;
    }
  }
  
  if ((vl2 = LP_Triangulate2D(vl)) == NULL)
    goto err3;
  
  LP_VertexList_Free(vl);
  LCW_Free(poly);
  return vl2;

 err3:
  LP_VertexList_Free(vl);
 err2:
  LCW_Free(poly);
 err:
  return NULL;
}

size_t LCW_NumSegs(const struct lcw_wire *wire, int slot) {
  if (slot < 0 || slot >= 2)
    return 0;
  
  return wire->num_seg[slot];
}

const struct lcw_seg *LCW_Segs(const struct lcw_wire *wire, int slot) {
  if (slot < 0 || slot >= 2)
    return NULL;
  
  return wire->seg[slot];
}

struct lcw_list *LCW_ListNew(void) {
  struct lcw_list *list;
  
  if ((list = malloc(sizeof(*list))) == NULL)
    goto err;
  memset(list, 0, sizeof(*list));

  return list;

 err:
  return NULL;
}

void LCW_ListFree(struct lcw_list *list, int free_wires) {
  struct lcw_list *next;
  
  while (list != NULL) {
    next = list->next;
    if (free_wires)
      LCW_Free(list->wire);
    free(list);
    list = next;
  }
}
