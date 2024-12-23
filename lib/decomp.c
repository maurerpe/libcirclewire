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

static void TransPt(float *dst, const float *pt, const float *base, const float *rot) {
  float delta[2];

  delta[0] = pt[0] - base[0];
  delta[1] = pt[1] - base[1];
  
  dst[0] = delta[0] * rot[0] + delta[1] * rot[1];
  dst[1] = delta[1] * rot[0] - delta[0] * rot[1];
}

static float PtError(const float *pt, float len) {
  float val[2];
  
  if (pt[0] < 0) {
    return Norm(pt);
  } else if (pt[0] > len) {
    val[0] = pt[0] - len;
    val[1] = pt[1];
    return Norm(val);
  }
  
  return fabsf(pt[1]);
}

static void CalcSegMaxError(float *maxerr, float *xloc, const float *a1, const float *b1, float alpha1, const float *a2, const float *b2) {
  float rot[2], len, aa[2], bb[2], err, cc[2], rad, pt[2], dir[2], tt;
  
  rot[0] = b2[0] - a2[0];
  rot[1] = b2[1] - a2[1];
  len = Norm(rot);
  if (len == 0)
    return;
  rot[0] /= len;
  rot[1] /= len;
  
  TransPt(aa, a1, a2, rot);
  TransPt(bb, b1, a2, rot);

  err = PtError(aa, len);
  if (err > *maxerr) {
    *maxerr = err;
    *xloc = a1[0];
  }

  err = PtError(bb, len);
  if (err > *maxerr) {
    *maxerr = err;
    *xloc = b1[0];
  }
  
  if (alpha1 < 0) {
    ArcCenter(cc, aa, bb, alpha1);
    rad = Dist(aa, bb) / (2 * fabsf(alpha1));
    dir[0] = 0;
    dir[1] = cc[1];
    if (cc[0] < 0) {
      dir[0] = cc[0];
    } else if (cc[0] > len) {
      dir[0] = cc[0] - len;
    }
    Normalize(dir);
    tt = TatX(aa, bb, alpha1, cc[0] + dir[0] * rad);
    if (tt > 0 && tt < 1) {
      EvalSeg(pt, aa, bb, alpha1, tt);
      err = PtError(pt, len);
      if (err > *maxerr) {
	*maxerr = err;
	EvalSeg(pt, a1, b1, alpha1, tt);
	*xloc = pt[0];
      }
    }
  }
}

static void SegMaxError(float *maxerr, float *xloc, const float *aa, const float *bb, float alpha, const struct lcw_wire *hull, size_t *idx, int slot) {
  struct lcw_seg *seg;
  float tt, split[2];
  
  while (bb[0] > hull->seg[slot][*idx].pt[0] + hull->comb_tol) {
    if (*idx >= hull->num_seg[slot] - 1)
      return;
    (*idx)++;
  }

  seg = &hull->seg[slot][*idx];
  if (seg->alpha != 0)
    return;
  
  tt = TatX(aa, bb, alpha, seg[-1].pt[0]);
  if (tt > 0) {
    if (tt >= 1)
      return;
    SplitSeg(split, NULL, &alpha, aa, bb, alpha, tt);
    aa = split;
  }
  
  if (fabsf(aa[0] - bb[0]) < hull->comb_tol)
    return;
  
  CalcSegMaxError(maxerr, xloc, aa, bb, alpha, seg[-1].pt, seg->pt);
}

static void CalcMaxError(float *maxerr, float *xloc, struct vatti_state *vs, struct lcw_wire *hull) {
  float aa[2], bb[2];
  size_t idx, slot;

  for (slot = 0; slot < 2; slot++) {
    idx = 1;
    while (VattiAdv(vs, slot)) {
      if (slot > 0) {
	aa[0] =  vs->pt[1][0];
	aa[1] = -vs->pt[1][1];
	bb[0] =  vs->seg[1]->pt[0];
	bb[1] = -vs->seg[1]->pt[1];
	SegMaxError(maxerr, xloc,
		    aa, bb, -vs->seg[1]->alpha,
		    hull, &idx, slot);
      } else {
	SegMaxError(maxerr, xloc,
		    vs->pt[0], vs->seg[0]->pt, vs->seg[0]->alpha,
		    hull, &idx, slot);
      }
    }
  }
}

static int Decomp(struct lcw_list *list, float tol) {
  struct lcw_wire *hull, *cp;
  struct vatti_state vs;
  struct lcw_list *cut;
  float maxerr = 0, xloc = 0;
  
  if ((hull = LCW_ConvexHull(list->wire)) == NULL)
    goto err;
  
  if (LCW_VattiClip(hull) != LCW_NO_ERROR)
    goto err2;
  
  if ((cp = LCW_Copy(hull)) == NULL)
    goto err2;
  
  if (VattiInit(&vs, list->wire) < 0)
    goto err3;
  
  MirrorYTop(cp);
  CalcMaxError(&maxerr, &xloc, &vs, cp);
  
  if (maxerr > tol) {
    if (xloc >= list->wire->seg[0][list->wire->num_seg[0]-1].pt[0] - list->wire->comb_tol)
      xloc -= tol / 2;
    if ((cut = LCW_CutAtX(list->wire, xloc)) == NULL)
      goto err4;
    
    if (cut->next == NULL)
      goto err5;
    
    if (LCW_VattiClip(cut->wire) != LCW_NO_ERROR)
      goto err5;
    
    if (LCW_VattiClip(cut->next->wire) != LCW_NO_ERROR)
      goto err5;
    
    LCW_Swap(list->wire, cut->wire);
    cut->next->next = list->next;
    list->next = cut->next;
    cut->next = NULL;
    LCW_ListFree(cut, 1);
    
    if (Decomp(list->next, tol) < 0)
      goto err5;
    
    if (Decomp(list, tol) < 0)
      goto err5;
  } else {
    LCW_Swap(list->wire, hull);
  }
  
  VattiDestroy(&vs);
  LCW_Free(cp);
  LCW_Free(hull);
  return 0;
  
 err5:
  LCW_ListFree(cut, 1);
 err4:
  VattiDestroy(&vs);
 err3:
  LCW_Free(cp);
 err2:
  LCW_Free(hull);
 err:
  return -1;
}

struct lcw_list *LCW_ConvexDecomp(const struct lcw_wire *wire, float tol) {
  struct lcw_list *head;
  
  if ((head = LCW_ListNew()) == NULL)
    goto err;
  
  if ((head->wire = LCW_Copy(wire)) == NULL)
    goto err2;
  
  if (LCW_VattiClip(head->wire) != LCW_NO_ERROR)
    goto err2;
  
  if (Decomp(head, tol) < 0)
    goto err2;
  
  return head;

 err2:
  LCW_ListFree(head, 1);
 err:
  return NULL;
}
