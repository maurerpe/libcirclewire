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

int Cut(struct lcw_wire *out1, struct lcw_wire *out2, struct vatti_state *vs, float xval) {
  float tt, pt[2], a1, a2;
  int ret, slot;

  Reset(out1, vs->seg[0]->pt, 1);
  for (slot = 0; slot < 2; slot++) {
    if (!VattiAdv(vs, slot))
      return LCW_INTERNAL_ERROR;
    while (vs->seg[slot]->pt[0] < xval) {
      if ((ret = To(out1, vs->seg[slot]->pt, vs->seg[slot]->alpha, slot)) < 0)
	return ret;
      
      if (!VattiAdv(vs, slot))
	return LCW_INTERNAL_ERROR;
    }
  }

  for (slot = 0; slot < 2; slot++) {
    tt = TatX(vs->pt[slot], vs->seg[slot]->pt, vs->seg[slot]->alpha, xval);
    SplitSeg(pt, &a1, &a2, vs->pt[slot], vs->seg[slot]->pt, vs->seg[slot]->alpha, tt);
    if ((ret = To(out1, pt, a1, slot)) < 0)
      return ret;
    if (slot == 1 && (ret = To(out1, pt, 0, 0)) < 0)
      return ret;
  }
  
  for (slot = 0; slot < 2; slot++) {
    while (vs->seg[slot]->pt[0] <= xval) {
      if (!VattiAdv(vs, slot))
	return LCW_INTERNAL_ERROR;
    }
  }
  
  for (slot = 0; slot < 2; slot++) {
    tt = TatX(vs->pt[slot], vs->seg[slot]->pt, vs->seg[slot]->alpha, xval);
    SplitSeg(pt, &a1, &a2, vs->pt[slot], vs->seg[slot]->pt, vs->seg[slot]->alpha, tt);
    if (slot == 0) {
      Reset(out2, pt, 1);
    } else if ((ret = To(out2, pt, 0, 1)) < 0) {
      return ret;
    }
    if ((ret = To(out2, vs->seg[slot]->pt, a2, slot)) < 0)
      return ret;
  }
  
  for (slot = 0; slot < 2; slot++) {
    while (VattiAdv(vs, slot)) {
      if ((ret = To(out2, vs->seg[slot]->pt, vs->seg[slot]->alpha, slot)) < 0)
	return ret;
    }
  }
  
  return LCW_NO_ERROR;
}

struct lcw_list *LCW_CutAtX(const struct lcw_wire *wire, float xval) {
  struct vatti_state vs;
  struct lcw_list *ret = NULL;
  struct lcw_wire *ww[2];
  int idx;

  if (VattiInit(&vs, wire) < 0)
    goto err;
  
  if ((ret = LCW_ListNew()) == NULL)
    goto err2;
  
  if (xval <= vs.wire->seg[0][0].pt[0] + wire->comb_tol || xval >= vs.wire->seg[0][vs.wire->num_seg[0]-1].pt[0] - wire->comb_tol) {
    if ((ret->wire = LCW_Copy(vs.wire)) == NULL)
      goto err3;
    return ret;
  }
  
  if ((ret->next = LCW_ListNew()) == NULL)
    goto err3;
  
  for (idx = 0; idx < 2; idx++)
    if ((ww[idx] = LCW_New(vs.wire->comb_tol)) == NULL)
      goto err4;
  
  if (Cut(ww[0], ww[1], &vs, xval) < 0)
    goto err4;
  
  VattiDestroy(&vs);
  ret->wire = ww[0];
  ret->next->wire = ww[1];
  return ret;
  
 err4:
  while (--idx >= 0)
    LCW_Free(ww[idx]);
 err3:
  LCW_ListFree(ret, 1);
 err2:
  VattiDestroy(&vs);
 err:
  return NULL;
}
