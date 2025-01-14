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
#include "wire.h"

static void Check_Seg(float *rad, const float *aa, const float *bb, float alpha, const float *center) {
  float tt, pt[2], dd;
  int cnt;
  
  for (cnt = 1; cnt < 100; cnt++) {
    tt = cnt / 100.0f;
    EvalSeg(pt, aa, bb, alpha, tt);
    dd = Dist(pt, center);
    if (dd > *rad)
      *rad = dd;
  }
}

static float Calc_Bound(const struct lcw_wire *wire, const float *center) {
  struct lcw_seg *seg;
  size_t cnt, ns;
  float rad, dd;
  
  seg = wire->seg[0];
  ns = wire->num_seg[0] - 1;
  rad = Dist(seg->pt, center);
  for (cnt = 0; cnt < ns; cnt++, seg++) {
    dd = Dist(seg[1].pt, center);
    if (dd > rad)
      rad = dd;
    if (seg[1].alpha != 0)
      Check_Seg(&rad, seg[0].pt, seg[1].pt, seg[1].alpha, center);
  }

  return rad;
}

float LCW_BoundingCircle(const struct lcw_wire *wire, const float *center) {
  struct lcw_wire *cp;
  float ret;
  
  if (!IS_VATTI(wire))
    return Calc_Bound(wire, center);

  if ((cp = LCW_Copy(wire)) == NULL)
    goto err;
  if (LCW_WireLoop(cp) != LCW_NO_ERROR)
    goto err2;
  ret = Calc_Bound(cp, center);
  LCW_Free(cp);
  return ret;

 err2:
  LCW_Free(cp);
 err:
  return INFINITY;
}
