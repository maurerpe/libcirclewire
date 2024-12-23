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

#include "str.h"
#include "util.h"
#include "wire.h"

const char *colors[] = {"black", "red", "green", "blue"};
#define NUM_COLORS (sizeof(colors) / sizeof(*colors))

static int Str_AddM(struct str *str, const float *pt) {
  char buf[256];
  int ret;
  
  if ((ret = snprintf(buf, sizeof(buf), "M%g %g", pt[0], pt[1])) < 0 || ret >= (int)sizeof(buf))
    return LCW_INTERNAL_ERROR;
  
  if ((ret = Str_AddStr(str, buf)) != LCW_NO_ERROR)
    return ret;
  
  return LCW_NO_ERROR;
}

static int Str_AddL(struct str *str, const float *pt) {
  char buf[256];
  int ret;
  
  if ((ret = snprintf(buf, sizeof(buf), " L%g %g", pt[0], pt[1])) < 0 || ret >= (int)sizeof(buf))
    return LCW_INTERNAL_ERROR;
  
  if ((ret = Str_AddStr(str, buf)) != LCW_NO_ERROR)
    return ret;
  
  return LCW_NO_ERROR;
}

static int Str_AddA(struct str *str, const float *prev_pt, const struct lcw_seg *seg) {
  char buf[1024];
  float rad;
  int dir, ret;
  
  rad = Dist(prev_pt, seg->pt) / fabsf(2 * seg->alpha);
  dir = seg->alpha > 0;
  
  if ((ret = snprintf(buf, sizeof(buf), " A%g %g 0 0 %d %g %g", rad, rad, dir, seg->pt[0], seg->pt[1])) < 0 || ret >= (int)sizeof(buf))
    return LCW_INTERNAL_ERROR;
  
  if ((ret = Str_AddStr(str, buf)) != LCW_NO_ERROR)
    return ret;
  
  return LCW_NO_ERROR;
}

static int Str_AddSeg(struct str *str, const struct lcw_seg *seg) {
  if (seg->alpha == 0)
    return Str_AddL(str, seg->pt);
  return Str_AddA(str, seg[-1].pt, seg);
}

static int ToSvgStr(struct str *str, const struct lcw_wire *wire) {
  struct lcw_seg *seg;
  size_t idx;
  int ret;
  
  if (wire->num_seg[0] == 0)
    return LCW_INTERNAL_ERROR;
  
  Str_AddM(str, wire->seg[0][0].pt);
  
  for (idx = 1; idx < wire->num_seg[0]; idx++) {
    seg = &wire->seg[0][idx];
    if ((ret = Str_AddSeg(str, seg)) != LCW_NO_ERROR)
      return ret;
  }
  
  return LCW_NO_ERROR;
}

static int ToSvgPath(char **str_out, const struct lcw_wire *wire) {
  struct str *str;
  struct lcw_wire *cp;
  int ret = LCW_NO_ERROR;
  
  if ((str = Str_New()) == NULL) {
    ret = LCW_OUT_OF_MEMORY;
    goto err;
  }
  
  if (!IS_VATTI(wire)) {
    if ((ret = ToSvgStr(str, wire)) != LCW_NO_ERROR)
      goto err2;
  } else {
    if ((cp = LCW_Copy(wire)) == NULL) {
      ret = LCW_OUT_OF_MEMORY;
      goto err2;
    }
    
    if ((ret = LCW_WireLoop(cp)) != LCW_NO_ERROR)
      goto err3;
  
    if ((ret = ToSvgStr(str, cp)) != LCW_NO_ERROR)
      goto err3;
    
    LCW_Free(cp);
  }
  
  *str_out = Str_FreeReturnStr(str);
  return LCW_NO_ERROR;

 err3:
  LCW_Free(cp);
 err2:
  Str_Free(str);
 err:
  return ret;
}

char *LCW_ToSvgPath(const struct lcw_wire *wire) {
  char *path;
  
  if (ToSvgPath(&path, wire) != LCW_NO_ERROR)
    return NULL;
  
  return path;
}

void LCW_PathFree(char *path) {
  free(path);
}

int LCW_ToSvg(const char *filename, const struct lcw_wire *wire, const char *unit, int flip_y) {
  struct lcw_list list;

  memset(&list, 0, sizeof(list));
  list.wire = (struct lcw_wire *) wire;

  return LCW_ListToSvg(filename, &list, unit, flip_y);
}

int LCW_ListToSvg(const char *filename, const struct lcw_list *list, const char *unit, int flip_y) {
  const struct lcw_list *cur;
  struct lcw_properties prop;
  char **path, **pp;
  FILE *out;
  float min[2], max[2], tmp;
  size_t num, pidx;
  int ret = LCW_NO_ERROR, idx;
  
  min[0] = INFINITY;
  min[1] = INFINITY;
  max[0] = -INFINITY;
  max[1] = -INFINITY;
  num = 0;
  for (cur = list; cur != NULL; cur = cur->next)
    num++;
  if ((path = calloc(num, sizeof(*path))) == NULL) {
    ret = LCW_OUT_OF_MEMORY;
    goto err;
  }
  pp = path;
  for (cur = list; cur != NULL; cur = cur->next, pp++) {
    if ((ret = ToSvgPath(pp, cur->wire)) != LCW_NO_ERROR)
      goto err2;
    
    if ((ret = LCW_Properties(&prop, cur->wire)) != LCW_NO_ERROR)
      goto err2;
    
    for (idx = 0; idx < 2; idx++) {
      if (prop.min[idx] < min[idx])
	min[idx] = prop.min[idx];
      if (prop.max[idx] > max[idx])
	max[idx] = prop.max[idx];
    }
  }
  if (list == NULL) {
    min[0] = 0;
    min[1] = 0;
    max[0] = 1;
    max[1] = 1;
  } else if (flip_y) {
    tmp = max[1];
    max[1] = -min[1];
    min[1] = -tmp;
  }
  if ((out = fopen(filename, "w")) == NULL) {
    ret = LCW_COULD_NOT_OPEN_FILE;
    goto err2;
  }
  
  if (fprintf(out, "<?xml version=\"1.0\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n<svg width=\"%g%s\" height=\"%g%s\" viewBox=\"%g %g %g %g\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n<g%s>\n",
	      max[0] - min[0],
	      unit,
	      max[1] - min[1],
	      unit,
	      min[0],
	      min[1],
	      max[0] - min[0],
	      max[1] - min[1],
	      flip_y ? " transform=\"scale(1,-1)\"" : "") < 0) {
    ret = LCW_COULD_NOT_WRITE_FILE;
    goto err3;
  }
  
  for (pidx = 0; pidx < num; pidx++, pp++) {
    if (fprintf(out, "  <path d=\"%s\" stroke=\"%s\" stroke-width=\"0.01\" fill=\"none\"/>\n", path[pidx], colors[pidx % NUM_COLORS]) < 0) {
      ret = LCW_COULD_NOT_WRITE_FILE;
      goto err3;
    }
  }

  if (fprintf(out, "</g>\n</svg>") < 0) {
    ret = LCW_COULD_NOT_WRITE_FILE;
    goto err3;
  }
  
  fclose(out);
  for (idx = 0; idx < num; idx++)
    LCW_PathFree(path[idx]);
  free(path);
  return LCW_NO_ERROR;
  
 err3:
  fclose(out);
 err2:
  for (idx = 0; idx < num; idx++)
    LCW_PathFree(path[idx]);
  free(path);
 err:
  return ret;
}
