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

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "wire.h"

static int ParseFloat(const char **str, float *val) {
  char *endptr;
  float vv;

  vv = strtof(*str, &endptr);
  
  if (endptr == NULL || (vv == HUGE_VALF && errno == ERANGE))
    return LCW_PARSE_ERROR;

  *str = endptr;
  *val = vv;
  return LCW_NO_ERROR;
}

static int ParseSep(const char **str) {
  const char *ss = *str;
  int had_comma = 0;
  
  while (isspace(*ss) || (!had_comma && *ss == ',')) {
    if (*ss == ',')
      had_comma = 1;
    ss++;
  }

  if (ss == *str)
    return LCW_PARSE_ERROR;
  
  *str = ss;
  return LCW_NO_ERROR;
}

static void ParseWhite(const char **str) {
  const char *ss = *str;

  while (isspace(*ss))
    ss++;
  
  *str = ss;
}

static int ParseFloats(const char **str, float *val, size_t num) {
  size_t cnt;
  int ret;
  
  ParseWhite(str);
  for (cnt = 0; cnt < num; cnt++) {
    if (cnt != 0 && (ret = ParseSep(str)) != LCW_NO_ERROR)
      return ret;
    if ((ret = ParseFloat(str, &val[cnt])) != LCW_NO_ERROR)
      return ret;
  }

  return LCW_NO_ERROR;
}

static void UpdatePos(float *cur, const float *val, char cmd) {
  if (cmd >= 'a' && cmd <= 'z') {
    cur[0] += val[0];
    cur[1] += val[1];
  } else {
    cur[0] = val[0];
    cur[1] = val[1];
  }
}

static int ParsePath(struct lcw_wire *wire, const char *str) {
  char cmd;
  float start[2], cur[2], val[7];
  int ret, closed = 0;
  
  ParseWhite(&str);
  cmd = *str++;
  if (cmd != 'M' && cmd != 'm')
    return LCW_PARSE_ERROR;
  if ((ret = ParseFloats(&str, val, 2)) != LCW_NO_ERROR)
    return ret;
  start[0] = cur[0] = val[0];
  start[1] = cur[1] = val[1];
  LCW_Reset(wire, cur);
  
  while (*str != '\0') {
    if ((ret = ParseSep(&str)) != LCW_NO_ERROR)
      return ret;
    if (*str == '\0')
      break;
    if (closed)
      return LCW_PARSE_ERROR;
    switch ((cmd = *str++)) {
    case 'L':
    case 'l':
      if ((ret = ParseFloats(&str, val, 2)) != LCW_NO_ERROR)
	return ret;
      UpdatePos(cur, val, cmd);
      if ((ret = LCW_LineTo(wire, cur)) != LCW_NO_ERROR)
	return ret;
      break;

    case 'H':
    case 'h':
      if ((ret = ParseFloats(&str, val, 1)) != LCW_NO_ERROR)
	return ret;
      if (cmd == 'h')
	cur[0] += val[0];
      else
	cur[0] = val[0];
      if ((ret = LCW_LineTo(wire, cur)) != LCW_NO_ERROR)
	return ret;
      break;

    case 'V':
    case 'v':
      if ((ret = ParseFloats(&str, val, 1)) != LCW_NO_ERROR)
	return ret;
      if (cmd == 'v')
	cur[1] += val[0];
      else
	cur[1] = val[0];
      if ((ret = LCW_LineTo(wire, cur)) != LCW_NO_ERROR)
	return ret;
      break;

    case 'A':
    case 'a':
      if ((ret = ParseFloats(&str, val, 7)) != LCW_NO_ERROR)
	return ret;
      if (val[0] != val[1])
	return LCW_ELLIPSE_NOT_SUPPORTED;
      UpdatePos(cur, &val[5], cmd);
      if ((ret = LCW_CircleRadiusTo(wire, cur, val[0], val[3] != 0, val[4] == 0)) != LCW_NO_ERROR)
	return ret;
      break;

    case 'Z':
    case 'z':
      closed = 1;
      if ((ret = LCW_LineTo(wire, start)) != LCW_NO_ERROR)
	return ret;
      break;

    default:
      return LCW_PARSE_ERROR;
    }
  }

  return LCW_NO_ERROR;
}

struct lcw_wire *LCW_FromSvgPath(const char *d_str, int flip_y, float seg_combine_tol) {
  struct lcw_wire *wire;
  
  if ((wire = LCW_New(seg_combine_tol)) == NULL)
    goto err;
  
  if (ParsePath(wire, d_str) != LCW_NO_ERROR)
    goto err2;
  
  if (flip_y)
    LCW_MirrorY(wire);
  
  return wire;
  
 err2:
  LCW_Free(wire);
 err:
  return NULL;
}
