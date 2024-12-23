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
#include "wire.h"

struct str {
  char *str;
  size_t loc;
  size_t alloc;
};

struct str *Str_New(void) {
  struct str *str;

  if ((str = malloc(sizeof(*str))) == NULL)
    goto err;
  memset(str, 0, sizeof(*str));

  str->alloc = 128;
  if ((str->str = malloc(str->alloc)) == NULL)
    goto err2;

  return str;

 err2:
  free(str);
 err:
  return NULL;
}

char *Str_FreeReturnStr(struct str *str) {
  if (str == NULL)
    return NULL;
  
  char *ss = str->str;
  
  free(str);
  return ss;
}

void Str_Free(struct str *str) {
  LCW_PathFree(Str_FreeReturnStr(str));
}

int Str_AddStr(struct str *str, const char *ss) {
  size_t len, target, alloc;
  char *nn;
  
  len = strlen(ss) + 1;
  target = str->loc + len;
  if (target > str->alloc) {
    for (alloc = str->alloc << 1; alloc < target; alloc <<= 1)
      ;
    
    if ((nn = realloc(str->str, alloc)) == NULL)
      return LCW_OUT_OF_MEMORY;
    
    str->str = nn;
    str->alloc = alloc;
  }
  
  memcpy(str->str + str->loc, ss, len);
  str->loc += len - 1;
  return LCW_NO_ERROR;
}
