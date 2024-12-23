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

/* Library to process 2D wire loops consisting of straight lines and
   circular arcs.  Segments are described by the parameter alpha.
     alpha = dist / (2 * radius)
   where dist is the distance between the endpoints of the segment.
   alpha ranges from -1 to 1, with positive being CCW, negative being
   CW and zero being a straight line.

   Uses Vatti clipping, but only supports a single top and single
   bottom wire.  Scans a vertical line from -x to +x.  Wires convertable
   to this form are refered to as "Vatti singular".
*/

#ifndef LIBCIRCLEWIRE_H
#define LIBCIRCLEWIRE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/*************************************************************************/
/* Data types                                                            */
/*************************************************************************/
struct lcw_seg {
  float pt[2];
  float alpha; /* alpha = dist / (2 * radius) */
  uint32_t flags; /* Reserved */
};

struct lcw_wire;

#define LCW_SLOT_MAIN 0
#define LcW_SLOT_BOT  0
#define LCW_SLOT_TOP  1

struct lcw_list {
  struct lcw_wire *wire;
  struct lcw_list *next;
};

#define LCW_X  0
#define LCW_Y  1
#define LCW_XY 2

struct lcw_properties {
  float area;
  float center_of_mass[2];
  float second_moment[3];
  float third_moment_x[3];
  float min[2];
  float max[2];
};

/*************************************************************************/
/* Functions that work on any wire                                       */
/*************************************************************************/
struct lcw_wire *LCW_New(float seg_combine_tol);
void LCW_Free(struct lcw_wire *wire);
struct lcw_list *LCW_ListNew(void);
void LCW_ListFree(struct lcw_list *list, int free_wires);

size_t LCW_NumSegs(const struct lcw_wire *wire, int slot);
const struct lcw_seg *LCW_Segs(const struct lcw_wire *wire, int slot);

struct lcw_wire *LCW_Copy(const struct lcw_wire *wire);
void LCW_Swap(struct lcw_wire *a, struct lcw_wire *b);
void LCW_Reset(struct lcw_wire *wire, const float *start);
int LCW_To(struct lcw_wire *wire, const float *dst, float alpha);
int LCW_LineTo(struct lcw_wire *wire, const float *dst);
int LCW_CircleRadiusTo(struct lcw_wire *wire, const float *dst, float radius, int large, int cw);
int LCW_CircleCenterTo(struct lcw_wire *wire, const float *dst, const float *center, int cw);
int LCW_Close(struct lcw_wire *wire);

void LCW_Translate(struct lcw_wire *wire, float dx, float dy);
void LCW_Scale(struct lcw_wire *wire, float cx, float cy, float scale);
int LCW_Rotate(struct lcw_wire *wire, float cx, float cy, float ccw_rad);
int LCW_MirrorX(struct lcw_wire *wire);
void LCW_MirrorY(struct lcw_wire *wire);

struct lcw_wire *LCW_FromSvgPath(const char *d_str, int flip_y, float seg_combine_tol);
char *LCW_ToSvgPath(const struct lcw_wire *wire);
void LCW_PathFree(char *path);
int LCW_ToSvg(const char *filename, const struct lcw_wire *wire, const char *unit, int flip_y);
int LCW_ListToSvg(const char *filename, const struct lcw_list *list, const char *unit, int flip_y);

int LCW_IsClosed(const struct lcw_wire *wire);
int LCW_WireLoop(struct lcw_wire *wire); /* Undo VattiClip */
int LCW_IsVatti(const struct lcw_wire *wire);
struct lcw_wire *LCW_ToPolygon(const struct lcw_wire *wire, float tol);
int LCW_Properties(struct lcw_properties *prop, const struct lcw_wire *wire);

/*************************************************************************/
/* Funcitons that require the wire to be Vatti singular                  */
/*************************************************************************/
int LCW_VattiClip(struct lcw_wire *wire);
struct lcw_wire *LCW_ConvexHull(const struct lcw_wire *wire);
struct lcw_list *LCW_CutAtX(const struct lcw_wire *wire, float xval);
struct lcw_list *LCW_ConvexDecomp(const struct lcw_wire *wire, float tol);

/*************************************************************************/
/* Error functions                                                       */
/*************************************************************************/
#define LCW_NO_ERROR               0
#define LCW_OUT_OF_MEMORY         -1
#define LCW_ALPHA_OUT_OF_RANGE    -2
#define LCW_START_AND_END_EQUAL   -3
#define LCW_RADIUS_TO_SMALL       -4
#define LCW_WIRE_NOT_CLOSED       -5
#define LCW_WIRE_NOT_VATTI_SINGLE -6
#define LCW_WIRE_TOO_FEW_SEGMENTS -7
#define LCW_PARSE_ERROR           -8
#define LCW_ELLIPSE_NOT_SUPPORTED -9
#define LCW_COULD_NOT_OPEN_FILE   -10
#define LCW_COULD_NOT_WRITE_FILE  -11
#define LCW_INTERNAL_ERROR        -1000

#ifdef __cplusplus
}
#endif

#endif
