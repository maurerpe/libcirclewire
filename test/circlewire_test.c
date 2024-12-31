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

#include "libcirclewire.h"

#define POLY_TOL 1e-4

const char first_lines_path[] = "M10,10 H20 V20 H10 Z";

const struct lcw_seg first_lines[] =
  {{{10, 10}, 0, 0},
   {{20, 10}, 0, 0},
   {{20, 20}, 0, 0},
   {{10, 20}, 0, 0},
   {{10, 10}, 0, 0}};

const struct lcw_seg first_lines_vt[] =
  {{{10, 10}, 0, 0},
   {{10, 20}, 0, 0},
   {{20, 20}, 0, 0}};

const struct lcw_seg first_lines_vb[] =
  {{{10, 10}, 0, 0},
   {{20, 10}, 0, 0},
   {{20, 20}, 0, 0}};

const char first_arc_path[] = "M10,10 A100,100,0,0,0,20,10 A5,5,0,0,1,20,20 Z";

const struct lcw_seg first_arc[] =
  {{{10, 10},  0,    0},
   {{20, 10}, -0.05, 0},
   {{20, 20},  1,    0},
   {{10, 10},  0,    0}};

const struct lcw_seg first_arc_vt[] =
  {{{10, 10},  0,                   0},
   {{20, 20},  0,                   0},
   {{25, 15}, -0.7071067811865475,  0}};

const struct lcw_seg first_arc_vb[] =
  {{{10, 10},  0,                   0},
   {{20, 10}, -0.05,                0},
   {{25, 15},  0.7071067811865475,  0}};

const char lines_path[] = "M0,0 V10 H10 V8 V15 V10 H18 L20,8 V0 V-5 H0 Z";

const struct lcw_seg lines[] =
  {{{ 0,  0}, 0, 0},
   {{ 0, 10}, 0, 0},
   {{10, 10}, 0, 0},
   {{10,  8}, 0, 0},
   {{10, 15}, 0, 0},
   {{10, 10}, 0, 0},
   {{18, 10}, 0, 0},
   {{20,  8}, 0, 0},
   {{20, -5}, 0, 0},
   {{0,  -5}, 0, 0},
   {{0,   0}, 0, 0}};

const struct lcw_seg lines_vt[] =
  {{{ 0, -5}, 0, 0},
   {{ 0, 10}, 0, 0},
   {{10, 10}, 0, 0},
   {{10,  8}, 0, 0},
   {{10, 15}, 0, 0},
   {{10, 10}, 0, 0},
   {{18, 10}, 0, 0},
   {{20,  8}, 0, 0}};

const struct lcw_seg lines_vb[] =
  {{{ 0, -5}, 0, 0},
   {{20, -5}, 0, 0},
   {{20,  8}, 0, 0}};

const char arc_path[] = "M5.18748,0 L15,6 A5,5,0,0,1,14,15 A4,4,0,0,0,8,12.9543 A4,4,0,0,0,6.22288,16 L4,17 H-4 A10.9826,10.9826,0,1,1,5.18748,0";

const struct lcw_seg arc[] =
  {{{  5.18748,  0},        0,        0},
   {{ 15,        6},        0,        0},
   {{ 14,       15},        0.905539, 0},
   {{  6.22288, 16},       -0.980145, 0},
   {{  4,       17},        0,        0},
   {{ -4,       17},        0,        0},
   {{-13.6618,   0.795738}, 0.858907, 0},
   {{ 5.18748,   0},        0.858907, 0}};

const struct lcw_seg arc_vt[] =
  {{{-14.9826,  6.01736},  0,        0},
   {{-4,       17},       -0.707104, 0},
   {{ 4,       17},        0,        0},
   {{ 6.22288, 16},        0,        0},
   {{14,       15},        0.980145, 0},
   {{17.3917,  10.2657},  -0.582379, 0}};

const struct lcw_seg arc_vb[] =
  {{{-14.9826,  6.01736},  0,        0},
   {{ 5.18748,  0},        0.958265, 0},
   {{15,        6},        0,        0},
   {{17.3917,  10.2657},   0.489046, 0}};

const float ext_param[] = {40, 0.5, M_PI/6};
const float rev_param[] = {0.8};

#define CHECK_ABS(a, b, t) (!(fabsf((a) - (b)) <= (t)))
#define CHECK(a, b, t, at) (!(fabsf((a) - (b)) <= (t) * fabsf(b)) && CHECK_ABS(a, b, at))

int Check_CutAtX(const char *msg, float xval, float alpha) {
  struct lcw_wire *wire;
  struct lcw_list *list;
  const struct lcw_seg *seg;
  int ret = 0;
  size_t ns;
  float pt[2];
  
  if ((wire = LCW_New(1e-12)) == NULL) {
    fprintf(stderr, "%s: Could not allocate wire\n", msg);
    ret = -1;
    goto err;
  }

  pt[0] = 10;
  pt[1] = 30;
  LCW_Reset(wire, pt);

  pt[0] = -10;
  pt[1] = 35;
  if (LCW_To(wire, pt, alpha) < 0) {
    fprintf(stderr, "%s: Could not create segment 0\n", msg);
    ret = -1;
    goto err2;
  }
  
  pt[0] = -10;
  pt[1] = -35;
  if (LCW_To(wire, pt, 0) < 0) {
    fprintf(stderr, "%s: Could not create segment 1\n", msg);
    ret = -1;
    goto err2;
  }
  
  pt[0] = 10;
  pt[1] = -30;
  if (LCW_To(wire, pt, alpha) < 0) {
    fprintf(stderr, "%s: Could not create segment 2\n", msg);
    ret = -1;
    goto err2;
  }
  
  if (LCW_Close(wire) < 0) {
    fprintf(stderr, "%s: Could not close wire\n", msg);
    ret = -1;
    goto err2;
  }
  
  if ((list = LCW_CutAtX(wire, xval)) == NULL) {
    fprintf(stderr, "%s: Could not cut wire at %g\n", msg, xval);
    ret = -1;
    goto err2;
  }
  
  if (list->next == NULL) {
    fprintf(stderr, "%s: Only one wire after cut\n", msg);
    ret = -1;
    goto err3;
  }
  
  seg = LCW_Segs(list->wire, 0);
  ns  = LCW_NumSegs(list->wire, 0);
  if (CHECK_ABS(seg[ns-1].pt[0], xval, 1e-5)) {
    fprintf(stderr, "%s: Incorrect point 1: expected %g, found %g\n", msg, xval, seg[ns-1].pt[0]);
    ret = -1;
  }
  
  seg = LCW_Segs(list->wire, 1);
  ns  = LCW_NumSegs(list->wire, 1);
  if (CHECK_ABS(seg[ns-1].pt[0], xval, 1e-5)) {
    fprintf(stderr, "%s: Incorrect point 2: expected %g, found %g\n", msg, xval, seg[ns-2].pt[0]);
    ret = -1;
  }
  
 err3:
  LCW_ListFree(list, 1);
 err2:
  LCW_Free(wire);
 err:
  return ret;
}

int Verify_Segs(const char *msg, const struct lcw_wire *wire, const struct lcw_seg **seg, size_t *num_seg, int check_path) {
  size_t num_wire, idx;
  const struct lcw_seg *ws;
  struct lcw_wire *wp;
  char *path, buf[1024];
  int ret = 0, slot;

  for (slot = 0; slot < 2; slot++) {
    num_wire = LCW_NumSegs(wire, slot);
    if (num_wire != num_seg[slot]) {
      fprintf(stderr, "%s: Slot %d: Wrong number of segments, expected %zu, found %zu\n",
	      msg, slot, num_seg[slot], num_wire);
      if (num_seg[slot] > num_wire)
	num_seg[slot] = num_wire;
      ret = -1;
    }

    ws = LCW_Segs(wire, slot);
    for (idx = 0; idx < num_seg[slot]; idx++) {
      if (CHECK(ws[idx].pt[0], seg[slot][idx].pt[0], 1e-5, 1e-5)) {
	fprintf(stderr, "%s: Slot %d: Segment %zu X mismatch, expected %g, found %g\n",
		msg, slot, idx, seg[slot][idx].pt[0], ws[idx].pt[0]);
	ret = -1;
      }
      if (CHECK(ws[idx].pt[1], seg[slot][idx].pt[1], 1e-5, 1e-5)) {
	fprintf(stderr, "%s: Slot %d: Segment %zu Y mismatch, expected %g, found %g\n",
		msg, slot, idx, seg[slot][idx].pt[1], ws[idx].pt[1]);
	ret = -1;
      }
      if (CHECK(ws[idx].alpha, seg[slot][idx].alpha, 1e-5, 1e-5)) {
	fprintf(stderr, "%s: Slot %d: Segment %zu alpha mismatch, expected %g, found %g\n",
		msg, slot, idx, seg[slot][idx].alpha, ws[idx].alpha);
	ret = -1;
      }
      if (ws[idx].flags != seg[slot][idx].flags) {
	fprintf(stderr, "%s: Slot %d: Segment %zu flags mismatch, expected 0x%lx, found 0x%lx\n",
		msg, slot, idx, (long) seg[slot][idx].flags, (long) ws[idx].flags);
	ret = -1;
      }
    }
  }
  
  if (ret < 0 || !check_path)
    return ret;
  
  if ((path = LCW_ToSvgPath(wire)) == NULL) {
    fprintf(stderr, "%s: Could not create path from wire\n", msg);
    return -1;
  }
  
  wp = LCW_FromSvgPath(path, 0, 1e-6);
  LCW_PathFree(path);
  if (wp == NULL) {
    fprintf(stderr, "%s: Could not convert back to wire from path\n", msg);
    return -1;
  }

  snprintf(buf, sizeof(buf), "%s: To path and back", msg);
  buf[sizeof(buf) - 1] = '\0';
  ret = Verify_Segs(buf, wp, seg, num_seg, 0);
  LCW_Free(wp);
  return ret;
}

#define CHECK_PROP(pp, at)						\
  do {									\
    if (CHECK(prop->pp, exp->pp, tol, (at))) {				\
      fprintf(stderr, "%s: " #pp " mismatch, expected %g, found %g\n",	\
	      msg, exp->pp, prop->pp);					\
      ret = -1;								\
    }									\
  } while(0)

int Check_Prop(const char *msg, const struct lcw_properties *prop, const struct lcw_properties *exp, float tol) {
  int ret = 0;

  CHECK_PROP(area, tol);
  CHECK_PROP(center_of_mass[0], tol);
  CHECK_PROP(center_of_mass[1], tol);
  CHECK_PROP(second_moment[0], tol);
  CHECK_PROP(second_moment[1], tol);
  CHECK_PROP(second_moment[2], tol);
  CHECK_PROP(third_moment_x[0], tol);
  CHECK_PROP(third_moment_x[1], tol);
  CHECK_PROP(third_moment_x[2], tol);
  CHECK_PROP(min[0], tol);
  CHECK_PROP(min[1], tol);
  CHECK_PROP(max[0], tol);
  CHECK_PROP(max[1], tol);
  
  return ret;
}

int Check_SolidProp(const char *msg, const struct lp_mass_properties *prop, const struct lp_mass_properties *exp, float tol) {
  int ret = 0;
  float tr = 1e-3 * (exp->inertia_tensor[0] + exp->inertia_tensor[4] + exp->inertia_tensor[8]);

  CHECK_PROP(volume, tol);
  CHECK_PROP(center_of_mass[0], tol);
  CHECK_PROP(center_of_mass[1], tol);
  CHECK_PROP(center_of_mass[2], tol);
  CHECK_PROP(inertia_tensor[0], tol * tr);
  CHECK_PROP(inertia_tensor[1], tol * tr);
  CHECK_PROP(inertia_tensor[2], tol * tr);
  CHECK_PROP(inertia_tensor[3], tol * tr);
  CHECK_PROP(inertia_tensor[4], tol * tr);
  CHECK_PROP(inertia_tensor[5], tol * tr);
  CHECK_PROP(inertia_tensor[6], tol * tr);
  CHECK_PROP(inertia_tensor[7], tol * tr);
  CHECK_PROP(inertia_tensor[8], tol * tr);
  
  return ret;
}

struct lcw_wire *Test_Path(const char *msg, const char *path, const struct lcw_seg **seg, size_t *num_seg, int check_path) {
  struct lcw_wire *wire;

  if ((wire = LCW_FromSvgPath(path, 0, 1e-4)) == NULL) {
    fprintf(stderr, "%s: Could not create wire from path\n", msg);
    return NULL;
  }
  
  if (Verify_Segs(msg, wire, seg, num_seg, check_path) < 0) {
    LCW_Free(wire);
    return NULL;
  }
  
  return wire;
}

int Full_Test(const char *msg, const char *path, const struct lcw_seg *seg, size_t num_seg, int check_path, int expect_closed, const struct lcw_seg *seg_vt, size_t num_seg_vt, const struct lcw_seg *seg_vb, size_t num_seg_vb, const char *svg_filename) {
  struct lcw_wire *wire = NULL, *poly = NULL, *vatti = NULL, *vpoly = NULL;
  struct lcw_properties prop[4];
  struct lcw_list list[4], *cut = NULL;
  const struct lcw_seg *seg_arr[] = {seg, NULL}, *seg_v_arr[] = {seg_vb, seg_vt};
  size_t num_seg_arr[] = {num_seg, 0}, num_seg_v_arr[] = {num_seg_vb, num_seg_vt};
  float mid;
  char buf[1024];
  int err = 0, ret = 0, idx = 0;

  memset(prop, 0, sizeof(prop));
  
  if ((wire = Test_Path(msg, path, seg_arr, num_seg_arr, check_path)) == NULL)
    return -1;
  if (LCW_IsClosed(wire) != expect_closed) {
    fprintf(stderr, "%s: Wire not %s\n", msg, expect_closed ? "closed" : "open");
    err = -1;
  }
  if ((poly = LCW_ToPolygon(wire, POLY_TOL)) == NULL) {
    fprintf(stderr, "%s: Could not convert wire to polygon\n", msg);
    err = -1;
  }
  if (poly && LCW_IsClosed(poly) != expect_closed) {
    fprintf(stderr, "%s: Polygon not %s\n", msg, expect_closed ? "closed" : "open");
    err = -1;
  }
  if ((ret = LCW_Properties(&prop[0], wire)) != LCW_NO_ERROR) {
    fprintf(stderr, "%s: Could not get properties: %d\n", msg, ret);
    err = -1;
  }
  if (poly && (ret = LCW_Properties(&prop[1], poly)) != LCW_NO_ERROR) {
    fprintf(stderr, "%s: Could not get polygon properties: %d\n", msg, ret);
    err = -1;
  }
  if (poly && (ret = Check_Prop(msg, &prop[0], &prop[1], POLY_TOL)) < 0)
    err = -1;
  if (seg_vt) {
    if ((vatti = LCW_Copy(wire)) == NULL) {
      fprintf(stderr, "%s: Could copy wire for Vatti clip\n", msg);
      err = -1;
    }
    if (vatti && (ret = LCW_VattiClip(vatti)) != LCW_NO_ERROR) {
      fprintf(stderr, "%s: Could not Vatti clip: %d\n", msg, ret);
      err = -1;
    }
    snprintf(buf, sizeof(buf), "%s: Vatti clip", msg);
    buf[sizeof(buf) - 1] = '\0';
    if (vatti && Verify_Segs(buf, vatti, seg_v_arr, num_seg_v_arr, 0) < 0)
      err = -1;
    if (vatti && (vpoly = LCW_ToPolygon(vatti, POLY_TOL)) == NULL) {
      fprintf(stderr, "%s: Could not convert vatti to polygon\n", buf);
      err = -1;
    }
    if (vatti && (ret = LCW_Properties(&prop[2], vatti)) != LCW_NO_ERROR) {
      fprintf(stderr, "%s: Could not get Vatti properties: %d\n", buf, ret);
      err = -1;
    }
    if (prop[0].area < 0) {
      prop[0].area *= -1;
      prop[0].second_moment[0] *= -1;
      prop[0].second_moment[1] *= -1;
      prop[0].second_moment[2] *= -1;
      prop[0].third_moment_x[0] *= -1;
      prop[0].third_moment_x[1] *= -1;
      prop[0].third_moment_x[2] *= -1;
    }
    if (vpoly && (ret = Check_Prop(buf, &prop[2], &prop[0], POLY_TOL)) < 0)
      err = -1;
    if (vpoly && (ret = LCW_Properties(&prop[3], vpoly)) != LCW_NO_ERROR) {
      fprintf(stderr, "%s: Could not get Vatti polygon properties: %d\n", buf, ret);
      err = -1;
    }
    snprintf(buf, sizeof(buf), "%s: Vatti clip polygon", msg);
    buf[sizeof(buf) - 1] = '\0';
    if (vpoly && (ret = Check_Prop(buf, &prop[3], &prop[2], POLY_TOL)) < 0)
      err = -1;
    if (num_seg_vt > 0) {
      mid = (seg_vt[0].pt[0] + seg_vt[num_seg_vt - 1].pt[0]) / 2;
      if ((cut = LCW_CutAtX(vatti, mid)) == NULL) {
	fprintf(stderr, "%s: Vatti Clip: Could not cut at %g\n", msg, mid);
	err = -1;
      }
    }
  }
  memset(list, 0, sizeof(list));
  list[0].wire = wire;
  if (poly) {
    list[idx].next = &list[idx+1];
    list[idx+1].wire = poly;
    idx++;
  }
  if (vatti) {
    list[idx].next = &list[idx+1];
    list[idx+1].wire = vatti;
    idx++;
  }
  if (vpoly) {
    list[idx].next = &list[idx+1];
    list[idx+1].wire = vpoly;
    idx++;
  }
  list[idx].next = cut;
  if (svg_filename && (ret = LCW_ListToSvg(svg_filename, list, "mm", 1)) < 0) {
    fprintf(stderr, "%s: Could not write svg to file: %d\n", msg, ret);
    err = -1;
  }

  LCW_ListFree(cut, 1);
  LCW_Free(vpoly);
  LCW_Free(vatti);
  LCW_Free(poly);
  LCW_Free(wire);
  return err;
}

static int DecompTest(const char *msg, const char *path, const char *svg_filename) {
  struct lcw_wire *wire, *hull;
  struct lcw_list list[2], *decomp;
  int ret_val = -1, ret;
  
  if ((wire = LCW_FromSvgPath(path, 0, 1e-4)) == NULL) {
    fprintf(stderr, "%s: Could not create wire from path\n", msg);
    goto err;
  }
  
  if ((ret = LCW_VattiClip(wire)) != LCW_NO_ERROR) {
    fprintf(stderr, "%s: Could vatti clip path: %d\n", msg, ret);
    goto err2;
  }
  
  if ((hull = LCW_ConvexHull(wire)) == NULL) {
    fprintf(stderr, "%s: Could not find convex hull\n", msg);
    goto err2;
  }

  if ((decomp = LCW_ConvexDecomp(wire, 1e-2)) == NULL) {
    fprintf(stderr, "%s: Could not find convex decomposition\n", msg);
    goto err3;
  }
  
  memset(list, 0, sizeof(list));
  list[0].wire = wire;
  list[0].next = &list[1];
  list[1].wire = hull;
  list[1].next = decomp;
  if (svg_filename && (ret = LCW_ListToSvg(svg_filename, list, "mm", 1)) < 0) {
    fprintf(stderr, "%s: Could not write svg to file: %d\n", msg, ret);
    goto err4;
  }
  
  ret_val = 0;
  /* Fall through */

 err4:
  LCW_ListFree(decomp, 1);
 err3:
  LCW_Free(hull);
 err2:
  LCW_Free(wire);
 err:
  return ret_val;
}

static int SolidTest(const char *msg, const char *path, enum lcw_solid_type type, const float *param, size_t num_param, const char *obj_path) {
  struct lcw_wire *wire;
  struct lcw_solid *solid;
  struct lp_mass_properties sprop, mprop;
  struct lp_vertex_list *mesh;
  struct lp_vl_list *decomp, *out, *mesh_list;
  int ret_val = 0;
  
  if ((wire = LCW_FromSvgPath(path, 0, 1e-5)) == NULL) {
    fprintf(stderr, "%s: Could not create wire from path\n", msg);
    ret_val = -1;
    goto err;
  }
  
  if ((solid = LCW_Solid(wire, type, param, num_param)) == NULL) {
    fprintf(stderr, "%s: Could not create solid\n", msg);
    ret_val = -1;
    goto err2;
  }
  
  if (LCW_SolidProperties(&sprop, solid) != LCW_NO_ERROR) {
    fprintf(stderr, "%s: Could not get solid properties\n", msg);
    ret_val = -1;
    goto err3;
  }
  
  if ((mesh = LCW_SolidMesh(solid, 1e-4)) == NULL) {
    fprintf(stderr, "%s: Could not mesh solid\n", msg);
    ret_val = -1;
  }
  
  LP_MassProperties(mesh, &mprop);
  if (Check_SolidProp(msg, &sprop, &mprop, POLY_TOL) < 0)
    ret_val = -1;
  
  if ((decomp = LCW_SolidConvexDecomp(solid, 1e-2, 1e-4)) == NULL) {
    fprintf(stderr, "%s: Could not perform solid convex decomposition\n", msg);
    ret_val = -1;
  }
  
  if (mesh && (mesh_list = LP_VertexList_ListAppend(NULL, mesh)) == NULL) {
    fprintf(stderr, "%s: Could not create list from mesh\n", msg);
    ret_val = -1;
  }
  if (mesh_list)
    mesh = NULL;

  if ((out = LP_VertexList_ListJoin(mesh_list, decomp)) == NULL) {
    fprintf(stderr, "%s: Could not join lists\n", msg);
    ret_val = -1;
    goto err4;
  }
  if (out) {
    mesh_list = NULL;
    decomp = NULL;
  }
  
  if (obj_path && LP_VertexList_Write(obj_path, out, 1) < 0) {
    fprintf(stderr, "%s: Could not write out obj to '%s'\n", msg, obj_path);
    ret_val = -1;
  }
  
  /* Fall through */

  LP_VertexList_ListFree(out);
 err4:
  LP_VertexList_ListFree(mesh_list);
  LP_VertexList_ListFree(decomp);
  LP_VertexList_Free(mesh);
 err3:
  LCW_SolidFree(solid);
 err2:
  LCW_Free(wire);
 err:
  return ret_val;
}

int main(void) {
  if (Check_CutAtX("Cut 1", 5, 0.05) < 0)
    exit(1);
  
  if (Check_CutAtX("Cut 2", 5, 0.2) < 0)
    exit(1);
  
  if (Check_CutAtX("Cut 3", -5, 0.05) < 0)
    exit(1);
  
  if (Check_CutAtX("Cut 4", -5, 0.2) < 0)
    exit(1);
  
  if (Full_Test("First lines",
		first_lines_path,
		first_lines, sizeof(first_lines) / sizeof(*first_lines),
		1,
		1,
		first_lines_vt, sizeof(first_lines_vt) / sizeof(*first_lines_vt),
		first_lines_vb, sizeof(first_lines_vb) / sizeof(*first_lines_vb),
		NULL) < 0)
    exit(1);
  
  if (Full_Test("First arc",
		first_arc_path,
		first_arc, sizeof(first_arc) / sizeof(*first_arc),
		1,
		1,
		first_arc_vt, sizeof(first_arc_vt) / sizeof(*first_arc_vt),
		first_arc_vb, sizeof(first_arc_vb) / sizeof(*first_arc_vb),
		"first_arc.svg") < 0)
    exit(1);
  
  if (Full_Test("Lines",
		lines_path,
		lines, sizeof(lines) / sizeof(*lines),
		1,
		1,
		lines_vt, sizeof(lines_vt) / sizeof(*lines_vt),
		lines_vb, sizeof(lines_vb) / sizeof(*lines_vb),
		"lines.svg") < 0)
    exit(1);
  
  if (Full_Test("Arc",
		arc_path,
		arc, sizeof(arc) / sizeof(*arc),
		1,
		1,
		arc_vt, sizeof(arc_vt) / sizeof(*arc_vt),
		arc_vb, sizeof(arc_vb) / sizeof(*arc_vb),
		"arc.svg") < 0)
    exit(1);

  if (DecompTest("First arc decomp", first_arc_path, "first_arc_decomp.svg") < 0)
    exit(1);
  
  if (DecompTest("Lines decomp", lines_path, "lines_decomp.svg") < 0)
    exit(1);
  
  if (DecompTest("Arc decomp", arc_path, "arc_decomp.svg") < 0)
    exit(1);
  
  if (SolidTest("First arc extrude",
		first_arc_path,
		lcw_extrude,
		ext_param,
		sizeof(ext_param) / sizeof(*ext_param),
		"first_arc_extrude.obj") < 0)
    exit(1);
  
  if (SolidTest("First arc revolve",
		first_arc_path,
		lcw_revolve,
		rev_param,
		sizeof(rev_param) / sizeof(*rev_param),
		"first_arc_revolve.obj") < 0)
    exit(1);
  
  if (SolidTest("Arc extrude",
		arc_path,
		lcw_extrude,
		ext_param,
		sizeof(ext_param) / sizeof(*ext_param),
		"arc_extrude.obj") < 0)
    exit(1);
  
  fprintf(stderr, "All tests passed\n");
  exit(0);
}
