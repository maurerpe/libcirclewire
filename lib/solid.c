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

struct lcw_solid {
  struct lcw_wire *wire;
  enum lcw_solid_type type;
  float param[3];
};

struct lcw_solid *LCW_Solid(const struct lcw_wire *wire, enum lcw_solid_type type, const float *param, size_t num_params) {
  struct lcw_solid *solid;
  struct lcw_properties prop;
  float dx = 0, dy = 0;
  
  if ((solid = malloc(sizeof(*solid))) == NULL)
    goto err;
  memset(solid, 0, sizeof(*solid));

  solid->type = type;
  switch (type) {
  case lcw_extrude:
    LCW_Properties(&prop, wire);
    dx = -prop.center_of_mass[0];
    dy = -prop.center_of_mass[1];
    
    if (num_params < 1 || num_params > 3)
      goto err2;
    solid->param[0] = param[0];

    if (num_params < 2)
      solid->param[1] = 1;
    else if (param[1] > 0 && param[1] < 2)
      solid->param[1] = param[1];
    else
      goto err2;

    if (num_params < 3)
      solid->param[2] = 0;
    else if (param[2] > -M_PI && param[2] < M_PI)
      solid->param[2] = param[2];
    else
      goto err2;
    break;
    
  case lcw_revolve:
    LCW_Properties(&prop, wire);
    if (prop.min[1] < 0)
      goto err2;
    
    if (num_params > 1)
      goto err2;

    if (num_params < 1)
      solid->param[0] = M_PI;
    else if (param[0] > 0 && param[0] <= M_PI)
      solid->param[0] = param[0];
    else
      goto err2;
    break;

  default:
    goto err2;
  }
  
  if ((solid->wire = LCW_Copy(wire)) == NULL)
    goto err2;
  if (dx != 0 || dy != 0)
    LCW_Translate(solid->wire, dx, dy);
  if (LCW_VattiClip(solid->wire) != LCW_NO_ERROR)
    goto err3;
  
  return solid;

 err3:
  LCW_Free(solid->wire);
 err2:
  free(solid);
 err:
  return NULL;
}

void LCW_SolidFree(struct lcw_solid *solid) {
  if (solid == NULL)
    return;
  
  LCW_Free(solid->wire);
  free(solid);
}

enum lcw_solid_type LCW_SolidType(const struct lcw_solid *solid) {
  return solid->type;
}

const struct lcw_wire *LCW_SolidWire(const struct lcw_solid *solid) {
  return solid->wire;
}

const float *LCW_SolidParams(const struct lcw_solid *solid) {
  return solid->param;
}

size_t LCW_SolidNumParams(const struct lcw_solid *solid) {
  if (solid->type == lcw_extrude)
    return 3;
  return 1;
}

int LCW_SolidProperties(struct lp_mass_properties *prop, const struct lcw_solid *solid) {
  struct lcw_properties prop2d;
  float h, scale, tan_sw, h2, h3, s2, s4, ss4;
  float ang, cx, cy, dx, dy, rcyl, Ixxx, Ixyy, Ixxy, q, j, k, sin_a, sinc_a;
  int ret;

  memset(prop, 0, sizeof(*prop));
  
  if ((ret = LCW_Properties(&prop2d, solid->wire)) != LCW_NO_ERROR)
    return ret;

  if (solid->type == lcw_extrude) {
    h = 2 * solid->param[0];
    scale = solid->param[1];
    tan_sw = tan(solid->param[2]);

    h2 = h * h;
    h3 = h2 * h;
    s2 = (scale - 2) * scale + 4;
    s4 = ((scale - 4) * scale + 4) * scale * scale + 4;
    ss4 = (((scale - 4) * scale + 16) * scale - 24) * scale + 16;
    prop->volume = prop2d.area * h * s2 / 3;
    prop->center_of_mass[0] = prop2d.area * h2 * (scale - 1) * tan_sw / (3 * prop->volume);
    prop->center_of_mass[1] = 0;
    prop->center_of_mass[2] = prop2d.area * h2 * (scale - 1) / (3 * prop->volume);
    
    prop->inertia_tensor[0] = prop2d.second_moment[LCW_X] * h * ss4 / 5 + prop2d.area * h3 * s4 / (20 * s2);
    prop->inertia_tensor[4] = prop2d.second_moment[LCW_Y] * h * ss4 / 5 + prop2d.area * h3 * s4 * (tan_sw * tan_sw + 1) / (20 * s2);
    prop->inertia_tensor[8] = (prop2d.second_moment[LCW_X] + prop2d.second_moment[LCW_Y]) * h * ss4 / 5 + prop2d.area * h3 * s4 * tan_sw * tan_sw / (20 * s2);

    prop->inertia_tensor[1] = prop->inertia_tensor[3] = -prop2d.second_moment[LCW_XY] * h * ss4 / 5;
    prop->inertia_tensor[2] = prop->inertia_tensor[6] = -prop2d.area * h3 * s4 * tan_sw / (20 * s2);
    prop->inertia_tensor[5] = prop->inertia_tensor[7] = 0;
  } else {
    ang = solid->param[0];
    cx = prop2d.center_of_mass[LCW_X];
    cy = prop2d.center_of_mass[LCW_Y];
    dx = prop2d.second_moment[LCW_XY] / (prop2d.area * cy);
    dy = prop2d.second_moment[LCW_X] / (prop2d.area * cy);
    rcyl = cy + dy;
    Ixxx = prop2d.third_moment_x[LCW_X]  - prop2d.area * cy * dy * dy;
    Ixyy = prop2d.third_moment_x[LCW_Y]  - prop2d.area * cy * dx * dx;
    Ixxy = prop2d.third_moment_x[LCW_XY] - prop2d.area * cy * dx * dy;
    q = prop2d.area * cy * rcyl * rcyl;
    j = (q + Ixxx) * sin(2 * ang) / 2;
    k = (q + Ixxx + 2 * Ixyy) * ang;
    sin_a = sin(ang);
    if (ang == 0)
      sinc_a = 0;
    else
      sinc_a = sin_a / ang;
    prop->volume = prop2d.area * cy * 2 * ang;
    prop->center_of_mass[0] = cx + dx;
    prop->center_of_mass[1] = sinc_a * rcyl;
    prop->center_of_mass[2] = 0;

    prop->inertia_tensor[0] = 2 * (q + Ixxx) * ang - 2 * q * sin_a * sinc_a;
    prop->inertia_tensor[4] = k - j;
    prop->inertia_tensor[8] = k + j - 2 * q * sin_a * sinc_a;
    prop->inertia_tensor[1] = prop->inertia_tensor[3] = -2 * Ixxy * sin_a;
    prop->inertia_tensor[2] = prop->inertia_tensor[6] = 0;
    prop->inertia_tensor[5] = prop->inertia_tensor[7] = 0;
  }
  
  return LCW_NO_ERROR;
}

static int SamePt(float *pt1, float *pt2, float comb_tol) {
  return
    fabsf(pt1[0] - pt2[0]) < comb_tol &&
    fabsf(pt1[1] - pt2[1]) < comb_tol &&
    fabsf(pt1[2] - pt2[2]) < comb_tol;
}

static int AddTri(struct lp_vertex_list *vl, float *pt1, float *pt2, float *pt3, float comb_tol) {
  if (SamePt(pt1, pt2, comb_tol) ||
      SamePt(pt2, pt3, comb_tol) ||
      SamePt(pt1, pt3, comb_tol))
    return 0;
  
  if (LP_VertexList_Add(vl, pt1) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(vl, pt2) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(vl, pt3) == UINT_MAX)
    return -1;

  return 0;
}

static int AddQuad(struct lp_vertex_list *vl, float *vert, float comb_tol) {
  if (AddTri(vl, &vert[1 * 3], &vert[0 * 3], &vert[2 * 3], comb_tol) < 0)
    return -1;
  if (AddTri(vl, &vert[2 * 3], &vert[3 * 3], &vert[1 * 3], comb_tol) < 0)
    return -1;
  
  return 0;
}

struct lp_vertex_list *LCW_SolidMesh(const struct lcw_solid *solid, float tol) {
  struct lcw_wire *poly;
  struct lp_vertex_list *vl, *end;
  struct lcw_properties prop;
  float vert[4 * 3], *vbase, *vv, h, scale, tan_sw;
  float ang, csa, sna, ang1, ang2, cos_a1, sin_a1, cos_a2, sin_a2;
  unsigned int *ind;
  int vidx, off;
  size_t fpv, idx, cnt, num, inum;
  
  if ((poly = LCW_ToPolygon(solid->wire, tol)) == NULL)
    goto err;
  if (LCW_WireLoop(poly) != LCW_NO_ERROR)
    goto err2;
  
  if ((end = LCW_Mesh(poly, tol)) == NULL)
    goto err2;

  if ((vl = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
    goto err3;

  if (solid->type == lcw_extrude) {
    h = solid->param[0];
    scale = solid->param[1];
    tan_sw = tan(solid->param[2]);

    fpv   = LP_VertexList_FloatsPerVert(end);
    inum  = LP_VertexList_NumInd(end);
    vbase = LP_VertexList_GetVert(end);
    ind   = LP_VertexList_GetInd(end);
    vidx = 0;
    while (inum-- > 0) {
      off = 2 * (1 - vidx); // Switch ind 0 and 2
      vv = vbase + ind[off] * fpv;
      vert[0] = vv[0] * scale + tan_sw * h;
      vert[1] = vv[1] * scale;
      vert[2] = h;
      if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	goto err4;
      ind++;
      vidx = (vidx + 1) % 3;
    }
    inum = LP_VertexList_NumInd(end);
    ind  = LP_VertexList_GetInd(end);
    while (inum-- > 0) {
      vv = vbase + *ind * fpv;
      vert[0] = vv[0] * (2 - scale) - tan_sw * h;
      vert[1] = vv[1] * (2 - scale);;
      vert[2] = -h;
      if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	goto err4;
      ind++;
    }
    
    for (idx = 1; idx < poly->num_seg[0]; idx++) {
      vert[ 0] = poly->seg[0][idx-1].pt[0] * scale + tan_sw * h;
      vert[ 1] = poly->seg[0][idx-1].pt[1] * scale;
      vert[ 2] = h;
      vert[ 3] = poly->seg[0][idx  ].pt[0] * scale + tan_sw * h;
      vert[ 4] = poly->seg[0][idx  ].pt[1] * scale;
      vert[ 5] = h;
      vert[ 6] = poly->seg[0][idx-1].pt[0] * (2 - scale) - tan_sw * h;
      vert[ 7] = poly->seg[0][idx-1].pt[1] * (2 - scale);
      vert[ 8] = -h;
      vert[ 9] = poly->seg[0][idx  ].pt[0] * (2 - scale) - tan_sw * h;
      vert[10] = poly->seg[0][idx  ].pt[1] * (2 - scale);
      vert[11] = -h;
      if (AddQuad(vl, vert, poly->comb_tol) < 0)
	goto err4;
    }
  } else {
    ang = solid->param[0];

    LCW_Properties(&prop, poly);
    if (ang < M_PI - poly->comb_tol / prop.max[1]) {
      fpv   = LP_VertexList_FloatsPerVert(end);
      inum  = LP_VertexList_NumInd(end);
      ind   = LP_VertexList_GetInd(end);
      vbase = LP_VertexList_GetVert(end);
      csa   = cos(ang);
      sna   = sin(ang);
      vidx = 0;
      while (inum-- > 0) {
	off = 2 * (1 - vidx); // Switch ind 0 and 2
	vv = vbase + ind[off] * fpv;
	vert[0] = vv[0];
	vert[1] = vv[1] * csa;
	vert[2] = vv[1] * sna;
	if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	  goto err4;
	ind++;
	vidx = (vidx + 1) % 3;
      }
      inum = LP_VertexList_NumInd(end);
      ind  = LP_VertexList_GetInd(end);
      while (inum-- > 0) {
	vv = vbase + *ind * fpv;
	vert[0] = vv[0];
	vert[1] = vv[1] * csa;
	vert[2] = vv[1] * -sna;
	if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	  goto err4;
	ind++;
      }
    }

    num = RevRequiredSegs(2 * ang, prop.max[1], tol);
    for (cnt = 0; cnt < num; cnt++) {
      ang1 = ((float) cnt) / num * 2 * ang - ang;
      ang2 = ((float) (cnt + 1)) / num * 2 * ang - ang;
      cos_a1 = cos(ang1);
      sin_a1 = sin(ang1);
      cos_a2 = cos(ang2);
      sin_a2 = sin(ang2);
      for (idx = 1; idx < poly->num_seg[0]; idx++) {
	vert[ 0] = poly->seg[0][idx-1].pt[0];
	vert[ 1] = poly->seg[0][idx-1].pt[1] * cos_a2;
	vert[ 2] = poly->seg[0][idx-1].pt[1] * sin_a2;
	vert[ 3] = poly->seg[0][idx  ].pt[0];
	vert[ 4] = poly->seg[0][idx  ].pt[1] * cos_a2;
	vert[ 5] = poly->seg[0][idx  ].pt[1] * sin_a2;
	vert[ 6] = poly->seg[0][idx-1].pt[0];
	vert[ 7] = poly->seg[0][idx-1].pt[1] * cos_a1;
	vert[ 8] = poly->seg[0][idx-1].pt[1] * sin_a1;
	vert[ 9] = poly->seg[0][idx  ].pt[0];
	vert[10] = poly->seg[0][idx  ].pt[1] * cos_a1;
	vert[11] = poly->seg[0][idx  ].pt[1] * sin_a1;
	if (AddQuad(vl, vert, poly->comb_tol) < 0)
	  goto err4;
      }
    }
  }
  
  LP_VertexList_Free(end);
  LCW_Free(poly);
  return vl;
  
 err4:
  LP_VertexList_Free(vl);
 err3:
  LP_VertexList_Free(end);
 err2:
  LCW_Free(poly);
 err:
  return NULL;
}

static float MaxScale(const struct lcw_solid *solid) {
  float scale;
  
  if (solid->type != lcw_extrude)
    return 1;
  
  scale = solid->param[1];
  if (scale >= 1)
    return scale;
  return 2 - scale;
}

static int AddHull(struct lp_vl_list ***tail, const struct lp_vertex_list *vl) {
  struct lp_vertex_list *hull;
  struct lp_vl_list *nn;
  
  if ((hull = LP_ConvexHull(vl)) == NULL)
    goto err;
  
  if ((nn = LP_VertexList_ListAppend(NULL, hull)) == NULL)
    goto err2;

  **tail = nn;
  *tail = &nn->next;
  return 0;

 err2:
  LP_VertexList_Free(hull);
 err:
  return -1;
}

static int AddExtrude(struct lcw_wire *poly, const float *param, struct lp_vl_list ***tail) {
  struct lp_vertex_list *vl;
  const struct lcw_seg *seg;
  float h, scale, tan_sw, vert[3];
  size_t cnt, ns;

  if (LCW_WireLoop(poly) != LCW_NO_ERROR)
    goto err;
  
  if ((vl = LP_VertexList_New(3, lp_pt_point)) == NULL)
    goto err;
  h      = param[0];
  scale  = param[1];
  tan_sw = tan(param[2]);

  seg = poly->seg[0];
  ns = poly->num_seg[0];
  for (cnt = 1; cnt < ns; cnt++, seg++) {
    vert[0] = seg->pt[0] * scale + tan_sw * h;
    vert[1] = seg->pt[1] * scale;
    vert[2] = h;
    if (LP_VertexList_Add(vl, vert) == UINT_MAX)
      goto err2;
    vert[0] = seg->pt[0] * (2 - scale) - tan_sw * h;
    vert[1] = seg->pt[1] * (2 - scale);
    vert[2] = -h;
    if (LP_VertexList_Add(vl, vert) == UINT_MAX)
      goto err2;
  }
  
  if (AddHull(tail, vl) < 0)
    goto err2;

  LP_VertexList_Free(vl);
  return 0;
  
 err2:
  LP_VertexList_Free(vl);
 err:
  return -1;
}

static int AddRevolve(struct lcw_wire *poly, const float *param, struct lp_vl_list ***tail, float dtol, float ptol) {
  struct lp_vertex_list *vl;
  struct lcw_seg *seg;
  struct lcw_properties prop;
  float ang, max_inner, da, start, aa, ca, sa, vert[3];
  size_t cnt, ns, dnum, pnum, dcnt, pcnt;
  
  if (LCW_VattiClip(poly) != LCW_NO_ERROR)
    goto err;
  
  ang = param[0];
  
  max_inner = 0;
  seg = poly->seg[0];
  ns = poly->num_seg[0];
  for (cnt = 0; cnt < ns; cnt++, seg++) {
    if (cnt > 0 && cnt == ns - 1 && fabsf(seg[-1].pt[0] - seg[0].pt[0]) < poly->comb_tol)
      continue;
    if (seg->pt[1] > max_inner)
      max_inner = seg->pt[1];
  }

  if (LCW_WireLoop(poly) != LCW_NO_ERROR)
    goto err;

  if (LCW_Properties(&prop, poly) != LCW_NO_ERROR)
    goto err;
  dnum = RevRequiredSegs(2 * ang, max_inner, dtol);
  pnum = RevRequiredSegs(2 * ang / dnum, prop.max[1], ptol);

  da = 2 * ang / dnum;
  ns = poly->num_seg[0];
  for (dcnt = 0; dcnt < dnum; dcnt++) {
    if ((vl = LP_VertexList_New(3, lp_pt_point)) == NULL)
      goto err;
    start = 2 * ang * dcnt / dnum - ang;
    
    for (pcnt = 0; pcnt <= pnum; pcnt++) {
      aa = da * pcnt / pnum + start;
      ca = cos(aa);
      sa = sin(aa);
      seg = poly->seg[0];
      for (cnt = 0; cnt < ns; cnt++, seg++) {
	vert[0] = seg->pt[0];
	vert[1] = seg->pt[1] * ca;
	vert[2] = seg->pt[1] * sa;
	if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	  goto err2;
      }
    }
    
    if (AddHull(tail, vl) < 0)
      goto err2;
    LP_VertexList_Free(vl);
  }
  
  return 0;
  
 err2:
  LP_VertexList_Free(vl);
 err:
  return -1;
}

struct lp_vl_list *LCW_SolidConvexDecomp(const struct lcw_solid *solid, float dtol, float ptol) {
  struct lp_vl_list *head = NULL, **tail = &head;
  struct lcw_list *decomp, *cur;
  struct lcw_wire *poly;
  float max_scale;

  max_scale = MaxScale(solid);
  if ((decomp = LCW_ConvexDecomp(solid->wire, dtol / max_scale)) == NULL)
    goto err;
  
  for (cur = decomp; cur != NULL; cur = cur->next) {
    if ((poly = LCW_ToPolygon(cur->wire, ptol / max_scale)) == NULL)
      goto err2;

    if (solid->type == lcw_extrude) {
      if (AddExtrude(poly, solid->param, &tail) < 0)
	goto err3;
    } else if (AddRevolve(poly, solid->param, &tail, dtol, ptol) < 0) {
      goto err3;
    }
    
    LCW_Free(poly);
  }

  LCW_ListFree(decomp, 1);
  return head;

 err3:
  LCW_Free(poly);
 err2:
  LCW_ListFree(decomp, 1);
 err:
  LP_VertexList_ListFree(head);
  return NULL;
}
