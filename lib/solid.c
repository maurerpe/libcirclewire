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
  LCW_Properties(&prop, wire);
  
  solid->type = type;
  switch (type) {
  case lcw_extrude:
    dx = -prop.center_of_mass[LCW_X];
    dy = -prop.center_of_mass[LCW_Y];
    
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
    else if (param[2] > -M_PI + 1e-6 && param[2] < M_PI - 1e-6)
      solid->param[2] = param[2];
    else
      goto err2;
    break;
    
  case lcw_revolve:
    dx = -(prop.center_of_mass[LCW_X] + prop.second_moment[LCW_XY] / (prop.area * prop.center_of_mass[LCW_Y]));
    
    if (prop.min[1] < 0)
      goto err2;
    
    if (num_params > 1)
      goto err2;

    if (num_params < 1)
      solid->param[0] = M_PI;
    else if (param[0] > 0 && param[0] <= M_PI + 1e-6)
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

struct lcw_solid *LCW_SolidCopy(struct lcw_solid *solid) {
  return LCW_Solid(LCW_SolidWire(solid),
		   LCW_SolidType(solid),
		   LCW_SolidParams(solid),
		   LCW_SolidNumParams(solid));
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
    j = (q + Ixxx) * sinf(2 * ang) / 2;
    k = (q + Ixxx + 2 * Ixyy) * ang;
    sin_a = sinf(ang);
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

static void Check_Rad(float *rad, const struct lcw_wire *wire, float ang, const float *center) {
  float ca, sa, c2[2], rad2d, dist, dd;
  
  ca = cosf(ang);
  sa = sinf(ang);
  c2[0] = center[0];
  c2[1] = center[1] * ca + center[2] * sa;
  rad2d = LCW_BoundingCircle(wire, c2);
  dist = center[2] * ca - center[1] * sa;
  dd = sqrtf(rad2d * rad2d + dist * dist);
  if (dd > *rad)
    *rad = dd;
}

float LCW_SolidBoundingSphere(const struct lcw_solid *solid, const float *center) {
  float h, scale, tan_sw, dz, c2[2], rad2d, rad, dd, ang, aa;
  
  if (solid->type == lcw_extrude) {
    h = 2 * solid->param[0];
    scale = solid->param[1];
    tan_sw = tan(solid->param[2]);

    dz = h - center[2];
    c2[0] = (center[0] - tan_sw * h) / scale;
    c2[1] =  center[1]               / scale;
    rad2d = LCW_BoundingCircle(solid->wire, c2) * scale;
    rad = sqrtf(rad2d * rad2d + dz + dz);
    
    dz = -h - center[2];
    c2[0] = (center[0] + tan_sw * h) / (2 - scale);
    c2[1] =  center[1]               / (2 - scale);
    rad2d = LCW_BoundingCircle(solid->wire, c2) * (2 - scale);
    dd = sqrtf(rad2d * rad2d + dz * dz);
    if (dd > rad)
      rad = dd;
  } else {
    ang = solid->param[0];

    rad = 0;
    Check_Rad(&rad, solid->wire, ang, center);
    Check_Rad(&rad, solid->wire, -ang, center);
    aa = atan2(center[2], center[1]);
    if (fabsf(aa) < ang)
      Check_Rad(&rad, solid->wire, aa, center);
  }

  return rad;
}

static int SamePt(const float *pt1, const float *pt2, float comb_tol) {
  return
    fabsf(pt1[0] - pt2[0]) < comb_tol &&
    fabsf(pt1[1] - pt2[1]) < comb_tol &&
    fabsf(pt1[2] - pt2[2]) < comb_tol;
}

static void Cross(float *result, const float *a, const float *b) {
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
}

static float PlaneNorm(float *norm, const float *p1, const float *p2, const float *p3) {
  float v1[3], v2[3], mag;
  
  v1[0] = p2[0] - p1[0];
  v1[1] = p2[1] - p1[1];
  v1[2] = p2[2] - p1[2];
  
  v2[0] = p3[0] - p2[0];
  v2[1] = p3[1] - p2[1];
  v2[2] = p3[2] - p2[2];
  
  Cross(norm, v1, v2);
  mag = sqrtf(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
  if (mag > 0) {
    norm[0] /= mag;
    norm[1] /= mag;
    norm[2] /= mag;
  }
  return mag;
}

static int AddTri(struct lp_vertex_list *vl, float *pt1, float *pt2, float *pt3, float comb_tol) {
  float norm[3];
  
  if (SamePt(pt1, pt2, comb_tol) ||
      SamePt(pt2, pt3, comb_tol) ||
      SamePt(pt1, pt3, comb_tol))
    return 0;
  
  if (PlaneNorm(norm, pt1, pt2, pt3) == 0)
    return 0;
  
  pt1[3] = pt2[3] = pt3[3] = norm[0];
  pt1[4] = pt2[4] = pt3[4] = norm[1];
  pt1[5] = pt2[5] = pt3[5] = norm[2];
  
  if (LP_VertexList_Add(vl, pt1) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(vl, pt2) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(vl, pt3) == UINT_MAX)
    return -1;

  return 0;
}

static int AddQuad(struct lp_vertex_list *vl, float *vert, float comb_tol) {
  if (AddTri(vl, &vert[1 * 8], &vert[0 * 8], &vert[2 * 8], comb_tol) < 0)
    return -1;
  if (AddTri(vl, &vert[2 * 8], &vert[3 * 8], &vert[1 * 8], comb_tol) < 0)
    return -1;
  
  return 0;
}

struct lp_vertex_list *LCW_SolidMesh(const struct lcw_solid *solid, float tol, float scale, struct lcw_list **stencil_out) {
  struct lcw_wire *poly;
  struct lcw_list *sten[4];
  struct lcw_seg *seg;
  struct lp_vertex_list *vl, *end;
  struct lcw_properties prop;
  float vert[4 * 8], *vbase, *vv, h, escale, tan_sw, pt[2];
  float tarc, parc, carc, arc, xmid, xoff, ymid, yoff, xsc, ysc, asc;
  float ang, csa, sna, ang1, ang2, cos_a1, sin_a1, cos_a2, sin_a2;
  unsigned int *ind;
  int vidx, off, snum;
  size_t fpv, idx, cnt, num, inum;
  
  if ((poly = LCW_ToPolygon(solid->wire, tol)) == NULL)
    goto err;
  if (LCW_WireLoop(poly) != LCW_NO_ERROR)
    goto err2;
  tarc = LCW_TotalArcLen(solid->wire);
  parc = LCW_TotalArcLen(poly);
  
  if ((end = LCW_Mesh(poly, tol, 1)) == NULL)
    goto err2;

  if ((vl = LP_VertexList_New(8, lp_pt_triangle)) == NULL)
    goto err3;

  memset(sten, 0, sizeof(sten));
  if (stencil_out) {
    for (idx = 0; idx < 4; idx++) {
      if ((sten[idx] = LCW_ListNew()) == NULL)
	goto err4;
      if (idx > 0)
	sten[idx-1]->next = sten[idx];
    }
  }
  
  LCW_Properties(&prop, solid->wire);
  
  if (solid->type == lcw_extrude) {
    h = solid->param[0];
    escale = solid->param[1];
    tan_sw = tan(solid->param[2]);

    fpv   = LP_VertexList_FloatsPerVert(end);
    inum  = LP_VertexList_NumInd(end);
    vbase = LP_VertexList_GetVert(end);
    ind   = LP_VertexList_GetInd(end);
    xmid  = (prop.max[0] + prop.min[0]) / 2;
    xoff  = 2.15 * h + 1.5 * (prop.max[0] - prop.min[0]);
    ymid  = (prop.max[1] + prop.min[0]) / 2;
    yoff  = 0.55 * tarc;
    xsc   = 1.0f / (2 * (prop.max[0] - prop.min[0]) + 2.2 * h);
    ysc   = 1.0f / (1.1 * tarc);
    vidx = 0;
    while (inum-- > 0) {
      off = 2 * (1 - vidx); // Switch ind 0 and 2
      vv = vbase + ind[off] * fpv;
      vert[0] = scale * (vv[0] * escale + tan_sw * h);
      vert[1] = scale * (vv[1] * escale);
      vert[2] = scale * h;
      vert[3] = 0;
      vert[4] = 0;
      vert[5] = 1;
      vert[6] = ((xmid - vv[0]) + xoff) * xsc;
      vert[7] = ((ymid - vv[1]) + yoff) * ysc;
      if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	goto err4;
      ind++;
      vidx = (vidx + 1) % 3;
    }
    
    if (stencil_out) {
      if ((sten[0]->wire = LCW_New(solid->wire->comb_tol)) == NULL)
	goto err4;
      pt[0] = 0;
      pt[1] = 0;
      Reset(sten[0]->wire, pt, 0);
      pt[1] = 1.1 * tarc;
      if (To(sten[0]->wire, pt, 0, 0) < 0)
	goto err4;
      pt[0] = 2.2 * h + 2 * (prop.max[0] - prop.min[0]);
      if (To(sten[0]->wire, pt, 0, 0) < 0)
	goto err4;
      pt[1] = 0;
      if (To(sten[0]->wire, pt, 0, 0) < 0)
	goto err4;
      if (LCW_Close(sten[0]->wire) < 0)
	goto err4;
      
      if ((sten[1]->wire = LCW_Copy(solid->wire)) == NULL)
	goto err4;
      LCW_MirrorY(sten[1]->wire);
      if (LCW_MirrorX(sten[1]->wire))
	goto err4;
      LCW_Translate(sten[1]->wire, xmid + xoff, ymid + yoff);
    }
    inum = LP_VertexList_NumInd(end);
    ind  = LP_VertexList_GetInd(end);
    xoff = 0.05 * h + 0.5 * (prop.max[0] - prop.min[0]);
    while (inum-- > 0) {
      vv = vbase + *ind * fpv;
      vert[0] = scale * (vv[0] * (2 - escale) - tan_sw * h);
      vert[1] = scale * (vv[1] * (2 - escale));
      vert[2] = scale * -h;
      vert[3] = 0;
      vert[4] = 0;
      vert[5] = -1;
      vert[6] = ((vv[0] - xmid) + xoff) * xsc;
      vert[7] = ((ymid - vv[1]) + yoff) * ysc;
      if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	goto err4;
      ind++;
    }
    if (stencil_out) {
      if ((sten[3]->wire = LCW_Copy(solid->wire)) == NULL)
	goto err4;
      LCW_MirrorY(sten[3]->wire);
      LCW_Translate(sten[3]->wire, -xmid + xoff, ymid + yoff);
    }
    
    xoff = 1.1 * h + (prop.max[0] - prop.min[0]);
    ymid = 0.5 * tarc;
    yoff = 0.55 * tarc;
    asc  = tarc / parc;
    carc = -0.5 * tarc;
    seg  = &poly->seg[0][1];
    for (idx = 1; idx < poly->num_seg[0]; idx++, seg++) {
      arc = ArcLen(seg[-1].pt, seg[0].pt, seg[0].alpha) * asc;
      vert[ 0] = scale * (seg[-1].pt[0] * escale + tan_sw * h);
      vert[ 1] = scale * (seg[-1].pt[1] * escale);
      vert[ 2] = scale * h;
      vert[ 3] = 0;
      vert[ 4] = 0;
      vert[ 5] = 0;
      vert[ 6] = (h + xoff) * xsc;
      vert[ 7] = ((ymid - carc) + yoff) * ysc;
      vert[ 8] = scale * (seg[ 0].pt[0] * escale + tan_sw * h);
      vert[ 9] = scale * (seg[ 0].pt[1] * escale);
      vert[10] = scale * h;
      vert[11] = 0;
      vert[12] = 0;
      vert[13] = 0;
      vert[14] = (h + xoff) * xsc;
      vert[15] = ((ymid - (carc + arc)) + yoff) * ysc;
      vert[16] = scale * (seg[-1].pt[0] * (2 - escale) - tan_sw * h);
      vert[17] = scale * (seg[-1].pt[1] * (2 - escale));
      vert[18] = scale * -h;
      vert[19] = 0;
      vert[20] = 0;
      vert[21] = 0;
      vert[22] = (-h + xoff) * xsc;
      vert[23] = ((ymid - carc) + yoff) * ysc;
      vert[24] = scale * (seg[ 0].pt[0] * (2 - escale) - tan_sw * h);
      vert[25] = scale * (seg[ 0].pt[1] * (2 - escale));
      vert[26] = scale * -h;
      vert[27] = 0;
      vert[28] = 0;
      vert[29] = 0;
      vert[30] = (-h + xoff) * xsc;
      vert[31] = ((ymid - (carc + arc)) + yoff) * ysc;
      if (AddQuad(vl, vert, poly->comb_tol * scale) < 0)
	goto err4;
      carc += arc;
    }
    if (stencil_out) {
      if ((sten[2]->wire = LCW_New(solid->wire->comb_tol)) == NULL)
	goto err4;
      pt[0] = -h + xoff;
      pt[1] = 1.05 * tarc;
      Reset(sten[2]->wire, pt, 0);
      pt[0] = h + xoff;
      pt[1] = 1.05 * tarc;
      if (To(sten[2]->wire, pt, 0, 0) < 0)
	goto err4;
      pt[0] = h + xoff;
      pt[1] = 0.05 * tarc;
      if (To(sten[2]->wire, pt, 0, 0) < 0)
	goto err4;
      pt[0] = -h + xoff;
      pt[1] = 0.05 * tarc;
      if (To(sten[2]->wire, pt, 0, 0) < 0)
	goto err4;
      LCW_Close(sten[2]->wire);
    }
  } else {
    ang = solid->param[0];
    h   = ang * prop.max[1];
    ysc = 1.0f / (1.1 * tarc);
    if (ang < M_PI - poly->comb_tol / prop.max[1]) {
      fpv   = LP_VertexList_FloatsPerVert(end);
      inum  = LP_VertexList_NumInd(end);
      ind   = LP_VertexList_GetInd(end);
      vbase = LP_VertexList_GetVert(end);
      csa   = cosf(ang);
      sna   = sinf(ang);
      xmid  = (prop.max[0] + prop.min[0]) / 2;
      xoff  = 2.15 * h + 1.5 * (prop.max[0] - prop.min[0]);
      ymid  = (prop.max[1] + prop.min[0]) / 2;
      yoff  = 0.55 * tarc;
      xsc   = 1.0f / (2 * (prop.max[0] - prop.min[0]) + 2.2 * h);
      vidx = 0;
      while (inum-- > 0) {
	off = 2 * (1 - vidx); // Switch ind 0 and 2
	vv = vbase + ind[off] * fpv;
	vert[0] = scale * vv[0];
	vert[1] = scale * vv[1] * csa;
	vert[2] = scale * vv[1] * sna;
	vert[3] = 0;
	vert[4] = -sna;
	vert[5] = csa;
	vert[6] = ((xmid - vv[0]) + xoff) * xsc;
	vert[7] = ((ymid - vv[1]) + yoff) * ysc;
	if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	  goto err4;
	ind++;
	vidx = (vidx + 1) % 3;
      }
      if (stencil_out) {
	if ((sten[0]->wire = LCW_New(solid->wire->comb_tol)) == NULL)
	  goto err4;
	pt[0] = 0;
	pt[1] = 0;
	Reset(sten[0]->wire, pt, 0);
	pt[1] = 1.1 * tarc;
	if (To(sten[0]->wire, pt, 0, 0) < 0)
	  goto err4;
	pt[0] = 2.2 * h + 2 * (prop.max[0] - prop.min[0]);
	if (To(sten[0]->wire, pt, 0, 0) < 0)
	  goto err4;
	pt[1] = 0;
	if (To(sten[0]->wire, pt, 0, 0) < 0)
	  goto err4;
	if (LCW_Close(sten[0]->wire) < 0)
	  goto err4;
	
	if ((sten[1]->wire = LCW_Copy(solid->wire)) == NULL)
	  goto err4;
	LCW_MirrorY(sten[1]->wire);
	if (LCW_MirrorX(sten[1]->wire))
	  goto err4;
	LCW_Translate(sten[1]->wire, xmid + xoff, ymid + yoff);
      }
      inum = LP_VertexList_NumInd(end);
      ind  = LP_VertexList_GetInd(end);
      xoff = 0.05 * h + 0.5 * (prop.max[0] - prop.min[0]);
      while (inum-- > 0) {
	vv = vbase + *ind * fpv;
	vert[0] = scale * vv[0];
	vert[1] = scale * vv[1] * csa;
	vert[2] = scale * vv[1] * -sna;
	vert[3] = 0;
	vert[4] = -sna;
	vert[5] = -csa;
	vert[6] = ((vv[0] - xmid) + xoff) * xsc;
	vert[7] = ((ymid - vv[1]) + yoff) * ysc;
	if (LP_VertexList_Add(vl, vert) == UINT_MAX)
	  goto err4;
	ind++;
      }
      if (stencil_out) {
	if ((sten[3]->wire = LCW_Copy(solid->wire)) == NULL)
	  goto err4;
	LCW_MirrorY(sten[3]->wire);
	LCW_Translate(sten[3]->wire, -xmid + xoff, ymid + yoff);
      }
      xoff = 1.1 * h + (prop.max[0] - prop.min[0]);
      snum = 2;
    } else {
      xoff = 1.05 * h;
      xsc  = 1.0f / (2.1 * h);
      snum = 1;
      if (stencil_out) {
	sten[1]->next = NULL;
	LCW_ListFree(sten[2], 1);
	if ((sten[0]->wire = LCW_New(solid->wire->comb_tol)) == NULL)
	  goto err4;
	pt[0] = 0;
	pt[1] = 0;
	Reset(sten[0]->wire, pt, 0);
	pt[1] = 1.1 * tarc;
	if (To(sten[0]->wire, pt, 0, 0) < 0)
	  goto err4;
	pt[0] = 2.1 * h;
	if (To(sten[0]->wire, pt, 0, 0) < 0)
	  goto err4;
	pt[1] = 0;
	if (To(sten[0]->wire, pt, 0, 0) < 0)
	  goto err4;
	if (LCW_Close(sten[0]->wire) < 0)
	  goto err4;
      }
    }

    num = RevRequiredSegs(2 * ang, prop.max[1], tol);
    ymid = 0.5 * tarc;
    yoff = 0.55 * tarc;
    asc  = tarc / parc;
    for (cnt = 0; cnt < num; cnt++) {
      ang1 = ((float) cnt) / num * 2 * ang - ang;
      ang2 = ((float) (cnt + 1)) / num * 2 * ang - ang;
      cos_a1 = cosf(ang1);
      sin_a1 = sinf(ang1);
      cos_a2 = cosf(ang2);
      sin_a2 = sinf(ang2);
      carc = -0.5 * tarc;
      seg = &poly->seg[0][1];
      for (idx = 1; idx < poly->num_seg[0]; idx++, seg++) {
	arc = ArcLen(seg[-1].pt, seg[0].pt, seg[0].alpha) * asc;
	vert[ 0] = scale * seg[-1].pt[0];
	vert[ 1] = scale * seg[-1].pt[1] * cos_a2;
	vert[ 2] = scale * seg[-1].pt[1] * sin_a2;
	vert[ 3] = 0;
	vert[ 4] = 0;
	vert[ 5] = 0;
	vert[ 6] = (ang2 * seg[-1].pt[1] + xoff) * xsc;
	vert[ 7] = ((ymid - carc) + yoff) * ysc;
	vert[ 8] = scale * seg[ 0].pt[0];
	vert[ 9] = scale * seg[ 0].pt[1] * cos_a2;
	vert[10] = scale * seg[ 0].pt[1] * sin_a2;
	vert[11] = 0;
	vert[12] = 0;
	vert[13] = 0;
	vert[14] = (ang2 * seg[ 0].pt[1] + xoff) * xsc;
	vert[15] = ((ymid - (carc + arc)) + yoff) * ysc;
	vert[16] = scale * seg[-1].pt[0];
	vert[17] = scale * seg[-1].pt[1] * cos_a1;
	vert[18] = scale * seg[-1].pt[1] * sin_a1;
	vert[19] = 0;
	vert[20] = 0;
	vert[21] = 0;
	vert[22] = (ang1 * seg[-1].pt[1] + xoff) * xsc;
	vert[23] = ((ymid - carc) + yoff) * ysc;
	vert[24] = scale * seg[ 0].pt[0];
	vert[25] = scale * seg[ 0].pt[1] * cos_a1;
	vert[26] = scale * seg[ 0].pt[1] * sin_a1;
	vert[27] = 0;
	vert[28] = 0;
	vert[29] = 0;
	vert[30] = (ang1 * seg[ 0].pt[1] + xoff) * xsc;
	vert[31] = ((ymid - (carc + arc)) + yoff) * ysc;
	if (AddQuad(vl, vert, poly->comb_tol * scale) < 0)
	  goto err4;
	carc += arc;
      }
    }
    if (stencil_out) {
      if ((sten[snum]->wire = LCW_New(solid->wire->comb_tol)) == NULL)
	goto err4;
      seg = &poly->seg[0][1];
      pt[0] = ang * seg[-1].pt[1] + xoff;
      pt[1] = 1.05 * tarc;
      Reset(sten[snum]->wire, pt, 0);
      carc = 0;
      for (idx = 1; idx < poly->num_seg[0]; idx++, seg++) {
	arc = ArcLen(seg[-1].pt, seg[0].pt, seg[0].alpha) * asc;
	carc += arc;
	pt[0] = ang * seg->pt[1] + xoff;
	pt[1] = (ymid - carc) + yoff;
	if (To(sten[snum]->wire, pt, 0, 0) < 0)
	  goto err4;
      }
      seg = &poly->seg[0][poly->num_seg[0] - 2];
      pt[0] = -ang * seg[1].pt[1] + xoff;
      pt[1] = 0.05 * tarc;
      if (To(sten[snum]->wire, pt, 0, 0) < 0)
	goto err4;
      carc = 0;
      for (idx = 1; idx < poly->num_seg[0]; idx++, seg--) {
	arc = ArcLen(seg[0].pt, seg[1].pt, seg[1].alpha) * asc;
	carc += arc;
	pt[0] = -ang * seg->pt[1] + xoff;
	pt[1] = (carc - ymid) + yoff;
	if (To(sten[snum]->wire, pt, 0, 0) < 0)
	  goto err4;
      }
      if (LCW_Close(sten[snum]->wire) < 0)
	goto err4;
    }
  }

  if (stencil_out)
    *stencil_out = sten[0];
  LP_VertexList_Free(end);
  LCW_Free(poly);
  return vl;
  
 err4:
  LCW_ListFree(sten[0], 1);
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

static int AddExtrude(struct lcw_wire *poly, const float *param, struct lp_vl_list ***tail, float scale) {
  struct lp_vertex_list *vl;
  const struct lcw_seg *seg;
  float h, escale, tan_sw, vert[3];
  size_t cnt, ns;

  if (LCW_WireLoop(poly) != LCW_NO_ERROR)
    goto err;
  
  if ((vl = LP_VertexList_New(3, lp_pt_point)) == NULL)
    goto err;
  h      = param[0];
  escale = param[1];
  tan_sw = tan(param[2]);

  seg = poly->seg[0];
  ns = poly->num_seg[0];
  for (cnt = 1; cnt < ns; cnt++, seg++) {
    vert[0] = scale * (seg->pt[0] * escale + tan_sw * h);
    vert[1] = scale * (seg->pt[1] * escale);
    vert[2] = scale * h;
    if (LP_VertexList_Add(vl, vert) == UINT_MAX)
      goto err2;
    vert[0] = scale * (seg->pt[0] * (2 - escale) - tan_sw * h);
    vert[1] = scale * (seg->pt[1] * (2 - escale));
    vert[2] = scale * -h;
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

static int AddRevolve(struct lcw_wire *poly, const float *param, struct lp_vl_list ***tail, float dtol, float ptol, float scale) {
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
      ca = cosf(aa);
      sa = sinf(aa);
      seg = poly->seg[0];
      for (cnt = 0; cnt < ns; cnt++, seg++) {
	vert[0] = scale * seg->pt[0];
	vert[1] = scale * seg->pt[1] * ca;
	vert[2] = scale * seg->pt[1] * sa;
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

struct lp_vl_list *LCW_SolidConvexDecomp(const struct lcw_solid *solid, float dtol, float ptol, float scale) {
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
      if (AddExtrude(poly, solid->param, &tail, scale) < 0)
	goto err3;
    } else if (AddRevolve(poly, solid->param, &tail, dtol, ptol, scale) < 0) {
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
