/*
 * sourcepos.cpp
 */

#include "geo2d.h"

char geo2d::sourcepos(int isrc, int *ix, int *iz)
/*
 * Locates position of source no. isrc within the grid. Returns 
 *  'i' if source is located within cell,
 *  'h' or 'v' if source is located on a horizontal or vertical
 *             edge between two cells, 
 *  'n' if source is located on a grid node of 4 adjacent cells.
 * Indices of corresponding conductivity cell are returned in 
 * ix and iz.
 */
{
  double xs = xsrc[isrc], zs = zsrc[isrc];
  int i;
  char pos;
  
  // first, search x-direction
  i = locate(xgrd, nxgrd, xs);
  i = max(i,0);
  i = min(i,nxgrd-2);
  if (fabs(xs-xgrd[i+1]) < 1.0e-15){    // right edge of cell
    *ix = config.NoEquiX + config.NoProlX + (config.NoRefineX + 1) * (i + 1);
    pos = 'v';
  }
  else {
    if (fabs(xs-xgrd[i]) < 1.0e-15){  // left edge of cell
      *ix = config.NoEquiX + config.NoProlX + (config.NoRefineX + 1) * i;
      pos = 'v';
    }
    else {                // interiour of cell
      *ix = config.NoEquiX + config.NoProlX + (config.NoRefineX + 1) * i;
      pos = 'i';
    }
  }

  // second, search z-direction
  i = locate(zgrd, nzgrd, zs);
  i = max(i,0);
  i = min(i,nzgrd-2);
  if (fabs(zs-zgrd[i]) < 1.0e-15){    // top edge of cell
    *iz = (config.NoRefineZ + 1) * i;
    if (pos == 'v'){
      pos = 'n';
    }
    else {
      pos = 'h';
    }
  }
  else {
    if (fabs(zs-zgrd[i+1]) < 1.0e-15){  // bottom edge of cell
      *iz = (config.NoRefineZ + 1) * (i + 1);
      if (pos == 'v') {
        pos = 'n';
      }
      else {
        pos = 'h';
      }
    }
    else {      // vertical edge or interiour of cell
      *iz = (config.NoRefineZ + 1) * i;
    }
  }
  return(pos);
}

