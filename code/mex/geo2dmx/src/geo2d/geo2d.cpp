/*
 * geo2d.cpp
 */

#include "geo2d.h"

geo2d::geo2d()
/* Standard constructor */
{
   config.FacProlX = 1.0;
   config.FacProlZ = 1.0;
   config.NoProlX = 0;
   config.NoEquiX = 0;
   config.NoRefineX = 0;
   config.NoProlZ = 0;
   config.NoEquiZ = 0;
   config.NoRefineZ = 0;
   config.NoLegendre = 8;
   config.NoLaguerre = 8;
   config.ModProlong = 'b';
   config.Background = 's';
   config.Boundary = 'm';

   readconfiguration("geo2d.config");

   xcell = NULL; zcell = NULL;
   xnode = NULL; znode = NULL;
   cnorm = NULL; cdiff = NULL;
   cne = NULL; cns = NULL; cna = NULL;
   cae = NULL; cas = NULL; caa = NULL; cap = NULL;
   wavnum = NULL; weight = NULL;
   C = NULL;	b = NULL; phi = NULL;
   cback = NULL; iback = NULL;
   nx = 0; nz = 0;
   cx = 0; cz = 0;
   nw = 0;
}

geo2d::~geo2d()
/* Destructor */
{
   if (xcell != NULL) FREE(xcell);
   if (zcell != NULL) FREE(zcell);
   if (xnode != NULL) FREE(xnode);
   if (znode != NULL) FREE(znode);
   if (cnorm != NULL) FREE(cnorm);
   if (cdiff != NULL) FREE(cdiff);
   if (cback != NULL)
   {
      for (int i = 0; i < nsrc; i++)
         if (cback[i] != NULL) FREE(cback[i]);
      FREE(cback);
   }
   if (iback != NULL)
   {
      for (int i = 0; i < nsrc; i++)
         if (iback[i] != NULL) FREE(iback[i]);
      FREE(iback);
   }
   if (cne != NULL) FREE(cne);
   if (cns != NULL) FREE(cns);
   if (cna != NULL) FREE(cna);
   if (cae != NULL) FREE(cae);
   if (cas != NULL) FREE(cas);
   if (caa != NULL) FREE(caa);
   if (cap != NULL) FREE(cap);
   if (wavnum != NULL) FREE(wavnum);
   if (weight != NULL) FREE(weight);
   if (C != NULL) FREE(C);
   if (b != NULL) FREE(b);
   if (phi != NULL) FREE(phi);
}

