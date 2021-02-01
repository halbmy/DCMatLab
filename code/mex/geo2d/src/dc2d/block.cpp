#include "dc2d.h"
void heapsort(unsigned long n, double *ra);

void application::findblocks(int nx, int nz, int *C)
{
  int iblock, icell, jcell,
    iblockcells, nblockcells, ineighbours, nneighbours;
  int *blockcells, *neighbours, *blocknum, *blockval; 
  int i,j,k;
  int *itop, **top;
  
  blockcells = new int [nx*nz];
  neighbours = new int [nx*nz];
  blocknum = new int [nx*nz];
  blockval = new int [nx*nz];

  for (icell = 0; icell < nx*nz; icell++)
  {
    blocknum[icell] = 0;
  }
  iblock = 0;
  ineighbours = 0;
  nneighbours = 1;
  neighbours[0] = 0;
  while (ineighbours < nneighbours)
  {
    icell = neighbours[ineighbours++];
    if (!blocknum[icell])
    {
      iblock++;
      iblockcells = 0;
      nblockcells = 1;
      blockcells[0] = icell;
      while (iblockcells < nblockcells)
      {
        icell = blockcells[iblockcells++];
        blocknum[icell] = iblock;
        blockval[iblock-1] = C[icell];
        jcell = icell+1;
        if (jcell % nz)
        {
          if (!blocknum[jcell])
          {
            if (C[icell] == C[jcell])
            {
              blockcells[nblockcells++] = jcell;
              blocknum[jcell] = iblock;
            }
            else
            {
              neighbours[nneighbours++] = jcell;
            }
          }
        }
        jcell = icell+nz;
        if (jcell < nx*nz)
        {
          if (!blocknum[jcell])
          {
            if (C[icell] == C[jcell])
            {
              blockcells[nblockcells++] = jcell;
              blocknum[jcell] = iblock;
            }
            else
            {
              neighbours[nneighbours++] = jcell;
            }
          }
        }
        jcell = icell-1;
        if (icell % nz)
        {
          if (!blocknum[jcell])
          {
            if (C[icell] == C[jcell])
            {
              blockcells[nblockcells++] = jcell;
              blocknum[jcell] = iblock;
            }
            else
            {
              neighbours[nneighbours++] = jcell;
            }
          }
        }
        jcell = icell-nz;
        if (jcell >= 0)
        {
          if (!blocknum[jcell])
          {
            if (C[icell] == C[jcell])
            {
              blockcells[nblockcells++] = jcell;
              blocknum[jcell] = iblock;
            }
            else
            {
              neighbours[nneighbours++] = jcell;
            }
          }
        }
      }
    }
  }
  
  itop = new int [iblock];
  top = new int * [iblock];
  for (i = 0; i < iblock; i++)
  {
    itop[i] = 0;
    top[i] = new int [iblock];
  }
  for (icell = 0; icell < nx*nz; icell += nz)
  {
    for (jcell = 1; jcell < nz; jcell++)
    {
      i = blocknum[icell+jcell-1]-1;
      j = blocknum[icell+jcell]-1;
      if (i != j)
      {
        for (k = 0; k < itop[i]; k++)
        {
          if (top[i][k] == j) break;
        }
        if (k == itop[i])
        {
          top[i][itop[i]++] = j;
          top[j][itop[j]++] = i;
        }
      }
    }
  }
  for (icell = nz; icell < nx*nz; icell += nz)
  {
    for (jcell = 0; jcell < nz; jcell++)
    {
      i = blocknum[icell+jcell-nz]-1;
      j = blocknum[icell+jcell]-1;
      if (i != j)
      {
        for (k = 0; k < itop[i]; k++)
        {
          if (top[i][k] == j) break;
        }
        if (k == itop[i])
        {
          top[i][itop[i]++] = j;
          top[j][itop[j]++] = i;
        }
      }
    }
  }

  block.n = iblock;
  block.cnd = new int [iblock];
  block.nnb = new int [iblock];
  block.inb = new int * [iblock];
  for (i = 0; i < iblock; i++)
  {
    block.cnd[i] = blockval[i];
    block.nnb[i] = itop[i];
    block.inb[i] = new int [itop[i]];
    for (j = 0; j < itop[i]; j++)
    {
      block.inb[i][j] = top[i][j];
    }
    cout << endl;
    delete top[i];
  }
  block.flag = new int [iblock];
  block.list = new int [iblock];
  delete top;
  delete itop;
  delete blockval;
  delete blocknum;
  delete neighbours;
  delete blockcells;
}

