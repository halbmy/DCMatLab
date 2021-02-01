#include "dc2d.h"

/* constructor */
application::application()
{
   geo2dobj = NULL;
   DataVal = NULL;
   DataStd = NULL;
   DataSyn = NULL;
   NoData = 0;
   NXGrid = 0;
   NZGrid = 0;
   XGrid = NULL;
   ZGrid = NULL;
   Conductivity = NULL;
   IConductivity = NULL;
   obj1 = 0.0;
   obj2 = 0;
   block.n = 0;
   block.cnd = NULL;
   block.inb = NULL;
   block.nnb = NULL;
   block.flag = NULL;
   block.list = NULL;
   CXM = NULL;
   CXP = NULL;
   CZM = NULL;
   CZP = NULL;
}

/* destructor */
application::~application()
{
   if (geo2dobj != NULL)
      delete geo2dobj;
   if (Conductivity != NULL) 
      delete [] Conductivity;
   if (ZGrid != NULL) 
      delete [] ZGrid;
   if (XGrid != NULL) 
      delete [] XGrid;
   if (IConductivity != NULL) 
      delete [] IConductivity;
   if (DataSyn != NULL) 
      delete [] DataSyn;
   if (DataStd != NULL) 
      delete [] DataStd;
   if (DataVal != NULL) 
      delete [] DataVal;
   if (CXM != NULL) 
      delete [] CXM;
   if (CXP != NULL) 
      delete [] CXP;
   if (CZM != NULL) 
      delete [] CZM;
   if (CZP != NULL) 
      delete [] CZP;
   if (block.list != NULL) 
      delete [] block.list;
   if (block.flag != NULL) 
      delete [] block.flag;
   if (block.inb != NULL)
   {
      while (block.n > 0)
      {
         block.n--;
         if (block.inb[block.n] != NULL) 
            delete [] block.inb[block.n];
      }
      delete [] block.inb;
   }
   if (block.nnb != NULL) 
      delete [] block.nnb;
   if (block.cnd != NULL) 
      delete [] block.cnd;
}
