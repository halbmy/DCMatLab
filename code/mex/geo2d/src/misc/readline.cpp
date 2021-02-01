#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>


int readline(ifstream *infile, double **x)
{
   char *string1, *string2;
	double *d, *t;
	int i, j, n;
	const int lstr = 64, inc = 16;

	string1 = new char [lstr];
	string2 = string1;
	n = inc;
	d = new double [n];

	*string1 = '\0';
	i = 1;
	while (infile->peek() != '\n')
	{
		i--;
		infile->get(string2,lstr-strlen(string1));
		if (i == n)
		{
			t = d;
			n += inc;
			d = new double [n];
			for (j = 0; j < i; j++) d[j] = t[j];
			delete t;
		}
		d[i++] = strtod(string1, &string2);
      while (*string2 == ' ' || *string2 == '\t') // eat white characters
         string2++;
		while (strlen(string2))
		{
			strcpy(string1, string2);
			if (i == n)
			{
				t = d;
				n += inc;
				d = new double [n];
				for (j = 0; j < i; j++) d[j] = t[j];
				delete t;
			}
			d[i++] = strtod(string1, &string2);
         while (*string2 == ' ' || *string2 == '\t') // eat white characters
            string2++;
		}
	}
	infile->ignore(); // skip eol

	t = d;
	d = new double [i];
	for (j = 0; j < i; j++) d[j] = t[j];
	delete t;
	delete string1;

	*x = d;
	return(i);
}
