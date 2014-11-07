#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NUMATOMS  13500  // 单原子数 Ni  Ti
#define CUTOFF    4.5    // 截断距离
#define STEP	100000     // 步长
#define FINALSTEP 500000
void quicksort(float c[], int low, int high)
{
    int i = low;
    int j = high;
    float temp = c[i];
    if(low < high)
    {
        while(i < j)
        {
            while((c[j] >= temp) && (i < j))
            {
                j--;
            }
            c[i] = c[j];
            while((c[i] <= temp) && (i < j))
            {
                i++;
            }
            c[j] = c[i];
        }
    c[i] = temp;
    quicksort(c, low, i-1);
    quicksort(c, j+1, high);
    }

}


int main(int argc, char* argv[])
{
    int step;

    for(step = 100000; step <= FINALSTEP; step = step + STEP)
    {
	   		int a, b, d;
    		int e[6] = {0,0,0,0,0,0};  // 用于统计严重偏离平衡位置的原子
    		float esum;
    		float f0, f1, f2;
    		float x, y, z;
    		float rx, ry, rz;
    		float x_cut, y_cut, z_cut;
    		float c[200];
    		float xx, yy, zz, xy, yz, xz;
    		char  buf1[256], buf2[256];
    		float Ni[NUMATOMS*3];
    		float Ti[NUMATOMS*3];

        FILE *fin;
        FILE *pin;

        sprintf(buf1, "%d.cfg", step);
        if((fin=fopen(buf1, "r"))==NULL)  exit(0);
        fscanf(fin, "%*[^\n]%*c");
        fscanf(fin, "%*[^\n]%*c");
        fscanf(fin, "%*7c\n%*c%f%*[^\n]%*c", &xx);
        fscanf(fin, "%*[^\n]%*c");
        fscanf(fin, "%*[^\n]%*c");
        fscanf(fin, "%*7c\n%*c%f%*[^\n]%*c", &xy);
        fscanf(fin, "%*7c\n%*c%f%*[^\n]%*c", &yy);
        fscanf(fin, "%*[^\n]%*c");
        fscanf(fin, "%*7c\n%*c%f%*[^\n]%*c", &xz);
        fscanf(fin, "%*7c\n%*c%f%*[^\n]%*c", &yz);
        fscanf(fin, "%*7c\n%*c%f%*[^\n]%*c", &zz);

        for(a = 0; a < 4; a++)  fscanf(fin, "%*[^\n]%*c");
        for(a = 0; a < NUMATOMS; a++)
        {
            fscanf(fin, "%f\n%f\n%f%*[^\n]%*c",&Ni[a*3], &Ni[a*3+1], &Ni[a*3+2]);
        }
        for(a = 0; a < 2; a++)  fscanf(fin,"%*c%*[^\n]%*c");
        for(a = 0; a < NUMATOMS; a++)
        {
            fscanf(fin, "%f\n%f\n%f%*[^\n]%*c", &Ti[a*3], &Ti[a*3+1], &Ti[a*3+2]);
        }
        fclose(fin);

        x_cut = xx - CUTOFF;
        y_cut = yy - CUTOFF;
        z_cut = zz - CUTOFF;
        sprintf(buf2, "op_%d.txt", step);
        if((pin = fopen(buf2, "w")) == NULL)  exit(0);
        for(a = 0; a < NUMATOMS; a++)
        {
          	d = 0;
          	for(b=0; b < NUMATOMS; b++)
          	{
              	rx = (Ti[b*3] - Ni[a*3]);
              	ry = (Ti[b*3+1] - Ni[a*3+1]);
              	rz = (Ti[b*3+2] - Ni[a*3+2]);
              	x = fabs(rx * xx + ry * xy + rz * xz);
              	if(x >= x_cut) x = xx - x;
								else if ( x <= CUTOFF) x = x;
              	else continue;
              	y = fabs(ry * yy + rz * yz);
              	if(y >= y_cut)y = yy - y;
								else if ( y <= CUTOFF) y = y;
              	else continue;
              	z = fabs(rz*zz);
              	if(z >= z_cut) z = zz - z;
								else if ( z <= CUTOFF) z = z;
              	else continue;
              	c[d] = sqrt(x*x + y*y + z*z);
              	d++;
          	}
          	quicksort(c, 0, d-1);
          	f0 = (c[0] + c[1] + c[6] + c[7]) / 2.0;
          	f1 = (c[2] + c[5] + c[3] + c[4]) / 2.0;
          	f2 = (f1 * (5.22 + 5.37) - f0 * (5.22 + 5.02)) / (5.22 * (5.02 - 5.37));
						//阶梯函数。。 
          	if(f2 <- 1.5)
						{
              	f2 = -2.0;
              	e[0]++;
          	}
        		else if(f2 >= -1.5 && f2 <- 0.2)
          	{
              	f2 = -1.0;
              	e[1]++;
          	}
        		else if(f2 >= -0.2 && f2 < 0.1)
          	{
              	f2 = 0.2;
              	e[2]++;
          	}
        		else if(f2 >= 0.1 && f2 < 1.8115)
          	{
              	f2 = 1.0;
              	e[3]++;
          	}
          	else if(f2 >= 1.8115 && f2 < 2.9)
          	{
              	f2 = 2.6133;
              	e[4]++;
          	}
          	else
          	{
              	f2 = 3.5;
              	e[5]++;
          	}
          	fprintf(pin, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.4f\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], f2);
        }
        fprintf(pin, "\n");
        fprintf(pin, "\n");

        for(a = 0; a < NUMATOMS; a++)
        {
            d = 0;
            for(b = 0; b < NUMATOMS; b++)
            {
                rx = (Ni[b*3] - Ti[a*3]);
                ry = (Ni[b*3+1] - Ti[a*3+1]);
                rz = (Ni[b*3+2] - Ti[a*3+2]);
								if(x >= x_cut) x = xx - x;
								else if ( x <= CUTOFF) x = x;
              	else continue;
              	y = fabs(ry * yy + rz * yz);
              	if(y >= y_cut)y = yy - y;
								else if ( y <= CUTOFF) y = y;
              	else continue;
              	z = fabs(rz*zz);
              	if(z >= z_cut) z = zz - z;
								else if ( z <= CUTOFF) z = z;
              	else continue;
               c[d] = sqrt(x*x + y*y + z*z);
                d++;
            }
            quicksort(c, 0, d-1);
            f0 = (c[0] + c[7] + c[1] + c[6])/2.0;
            f1 = (c[2] + c[5] + c[3] + c[4])/2.0;
            f2 = (f1 * (5.22 + 5.37) - f0 * (5.22 + 5.02)) / (5.22* (5.02 - 5.37));
            if(f2 <- 1.5)
	  				{
                f2 = -2.0;
                e[0]++;
            }
          	else if(f2 >= -1.5 && f2 <- 0.2)
            {
                f2 = -1.0;
                e[1]++;
            }
          	else if(f2 >= -0.2 && f2 < 0.1)
            {
                f2 = 0.2;
                e[2]++;
            }
          	else if(f2 >= 0.1 && f2 < 1.8115)
            {
                f2 = 1.0;
                e[3]++;
            }
            else if(f2 >= 1.8115 && f2 < 2.9)
            {
                f2 = 2.6133;
                e[4]++;
            }
            else
            {
                f2 = 3.5;
                e[5]++;
            }
            fprintf(pin,"%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.4f\n",c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],f2);
        }
        esum = (e[0] + e[1] + e[2] + e[3] + e[4] + e[5])*1.0;
        fprintf(pin,"B2: %.4f B19: %.4f B19': %.4f BCO: %.4f defect: %.4f\n",e[1]/esum, e[2]/esum, e[3]/esum, e[4]/esum, (e[0]+e[5])/esum);
        fclose(pin);
   }
  return 0;
}



