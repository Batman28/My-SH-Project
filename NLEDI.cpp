#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <math.h>
#include "memalloc.h"
#include "NLEDI.h"
//#include "image.h"
/* author: Xinfeng Zhang, 2009PCM paper*/
double thrld=8;
int imbilinear(double* imLR,double* imHR, int oldwidth, int oldheight, int width, int height)
{
	int i,j,x,y;
	double Xstep,Ystep,u,v;
	double** image2D;
	get_mem2Ddouble(&image2D,oldheight,oldwidth);
	Xstep=(double)oldheight/(double)height;
    Ystep=(double)oldwidth/(double)width;//compute the interpolating step
	for (j=0;j<oldheight;j++)//translate into 2D array
	{
		for (i=0;i<oldwidth;i++)
		{
			image2D[j][i]=imLR[j*oldwidth+i];
		}
	}
	for (j=0;j<height;j++)
	{
		for (i=0;i<width;i++)
		{
			u=j*Xstep;//row number of the new picture current point
			v=i*Ystep;//column number of the new picture current point
			x=(int)u;//row number of the LR image point  
			y=(int)v;//column number of the LR image point
			imHR[j*width+i]=image2D[min(x,oldheight-1)][min(y,oldwidth-1)]*(1+y-v)*(1+x-u)+image2D[min(x,oldheight-1)][min(y+1,oldwidth-1)]*(v-y)*(1+x-u)+image2D[min(x+1,oldheight-1)][min(y,oldwidth-1)]*(1+y-v)*(u-x)+image2D[min(x+1,oldheight-1)][min(y+1,oldwidth-1)]*(v-y)*(u-x);
		}
	}	
	free_mem2Ddouble(image2D);
	return 0;
}
int NonLocalEDI(double* imLR,double* imHR_out,int xs_LR,int ys_LR, int xs,int ys)
{
	int j,i,k,jj,ii,u,v,uu,vv,m,n,win_rad,swin_rad,totalsample,icheck;
// 	int cen_pos;
	int nx[64]={-2,2,-2,2},ny[64]={-2,-2,2,2},my[4]={0,0,1,-1},mx[4]={-1,1,0,0};
	double **C,**R,*P,maxP,**nExX,**nEXX,**C_T,*A,*b,sum;
	double var,avg,is;
	double **center,**sample,**imHR;
	imbilinear(imLR,imHR_out,xs_LR,ys_LR,xs,ys);
	get_mem2Ddouble(&imHR,ys,xs);
	for (j = 0; j < ys; j ++)
	{
		for (i = 0; i< xs; i++)
		{
			imHR[j][i] = imHR_out[j*xs + i];
		}
	}
	win_rad=4;swin_rad=3;
	get_mem2Ddouble(&C,2*win_rad*2*win_rad,4);
	get_mem2Ddouble(&C_T,4,2*win_rad*2*win_rad);
	get_mem2Ddouble(&R,2*win_rad*2*win_rad,1);
	get_mem1Ddouble(&P,2*win_rad*2*win_rad);
	get_mem2Ddouble(&nExX,4,1);
	get_mem2Ddouble(&nEXX,4,4);
	get_mem1Ddouble(&A,4*4);
	get_mem1Ddouble(&b,4*1);
	get_mem2Ddouble(&center,2*swin_rad+1,2*swin_rad+1);
	get_mem2Ddouble(&sample,2*swin_rad+1,2*swin_rad+1);
	
	//interpolate the middle point
	for (j=1;j<ys-1;j=j+2)
	{

		for (i=1;i<xs-1;i=i+2)
		{
			avg=(imHR[j-1][i-1]+imHR[j-1][i+1]+imHR[j+1][i-1]+imHR[j+1][i+1])/4;
			var=sqrt(pow(imHR[j-1][i-1]-avg,2)+pow(imHR[j-1][i+1]-avg,2)+pow(imHR[j+1][i-1]-avg,2)+pow(imHR[j+1][i+1]-avg,2));
			if (var>thrld)
			{
			
			//////////////////////////////////////////////////////////////////////////
			for (jj=-swin_rad;jj<=swin_rad;jj++)
			{
				u=iClip3(0,ys-1,j+jj);
				for (ii=-swin_rad;ii<=swin_rad;ii++)
				{					
					v=iClip3(0,xs-1,i+ii);
					center[jj+swin_rad][ii+swin_rad]=imHR[u][v];
				}
			}
			//////////////////////////////////////////////////////////////////////////
			k=0;maxP=0;
			for (jj=-2*win_rad+2;jj<=2*win_rad;jj=jj+2)
			{
				u=iClip3(0,ys-1,j-1+jj);
				for (ii=-2*win_rad+2;ii<=2*win_rad;ii=ii+2)
				{
					v=iClip3(0,xs-1,i-1+ii);//u,v here are the coordinate of sample
					//*************pixel in similarity window
					for (uu=-2*swin_rad;uu<=2*swin_rad;uu=uu+2)
					{
						m=iClip3(0,ys-1,u+uu);
						for (vv=-2*swin_rad;vv<=2*swin_rad;vv=vv+2)
						{					
							n=iClip3(0,xs-1,v+vv);
							sample[(uu+2*swin_rad)/2][(vv+2*swin_rad)/2]=imHR[m][n];
						}
					}
					//***************

					for (uu=0;uu<4;uu++)
					{
						C[k][uu]=imHR[iClip3(0,ys-1,u+ny[uu])][iClip3(0,xs-1,v+nx[uu])];
					}
					R[k][0]=imHR[iClip3(0,ys-1,u)][iClip3(0,xs-1,v)];
					ComputWeight(center,sample, swin_rad, k,P);
					k++;

				}
			}
// 			P[cen_pos]=maxP;
			totalsample=k;
			/////////////////////////////normalization///////////////////////////////////////
			sum=0;
			for (jj=0;jj<totalsample;jj++)
			{
				sum+=P[jj];
			}
			for (jj=0;jj<totalsample;jj++)
			{
				P[jj]=P[jj]/sum;
			}
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////
			for (jj=0;jj<totalsample;jj++)
			{
				R[jj][0]=R[jj][0]*P[jj];
			}
			MatrixTrans_NLEDI(C,totalsample,4,C_T);
			MatrixMul_NLEDI(C_T,R,nExX,4,totalsample,1);
			for (jj=0;jj<4;jj++)
			{
				for (ii=0;ii<totalsample;ii++)
				{
					C_T[jj][ii]=C_T[jj][ii]*P[ii];
				}

			}
			MatrixMul_NLEDI(C_T,C,nEXX,4,totalsample,4);
			for (jj=0;jj<4;jj++)
			{
				for (ii=0;ii<4;ii++)
				{
					A[jj*4+ii]=nEXX[jj][ii];
				}
				b[jj]=nExX[jj][0];
			}
			icheck=agaus(A,b,4);
			
			
			if (icheck!=0)
			{
				is=b[0]*imHR[j-1][i-1]+b[1]*imHR[j-1][i+1]+b[2]*imHR[j+1][i-1]+b[3]*imHR[j+1][i+1];
				if (is>0&&is<255)
				{
					imHR[j][i]=is;
				}
			}
			
			}
		}
	}


	//interpolate vertical edge point
	k=0;
	for (j=-(win_rad-1);j<=win_rad;j++)
	{
		for (i=-(win_rad-1);i<=win_rad;i++)
		{
			nx[k]=j-i;ny[k]=j+i-1;k++;
		}
	}

	for (j=1;j<ys-1;j=j+2)
	{
		for (i=2;i<xs-1;i=i+2)
		{
			avg=(imHR[j][i-1]+imHR[j][i+1]+imHR[j+1][i]+imHR[j-1][i])/4;
			var=sqrt(pow(imHR[j][i-1]-avg,2)+pow(imHR[j][i+1]-avg,2)+pow(imHR[j+1][i]-avg,2)+pow(imHR[j-1][i]-avg,2));
			if (var>thrld)
			{
			//////////////////////////////////////////////////////////////////////////
			for (jj=-swin_rad;jj<=swin_rad;jj++)
			{
				u=iClip3(0,ys-1,j+jj);
				for (ii=-swin_rad;ii<=swin_rad;ii++)
				{					
					v=iClip3(0,xs-1,i+ii);
					center[jj+swin_rad][ii+swin_rad]=imHR[u][v];
				}
			}
			//////////////////////////////////////////////////////////////////////////

			for (k=0;k<(2*win_rad*2*win_rad)-1;k++)
			{
				maxP=0;
				R[k][0]=imHR[iClip3(0,ys-1,j+ny[k])][iClip3(0,xs-1,i+nx[k])];
				for (jj=0;jj<4;jj++)
				{
					C[k][jj]=imHR[iClip3(0,ys-1,j+ny[k]+2*my[jj])][iClip3(0,xs-1,i+nx[k]+2*mx[jj])];
				}
				//*************pixel in similarity window
				for (uu=-2*swin_rad;uu<=2*swin_rad;uu=uu+2)
				{
					m=iClip3(0,ys-1,j+ny[k]+uu);
					for (vv=-2*swin_rad;vv<=2*swin_rad;vv=vv+2)
					{					
						n=iClip3(0,xs-1,i+nx[k]+vv);
						sample[(uu+2*swin_rad)/2][(vv+2*swin_rad)/2]=imHR[m][n];
					}
				}
				//***************
				ComputWeight(center,sample, swin_rad, k,P);
			
			}
			totalsample=k;
			/////////////////////////////normalization///////////////////////////////////////
			sum=0;
			for (jj=0;jj<totalsample;jj++)
			{
				sum+=P[jj];
			}
			for (jj=0;jj<totalsample;jj++)
			{
				P[jj]=P[jj]/sum;
			}
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////
			for (jj=0;jj<totalsample;jj++)
			{
				R[jj][0]=R[jj][0]*P[jj];
			}
			MatrixTrans_NLEDI(C,totalsample,4,C_T);
			MatrixMul_NLEDI(C_T,R,nExX,4,totalsample,1);
			for (jj=0;jj<4;jj++)
			{
				for (ii=0;ii<totalsample;ii++)
				{
					C_T[jj][ii]=C_T[jj][ii]*P[ii];
				}

			}
			MatrixMul_NLEDI(C_T,C,nEXX,4,totalsample,4);
			for (jj=0;jj<4;jj++)
			{
				for (ii=0;ii<4;ii++)
				{
					A[jj*4+ii]=nEXX[jj][ii];
				}
				b[jj]=nExX[jj][0];
			}
			icheck=agaus(A,b,4);


			if (icheck!=0)
			{
				is=b[0]*imHR[j][i-1]+b[1]*imHR[j][i+1]+b[2]*imHR[j+1][i]+b[3]*imHR[j-1][i];
				if (is>0&&is<255)
				{
					imHR[j][i]=is;
				}
			}
		}
		}
	}



	//interpolate horizontal edge point
	for (j=2;j<ys-1;j=j+2)
	{
		for (i=1;i<xs-1;i=i+2)
		{
			avg=(imHR[j][i-1]+imHR[j][i+1]+imHR[j+1][i]+imHR[j-1][i])/4;
			var=sqrt(pow(imHR[j][i-1]-avg,2)+pow(imHR[j][i+1]-avg,2)+pow(imHR[j+1][i]-avg,2)+pow(imHR[j-1][i]-avg,2));
			if (var>thrld)
			{
			
			//////////////////////////////////////////////////////////////////////////
			for (jj=-swin_rad;jj<=swin_rad;jj++)
			{
				u=iClip3(0,ys-1,j+jj);
				for (ii=-swin_rad;ii<=swin_rad;ii++)
				{					
					v=iClip3(0,xs-1,i+ii);
					center[jj+swin_rad][ii+swin_rad]=imHR[u][v];
				}
			}
			//////////////////////////////////////////////////////////////////////////

			for (k=0;k<(2*win_rad*2*win_rad)-1;k++)
			{
				maxP=0;
				R[k][0]=imHR[iClip3(0,ys-1,j+ny[k])][iClip3(0,xs-1,i+nx[k])];
				for (jj=0;jj<4;jj++)
				{
					C[k][jj]=imHR[iClip3(0,ys-1,j+ny[k]+2*my[jj])][iClip3(0,xs-1,i+nx[k]+2*mx[jj])];
				}
				//*************pixel in similarity window
				for (uu=-2*swin_rad;uu<=2*swin_rad;uu=uu+2)
				{
					m=iClip3(0,ys-1,j+ny[k]+uu);
					for (vv=-2*swin_rad;vv<=2*swin_rad;vv=vv+2)
					{					
						n=iClip3(0,xs-1,i+nx[k]+vv);
						sample[(uu+2*swin_rad)/2][(vv+2*swin_rad)/2]=imHR[m][n];
					}
				}
				//***************
				ComputWeight(center,sample, swin_rad, k,P);
			
			}
			totalsample=k;
			/////////////////////////////normalization///////////////////////////////////////
			sum=0;
			for (jj=0;jj<totalsample;jj++)
			{
				sum+=P[jj];
			}
			for (jj=0;jj<totalsample;jj++)
			{
				P[jj]=P[jj]/sum;
			}
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////
			for (jj=0;jj<totalsample;jj++)
			{
				R[jj][0]=R[jj][0]*P[jj];
			}
			MatrixTrans_NLEDI(C,totalsample,4,C_T);
			MatrixMul_NLEDI(C_T,R,nExX,4,totalsample,1);
			for (jj=0;jj<4;jj++)
			{
				for (ii=0;ii<totalsample;ii++)
				{
					C_T[jj][ii]=C_T[jj][ii]*P[ii];
				}

			}
			MatrixMul_NLEDI(C_T,C,nEXX,4,totalsample,4);
			for (jj=0;jj<4;jj++)
			{
				for (ii=0;ii<4;ii++)
				{
					A[jj*4+ii]=nEXX[jj][ii];
				}
				b[jj]=nExX[jj][0];
			}
			icheck=agaus(A,b,4);

			
			if (icheck!=0)
			{
				is=b[0]*imHR[j][i-1]+b[1]*imHR[j][i+1]+b[2]*imHR[j+1][i]+b[3]*imHR[j-1][i];
				if (is>=0&&is<=255)
				{
					imHR[j][i]=is;
				}
			}
			}
		}
	}
	for (j = 0; j < ys; j ++)
	{
		for (i = 0; i< xs; i++)
		{
			imHR_out[j*xs + i] = imHR[j][i];
		}
	}
	free_mem2Ddouble(C);
	free_mem2Ddouble(C_T);
	free_mem2Ddouble(R);
	free_mem1Ddouble(P);
	free_mem2Ddouble(nExX);
	free_mem2Ddouble(nEXX);
	free_mem1Ddouble(A);
	free_mem1Ddouble(b);
	free_mem2Ddouble(center);
	free_mem2Ddouble(sample);
	free_mem2Ddouble(imHR);
	return 0;
}


/*static inline */int iClip3(int low, int high, int x)
{
  x = max(x, low);
  x = min(x, high);
  return x;
}


void ComputWeight(double **center,double **sample, int swin_rad,int k, double *P)
{
	int j,i,m;
	double d=0,value1;
	double **DistantKenerl;
	get_mem2Ddouble(&DistantKenerl,2*swin_rad+1,2*swin_rad+1);
	for (m=0;m<2*swin_rad+1;m++)
	{
		memset(DistantKenerl[m],0,2*swin_rad+1);
	}
	for (m=0;m<=swin_rad;m++)
	{
		value1=1/(double)((2*(m+1)+1)*(2*(m+1)+1));
		for (j=-m;j<=m;j++)
		{
			for (i=-m;i<=m;i++)
			{
				DistantKenerl[swin_rad-j][swin_rad-i]=DistantKenerl[swin_rad-j][swin_rad-i]+value1;
			}
		}
	}
	value1=0;
	for (j=0;j<2*swin_rad+1;j++)
	{
		for (i=0;i<2*swin_rad+1;i++)
		{
			value1+=DistantKenerl[j][i];	
		}
	}
	for (j=0;j<2*swin_rad+1;j++)
	{
		for (i=0;i<2*swin_rad+1;i++)
		{
			DistantKenerl[j][i]=DistantKenerl[j][i]/value1;	
		}
	}
	for (j=0;j<2*swin_rad+1;j++)
	{
		for (i=0;i<2*swin_rad+1;i++)
		{
			d+=(center[j][i]-sample[j][i])*DistantKenerl[j][i]*(center[j][i]-sample[j][i])*DistantKenerl[j][i]*0.01;
		}
	}
	d=exp(-sqrt(d));
	if (d<0.0000001)
	{
		P[k]=0.0;
	} 
	else
	{
		P[k]=d;
	}
	free_mem2Ddouble(DistantKenerl);
}

/**********************************
* C=A*B
*
*
**********************************/
void MatrixMul_NLEDI(double** A,double** B,double **C,int row,int cn,int col)
{
	int i,j,k;
	double d=0;
	for (j=0;j<row;j++)
	{
		for (i=0;i<col;i++)
		{
			d=0;
			for (k=0;k<cn;k++)
			{
				d+=A[j][k]*B[k][i];
			}
			C[j][i]=d;
		}
	}
}
/************************************************************************/
/* matrix transpose                                                      */
/************************************************************************/
void MatrixTrans_NLEDI(double **A,int row, int col,double **A_T)
{
	int j,i;
	for (j=0;j<row;j++)
	{
		for (i=0;i<col;i++)
		{
			A_T[i][j]=A[j][i];
		}
	}
}
/************************************************************************/
/* solve the equation:ax=b                                               */
/************************************************************************/
 int agaus(double*a,double*b,int n)
  { int *js,l,k,i,j,is,p,q;
    double d,t;
    js=(int*)malloc(n*sizeof(int));
    l=1;
    for (k=0;k<=n-2;k++)
      { d=0.0;
        for (i=k;i<=n-1;i++)
          for (j=k;j<=n-1;j++)
            { t=fabs(a[i*n+j]);
              if (t>d) { d=t; js[k]=j; is=i;}
            }
        if (d+1.0==1.0) l=0;
        else
          { if (js[k]!=k)
              for (i=0;i<=n-1;i++)
                { p=i*n+k; q=i*n+js[k];
                  t=a[p]; a[p]=a[q]; a[q]=t;
                }
            if (is!=k)
              { for (j=k;j<=n-1;j++)
                  { p=k*n+j; q=is*n+j;
                    t=a[p]; a[p]=a[q]; a[q]=t;
                  }
                t=b[k]; b[k]=b[is]; b[is]=t;
              }
          }
        if (l==0)
          { free(js); //printf("fail\n");
            return(0);
          }
        d=a[k*n+k];
        for (j=k+1;j<=n-1;j++)
          { p=k*n+j; a[p]=a[p]/d;}
        b[k]=b[k]/d;
        for (i=k+1;i<=n-1;i++)
          { for (j=k+1;j<=n-1;j++)
              { p=i*n+j;
                a[p]=a[p]-a[i*n+k]*a[k*n+j];
              }
            b[i]=b[i]-a[i*n+k]*b[k];
          }
      }
    d=a[(n-1)*n+n-1];
    if (fabs(d)+1.0==1.0)
      { free(js); //printf("fail\n");
        return(0);
      }
    b[n-1]=b[n-1]/d;
    for (i=n-2;i>=0;i--)
      { t=0.0;
        for (j=i+1;j<=n-1;j++)
          t=t+a[i*n+j]*b[j];
        b[i]=b[i]-t;
      }
    js[n-1]=n-1;
    for (k=n-1;k>=0;k--)
      if (js[k]!=k)
        { t=b[k]; b[k]=b[js[k]]; b[js[k]]=t;}
    free(js);
    return(1);
  }