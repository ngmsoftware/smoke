#include <stdio.h>
#include <math.h>
#include "mex.h"

/* To compile:

mex -Dchar16_t=uint16_T exemplo.c

*/

#define NDIFUSE 10

/* Input */

#define	k_IN	prhs[0]
#define	pho_in_IN	prhs[1]
#define	nu_IN	prhs[2]
#define	vel_x_in_IN	prhs[3]
#define	vel_y_in_IN	prhs[4]
#define	dt_IN	prhs[5]
#define	srcPos_IN	prhs[6]

/* Output */

#define	pho_out_OUT	plhs[0]

#define pho_out(l,c) pho_out[(l)+pho_out_rows*(c)] // same
#define	vel_x_out_OUT	plhs[1]

#define vel_x_out(l,c) vel_x_out[(l)+vel_x_out_rows*(c)] // same
#define	vel_y_out_OUT	plhs[2]

#define vel_y_out(l,c) vel_y_out[(l)+vel_y_out_rows*(c)] // same


/* Matrix access */

#define k(l,c) k[(l)+k_rows*(c)] // k_rows = number of lines in k
#define pho_in(l,c) pho_in[(l)+pho_in_rows*(c)] // pho_in_rows = number of lines in pho_in
#define nu(l,c) nu[(l)+nu_rows*(c)] // nu_rows = number of lines in nu
#define vel_x_in(l,c) vel_x_in[(l)+vel_x_in_rows*(c)] // vel_x_in_rows = number of lines in vel_x_in
#define vel_y_in(l,c) vel_y_in[(l)+vel_y_in_rows*(c)] // vel_y_in_rows = number of lines in vel_y_in
#define dt(l,c) dt[(l)+dt_rows*(c)] // dt_rows = number of lines in dt
#define x(i, j) x[(i)+(j)*Ny]
#define d(i, j) d[(i)+(j)*Ny]
#define x0(i, j) x0[(i)+(j)*Ny]
#define d0(i, j) d0[(i)+(j)*Ny]
#define vx(i, j) vx[(i)+(j)*Ny]
#define vy(i, j) vy[(i)+(j)*Ny]
#define vx0(i, j) vx0[(i)+(j)*Ny]
#define vy0(i, j) vy0[(i)+(j)*Ny]

void copy(double *x_src, double *x_dst, int Nx, int Ny) {
    int i;
    
    for (i=0; i<Nx*Ny; i++) {
        x_dst[i] = x_src[i];
    }
}


void addSource(double *x, int Nx, int Ny, double srcPosX, double srcPosY) {
    x((int)(srcPosX*Nx),(int)(srcPosY*Ny)) = 1.0;
}


void setBoundary(double *x, int Nx, int Ny) {
    int i,j;
    
    j = 0;
    for (i=0; i<Nx; i++) {
        //x[i+j*Ny] = -x[i+(j+1)*Ny];
        //x[i+j*Ny] = x[i+(j+1)*Ny];
        x[i+j*Ny] = 0.0;
    }
    j = Ny-1;
    for (i=0; i<Nx; i++) {
        //x[i+j*Ny] = -x[i+(j-1)*Ny];
        //x[i+j*Ny] = x[i+(j-1)*Ny];
        x[i+j*Ny] = 0.0;
    }

    i = 0;
    for (j=0; j<Ny; j++) {
        //x[i+j*Ny] = -x[(i+1)+j*Ny];
        //x[i+j*Ny] = x[(i+1)+j*Ny];
        x[i+j*Ny] = 0.0;
    }
    i = Nx-1;
    for (j=0; j<Ny; j++) {
        //x[i+j*Ny] = -x[(i-1)+j*Ny];;
        //x[i+j*Ny] = x[(i-1)+j*Ny];;
        x[i+j*Ny] = 0.0;
    }
}

void difuse( double *x, double k, double dt, int Nx, int Ny, int maxIter, double *x0 ) {
    int i, j, iter;
    double a = dt*k*Nx*Ny;
    double si,so;
  
    /*
    si = 0.0;
    for (i=1; i<Nx-1; i++) {
        for (j=1; j<Ny-1; j++) {
            si += x0(i,j);
        }
    }
    si /= Nx*Ny;
    */
    
    for (iter=0; iter<maxIter; iter++) {
        for (i=1; i<Nx-1; i++) {
            for (j=1; j<Ny-1; j++) {
                x(i, j) = (x0(i,j) + a*( x(i-1,j) + x(i+1,j) + x(i,j-1) + x(i,j+1) ))/(1.0+4.0*a);
            }
        }
    }    
    
    /*
    so = 0.0;
    for (i=1; i<Nx-1; i++) {
        for (j=1; j<Ny-1; j++) {
            so += x(i,j);
        }
    }
    so /= Nx*Ny;
    
    for (i=1; i<Nx-1; i++) {
        for (j=1; j<Ny-1; j++) {
            x(i,j) *= si/so;
        }
    }    
    */
    
    setBoundary(x, Nx, Ny);
}




// void advert2( double *d, double *vx, double *vy, double dt, int Nx, int Ny, double *d0) {
//     int i, j, i0, j0, i1, j1;
//     double x, y, s0, t0, s1, t1;
//     
//     for (i=1; i<Nx-1; i++) {
//         for (j=1; j<Ny-1; j++) {
//             x = i-dt*vx(i,j);
//             y = j-dt*vy(i,j);
//             if (x<0.5) x=0.5; 
//             if (x>Nx+0.5) x=Nx+0.5;
//             if (y<0.5) y=0.5; 
//             if (y>Ny+0.5) y=Ny+0.5;
//             i0 = (int)x;
//             j0 = (int)y;
//             i1 = i0+1;
//             j1 = j0+1;
//             s1 = x-i0;
//             s0 = 1.0-s1;
//             t1 = y-j0;
//             t0 = 1.0-t1;
//             
//             d(i,j) = s0*( t0*d0(i0,j0) + t1*d0(i0,j1) ) + s1*( t0*d0(i1,j0) + t1*d0(i1,j1) );
//             
//         }
//     }
// 
//     setBoundary(d, Nx, Ny);
// }


void advect( double *d, double *vx, double *vy, double dt, int Nx, int Ny, double *d0) {
    int i, j;
    double dx, dy, dd_dx, dd_dy;
    double si, so;

    si = 0.0;
    for (i=1; i<Nx-1; i++) {
        for (j=1; j<Ny-1; j++) {
            si += d0(i,j);
        }
    }
    si /= Nx*Ny;
    
    for (i=1; i<Nx-1; i++) {
        for (j=1; j<Ny-1; j++) {
            dx = 1.0/(double)Nx;
            dy = 1.0/(double)Ny;
            
            if (j==1) {
                dd_dx = (d0(i,j+1)-d0(i,j-1))/(2.0*dx);
            } else {
                dd_dx = (-d0(i,j+2)+8*d0(i,j+1)-8*d0(i,j-1)+d0(i,j-2))/(12.0*dx);
            }
            if (i==1) {
                dd_dy = (d0(i+1,j)-d0(i-1,j))/(2.0*dy);
            } else {
                dd_dy = (-d0(i+2,j)+8*d0(i+1,j)-8*d0(i-1,j)+d0(i-2,j))/(12.0*dx);
            }
            
            d(i,j) = d0(i,j) - (vx(i,j)*dd_dx + vy(i,j)*dd_dy)*dt;            
        }
    }

    so = 0.0;
    for (i=1; i<Nx-1; i++) {
        for (j=1; j<Ny-1; j++) {
            so += d(i,j);
        }
    }
    so /= Nx*Ny;
    

    for (i=1; i<Nx-1; i++) {
        for (j=1; j<Ny-1; j++) {
            d(i,j) *= si/(so+1e-10);
        }
    }
    
    setBoundary(d, Nx, Ny);
}

















void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )
{
    int i,j;
    

    double *k;
    int *D_k, k_rows, k_cols;
    
    double *pho_in;
    int *D_pho_in, pho_in_rows, pho_in_cols;
    
    double *nu;
    int *D_nu, nu_rows, nu_cols;
    
    double *vel_x_in;
    int *D_vel_x_in, vel_x_in_rows, vel_x_in_cols;
    
    double *vel_y_in;
    int *D_vel_y_in, vel_y_in_rows, vel_y_in_cols;
    
    double *dt;
    int *D_dt, dt_rows, dt_cols;    

    double *srcPos;
    int *D_srcPos, srcPos_rows, srcPos_cols;    
        
    double *pho_out;
    int pho_out_rows, pho_out_cols;

    double *vel_x_out;
    int vel_x_out_rows, vel_x_out_cols;

    double *vel_y_out;
    int vel_y_out_rows, vel_y_out_cols;    
    
    
// Transforms Matlab variables mx_array into
// *double variables

    k = mxGetPr(k_IN);
    pho_in = mxGetPr(pho_in_IN);
    nu = mxGetPr(nu_IN);
    vel_x_in = mxGetPr(vel_x_in_IN);
    vel_y_in = mxGetPr(vel_y_in_IN);
    dt = mxGetPr(dt_IN);
    srcPos = mxGetPr(srcPos_IN);

// Size of the input matrices

    D_k = mxGetDimensions(k_IN);
    k_rows = D_k[0];
    k_cols = D_k[1];
    
    D_pho_in = mxGetDimensions(pho_in_IN);
    pho_in_rows = D_pho_in[0];
    pho_in_cols = D_pho_in[1];
    
    D_nu = mxGetDimensions(nu_IN);
    nu_rows = D_nu[0];
    nu_cols = D_nu[1];
    
    D_vel_x_in = mxGetDimensions(vel_x_in_IN);
    vel_x_in_rows = D_vel_x_in[0];
    vel_x_in_cols = D_vel_x_in[1];
    
    D_vel_y_in = mxGetDimensions(vel_y_in_IN);
    vel_y_in_rows = D_vel_y_in[0];
    vel_y_in_cols = D_vel_y_in[1];
    
    D_dt = mxGetDimensions(dt_IN);
    dt_rows = D_dt[0];
    dt_cols = D_dt[1];
    
    D_srcPos = mxGetDimensions(srcPos_IN);
    srcPos_rows = D_srcPos[0];
    srcPos_cols = D_srcPos[1];

    
    pho_out_rows = pho_in_rows;
    pho_out_cols = pho_in_rows;
    vel_x_out_rows = vel_x_in_rows;
    vel_x_out_cols = vel_x_in_cols;
    vel_y_out_rows = vel_y_in_rows;
    vel_y_out_cols = vel_y_in_cols;
        
// create output matrices (in this example, the same size inputs)

    pho_out_OUT = mxCreateDoubleMatrix(pho_out_rows,pho_out_cols,mxREAL);

    pho_out = mxGetPr(pho_out_OUT);
    vel_x_out_OUT = mxCreateDoubleMatrix(vel_x_out_rows,vel_x_out_cols,mxREAL);

    vel_x_out = mxGetPr(vel_x_out_OUT);
    vel_y_out_OUT = mxCreateDoubleMatrix(vel_y_out_rows,vel_y_out_cols,mxREAL);

    vel_y_out = mxGetPr(vel_y_out_OUT);

// here we go   
    
    // pho, vel_x and vel_y   MUST   be the same size
    int Ncols = pho_out_cols;
    int Nrows = pho_out_rows;

    for(i = 0; i<Nrows; i++) {
        for(j = 0; j<Ncols; j++) {
            pho_out(i,j) = pho_in(i,j);
        }
    }
    for(i = 0; i<Nrows; i++) {
        for(j = 0; j<Ncols; j++) {
            vel_x_out(i,j) = vel_x_in(i,j);
        }
    }
    for(i = 0; i<Nrows; i++) {
        for(j = 0; j<Ncols; j++) {
            vel_y_out(i,j) = vel_y_in(i,j);
        }
    }
        
    
    //addSource(pho_in, pho_in_rows, pho_in_cols, srcPos[0], srcPos[1]);
    difuse( pho_out, k[0], dt[0], Ncols, Nrows, NDIFUSE, pho_in);
    copy(pho_out, pho_in, Ncols, Nrows);
    advect(pho_out, vel_x_in, vel_y_in, dt[0], Ncols, Nrows, pho_in);
    
    difuse(vel_x_out, nu[0], dt[0], Ncols, Nrows, NDIFUSE, vel_x_in);
    difuse(vel_y_out, nu[0], dt[0], Ncols, Nrows, NDIFUSE, vel_y_in);
    copy(vel_x_out, vel_x_in, Ncols, Nrows);
    copy(vel_y_out, vel_y_in, Ncols, Nrows);
    advect(vel_x_out, vel_x_in, vel_y_in, dt[0], Ncols, Nrows, vel_x_in);
    advect(vel_y_out, vel_x_in, vel_y_in, dt[0], Ncols, Nrows, vel_y_in);
        
        
    
// here we went   
   
}

