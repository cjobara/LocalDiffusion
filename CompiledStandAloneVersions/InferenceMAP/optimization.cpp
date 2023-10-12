/*
 *
 *  InferenceMAP v1.0
 *  4/3/2015
 *
 *  Author: Mohamed El Beheiry, Physico-Chimie Curie, Institut Curie
 *  		Jean-Baptiste Masson, Physics of Biological Systems, Institut Pasteur
 *  Contact e-mail: mohamed.elbeheiry@gmail.com
 *  Copyright (c) 2015, Mohamed El Beheiry, Jean-Baptiste Masson, Institut Curie, Institut Pasteur
 *  All rights Reserved.
 *
 *  InferenceMAP is released under an "academic use only" licence.
 *  Details of which are provided in the InferenceMAP_License.doc file.
 *  Usage of InferenceMAP requires acceptance of this license.
 *
 *  User instructions for using InferenceMAP are provided in the InferenceMAP User Manual.
 *
 */

#include "stdafx.h"

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "globals.h"
#include "optimization.h"
#include "file.h"

extern Globals *iMAP;
extern File *file;

void dfpminOriginal(double p[],int ndim,  double gtol, int iter[], double fret[], double (func)(double [] ), void (dfunc)(double (func)(double []),double [], double []) ) {

	int check, i, ii, its, j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg, sumxi,temp,test;
	double* dg = new double[ndim];
	double* g = new double[ndim];
	double* hdg = new double[ndim];
	double* pnew = new double[ndim];
	double* xi = new double[ndim];

	double** hessin = new double*[ndim];
	for (int rr = 0; rr < ndim; rr++) {
		hessin[rr] = new double[ndim];
	}

	for (ii=0; ii < ndim; ii++) {
		dg[ii] = 0.0;
		g[ii] = 0.0;
		hdg[ii] = 0.0;
		pnew[ii] = 0.0;
		xi[ii] = 0.0;
	}

	fp = (func)(p) ;

	(dfunc)((func),p,g);

	for (i=0;i<ndim;i++){
		for (j=0; j<ndim; j++) {
			hessin[i][j]=0.0;
		}
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}

	stpmax = STPMX*FMAX(sqrt(sum), (double)ndim);

	for ( its=0; its<ITMAX; its++){
		iter[0]=its;

		fprintf(stderr,"Iteration %i\n",iter[0]+1);

		lnsrch(ndim,p,fp,g,xi,pnew,fret,stpmax,&check,(func));
		fp = fret[0];

		for (i=0; i<ndim; i++){
			xi[i] = pnew[i]-p[i];
			p[i] = pnew[i];
		}
		test =0.0;
		for (i=0 ; i<ndim ; i++){
			temp = fabs(xi[i])/FMAX(fabs(p[i]), 1.0);
			if (temp > test) {test = temp;}
		}
		if (test<TOLY){
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;
			for (int rr = 0; rr < ndim; rr++) {
				delete [] hessin[rr];
			}
			delete [] hessin;

			return;
		}

		for (i=0 ; i<ndim ; i++){ dg[i]=g[i];}

		(dfunc)((func),p,g);

		test = 0.0;
		den = FMAX(fret[0],1.0);
		for (i=0; i<ndim ;  i++){
			temp = fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) {test =temp;}
		}

		if (test < gtol){
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;
			for (int rr = 0; rr < ndim; rr++) {
				delete [] hessin[rr];
			}
			delete [] hessin;
			return;
		}

		for (i=0; i<ndim; i++) {dg[i] = g[i]-dg[i];}
		for (i=0; i<ndim; i++) {
			hdg[i] = 0.0;
			for (j=0; j<ndim; j++){ hdg[i] += hessin[i][j]*dg[j];}
		}
		fac =fae = sumdg = sumxi = 0.0;
		for (i=0; i<ndim ; i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}

		if (fac > sqrt(EPS*sumdg*sumxi)){
			fac = 1.0/fac;
			fad = 1.0/fae;
			for (i=0; i<ndim; i++ ){ dg[i]= fac*xi[i]-fad*hdg[i];}
			for (i=0;i<ndim ;i++){
				for (j=i ; j<ndim ; j++ ){
					hessin[i][j] += fac*xi[i]*xi[j]
					- fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
					hessin[j][i] = hessin[i][j];
				}
			}
		}
		for (i=0; i<ndim ; i++) {
			xi[i] =0.0;
			for (j=0; j<ndim; j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	
	delete [] dg;
	delete [] g;
	delete [] hdg;
	delete [] pnew;
	delete [] xi;
	for (int rr = 0; rr < ndim; rr++) {
		delete [] hessin[rr];
	}
	delete [] hessin;

	fprintf(stderr,"dfpmin fail\n"); return; exit(-1);
}

//void dfpmin(double p[],int ndim,  double gtol, int iter[], double fret[], double (func)(double [] ), void (dfunc)(double (func)(double []),double [], double []) ) {
void dfpmin(double *p,int ndim,  double gtol, int *iter, double *fret, double (func)(double *), void (dfunc)(double (func)(double *),double *, double *) ) {
	int check, i, ii, its, j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg, sumxi,temp,test;
	double* dg = new double[ndim];
	double* g = new double[ndim];
	double* hdg = new double[ndim];
	double* pnew = new double[ndim];
	double* xi = new double[ndim];

	double **hessin;
	char progress[100];

	// stack problems occur if n x n array declared inside function
	switch(file->meshType) {
	case 0: // square meshing
		hessin = file->squareMesh->hessian;
		break;
	case 1: // voroonoi meshing
		hessin = file->voronoiMesh->hessian;
		break;
	case 2: // quad-tree meshing
		hessin = file->treeMesh->hessian;
		break;
	default:
		// custom selection
		// single trajectory
		// batch trajectory
		hessin = new double*[3];
		for(int i = 0; i < 3; i++) {
			hessin[i] = new double[3];
		}
		break;
	}

	for (ii=0; ii < ndim; ii++) {
		dg[ii] = 0.0;
		g[ii] = 0.0;
		hdg[ii] = 0.0;
		pnew[ii] = 0.0;
		xi[ii] = 0.0;
	}

	fp = (func)(p) ;

	(dfunc)((func),p,g);

	for (i = 0; i < ndim; i++){
		if (iMAP->pauseCalculation == false) {
			for (j=0; j<ndim; j++) {
				hessin[i][j]=0.0;
			}
			hessin[i][i]=1.0;
		}
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	iMAP->pauseCalculation = false;

	stpmax = STPMX*FMAX(sqrt(sum), (double)ndim);

	for ( its=0; its<ITMAX; its++) {
		iter[0]=its;

		if ( (iMAP->updateDisplay && file->optimizationMode == 3)||
				(iMAP->updateDisplay && iMAP->updatePotentialCalculation)||
				(file->meshType == -1)||
				(file->optimizationMode == 4 && iMAP->updateDisplay)||
				(iMAP->updateDisplay && iMAP->smoothingPriorActive)) {

			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
			sprintf(progress,"Iteration %i:\t%.3f\n",iter[0]+1,fp);
			textDisplayUpdate(progress);
		}
		else { fprintf(stderr,"Iteration %i:\t%.3f\n",iter[0]+1,fp); }

		if (iMAP->pauseCalculation || iMAP->stopCalculation) { return; }

		lnsrch(ndim,p,fp,g,xi,pnew,fret,stpmax,&check,(func));

		fp = fret[0];

		for (i=0; i<ndim; i++){
			xi[i] = pnew[i]-p[i];
			p[i] = pnew[i];
		}
		test = 0.0;
		for (i=0 ; i<ndim ; i++){
			temp = fabs(xi[i])/FMAX(fabs(p[i]), 1.0);
			if (temp > test) {
				test = temp;
			}
		}

		if (test<TOLY){
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;

			return;
		}

		for (i=0 ; i<ndim ; i++){ dg[i]=g[i];}
		(dfunc)((func),p,g);
		test = 0.0;

		den = FMAX(fret[0],1.0);
		for (i=0; i<ndim ;  i++){
			temp = fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) {test =temp;}
		}

		if (test < gtol){ 
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;
			
			return;
		}

		for (i=0; i<ndim; i++) {dg[i] = g[i]-dg[i];}

		for (i=0; i<ndim; i++) {
			hdg[i] = 0.0;
			for (j=0; j<ndim; j++){ hdg[i] += hessin[i][j]*dg[j];}
		}

		fac =fae = sumdg = sumxi = 0.0;
		for (i=0; i<ndim ; i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}

		if (fac > sqrt(EPS*sumdg*sumxi)){
			fac = 1.0/fac;
			fad = 1.0/fae;
			for (i=0; i<ndim; i++ ){ dg[i]= fac*xi[i]-fad*hdg[i];}
			for (i=0; i<ndim ;i++){
				for (j=i ; j<ndim ; j++ ){
					hessin[i][j] += fac*xi[i]*xi[j]
					- fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
					hessin[j][i] = hessin[i][j];
				}
			}
		}
		for (i=0; i<ndim ; i++) {
			xi[i] =0.0;
			for (j=0; j<ndim; j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	delete [] dg;
	delete [] g;
	delete [] hdg;
	delete [] pnew;
	delete [] xi;

	fprintf(stderr,"dfpmin fail\n"); return; exit(-1);
}

void dfpminRO(double p[],int ndim,  double gtol, int iter[], double fret[], double (func)(double [] ), void (dfunc)(double (func)(double []),double [], double []), int itmax, int roIteration) {
	int check, i, ii, its, j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg, sumxi,temp,test;
	double* dg = new double[ndim];
	double* g = new double[ndim];
	double* hdg = new double[ndim];
	double* pnew = new double[ndim];
	double* xi = new double[ndim];

	double **hessin;
	char progress[60];

	// stack problems occur if n x n array declared inside function
	switch(file->meshType) {
	case 0: // regular meshing
		hessin = file->squareMesh->hessian;
		break;
	case 1: // irregular meshing
		hessin = file->voronoiMesh->hessian;
		break;
	case 2: // quad-tree meshing
		hessin = file->treeMesh->hessian;
		break;
	default:
		// custom selection
		// single trajectory
		// batch trajectory
		hessin = new double*[3];
		for(int i = 0; i < 3; i++) {
			hessin[i] = new double[3];
		}
		break;
	}

	for (ii=0; ii < ndim; ii++) {
		dg[ii] = 0.0;
		g[ii] = 0.0;
		hdg[ii] = 0.0;
		pnew[ii] = 0.0;
		xi[ii] = 0.0;
	}

	fp = (func)(p) ;

	(dfunc)((func),p,g);

	for (i=0;i<ndim;i++){
		for (j=0; j<ndim; j++) {
			hessin[i][j]=0.0;
		}
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}

	stpmax = STPMX*FMAX(sqrt(sum), (double)ndim);

	for ( its=0; its<itmax; its++){
		iter[0]=its;

		if (iMAP->updateDisplay) {
			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
			sprintf(progress,"Iteration %i : %i / %i\n",roIteration, iter[0]+1, itmax);
			textDisplayUpdate(progress);
		}
		else { fprintf(stderr,"Iteration %i\n",iter[0]+1); }

		if (iMAP->pauseCalculation || iMAP->stopCalculation) { return; }

		lnsrch(ndim,p,fp,g,xi,pnew,fret,stpmax,&check,(func));

		fp = fret[0];

		for (i=0; i<ndim; i++){
			xi[i] = pnew[i]-p[i];
			p[i] = pnew[i];
		}
		test =0.0;
		for (i=0 ; i<ndim ; i++){
			temp = fabs(xi[i])/FMAX(fabs(p[i]), 1.0);
			if (temp > test) {test = temp;}
		}
		if (test<TOLY){
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;
			return;
		}

		for (i=0 ; i<ndim ; i++){ dg[i]=g[i];}

		(dfunc)((func),p,g);

		test = 0.0;
		den = FMAX(fret[0],1.0);
		for (i=0; i<ndim ;  i++){
			temp = fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) {test =temp;}
		}

		if (test < gtol) {
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;
			return;
		}

		for (i=0; i<ndim; i++) {dg[i] = g[i]-dg[i];}

		for (i=0; i<ndim; i++) {
			hdg[i] = 0.0;
			for (j=0; j<ndim; j++){ hdg[i] += hessin[i][j]*dg[j];}
		}

		fac =fae = sumdg = sumxi = 0.0;
		for (i=0; i<ndim ; i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}

		if (fac > sqrt(EPS*sumdg*sumxi)){
			fac = 1.0/fac;
			fad = 1.0/fae;
			for (i=0; i<ndim; i++ ){ dg[i]= fac*xi[i]-fad*hdg[i];}
			for (i=0;i<ndim ;i++){
				for (j=i ; j<ndim ; j++ ){
					hessin[i][j] += fac*xi[i]*xi[j]
					- fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
					hessin[j][i] = hessin[i][j];
				}
			}
		}
		for (i=0; i<ndim ; i++) {
			xi[i] =0.0;
			for (j=0; j<ndim; j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	delete [] dg;
	delete [] g;
	delete [] hdg;
	delete [] pnew;
	delete [] xi;
}

void lnsrch(int ndim,double xold[], double fold, double g[] , double p[], double xx[] , double f[], double stpmax, int check[],  double (func)(double [])  ){

	int i;
	double a,alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp, test, tmplam;

	check[0]=0;
	for (sum=0.0, i=0; i<ndim; i++){ sum += p[i]*p[i];}
	sum = sqrt(sum);

	if (sum > stpmax)
		for (i=0; i<ndim;i++) p[i] *= stpmax/sum;
	for (slope=0.0, i=0;i<ndim;i++ ) {
		slope += g[i]*p[i];
	}
	
	if (slope >= 0.0) {

		iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
		textDisplayUpdate("Line Search Failure\n");

		return;
		exit(1);
	}
	test = 0.0;
	for (i=0; i<ndim; i++){
		temp = fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test = temp;
	}
	alamin = TOLX/test;
	alam = 1.0;
	for (;;) {
		for (i=0;i<ndim;i++) xx[i]=xold[i]+alam*p[i];

		f[0] =   (func)( xx );


		if (alam < alamin) {
			for (i=0; i<ndim ; i++){ xx[i] = xold[i];}
			check[0]=1;
			return ;
		} else if (f[0] <= fold +ALF*alam*slope) {
			return;
		}
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f[0] - fold -slope));
			else {
				rhs1 = f[0] -fold -alam*slope;
				rhs2 = f2 - fold -alam2*slope;
				a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2);
				b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0){ tmplam = -slope/(2.0*b);}
				else {
					disc = b*b-3.0*a*slope;
					if (disc < 0.0){ tmplam = 0.5*alam;}
					else if (b <=0.0){ tmplam = (-b +sqrt(disc))/(3.0*a);}
					else {tmplam = -slope/(b+sqrt(disc));}
				}

				if (tmplam > 0.5 *alam)
					tmplam = 0.5 *alam;
			}
		}
		alam2 = alam;
		f2 =f[0];
		alam = FMAX(tmplam, 0.1*alam);
	}

}

void dfunc(double (func)(double *), double *yy, double *ans){

	double h=1.e-8, ERR;
	int i,j, indice, indice2;
	double errt,fac,hh;

	int ndim;

	switch(file->meshType) {
		case 0: // square mesh
			ndim = file->squareMesh->getDimensions();
//			printf("ndim = %i\n",ndim);
			break;
		case 1: // voronoi mesh
			ndim = file->voronoiMesh->getDimensions();
			break;
		case 2: // quad-tree mesh
			ndim = file->treeMesh->getDimensions();
			break;
		default:
			ndim = 3;
			break;
	}

	double* veclocal1 = new double[ndim];
	double* veclocal2 = new double[ndim];

	double a[NTAB][NTAB];

	for (int aa = 0; aa < NTAB; aa++) {
		for (int bb = 0; bb < NTAB; bb++) {
			a[aa][bb] = 0.0;
		}
	}

	for (indice=0; indice<ndim ; indice++){

		hh=h;
		for (indice2=0; indice2<ndim ; indice2++ ){
			if (indice2 == indice){
				veclocal1[indice2]= yy[indice2]+hh;
				veclocal2[indice2]= yy[indice2]-hh;
			} else {
				veclocal1[indice2] = yy[indice2];
				veclocal2[indice2] = yy[indice2];
			}
		}

		a[0][0]=( (func)(veclocal1 ) - (func)(veclocal2)  )/(2.0*hh);

		ERR=BIG;

		for(i=1 ; i<NTAB ;  i++ ){
			hh /= CON;

			for (indice2=0; indice2<ndim ; indice2++ ){
				if (indice2 == indice){
					veclocal1[indice2]= yy[indice2]+hh;
					veclocal2[indice2]= yy[indice2]-hh;
				} else {
					veclocal1[indice2] = yy[indice2];
					veclocal2[indice2] = yy[indice2];
				}
			}

			a[0][i]=( (func)(veclocal1) - (func)(veclocal2)  )/(2.0*hh);

			fac = CON2;
			for (j=1 ; j<i; j++){
				a[j][i]= (a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
				fac = CON2*fac;
				errt=FMAX( fabs(a[j][i]-a[j-1][i]) , fabs(a[j][i]-a[j-1][i-1])  );
				if (errt <= ERR){
					ERR = errt;
					ans[indice] = a[j][i];
				}
			}
			if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(ERR)) break;
		}

	}

	delete [] veclocal1;
	delete [] veclocal2;

	return ;
}

void vander(double *xxx, double *cof, double *yyy, int nn){
	// resout x^{k}*w=V input x et V et cela donne w
	int k,j,i;
	double phi,ff,b;

	double* s = new double[nn+1];
	for (i=0;i<=nn;i++) s[i]=cof[i]=0.0;
	s[nn] = -xxx[0];
	for (i=1;i<=nn;i++) {
		for (j=nn-i;j<=nn-1;j++)
			s[j] -= xxx[i]*s[j+1];
		s[nn] -= xxx[i];
	}
	for (j=0;j<=nn;j++) {
		phi=nn+1;
		for (k=nn;k>=1;k--)
			phi=k*s[k]+xxx[j]*phi;
		ff=yyy[j]/phi;
		b=1.0;
		for (k=nn;k>=0;k--) {
			cof[k] += b*ff;
			b=s[k]+xxx[j]*b;
		}
	}
	delete [] s;
}

double amotry(double **p, int ndim, double y[], double psum[], double funk(double []), int ihi, double fac){

	int j,jj;
	double fac1, fac2, ytry;
	double *ptry;
	// double ptry[2];

	ptry  =  (double *)malloc((size_t) ndim*sizeof(double));

	fac1 = (1.0-fac)/ndim;
	fac2 = fac1 - fac;
	for (j = 0; j < ndim; j++) { ptry[j] = psum[j]*fac1-p[ihi][j]*fac2; }
	ytry = funk(ptry);
	if (ytry < y[ihi]) {
	y[ihi] = ytry;
		for (j=0; j<ndim; j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}

	for (jj=0;jj<ndim; jj++) {
	ptry[jj]=0.0;
	}
	return ytry;

}

int amoeba(double **p,int ndim,  double y[], double ftol, double funk(double []), int *nfunk, int maxi) {
	double amotry(double **p, int ndim, double y[], double psum[] , double funk(double []), int ihi, double fac);

	int i,ii,ihi,ilo,inhi,j,mpts=ndim+1;

	double rtol,sum,swap,ysave,ytry;
	double *psum;
	psum =    (double *)malloc((size_t) ndim*sizeof(double));
	*nfunk = 0;
	GET_PSUM

	for (j=0; j< ndim ; j++) {
		for (sum=0.0, i=0; i < mpts; i++ ) {
			sum += p[i][j];
		}
		psum[j]=sum;
	}

	for (;;) {
		ilo = 1;
		ihi = y[0] > y[1] ? (inhi=1,0) : (inhi=0,1);

		for (i=0 ; i< mpts ; i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi] ) {
				inhi = ihi;
				ihi = i;
			} else if (y[i] > y[inhi] && i != ihi) { inhi = i; }
		}

		rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
		if (rtol < ftol) {
			SWAP(y[0],y[ilo]);
			for (i = 0; i < ndim ; i++) { SWAP(p[0][i], p[ilo][i]); }
			fprintf(stderr,"Optimized\n") ;
			return 1;
		}

		if (*nfunk >= maxi) {
			fprintf(stderr,"Maximum Iterations\n") ;
			return 0;
		}

		*nfunk += 2;

		ytry = amotry(p, ndim, y, psum, funk, ihi, -1.0);

		if (ytry <= y[ilo]) {
			ytry = amotry(p, ndim, y, psum, funk, ihi, 2.0);
		}
		else if (ytry >= y[inhi]) {
			ysave = y[ihi];
			ytry = amotry(p, ndim, y, psum, funk, ihi, 0.5);

			if (ytry >= ysave ) {
				for (i=0; i < mpts; i++) {
					if (i != ilo) {
						for (j=0; j<ndim ; j++) {
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						}
						y[i]= funk(psum);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}

		} else --(*nfunk);
	}
	for (ii=0; ii<ndim; i++ ) {
		psum[ii]=0.0;
	}
	return 0;
}

double *nrvector(int nl, int nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) fprintf(stderr,"allocation failure in vector()\n");
	return v-nl+NR_END;
}

void free_vector(double *v, int nl, int nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double ran1(long *idum){
	/*
	 * "Minimal" random number generator of Park and Miller with Bays-Durham shu e and added
	 * safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
	 * values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
	 * successive deviates in a sequence. RNMX should approximate the largest floating value that is
	 * less than 1.
	 */

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB_bis];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB_bis+7;j>=0;j--) {
			k=(*idum)/IQ_bis;
			*idum=IA_bis*(*idum-k*IQ_bis)-IR_bis*k;
			if (*idum < 0) *idum += IM_bis;
			if (j < NTAB_bis) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ_bis;
	*idum=IA_bis*(*idum-k*IQ_bis)-IR_bis*k;
	if  (*idum < 0)  *idum += IM_bis;
	j=iy/NDIV_bis;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM_bis*iy) > RNMX_bis) return RNMX_bis;
	else return temp;
}

double gasdev(long *idum) {
	double ran1(long *idum);
	static int iset=0;
	static float gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

