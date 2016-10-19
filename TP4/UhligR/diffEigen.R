#
#   diffEigen package
#   Copyright (C) 2007  Jan de Leeuw <deleeuw@stat.ucla.edu>
#   UCLA Department of Statistics, Box 951554, Los Angeles, CA 90095-1554
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################
#
# version 0.0.1, 2007-12-15   Initial Alpha Release
# version 0.0.2, 2007-12-16   Added scaled versions of vectors
#


#  gevd compute the generalized eigenvalue 
# decomposition for (a,b)

gevd<-function(a,b=diag(nrow(a))) {
	bs<-mfunc(b,function(x) ginvx(sqrt(x)))
	ev<-eigen(bs%*%a%*%bs)
	return(list(gvalues=ev$values,gvectors=bs%*%ev$vectors))
}

# gsvd computes the generalized singular value 
# decomposition for (r,p,q)

gsvd<-function(r,p=diag(nrow(a)),q=diag(ncol(a))) {
	ps<-mfunc(p,function(x) ginvx(sqrt(x)))
	qs<-mfunc(q,function(x) ginvx(sqrt(x)))
	sv<-svd(ps%*%r%*%qs)
	return(list(gd=sv$d,gu=ps%*%sv$u,gv=qs%*%sv$v))
}

# ginvgevd: the (1,2) inverse needed for derivatives of
# generalized eigenvectors

ginvgevd<-function(ge,ind) {
	y<-ge$gvectors; v<-ge$gvalues-(ge$gvalues[ind])
	return(tcrossprod(y%*%diag(ginvx(v)),y))
}

# ginvgsvd: the (1,2) inverse needed for derivatives of
# generalized singular vectors

ginvgsvd<-function(gs,p,ind) {
	x<-gs$gu; z<-gs$gv; gm<-gs$gd; gi<-gm[ind]
	iv1<-(ginvx(gm-gi)+ginvx(-(gm+gi)))/2
	iv2<-(ginvx(gm-gi)-ginvx(-(gm+gi)))/2
	n<-nrow(x); m<-nrow(z)
	a<-matrix(0,n+m,n+m)
	b<-tcrossprod(z%*%diag(iv2),x)
	a[1:n,n+(1:m)]<-t(b); a[n+(1:m),1:n]<-b
	a[n+(1:m),n+(1:m)]<-tcrossprod(z%*%diag(iv1),z)
	a[1:n,1:n]<-tcrossprod(x%*%diag(iv1),x)
	a[1:n,1:n]<-a[1:n,1:n]-(solve(p)-tcrossprod(x))/gi
	return(a)
}

# gevdDer: generalized eigenvalue decomposition plus derivatives
# needs four functions to compute a, b, da, and db

gevdDer<-function(par,aPar,bPar,daPar,dbPar,ind=1) {
	a<-aPar(par); b<-bPar(par); da<-daPar(par); db<-dbPar(par)
	nord<-nrow(a); neval<-length(ind); npars<-length(par)
	ge<-gevd(a,b); gv<-ge$gvectors; gd<-ge$gvalues
	dl<-matrix(0,neval,npars); dy<-array(0,c(nord,neval,npars))
	for (i in 1:neval) {
		j<-ind[i]; y<-gv[,j]; lb<-gd[j]
		dl[i,]<-apply(da-lb*db,3,function(x) y%*%x%*%y)
		aux0<-ginvgevd(ge,j)
		aux1<-apply(da-lb*db,3,function(x) aux0%*%x%*%y)
		aux2<-outer(y,apply(db,3,function(x) y%*%x%*%y))/2
		dy[,i,]<--(aux1+aux2)
		}
	return(list(gd=gd,gv=gv,dl=dl,dy=dy,ind=ind))
}

# gsvdDer: generalized singular value decomposition plus derivatives
# needs six functions to compute p, q, r, dp, dq, and dr


gsvdDer<-function(par,pPar,qPar,rPar,dpPar,dqPar,drPar,ind=1) {
	p<-pPar(par); q<-qPar(par); r<-rPar(par)
	dp<-dpPar(par); dq<-dqPar(par); dr<-drPar(par)
	nrows<-nrow(p); ncols<-nrow(q); neval<-length(ind); npars<-length(par)
	gs<-gsvd(r,p,q); gu<-gs$gu; gv<-gs$gv; gd<-gs$gd
	dl<-matrix(0,neval,npars)
	dx<-array(0,c(nrows,neval,npars)); dz<-array(0,c(ncols,neval,npars))
	for (i in 1:neval) {
		j<-ind[i]; x<-gu[,j]; z<-gv[,j]; gi<-gd[j]
		xz<-rbind(cbind(x),cbind(z))
		dpx<-apply(dp,3,function(d) x%*%d%*%x) 
		dqz<-apply(dq,3,function(d) z%*%d%*%z)
		dl[i,]<-apply(dr,3,function(d) x%*%d%*%z)-gi*(dpx+dqz)/2
		aux0<-ginvgsvd(gs,p,j)
		kv<-array(0,c(nrows+ncols,nrows+ncols,npars))
		kv[1:nrows,1:nrows,]<--gi*dp
		kv[1:nrows,nrows+(1:ncols),]<-dr
		kv[nrows+(1:ncols),1:nrows,]<-aperm(dr,c(2,1,3))
		kv[nrows+(1:ncols),nrows+(1:ncols),]<--gi*dq 
		aux1<-apply(kv,3,function(d) aux0%*%d%*%xz)
		aux2<-drop(outer(xz,(dpx+dqz)/4))
		dxz<--(aux1+aux2)
		dx[,i,]<-dxz[1:nrows,]
		dz[,i,]<-dxz[nrows+(1:ncols),]
		}
	return(list(gd=gd,gu=gu,gv=gv,dl=dl,dx=dx,dz=dz,ind=ind))	
}

# gevdScal: generalized eigen value decomposition with eigen vector
# length scaled to generalized eigenvalue

gevdScal<-function(ge) {
	gd<-ge$gd; gv<-ge$gv; dl<-ge$dl; dy<-ge$dy
	dys<-dy; gvs<-gv; ind<-ge$ind; neval<-length(ind)
	for (i in 1:ncol(gv))
		gvs[,i]<-gv[,i]*sqrt(gd[i])
	for (i in 1:neval) {
		j<-ind[i]
		dys[,i,]<-1/(2*sqrt(gd[j]))*outer(gv[,j],dl[i,])+sqrt(gd[j])*dy[,i,]
		}
return(list(gd=gd,gv=gvs,dl=dl,dy=dys,ind=ind))	
}

# gsvdScal: generalized singular value decomposition with singular vector
# length scaled to generalized singular value (in four ways)

gsvdScal<-function(gs,scal="be") {
	gd<-gs$gd; gu<-gs$gu; gv<-gs$gv
	dl<-gs$dl; dx<-gs$dx; dz<-gs$dz
	dxs<-dx; dzs<-dz; gus<-gu; gvs<-gv
	ind<-gs$ind; neval<-length(ind)
	if (scal=="be") {
		for (i in length(gd)) {
			gus[,i]<-gu[,i]*gd[i]
			gvs[,i]<-gv[,i]*gd[i]
			}
		for (i in 1:neval) {
			j<-ind[i]
			dxs[,i,]<-outer(gu[,j],dl[i,])+gd[j]*dx[,i,]
			dzs[,i,]<-outer(gv[,j],dl[i,])+gd[j]*dz[,i,]
			}
		}	
	if (scal=="go") {
		for (i in length(gd)) {
			gus[,i]<-gu[,i]*sqrt(gd[i])
			gvs[,i]<-gv[,i]*sqrt(gd[i])
			}
		for (i in 1:neval) {
			j<-ind[i]
			dxs[,i,]<-1/(2*sqrt(gd[j]))*outer(gu[,j],dl[i,])+sqrt(gd[j])*dx[,i,]
			dzs[,i,]<-1/(2*sqrt(gd[j]))*outer(gv[,j],dl[i,])+sqrt(gd[j])*dz[,i,]
			}
		}	
	if (scal=="rc") {
		for (i in length(gd)) {
			gvs[,i]<-gv[,i]*gd[i]
			}
		for (i in 1:neval) {
			j<-ind[i]
			dzs[,i,]<-outer(gv[,j],dl[i,])+gd[j]*dz[,i,]
			}
		}	
	if (scal=="cr") {
		for (i in length(gd)) {
			gus[,i]<-gu[,i]*gd[i]
			}
		for (i in 1:neval) {
			j<-ind[i]
			dxs[,i,]<-outer(gu[,j],dl[i,])+gd[j]*dx[,i,]
			}
		}	
return(list(gd=gd,gu=gus,gv=gvs,dl=dl,dx=dxs,dz=dzs,ind=ind))	
}

# ginvx is a helper to compute reciprocals

ginvx<-function(x) {ifelse(x==0,0,1/x)}

# mfunc is a helper to compute matrix functions

mfunc<-function(a,fn=sqrt) {
	e<-eigen(a); y<-e$vectors; v<-e$values
	return(tcrossprod(y%*%diag(fn(v)),y))
}

