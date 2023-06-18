C***********************************************************************
	subroutine ctrsc(x,n,p,xbar,sdv,mvcode)
C Centers and scales the data matrix so that the observed data in every
C column have mean zero and variance one. If a column has zero variance
C or less than 2 observations then the data are centered (set equal 
C to zero)
	integer n,p,count
C real x(n,p),mvcode was replaced by line below
        double precision x(n,p),mvcode
	double precision xbar(p),sdv(p),sum1,sum2
	do 10 j=1,p
           sum1=0
           sum2=0
           count=0
           do 5 i=1,n
              if(x(i,j).ne.mvcode) then
                 count=count+1
                 sum1=sum1+x(i,j)
                 sum2=sum2+x(i,j)**2.
              endif
5	   continue
           if(count.gt.0) then
              xbar(j)=sum1/count
              sdv(j)=sqrt((sum2-(sum1**2.)/count)/count)
              do 7 i=1,n
                 if(x(i,j).ne.mvcode) x(i,j)=x(i,j)-xbar(j)
7             continue
              if(sdv(j).gt.0.) then
                 do 9 i=1,n
                    if(x(i,j).ne.mvcode) x(i,j)=x(i,j)/sdv(j)
9                continue
               else 
                 sdv(j)=1.
              endif
           else
              sdv(j)=1.
           endif
10	continue
	return
	end
C***********************************************************************
        subroutine mkpsi(p,psi)
C Generates a symmetric matrix of integers indicating the linear
C position in packed storage of the matrix elements
        integer p,psi(0:p,0:p),posn
        posn=0
        do 10 j=0,p
           posn=posn+1
           psi(j,j)=posn
           do 5 k=j+1,p
              posn=posn+1
              psi(j,k)=posn
              psi(k,j)=posn
5          continue
10      continue
        return
        end
C***********************************************************************
        subroutine swp(d,theta,pivot,p,psi,submat,dir)
C Performs sweep on a symmetric matrix in packed storage.
C Sweeps on pivot position. Sweeps only the (0:submat,0:submat)
C submatrix.
C If dir=1, performs ordinary sweep. If dir=-1, performs reverse sweep.
        integer d,p,pivot,psi(0:p,0:p),submat,dir
        double precision theta(d),a,b,c
        a=theta(psi(pivot,pivot))
        theta(psi(pivot,pivot))=-1./a
        do 10 j=0,submat
           if(j.ne.pivot) theta(psi(j,pivot))=theta(psi(j,pivot))/a*dir
10      continue
        do 30 i=0,submat
           do 20 j=i,submat
              if((i.ne.pivot).and.(j.ne.pivot))then
                 b=theta(psi(i,pivot))
                 c=theta(psi(j,pivot))
                 theta(psi(i,j))=theta(psi(i,j))-a*b*c
              endif
20         continue
30      continue
        return
        end
C***********************************************************************
        subroutine initn(d,theta)
C Initializes theta
        integer d
        double precision theta(d)
        theta(1)=1.
        do 1 i=2,d
           theta(i)=0.
1       continue
        return
        end
C***********************************************************************
        subroutine stvaln(d,theta,p,psi)
C Gets starting value of theta: mu=0 and sigma=I
        integer d,p,psi(0:p,0:p)
        double precision theta(d)
        call initn(d,theta)
        theta(1)=-1.
        do 5 j=1,p
           theta(psi(j,j))=1.
5       continue
        return
        end
C***********************************************************************
        subroutine tobsn(d,tobs,p,psi,n,x,npatt,r,mdpst,nmdp,oc)
C Tabulates the known part of the sscp matrix for all missingness
C patterns.
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,patt
C real x(n,p) was changed to the line below
        double precision x(n,p)
        double precision tobs(d)
        call initn(d,tobs)
        do 40 patt=1,npatt
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 30 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 20 j=1,noc
                 tobs(psi(0,oc(j)))=tobs(psi(0,oc(j)))+x(i,oc(j))
                 do 10 k=j,noc
                    tobs(psi(oc(j),oc(k)))=tobs(psi(oc(j),oc(k)))+
     /                   x(i,oc(j))*x(i,oc(k))
10               continue
20            continue
30         continue
40      continue
        return
        end
C************************************************************************
        subroutine emn(d,theta,t,tobs,p,psi,n,x,npatt,r,mdpst,nmdp,
     /      oc,mc,c,mle,tau,m,mu,lmbinv)
C Performs one step of em. Theta must be in sweep(0) condition. 
C After execution, the new parameter value is contained in t, and
C theta is left swept on columns corresponding to observed variables
C in the first missingness pattern. 
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,mc(p),nmc,patt,mle
C real x(n,p) was changed to the line below
        double precision x(n,p)
        double precision theta(d),t(d),tobs(d),c(p),tau,m,mu(p)
        double precision lmbinv(p,p)
        do 1 i=1,d
           t(i)=tobs(i)
1       continue
        do 200 patt=1,npatt
           call swpobs(d,theta,p,psi,npatt,r,patt)
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 150 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 50 j=1,nmc
                 c(mc(j))=theta(psi(0,mc(j)))
                 do 40 k=1,noc
                    c(mc(j))=c(mc(j))+theta(psi(oc(k),mc(j)))*x(i,oc(k))
40               continue
50            continue
              do 100 j=1,nmc
                 t(psi(0,mc(j)))=t(psi(0,mc(j)))+c(mc(j))
                 do 70 k=1,noc
                    t(psi(oc(k),mc(j)))=t(psi(oc(k),mc(j)))+
     /                   x(i,oc(k))*c(mc(j))
70               continue
                 do 80 k=j,nmc
                    t(psi(mc(k),mc(j)))=t(psi(mc(k),mc(j)))+
     /                c(mc(k))*c(mc(j))+theta(psi(mc(k),mc(j)))
80               continue
100           continue
150        continue
200     continue
        if(mle.eq.0) then
           call moden(d,t,p,psi,n,tau,m,mu,lmbinv)
        endif
        do 210 i=2,d
           t(i)=t(i)/dble(n)
210     continue
        call swp(d,t,0,p,psi,p,1)
        return
        end
C************************************************************************
        subroutine moden(d,t,p,psi,n,tau,m,mu,lmbinv)
C Alters the sufficient statistics to yield a posterior mode.
        integer d,p,psi(0:p,0:p),n
        double precision t(d),tau,m,mu(p),lmbinv(p,p),c,e
        do 5 j=1,p
           mu(j)=mu(j)*dble(n)
5       continue
        c=tau/(dble(n)*(tau+dble(n)))
        e=dble(n)/(dble(n)+m+dble(p)+2.)
        do 20 j=1,p
           do 10 k=j,p
              t(psi(j,k))=t(psi(j,k))+lmbinv(j,k)-(t(psi(0,j))*
     /           t(psi(0,k)))/dble(n)
              t(psi(j,k))=t(psi(j,k))+c*(t(psi(0,j))-mu(j))*
     /           (t(psi(0,k))-mu(k))
              t(psi(j,k))=t(psi(j,k))*e
10         continue
20      continue
        c=dble(n)/(tau+dble(n))
        e=1.-c
        do 30 j=1,p
           t(psi(0,j))=c*t(psi(0,j))+e*mu(j)
30      continue
        do 50 j=1,p
           do 40 k=j,p
              t(psi(j,k))=t(psi(j,k))+(t(psi(0,j))*
     /           t(psi(0,k)))/dble(n)
40         continue
50      continue
        return
        end
C************************************************************************
        subroutine lprin(d,theta,p,psi,c,tau,m,mu,lmbinv,logpri)
C Evaluates log-prior at theta
        integer d,p,psi(0:p,0:p)
        double precision theta(d),c(p),logpri,logdet,trace
        double precision tau,m,mu(p),lmbinv(p,p)
        logdet=dble(0)  
        do 5 j=1,p
           c(j)=theta(psi(0,j))-mu(j)
5       continue
        do 10 j=1,p
           logdet=logdet+log(theta(psi(j,j)))
           call swp(d,theta,j,p,psi,p,1)
10      continue
        logpri=-logdet*(m+dble(p)+2.)/2.
        trace=dble(0)
        do 120 j=1,p
           do 110 k=1,p
              trace=trace-theta(psi(j,k))*(lmbinv(j,k)+tau*c(j)*c(k))
110        continue
120     continue
        logpri=logpri-trace/2.
        return
        end
C************************************************************************
        subroutine lobsn(d,theta,t,p,psi,n,x,npatt,r,mdpst,nmdp,oc,
     /    c,loglik)
C Evaluates observed-data loglikelihood at theta
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,patt
C real x(n,p) was changed to the line below
        double precision x(n,p)
        double precision theta(d),t(d),c(p),loglik,logdet,trace
        loglik=dble(0)
        logdet=dble(0)  
        do 5 j=1,p
           c(j)=theta(psi(0,j))
5       continue
        do 200 patt=1,npatt
           call initn(d,t)
           do 10 j=1,p
              if((r(patt,j).eq.1).and.(theta(psi(j,j)).gt.0.))then
                 logdet=logdet+log(theta(psi(j,j)))
                 call swp(d,theta,j,p,psi,p,1)
              elseif((r(patt,j).eq.0).and.(theta(psi(j,j)).lt.0.))then
                 call swp(d,theta,j,p,psi,p,-1)
                 logdet=logdet-log(theta(psi(j,j)))
              endif
10         continue
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 100 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 30 j=1,noc
                 t(psi(0,oc(j)))=x(i,oc(j))-c(oc(j))
30            continue
              do 50 j=1,noc
                 do 40 k=j,noc
                    t(psi(oc(j),oc(k)))=t(psi(oc(j),oc(k)))+
     /                 t(psi(0,oc(j)))*t(psi(0,oc(k)))
40               continue
50            continue
100        continue
           trace=dble(0)
           do 120 j=1,noc
              do 110 k=1,noc
                 trace=trace-(theta(psi(oc(j),oc(k)))*
     /              t(psi(oc(j),oc(k))))
110           continue
120        continue
           loglik=loglik-(dble(nmdp(patt))*logdet/2.)-(trace/2.)
200     continue
        return
        end
C************************************************************************
        subroutine swpobs(d,theta,p,psi,npatt,r,patt)
C Sweeps theta to condition on the observed variables
        integer d,p,psi(0:p,0:p),npatt,r(npatt,p),patt
        double precision theta(d)
        do 10 j=1,p
           if((r(patt,j).eq.1).and.(theta(psi(j,j)).gt.0.))then
              call swp(d,theta,j,p,psi,p,1)
           elseif((r(patt,j).eq.0).and.(theta(psi(j,j)).lt.0.))then
              call swp(d,theta,j,p,psi,p,-1)
           endif
10      continue
        return
        end
C************************************************************************
        subroutine gtmc(p,npatt,r,patt,mc,nmc,last)
C Finds the column numbers of the missing variables, and stores them
C in the first nmc elements of mc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,mc(p),nmc,last
        nmc=0
        do 10 j=1,last
           if(r(patt,j).eq.0)then
              nmc=nmc+1
              mc(nmc)=j
           endif
10      continue
        return
        end
C************************************************************************
        subroutine gtoc(p,npatt,r,patt,oc,noc,last)
C Finds the column numbers of the observed variables, and stores them
C in the first noc elements of oc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,oc(p),noc,last
        noc=0
        do 10 j=1,last
           if(r(patt,j).eq.1)then
              noc=noc+1
              oc(noc)=j
           endif
10      continue
        return
        end
C************************************************************************
        function rangen(init)
        integer a,p,ix,b15,b16,xhi,xalo,leftflo,fhi,k,init
        data a/16807/,b15/32768/,b16/65536/,p/2147483647/
        save ix
        if(init.ne.0) ix=init
        xhi=ix/b16
        xalo=(ix-xhi*b16)*a
        leftflo=xalo/b16
        fhi=xhi*a+leftflo
        k=fhi/b15
        ix=(((xalo-leftflo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix.lt.0)ix=ix+p
        rangen=float(ix)*4.656612875E-10
        return
        end
C***********************************************************************
        subroutine rngs(seed)
C initializes rangen with seed
        integer seed
        tmp=rangen(seed)
        return
        end
C***********************************************************************
        function chisq(df)
C Generates a chisquare deviate with df degrees of freedom
C real df was changed to the line below
	double precision df
        chisq=2.*gamm(df/2.)
        return
        end
C***********************************************************************
        function gamm(a)
C Generates a random gamma(a) deviate. If a>=1, uses the method of 
C Fishman (1976); if 0<a<1, the method of Ahrens (1974)
C real a,u,y,q,e,b,p,u1
        double precision a,u,y,q,e,b,p,u1
        data e/2.718282/
        if(a.ge.1)then
1          continue
           u=rangen(0)
           y=-log(rangen(0))
           q=(y/exp(y-1))**(a-1)
           if(u.le.q)then
              gamm=a*y
           else
              goto 1
           endif
        else
2          continue
           u=rangen(0)
           b=(e+a)/e
           p=b*u
           if(p.gt.1) goto 4
3          continue
           x=p**(1/a)
           u1=rangen(0)
           if(u1.gt.(e**(-x)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
4          continue
           x=-log((b-p)/a)
           u1=rangen(0)
           if(u1.gt.(x**(a-1)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
        endif
10      continue
        return
        end
C***********************************************************************
        subroutine is1n(d,theta,t,tobs,p,psi,n,x,npatt,r,mdpst,nmdp,
     /    oc,mc,z,c)
C Performs I-step of data augmentation. Randomly draws missing data Xmis
C given theta, and stores the sufficient statistics in t.
C Theta must be in sweep(0) condition. Answer is returned in
C unswept condition.
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,mc(p),nmc,patt
C real x(n,p),z(p),junk
	double precision x(n,p),z(p),junk
        double precision theta(d),t(d),tobs(d),c(d)
        junk=gauss()
        do 1 i=1,d
           t(i)=tobs(i)
1       continue
        do 200 patt=npatt,1,-1
           call swpobs(d,theta,p,psi,npatt,r,patt)
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           call gtoc(p,npatt,r,patt,oc,noc,p)
           call sigex(d,theta,c,p,psi,mc,nmc)
           call chols(d,c,p,psi,mc,nmc)
           do 150 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 100 j=1,nmc
                 x(i,mc(j))=theta(psi(0,mc(j)))
                 do 40 k=1,noc
                    x(i,mc(j))=x(i,mc(j))+theta(psi(oc(k),
     /                 mc(j)))*x(i,oc(k))
40               continue
                 z(mc(j))=gauss()
                 do 50 k=1,j
                    x(i,mc(j))=x(i,mc(j))+z(mc(k))*c(psi(mc(j),mc(k)))
50               continue
                 t(psi(0,mc(j)))=t(psi(0,mc(j)))+x(i,mc(j))
                 do 60 k=1,noc
                    t(psi(mc(j),oc(k)))=t(psi(mc(j),oc(k)))+
     /                   x(i,mc(j))*x(i,oc(k))
60               continue
                 do 70 k=1,j
                    t(psi(mc(j),mc(k)))=t(psi(mc(j),mc(k)))+
     /                   x(i,mc(j))*x(i,mc(k))
70               continue
100           continue
150        continue
200     continue
C Divide t by n
        do 210 i=2,d
           t(i)=t(i)/dble(n)
210     continue
        return
        end
C***********************************************************************
        subroutine sigex(d,theta,extr,p,psi,mc,nmc)
C Extracts submatrix of theta corresponding to the columns of mc
        integer d,p,psi(0:p,0:p),mc(p),nmc
        double precision theta(d),extr(d)
        do 2 j=1,nmc
           do 1 k=j,nmc
              extr(psi(mc(j),mc(k)))=theta(psi(mc(j),mc(k)))
1          continue
2       continue
        return
        end
C***********************************************************************
	subroutine chols(d,theta,p,psi,mc,nmc)
	integer d,p,psi(0:p,0:p),mc(p),nmc
	double precision theta(d),tmp
	do 40 i=1,nmc
	  tmp=0.
	  do 10 k=1,i-1
	    tmp=tmp+theta(psi(mc(k),mc(i)))**2
10	  continue
	  theta(psi(mc(i),mc(i)))=sqrt(theta(psi(mc(i),mc(i)))-tmp)
	  do 30 j=i+1,nmc
	    tmp=0.
	    do 20 k=1,i-1
	      tmp=tmp+theta(psi(mc(k),mc(i)))*theta(psi(mc(k),mc(j)))
20	    continue
	    theta(psi(mc(i),mc(j)))=(theta(psi(mc(i),mc(j)))-tmp)
     /             /theta(psi(mc(i),mc(i)))
30	  continue
40	continue
        return
	end
C***********************************************************************
	function gauss()
        save alt,next
	integer alt
C real next
	double precision next
	data pi/3.141593/
        if((alt.ne.0).and.(alt.ne.1)) alt=0
	if(alt.eq.0)then
	  u1=rangen(0)
	  u2=rangen(0)
	  gauss=sqrt(-2*log(u1))*cos(2*pi*u2)
	  next=sqrt(-2*log(u1))*sin(2*pi*u2)
	  alt=1
	else
	  gauss=next
	  alt=0
	endif
	return
	end
C***********************************************************************
        subroutine ps1n(d,t,m,tau,theta,p,psi,n,mat,z,b,c)
C Performs P-step of data augmentation. Theta initially contains mu0 
C and lambdainverse in packed storage, and sufficient statistics are in
C t. The new parameter is written to theta. Requires workspaces 
C mat(p,p), z(p), b(d), and c(p).
        integer d,p,psi(0:p,0:p),n,c(p)
C real z(p)
        double precision z(p)
        double precision m,tau,t(d),theta(d),mat(p,p),b(d)
        call swp(d,t,0,p,psi,p,1)
        do 20 j=1,p
           do 10 k=j,p
              theta(psi(j,k))=theta(psi(j,k))+dble(n)*t(psi(j,k))
     /     + tau*dble(n)/(tau+dble(n))*
     /     (t(psi(0,j))-theta(psi(0,j)))*(t(psi(0,k))-theta(psi(0,k)))
10         continue
20      continue
        do 30 j=1,p
           theta(psi(0,j))=dble(n)*t(psi(0,j))/(dble(n)+tau)
     /     + tau*theta(psi(0,j))/(dble(n)+tau)
30      continue
        tau=tau+dble(n)
        m=m+dble(n)
        call ninvwn(d,theta,tau,m,p,psi,mat,z,b,c)
        return
        end
C***********************************************************************
        subroutine invtrn(d,t,p,psi)
C Inverts triangular matrix in packed storage
        integer d,p,psi(0:p,0:p)
        double precision t(d),sum
        t(psi(1,1))=1./t(psi(1,1))
        do 10 k=2,p
           t(psi(k,k))=1./t(psi(k,k))
           do 5 j=1,k-1
              sum=0
              do 3 i=j,k-1
                 sum=sum+t(psi(i,j))*t(psi(i,k))
3             continue
              t(psi(j,k))=-sum*t(psi(k,k))
5          continue
10      continue
        return
        end
C***********************************************************************
        subroutine bfac(d,b,p,psi,m)
C draws triangular square-root of a Wishart(m,I) using Bartlett
C decomposition, putting result into packed storage
        integer d,p,psi(0:p,0:p)
        double precision b(d),m
	do 10 j=1,p
	   b(psi(j,j))=dble(sqrt(2.*gamm((dble(m)-float(j)+1.)/2.)))
10      continue
        do 30 j=1,p-1
	   do 20 k=j+1,p
             b(psi(j,k))=dble(gauss())
20         continue
30      continue
        return
        end
C***********************************************************************
        subroutine mmn(d,l,u,p,psi,m)
C Multiplies lower triangular matrix l by upper-triangular matrix u,
C both in packed storage, and puts result into m which is unpacked
        integer d,p,psi(0:p,0:p)
        double precision l(d),u(d),sum,m(p,p)
        do 10 i=1,p
	   do 5 j=1,p
	     sum=0
             do 2 k=1,min(i,j)
	        sum=sum+l(psi(i,k))*u(psi(k,j))
2            continue
             m(i,j)=sum
5          continue
10      continue
        return
        end
C***********************************************************************
        subroutine ninvwn(d,pars,tau,m,p,psi,mat,z,b,c)
C Draws from normal-inverted Wishart distribution with parameters 
C given by pars (contains mu0 and lambdainverse, in packed storage),
C tau, and m. Needs mat(p,p), z(p), b(d), c(p) as workspaces. 
C Overwrites pars with result.
        integer d,p,psi(0:p,0:p),c(p)
        double precision pars(d),tau,m,mat(p,p),b(d),sum
C real junk, z(p)
        double precision junk, z(p)
        junk=gauss()
        do 2 j=1,p
           c(j)=j
2       continue
        call chols(d,pars,p,psi,c,p)
        call bfac(d,b,p,psi,m)
        call invtrn(d,b,p,psi)
        call mmn(d,b,pars,p,psi,mat)
C add normal errors to mu0
        do 5 j=1,p
           z(j)=gauss()
5       continue
        do 20 j=1,p
           sum=0
           do 10 k=1,p
              sum=sum+mat(k,j)*z(k)
10         continue
           pars(psi(0,j))=pars(psi(0,j))+sum/sqrt(tau)
20      continue
C put sigma=mat^t mat into pars
        do 60 i=1,p
           do 50 j=i,p
              sum=0
              do 40 k=1,p
                 sum=sum+mat(k,i)*mat(k,j)
40            continue
              pars(psi(i,j))=sum
50         continue
60      continue
        pars(psi(0,0))=-1.
        return
        end
C***********************************************************************
        subroutine ph2thn(d,theta,p,psi)
C Converts phi to theta
        integer d,p,psi(0:p,0:p)
        double precision theta(d)
        do 1 j=1,p-1
           call swp(d,theta,j,p,psi,j,1)
1       continue
        do 2 j=1,p-1
           call swp(d,theta,j,p,psi,p,-1)
2       continue
        return
        end
C***********************************************************************
        subroutine sjn(p,npatt,r,sj)
C computes s_j, the number of the last missingness pattern for which
C the jth variable needs to be imputed to complete the monotone pattern
        integer p,npatt,r(npatt,p),sj(p),patt,tmp
        do 30 j=1,p
           patt=npatt+1
10         continue
           patt=patt-1
           if(patt.ge.1)then
              if(r(patt,j).eq.0)goto 10
           endif
           sj(j)=patt
30      continue
        tmp=sj(p)
        do 40 j=p-1,1,-1
           sj(j)=max0(sj(j),tmp)
           tmp=sj(j)
40      continue
        return
        end
C***********************************************************************
        subroutine nmons(p,npatt,r,nmdp,sj,nmon)
C Computes the number of observations in (Xobs,Xmis*) for each column
C (called N_j in Figure 6.11 of book)
        integer p,npatt,r(npatt,p),nmdp(npatt),sj(p),nmon(p),patt
        do 40 j=1,p
           nmon(j)=0
           do 20 patt=1,sj(j)
              nmon(j)=nmon(j)+nmdp(patt)
20         continue
40      continue
        return
        end
C***********************************************************************
        subroutine lasts(p,npatt,sj,last)
C Finds last variable in each missingness pattern
C to complete a monotone pattern
        integer p,npatt,sj(p),last(npatt),patt,start
        do 50 j=p,1,-1
           if(j.eq.p)then
             start=1
           else
              start=sj(j+1)+1
           endif
           do 40 patt=start,sj(j)
              last(patt)=j
40         continue
50      continue
        return
        end
C***********************************************************************
        subroutine is2n(d,theta,p,psi,n,x,npatt,r,mdpst,nmdp,sj,last,
     /    oc,mc,z,c)
C Performs I-step of monotone data augmentation. Randomly draws missing 
C data Xmis* given theta to complete the monotone pattern.
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,mc(p),nmc,patt,sj(p),last(npatt)
C real x(n,p),z(p),junk
        double precision x(n,p),z(p),junk
        double precision theta(d),c(d)
        junk=gauss()
        do 200 patt=1,npatt
           call swpobs(d,theta,p,psi,npatt,r,patt)
           call gtmc(p,npatt,r,patt,mc,nmc,last(patt))
           call gtoc(p,npatt,r,patt,oc,noc,last(patt))
           call sigex(d,theta,c,p,psi,mc,nmc)
           call chols(d,c,p,psi,mc,nmc)
           do 150 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 100 j=1,nmc
                 x(i,mc(j))=theta(psi(0,mc(j)))
                 do 40 k=1,noc
                    x(i,mc(j))=x(i,mc(j))+theta(psi(oc(k),
     /                 mc(j)))*x(i,oc(k))
40               continue
                 z(mc(j))=gauss()
                 do 50 k=1,j
                    x(i,mc(j))=x(i,mc(j))+z(mc(k))*c(psi(mc(j),mc(k)))
50               continue
100           continue
150        continue
200     continue
        return
        end
C***********************************************************************
        subroutine layers(p,sj,layer,nlayer)
C Finds layers for observed parts of the sufficient statistics
        integer p,sj(p),layer(p),nlayer
        nlayer=0
        do 10 j=p,1,-1
           if(j.eq.p)then
              if(sj(j).gt.0) nlayer=nlayer+1
           else
              if(sj(j).gt.sj(j+1)) nlayer=nlayer+1
           endif
           layer(j)=nlayer
10      continue
        return
        end
C***********************************************************************
        subroutine tobsmn(p,psi,n,x,npatt,r,mdpst,nmdp,last,oc,sj,
     /     layer,nlayer,d,tobs)
C Calculates known parts of the sufficient statistics for monotone
C data augmentation
        integer p,psi(0:p,0:p),d,n,npatt,r(npatt,p),mdpst(npatt),
     /     nmdp(npatt),last(npatt),oc(p),noc,layer(p),nlayer,l,patt,
     /     fpatt,lpatt,sj(p)
C real x(n,p)
        double precision x(n,p)
        double precision tobs(nlayer,d)
        do 2 l=1,nlayer
           do 1 i=1,d
              tobs(l,i)=dble(0.)
1          continue
2       continue
        lpatt=0
        do 100 l=1,nlayer
           fpatt=lpatt+1
           j=p+1
5          continue
             j=j-1
             if(layer(j).ne.l)goto 5
           lpatt=sj(j)
           do 80 patt=fpatt,lpatt
              call gtoc(p,npatt,r,patt,oc,noc,last(patt))
              do 30 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
                 tobs(l,psi(0,0))=tobs(l,psi(0,0))+dble(1.)
                 do 20 j=1,noc
                    tobs(l,psi(0,oc(j)))=tobs(l,psi(0,oc(j)))
     /                 +x(i,oc(j))
                    do 10 k=j,noc
                       tobs(l,psi(oc(j),oc(k)))=tobs(l,psi(oc(j),
     /                    oc(k)))+x(i,oc(j))*x(i,oc(k))
10                  continue
20               continue
30            continue
80         continue
100     continue
        return
        end
C***********************************************************************
	subroutine chol2(d,theta,p,psi,last)
	integer d,p,psi(0:p,0:p),last
	double precision theta(d),tmp
	do 40 i=0,last
	  tmp=0.
	  do 10 k=0,i-1
	    tmp=tmp+theta(psi(k,i))**2
10	  continue
	  theta(psi(i,i))=sqrt(theta(psi(i,i))-tmp)
	  do 30 j=i+1,last
	    tmp=0.
	    do 20 k=0,i-1
	      tmp=tmp+theta(psi(k,i))*theta(psi(k,j))
20	    continue
	    theta(psi(i,j))=(theta(psi(i,j))-tmp)/theta(psi(i,i))
30	  continue
40	continue
        return
	end
C***********************************************************************
        subroutine ps2n(p,psi,n,x,npatt,r,mdpst,nmdp,oc,mc,nmon,
     /     sj,nlayer,d,tobs,t,c,v,theta)
C P-step of monotone data augmentation
        integer p,psi(0:p,0:p),d,n,npatt,r(npatt,p),mdpst(npatt),
     /     nmdp(npatt),oc(p),noc,mc(p),nmc,nlayer,patt,nmon(p),
     /     sj(p),lsj,l
C real x(n,p),df
        double precision x(n,p),df
        double precision tobs(nlayer,d),t(d),c(d),v(0:p),theta(d)
        do 1 i=1,d
           t(i)=dble(0.)
1       continue
        lsj=0
        l=0
        do 100 j=p,1,-1
           if(sj(j).gt.lsj) then
              l=l+1
              t(psi(0,0))=t(psi(0,0))+tobs(l,psi(0,0))
              do 10 k=1,j
                 t(psi(0,k))=t(psi(0,k))+tobs(l,psi(0,k))
                 do 5 m=k,j
                    t(psi(m,k))=t(psi(m,k))+tobs(l,psi(m,k))
5                continue
10            continue
              do 20 patt=(lsj+1),sj(j)
                 call gtmc(p,npatt,r,patt,mc,nmc,j)
                 call gtoc(p,npatt,r,patt,oc,noc,j)
                 do 18 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
                    do 15 k=1,nmc
                       t(psi(0,mc(k)))=t(psi(0,mc(k)))+x(i,mc(k))
                       do 12 m=1,noc
                          t(psi(mc(k),oc(m)))=t(psi(mc(k),oc(m)))+
     /                       x(i,mc(k))*x(i,oc(m))
12                     continue
                       do 13 m=1,k
                          t(psi(mc(k),mc(m)))=t(psi(mc(k),mc(m)))+
     /                         x(i,mc(k))*x(i,mc(m))
13                     continue
15                  continue
18               continue
20            continue
           endif
           if(sj(j).gt.lsj) then
              do 25 k=0,j-1
                 call swp(d,t,k,p,psi,j,1)
25            continue
           endif
           df=float(nmon(j)+3*(p-j)-1)
           theta(psi(j,j))=t(psi(j,j))/dble(chisq(df))
           do 30 k=0,j-1
              do 28 m=k,j-1
                 c(psi(k,m))=-theta(psi(j,j))*t(psi(k,m))
28            continue
30         continue
           call chol2(d,c,p,psi,j-1)
           do 50 k=0,j-1
              v(k)=dble(gauss())
              theta(psi(k,j))=t(psi(k,j))
              do 40 m=0,k
                 theta(psi(k,j))=theta(psi(k,j))+c(psi(m,k))*v(m)
40            continue
50         continue
           if(j.gt.1)then
              if(sj(j-1).gt.sj(j))then
                 do 60 k=0,j-1
                    call swp(d,t,k,p,psi,j-1,-1)
60               continue
              elseif(sj(j-1).eq.sj(j))then
                 call swp(d,t,j-1,p,psi,j-1,-1)
              endif
           endif
           lsj=sj(j)
100     continue
        theta(psi(0,0))=dble(-1.)
        call ph2thn(d,theta,p,psi)
        return
        end
C***********************************************************************
