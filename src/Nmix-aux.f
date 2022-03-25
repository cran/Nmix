	subroutine stdalloc(y,n,wt,mu,ssq,ncmax,start,leng,next,pw,inext,
     &		first,delta,qprior)

	integer first,start
	real mu
	logical qprior
	dimension start(ncmax),leng(ncmax),inext(n),next(ncmax),
     &		pw(ncmax),y(n),wt(ncmax),mu(ncmax),ssq(ncmax)

	j = first
	do while(j.ne.0)
 	    start(j) = 0
	    leng(j) = 0
	    j = next(j)
	end do

	do i = 1,n

	pwsum = 0.0

	j = first
	do while(j.ne.0)
	if(qprior) then
	    pw(j) = wt(j)
	else
	    pw(j) = wt(j)
     &		*exp(max(-20.0,-0.5*(y(i)-mu(j))**2/ssq(j)))/sqrt(ssq(j))
	end if
	    pwsum = pwsum+pw(j)
	    j = next(j)
	end do
	usd = sdrand(u)
	u = pwsum*usd

	j = first
	do while(j.ne.0)
	    u = u-pw(j)
	    if(u.lt.0.0) go to 43
	    j = next(j)
	end do

	j = first

43	inext(i) = start(j)
	start(j) = i
	leng(j) = leng(j)+1

	end do

	return

	end

c---------------------------------------------------

	subroutine reorder(mu,ncmax,next,prev,first)

c.. check ordering of mu's and correct if necessary

	integer prev,first
	real mu

	dimension mu(ncmax),prev(ncmax),next(ncmax)

	jp = first
73	jnext = next(jp)
	jq = prev(jp)
74	if(jq.eq.0) go to 75
	if(mu(jq).lt.mu(jp)) go to 75
	jq = prev(jq)
	go to 74
75	if(prev(jp).ne.jq) then
		if(prev(jp).le.0) return
c.. remove jp from list
		next(prev(jp)) = jnext
		if(jnext.ne.0) prev(jnext) = prev(jp)
c.. re-insert jp, after jq
		prev(jp) = jq
		if(jq.ne.0) then
			next(jp) = next(jq)
			if(next(jq).le.0) return
			prev(next(jq)) = jp
			next(jq) = jp
		else
			next(jp) = first
			prev(first) = jp
			first = jp
		end if
	end if
	jp = jnext
	if(jp.ne.0) go to 73

	return

	end

c---------------------------------------------------

      subroutine sort(na,llo,lhi,q,nq)
      dimension na(lhi),q(nq)
c
c    sorts array na(l), l = llo to lhi
c    according to values of q(na(l))
c    so if llo = 1 and na(l) = l on entry,
c    obtain na(l) = index of l-th smallest q on exit.
c

      mesh = lhi-llo+1
      do while(mesh.gt.1)
         mesh = (mesh+1)/3
         lst = llo+mesh
         do l1 = lst,lhi
      	nal = na(l1)
      	ql = q(nal)
      	lnow = l1
      	do l2 = lst,l1,mesh
      	   lnext = lnow-mesh
      	   nan = na(lnext)
      	   if(ql.ge.q(nan)) exit
      	   na(lnow) = nan
      	   lnow = lnext
		end do
         na(lnow) = nal
	   end do
      end do
      return
      end

c---------------------------------------------------

	function entropy(n,leng,ncmax,first,next)

c.. calculates ent = -\sum_j (n_j/n) \log (n_j/n)

	integer first,leng(ncmax),next(ncmax)

	entropy = log(real(n))

	j = first
	do while(j.ne.0)
		if(leng(j).ne.0) then
		    temp = real(leng(j))
		    entropy = entropy-(temp/n)*log(temp)
		end if
		j = next(j)
	end do

	return

	end

c---------------------------------------------------

	subroutine rmu(mu,y,xi,kappa,alpha,beta)

c		mu = rv from f propto
c                   exp(-kappa(mu-xi)^2/2)/(beta+.5*(yi-mu)^2)^(alpha+.5)

	real mu,kappa

1	xx = sdrand(u)
	yy = xx*y+0.2*(sdrand(u)-0.5)
	mu = yy/xx
	if(xx**2.gt.exp(-0.5*kappa*(mu-xi)**2)/
     &	    (1.0+(0.5/beta)*(y-mu)**2)**(alpha+0.5)) go to 1

	return

	end

c---------------------------------------------------

	subroutine rbeta2(a1,a2,z)
	external rgamma
c	x = rgamma(a1)
c	y = rgamma(a2)
	call rgamma2(a1,x)
	call rgamma2(a2,y)
	z = x/(x+y)
	return
	end

c---------------------------------------------------

	real function logbeta(x1,x2)
	integer, parameter :: dbl = kind(1.0d0)
	real(kind=dbl) :: dlgama
	logbeta = sngl(dlgama(dble(x1))+dlgama(dble(x2))
     &		-dlgama(dble(x1+x2)))
	return
	end
