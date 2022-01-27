	subroutine nmixsub(y,n,pw,den,avn,ncmax1,ncd1,ngrid1,ncmax21,
     &	wtav,muav,sigav,acctry,enttr,ktr,off,partr,nparsamp,devtr,
c.. input, dimensions (6)
     &	nmax,ncmax,ncmax2,ncd,ngrid,kkzz,
c.. input, read only:
c.. sampling pars (5): 
     &	nsweep,nburnin,imoves,nstep,kinit,
c.. switches (13): 
     &	qempty,qprior,qunif,qfix,qfull,qrkpos,qrange,qkappa,
     &	qbeta,qpalloc,qpwms,qpclass,qkreg,
c.. model pars (11): 
     &	alpha,delta,eee,fff,ggg,hhh,unhw,     
     &	kappa,lambda,xi,sp,
c.. output control (4): 
     &	nspace,nsamp,nskdel,idebug,
c.. input, but not read-only (2):
     &	beta,qdebug)

c  arguments not (even implicitly) mentioned on Nmix webpage
c  ncmax,ncmax2,ncd,ngrid,kkzz,
c  qunif,unhw,
c  qrkpos,
c  qkreg,nskdel,
c  idebug,qdebug

	character*6 num,num2
	character*80 line
	character moveselect,sep
	character*5 movnam(8)
	character*8 movabb

	logical qbeta,qkappa,qq,qdebug,qprior,qrange,
     &		qfull,qfix,qpalloc,
     &		qunif,qempty,qstop,qpwms,qpclass,
     &		qkreg,qrkpos

	integer hwm,first,free,st1,st2,sp,ckrep
	integer imoves(nstep),acctry(2,4)

	real kappa,lambda,mu1,mu2,muc,ms,logratio,logbeta,lprob
	real kwas,loglikrat,mustar(1)

	INTEGER, PARAMETER :: dbl = KIND(1.0d0)
	REAL(KIND=dbl) :: dlgama

	real y(nmax),den(ncd,ngrid),pw(ncmax)
	real wtav(ncmax2,ncmax2),muav(ncmax2,ncmax2),sigav(ncmax2,ncmax2)
	real avn(ncmax21,ncmax21)
	real enttr(nsweep/nspace),devtr(nsweep/nspace)
	real partr(3,nparsamp)

	integer ktr(nsweep/nspace),off(nsweep/nspace),off0

	integer split,combine,allocate,weights,parameters,hyper,
     &		birth,death
	data split,combine,allocate,weights,parameters,hyper,
     &		birth,death/1,2,3,4,5,6,7,8/
	data movnam/'split','comb ','alloc','wts  ','param','hyper',
     &		'birth','death'/
	data movabb/'scawphbd'/

	integer, dimension(:), allocatable :: start,leng,prev,next,count,
     &	countpos
	integer, dimension(:), allocatable :: inext,na,z
	integer, dimension(:), allocatable :: nfkemp,countden

	real, dimension(:), allocatable :: lp,yv
	real, dimension(:), allocatable :: wt,mu,mun,bf,avdev,avdevc,pf
	real, dimension(:), allocatable :: ssq
	real, dimension(:,:), allocatable :: pclass,pz
	real, dimension(:), allocatable :: ggh,b,d

	logical, dimension(:), allocatable :: fileopen

	pi = 4.0*atan(1.0)
	rr2pi = 1.0/sqrt(2.0*pi)

	allocate(start(ncmax),leng(ncmax))
	allocate(prev(ncmax),next(ncmax))
	allocate(count(ncmax),countpos(ncmax))
	allocate(inext(nmax))
	allocate(na(nmax),z(nmax))
	allocate(nfkemp(0:9))
	allocate(countden(ncd1))

	allocate(lp(ncmax))
	allocate(yv(nmax))
	allocate(wt(ncmax))
	allocate(mu(ncmax),mun(ncmax))
	allocate(bf(ncmax))
	allocate(avdev(ncmax1),avdevc(ncmax1))
	allocate(pf(ncmax))
	allocate(ssq(ncmax))
	allocate(fileopen(ncmax))


	allocate(pclass(kkzz,ngrid),pz(kkzz,nmax))
	allocate(ggh(ncd1))
	allocate(b(ncmax),d(ncmax))

c	call realpr('y',1,y,n)

c.. sort out filename stuff

	sep = "/"
	ijk = 1

1	istd = 6-int(log10(0.5+ijk))
c	write(num,'(i6)') ijk

c	inquire(file=base(1:ncbase)//sep//num(istd:6)//".log",exist=qq)
c	if(qq) then
c		ijk = ijk+1
c		go to 1
c	end if

c	prfx = base(1:ncbase)//sep//num(istd:6)
c	npf = ncbase+8-istd

c	open(4,file=prfx(1:npf)//".out",status='unknown')
c	open(7,file=prfx(1:npf)//".log",status='unknown')
c	if(.not.qfix) open(8,file=prfx(1:npf)//".bdlog",status='unknown')

c	if(qfull) then
c	    open(9,file=prfx(1:npf)//".dev",status='unknown')
c	    if(qpclass) then
c		open(14,file=prfx(1:npf)//".pcl",status='unknown')
c		open(15,file=prfx(1:npf)//".scl",status='unknown')
c	    end if
c	end if
c	if(qbeta.or.qkappa)
c     &	    open(10,file=prfx(1:npf)//".bk",status='unknown')
c	open(12,file=prfx(1:npf)//".k",status='unknown')
c	if(qfull)
c     &		open(13,file=prfx(1:npf)//".den",status='unknown')
c	open(17,file=prfx(1:npf)//".ent",status='unknown')
c	open(11,file=prfx(1:npf)//".pe",status='unknown')
c	if(qfull) then
c	    open(16,file=prfx(1:npf)//".avn",status='unknown')
c	end if

c	if(idebug.le.nsweep) then
c	    open(18,file=prfx(1:npf)//".db",status='unknown')
c	end if

c.. proposal parameters

	ws = 2.0
	ms = 2.0
	ss = 1.0

	do k = 2,ncmax-1
		b(k) = 0.5
		d(k) = 0.5
	end do
	b(ncmax) = 0.0
	d(ncmax) = 1.0
	b(1) = 1.0
	d(1) = 0.0

c.. basic statistics

	ymin = y(1)
	ymax = ymin
	do i = 2,n
		ymin = min(ymin,y(i))
		ymax = max(ymax,y(i))
	end do

	ymid = 0.5*(ymin+ymax)
	yrange = ymax-ymin
	yd0 = ymin-0.05*yrange
	ydinc = 1.1*yrange/ngrid

	ysum = 0.0
	do i = 1,n
		ysum = ysum+y(i)
	end do
	ysum = ysum/n
	ssd  = 0.0
	do i = 1,n
		ssd = ssd+(y(i)-ysum)**2
	end do

c	call realpr('ssd',1,ssd,1)

c.. adjust hyperparameters if range-based

	if(qrange) then
		xiwas = xi
		xi = ymid+xi*0.5*yrange
		kwas = kappa
		kappa = kappa/yrange**2
		unhwwas = unhw
		unhw = unhw*0.5*yrange
		bwas = beta
		beta = beta*yrange**2
		fwas = fff
		fff = fff*yrange**2
		hwas = hhh
		hhh = hhh/yrange**2
	end if

c.. set up prior on k

	    if(lambda.gt.0.0) then
		temp = -log(exp(lambda)-1.0)
		do k = 1,ncmax
			temp = temp+log(lambda/k)
			lp(k) = temp
		end do
	    else if(lambda.lt.0.0) then
		do k = 1,ncmax
			lp(k) = 0.0
		end do
	    else
		do k = 1,ncmax
			lp(k) = log(1.0/k)
		end do
	    end if

c.. log parameter values

c	write(0,*) "Run: ",prfx(1:npf)," : ",line(1:ncline)
c	write(7,*) "Run: ",prfx(1:npf)," : ",line(1:ncline)
c	write(7,'("random number seed:",i12)') iseed
c	if(qprior) write(7,'("PRIOR simulation")')
c	if(qunif) write(7,'("Uniform prior for component means")')
c	if(qfix) write(7,'("fixed k")')
c	if(qfull) write(7,'("full output")')
c	if(qrkpos) write(7,'("reporting kpos option")')
c	if(qrange) then
c		write(7,
c     &	 '("xi      =   ",g11.4,"*R/2+ymid =",g11.4)') xiwas,xi
c	else
c		write(7,'("xi      =   ",g11.4)') xi
c	end if
c	if(qunif) then
c		if(qrange) then
c			write(7,
c     &	'(" unhw = ",g11.4,"*R/2 =",g11.4)') unhwwas,unhw
c		else
c			write(7,'(" unhw = ",g11.4)') unhw
c		end if
c	else if(qkappa) then
c		write(7,'("kappa is Gamma(e,f) with")')
c		write(7,'("  e     =   ",g11.4)') eee
c		if(qrange) then
c			write(7,
c     &	 '("  f     =   ",g11.4,"*R^2 =",g11.4)') fwas,fff
c		else
c			write(7,'("  f     =   ",g11.4)') fff
c		end if
c	else
c		if(qrange) then
c			write(7,
c     &	 '("kappa   =   ",g11.4,"/R^2 =",g11.4)') kwas,kappa
c		else
c			write(7,'("kappa   =   ",g11.4)') kappa
c		end if
c	end if
c	if(sp.gt.1) write(7,'(" spacing prior, sp=",i3)') sp
c	write(7,'("alpha   =   ",g11.4)') alpha
c	if(qbeta) then
c		write(7,'("beta is Gamma(g,h) with")')
c		write(7,'("  g     =   ",g11.4)') ggg
c		if(qrange) then
c			write(7,
c     &	 '("  h     =   ",g11.4,"/R^2 =",g11.4)') hwas,hhh
c		else
c			write(7,'("  h     =   ",g11.4)') hhh
c		end if
c	else
c		if(qrange) then
c			write(7,
c     &	 '("beta   =   ",g11.4,"*R^2 =",g11.4)') bwas,beta
c		else
c			write(7,'("beta    =   ",g11.4)') beta
c		end if
c	end if
c	write(7,'("delta   =   ",g11.4)') delta
c	write(7,'("lambda  =   ",g11.4)') lambda
c	write(7,'("nsweep  =   ",i8)') nsweep
c	write(7,'("nstep   =   ",i8)') nstep
c	write(7,'("nburnin =   ",i8)') nburnin
c	write(7,'("n       =   ",i8)') n
c	write(7,'("ymin    =   ",g11.4)') ymin
c	write(7,'("ymax    =   ",g11.4)') ymax

c	write(7,*) "move schedule: ",moves(1:nstep)

c.. initialise dynamic variables
c.. initialise number of components, k, and allocation

c	call intpr('kinit',5,kinit,1)

	if(kinit.eq.0) then
c		open(3,file=base(1:ncbase)//".zin",status='unknown')
c		read(3,*) (na(i),i=1,n)
c		close(3)
		do i = 1,n
			if(na(i).lt.1.or.na(i).gt.ncmax) return
		end do
		do j = 1,ncmax
			start(j) = 0
			leng(j) = 0
		end do
		do i = 1,n
			j = na(i)
			inext(i) = start(j)
			start(j) = i
			leng(j) = leng(j)+1
		end do
		k = 0
		first = 0
		free = 0
		do j = 1,ncmax
		    next(j) = 0
		    if(leng(j).gt.0) then
			hwm = j
			k = k+1
			if(first.eq.0) then
				first = j
				prev(j) = 0
			else
				prev(j) = jp
				next(jp) = j
			end if
			jp = j
		    else
			if(free.eq.0) then
				free = j
			else
				next(jq) = j
			end if
			prev(j) = -1
			jq = j
		    end if
		end do
c	write(7,'("initial allocation read from .zin file: k =",i3)') k
	else
c	write(7,'("initial number of components:",i4)') kinit
	if(kinit.eq.1) then
	    start(1) = 1
	    leng(1) = n
	    do i = 1,n
		inext(i) = i+1
	    end do
	    inext(n) = 0
	else
	    do i = 1,n
		na(i) = i
	    end do
	    call sort(na,1,n,y,n)
	    do j = 1,kinit
		l1 = (n*(j-1))/kinit+1
		l2 = (n*j)/kinit
		leng(j) = l2-l1+1
		if(leng(j).gt.0) then
			start(j) = na(l1)
		else
			start(j) = 0
		end if
		do l = l1,l2-1
			inext(na(l)) = na(l+1)
		end do
		if(l2.gt.0) inext(na(l2)) = 0
	    end do
	end if
	    first = 1
	    do j = 1,kinit
		prev(j) = j-1
		next(j) = j+1
	    end do
	    next(kinit) = 0
	    free = kinit+1
	    do j = kinit+1,ncmax
		prev(j) = -1
		next(j) = j+1
	    end do
	    next(ncmax) = 0
	    k = kinit
	    hwm = kinit
	end if

c	call intpr('hwm',3,hwm,1)

	kemp = 0
	j = first
	do while(j.ne.0)
	    if(leng(j).eq.0) kemp = kemp+1
	    j = next(j)
	end do

c	call intpr('kemp',4,kemp,1)

c.. initialise weights

	wtsum = 0.0
	j = first
	do while(j.ne.0)
c	    wt(j) = rgamma(delta+leng(j))
		call rgamma2(delta+leng(j),wt(j))
c	call realpr('wt(j)',5,wt(j),1)
	    wtsum = wtsum+wt(j)
	    j = next(j)
	end do
	j = first
	do while(j.ne.0)
	    wt(j) = wt(j)/wtsum
	    j = next(j)
	end do

c	call realpr('wtsum',5,wtsum,1)

c.. initialise means

	call gauss(mu,hwm)

c	call realpr('mu',2,mu,hwm)

c	if(qdebug) then
c		write(18,*) 'mu',(mu(j),j=1,hwm)
c	end if
	ssq0 = beta/alpha
	j = first
	do while(j.ne.0)
	if(qprior) then
	    if(qunif) then
		mu(j) = xi+unhw*(2.0*sdrand(u)-1.0)
	    else
		mu(j) = xi+mu(j)/sqrt(kappa)
	    end if
	else
	    ysum = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ysum = ysum+y(i)
		i = inext(i)
	    end do
	    if(qunif) then
		if(leng(j).eq.0) then
		    mu(j) = xi+unhw*(2.0*sdrand(u)-1.0)
		else
44		    mu(j) = ysum/leng(j)+sqrt(ssq0/leng(j))*mu(j)
		    if(abs((mu(j)-xi)/unhw).gt.1.0) then
			call gauss(mu(j),1)
			go to 44
		    end if
		end if
	    else
		con = 1.0/(leng(j)/ssq0+kappa)
		mu(j) = con*(ysum/ssq0+xi*kappa)+sqrt(con)*mu(j)
	    end if
c	call realpr('mu(j)',5,mu(j),1)
	end if
	    j = next(j)
	end do

c.. initialise variances

	j = first
	do while(j.ne.0)
	if(qprior) then
c	    ssq(j) = beta/rgamma(alpha)
		call rgamma2(alpha,temp)
		ssq(j) = beta/temp
	else
	    ssd = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ssd = ssd+(y(i)-mu(j))**2
		i = inext(i)
	    end do
c	    ssq(j) = (beta+0.5*ssd)/rgamma(alpha+0.5*leng(j))
		call rgamma2(alpha+0.5*leng(j),temp)
		ssq(j) = (beta+0.5*ssd)/temp
	end if
c	call realpr('ssq(j)',6,ssq(j),1)

	j = next(j)
	end do

c.. now check ordering of mu, and correct

	call reorder(mu,ncmax,next,prev,first)
c	call realpr('mu',2,mu,hwm)

c.. initialise v
	do i = 1,n
		yv(i) = y(i)
	end do
c	call realpr('yv',2,yv,n)

c.. initialise accumulators

	do j = 1,ncmax
		count(j) = 0
		countpos(j) = 0
		fileopen(j) = .false.
		avdev(j) = 0.0
		avdevc(j) = 0.0
	end do

	do j = 1,ncmax2
	do ij = 1,j
		wtav(j,ij) = 0.0
		muav(j,ij) = 0.0
		sigav(j,ij) = 0.0
	end do
	end do

	if(qfull) then
		do j = 1,ncmax2
		do ij = 1,j
			avn(j,ij) = 0.0
		end do
		end do
		do krec = 1,ncd
	    		countden(krec) = 0
	    		do iyd = 1,ngrid
			den(krec,iyd) = 0.0
	    		end do
		end do
	end if

	nkemp = 0
	avkemp = 0.0
	ppkemp = 0.0
	do j = 0,9
		nfkemp(j) = 0
	end do

	if(qpclass) then
	    do iyd = 1,ngrid
	    do j = 1,kkzz
		pclass(j,iyd) = 0.0
	    end do
	    end do

	    do i = 1,n
	    do j = 1,kkzz
		pz(j,i) = 0.0
	    end do
	    end do
	end if

	ntrys = 0
	ntryc = 0
	ntryb = 0
	ntryd = 0
	naccs = 0
	naccc = 0
	naccb = 0
	naccd = 0
	nrejr = 0

	if(qunif) then
	    const1 = -sngl(dlgama(dble(alpha)))
	else
	    const1 = 0.5*log(0.5/pi)-sngl(dlgama(dble(alpha)))
	end if
	const2 = logbeta(ws,ws)+logbeta(ms,ms)+logbeta(ss,ss)
c	call realpr('const2',6,const2,1)
c	if(idebug.lt.0) write(18,'("consts",2g11.5)') const1,const2

	do j = 1,ncmax
	temp = sp*(j+1)
	    do jj = 1,sp-1
		temp = temp*real(sp*(j+1)+jj)/real(jj)
	    end do
	pf(j) = temp
	end do
c	call realpr('pf',2,pf,ncmax)

	off0 = 0

c	if(qdebug) then
c		write(18,'("initialised")')
c		write(18,*) 'k',k
c		j = first
c		do while(j.ne.0)
c		    write(18,*) 'j,leng,..',j,leng(j),mu(j),ssq(j),wt(j)
c		    j = next(j)
c		end do
c	end if

c.. main loop: MCMC iteration

	do isweep = 1-nburnin,nsweep

c	if(mod(nsweep-isweep,nsweep/10).eq.0) then
c	    write(0,'(i3,$)') (nsweep-isweep)/(nsweep/10)
c	end if

	qdebug = isweep.ge.idebug

	do istep = 0,nstep-1
c	call intpr('istep',5,istep,1)
c.. select next move type

c	moveselect = moves(istep+1:istep+1)
	im = imoves(istep+1)
	moveselect = movabb(im:im)
	if(moveselect.eq.'s') then
	    usw = sdrand(u)
c	    if(qdebug) write(18,*) 'isweep,usw',isweep,usw
	    if(usw.lt.b(k)) then
		move = split
	    else
		move = combine
	    end if
	else if(moveselect.eq.'w') then
		move = weights
	else if(moveselect.eq.'p') then
		move = parameters
	else if(moveselect.eq.'a') then
		move = allocate
	else if(moveselect.eq.'h') then
		move = hyper
	else if(moveselect.eq.'b') then
	    usw = sdrand(u)
	    if(usw.lt.b(k)) then
		move = birth
	    else
		move = death
	    end if
	end if

	if(move.eq.split) then

c	   write(18,*)'split'

c.. split move ----------------------------------------------

	ntrys = ntrys+1

20	j = 1+int(hwm*sdrand(u))
	if(prev(j).lt.0) go to 20

c	if(qdebug) write(18,*) 'split:isweep,j',isweep,j

	wtc = wt(j)
	muc = mu(j)
	ssqc = ssq(j)

c	u1 = rbeta(ws,ws)
	call rbeta2(ws,ws,u1)
	cu1 = 1.0-u1
c	u2 = rbeta(ms,ms)
	call rbeta2(ms,ms,u2)
c	if(qdebug) write(18,*) 'split:isweep,u1,u2',isweep,u1,u2
	wt1 = wtc*u1
	wt2 = wtc-wt1
	mu1 = muc-sqrt(ssqc*wt1*wt2)*u2/wt1
	mu2 = (wtc*muc-wt1*mu1)/wt2

c.. check mu in order: else reject

	j1 = prev(j)
	if(j1.ne.0) then
c		if(qdebug) write(18,*) 'mu1,mu(j1)',mu1,mu(j1)
		if(mu(j1).gt.mu1) then
			nrejr = nrejr+1
			go to 24
		end if
	else if(qunif) then
		if(mu1.lt.xi-unhw) go to 24
	end if
	j1 = next(j)
	if(j1.ne.0) then
c		if(qdebug) write(18,*) 'mu2,mu(j1)',mu2,mu(j1)
		if(mu(j1).lt.mu2) then
			nrejr = nrejr+1
			go to 24
		end if
	else if(qunif) then
		if(mu2.gt.xi+unhw) go to 24
	end if

c	u3 = rbeta(ss,ss)
	call rbeta2(ss,ss,u3)
c	if(qdebug) write(18,*) 'isweep,u3',isweep,u3
	cu3 = 1.0-u3
	temp = wtc*ssqc*(1.0-u2**2)
	ssq1 = temp*u3/wt1
	ssq2 = temp*cu3/wt2

c.. allocate observations in component to be split

	lprob = 0.0
	st1 = 0
	st2 = 0
	l1 = 0
	l2 = 0
	i = start(j)
	do while(i.ne.0)
	if(qprior) then
	    p1 = wt1
	    p2 = wt2
	else
	    p1 = wt1*exp(max(-20.0,-0.5*(yv(i)-mu1)**2/ssq1))/sqrt(ssq1)
	    p2 = wt2*exp(max(-20.0,-0.5*(yv(i)-mu2)**2/ssq2))/sqrt(ssq2)
	end if
	    in = inext(i)
	    if(sdrand(u).lt.p1/(p1+p2)) then
		inext(i) = st1
		if(st1.eq.0) ilast = i
		st1 = i
		l1 = l1+1
		lprob = lprob+log(p1)-log(p1+p2)
	    else
		inext(i) = st2
		st2 = i
		l2 = l2+1
		lprob = lprob+log(p2)-log(p1+p2)
	    end if
	    i = in
	end do

c	if(qdebug) write(18,*) 'isweep,l1,l2',isweep,l1,l2

	klow = k

c.. compute ratio and decide

	logratio = 0.0
	if(.not.qprior) then
	i = st1
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq1/ssqc)
     &		-0.5*((yv(i)-mu1)**2/ssq1-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	i = st2
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq2/ssqc)
     &		-0.5*((yv(i)-mu2)**2/ssq2-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	end if

c	if(qdebug) then
c		write(18,*) 'split1:isw,move,logr',isweep,move,logratio
c	end if

c.. p(k,w,z) terms

	logratio = logratio+(lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wt1*wt2/wtc)-logbeta(delta,klow*delta)
     &		+l1*log(wt1/wtc)+l2*log(wt2/wtc)

c	if(qdebug) then
c		write(18,*) 'split2:isw,move,logr',isweep,move,logratio
c	end if

c.. p(theta) terms

	if(qunif) then
	    logratio = logratio+const1+alpha*log(beta)
     &		-log(2.0*unhw)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
	else
	    logratio = logratio+const1+alpha*log(beta)+0.5*log(kappa)
     &		-0.5*kappa*((mu1-xi)**2+(mu2-xi)**2-(muc-xi)**2)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
c	if(qdebug) then
c		write(18,*) const1,alpha,beta,kappa,ssq1,ssq2,ssqc
c		write(18,*) const1+alpha*log(beta)+0.5*log(kappa)
c		write(18,*) -0.5*kappa*((mu1-xi)**2+(mu2-xi)**2-(muc-xi)**2)
c		write(18,*) -(alpha+1.0)*log(ssq1*ssq2/ssqc)
c		write(18,*) -beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
c		write(18,*) +log(pf(klow))
c	end if
	    if(sp.gt.1) then
		sqrk = sqrt(kappa)
		jplus = next(j)
		jminus = prev(j)
		if(jminus.eq.0) then
		    temp = log(pnorm((mu1-xi)*sqrk)/
     &			pnorm((muc-xi)*sqrk))
		else
		    temp = log((pnorm((mu1-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk))
     &			/(pnorm((muc-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk)))
		end if
		temp = temp+log(pnorm((mu2-xi)*sqrk)-pnorm((mu1-xi)*sqrk))
		if(jplus.eq.0) then
		    temp = temp+log((1.0-pnorm((mu2-xi)*sqrk))/
     &			(1.0-pnorm((muc-xi)*sqrk)))
		else
		    temp = temp+log((pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((mu2-xi)*sqrk))
     &			/(pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((muc-xi)*sqrk)))
		end if

		logratio = logratio+(sp-1)*temp
	    end if
	end if

c	if(qdebug) then
c		write(18,*) 'split3:isw,move,logr',isweep,move,logratio
c	end if

c.. proposal terms


	logratio = logratio+const2
     &		+log(d(klow+1)/b(klow))-lprob
     &		-(ws-1.0)*log(u1*cu1)
     &		-(ms-1.0)*log(u2*(1.0-u2))
     &		-(ss-1.0)*log(u3*cu3)

c	if(qdebug) then
c		write(18,*) 'split4:isw,move,logr',isweep,move,logratio
c	end if

c.. Jacobian terms

	logratio = logratio
     &		+log(wtc*abs(mu1-mu2)*ssq1*ssq2/ssqc)
     &		-log(u2*(1.0-u2**2)*u3*cu3)

c	if(qdebug) then
c		write(18,*) 'split5:isw,move,logr',isweep,move,logratio
c	end if

	logratio = max(-20.0,min(20.0,logratio))

      if(sdrand(u).lt.exp(logratio)) then

c.. accept split

c		if(qdebug) then
c		    write(18,'("split accepted")')
c		end if
		naccs = naccs+1

		if(free.eq.0) return
		jnew = free
		free = next(jnew)

		jnext = next(j)
		next(jnew) = jnext
		if(jnext.ne.0) prev(jnext) = jnew
		next(j) = jnew
		prev(jnew) = j

		hwm = max(hwm,jnew)

		start(j) = st1
		start(jnew) = st2
		leng(j) = l1
		leng(jnew) = l2
		wt(j) = wt1
		wt(jnew) = wt2
		mu(j) = mu1
		mu(jnew) = mu2
		ssq(j) = ssq1
		ssq(jnew) = ssq2

		k = k+1
		if(l1*l2.eq.0) kemp = kemp+1
		kpos = k-kemp

c		if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

	else

	    if(l1.ne.0) then
		start(j) = st1
		inext(ilast) = st2
	    else
		start(j) = st2
	    end if

	end if

24	continue

	else if(move.eq.combine) then

c	   write(18,*)'combine'

c.. combine move ----------------------------------------------

	ntryc = ntryc+1


30	j1 = 1+int(hwm*sdrand(u))
	if(prev(j1).lt.0) go to 30
	j2 = next(j1)
	if(j2.eq.0) go to 30


	st1 = start(j1)
	st2 = start(j2)
	l1 = leng(j1)
	l2 = leng(j2)

	wt1 = wt(j1)
	wt2 = wt(j2)
	mu1 = mu(j1)
	mu2 = mu(j2)
	ssq1 = ssq(j1)
	ssq2 = ssq(j2)
	wtc = wt1+wt2
	muc = (wt1*mu1+wt2*mu2)/wtc
c	ssqc = (wt1*(mu1**2+ssq1)
c     &		+wt2*(mu2**2+ssq2))/wtc-muc**2
	ssqc = wt1*wt2*((mu1-mu2)/wtc)**2 +(wt1*ssq1+wt2*ssq2)/wtc
	u1 = wt1/wtc
	cu1 = wt2/wtc
c	u2 = (muc-mu1)*wt1/sqrt(ssqc*wt1*wt2)
	u2 = (mu2-mu1)*sqrt(wt1*wt2/ssqc)/wtc
	u2 = max(u2,1e-12)
	u2 = min(u2,1.0-1e-4)
	u3 = wt1*ssq1/(wt1*ssq1+wt2*ssq2)
	cu3 = wt2*ssq2/(wt1*ssq1+wt2*ssq2)

c	write(18,'(3f8.5)')wt1,wt2,wtc,mu1,mu2,muc,ssq1,ssq2,ssqc,u1,u2,u3
c	write(18,*)ssq1,ssq2,ssqc,u2

	lprob = 0.0
	i = st1
	do while(i.ne.0)
	if(qprior) then
	    p1 = wt1
	    p2 = wt2
	else
	    p1 = wt1*exp(max(-20.0,-0.5*(yv(i)-mu1)**2/ssq1))/sqrt(ssq1)
	    p2 = wt2*exp(max(-20.0,-0.5*(yv(i)-mu2)**2/ssq2))/sqrt(ssq2)
	end if
	    lprob = lprob+log(p1)-log(p1+p2)
	    ilast = i
	    i = inext(i)
	end do
	i = st2
	do while(i.ne.0)
	if(qprior) then
	    p1 = wt1
	    p2 = wt2
	else
	    p1 = wt1*exp(max(-20.0,-0.5*(yv(i)-mu1)**2/ssq1))/sqrt(ssq1)
	    p2 = wt2*exp(max(-20.0,-0.5*(yv(i)-mu2)**2/ssq2))/sqrt(ssq2)
	end if
	    lprob = lprob+log(p2)-log(p1+p2)
	    i = inext(i)
	end do

	klow = k-1

c.. compute ratio and decide

	logratio = 0.0
	if(.not.qprior) then
	i = st1
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq1/ssqc)
     &		-0.5*((yv(i)-mu1)**2/ssq1-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	i = st2
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq2/ssqc)
     &		-0.5*((yv(i)-mu2)**2/ssq2-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	end if

	loglikrat = logratio

c	if(qdebug) then
c		write(18,'(2i6,1pe16.7)') isweep,move,logratio
c	end if

c.. p(k,w,z) terms

	logratio = logratio+(lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wt1*wt2/wtc)-logbeta(delta,klow*delta)
     &		+l1*log(wt1/wtc)+l2*log(wt2/wtc)

c	if(qdebug) then
c		write(18,'(2i6,1pe16.7)') isweep,move,logratio
c	end if

c.. p(theta) terms

	if(qunif) then
	    logratio = logratio+const1+alpha*log(beta)
     &		-log(2.0*unhw)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
	else
	    logratio = logratio+const1+alpha*log(beta)+0.5*log(kappa)
     &		-0.5*kappa*((mu1-xi)**2+(mu2-xi)**2-(muc-xi)**2)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
	    if(sp.gt.1) then
		sqrk = sqrt(kappa)
		jplus = next(j2)
		jminus = prev(j1)
		if(jminus.eq.0) then
		    temp = log(pnorm((mu1-xi)*sqrk)/
     &			pnorm((muc-xi)*sqrk))
		else
		    temp = log((pnorm((mu1-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk))
     &			/(pnorm((muc-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk)))
		end if
		temp = temp+log(pnorm((mu2-xi)*sqrk)-pnorm((mu1-xi)*sqrk))
		if(jplus.eq.0) then
		    temp = temp+log((1.0-pnorm((mu2-xi)*sqrk))/
     &			(1.0-pnorm((muc-xi)*sqrk)))
		else
		    temp = temp+log((pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((mu2-xi)*sqrk))
     &			/(pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((muc-xi)*sqrk)))
		end if

		logratio = logratio+(sp-1)*temp
	    end if
	end if

c	if(qdebug) then
c		write(18,'(2i6,1pe16.7)') isweep,move,logratio
c	end if


c.. proposal terms

	logratio = logratio+const2
     &		+log(d(klow+1)/b(klow))-lprob
     &		-(ws-1.0)*log(u1*cu1)
     &		-(ms-1.0)*log(u2*(1.0-u2))
     &		-(ss-1.0)*log(u3*cu3)

c	if(qdebug) then
c		write(18,'(2i6,1pe16.7)') isweep,move,logratio
c	end if


c.. Jacobian terms

	logratio = logratio
     &		+log(wtc*abs(mu1-mu2)*ssq1*ssq2/ssqc)
     &		-log(u2*(1.0-u2**2)*u3*cu3)

c	if(qdebug) then
c		write(18,'(2i6,1pe16.7)') isweep,move,logratio
c	end if

	logratio = max(-20.0,min(20.0,logratio))


	if(sdrand(u).lt.exp(-logratio)) then

c.. accept combine

c		if(qdebug) then
c		    write(18,'("combine accepted")')
c		end if
	    naccc = naccc+1

c 	    next(j1) = next(j2)
c	    if(next(j1).ne.0) prev(next(j1)) = j1
	if(prev(j2).ne.0) then
		next(prev(j2)) = next(j2)
	else
		first = next(j2)
	end if
	if(next(j2).ne.0) prev(next(j2)) = prev(j2)
	    next(j2) = free
	    free = j2
	    prev(j2) = -1

	    leng(j1) = l1+l2
	    wt(j1) = wtc
	    mu(j1) = muc
	    ssq(j1) = ssqc
	    if(l1.ne.0) then
		inext(ilast) = st2
	    else
		start(j1) = st2
	    end if

	    k = k-1
	    if(l1*l2.eq.0) kemp = kemp-1
	    kpos = k-kemp

c		if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

	end if


	else if(move.eq.allocate) then

c	   write(18,*)'alloc'

c.. Gibbs move for allocation -------------------------------

	if(k.gt.1) then

	   call stdalloc(yv,n,wt,mu,ssq,ncmax,start,leng,next,pw,inext,
     &		first,delta,qprior)
	   kemp = 0
	   j = first
	   do while(j.ne.0)
	      if(leng(j).eq.0) kemp = kemp+1
	      j = next(j)
	   end do

c	   if(qdebug) then
c	      write(18,'("allocation move")')
c	      j = first
c	      do while(j.ne.0)
c		 write(18,'(2i6)') j,leng(j)
c		 j = next(j)
c	      end do
c	   end if
	   

	   
	end if
	
	else if(move.eq.weights) then

c	   write(18,*)'weights'

c.. Gibbs move for weights  ---------------------------------

c	if(qdebug) then
c		write(18,'("weights move")')
c	end if
	wtsum = 0.0
	j = first
	do while(j.ne.0)
c	    wt(j) = rgamma(delta+leng(j))
	    call rgamma2(delta+leng(j),wt(j))
	    wtsum = wtsum+wt(j)
	    j = next(j)
	end do
	j = first
	do while(j.ne.0)
	    wt(j) = wt(j)/wtsum
c	    if(qdebug) write(18,'(i4,f10.4)') j,wt(j)
	    j = next(j)
	end do


	else if(move.eq.parameters) then

c	   write(18,*)'params'

c.. Metropolis and/or Gibbs moves for component parameters --


	call gauss(mun,hwm)
	j = first
	do while(j.ne.0)
	if(qprior) then
	    if(qunif) then
		mun(j) = xi+unhw*(2.0*sdrand(u)-1.0)
	    else
		mun(j) = xi+mun(j)/sqrt(kappa)
	    end if
	else
	    ysum = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ysum = ysum+yv(i)
		i = inext(i)
	    end do
	    if(qunif) then
		if(leng(j).eq.0) then
		    mun(j) = xi+unhw*(2.0*sdrand(u)-1.0)
		else
		    mun(j) = ysum/leng(j)+sqrt(ssq(j)/leng(j))*mun(j)
		end if
	    else
		con = 1.0/(leng(j)/ssq(j)+kappa)
		mun(j) = con*(ysum/ssq(j)+xi*kappa)+sqrt(con)*mun(j)
	    end if
	end if

	j = next(j)
	end do

c.. check order first
	jp = first
	jq = next(jp)
	do while(jq.ne.0)
	   if(mun(jp).ge.mun(jq)) go to 66
	   jp = jq
	   jq = next(jq)
	end do
	
c.. for unif option, have to check in range
	if(qunif) then
	   if(mun(first).lt.xi-unhw.or.mun(jp).gt.xi+unhw)
     &			go to 66
	end if

c.. if sp>1, calc acceptance prob, and reject if necessary
	if(sp.gt.1) then
	   sqrk = sqrt(kappa)
	   jp = first
	   temp = log(pnorm((mun(jp)-xi)*sqrk)/pnorm((mu(jp)-xi)*sqrk))
	   jq = next(jp)
	   do while(jq.ne.0)
	      temp = temp+log((pnorm((mun(jq)-xi)*sqrk)
     &		   -pnorm((mun(jp)-xi)*sqrk))
     &			/(pnorm((mu(jq)-xi)*sqrk)
     &			-pnorm((mu(jp)-xi)*sqrk)))
	      jp = jq
	      jq = next(jq)
	   end do
	   temp = temp+log((1.0-pnorm((mun(jp)-xi)*sqrk))
     &			/(1.0-pnorm((mu(jp)-xi)*sqrk)))
	   logratio = (sp-1)*temp
	   logratio = max(-20.0,min(0.0,logratio))
	   if(sdrand(u).gt.exp(logratio)) go to 66
	end if

c.. copy over if accepted
	j = first
	do while(j.ne.0)
	   mu(j) = mun(j)
	   j = next(j)
	end do
	

 66	if(qdebug) then
c		write(18,'("means updated")')
c		j = first
c		do while(j.ne.0)
c		    write(18,'(i4,f10.4)') j,mu(j)
c		    j = next(j)
c		end do
	end if

	j = first
	do while(j.ne.0)
	if(qprior) then
c	    ssq(j) = beta/rgamma(alpha)
		call rgamma2(alpha,temp)
		ssq(j) = beta/temp
	else
	    ssd = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ssd = ssd+(yv(i)-mu(j))**2
		i = inext(i)
	    end do
c	    ssq(j) = (beta+0.5*ssd)/rgamma(alpha+0.5*leng(j))
		call rgamma2(alpha+0.5*leng(j),temp)
		ssq(j) = (beta+0.5*ssd)/temp
	end if

	j = next(j)
	end do

c	if(qdebug) then
c		write(18,'("variances updated")')
c		j = first
c		do while(j.ne.0)
c		    write(18,'(i4,f10.4)') j,ssq(j)
c		    j = next(j)
c		end do
c	end if

	else if(move.eq.hyper) then

c	   write(18,*)'hyper'

c.. Gibbs move for hyperparameters --------------------------

	if(qbeta) then
		j = first
		sum = 0.0
		do while(j.ne.0)
			sum = sum+1.0/ssq(j)
			j = next(j)
		end do
c		beta = rgamma(ggg+k*alpha)/(hhh+sum)
		call rgamma2(ggg+k*alpha,temp)
		beta = temp/(hhh+sum)
	end if
	if(qkappa) then
		j = first
		sum = 0.0
		do while(j.ne.0)
			sum = sum+(mu(j)-xi)**2
			j = next(j)
		end do
c		kappa = rgamma(eee+0.5*k)/(fff+0.5*sum)
		call rgamma2(eee+0.5*k,temp)
		kappa = temp/(fff+0.5*sum)
	end if


c.. Birth of an empty component

	else if(move.eq.birth) then

c	   write(18,*)'birth'

	    ntryb = ntryb+1

c.. compute logratio

	    klow = k
	    kemplow = kemp

c	    wtstar = rbeta(1.0,real(klow))
		call rbeta2(1.0,real(klow),wtstar)

c	    ssqstar = beta/rgamma(alpha)
		call rgamma2(alpha,temp)
		ssqstar = beta/temp

	    if(qunif) then
		mustar(1) = xi+unhw*(2.0*sdrand(u)-1.0)
	    else
	    call gauss(mustar,1)
	    mustar(1) = xi+mustar(1)/sqrt(kappa)
	    end if

c.. p(k,w,z) terms

	    logratio = (lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wtstar)-logbeta(delta,klow*delta)
     &		+(n+klow*(delta-1.0))*log(1.0-wtstar)

c.. p(theta) and proposal terms

	    logratio = logratio
     &		+log((klow+1)*d(klow+1)/(b(klow)*(kemplow+1)))
     &		-(klow-1)*log(1.0-wtstar)+logbeta(1.0,real(klow))

c.. Jacobian terms

	    logratio = logratio+(klow-1)*log(1.0-wtstar)

	    logratio = max(-20.0,min(20.0,logratio))

	    if(sdrand(u).lt.exp(logratio)) then

		scale = 1.0-wtstar
		j = first
		do while(j.ne.0)
		    wt(j) = scale*wt(j)
		    j = next(j)
		end do

		if(free.eq.0) return
		jnew = free
		free = next(jnew)

		jprev = 0
		j = first
		do while(j.ne.0)
			if(mu(j).gt.mustar(1)) go to 81
			jprev = j
			j = next(j)
		end do

81		next(jnew) = j
		prev(jnew) = jprev
		if(jprev.ne.0) then
			next(jprev) = jnew
		else
			first = jnew
		end if
		if(j.ne.0) prev(j) = jnew

		hwm = max(hwm,jnew)

		start(jnew) = 0
		leng(jnew) = 0
		mu(jnew) = mustar(1)
		ssq(jnew) = ssqstar
		wt(jnew) = wtstar

		k = k+1
		kemp = kemp+1
		kpos = k-kemp

c		if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

		naccb = naccb+1

	    end if

c.. Death of an empty component

	else if(move.eq.death) then

c	   write(18,*)'death'

	    ntryd = ntryd+1

	    if(kemp.gt.0) then
80		jkill = 1+int(hwm*sdrand(u))
		if(prev(jkill).lt.0) go to 80
		if(leng(jkill).ne.0) go to 80

c.. compute logratio
		wtstar = wt(jkill)
		mustar(1) = mu(jkill)
		ssqstar = ssq(jkill)
		klow = k-1
		kemplow = kemp-1

c.. p(k,w,z) terms

		logratio = (lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wtstar)-logbeta(delta,klow*delta)
     &		+(n+klow*(delta-1.0))*log(1.0-wtstar)

c.. p(theta) and proposal terms

		logratio = logratio
     &		+log((klow+1)*d(klow+1)/(b(klow)*(kemplow+1)))
     &		-(klow-1)*log(1.0-wtstar)+logbeta(1.0,real(klow))

c.. Jacobian terms

		logratio = logratio+(klow-1)*log(1.0-wtstar)

		logratio = max(-20.0,min(20.0,logratio))

		if(sdrand(u).lt.exp(-logratio)) then

			if(prev(jkill).ne.0) then
				next(prev(jkill)) = next(jkill)
			else
				first = next(jkill)
			end if
			if(next(jkill).ne.0)
     &			    prev(next(jkill)) = prev(jkill)
			next(jkill) = free
			free = jkill
			prev(jkill) = -1

			scale = 1.0/(1.0-wtstar)
			j = first
			do while(j.ne.0)
			    wt(j) = scale*wt(j)
			    j = next(j)
			end do

			k = k-1
			kemp = kemp-1
			kpos = k-kemp

c			if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

			naccd = naccd+1

		end if

	   end if


c.. Trap for illegal move

	else

		return

	end if



c.. end of 'istep' loop:
	end do

c.. logging stage

	kpos = k-kemp

c	if(qkreg.and..not.qfix) then
c	    if(mod(isweep,nskdel).eq.0) then
c		write(8,'(2i3,i8)') k,kpos,isweep
c	    end if
c	end if

	if(isweep.gt.0) then

	count(k) = count(k)+1
	countpos(kpos) = countpos(kpos)+1

c.. updates complete: record information

	if(qrkpos) then
		krep = kpos
	else
		krep = k
	end if
	if(krep.le.10) then
		j = first
		ij = 1
		do while(j.ne.0)
		    if((.not.qrkpos).or.leng(j).gt.0) then
			wtav(krep,ij) = wtav(krep,ij)+wt(j)
			muav(krep,ij) = muav(krep,ij)+mu(j)
			sigav(krep,ij) = sigav(krep,ij)+sqrt(ssq(j))
			if(qfull) avn(krep,ij) = avn(krep,ij)+leng(j)
			ij = ij+1
		    end if
		    j = next(j)
		end do
	end if

	if((.not.qprior).and.qfull) then

	    dev = 0.0
	    devc = 0.0
	    do i = 1,n
	       dvdy = 1.0
		dens = 0.0
		dens2 = 0.0
		j = first
		do while(j.ne.0)
		    temp = rr2pi*
     &			exp(max(-20.0,-0.5*(yv(i)-mu(j))**2/ssq(j)))/
     &			sqrt(ssq(j))
		    dens = dens+wt(j)*temp
		    dens2 = dens2+((leng(j)+delta)/(n+k*delta))*temp
		    j = next(j)
		end do
		dev = dev+log(dens*dvdy)
		devc = devc+log(dens2*dvdy)
	    end do
	    dev = -2.0*dev
	    devc = -2.0*devc
c	    devc = dev+devcorrection

	    avdev(krep) = avdev(krep)+dev
	    avdevc(krep) = avdevc(krep)+devc

	end if

	if(mod(isweep,nspace).eq.0) then

	ent = entropy(n,leng,ncmax,first,next)

c	if(qbeta.or.qkappa) write(10,*) beta,kappa,k,kpos

c	if(qfull) write(9,'(2i4,2f8.2)') k,kpos,dev,devc

	if(krep.le.10.and.qfull.and.(qpwms.or.qpalloc)) then
c	    if(.not.fileopen(krep)) then
c		istd2 = 6-int(log10(0.5+krep))
c		write(num2,'(i6)') krep
c		if(qpwms) then
c		open(25+krep,file=
c     &			prfx(1:npf)//".wms."//num2(istd2:6),
c     &                  status='unknown')
c		end if
c		if(qpalloc) then
c		open(25+ncmax+krep,file=
c     &			prfx(1:npf)//".z."//num2(istd2:6),
c     &                  status='unknown')
c		end if
c		fileopen(krep) = .true.
c	    end if
	    if(qpwms) then
	        j = first
		do while(j.ne.0)
c		    if((.not.qrkpos).or.leng(j).gt.0) then
c			write(25+krep,'(3f12.5)') wt(j),mu(j),sqrt(ssq(j))
c		    end if
		    j = next(j)
		end do
	    end if
	    if(qpalloc) then
		j = first
		ij = 1
		do while(j.ne.0)
			i = start(j)
			do while(i.ne.0)
				z(i) = ij
				i = inext(i)
			end do
			j = next(j)
			ij = ij+1
		end do
c		write(25+ncmax+krep,'(20i3)') (z(i),i=1,n)
	    end if
	end if

	end if

	nkemp = nkemp+1
	avkemp = avkemp+kemp
	if(kemp.gt.0) ppkemp = ppkemp+1.0
	if(kemp.le.9) nfkemp(kemp) = nfkemp(kemp)+1

	if(mod(isweep,nspace).eq.0) then
c		write(17,*) k,kpos,ent
		enttr(isweep/nspace) = ent
		ktr(isweep/nspace) = k
		off(isweep/nspace) = off0
		if(off0+k.le.nparsamp) then
			j = first
			do while(j.ne.0)
				off0 = off0+1
				partr(1,off0) = wt(j)
				partr(2,off0) = mu(j)
				partr(3,off0) = sqrt(ssq(j))
				j = next(j)
			end do
		end if
		if(qfull) devtr(isweep/nspace) = dev
	end if

	if(mod(isweep,max(1,nsweep/nsamp)).eq.0) then

c	write(4,'(2i3)') k,kpos
	j = first
	do while(j.ne.0)
	    stdev = sqrt(ssq(j))
c	    write(4,'(3f12.5,i4)') wt(j),mu(j),stdev,leng(j)
	    j = next(j)
	end do
	end if

c.. accumulate mean densities

	if(qfull) then
	    krec = min(k,ncd)
	    countden(krec) = countden(krec)+1
	    do iyd = 1,ngrid
	       yvg = yd0+iyd*ydinc
	       dvdy = 1.0
		pwsum = 0.0
		j = first
		ij = 1
		do while(j.ne.0)
		    pw(ij) = rr2pi*wt(j)*
     &			exp(max(-20.0,-0.5*(yvg-mu(j))**2/
     &			ssq(j)))/sqrt(ssq(j))
		if(pw(ij).gt.huge(1.2)) then
			call intpr('krec',4,krec,1)
			call intpr('isweep',6,isweep,1)
			call realpr('ssq(j)',6,ssq(j),1)
			pw(ij) = 0.0
		end if			
		    den(krec,iyd) = den(krec,iyd)+pw(ij)*dvdy
		    pwsum = pwsum+pw(ij)
		    j = next(j)
		    ij = ij+1
		end do
		if(qpclass) then
		    if(k.ge.k1.and.k.le.k2) then
		    ko = ((k+k1-3)*(k-k1))/2
		    do ij = 1,k-1
			pclass(ko+ij,iyd) = pclass(ko+ij,iyd)
     &				+pw(ij)/pwsum
		    end do
		    end if
		end if
	    end do

		if(qpclass) then
		if(k.ge.k1.and.k.le.k2) then
		    ko = ((k+k1-3)*(k-k1))/2
		    j = first
		    ij = 1
		    do while(ij.lt.k)
			i = start(j)
			do while(i.ne.0)
			    pz(ko+ij,i) = pz(ko+ij,i)+1
			    i = inext(i)
			end do
			j = next(j)
			ij = ij+1
		    end do
		end if
		end if
	end if

	if(qstop) return

	end if

c.. end of main loop

	if(mod(isweep,100).eq.0) then
		hwm = 0
		j = first
		do while(j.ne.0)
			hwm = max(hwm,j)
			j = next(j)
		end do
	end if

c	if(qdebug.and.mod(isweep,nsweep/50).eq.0) write(18,*) isweep

	end do

c	write(0,*)

c	write(7,'("splits:  ",i6," out of",i7)') naccs,ntrys
c	write(7,'("split rej r:  ",i6)') nrejr
c	write(7,'("combines:",i6," out of",i7)') naccc,ntryc
c	if(qempty) then
c		write(7,'("births:  ",i6," out of",i7)') naccb,ntryb
c		write(7,'("deaths:  ",i6," out of",i7)') naccd,ntryd
c	end if

	acctry(1,1) = naccs
	acctry(2,1) = ntrys
	acctry(1,2) = naccc
	acctry(2,2) = ntryc
	acctry(1,3) = naccb
	acctry(2,3) = ntryb
	acctry(1,4) = naccd
	acctry(2,4) = ntryd

c.. output posterior probs for k, and Bayes factors

	kbase = 0
	do k = 1,ncmax
		pw(k) = real(count(k))/nsweep
		if(kbase.eq.0.and.count(k).ne.0) kbase = k
	end do
c	write(12,'(5f16.6)') pw

	do k = 1,ncmax
	    bf(k) = (pw(k)/pw(kbase))/exp(lp(k)-lp(kbase))
	end do

c	write(12,'(1x)')
c	write(12,'(5f16.6)') bf

c.. write prior probs p(k)

	do k = 1,ncmax
	    pw(k) = exp(lp(k))
	end do
c	write(12,'(1x)')
c	write(12,'(5f16.6)') pw

c.. write posterior probs for number of nonempty components

	do k = 1,ncmax
		   pw(k) = real(countpos(k))/nsweep
	end do
c	write(12,'(1x)')
c	write(12,'(5f16.6)') pw

	avkemp = avkemp/nkemp
	ppkemp = ppkemp/nkemp
c	write(12,'(2f10.4)') avkemp,ppkemp
c	write(12,'(10i8)') nfkemp

	if(qfull) then

	    avdevall = 0.0
	    avdevcall = 0.0
	    do k = 1,ncmax
		avdevall = avdevall+avdev(k)
		avdevcall = avdevcall+avdevc(k)
	    end do
	    avdevall = avdevall/nsweep
	    avdevcall = avdevcall/nsweep

	    if(qrkpos) then
		do k = 1,ncmax
		    avdev(k) = avdev(k)/max(1,countpos(k))
		    avdevc(k) = avdevc(k)/max(1,countpos(k))
		end do
	    else
		do k = 1,ncmax
		    avdev(k) = avdev(k)/max(1,count(k))
		    avdevc(k) = avdevc(k)/max(1,count(k))
		end do
	    end if
c	    write(12,'(1x)')
c	    write(12,'(5f16.6)') avdev
c	    write(12,'(1x)')
c	    write(12,'(5f16.6)') avdevc
	end if

c.. output posterior expectations of parameters

	do k = 1,ncmax2
	if(qrkpos) then
		ckrep = countpos(k)
	else
		ckrep = count(k)
	end if
c	write(11,'(i3,i8)') k,ckrep
	ckrep = max(1,ckrep)
	do ij = 1,k
		wtav(k,ij) = wtav(k,ij)/ckrep
		muav(k,ij) = muav(k,ij)/ckrep
		sigav(k,ij) = sigav(k,ij)/ckrep
c		write(11,'(3f12.5)') wtav(k,ij),muav(k,ij),sigav(k,ij)
	end do
	end do

c.. output average group sizes

	if(qfull) then
	    do k = 1,ncmax2
		if(qrkpos) then
		    ckrep = countpos(k)
		else
		    ckrep = count(k)
		end if
c		write(16,'(2i8)') k,ckrep
		ckrep = max(1,ckrep)
		do ij = 1,k
		    avn(k,ij) = avn(k,ij)/ckrep
		end do
c		write(16,'(10f8.2)') (avn(k,ij),ij=1,k)
	    end do
	end if
c.. output densities

c	do iyd=1,6
c	call realpr('den(,iyd)',9,den(1,iyd),7)
c	end do

	if(qfull) then
	    do k = 1,ncd-1
	    countden(ncd) = countden(ncd)+countden(k)
	    do iyd = 1,ngrid
	    den(ncd,iyd) = den(ncd,iyd)+den(k,iyd)
	    end do
	    end do
c	    if(qfull)then
c		write(13,*) (yd0+iyd*ydinc,iyd=1,ngrid)
c	    end if
	    do k = 1,ncd
	    do iyd = 1,ngrid
	    den(k,iyd) = den(k,iyd)/max(countden(k),1)
	    end do
c	    if(qfull)then
c		write(13,*) (den(k,iyd),iyd=1,ngrid)
c	    end if
	    end do

c	    write(12,'("     D-hat:",g13.5)') avdevall
c	    write(12,'(" new D-hat:",g13.5)') avdevcall
	    do k = 1,ncd
		ggh(k) = 0.0
		if(countden(k).gt.0) then
		  do i = 1,n
		    iyd = max(1,min(ngrid,int(0.5+(y(i)-yd0)/ydinc)))
		    ggh(k) = ggh(k)+log(den(k,iyd))
		  end do
		  ggh(k) = -2.0*ggh(k)
		end if
	    end do
c	    write(12,'("  G(g-hat):",g13.5)') ggh(ncd)
c	    write(12,'("G(g-hat_k):",3g13.5/(11x,3g13.5))')
c     &		(ggh(k),k=1,ncd-1)

	end if

c.. output predictive classification

	if(qfull.and.qpclass) then

c	write(14,*) (yd0+iyd*ydinc,iyd=1,ngrid)
	do k = k1,k2
	    ko = ((k+k1-3)*(k-k1))/2
c	    write(14,'(i4)') k
	    do ij = 1,k-1
		do iyd = 1,ngrid
		    pclass(ko+ij,iyd) = pclass(ko+ij,iyd)/max(1,count(k))
		end do
		if(ij.ne.1) then
		    do iyd = 1,ngrid
			pclass(ko+ij,iyd) = pclass(ko+ij,iyd)+
     &				pclass(ko+ij-1,iyd)
		    end do
		end if
c	    write(14,'(10f8.4)') (pclass(ko+ij,iyd),iyd=1,ngrid)
	    end do
	end do

c.. output within-sample classification

	do k = k1,k2
	    ko = ((k+k1-3)*(k-k1))/2
c	    write(15,'(i4)') k
	    do ij = 1,k-1
		do i = 1,n
		    pz(ko+ij,i) = pz(ko+ij,i)/max(1,count(k))
		end do
		if(ij.ne.1) then
		    do i = 1,n
			pz(ko+ij,i) = pz(ko+ij,i)+pz(ko+ij-1,i)
		    end do
		end if
c		write(15,'(10f8.4)') (pz(ko+ij,i),i=1,n)
	    end do
	end do

	end if

	return
	end

