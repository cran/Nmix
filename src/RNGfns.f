	function sdrand(u)
	real u
	double precision unifrnd

	call rndstart()
	u = unifrnd()
	sdrand = u
	call rndend()
	return
	end 

	subroutine gauss(z,n)

c..    generates an n-vector of i.i.d. N(0,1) r.v.'s z

	real z(n)
	double precision normrnd

	call rndstart()
	do i = 1,n
		z(i) = normrnd()
	end do
	call rndend()

	return
	end
