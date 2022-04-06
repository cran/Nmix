      subroutine rgamma2(s,res)
      data e/2.71828182/

c.. Returns one value from the gamma distribution
c   GAMMA(s,1), i.e. with pdf = y^(s-1) exp(-y) / GAMMA(s)

      if(s.lt.1) then
            b = (s+e)/e
            c1 = 1.0/s
            do
                  bu = b*sdrand(u)
                  if(bu.le.1.0) then
                        res = exp(max(-30.0,c1*log(bu)))
                        if(sdrand(u).lt.exp(-res)) exit
                  else
                        res = -log((b-bu)/s)
                        if(sdrand(u).lt.res**(s-1)) exit
                  end if
            end do

      else if(s.eq.1) then
            u = sdrand(u)
            res = -log(u)

      else
            c1 = s-1.0
            c2 = (s-1.0/(6.0*s))/c1
            c3 = 2.0/c1
            c4 = c3+2.0
            c5 = 1.0/sqrt(s)
            do
                  u1 = sdrand(u)
                  u2 = sdrand(u)
                  if(s.gt.2.5) u1 = u2+c5*(1.0-1.86*u1)
                  if(u1.gt.0.0.and.u1.lt.1.0) then
                        w = c2*u2/u1
                        if(c3*u1+w+1.0/w.le.c4) exit
                        if(c3*log(u1)-log(w)+w.lt.1.0) exit
                  end if
            end do
            res = c1*w
      end if
      
      return

      end
