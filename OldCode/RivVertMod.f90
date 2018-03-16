	MODULE RivVertMod
	USE RivVarVertMod
	USE CalcCond
	IMPLICIT NONE
	CONTAINS
	
	SUBROUTINE ZOTOCDTWO()
	IMPLICIT NONE
	INTEGER :: i,j
	REAL :: tmp1, tmp2
	REAL :: vck = 0.4
	    do i = 1,ns
	        do j = 1,nn
	            tmp1 = log(hl(i,j)/znaught(i,j)) - 1
!	            tmp2 = 
	            cd(i,j) = ((1/vck)*tmp1)**-2
	            if(cd(i,j) < cdmin)then  
		            cd(i,j) = cdmin
		        endif
		        if(cd(i,j) > cdmax) then
		            cd(i,j) = cdmax
		        endif
		        totcd(i,j) = cd(i,j)

	        enddo
	    enddo
    END SUBROUTINE 

	
	
	SUBROUTINE Z0TOCD()
	INTEGER :: i, j, k
!	cdmin = 0.001
!	cdmax = 0.1
		call alloc_vert()
         pi=acos(-1.)
        vkc=0.4
        beta=6.24
        do 10 i=1,ns
        do 10 j=1,nn
		if(hl(i,j) > 0) then
	        zeta0(i,j)= znaught(i,j)/hl(i,j)
		else
			zeta0(i,j)= znaught(i,j)/.01
		endif
        sczeta=(log(1./zeta0(i,j)))/(nz-1)
        do 4 k=1,nz
        zeta(k)=zeta0(i,j)*exp((k-1)*sczeta)
        if(zeta(k).lt.0.2) then
         ecoef(i,j,k)=vkc*zeta(k)*(1.-zeta(k))
        else
         ecoef(i,j,k)=vkc/beta
        endif
        if(k.eq.1) then
         f1(i,j,k)=0.
         fone(i,j)=0.
         f3(i,j,k)=0.
        else
         dzeta(k)=zeta(k)-zeta(k-1)
         finc1=(1.-zeta(k))/ecoef(i,j,k)
         finc2=(1.-zeta(k-1))/ecoef(i,j,k-1)
         finc=(finc1+finc2)
         dzp=0.5*dzeta(k)
         f1(i,j,k)=f1(i,j,k-1)+finc*dzp
         fone(i,j)=fone(i,j)+(f1(i,j,k)+f1(i,j,k-1))*dzp
         f3(i,j,k)=f3(i,j,k-1)+(f1(i,j,k)**2.+f1(i,j,k-1)**2.)*dzp
        endif
4       continue
        
		 cd(i,j) = 1/(fone(i,j)**2.)
		 if(cd(i,j) > cdmax) then 
		    cd(i,j) = cdmax
		 endif
		 if(cd(i,j) < cdmin) then
		    cd(i,j) = cdmin
		 endif
		 if(ibc(i,j).eq.0) then
		    cd(i,j) = 0
		 endif
		 totcd(i,j) = cd(i,j)

!		 write(6,*) i, j, fone(i,j), cd(i,j)
10		continue

	call dealloc_vert()


	END SUBROUTINE

    SUBROUTINE vert()
    IMPLICIT NONE
	INTEGER :: i, j, k, ismoo
	REAL :: tmp
		call alloc_vert()
         pi=acos(-1.)
        vkc=0.4
        beta=6.24
        do 10 i=1,ns
        do 10 j=1,nn
			if(roughnesstype == 1) then
				if(hl(i,j) > 0) then
					zeta0(i,j)= znaught(i,j)/hl(i,j)
				else	
					zeta0(i,j)= .001
				endif
					
			else
				if(hl(i,j) > 0) then 
!					zeta0(i,j) = exp(-1.*(vkc/sqrt(cd(i,j)))-1.)
!					zeta0(i,j) = (1./(exp(vkc/sqrt(totcd(i,j)))))/hl(i,j)
					zeta0(i,j) = 1.-(1./exp(vkc*totcd(i,j)))
				else
					zeta0(i,j) = 0.001
				endif
			endif
        sczeta=(log(1./zeta0(i,j)))/(nz-1)
        
        do 4 k=1,nz
        zeta(k)=zeta0(i,j)*exp((k-1)*sczeta)
!		if(ibc(i,j) .ne. -1) then
!			zz(i, j, k) = eta(i,j)-100. + (100.*zeta(k))
!		else 
			zz(i, j, k) = eta(i,j) + (hl(i,j)*zeta(k))
!		endif
!		if(ibc(i,j) == 0) then
!			zz(i, j, k) = eta(i,j)-100. + (100.*zeta(k))
!		else
!			zz(i, j, k) = eta(i,j) + (hl(i,j)*zeta(k))
!		endif					
        if(zeta(k).lt.0.2) then
         ecoef(i,j,k)=vkc*zeta(k)*(1.-zeta(k))
        else
         ecoef(i,j,k)=vkc/beta
        endif
        if(k.eq.1) then
         f1(i,j,k)=0.
         fone(i,j)=0.
         f3(i,j,k)=0.
        else
         dzeta(k)=zeta(k)-zeta(k-1)
         finc1=(1.-zeta(k))/ecoef(i,j,k)
         finc2=(1.-zeta(k-1))/ecoef(i,j,k-1)
         finc=(finc1+finc2)
         dzp=0.5*dzeta(k)
         f1(i,j,k)=f1(i,j,k-1)+finc*dzp
         fone(i,j)=fone(i,j)+(f1(i,j,k)+f1(i,j,k-1))*dzp
         f3(i,j,k)=f3(i,j,k-1)+(f1(i,j,k)**2.+f1(i,j,k-1)**2.)*dzp
        endif
4       continue
        do 6 k=1,nz
        f3(i,j,k)=f3(i,j,nz)-f3(i,j,k)
        if(k.eq.1) then
         f2(i,j,k)=0.
         ftwo(i,j)=0.
        else
         finc1=f3(i,j,k)/ecoef(i,j,k)
         finc2=f3(i,j,k-1)/ecoef(i,j,k-1)
         f2(i,j,k)=f2(i,j,k-1)+(finc1+finc2)*0.5*dzeta(k)
         ftwo(i,j)=ftwo(i,j)+(f2(i,j,k)+f2(i,j,k-1))*0.5*dzeta(k)
        endif 
6       continue       
        if(taus(i,j).eq.0.and.taun(i,j).eq.0) then
         theta(i,j)=1000.
        else
         theta(i,j)=atan2(taun(i,j),taus(i,j))
        endif
        ustr2(i,j)=(taus(i,j)**2.+taun(i,j)**2.)**0.5
10      continue
        do 50 i=1,ns
        do 50 j=1,nn
        if(i.eq.1.or.theta(i-1,j).eq.1000.) then
         dth=theta(i+1,j)-theta(i,j)
         scn=1.
        elseif(i.eq.ns.or.theta(i+1,j).eq.1000.) then
         dth=theta(i,j)-theta(i-1,j)
         scn=1.
        else
         dth=(theta(i+1,j)-theta(i-1,j))
         scn=2.
        endif
        if(dth.gt.pi) then
         dth=dth-2.*pi
        elseif(dth.lt.(-1.*pi)) then
         dth=2.*pi+dth
        endif
        dthds=dth/(scn*ds*rn(i,j))
        if(j.eq.1.or.theta(i,j-1).eq.1000.) then
         dth=theta(i,j+1)-theta(i,j)
         scn=1.
        elseif(j.eq.nn.or.theta(i,j+1).eq.1000.) then
         dth=theta(i,j)-theta(i,j-1)
         scn=1.
        else
         dth=(theta(i,j+1)-theta(i,j-1))
         scn=2.
        endif
        if(dth.gt.pi) then
         dth=dth-2.*pi
        elseif(dth.lt.(-1.*pi)) then
         dth=2.*pi+dth
        endif
        dthdn=dth/(scn*dn)
        if(ustr2(i,j).eq.0.or.ibc(i,j).ne.-1) then 
         rs(i,j)=10**8.
        else
!		Adjusted for Kootenai 7/22/03 jmn & rmcd
!		Center line and stream line curvature accounted for
        If(CALCQUASI3DRS) Then
         curvs=(dthds+(1./r(i)))*taus(i,j)/ustr2(i,j)
         curvn=dthdn*taun(i,j)/ustr2(i,j)
        ELSE
         curvs=((1./r(i)))*taus(i,j)/ustr2(i,j)
		 curvn = 0
        ENDIF
        
!         curvs=(dthds+(1./r(i)))*taus(i,j)/ustr2(i,j)
!         curvn=dthdn*taun(i,j)/ustr2(i,j)
         curv=curvs+curvn
         if(curv.eq.0.) then
          rs(i,j)=10**8.
         else 
          rs(i,j)=1./curv
         endif
        endif
        if(abs(rs(i,j)).lt.MinRS) then
         rs(i,j)=MinRS*sign(1.,rs(i,j))
        endif
        taux=ustr2(i,j)*hl(i,j)/(rn(i,j)*rs(i,j)) 
        vs=((ftwo(i,j)/fone(i,j))-f3(i,j,1))
        taux=taux*vs
        do 40 k=1,nz
        uu1=(ustr2(i,j)**.5)*cos(theta(i,j))*f1(i,j,k)
        vv1=(ustr2(i,j)**.5)*sin(theta(i,j))*f1(i,j,k)
        g1=(ustr2(i,j)**.5)*hl(i,j)/(rn(i,j)*rs(i,j))
        g2=((ftwo(i,j)/fone(i,j))*f1(i,j,k))-f2(i,j,k)
        uu2=-1.*g1*g2*sin(theta(i,j))
        vv2=g1*g2*cos(theta(i,j)) 
        uz(i,j,k)=uu1+uu2
40      vz(i,j,k)=vv1+vv2
        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
        taun(i,j)=taun(i,j)+taux*cos(theta(i,j))        
!        if(i.eq.1) then
!         dth=theta(i+1,j)-theta(i,j)
!         scn=1.
!        elseif(i.ne.1.and.theta(i-1,j).eq.1000.) then
!         dth = theta(i+1,j)-theta(i,j)
!         scn=1.
!        elseif(i.eq.ns) then
!         dth=theta(i,j)-theta(i-1,j)
!         scn=1.
!        elseif(i.ne.ns.and.theta(i+1,j).eq.1000.) then
!         dth=theta(i,j)-theta(i-1,j)
!         scn=1.
!        else
!!         dth=(theta(i+1,j)-theta(i-1,j))
!!         scn=2.
!         dth=(theta(i+1,j)-theta(i,j))
!         scn=1.
!       endif
!        if(dth.gt.pi) then
!         dth=dth-2.*pi
!        elseif(dth.lt.(-1.*pi)) then
!         dth=2.*pi+dth
!        endif
!        dthds=dth/(scn*ds*rn(i,j))
!        
!       
!        if(j.eq.1) then
!        dth=theta(i,j+1)-theta(i,j)
!        scn=1.
!        elseif(j.eq.nn) then
!        dth=theta(i,j)-theta(i,j-1)
!        scn=1.
!        else if(j.ne.1.and.theta(i,j-1).eq.1000.) then
!        dth = theta(i,j+1)-theta(i,j)
!        scn=1.
!        elseif(j.ne.nn.and.theta(i,j+1).eq.1000.) then
!        dth=theta(i,j)-theta(i,j-1)
!        scn=1.
!        else
!!        dth=(theta(i,j+1)-theta(i,j-1))
!!        scn=2.
!        dth=(theta(i,j+1)-theta(i,j))
!        scn=1.
!        endif
!        if(dth.gt.pi) then
!         dth=dth-2.*pi
!        elseif(dth.lt.(-1.*pi)) then
!         dth=2.*pi+dth
!        endif
!        dthdn=dth/(scn*dn)
!        
!        if(ustr2(i,j).eq.0.or.ibc(i,j).ne.-1) then 
!         rs(i,j)=10**8.
!        else
!!		Adjusted for Kootenai 7/22/03 jmn & rmcd
!!		Center line and stream line curvature accounted for
!        If(CALCQUASI3DRS) Then
!         curvs=(dthds+(1./r(i)))*taus(i,j)/ustr2(i,j)
!         curvn=dthdn*taun(i,j)/ustr2(i,j)
!        ELSE
!         curvs=((1./r(i)))*taus(i,j)/ustr2(i,j)
!		 curvn = 0
!        ENDIF
!         curv=curvs+curvn
!        if(curv.eq.0) then
!          rs(i,j)=10**8.
!        else 
!          rs(i,j)=1./curv
!         endif
!        endif
!        if(abs(rs(i,j)).lt.MinRS) then
!         rs(i,j)=MinRS*sign(1.,rs(i,j))
!        endif
!        taux=ustr2(i,j)*hl(i,j)/(rn(i,j)*rs(i,j)) 
!        vs=((ftwo(i,j)/fone(i,j))-f3(i,j,1))
!        taux=taux*vs
!        do 40 k=1,nz
!        uu1=(ustr2(i,j)**.5)*cos(theta(i,j))*f1(i,j,k)
!        vv1=(ustr2(i,j)**.5)*sin(theta(i,j))*f1(i,j,k)
!        g1=(ustr2(i,j)**.5)*hl(i,j)/(rn(i,j)*rs(i,j))
!        g2=((ftwo(i,j)/fone(i,j))*f1(i,j,k))-f2(i,j,k)
!        uu2=-1.*g1*g2*sin(theta(i,j))
!        vv2=g1*g2*cos(theta(i,j)) 
!        tmp = uu1+uu2
!        uz(i,j,k)=uu1+uu2
!40      vz(i,j,k)=vv1+vv2
!!		uz(i,j,k) = (ustr2(i,j)**.5)*f1(i,j,k)
!!40		vz(i,j,k) = g1*g2
!        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
!        taun(i,j)=taun(i,j)+taux*cos(theta(i,j))
50      continue
      IF(RSSMOO.gt.0) THEN
      DO  ISMOO=1,RSSMOO
          DO  I=2,ns-1
              DO  J=1,nn
              if(j.eq.1) then
               DUM1a(I,1)=(rs(I-1,1)+rs(I,1)+rs(I+1,1)+rs(I,j+1))/4.
              elseif(j.eq.nn) then
               DUM1a(I,25)=(rs(I-1,nn)+rs(I,nn)+rs(I+1,nn)+rs(i,nn-1))/4.
              else
               DUM1a(I,J)=rs(I,J)+SEDSMOOWGHT*rs(I+1,J)+SEDSMOOWGHT*rs(I-1,J)
               DUM1a(I,J)=DUM1a(I,J)+SEDSMOOWGHT*rs(I,J-1)+SEDSMOOWGHT*rs(I,J+1)
               DUM1a(I,J)=DUM1a(I,J)/5.
              endif
              ENDDO
          ENDDO
      DO 670 I=2,ns-1
      DO 670 J=1,nn
670    rs(I,J)=DUM1a(I,J)
      ENDDO


        do 750 i=1,ns
        do 750 j=1,nn
        if(abs(rs(i,j)).lt.MinRS) then
         rs(i,j)=MinRS*sign(1.,rs(i,j))
        endif
        taux=ustr2(i,j)*hl(i,j)/(rn(i,j)*rs(i,j)) 
        vs=((ftwo(i,j)/fone(i,j))-f3(i,j,1))
        taux=taux*vs
        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
750     taun(i,j)=taun(i,j)+taux*cos(theta(i,j))

     ENDIF

!! Changes below made 9/25/06.  Program was crashing in calculation
!! of dth. Previous version allowed array bounds to be exceeded.
!        do 50 i=1,ns
!        do 50 j=1,nn
!        	if(i.ne.1.and.i.ne.ns) then
!        		if(theta(i-1,j).eq.1000.) then
!				 dth=theta(i+1,j)-theta(i,j)
!				 scn=1.
!				elseif(theta(i+1,j).eq.1000.) then
!				 dth=theta(i,j)-theta(i-1,j)
!				 scn=1.
!				else
!				 dth=(theta(i+1,j)-theta(i-1,j))
!				 scn=2.
!				endif
!			elseif(i.eq.1) then
!				 dth=theta(i+1,j)-theta(i,j)
!				 scn=1.
!			else
!				 dth=theta(i,j)-theta(i-1,j)
!				 scn=1.
!			endif			
!
!        if(dth.gt.pi) then
!         dth=dth-2.*pi
!        elseif(dth.lt.(-1.*pi)) then
!         dth=2.*pi+dth
!        endif
!        dthds=dth/(scn*ds*rn(i,j))
!        	if(j.ne.1.and.j.ne.nn) then
!        		if(theta(i,j-1).eq.1000.) then
!					dth=theta(i,j+1)-theta(i,j)
!					scn=1.
!				elseif(theta(i,j+1).eq.1000.) then
!					dth=theta(i,j)-theta(i,j-1)
!					scn=1.
!				else
!					dth=(theta(i,j+1)-theta(i,j-1))
!					scn=2.
!				endif
!			elseif(j.eq.1) then
!				dth=theta(i,j+1)-theta(i,j)
!				scn=1.
!			else        	
!				dth=theta(i,j)-theta(i,j-1)
!				scn=1.
!			endif        	
!        if(dth.gt.pi) then
!         dth=dth-2.*pi
!        elseif(dth.lt.(-1.*pi)) then
!         dth=2.*pi+dth
!        endif
!        dthdn=dth/(scn*dn)
!        if(ustr2(i,j).eq.0.or.ibc(i,j).ne.-1) then 
!         rs(i,j)=10**8.
!        else
!!		Adjusted for Kootenai 7/22/03 jmn & rmcd
!!		Center line and stream line curvature accounted for
!        If(CALCQUASI3DRS == 1) Then
!         curvs=(dthds+(1./r(i)))*taus(i,j)/ustr2(i,j)
!         curvn=dthdn*taun(i,j)/ustr2(i,j)
!        ELSE
!         curvs=((1./r(i)))*taus(i,j)/ustr2(i,j)
!		 curvn = 0
!        ENDIF
!         curv=curvs+curvn
!        if(curv.eq.0) then
!          rs(i,j)=10**8.
!        else 
!          rs(i,j)=1./curv
!         endif
!        endif
!        if(abs(rs(i,j)).lt.MinRS) then
!         rs(i,j)=MinRS*sign(1.,rs(i,j))
!        endif
!        taux=ustr2(i,j)*hl(i,j)/(rn(i,j)*rs(i,j)) 
!        vs=((ftwo(i,j)/fone(i,j))-f3(i,j,1))
!        taux=taux*vs
!        do 40 k=1,nz
!        uu1=(ustr2(i,j)**.5)*cos(theta(i,j))*f1(i,j,k)
!        vv1=(ustr2(i,j)**.5)*sin(theta(i,j))*f1(i,j,k)
!        g1=(ustr2(i,j)**.5)*hl(i,j)/(rn(i,j)*rs(i,j))
!        g2=((ftwo(i,j)/fone(i,j))*f1(i,j,k))-f2(i,j,k)
!        uu2=-1.*g1*g2*sin(theta(i,j))
!        vv2=g1*g2*cos(theta(i,j)) 
!        tmp = uu1+uu2
!        uz(i,j,k)=uu1+uu2
!40      vz(i,j,k)=vv1+vv2
!!		uz(i,j,k) = (ustr2(i,j)**.5)*f1(i,j,k)
!!40		vz(i,j,k) = g1*g2
!        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
!        taun(i,j)=taun(i,j)+taux*cos(theta(i,j))
!50      continue

!      DO 670 ISMOO=1,1
!      DO 660 I=2,ns-1
!      DO 660 J=1,nn
!      if(j.eq.1) then
!       DUM1a(I,1)=(rs(I-1,1)+rs(I,1)+rs(I+1,1)+rs(I,j+1))/4.
!      elseif(j.eq.nn) then
!       DUM1a(I,25)=(rs(I-1,nn)+rs(I,nn)+rs(I+1,nn)+rs(i,nn-1))/4.
!      else
!       DUM1a(I,J)=rs(I,J)+rs(I+1,J)+rs(I-1,J)
!       DUM1a(I,J)=DUM1a(I,J)+rs(I,J-1)+rs(I,J+1)
!       DUM1a(I,J)=DUM1a(I,J)/5.
!      endif
!660   continue
!      DO 670 I=2,ns-1
!      DO 670 J=1,nn
!670    rs(I,J)=DUM1a(I,J)
!        do 750 i=1,ns
!        do 750 j=1,nn
!        if(abs(rs(i,j)).lt.MinRS) then
!         rs(i,j)=3000.*sign(1.,rs(i,j))
!        endif
!        taux=ustr2(i,j)*hl(i,j)/(rn(i,j)*rs(i,j)) 
!        vs=((ftwo(i,j)/fone(i,j))-f3(i,j,1))
!        taux=taux*vs
!        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
!750     taun(i,j)=taun(i,j)+taux*cos(theta(i,j))
!	call dealloc_vert()
!      return
      END SUBROUTINE vert

	  END MODULE RivVertMod
