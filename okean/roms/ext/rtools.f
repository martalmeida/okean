c f2py -c rtools.f -m rtools --fcompiler=gnu95
c f2py -c rtools.f -m rtools --f90flags==-32
c f2py -c rtools.f -m rtools --fcompiler=pg



      subroutine roms_slicez(s,v,h,Zeta,hc,tts,ttb,lev,SN,Nr,N,M,L)

      integer N,M,L, k,j,i,Nr
      real z_w(0:Nr,M,L),h(M,L),Zeta(M,L),z_r(Nr,M,L),hc,tts,
     &     ttb, v(N,M,L), s(M,L), lev,a,b
      logical SN ! nans near syrface

Cf2py intent(out) s
Cf2py integer optional, intent(in) :: n=shape(v,0)
Cf2py integer optional, intent(in) :: m=shape(v,1)
Cf2py integer optional, intent(in) :: l=shape(v,2)

      call s_levels(h,Zeta,z_r,z_w,hc,tts,ttb,Nr,M,L)

      if (N.eq.Nr) then

        do i=1,L
          do j=1,M
!            if (z_r(1,j,i).gt.lev .or. z_r(N,j,i).lt.lev) then
!              s(j,i)=-99.

            if (z_r(1,j,i).gt.lev) then
              s(j,i)=-99.
            elseif (z_r(N,j,i).lt.lev) then
              if (SN) then
                s(j,i)=-99.
              else
                s(j,i)=v(N,j,i)
              endif

            else
              k=2
              do while (z_r(k,j,i).lt.lev)
                k=k+1
              enddo
              a=z_r(k,j,i)-lev
              b=lev-z_r(k-1,j,i)
              s(j,i)=(v(k,j,i)*b+v(k-1,j,i)*a)/(a+b)
            endif
          enddo
        enddo

      else
        do i=1,L
          do j=1,M
!            if (z_w(0,j,i).gt.lev .or. z_w(N,j,i).lt.lev) then
!              s(j,i)=-99.

            if (z_w(0,j,i).gt.lev) then
              s(j,i)=-99.
            elseif (z_w(N,j,i).lt.lev) then
              if (SN) then
                s(j,i)=-99.
              else
                s(j,i)=v(N,j,i)
              endif

            else
              k=1
              do while (z_w(k,j,i).lt.lev)
                k=k+1
              enddo
              a=z_w(k,j,i)-lev
              b=lev-z_w(k-1,j,i)
              s(j,i)=(v(k,j,i)*b+v(k-1,j,i)*a)/(a+b)
            endif
          enddo
        enddo
      endif

      end


      subroutine s_levels(h,Zeta,z_r,z_w,hc,theta_s,theta_b,N,M,L)
      integer N,M,L, k,j,i
      real z_w(0:N,M,L),h(M,L),cff_w,cff1_w,cff2_w,cff_r,cff1_r,cff2_r,
     &     sc_w(0:N), Cs_w(0:N), sc_r(N),Cs_r(N), z_w0, z_r0,
     &     Zeta(M,L), hinv(M,L),
     &     hc, cff1, cff2, theta_s, theta_b,
     &     z_r(N,M,L)

Cf2py intent(out) z_r
Cf2py intent(out) z_w

      cff1=1./sinh(theta_s)
      cff2=0.5/tanh(0.5*theta_s)

      sc_w(0)=-1.0
      Cs_w(0)=-1.0

      cff=1./float(N)
      do k=1,N,+1
        sc_w(k)=cff*float(k-N)
        Cs_w(k)=(1.-theta_b)*cff1*sinh(theta_s*sc_w(k))
     &             +theta_b*(cff2*tanh(theta_s*(sc_w(k)+0.5))-0.5)

        sc_r(k)=cff*(float(k-N)-0.5)
        Cs_r(k)=(1.-theta_b)*cff1*sinh(theta_s*sc_r(k))
     &             +theta_b*(cff2*tanh(theta_s*(sc_r(k)+0.5))-0.5)
      enddo

      do j=1,M
        do i=1,L
          z_w(0,j,i)=-h(j,i)
        enddo
        do k=1,N
          cff_w=hc*(sc_w(k)-Cs_w(k))
          cff1_w=Cs_w(k)
          cff2_w=sc_w(k)+1.

          cff_r=hc*(sc_r(k)-Cs_r(k))
          cff1_r=Cs_r(k)
          cff2_r=sc_r(k)+1.

          do i=1,L
            hinv(j,i)=1./h(j,i)

            z_w0=cff_w+cff1_w*h(j,i)
            z_w(k,j,i)=z_w0+Zeta(j,i)*(1.+z_w0*hinv(j,i))

            z_r0=cff_r+cff1_r*h(j,i)
            z_r(k,j,i)=z_r0+Zeta(j,i)*(1.+z_r0*hinv(j,i))
          enddo
        enddo
      enddo

      end



      subroutine  depthof(v,z,mask,val,S,M,N,depth)

      integer S, mask(M,N),dir
      real v(S,M,N),z(S,M,N),depth(M,N),val

Cf2py intent(out) depth
! Cf2py integer optional, intent(in) :: S=shape(v,0)
! Cf2py integer optional, intent(in) :: M=shape(v,1)
! Cf2py integer optional, intent(in) :: N=shape(v,2)

      ! 999 --> surface higher than val
      ! 9999 --> all whater column lower than val


      ! find if v increases or decreases with depth:
      dir=0
      do i=1,M
        do j=1,N
          if (mask(i,j).eq.1) then
            if (v(1,i,j).gt.v(S,i,j)) then
              dir=1
            endif
            exit
          endif
        enddo
      enddo

      if (dir.eq.1) then ! increases:
         do i=1,M
           do j=1,N
             depth(i,j)=999.
             if ((mask(i,j).eq.1).and.(v(S,i,j).le.val)) then
               do k=S,2,-1
                 if (v(k,i,j).gt.val) then
                   a=v(k,i,j)-val
                   b=val-v(k+1,i,j)
                   depth(i,j)=(z(k,i,j)*b+z(k+1,i,j)*a)/(a+b)
                   exit
                 endif
               enddo
               if (depth(i,j).eq.999.) then
                 depth(i,j)=9999.
               endif
             endif
           enddo
         enddo

      else !  decrease
         do i=1,M
           do j=1,N
             depth(i,j)=999.
             if ((mask(i,j).eq.1).and.(v(S,i,j).ge.val)) then
               do k=S,2,-1
                 if (v(k,i,j).lt.val) then
                   a=v(k,i,j)-val
                   b=val-v(k+1,i,j)
                   depth(i,j)=(z(k,i,j)*b+z(k+1,i,j)*a)/(a+b)
                   exit
                 endif
               enddo
               if (depth(i,j).eq.999.) then
                 depth(i,j)=9999.
               endif
             endif
           enddo
         enddo


      endif

      end
