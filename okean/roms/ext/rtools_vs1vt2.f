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

      call s_levels(h,Zeta,hc,tts,ttb,Nr,M,L,z_r,z_w)

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


      subroutine scoord(theta_s,theta_b,N,sc_r,Cs_r,sc_w,Cs_w)
      integer N,k
      real ds,sc_w(0:N),Cs_w(0:N),
     &     theta_s,theta_b, sc_r(N),
     &     Cs_r(N)
      integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-bit

Cf2py intent(out) sc_r
Cf2py intent(out) Cs_r
Cf2py intent(out) sc_w
Cf2py intent(out) Cs_w

        IF (theta_s.ne.0.0_r8) THEN
          cff1=1.0_r8/SINH(theta_s)
          cff2=0.5_r8/TANH(0.5_r8*theta_s)
        END IF
        sc_w(0)=-1.0_r8
        Cs_w(0)=-1.0_r8
        ds=1.0_r8/REAL(N,r8)
        DO k=1,N
          sc_w(k)=ds*REAL(k-N,r8)
          sc_r(k)=ds*(REAL(k-N,r8)-0.5_r8)
          IF (theta_s.ne.0.0_r8) THEN
            Cs_w(k)=(1.0_r8-theta_b)*cff1*SINH(theta_s*sc_w(k))+
     &              theta_b*(cff2*TANH(theta_s*(sc_w(k)+0.5_r8))-
     &              0.5_r8)

            Cs_r(k)=(1.0_r8-theta_b)*cff1*SINH(theta_s*sc_r(k))+
     &               theta_b*(cff2*TANH(theta_s*(sc_r(k)+0.5_r8))-
     &               0.5_r8)
          ELSE
            Cs_w(k)=sc_w(k)
            Cs_r(k)=sc_r(k)
          END IF
        END DO
      end


      subroutine s_levels(h,zeta,hc,theta_s,theta_b,N,M,L,z_r,z_w)
      integer N,M,L, k,j,i
      real z_w(0:N,M,L),h(M,L),cff_w,cff1_w,cff_r,cff1_r,cff2_w,
     &     cff2_r, sc_w(0:N), Cs_w(0:N), sc_r(N),Cs_r(N),
     &     zeta(M,L), hinv, hc, theta_s, theta_b, z_r(N,M,L), hwater

Cf2py intent(out) z_r
Cf2py intent(out) z_w

      call scoord(theta_s,theta_b,N,sc_r,Cs_r,sc_w,Cs_w)


        DO j=1,M
          DO i=1,L
            z_w(0,j,i)=-h(j,i)
          END DO
          DO k=1,N
            cff_r=hc*sc_r(k)
            cff_w=hc*sc_w(k)
            cff1_r=Cs_r(k)
            cff1_w=Cs_w(k)
            DO i=1,L
              hwater=h(j,i)
              hinv=1.0/(hc+hwater)
              cff2_r=(cff_r+cff1_r*hwater)*hinv
              cff2_w=(cff_w+cff1_w*hwater)*hinv

              z_w(k,j,i)=zeta(j,i)+(zeta(j,i)+hwater)*cff2_w
              z_r(k,j,i)=zeta(j,i)+(zeta(j,i)+hwater)*cff2_r
            END DO
          END DO
        END DO
      end

