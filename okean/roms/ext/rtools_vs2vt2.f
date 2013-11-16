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
      real Aweight,Bweight,ds,sc_w(0:N),Cs_w(0:N),sc__w,
     &     theta_s,theta_b,Csur,Cbot,Cweight, sc_r(N),
     &     Cs_r(N),sc__r

Cf2py intent(out) sc_r
Cf2py intent(out) Cs_r
Cf2py intent(out) sc_w
Cf2py intent(out) Cs_w

        Aweight=1.0
        Bweight=1.0
        ds=1.0/N

        sc_w(N)=0.0
        Cs_w(N)=0.0
        DO k=N-1,1,-1
          sc__w=ds*(k-N)
          sc_w(k)=sc__w
          IF (theta_s.gt.0.0) THEN
            Csur=(1.0-COSH(theta_s*sc__w))/
     &           (COSH(theta_s)-1.0)
            IF (theta_b.gt.0.0) THEN
              Cbot=SINH(theta_b*(sc__w+1.0))/
     &             SINH(theta_b)-1.0
              Cweight=(sc__w+1.0)**Aweight*
     &                (1.0+(Aweight/Bweight)*
     &                        (1.0-(sc__w+1.0)**Bweight))
              Cs_w(k)=Cweight*Csur+(1.0-Cweight)*Cbot
            ELSE
              Cs_w(k)=Csur
            END IF
          ELSE
            Cs_w(k)=sc__w
          END IF
        END DO
        sc_w(0)=-1.0
        Cs_w(0)=-1.0

        DO k=1,N
          sc__r=ds*(k-N-0.5)
          sc_r(k)=sc__r
          IF (theta_s.gt.0.0) THEN
            Csur=(1.0-COSH(theta_s*sc__r))/
     &           (COSH(theta_s)-1.0)
            IF (theta_b.gt.0.0) THEN
              Cbot=SINH(theta_b*(sc__r+1.0))/
     &             SINH(theta_b)-1.0
              Cweight=(sc__r+1.0)**Aweight*
     &                (1.0+(Aweight/Bweight)*
     &                        (1.0-(sc__r+1.0)**Bweight))
              Cs_r(k)=Cweight*Csur+(1.0-Cweight)*Cbot
            ELSE
              Cs_r(k)=Csur
            END IF
          ELSE
            Cs_r(k)=sc__r
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

