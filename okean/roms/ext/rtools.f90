! f2py -c rtools.f -m rtools --fcompiler=gnu95
! f2py -c rtools.f -m rtools --f90flags==-32
! f2py -c rtools.f -m rtools --fcompiler=pg


      function interp_spl(z,v,lev,N)
      integer :: N
      real ( kind = 8 ) :: z(N), v(N),c(4,N), lev, interp_spl, ppvalu
      c(1,:)=v
      call cubspl(z,c,N,0,0)
      interp_spl=ppvalu(z,c(:,1:N-1),N-1,4,lev,0)
      return
      end


      function interp_lin(z,v,lev,N)
      ! Nos extremos extrapola linearmente usado o primeiro ou ultimo declives
      integer :: k, N
      real ( kind = 8 ) :: z(N), v(N), lev, interp_lin
      real (kind=8) :: a,b

      if (lev.le.z(1)) then
        k=2
      else if (lev.gt.z(N)) then
        k=N
      else
        do k=2,N
          if (lev.le.z(k)) then
            exit
          endif
        enddo
      endif

      a=z(k)-lev
      b=lev-z(k-1)
      interp_lin=(v(k)*b+v(k-1)*a)/(a+b)
      return
      end

      subroutine roms_slicez(s,v,h,zeta,tts,ttb,hc,Nr, &
      Vtransform,Vstretching,lev,SN,SPL, &
      N,Ny,Nx)

      integer :: N,Ny,Nx,j,i,Nr, Vtransform,Vstretching
      real ( kind=8 ) :: z_w(0:Nr,Ny,Nx),h(Ny,Nx),zeta(Ny,Nx), &
          z_r(Nr,Ny,Nx),hc,tts,ttb, v(N,Ny,Nx), s(Ny,Nx), &
          lev(Ny,Nx), interp_lin,interp_spl
      logical :: SN  ! masked near surface
      logical :: SPL ! use spline instead of linear interpolation

!f2py intent(out) s
!Cf2py integer optional, intent(in) :: n=shape(v,0)
!Cf2py integer optional, intent(in) :: m=shape(v,1)
!Cf2py integer optional, intent(in) :: l=shape(v,2)

      call s_levels(h,zeta,tts,ttb,hc,Nr,Ny,Nx,z_r,z_w, &
      Vtransform,Vstretching)

      if (N.eq.Nr) then

        do j=1,Nx
          do i=1,Ny
!            if (z_r(1,i,j).gt.lev .or. z_r(N,i,j).lt.lev) then
!              s(i,j)=-99.

            if (z_r(1,i,j).gt.lev(i,j)) then
              s(i,j)=-99.
            elseif (z_r(N,i,j).lt.lev(i,j)) then
              if (SN) then
                s(i,j)=-99.
              else
                s(i,j)=v(N,i,j)
              endif

            else
              if (SPL) then
                s(i,j)=interp_spl(z_r(:,i,j),v(:,i,j),lev(i,j),N)
              else
                s(i,j)=interp_lin(z_r(:,i,j),v(:,i,j),lev(i,j),N)
              endif

            endif
          enddo
        enddo

      else
        do j=1,Nx
          do i=1,Ny
!            if (z_w(0,i,j).gt.lev .or. z_w(N,i,j).lt.lev) then
!              s(i,j)=-99.

            if (z_w(0,i,j).gt.lev(i,j)) then
              s(i,j)=-99.
            elseif (z_w(Nr,i,j).lt.lev(i,j)) then
              if (SN) then
                s(i,j)=-99.
              else
                s(i,j)=v(N,i,j)
              endif

            else
              if (SPL) then
                s(i,j)=interp_spl(z_w(:,i,j),v(:,i,j),lev(i,j),N)
              else
                s(i,j)=interp_lin(z_w(:,i,j),v(:,i,j),lev(i,j),N)
              endif

            endif
          enddo
        enddo
      endif

      end


      subroutine generic_slicez(u,v,z,mask,lev,SPL,S,M,N)

      integer :: S,M,N,i,j,k
      real ( kind=8 ) :: z(S,M,N),v(S,M,N), mask(S,M,N), &
                         u(M,N),lev(M,N), interp_lin,interp_spl, &
                         kini(M,N),kend(M,N)
      logical :: SPL ! use spline instead of linear interpolation

!f2py intent(out) u


      ! 999 --> k0 (initial level not masked) higher than val
      ! 9999 --> all water column lower than val
      ! -999 --> all masked


      ! assume q z aumenta de 1 para S

      nDec=0
      do i=1,M
        do j=1,N
          ! assumir q existe mascara so num (ou ambos) os extremos
          ! achar 1o ultimos elementos nao mascarado
          k0=0
          k1=0
          do k=1,S
            if ((k0.eq.0).and.(mask(k,i,j).eq.1)) then
              k0=k
            else if ((k0.gt.0).and.(k1.eq.0).and.(mask(k,i,j).eq.0)) then
              k1=k-1
            endif
          enddo
          if ((k0.gt.0).and.(k1.eq.0)) then
            k1=S
          endif

          kini(i,j)=k0
          kend(i,j)=k1
        enddo
      enddo


      do i=1,M
        do j=1,N
          k0=kini(i,j)
          k1=kend(i,j)

          if ((k0.gt.0).and.(k1.gt.0)) then
            if (k0.eq.k1) then ! se so tem um avalor nao mascarado
              if (z(k0,i,j).eq.lev(i,j)) then
                u(i,j)=z(k0,i,j)
              else
                u(i,j)=999. ! como se primeiro valor ja fora do rango
              endif
            else
              if (z(k0,i,j).le.lev(i,j)) then
                if (z(k1,i,j).ge.lev(i,j)) then

                  if ((i.eq.101).and.(j.eq.101)) then
                    print*,k0,k1,z(k0:k1,i,j),lev(i,j)
                    print*,k0,k1,v(k0:k1,i,j),lev(i,j)
                  endif

                  if (SPL) then
                    u(i,j)=interp_spl(z(k0:k1,i,j),v(k0:k1,i,j),lev(i,j),k1-k0+1)
                  else
                    u(i,j)=interp_lin(z(k0:k1,i,j),v(k0:k1,i,j),lev(i,j),k1-k0+1)
                  endif

                else
                   ! primeiro valor ok, mas lev nao encontado no rango
                   u(i,j)=9999.
                endif
              else
                u(i,j)=999. ! primeiro valor ja esta fora do rango
              endif

            endif ! k0 eq k1
          else
             u(i,j)=-999. ! tudo mascarado
          endif
        enddo
      enddo

      end



      subroutine scoord(theta_s,theta_b,N,sc_r,Cs_r,sc_w,Cs_w, &
      Vstretching)
      implicit none
      integer :: N,k,Vstretching
!      integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-bit
!      not a good idea ... check:
!      https://sysbio.ioc.ee/projects/f2py2e/FAQ.html#q-what-if-fortran-90-code-uses-type-spec-kind-kind
!      thus:
      integer, parameter :: r8 = 8
      real (kind=r8) :: theta_s,theta_b,sc_r(N),Cs_r(N),sc_w(0:N), &
      Cs_w(0:N)
      real (kind=r8) :: cff1,cff2,ds,sc__w,sc__r
      real (kind=r8) :: Aweight,Bweight,Csur,Cbot,Cweight

!f2py intent(out) sc_r
!f2py intent(out) Cs_r
!f2py intent(out) sc_w
!f2py intent(out) Cs_w

!     -----------------------------------------------------------------------
!     Original vertical strectching function, Song and Haidvogel (1994).
!     -----------------------------------------------------------------------
      if (Vstretching.eq.1) then
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
            Cs_w(k)=(1.0_r8-theta_b)*cff1*SINH(theta_s*sc_w(k))+ &
                   theta_b*(cff2*TANH(theta_s*(sc_w(k)+0.5_r8))- &
                   0.5_r8)

            Cs_r(k)=(1.0_r8-theta_b)*cff1*SINH(theta_s*sc_r(k))+ &
                    theta_b*(cff2*TANH(theta_s*(sc_r(k)+0.5_r8))- &
                    0.5_r8)
          ELSE
            Cs_w(k)=sc_w(k)
            Cs_r(k)=sc_r(k)
          END IF
        END DO

!     -----------------------------------------------------------------------
!     A. Shchepetkin vertical stretching function. This function was
!     improved further to allow bottom refiment (see Vstretching=4).
!     -----------------------------------------------------------------------
      else if (Vstretching.eq.2) then
        Aweight=1.0_r8
        Bweight=1.0_r8
        ds=1.0_r8/REAL(N,r8)

        sc_w(N)=0.0_r8
        Cs_w(N)=0.0_r8
        DO k=N-1,1,-1
          sc__w=ds*real(k-N,r8)
          sc_w(k)=sc__w
          IF (theta_s.gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s*sc__w))/ &
                (COSH(theta_s)-1.0_r8)
            IF (theta_b.gt.0.0_r8) THEN
              Cbot=SINH(theta_b*(sc__w+1.0_r8))/ &
                  SINH(theta_b)-1.0_r8
              Cweight=(sc__w+1.0_r8)**Aweight* &
                     (1.0_r8+(Aweight/Bweight)* &
                             (1.0_r8-(sc__w+1.0_r8)**Bweight))
              Cs_w(k)=Cweight*Csur+(1.0_r8-Cweight)*Cbot
            ELSE
              Cs_w(k)=Csur
            END IF
          ELSE
            Cs_w(k)=sc__w
          END IF
        END DO
        sc_w(0)=-1.0_r8
        Cs_w(0)=-1.0_r8

        DO k=1,N
          sc__r=ds*(k-N-0.5_r8)
          sc_r(k)=sc__r
          IF (theta_s.gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s*sc__r))/ &
                (COSH(theta_s)-1.0_r8)
            IF (theta_b.gt.0.0_r8) THEN
              Cbot=SINH(theta_b*(sc__r+1.0_r8))/ &
                  SINH(theta_b)-1.0_r8
              Cweight=(sc__r+1.0_r8)**Aweight* &
                     (1.0_r8+(Aweight/Bweight)* &
                             (1.0_r8-(sc__r+1.0_r8)**Bweight))
              Cs_r(k)=Cweight*Csur+(1.0_r8-Cweight)*Cbot
            ELSE
              Cs_r(k)=Csur
            END IF
          ELSE
            Cs_r(k)=sc__r
          END IF
        END DO

!     -----------------------------------------------------------------------
!     R. Geyer stretching function for high bottom boundary layer
!     resolution.
!     -----------------------------------------------------------------------
      else if (Vstretching.eq.3) then
        print*, 'TODO ...'

!     -----------------------------------------------------------------------
!     A. Shchepetkin improved double vertical stretching functions with
!     bottom refiment.
!     -----------------------------------------------------------------------
      else if (Vstretching.eq.4) then
        ds=1.0_r8/real(N,r8)
!
        sc_w(N)=0.0_r8
        Cs_w(N)=0.0_r8
        DO k=N-1,1,-1
          sc__w=ds*(k-N)
          sc_w(k)=sc__w
          IF (theta_s.gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s*sc__w))/ &
                (COSH(theta_s)-1.0_r8)
          ELSE
            Csur=-sc__w**2
          END IF
          IF (theta_b.gt.0.0_r8) THEN
            Cbot=(EXP(theta_b*Csur)-1.0_r8)/ &
                (1.0_r8-EXP(-theta_b))
            Cs_w(k)=Cbot
          ELSE
            Cs_w(k)=Csur
          END IF
        END DO
        sc_w(0)=-1.0_r8
        Cs_w(0)=-1.0_r8
!
        DO k=1,N
          sc__r=ds*(k-N-0.5_r8)
          sc_r(k)=sc__r
          IF (theta_s.gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s*sc__r))/ &
                (COSH(theta_s)-1.0_r8)
          ELSE
            Csur=-sc__r**2
          END IF
          IF (theta_b.gt.0.0_r8) THEN
            Cbot=(EXP(theta_b*Csur)-1.0_r8)/ &
                (1.0_r8-EXP(-theta_b))
            Cs_r(k)=Cbot
          ELSE
            Cs_r(k)=Csur
          END IF
        END DO

      endif
      end



      subroutine s_levels(h,zeta,theta_s,theta_b,hc,N,Ny,Nx,z_r,z_w, &
      Vtransform,Vstretching)
      implicit none
      integer :: N,Ny,Nx, k,j,i, Vtransform,Vstretching
!      integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-bit
      integer, parameter :: r8 = 8
      real (kind=r8) :: h(Ny,Nx),zeta(Ny,Nx),hc, theta_s, theta_b
      real (kind=r8) :: z_r(N,Ny,Nx), z_w(0:N,Ny,Nx)

      real (kind=r8) :: cff_w,cff1_w,cff_r,cff1_r,cff2_w, &
          cff2_r, sc_w(0:N), Cs_w(0:N), sc_r(N),Cs_r(N), &
          hinv, hwater, z_w0, z_r0

!f2py intent(out) z_r
!f2py intent(out) z_w

      call scoord(theta_s,theta_b,N,sc_r,Cs_r,sc_w,Cs_w,Vstretching)

!     -----------------------------------------------------------------------
!     Original formulation: Compute vertical depths (meters, negative) at
!                           RHO- and W-points, and vertical grid
!     -----------------------------------------------------------------------
      if (Vtransform.eq.1) then

        DO j=1,Nx
          DO i=1,Ny
            z_w(0,i,j)=-h(i,j)
           END DO
           DO k=1,N
             cff_r=hc*(sc_r(k)-Cs_r(k))
             cff_w=hc*(sc_w(k)-Cs_w(k))
             cff1_r=Cs_r(k)
             cff1_w=Cs_w(k)
             DO i=1,Ny
              hwater=h(i,j)
              hinv=1.0_r8/hwater

              z_w0=cff_w+cff1_w*hwater
              z_w(k,i,j)=z_w0+zeta(i,j)*(1.0_r8+z_w0*hinv)
              z_r0=cff_r+cff1_r*hwater
              z_r(k,i,j)=z_r0+zeta(i,j)*(1.0_r8+z_r0*hinv)
            END DO
          END DO
        END DO


!     -----------------------------------------------------------------------
!     New formulation: Compute vertical depths (meters, negative) at
!                      RHO- and W-points, and vertical grid thicknesses.
!     -----------------------------------------------------------------------
      else if (Vtransform.eq.2) then

        DO j=1,Nx
          DO i=1,Ny
            z_w(0,i,j)=-h(i,j)
          END DO
          DO k=1,N
            cff_r=hc*sc_r(k)
            cff_w=hc*sc_w(k)
            cff1_r=Cs_r(k)
            cff1_w=Cs_w(k)
            DO i=1,Ny
              hwater=h(i,j)
              hinv=1.0_r8/(hc+hwater)
              cff2_r=(cff_r+cff1_r*hwater)*hinv
              cff2_w=(cff_w+cff1_w*hwater)*hinv

              z_w(k,i,j)=zeta(i,j)+(zeta(i,j)+hwater)*cff2_w
              z_r(k,i,j)=zeta(i,j)+(zeta(i,j)+hwater)*cff2_r
            END DO
          END DO
        END DO

      endif
      end


      subroutine  depthof(v,z,mask,val,S,M,N,depth)

      integer S, mask(M,N),dir, nInc, nDec
      real (kind=8) :: v(S,M,N),z(S,M,N),depth(M,N),val(M,N)
      real (kind=8) :: a,b

!f2py intent(out) depth
! Cf2py integer optional, intent(in) :: S=shape(v,0)
! Cf2py integer optional, intent(in) :: M=shape(v,1)
! Cf2py integer optional, intent(in) :: N=shape(v,2)

      ! 999 --> surface higher than val
      ! 9999 --> all whater column lower than val


      ! find if v increases or decreases with depth:
!      dir=0
!      do i=1,M
!        do j=1,N
!          if (mask(i,j).eq.1) then
!            if (v(1,i,j).gt.v(S,i,j)) then
!              dir=1
!            endif
!            exit
!          endif
!        enddo
!      enddo

! can't just test the 1st non masked element. It can be on mask
! with some extrapolated data, or it may not represent the whole
! dataset. So, find the number or points where v increases and
! decreases:

      dir=0
      nInc=0
      nDec=0
      do i=1,M
        do j=1,N
          if (mask(i,j).eq.1) then
            if (v(1,i,j).gt.v(S,i,j)) then
              nInc=nInc+1 ! aumenta da superfiie para o fundo
            else
              nDec=nDec+1
            endif
          endif
        enddo
      enddo
      if (nInc.gt.nDec) then
        dir=1
      else
        dir=0
      endif
      !print*, 'DIR=',dir, nInc, nDec

      if (dir.eq.1) then ! increases with depth (from k=S to k=1)
         do i=1,M
           do j=1,N
             depth(i,j)=999.
             if ((mask(i,j).eq.1).and.(v(S,i,j).le.val(i,j))) then
               do k=S-1,1,-1
                 if (v(k,i,j).gt.val(i,j)) then
                   a=v(k+1,i,j)-val(i,j)
                   b=val(i,j)-v(k,i,j)
                   depth(i,j)=(z(k,i,j)*a+z(k+1,i,j)*b)/(a+b)
                   exit
                 endif
               enddo
               if (depth(i,j).eq.999.) then
                 depth(i,j)=9999.
               endif
             endif
           enddo
         enddo

      else !  decrease with depth (from k=S to k=1)
         do i=1,M
           do j=1,N
             depth(i,j)=999.
             if ((mask(i,j).eq.1).and.(v(S,i,j).ge.val(i,j))) then
               do k=S-1,1,-1
                 if (v(k,i,j).lt.val(i,j)) then
                   a=v(k+1,i,j)-val(i,j)
                   b=val(i,j)-v(k,i,j)
                   depth(i,j)=(z(k,i,j)*a+z(k+1,i,j)*b)/(a+b)
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

      subroutine  depthof_mask3d(v,z,mask,val,S,M,N,depth)
      ! a busca comeca em k0=1 (ou proximo indice nao mascarado)
      ! ou seja, em geral, para oceano admite-se q k=1 seja a superficie

      integer :: S, mask(S,M,N),dir, nInc, nDec
      integer :: k0,k1,kini(M,N),kend(M,N)
      real (kind=8) :: v(S,M,N),z(S,M,N),depth(M,N),val(M,N)
      real (kind=8) :: a,b
      real (kind=8) :: interp_lin

!f2py intent(out) depth

      ! 999 --> surface higher than val
      ! 9999 --> all water column lower than val
      ! -999 --> all masked


      ! find if v increases or decreases with depth:

      dir=0
      nInc=0
      nDec=0
      do i=1,M
        do j=1,N
          ! assumir q existe mascara so num (ou ambos) os extremos
          ! achar 1o ultimos elementos nao mascarado
          k0=0
          k1=0
          do k=1,S
            if ((k0.eq.0).and.(mask(k,i,j).eq.1)) then
              k0=k
            else if ((k0.gt.0).and.(k1.eq.0).and.(mask(k,i,j).eq.0)) then
              k1=k-1
            endif
          enddo
          if ((k0.gt.0).and.(k1.eq.0)) then
            k1=S
          endif

          kini(i,j)=k0
          kend(i,j)=k1

          if ((k0.gt.0).and.(k1.gt.0).and.v(k0,i,j).gt.v(k1,i,j)) then
            nDec=nDec+1 ! diminiu de k0 para k1
          else
            nInc=nInc+1 ! aumenta de k0 para k1
          endif
        enddo
      enddo

      if (nInc.gt.nDec) then
        dir=1
      else
        dir=0
      endif
      !print*, 'DIR=',dir, nInc, nDec

      if (dir.eq.1) then ! increases de k0 para k1
         do i=1,M
           do j=1,N
             k0=kini(i,j)
             k1=kend(i,j)
             if ((k0.gt.0).and.(k1.gt.0)) then
               if (k0.eq.k1) then ! se so tem um avalor nao mascarado
                 if (v(k0,i,j).eq.val(i,j)) then
                   depth(i,j)=z(k0,i,j)
                 else
                   depth(i,j)=999. ! como se primeiro valor ja fora do rango
                 endif
               else
                 if (v(k0,i,j).le.val(i,j)) then
                   if (v(k1,i,j).ge.val(i,j)) then
                     !depth(i,j)=interp_lin(v(k0:k1,i,j),z(k0:k1,i,j),val(i,j),k1-k0+1)
                     ! interp_lin e interp_spl pode dar maus resultados
                     ! (incluindo extrapolacoes) pois v[k0:j1,i,j] nao cresce
                     ! monotonicamente. Ou seja, vou procurar nivel a nivel:
                     do k=k0+1,k1
                       if (v(k,i,j).gt.val(i,j)) then
                         a=v(k-1,i,j)-val(i,j)
                         b=val(i,j)-v(k,i,j)
                         depth(i,j)=(z(k,i,j)*a+z(k-1,i,j)*b)/(a+b)
                         exit
                       endif
                     enddo
                   else
                     ! primeiro valor ok, mas val nao encontado no rango
                     depth(i,j)=9999.
                   endif
                 else
                   depth(i,j)=999. ! primeiro valor ja esta fora do rango
                 endif
               endif ! k0 eq k1
             else
               depth(i,j)=-999. ! tudo mascarado
             endif

           enddo
         enddo

      else !  decrease de k0 para k1
         do i=1,M
           do j=1,N
             k0=kini(i,j)
             k1=kend(i,j)
             if ((k0.gt.0).and.(k1.gt.0)) then
               if (k0.eq.k1) then ! se so tem um avalor nao mascarado
                 if (v(k0,i,j).eq.val(i,j)) then
                   depth(i,j)=z(k0,i,j)
                 else
                   depth(i,j)=999. ! como se primeiro valor ja fora do rango
                 endif
               else
                 if (v(k0,i,j).ge.val(i,j)) then
                   if (v(k1,i,j).le.val(i,j)) then
                     do k=k0+1,k1
                       if (v(k,i,j).lt.val(i,j)) then
                         a=v(k-1,i,j)-val(i,j)
                         b=val(i,j)-v(k,i,j)
                         depth(i,j)=(z(k,i,j)*a+z(k-1,i,j)*b)/(a+b)
                         exit
                       endif
                     enddo
                   else
                     ! primeiro valor ok, mas val nao encontrado no rango
                     depth(i,j)=9999.
                   endif
                 else
                   depth(i,j)=999. ! primeiro valor ja esta fora do rango
                 endif
               endif ! k0 eq k1
             else
               depth(i,j)=-999. ! tudo mascarado
             endif ! k0>0 e k1>0

           enddo
         enddo


      endif

      end
