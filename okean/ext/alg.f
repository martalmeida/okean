      subroutine meetpoint(x1,y1,x2,y2,resu,n1,n2,Nu)
       implicit none
       integer :: n1, n2, i, j, N,Nu !,Nf,k,ii, Nex
       integer, parameter :: r8=8
!       integer :: mask(0:(n1-1)*(n2-1)+1)
       real (kind=r8) :: x1(n1), x2(n2), y1(n1), y2(n2), epsx, epsy
       real (kind=r8) :: Xa,Ya,xa_,ya_,Xb,Yb,xb_,yb_,dx1,dy1,dx2,dy2
       real (kind=r8) :: x,y,m1,m2, mx,my,mx1,my1,mx2,my2
       real (kind=r8) :: res(2,n1*n2), resu(2,n1*n2)
! , resex(2,n1*n2)
!       real (kind=r8) :: resf(2,n1*n2)
       logical cond1,cond2, found

! max size of res is n1*n2 and will only happen in n1=n2=2 and the two
! segments are equal!

Cf2py intent(out) resu,Nu
Cf2py intent(in) :: n1=shape(x1,0)
Cf2py intent(in) :: n2=shape(x2,0)

       ! find a scaled epsilon:
       ! mx=0.25*(maxval(x1)+minval(x1)+maxval(x2)+minval(x2))
       ! my=0.25*(maxval(y1)+minval(y1)+maxval(y2)+minval(y2))
       ! use max to avoid zero:
       mx1=max(maxval(x1),maxval(-x1))
       mx2=max(maxval(x2),maxval(-x2))
       mx=max(mx1,mx2)

       my1=max(maxval(y1),maxval(-y1))
       my2=max(maxval(y2),maxval(-y2))
       my=max(my1,my2)

       epsx=mx*epsilon(mx)*10
       epsy=my*epsilon(my)*10

       N=0
       do i=1,n1-1
         Xa=max(x1(i),x1(i+1))
         Ya=max(y1(i),y1(i+1))
         xa_=min(x1(i),x1(i+1))
         ya_=min(y1(i),y1(i+1))
         dx1=x1(i+1)-x1(i)
         dy1=y1(i+1)-y1(i)
         do j=1,n2-1

           if ( ((x1(i).eq.x2(j)).and.(y1(i).eq.(y2(j)) )).or.
     &          ((x1(i).eq.x2(j+1)).and.(y1(i).eq.(y2(j+1)) )) ) then
             N=N+1
             res(1,N)=x1(i)
             res(2,N)=y1(i)
           elseif ( ((x1(i+1).eq.x2(j)).and.(y1(i+1).eq.(y2(j)) )).or.
     &         ((x1(i+1).eq.x2(j+1)).and.(y1(i+1).eq.(y2(j+1)) )) ) then
             N=N+1
             res(1,N)=x1(i+1)
             res(2,N)=y1(i+1)

           else
             Xb=max(x2(j),x2(j+1))
             Yb=max(y2(j),y2(j+1))
             xb_=min(x2(j),x2(j+1))
             yb_=min(y2(j),y2(j+1))
             dx2=x2(j+1)-x2(j)
             dy2=y2(j+1)-y2(j)

             if ((Xa.ge.xb_).and.(xa_.le.Xb).and.
     &           (Ya.ge.yb_).and.(ya_.le.Yb)) then

               cond1=.true.
               if ((dx1.ne.0.).and.(dx2.ne.0)) then
                 m1 = dy1/dx1
                 m2 = dy2/dx2
!!!                 print*, i,j, m1, m2, m1.eq.m2, m1.eq.0.,m2.eq.0.
                 if (m1.eq.m2) then
                   cond1=.false.
                 else
                    y = 1./(m1-m2) *
     &                  ( m1*y2(j) - m2*y1(i) + m1*m2*(x1(i)-x2(j)))

                    ! if some slope is zero, fix y (prev calc has some
                    ! errors
                    if (m1.eq.0.) then
                      y=y1(i)
                    elseif (m2.eq.0.) then
                      y=y2(j)
                    endif

                    if (m1.ne.0.) then
                      x = (y-y1(i)+m1*x1(i))/m1
                    elseif (m2.ne.0.) then
                      x = (y-y2(j)+m2*x2(j))/m2
                    endif
                 endif

               elseif ((dx1.eq.0.).and.(dx2.ne.0.)) then
                 m2 = dy2/dx2
                 x=x1(i)
                 y = m2*(x-x2(j)) + y2(j)
               elseif ((dx1.ne.0.).and.(dx2.eq.0.)) then
                  m1 = dy1/dx1
                  x = x2(j)
                  y = m1*(x-x1(i)) + y1(i)
               else
                 cond1=.false.
               endif


          !             cond2=(x.ge.max(xa_,xb_)).and.(x.le.min(Xa,Xb)).and.
          !     &             (y.ge.max(ya_,yb_)).and.(y.le.min(Ya,Yb))
          !
          !
          !             print*, x.ge.max(xa_,xb_), x.le.min(Xa,Xb),
          !     &               y.ge.max(ya_,yb_), y.le.min(Ya,Yb)
          !

               cond2=(x.ge.(max(xa_,xb_)-epsx)).and.
     &               (x.le.(min(Xa,Xb)+epsx)).and.
     &               (y.ge.(max(ya_,yb_)-epsy)).and.
     &               (y.le.(min(Ya,Yb)+epsy))


               if (cond1.and.cond2) then
                 N=N+1
                 res(1,N)=x
                 res(2,N)=y
               endif

             endif
           endif

         enddo
       enddo

       ! eliminate repeated points -------------------:
       Nu=0
       if (N.ge.1) then
         resu(1,1)=res(1,1)
         resu(2,1)=res(2,1)
         Nu=1

         do j=2,N
           do i=1,j-1
             found=.false.
             if ( (res(1,j).eq.res(1,i)).and.
     &                  (res(2,j).eq.res(2,i)) ) then
               found=.true.
               exit
             endif
           enddo

          if (.not.found) then
            Nu=Nu+1
            resu(1,Nu)=res(1,j)
            resu(2,Nu)=res(2,j)
          endif
         enddo
       endif


      end


      subroutine pdist(x0,x1,y0,y1,x,y,d)
      ! point (x,y) to line ([x0,y0]->[x1,y1]) distance
      real x0,x1,y0,y1,x,y,d

      d=((y0-y1)*x+(x1-x0)*y+(x0*y1-x1*y0))
     &           /sqrt((x1-x0)**2+(y1-y0)**2)
      if (d.lt.0.) then
        d=-d
      endif
      end


      subroutine bilin(x,y,v,xi,yi,vi,mask,maski,L,M,Li,Mi)
      !2d interpolation
      ! x,y need not to be fully regular!

      integer L,M,Li,Mi, i,j, ii,jj,maski(Li,Mi)
      real x(L,M), y(L,M), v(L,M), xi(Li,Mi), yi(Li,Mi), vi(Li,Mi)
!     real ax,bx,ay,by
      real xp(5), yp(5)
      logical in_,mask(L,M)
      real dr,dl,dt,db,r1,r2

Cf2py intent(out) vi,maski

      do i=1,Li
        do j=1,Mi
          maski(i,j)=0 ! outside domain x,y
        enddo
      enddo

      do ii=1,Li
        do jj=1,Mi
          in_=.false.
          i=0
          do while ((.not.in_).and.(i.lt.L-1))
            i=i+1
            j=0
            do while ((.not.in_).and.(j.lt.M-1))
              j=j+1

              xp(1)=x(i,j)
              xp(2)=x(i+1,j)
              xp(3)=x(i+1,j+1)
              xp(4)=x(i,j+1)
              xp(5)=x(i,j)

              yp(1)=y(i,j)
              yp(2)=y(i+1,j)
              yp(3)=y(i+1,j+1)
              yp(4)=y(i,j+1)
              yp(5)=y(i,j)

              xp1=min(xp(1),xp(2),xp(3),xp(4),xp(5))
              yp1=min(yp(1),yp(2),yp(3),yp(4),yp(5))
              xp2=max(xp(1),xp(2),xp(3),xp(4),xp(5))
              yp2=max(yp(1),yp(2),yp(3),yp(4),yp(5))

              if (xi(ii,jj).ge.xp1 .and. xi(ii,jj).le.xp2
     &          .and. yi(ii,jj).ge.yp1 .and. yi(ii,jj).le.yp2) then
!                call in_polygon(xp,yp,xi(ii,jj),yi(ii,jj),5,in_)
                call PNPOLY(xi(ii,jj),yi(ii,jj),xp,yp,5,inn)
                if (inn.eq.-1) then
                  in_=.false.
                else
                  in_=.true.
                endif
              endif

            enddo
          enddo
          if (in_) then
!            ax=xi(ii,jj)-x(i,j)
!            ay=yi(ii,jj)-y(i,j)
!            bx=x(i+1,j+1)-xi(ii,jj)
!            by=y(i+1,j+1)-yi(ii,jj)
!            vi(ii,jj)=(bx*by*v(i,j)+
!     &                 ax*by*v(i,j+1)+
!     &                 ax*ay*v(i+1,j+1)+
!     &                 bx*ay*v(i+1,j))/(bx*by + ax*by +ax*ay +bx*ay)
!

            ! distance to right edge:
            call pdist(xp(2),xp(3),yp(2),yp(3),xi(ii,jj),yi(ii,jj),dr)
            ! distance to left edge:
            call pdist(xp(1),xp(4),yp(1),yp(4),xi(ii,jj),yi(ii,jj),dl)
            ! distance to bottom edge:
            call pdist(xp(1),xp(2),yp(1),yp(2),xi(ii,jj),yi(ii,jj),db)
            ! distance to top edge:
            call pdist(xp(3),xp(4),yp(3),yp(4),xi(ii,jj),yi(ii,jj),dt)

            if (mask(i,j).or.mask(i+1,j).or.mask(i+1,j+1)
     &         .or.mask(i,j+1)) then
              maski(ii,jj)=-1 ! x,y domain is also masked
            else
              r1=(dr*v(i,j)+dl*v(i+1,j))/(dr+dl)
              r2=(dr*v(i,j+1)+dl*v(i+1,j+1))/(dr+dl)
              vi(ii,jj)=(dt*r1+db*r2)/(dt+db)
              maski(ii,jj)=1 ! not masked
            endif

          endif

        enddo
      enddo

      end


      subroutine bilin_fast(v,vi,maski,coefs,inds,L,M,Li,Mi)
      integer L,M,Li,Mi, i,j,ii,jj,inds(0:2,Li,Mi),maski(Li,Mi)
      real v(L,M), vi(Li,Mi), coefs(4,Li,Mi)
      real r1,r2,dr,dl,db,dt

Cf2py intent(out) vi,maski
      do ii=1,Li
        do jj=1,Mi
          maski(ii,jj)=inds(0,ii,jj)
          if (inds(0,ii,jj).eq.1) then

            dr=coefs(1,ii,jj)
            dl=coefs(2,ii,jj)
            db=coefs(3,ii,jj)
            dt=coefs(4,ii,jj)

            i=inds(1,ii,jj)
            j=inds(2,ii,jj)

            r1=(dr*v(i,j)+dl*v(i+1,j))/(dr+dl)
            r2=(dr*v(i,j+1)+dl*v(i+1,j+1))/(dr+dl)
            vi(ii,jj)=(dt*r1+db*r2)/(dt+db)
          endif
        enddo
      enddo
      end


      subroutine bilin_coefs(x,y,xi,yi,inds,coefs,mask,L,M,Li,Mi)

      integer L,M,Li,Mi, i,j, ii,jj,k,inds(0:2,Li,Mi)
      real x(L,M), y(L,M), xi(Li,Mi), yi(Li,Mi)
      real xp(5), yp(5),coefs(4,Li,Mi)
      real dr,dl,dt,db
      logical in_,mask(L,M)

Cf2py intent(out) coefs,inds

      do i=1,Li
        do j=1,Mi
          do k=1,4
            coefs(k,i,j)=0.
          enddo
          inds(0,i,j)=0
        enddo
      enddo

      do ii=1,Li
        do jj=1,Mi
          in_=.false.
          i=0
          do while ((.not.in_).and.(i.lt.L-1))
            i=i+1
            j=0
            do while ((.not.in_).and.(j.lt.M-1))
              j=j+1

              xp(1)=x(i,j)
              xp(2)=x(i+1,j)
              xp(3)=x(i+1,j+1)
              xp(4)=x(i,j+1)
              xp(5)=x(i,j)

              yp(1)=y(i,j)
              yp(2)=y(i+1,j)
              yp(3)=y(i+1,j+1)
              yp(4)=y(i,j+1)
              yp(5)=y(i,j)

              xp1=min(xp(1),xp(2),xp(3),xp(4),xp(5))
              yp1=min(yp(1),yp(2),yp(3),yp(4),yp(5))
              xp2=max(xp(1),xp(2),xp(3),xp(4),xp(5))
              yp2=max(yp(1),yp(2),yp(3),yp(4),yp(5))

              if (xi(ii,jj).ge.xp1 .and. xi(ii,jj).le.xp2
     &          .and. yi(ii,jj).ge.yp1 .and. yi(ii,jj).le.yp2) then
!                call in_polygon(xp,yp,xi(ii,jj),yi(ii,jj),5,in_)
                call PNPOLY(xi(ii,jj),yi(ii,jj),xp,yp,5,inn)
                if (inn.eq.-1) then
                  in_=.false.
                else
                  in_=.true.
                endif
              endif

            enddo
          enddo
          if (in_) then
            ! distance to right edge:
            call pdist(xp(2),xp(3),yp(2),yp(3),xi(ii,jj),yi(ii,jj),dr)
            ! distance to left edge:
            call pdist(xp(1),xp(4),yp(1),yp(4),xi(ii,jj),yi(ii,jj),dl)
            ! distance to bottom edge:
            call pdist(xp(1),xp(2),yp(1),yp(2),xi(ii,jj),yi(ii,jj),db)
            ! distance to top edge:
            call pdist(xp(3),xp(4),yp(3),yp(4),xi(ii,jj),yi(ii,jj),dt)

            coefs(1,ii,jj)=dr
            coefs(2,ii,jj)=dl
            coefs(3,ii,jj)=db
            coefs(4,ii,jj)=dt

            inds(0,ii,jj)=1
            inds(1,ii,jj)=i
            inds(2,ii,jj)=j

            if (mask(i,j).or.mask(i+1,j).or.mask(i+1,j+1)
     &         .or.mask(i,j+1)) then
              inds(0,ii,jj)=-1
            endif

          endif

        enddo
      enddo

      end




! PNPOLY - Point Inclusion in Polygon Test
! W. Randolph Franklin (WRF)
! http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
C
C     ..................................................................
C
C        SUBROUTINE PNPOLY
C
C        PURPOSE
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON
C
C        USAGE
C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )
C
C        DESCRIPTION OF THE PARAMETERS
C           PX      - X-COORDINATE OF POINT IN QUESTION.
C           PY      - Y-COORDINATE OF POINT IN QUESTION.
C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF
C                     VERTICES OF POLYGON.
C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF
C                     VERTICES OF POLYGON.
C           N       - NUMBER OF VERTICES IN THE POLYGON.
C           INOUT   - THE SIGNAL RETURNED:
C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
C                      1 IF THE POINT IS INSIDE OF THE POLYGON.
C
C        REMARKS
C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.
C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY
C           OPTIONALLY BE INCREASED BY 1.
C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING
C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX
C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING
C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.
C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.
C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM
C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT
C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE
C           POINT IS INSIDE OF THE POLYGON.
C
C     ..................................................................
C
      SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)
      REAL X(200),Y(200),XX(N),YY(N)
      LOGICAL MX,MY,NX,NY
      INTEGER O

Cf2py intent(out) INOUT
Cf2py integer optional, intent(in) :: N=shape(XX,0)


C      OUTPUT UNIT FOR PRINTED MESSAGES
      DATA O/6/
      MAXDIM=200
      IF(N.LE.MAXDIM)GO TO 6
      WRITE(O,7)
7     FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY.
     1RESULTS INVALID')
      RETURN
6     DO 1 I=1,N
      X(I)=XX(I)-PX
1     Y(I)=YY(I)-PY
      INOUT=-1
      DO 2 I=1,N
      J=1+MOD(I,N)
      MX=X(I).GE.0.0
      NX=X(J).GE.0.0
      MY=Y(I).GE.0.0
      NY=Y(J).GE.0.0
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3
      INOUT=-INOUT
      GO TO 2
3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
4     INOUT=0
      RETURN
5     INOUT=-INOUT
2     CONTINUE
      RETURN
      END


      subroutine extrap2(x,y,v,vout,mask,dir,L,M)
      ! 2d extrapolation

      integer L,M
      real x(L,M), y(L,M),v(L,M),xyratio,dir(4),vout(L,M)
      logical mask(L,M)

!Cf2py intent(in,out) v
Cf2py intent(out) vout
      vout=v

      xyratio=(maxval(x)-minval(x))/(maxval(y)-minval(y))

      call extrap_once(x,y,vout,mask,xyratio,dir,L,M)
      if (any(mask)) then
        call extrap_once(x,y,vout,mask,xyratio,dir,L,M)
      endif
      end



      subroutine extrap_once(x,y,v,mask,xyratio,dir,L,M)
      ! called by extrap2
      ! extrap_once may return data with mask! if so, rerun extrap_once
      ! mma, feb 2008

      integer L,M,i,j
      real x(L,M), y(L,M),v(L,M), A(L,M),B(L,M),C(L,M),D(L,M)
      real xyratio,dir(4)
      logical Amask(L,M),Bmask(L,M),Cmask(L,M),Dmask(L,M),mask(L,M)

      A=v
      B=v
      C=v
      D=v
      Amask=mask
      Bmask=mask
      Cmask=mask
      Dmask=mask

      do i=1,L
        do j=1,M
          if (Amask(i,j)) then
            call extrap_aux(x,y,A,Amask,xyratio,i,j,L,M)
          endif
        enddo
      enddo

      do i=L,1,-1
        do j=1,M
          if (Bmask(i,j)) then
            call extrap_aux(x,y,B,Bmask,xyratio,i,j,L,M)
          endif
        enddo
      enddo

      do i=1,L
        do j=M,1,-1
          if (Cmask(i,j)) then
            call extrap_aux(x,y,C,Cmask,xyratio,i,j,L,M)
          endif
        enddo
      enddo

      do i=L,1,-1
        do j=M,1,-1
          if (Dmask(i,j)) then
            call extrap_aux(x,y,D,Dmask,xyratio,i,j,L,M)
          endif
        enddo
      enddo

      do i=1,L
        do j=1,M
          v(i,j)=(A(i,j)*dir(1)+B(i,j)*dir(2)+C(i,j)*dir(3)
     &            +D(i,j)*dir(4))/(dir(1)+dir(2)+dir(3)+dir(4))
          if (Amask(i,j).or.Bmask(i,j).or.Cmask(i,j).or.Dmask(i,j)) then
            mask(i,j)=.true.
          else
            mask(i,j)=.false.
          endif
        enddo
      enddo

      end



      subroutine extrap_aux(x,y,v,mask,xyratio,i0,j0,L,M)
      !called by extrap_once
      !extraps for a single point, v(i0,j0)
      !mma, feb 2008

      integer L,M,i,j,i0,j0
      real x(L,M), y(L,M),v(L,M),d(4),u(4),xyratio
      logical mask(L,M)

      do i=1,4
        d(i)=0.
        u(i)=0.
      enddo

c      ! find d,u 1 and 3 (left, right)
        j=j0-1
        do while (j.ge.1)
          if (.not.mask(i0,j)) then
            u(1)=v(i0,j)
            d(1)=1./sqrt((x(i0,j0)-x(i0,j))**2 + (y(i0,j0)-y(i0,j))**2 )
            j=1
          endif
          j=j-1
        enddo

        j=j0+1
        do while (j.le.M)
          if (.not.mask(i0,j)) then
            u(3)=v(i0,j)
            d(3)=1./sqrt((x(i0,j0)-x(i0,j))**2 + (y(i0,j0)-y(i0,j))**2 )
            j=M
          endif
          j=j+1
        enddo

c     ! find d,u 2 and 4 (top, bottom)
        i=i0-1
        do while (i.ge.1)
          if (.not.mask(i,j0)) then
            u(2)=v(i,j0)
            d(2)=1./sqrt((x(i0,j0)-x(i,j0))**2 + (y(i0,j0)-y(i,j0))**2 )
            i=1
          endif
          i=i-1
        enddo

        i=i0+1
        do while (i.le.L)
          if (.not.mask(i,j0)) then
            u(4)=v(i,j0)
            d(4)=1./sqrt((x(i0,j0)-x(i,j0))**2 + (y(i0,j0)-y(i,j0))**2 )
            i=L
          endif
          i=i+1
        enddo

      if (sum(d).ne.0.) then
        d(1)=d(1)*xyratio
        d(2)=d(2)*(1.-xyratio)
        d(3)=d(3)*xyratio
        d(4)=d(4)*(1.-xyratio)
        v(i0,j0)=(u(1)*d(1)+u(2)*d(2)+u(3)*d(3)+u(4)*d(4))/
     &           (d(1)+d(2)+d(3)+d(4))

        mask(i0,j0)=.false.

      endif

      end


! alternative point in polygon test
! pnpoly works fine... so this is not needed for now!
!
! mma jan 2014.
!
!
!      subroutine inpolygon(xp,yp,x,y,res,n,p)
!      integer n, i,p
!      logical in_, res(p)
!      real xp(n), yp(n), x(p), y(p)
!
!Cf2py intent(out) res
!Cf2py integer optional, intent(in) :: n=shape(xp,0)
!Cf2py integer optional, intent(in) :: p=shape(x,0)
!
!      do i=1,p
!        call in_polygon(xp,yp,x(i),y(i),n,in_)
!        res(i)=in_
!      enddo
!      end
!
!
!      subroutine in_polygon(xp,yp,x,y,n,in_)
!      ! checks if 2d point (x,y) is inside the polygon (xp,yp)
!      ! mma 12-2007
!
!      integer n, m, i
!      logical in_
!      real xp(n), yp(n), x, y, xi(n-1), yi(n-1), xx(2),xp1,xp2,yp1,yp2
!      parameter (mask=-99.)
!
!Cf2py intent(out) in_
!
!      xp1=xp(1)
!      xp2=xp(1)
!      yp1=yp(1)
!      yp2=yp(1)
!      do i=2,n
!        xp1=min(xp1,xp(i))
!        xp2=max(xp2,xp(i))
!        yp1=min(yp1,yp(i))
!        yp2=max(yp2,yp(i))
!      enddo
!
!      if (x.ge.xp1 .and. x.le.xp2 .and. y.ge.yp1 .and. y.le.yp2) then
!
!        xx(1)=x
!        xx(2)=xp(1)
!        do i=2,n
!          if (xp(i).gt.xx(2)) then
!            xx(2)=xp(i)
!          endif
!        enddo
!        xx(2)=xx(2)+1.
!
!        call hmp(xp,yp,xx,y,xi,yi,n)
!        m=0
!        do i=1,n-1
!          if (xi(i).ne.mask) then
!            m=m+1
!          endif
!        enddo
!
!        if (m/2.eq.m/2.) then
!          in_=.false.
!        else
!          in_=.true.
!        endif
!
!      else
!        in_=.false.
!      endif
!
!
!      end
!
!
!      subroutine hmp(xp,yp,x,y,xi,yi,n)
!      ! called by in_polygon
!      ! horizontal meet point
!
!      integer n,i
!      real xp(n), yp(n), x(2), y, m, xi(n-1), yi(n-1)
!      real x1,x2,xp1,xp2,yp1,yp2
!      parameter (mask=-99.)
!
!Cf2py intent(out) xi, yi
!
!      x1=min(x(1),x(2))
!      x2=max(x(1),x(2))
!      do i=1,n-1
!        xp1=min(xp(i),xp(i+1))
!        xp2=max(xp(i),xp(i+1))
!
!        yp1=min(yp(i),yp(i+1))
!        yp2=max(yp(i),yp(i+1))
!
!        if ( xp2.lt.x1 .or. x1.gt.x2 .or. y.gt.yp2 .or. y.lt.yp1) then 
!          xi(i)=mask
!          yi(i)=mask
!        else
!          if (xp(i+1).eq.xp(i)) then
!            xi(i)=xp(i)
!            yi(i)=y
!          elseif (yp(i+1).eq.yp(i)) then
!            xi(i)=mask
!            yi(i)=mask
!          else
!            m=(yp(i+1)-yp(i))/(xp(i+1)-xp(i))
!            xi(i)=(y-yp(i))/m +xp(i)
!            yi(i)=y
!          endif
!        endif
!
!        if (xi(i).ne.mask.and.(xi(i).lt.max(xp1,x1) .or.
!     &                         xi(i).gt.min(xp2,x2) )) then
!          xi(i)=mask
!          yi(i)=mask
!        endif
!
!      enddo
!
!      end
