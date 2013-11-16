      subroutine meetpoint(x1,y1,x2,y2,res,n1,n2,N)
       integer n1, n2, i, j, N
       parameter (inf=9999, mask=999)
       real x1(n1), x2(n2), y1(n1), y2(n2),
     &      Xa,Ya,xa_,ya_,Xb,Yb,xb_,yb_,dx1,dx2,
     &      dy1,dy2,m1,m2,x,y
       real res(2,(n1-1)*(n2-1))
       logical cond

Cf2py intent(out) res,N
Cf2py intent(in) :: n1=shape(x1,0)
Cf2py intent(in) :: n2=shape(x2,0)

       do i=1,((n1-1)*(n2-1))
         res(1,i)=mask
         res(2,i)=mask
       enddo

       N=0
       do i=1,n1-1
         Xa=max(x1(i),x1(i+1))
         Ya=max(y1(i),y1(i+1))
         xa_=min(x1(i),x1(i+1))
         ya_=min(y1(i),y1(i+1))
         dx1=x1(i+1)-x1(i)
         dy1=y1(i+1)-y1(i)
         do j=1,n2-1
           Xb=max(x2(j),x2(j+1))
           Yb=max(y2(j),y2(j+1))
           xb_=min(x2(j),x2(j+1))
           yb_=min(y2(j),y2(j+1))
           dx2=x2(j+1)-x2(j)
           dy2=y2(j+1)-y2(j)

           if ((Xa.ge.xb_).and.(xa_.le.Xb).and.
     &         (Ya.ge.yb_).and.(ya_.le.Yb)) then

             if ((dx1.ne.0.).and.(dx2.ne.0)) then
               m1 = dy1/dx1
               m2 = dy2/dx2
               if (m1.eq.m2) then
                 x=inf
                 y=inf
               else
                  y = 1./(m1-m2) *
     &                ( m1*y2(j) - m2*y1(i) + m1*m2*(x1(i)-x2(j)))
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
               x=inf
               y=inf
             endif

             cond=(x.ge.max(xa_,xb_)).and.(x.le.min(Xa,Xb)).and.
     &            (y.ge.max(ya_,yb_)).and.(y.le.min(Ya,Yb))

             if (cond) then
               N=N+1
               res(1,N)=x
               res(2,N)=y
             endif

           endif

         enddo
       enddo

      end
