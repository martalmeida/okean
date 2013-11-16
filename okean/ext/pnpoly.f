!
! PNPOLY - Point Inclusion in Polygon Test
! W. Randolph Franklin (WRF)
! http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
!

C>>>PNP1                                                                
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
