program ptcloud

!*************************************************************************
!
! Program generates a set of points and their corresponding outward 
! pointing normal vectors. Points are sampled in a roughly
! equidistant fashion from a parametric surface defined by the
! Gielis' superformula.
!
! Note: Parameters a,b,m1,n1-n3 are associated with R1 in 3D superformula
! and parameters c,d,m2,n4-n6 are associated with R2 in 3D superformula
!
!*************************************************************************



 integer (kind = 4), parameter :: numAngle = 10000   !at most you can generate (numAngle)^2 points
 integer (kind = 4) h,i,j,k,l,p,t

!Following parameters are associated with Gielis' superformula and determine the 3D surface generated
!****************************************************************************************************
 real (kind = 8),parameter :: m1=5.0d+00
 real (kind = 8),parameter :: m2=1.0d+00
 real (kind = 8),parameter :: n1=0.1d+00
 real (kind = 8),parameter :: n2=1.7d+00
 real (kind = 8),parameter :: n3=1.7d+00
 real (kind = 8),parameter :: n4=0.3d+00
 real (kind = 8),parameter :: n5=0.5d+00
 real (kind = 8),parameter :: n6=0.5d+00
 real (kind = 8),parameter :: a=1.0d+00
 real (kind = 8),parameter :: b=1.0d+00
 real (kind = 8),parameter :: c=1.0d+00
 real (kind = 8),parameter :: d=1.0d+00 
!****************************************************************************************************
   
 real (kind = 8),parameter :: epsilon_=0.01d+00     !infinitesimal amount necessary to avoid singularity at theta=phi=0
 real (kind = 8),dimension(numAngle) :: theta,phi 
 real (kind = 8) delta_theta 
 real (kind = 8) delta_phi 
 real (kind = 8) pi
 
 real (kind = 8),allocatable,dimension(:,:) :: ptCld_xyz,nVec_xyz

 real (kind = 8),parameter :: S_theta=0.094d+00   !Initialize arc length. smaller => more points generated
 real (kind = 8),parameter :: S_phi=0.094d+00
 real (kind = 8) A_
 real (kind = 8) R1,dR1,R2,dR2
 real (kind = 8) dR1_1,dR1_21,dR1_22
 real (kind = 8) dR2_1,dR2_21,dR2_22 
 

 
 pi = 4.0d+00*atan(1.0d+00)
 A_ = (abs(cos(m2*epsilon_/4.0d+00)/c)**n5 + abs(cos(m2*epsilon_/4.0d+00)/d)**n6)**(-1/n4) 
 
 
!Following loops generate the necessary values of theta and phi to sample approximately equidistant points on 3D surface

 p = 1
 theta(1) = epsilon_
 
   do i=1,numAngle
 
 
       R1 = (abs(cos(m1*theta(i)/4.0d+00)/a)**n2 + abs(cos(m1*theta(i)/4.0d+00)/b)**n3)**(-1/n1)       
       dR1_1 = (1/n1)*(abs(cos(m1*theta(i)/4.0d+00)/a)**n2 + abs(sin(m1*theta(i)/4.0d+00)/b)**n3)**(-1/n1-1)       
       dR1_21 = m1*n2*abs(cos(m1*theta(i)/4.0d+00)/a)**n2*sin(m1*theta(i)/4.0d+00)/(4.0d+00*cos(m1*theta(i)/4.0d+00))       
       dR1_22 = -(m1*n3*abs(sin(m1*theta(i)/4.0d+00)/b)**n2*cos(m1*theta(i)/4.0d+00)/a)/(4.0d+00*sin(m1*theta(i)/4.0d+00))
       dR1 = dR1_1*(dR1_21+dR1_22) 
       
 
       call deltaTheta(delta_theta,R1,dR1,A_,S_theta)
       
       theta(i+1) = theta(i) + delta_theta
              
       if( (theta(i+1) .eq. 2*pi/m1) .or. (theta(i+1) .eq. 4*pi/m1)) then
       theta(i+1) = theta(i+1) + epsilon_    
       end if
        
       if( theta(i+1) .ge. pi ) then
       GOTO 10
       else
       p=p+1
       end if
          
       
   end do       
         
10 continue  

 theta(p) = -1*epsilon_
 
 l = 1 
   do i=p,numAngle
  
       R1 = (abs(cos(m1*theta(i)/4.0d+00)/a)**n2 + abs(cos(m1*theta(i)/4.0d+00)/b)**n3)**(-1/n1)       
       dR1_1 = (1/n1)*(abs(cos(m1*theta(i)/4.0d+00)/a)**n2 + abs(sin(m1*theta(i)/4.0d+00)/b)**n3)**(-1/n1-1)       
       dR1_21 = m1*n2*abs(cos(m1*theta(i)/4.0d+00)/a)**n2*sin(m1*theta(i)/4.0d+00)/(4.0d+00*cos(m1*theta(i)/4.0d+00))       
       dR1_22 = -(m1*n3*abs(sin(m1*theta(i)/4.0d+00)/b)**n2*cos(m1*theta(i)/4.0d+00)/a)/(4.0d+00*sin(m1*theta(i)/4.0d+00))
       dR1 = dR1_1*(dR1_21+dR1_22) 
       
       call deltaTheta(delta_theta,R1,dR1,A_,S_theta)
       
       theta(i+1) = theta(i) - delta_theta
       
       if( (theta(i+1) .eq. -1*2*pi/m1) .or. (theta(i+1) .eq. -1*4*pi/m1)) then
       theta(i+1) = theta(i+1) - epsilon_    
       end if       
       
       if( theta(i+1) .le. -1*pi ) then
       GOTO 20
       else
       l = l+1
       end if            
       
   end do       

20 continue

 k = 1
 phi(1) = epsilon_
 
   do j=1,p+l
      do i=1,numAngle
      
       R1 = (abs(cos(m1*theta(j)/4.0d+00)/a)**n2 + abs(cos(m1*theta(j)/4.0d+00)/b)**n3)**(-1/n1)

       R2 = (abs(cos(m2*phi(i)/4.0d+00)/c)**n5 + abs(cos(m2*phi(i)/4.0d+00)/d)**n6)**(-1/n4)       
       dR2_1 = (1/n3)*(abs(cos(m2*phi(i)/4.0d+00)/c)**n5 + abs(sin(m2*phi(i)/4.0d+00)/d)**n6)**(-1/n4-1)       
       dR2_21 = m2*n5*abs(cos(m2*phi(i)/4.0d+00)/c)**n5*sin(m2*phi(i)/4.0d+00)/(4.0d+00*cos(m2*phi(i)/4.0d+00))       
       dR2_22 = -(m2*n6*abs(sin(m2*phi(i)/4.0d+00)/d)**n2*cos(m1*phi(i)/4.0d+00)/c)/(4.0d+00*sin(m2*phi(i)/4.0d+00))
       dR2 = dR2_1*(dR2_21+dR2_22) 


       call deltaPhi(delta_phi,R2,dR2,R1,S_phi,phi(i)) 
       
       phi(i+1) = phi(i) + delta_phi
       
       if( (phi(i+1) .eq. 2*pi/m1)) then
       phi(i+1) = phi(i+1) + epsilon_    
       end if
       
       if( phi(i+1) .ge. pi/2 ) then
       GOTO 30
       else
       k = k+1
       end if
       
       
       end do        
   end do  

30 continue
       

 t = 1
 phi(k) = -1*epsilon_
 
   do j=1,p+l
      do i=k,numAngle
      
       R1 = (abs(cos(m1*theta(j)/4.0d+00)/a)**n2 + abs(cos(m1*theta(j)/4.0d+00)/b)**n3)**(-1/n1)

       R2 = (abs(cos(m2*phi(i)/4.0d+00)/c)**n5 + abs(cos(m2*phi(i)/4.0d+00)/d)**n6)**(-1/n4)       
       dR2_1 = (1/n3)*(abs(cos(m2*phi(i)/4.0d+00)/c)**n5 + abs(sin(m2*phi(i)/4.0d+00)/d)**n6)**(-1/n4-1)       
       dR2_21 = m2*n5*abs(cos(m2*phi(i)/4.0d+00)/c)**n5*sin(m2*phi(i)/4.0d+00)/(4.0d+00*cos(m2*phi(i)/4.0d+00))       
       dR2_22 = -(m2*n6*abs(sin(m2*phi(i)/4.0d+00)/d)**n2*cos(m1*phi(i)/4.0d+00)/c)/(4.0d+00*sin(m2*phi(i)/4.0d+00)) 
       dR2 = dR2_1*(dR2_21+dR2_22) 


       call deltaPhi(delta_phi,R2,dR2,R1,S_phi,phi(i)) 
       
       phi(i+1) = phi(i) - delta_phi
       
       if( (phi(i+1) .eq. -1*2*pi/m1)) then
       phi(i+1) = phi(i+1) - epsilon_    
       end if       
             
       if( phi(i+1) .le. -1*pi/2 ) then
       GOTO 40
       else
       t = t+1
       end if
       
       
      end do     
   end do  
   

40 continue


allocate(ptCld_xyz(3,(k+l)*(p+t)))
allocate(nVec_xyz(3,(k+l)*(p+t)))

   call cldGen(ptCld_xyz,nVec_xyz,theta,phi,h,k,l,p,t,m1,m2,n1,n2,n3,n4,n5,n6,a,b,c,d)
   call writeClds(ptCld_xyz,nVec_xyz,h)


end program


!Calculates first order approximation to delta_theta. For surfaces with high curvature, higher order
!approximations may need to be added to this subroutine or special consideration given to areas near singularities
subroutine deltaTheta(delta_theta,R1,dR1,A_,S) 

!Inputs
 real (kind = 8) R1,dR1
 real (kind = 8) A_
 real (kind = 8) S 


!Input/Outputs
 real (kind = 8) delta_theta
 

 delta_theta = S/(A_*sqrt((dR1)**2 + (R1)**2))
 

end subroutine

!Calculates first order approximation to delta_phi. For surfaces with high curvature, higher order
!approximations may need to be added to this subroutine or special consideration given to areas near singularities
subroutine deltaPhi(delta_phi,R2,dR2,R1,S,phi) 

!Inputs
 real (kind = 8) R2,dR2
 real (kind = 8) R1
 real (kind = 8) S 
 real (kind = 8) phi


!Input/Outputs
 real (kind = 8) delta_phi
 

 delta_phi = S/(sqrt(R1*(dR2*cos(phi) - R2*sin(phi))**2 + (dR2*sin(phi) + R2*cos(phi))**2))

 
end subroutine


!Generates point cloud and corresponding normal vectors
subroutine cldGen(ptCld_xyz,nVec_xyz,theta,phi,h,k,l,p,t,m1,m2,n1,n2,n3,n4,n5,n6,a,b,c,d)


!Input
 integer (kind = 4) h,k,l,p,t

 real (kind = 8), dimension(350) :: theta,phi 
 real (kind = 8) a,b,c,d
 real (kind = 8) m1,m2
 real (kind = 8) n1,n2,n3,n4,n5,n6
  
!Local
 integer (kind = 4) i,j   
 
 real (kind = 8) R1,dR1,R2,dR2
 real (kind = 8) dR1_1,dR1_21,dR1_22
 real (kind = 8) dR2_1,dR2_21,dR2_22 
 
!Output
 real (kind = 8), dimension(3,(k+l)*(p+t)) :: ptCld_xyz
 real (kind = 8), dimension(3,(k+l)*(p+t)) :: nVec_xyz
 
     
 h = 1
   do i=1,p+l
      do j=1,k+t
      
       R1 = (abs(cos(m1*theta(i)/4.0d+00)/a)**n2 + abs(cos(m1*theta(i)/4.0d+00)/b)**n3)**(-1/n1)      
       dR1_1 = (1/n1)*(abs(cos(m1*theta(i)/4.0d+00)/a)**n2 + abs(sin(m1*theta(i)/4.0d+00)/b)**n3)**(-1/n1-1)      
       dR1_21 = m1*n2*abs(cos(m1*theta(i)/4.0d+00)/a)**n2*sin(m1*theta(i)/4.0d+00)/(4.0d+00*cos(m1*theta(i)/4.0d+00))       
       dR1_22 = -(m1*n3*abs(sin(m1*theta(i)/4.0d+00)/b)**n2*cos(m1*theta(i)/4.0d+00)/a)/(4.0d+00*sin(m1*theta(i)/4.0d+00))
       dR1 = dR1_1*(dR1_21+dR1_22) 
      
       R2 = (abs(cos(m2*phi(j)/4.0d+00)/c)**n5 + abs(cos(m2*phi(j)/4.0d+00)/d)**n6)**(-1/n4)       
       dR2_1 = (1/n3)*(abs(cos(m2*phi(j)/4.0d+00)/c)**n5 + abs(sin(m2*phi(j)/4.0d+00)/d)**n6)**(-1/n4-1)       
       dR2_21 = m2*n5*abs(cos(m2*phi(j)/4.0d+00)/c)**n5*sin(m2*phi(j)/4.0d+00)/(4.0d+00*cos(m2*phi(j)/4.0d+00))       
       dR2_22 = -(m2*n6*abs(sin(m2*phi(j)/4.0d+00)/d)**n2*cos(m1*phi(j)/4.0d+00)/c)/(4.0d+00*sin(m2*phi(j)/4.0d+00)) 
       dR2 = dR2_1*(dR2_21+dR2_22)
      
       ptCld_xyz(1,h) = R1*cos(theta(i))*R2*cos(phi(j))
       ptCld_xyz(2,h) = R1*sin(theta(i))*R2*cos(phi(j))
       ptCld_xyz(3,h) = R2*sin(phi(j))            
       
       nVec_xyz(1,h) = R2*cos(phi(j))*(dR1*sin(theta(i))+R1*cos(theta(i)))*(dR2*sin(phi(j))+R2*cos(phi(j)))
       nVec_xyz(2,h) = -1*R2*cos(phi(j))*(dR1*cos(theta(i))-R1*sin(theta(i)))*(dR2*sin(phi(j))+R2*cos(phi(j)))
       nVec_xyz(3,h) = R1*R2*sin(theta(i))*cos(phi(j))*(dR1*cos(theta(i))-R1*sin(theta(i)))*(-1*R2*sin(phi(j))+dR2*cos(phi(j)))&
       & - R1*R2*(cos(phi(j)))**2*(R1*cos(theta(i))+dR1*sin(theta(i)))*(-1*R2*sin(phi(j))+dR2*cos(phi(j)))
       
       nVec_xyz(:,h) = nVec_xyz(:,h)/sqrt((nVec_xyz(1,h))**2 + (nVec_xyz(2,h))**2 + (nVec_xyz(3,h))**2)
            
       h = h+1
      end do
   end do   
  
end subroutine


!Writes points and normals to files
subroutine writeClds(ptCld_xyz,nVec_xyz,h)


 integer (kind = 4) h,i
 real (kind = 8),dimension(3,h) :: ptCld_xyz,nVec_xyz
 
 
   open(25,file='ptcldx.txt')
   open(27,file='ptcldy.txt')
   open(28,file='ptcldz.txt')
   open(29,file='ptcld.txt')
   open(26,file='norm_ptcld.txt') 

   do i=1,h-1
       write(25,*)ptCld_xyz(1,i)
       write(27,*)ptCld_xyz(2,i)
       write(28,*)ptCld_xyz(3,i)
       write(26,*)nVec_xyz(1,i),nVec_xyz(2,i),nVec_xyz(3,i),i  
       write(29,*)ptCld_xyz(1,i),ptCld_xyz(2,i),ptCld_xyz(3,i),i   
   end do    
   
   close(25)
   close(26)
   close(27)
   close(28)
   close(29)
end subroutine