program mom3d

!*****************************************************************************80
!
!   mom3d is the main program for the computing of the field scattered by a 
!  3D object.
!
!  Discussion:
!    
!    mom3d is based on Rao et al paper :
!    
!    Electromagnetic Scattering by Surfaces of Arbitrary Shape 
!    IEEE TAP, VOL AP-30,NO 3, 409-418 pp,(1982)
!    
!     on the numerical implementation from Makarov's book :
!
!    Antenna and EM Modelling with Matlab (Wiley,2003)
!
!    and on some of the subroutines provided by John Burkardt
!    
!    http://people.sc.fsu.edu/~burkardt
!
!    This version considers only PEC surfaces. More specifically, it treats the benchmark problem of a
!    perfectly conducting plane plate.                         
!
!
!  Date:
!
!    03 february 2009
!
!  Modified: 
!
!    23/01/2012
!
!
!  Author:
!
!    Demetrio Macias
!
!  Usage:
!
!    Scattering simulations
!
  implicit none

  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), parameter :: nelemx = 10
  integer ( kind = 4 ), parameter :: nelemy = 10
  character ( len =4 ), parameter :: type_surf='esf'                ! type_surf='esf' for the sphere and type_surf='gri' for a plane grid
  character (len = 4), parameter ::type_material='dec'              !type_material='pec' for perfect conductor and 'dec' for dielectric 
  logical, parameter :: header = .true.
  integer ( kind = 4 ) i,j,m,n,r, node, node_num 
  ! (P for Makarov= node_num for Burkardt) is the number of nodes
  ! (N for Makarov=element_node for Burkardt) is the number of triangles
        
  integer ( kind = 4 ) EdgesTotal                                   ! total number of common edges
  integer ( kind = 4 ) TrianglesTotal 
  
  integer ( kind = 4 ), allocatable, dimension (:, :) :: t          !(element_node for Burkardt) this array contains the number of nodes for each triangle
  real    ( kind = 8 ), allocatable, dimension(:,:) :: p            ! (node_xy for Burkardt) this array contains the coordinates of each node
  real    ( kind = 8 ), allocatable, dimension(:,:) :: Vec1,Vec2,Vec1pVec2
  integer ( kind = 4 ), allocatable, dimension (:,:) :: Edge_,tempEdge_       ! it contains the edge first and second numbers
  integer ( kind = 4 ), allocatable, dimension (:) :: TrianglePlus  ! Plus Triangle number
  integer ( kind = 4 ), allocatable, dimension (:) :: TriangleMinus,tempTriangleMinus ! Minus Triangle number 
  real    ( kind = 8 ), allocatable, dimension (:) :: Edgelength,tempEdgelength    ! Edge length
  real    ( kind = 8 ), allocatable, dimension (:) :: Area          ! Area of each triangle
  real    ( kind = 8 ), allocatable, dimension (:,:) :: Center      ! Center of each triangle
  real     ( kind = 8 ), allocatable, dimension (:,:,:) :: Center_  ! 9 sub-triangle midpoints

!****************************************************************************
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: row
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable,  dimension (:)  :: ind
!****************************************************************************

!****************************************************************************
  integer  ( kind = 4 ) n1,n2,n3 
  real     ( kind = 8 ), dimension (3)  :: r1,r2,r3
  real     ( kind = 8 ), dimension (3)  :: r12,r23,r13
  real     ( kind = 8 ), dimension (3)  :: C1,C2,C3,C4,C5,C6
  real     ( kind = 8 ), dimension (3)  :: a1,a2,a3,a4,a5,a6,a7,a8,a9
  integer  ( kind = 4 ) NoPlus,NoMinus 
  real     ( kind = 8 ), allocatable, dimension (:,:) :: RHO_Plus
  real     ( kind = 8 ), allocatable, dimension (:,:) :: RHO_Minus
  real     ( kind = 8 ), allocatable, dimension (:) :: FreeVertex
  real     ( kind = 8 ), allocatable, dimension (:,:,:):: RHO__Plus,RHO__Minus

!***************************************************************************

!***************************************************************************

  real     ( kind = 8 )  f,epsilon_,epsilon_2,mu_,pi,theta,phi
  real     ( kind = 8 )  c_,c_2,eta_,Emtheta,Emphi
  real     ( kind = 8 )  omega,k,k2,Constant1,Factor,Constant3,Constant4
  real     ( kind = 8 ), allocatable, dimension (:,:,:):: RHO_P, RHO_M
  complex  ( kind = 8 ), allocatable, dimension (:,:) :: Factord
  complex  ( kind = 8 )   Kc,Kc2,im1,Constant2 ,im2
  complex  ( kind = 8 ), allocatable, dimension (:,:) ::  FactorA, FactorFi
  complex  ( kind = 8 ), allocatable, dimension (:,:) ::  FactorA2, FactorFi2
  complex  ( kind = 8 ), allocatable, dimension (:,:) ::  FactorA2_2, FactorFi_2,FactorFi2_2
  complex  ( kind = 8 ), allocatable, dimension(:,:):: Matimp, dMatimp 
  complex  ( kind = 8 ), allocatable, dimension(:,:):: DMatsub1
  complex  ( kind = 8 ), allocatable, dimension(:,:):: DMatsub2
!***************************************************************************
!***************************************************************************
  integer ( kind = 4 )   ::   INFO  
  integer ( kind = 4 ),allocatable,dimension(:) ::  IPIV,IPIVd
  complex ( kind = 8 ) :: ScalarProduct,ScalarPlus,ScalarMinus
  complex ( kind = 8 ), allocatable,dimension(:)::  V,VH,Vd   
  real    ( kind = 8 ), allocatable, dimension (:,:) ::transpose_center,transp_row
  integer ( kind = 4 ),dimension (1,3):: d,Pol,kv,kv2
  complex ( kind = 8 ), dimension(3,1):: EmPlus, EmMinus 
  complex ( kind = 8 ), dimension(3,1):: HPlus, HMinus
  real    ( kind = 8 ), dimension(3,1)::temp_RHO_Plus,temp_RHO_Minus
!***************************************************************************

!***************************************************************************

  real    ( kind = 8 ), allocatable,dimension(:)::CurrentNorm,JdCurrent,MdCurrent
  complex ( kind = 8 ) :: IE,IE2,ME
  complex ( kind = 8 ), allocatable,dimension(:):: compt  


  
!***************************************************************************

!****************************************************************************
  integer ( kind = 4 )   :: Kcompt, Index
  real    ( kind = 8 ) :: x0,x1,y0,y1
  real    ( kind = 8 ),allocatable, dimension(:)::x,grand_x
  real    ( kind = 8 ),allocatable, dimension(:,:)::temp_x, Dist
  real    ( kind = 8 ),allocatable, dimension(:)::y,grand_y
  real    ( kind = 8 ),allocatable, dimension(:,:)::temp_y
  real    ( kind = 8 ),dimension(7)::xi,gxi
  real    ( kind = 8 ),dimension(8)::yi,gyi
!****************************************************************************
  integer ( kind = 4 ),parameter :: N_pts = 20  
    
  complex ( kind = 8 ),allocatable, dimension(:,:,:):: nField 
  complex (kind = 8),dimension(3) :: fField



! BEGINS THE GENERATION OF THE GEOMETRY (TRIANGULATION)


 
  if (type_surf.eq.'gri') then 

!  How many elements are there?
!
     call grid_t3_element_num ( nelemx, nelemy,  TrianglesTotal )
     allocate ( t(element_order, TrianglesTotal) )
!
!  How many nodes are there?
!
     call grid_t3_node_num ( nelemx, nelemy, node_num )
     allocate ( p(3,node_num) )
!
!  Get the nodes that make up each element.
!
     call grid_t3_element ( nelemx, nelemy, t )
         
     call grid(p,nelemx,nelemy,node_num)
!
!  Write the elements and nodes to files.
!
     call r8table_write ( 'grid_t3_nodes.txt', 3, node_num, p,header )

     call itable_write ( 'grid_t3_elements.txt', element_order, TrianglesTotal,t, header )

  else

     call sphere_grid_t3_element_num ( nelemx, nelemy, TrianglesTotal )
     allocate ( t(1:element_order,1:TrianglesTotal) )

     call sphere_grid_t3_node_num ( nelemx, nelemy, node_num )
     allocate ( p(1:3,1:node_num) )

     write ( *, '(a)' ) ' '
     write ( *, '(a,i8)' ) '  Expected number of nodes =    ', node_num
     write ( *, '(a,i8)' ) '  Expected number of elements = ', TrianglesTotal!element_num

     call sphere_grid_t3_element ( nelemx, nelemy, t )



     call sphere_grid_t3_node_xyz ( nelemx, nelemy, p )


!
!  Write the elements and nodes to files.
!
     call r8table_write ( 'sphere_t3_nodes.txt', 3, node_num, p,header )

     call itable_write ( 'sphere_t3_elements.txt', element_order, TrianglesTotal,t, header )


  endif




  TrianglesTotal=size(t, DIM=2)

  allocate ( Vec1(3,1) )
  allocate ( Vec2(3,1) )
  allocate ( Vec1pVec2(3,1) )
  allocate(Area(TrianglesTotal))
  allocate(Center(3,TrianglesTotal))

  do m=1,TrianglesTotal

     Vec1(:,1)=p(:,t(1,m))-p(:,t(2,m))
     Vec2(:,1)=p(:,t(3,m))-p(:,t(2,m))
   
     Vec1pVec2(1,1)= Vec1(2,1)*Vec2(3,1)- Vec1(3,1)*Vec2(2,1)
     Vec1pVec2(2,1)=-( Vec1(1,1)*Vec2(3,1)- Vec1(3,1)*Vec2(1,1))
     Vec1pVec2(3,1)=Vec1(1,1)*Vec2(2,1)- Vec1(2,1)*Vec2(1,1)
   
     Area(m)=(0.5d+00)*sqrt(sum((Vec1pVec2(:,1)**2),DIM=1))
     Center(:,m)=(1.0d+00/3.0d+00)*(p(:,t(1,m))+p(:,t(2,m))+p(:,t(3,m)))
   
  enddo

  allocate ( triangle_neighbor(1:3,1:TrianglesTotal) )
  allocate ( row(3*TrianglesTotal,4) )

  call triangulation_order3_neighbor_triangles ( TrianglesTotal, &
       t, triangle_neighbor,row )
 
  r=0
  do m=1,3*TrianglesTotal
     do n=m+1,3*TrianglesTotal
        if (all(row(m,1:2).eq.row(n,1:2)).eqv..true.) then
           r=r+1
        endif
     enddo
  enddo
  EdgesTotal=r


  allocate (Edge_ (EdgesTotal,2)) 
  allocate (tempEdge_ (EdgesTotal,2))      
  allocate (TrianglePlus (EdgesTotal)) 
  allocate (TriangleMinus (EdgesTotal))
  allocate (tempTriangleMinus (EdgesTotal))  
  allocate (Edgelength (EdgesTotal)) 
  allocate (tempEdgelength (EdgesTotal)) 

  r=0
  do m=1,3*TrianglesTotal
     do n=m+1,3*TrianglesTotal
        if (all(row(m,1:2).eq.row(n,1:2)).eqv..true.) then
           r = r+1
           Edge_(r,1:2) = row(m,1:2)
           TrianglePlus(r) = row(m,4)
           TriangleMinus(r) = row(n,4)
           Edgelength(r) = sqrt(sum((p(:,Edge_(r,1))-p(:,Edge_(r,2)))**2, dim=1))
        endif
     enddo
  enddo
 
  allocate(ind(EdgesTotal))

  call quick_sort(TrianglePlus,ind,EdgesTotal)

  do m=1, EdgesTotal
     tempTriangleMinus(m)=TriangleMinus(ind(m))
     tempEdgelength(m)= Edgelength(ind(m))
     tempEdge_(m,:)=Edge_(ind(m),1:2)
  enddo
  TriangleMinus= tempTriangleMinus
  Edgelength=tempEdgelength
  Edge_=tempEdge_

 


  allocate(Center_ (3,9,TrianglesTotal)) 
 

  do m=1,TrianglesTotal
     n1             = t(1,m)
     n2             = t(2,m)
     n3             = t(3,m) 
     r1             = p(:,n1)
     r2             = p(:,n2)
     r3             = p(:,n3)
     r12            = r2-r1
     r23            = r3-r2
     r13            = r3-r1
     C1             = r1+(r12/3.0d+00)
     C2             = r1+(2*r12/3.0d+00)
     C3             = r2+(r23/3.0d+00)
     C4             = r2+(2.0d+00*r23/3.0d+00)
     C5             = r1+(r13/3.0d+00)
     C6             = r1+(2.0d+00*r13/3.0d+00)
     a1             = 1.0d+00/3.0d+00*(C1+C5+r1)
     a2             = 1.0d+00/3.0d+00*(C1+C2+Center(:,m))
     a3             = 1.0d+00/3.0d+00*(C2+C3+r2)
     a4             = 1.0d+00/3.0d+00*(C2+C3+Center(:,m))
     a5             = 1.0d+00/3.0d+00*(C3+C4+Center(:,m))
     a6             = 1.0d+00/3.0d+00*(C1+C5+Center(:,m))
     a7             = 1.0d+00/3.0d+00*(C5+C6+Center(:,m))
     a8             = 1.0d+00/3.0d+00*(C4+C6+Center(:,m))
     a9             = 1.0d+00/3.0d+00*(C4+C6+r3)
     Center_(:,1,m) = a1
     Center_(:,2,m) = a2
     Center_(:,3,m) = a3
     Center_(:,4,m) = a4
     Center_(:,5,m) = a5
     Center_(:,6,m) = a6
     Center_(:,7,m) = a7
     Center_(:,8,m) = a8
     Center_(:,9,m) = a9
enddo

allocate(RHO_Plus(3,EdgesTotal))
allocate(RHO_Minus(3,EdgesTotal))
allocate(FreeVertex(3))
allocate(RHO__Plus(3,9,EdgesTotal))
allocate(RHO__Minus(3,9,EdgesTotal))

!PLUS
do m=1,EdgesTotal
   NoPlus=TrianglePlus(m)
   n1=t(1,NoPlus)
   n2=t(2,NoPlus)
   n3=t(3,NoPlus) 
   if((n1.NE.Edge_(m,1)).and.(n1.NE.Edge_(m,2))) then 
      NODE=n1 
   endif
   if((n2.NE.Edge_(m,1)).and.(n2.NE.Edge_(m,2))) then 
      NODE=n2 
   endif
   if((n3.NE.Edge_(m,1)).and.(n3.NE.Edge_(m,2))) then 
      NODE=n3 
   endif
   FreeVertex=p(:,NODE)
    
   RHO_Plus(:,m)   = (Center(:,NoPlus)-FreeVertex)!/Area(NoPlus)
   !Nine rho's of the "plus" triangle
   RHO__Plus(:,:,m)  = +Center_(:,:,NoPlus)-spread(FreeVertex,Dim=2,NCOPIES=9)
enddo
!MINUS



do m=1,EdgesTotal
   NoMinus=TriangleMinus(m)
   n1=t(1,NoMinus)
   n2=t(2,NoMinus)
   n3=t(3,NoMinus) 
   if((n1.NE.Edge_(m,1)).AND.(n1.NE.Edge_(m,2))) then 
      NODE=n1
   endif
   if((n2.NE.Edge_(m,1)).AND.(n2.NE.Edge_(m,2))) then 
      NODE=n2
   endif
   if((n3.NE.Edge_(m,1)).AND.(n3.NE.Edge_(m,2))) then 
      NODE=n3 
   endif
   FreeVertex=p(:,NODE)
   
   RHO_Minus(:,m)   = (-Center(:,NoMinus) + FreeVertex)!/Area(NoMinus)
    !Nine rho's of the "minus" triangle
   RHO__Minus(:,:,m)  = -Center_(:,:,NoMinus) + spread(FreeVertex,Dim=2,NCOPIES=9)
enddo


!END OF GEOMETRY GENERATION


!BEGINS THE GENERATION OF THE IMPEDANCE MATRIX


!EM parameters (f=3e8 means that f=300 MHz) 
f           =3e8  
epsilon_    =8.854e-012
epsilon_2   =0.5*8.854e-012
mu_         =1.257e-006

!Speed of light
 c_=1/sqrt(epsilon_*mu_)
 c_2=1/sqrt(epsilon_2*mu_)

!Free-space impedance 
eta_=sqrt(mu_/epsilon_)

!Contemporary variables - introduced to speed up 
!the impedance matrix calculation
pi=4.0d+00*atan(1.0d+00)
im1=(0.0d+00,1.0d+00)
im2=(1.0d+00,0.0d+00)
omega       =2.0d+00*pi*f                                            
k           =omega/c_
k2          =omega/c_2
Kc          =im1*k
Kc2         =im1*k2 
Constant1   =mu_/(4.0d+00*pi)
Constant2   =1.0d+00/(im1*4.0d+00*pi*omega*epsilon_)
Constant3   =1.0d+00/(4.0d+00*pi)
Constant4   =1.0d+00/(im1*4.0d+00*pi*omega*epsilon_2)
Factor      =1.0d+00/9.0d+00    

allocate (FactorA (1,EdgesTotal))
allocate (FactorFi (1,EdgesTotal))
allocate (FactorFi_2 (1,EdgesTotal)) 
allocate (Factord (1,EdgesTotal))
allocate (FactorA2 (1,EdgesTotal))
allocate (FactorA2_2 (1,EdgesTotal))
allocate (FactorFi2 (1,EdgesTotal))
allocate (FactorFi2_2 (1,EdgesTotal))   

FactorA(1,:)      =Factor*(im1*omega*EdgeLength(:))/(4.0d+00)*Constant1
FactorFi(1,:)     =Factor*(EdgeLength(:))*Constant2
FactorFi_2(1,:)   =Factor*(EdgeLength(:))*Constant4
Factord(1,:)      =Factor*(im2*EdgeLength(:))/(4.0d+00)*Constant3
FactorA2(1,:)     =FactorA(1,:)*(epsilon_/mu_)
FactorA2_2(1,:)   =FactorA(1,:)*(epsilon_2/mu_)
FactorFi2(1,:)    =FactorFi(1,:)*(mu_/epsilon_)
FactorFi2_2(1,:)  =FactorFi_2(1,:)*(mu_/epsilon_2)

allocate(RHO_P(3,9,EdgesTotal))
allocate(RHO_M(3,9,EdgesTotal))

do m=1,EdgesTotal
   RHO_P(:,:,m)=spread(RHO_Plus(:,m),DIM=2,NCOPIES=9)   ![3 9 EdgesTotal]
   RHO_M(:,:,m)=spread(RHO_Minus(:,m),DIM=2,NCOPIES=9)  ![3 9 EdgesTotal]
enddo
   
allocate(Matimp(EdgesTotal,EdgesTotal))
allocate(dMatimp(2*EdgesTotal,2*EdgesTotal))
    
 call impmet( Matimp,dMatimp,EdgesTotal,TrianglesTotal,EdgeLength,Area,Kc,Kc2,Center,Center_,TrianglePlus,TriangleMinus,&
        & RHO_P,RHO_M,RHO__Plus,RHO__Minus,FactorA,FactorFi,FactorFi_2,Factord,FactorA2,FactorA2_2,FactorFi2,FactorFi2_2)   


!ENDS THE GENERATION OF THE IMPEDANCE MATRIX

! BEGINS GENERATION OF THE INCIDENT FIELD


!Example: d=[0 0 -1] means that the incident signal
! is in the -z direction. 

!Plate - normal incidence


! The 3D incident field with arbitrary polarization. The variables Emtheta and Emphi correspond to E superscript-i subscript-V and
!E superscript-i subscript-H, respectively, in Gedney's notes on pg 36. I think some additional work needs to be done in order to
!determine the values of these constants for given values of theta and phi.
theta = pi   !theta_incident
phi = pi     !phi_incident
Emtheta = 1.0d+00
Emphi = 0.0d+00
d(1,:)=(/sin(theta)*cos(phi), sin(theta)*sin(phi), cos(phi)/)     
Pol(1,:)=(/Emtheta*cos(theta)*cos(phi)-Emphi*sin(phi), Emtheta*cos(theta)*sin(phi)+Emphi*cos(phi),&
& -1.0d+00*Emtheta*sin(theta)/)      



k=omega/c_
k2=omega/c_2
kv=k*d
kv2=k2*d

allocate(transpose_center(TrianglesTotal,3))
allocate(transp_row(1,3))
allocate(V(EdgesTotal))
allocate(VH(EdgesTotal))
allocate(Vd(2*EdgesTotal)) 
allocate(IPIV(EdgesTotal))
allocate(IPIVd(2*EdgesTotal))

transpose_center=transpose(Center)


transp_row(1,:)=transpose_center(TrianglePlus(m),:) 

do m=1,EdgesTotal   
   transp_row(1,:)=transpose_center(TrianglePlus(m),:) 
   ScalarProduct=sum(kv*transp_row)
   EmPlus =transpose(Pol)*exp(-im1*ScalarProduct)      
   transp_row(1,:)=transpose_center(TriangleMinus(m),:) 
   ScalarProduct=sum(kv*transp_row)
   EmMinus=transpose(Pol)*exp(-im1*ScalarProduct)      

   temp_RHO_Plus(:,1)=RHO_Plus(:,m)
   ScalarPlus =sum(EmPlus*temp_RHO_Plus)
   temp_RHO_Minus(:,1)=RHO_Minus(:,m)
   ScalarMinus=sum(EmMinus*temp_RHO_Minus)


      V(m)=EdgeLength(m)*(ScalarPlus/2+ScalarMinus/2)   
enddo


!H field generation. Loop determines H field using the fact that H=(1/k*mu)(k X E)
if(type_material.eq.'dec') then
do m=1,EdgesTotal   
   transp_row(1,:)=transpose_center(TrianglePlus(m),:) 
   ScalarProduct=sum(kv*transp_row)
   EmPlus =transpose(Pol)*exp(-im1*ScalarProduct)
   HPlus(1,1) = (kv(1,2)*EmPlus(3,1)-kv(1,3)*EmPlus(2,1))
   HPlus(2,1) = (kv(1,1)*EmPlus(3,1)-kv(1,3)*EmPlus(1,1))   
   HPlus(3,1) = (kv(1,1)*EmPlus(2,1)-kv(1,2)*EmPlus(1,1))      
   transp_row(1,:)=transpose_center(TriangleMinus(m),:) 
   ScalarProduct=sum(kv*transp_row)

   EmMinus=transpose(Pol)*exp(-im1*ScalarProduct)
   HMinus(1,1) = (kv(1,2)*EmMinus(3,1)-kv(1,3)*EmMinus(2,1))
   HMinus(2,1) = (kv(1,1)*EmMinus(3,1)-kv(1,3)*EmMinus(1,1))   
   HMinus(3,1) = (kv(1,1)*EmMinus(2,1)-kv(1,2)*EmMinus(1,1))

   temp_RHO_Plus(:,1)=RHO_Plus(:,m)
   ScalarPlus =sum(HPlus*temp_RHO_Plus)
   temp_RHO_Minus(:,1)=RHO_Minus(:,m)
   ScalarMinus=sum(HMinus*temp_RHO_Minus)

      VH(m)=EdgeLength(m)*(ScalarPlus/2+ScalarMinus/2)

end do
end if

   VH=VH*(1.0d+00/(k*mu_))

 !combined E-field and H-field vector for 'dec' case

do m=1,EdgesTotal
   Vd(m)=V(m)
end do

do m=EdgesTotal+1,2*EdgesTotal
   Vd(m)=VH(m-EdgesTotal)
end do


!END THE GENERATION OF THE INCIDENT FIELD
   
!BEGINGS THE SOLUTION OF THE LINEAR EQUATIONS SYSTEM  

!if(type_material.eq.'pec') then
 call  ZGESV( EdgesTotal, 1, Matimp, EdgesTotal, IPIV, V, EdgesTotal, INFO )


!else
 !'dec' case
 call ZGESV( 2*EdgesTotal, 1, dMatimp, 2*EdgesTotal, IPIVd, Vd, 2*EdgesTotal, INFO) 

!end if


!ENDS THE SOLUTION OF THE SYSTEM



   allocate(compt(3)) 
   allocate(CurrentNorm(size(t(3,:))))
   allocate(JdCurrent(size(t(3,:))))
   allocate(MdCurrent(size(t(3,:)))) 




!if(type_material.eq.'pec') then
   !%FINDING OF THE CURRENT DENSITY FOR EVERY TRIANGLE
   do j=1,size(t(3,:))
    compt=0!   i=[0 0 0]'
       do m=1,EdgesTotal
           IE=V(m)*EdgeLength(m)
           if(TrianglePlus(m).eq.j) then
               compt=compt+IE*RHO_Plus(:,m)/(2*Area(TrianglePlus(m)))
           endif
           if(TriangleMinus(m).eq.j) then
               compt=compt+IE*RHO_Minus(:,m)/(2*Area(TriangleMinus(m)))
           endif
      
       enddo
 

       CurrentNorm(j)=real(sqrt(sum(compt*conjg(compt))))
 
   enddo

!else
   do j=1,size(t(3,:))
       compt=0!   i=[0 0 0]'
          do m=1,EdgesTotal
              IE2=Vd(m)*EdgeLength(m)
              if(TrianglePlus(m).eq.j) then
                  compt=compt+IE2*RHO_Plus(:,m)/(2*Area(TrianglePlus(m)))
              endif
              if(TriangleMinus(m).eq.j) then
                  compt=compt+IE2*RHO_Minus(:,m)/(2*Area(TriangleMinus(m)))
              endif
      
          enddo
 

          JdCurrent(j)=real(sqrt(sum(compt*conjg(compt))))
 
   enddo

   do j=1,size(t(3,:))
          compt=0!   i=[0 0 0]'
             do m=EdgesTotal+1,2*EdgesTotal
                 ME=Vd(m)*EdgeLength(m-EdgesTotal)
                 if(TrianglePlus(m-EdgesTotal).eq.j) then
                     compt=compt+ME*RHO_Plus(:,m-EdgesTotal)/(2*Area(TrianglePlus(m-EdgesTotal)))
                 endif
                 if(TriangleMinus(m-EdgesTotal).eq.j) then
                     compt=compt+ME*RHO_Minus(:,m-EdgesTotal)/(2*Area(TriangleMinus(m-EdgesTotal)))
                 endif
      
             enddo
 

             MdCurrent(j)=real(sqrt(sum(compt*conjg(compt))))
 
   enddo


!end if

 allocate(nField(3,1,N_pts**2))

 call nearfield( nField,EdgesTotal,TrianglesTotal,EdgeLength,Area,Kc,Center_,TrianglePlus,TriangleMinus,&
     & RHO_P,RHO_M,RHO__Plus,RHO__Minus,Vd,t,element_order,f,epsilon_,mu_,c_,omega,k,pi,N_pts)

 call farfield(fField,EdgesTotal,TrianglesTotal,EdgeLength,Area,Kc,Center_,TrianglePlus,TriangleMinus,&
     & RHO_P,RHO_M,RHO__Plus,RHO__Minus,Vd,t,element_order,f,epsilon_,mu_,c_,omega,k,pi)    


!%Comparison with respect to published data

Kcompt=15
x0=minval(p(1,:))
x1=maxval(p(1,:))
y0=minval(p(2,:))
y1=maxval(p(2,:))
allocate(x(Kcompt+1))
allocate(temp_x(3,1))
allocate(Dist(size(temp_x,dim=1),TrianglesTotal))
allocate(grand_x(Kcompt+1))
allocate(y(Kcompt+1))
allocate(temp_y(3,1))
allocate(grand_y(Kcompt+1))

do n=1,Kcompt+1
    x(n)=x0+(n-1)*(x1-x0)/Kcompt
    temp_x(1,1)=x(n)
    temp_x(2,1)=0.0d+00
    temp_x(3,1)=0.0d+00
    Dist=spread(temp_x(:,1),dim=2,ncopies=TrianglesTotal)-Center
 
    Index=minloc( sum(Dist*Dist,dim=1),dim=1)
    grand_x(n)=CurrentNorm(Index)*eta_     !%eta_ means normalization over Hinc
enddo

do n=1,Kcompt+1
    y(n)=y0+(n-1)*(y1-y0)/Kcompt
 temp_y(1,1)=0.0d+00
    temp_y(2,1)=y(n)
    temp_y(3,1)=0.0d+00
    Dist=spread(temp_y(:,1),dim=2,ncopies=TrianglesTotal)-Center
  
    Index=minloc( sum(Dist*Dist,dim=1),dim=1)
    grand_y(n)=CurrentNorm(Index)*eta_     !%eta_ means normalization over Hinc   
enddo

if (type_surf.eq.'gri') then 

   open(25,file='plate_scatt_jx.dat')
   do i=1,size(x,dim=1)
      write(25,*)x(i),grand_x(i)
   enddo
   close(25)
   open(25,file='plate_scatt_jy.dat')
   do i=1,size(y,dim=1)
      write(25,*)y(i),grand_y(i)
   enddo
   close(25)
!!%Data of Catedra et al.

   xi(1)=0.125d+00 
   xi(2)=0.250d+00 
   xi(3)=0.375d+00 
   xi(4)=0.500d+00 
   xi(5)=0.625d+00 
   xi(6)=0.750d+00 
   xi(7)=0.875d+00 

   gxi(1)=1.625d+00
   gxi(2)=2.083d+00
   gxi(3)=2.604d+00
   gxi(4)=2.896d+00
   gxi(5)=2.604d+00
   gxi(6)=2.083d+00
   gxi(7)=1.625d+00

   open(25,file='plate_catedra_jx.dat')
   do i=1,7
      write(25,*)xi(i)-0.5d00,gxi(i)
   enddo
   close(25)

   yi(1)=0.125d+00*(0.0d+00+0.5d+00)
   yi(2)=0.125d+00*(1.0d+00+0.5d+00) 
   yi(3)=0.125d+00*(2.0d+00+0.5d+00) 
   yi(4)=0.125d+00*(3.0d+00+0.5d+00) 
   yi(5)=0.125d+00*(4.0d+00+0.5d+00) 
   yi(6)=0.125d+00*(5.0d+00+0.5d+00)
   yi(7)=0.125d+00*(6.0d+00+0.5d+00) 
   yi(8)=0.125d+00*(7.0d+00+0.5d+00)

   gyi(1)=4.833d+00
   gyi(2)=2.6666d+00
   gyi(3)=2.75d+00
   gyi(4)=2.896d+00
   gyi(5)=2.896d+00
   gyi(6)=2.75d+00
   gyi(7)=2.6666d+00
   gyi(8)=4.833d+00

   open(25,file='plate_catedra_jy.dat')
   do i=1,8
      write(25,*)yi(i)-0.5d+00,gyi(i)
   enddo
   close(25)
else

   open(25,file='sphere_scatt_jx.dat')
   do i=1,size(x,dim=1)
      write(25,*)x(i),grand_x(i)
   enddo
   close(25)
   open(25,file='sphere_scatt_jy.dat')
   do i=1,size(y,dim=1)
      write(25,*)y(i),grand_y(i)
   enddo
   close(25)
endif


  write ( *, '(a)' ) '  Normal end of execution.'

deallocate(x)
deallocate(temp_x)
deallocate(Dist)
deallocate(grand_x)
deallocate(y)
deallocate(temp_y)
deallocate(grand_y)
deallocate(compt) 
deallocate(CurrentNorm) 
deallocate(IPIV)
deallocate(transpose_center)
deallocate(transp_row)
deallocate(V) 
deallocate ( Vec1 )
deallocate ( Vec2 )
deallocate ( Vec1pVec2 )
deallocate(Area)
deallocate(Center)
deallocate(RHO_Plus)
deallocate(RHO_Minus)
deallocate(FreeVertex)
deallocate(RHO__Plus)
deallocate(RHO__Minus)
deallocate (FactorA )
deallocate (FactorFi)
deallocate (Factord)
deallocate (FactorA2)
deallocate (FactorFi2)
deallocate (FactorFi2_2)
deallocate (FactorFi_2)
deallocate (FactorA2_2)  
deallocate(RHO_P)
deallocate(RHO_M)
deallocate(Matimp)
deallocate(dMatimp)
deallocate(nField)    

end program mom3d
!*************************************************************************************************************************************************************************************************
!Following subroutines implement a barycentric approximation to near and far field scattering calculations
! Written by: Murphy Griffin
! Summer of 2012
!
Include 'xscatt_mom3d.f90'
!Include 'zegsv.f'

subroutine impmet( Matimp,dMatimp,EdgesTotal,TrianglesTotal,EdgeLength,Area,Kc,Kc2,Center,Center_,TrianglePlus,TriangleMinus,&
     & RHO_P,RHO_M,RHO__Plus,RHO__Minus,FactorA,FactorFi,FactorFi_2,Factord,FactorA2,FactorA2_2,FactorFi2,FactorFi2_2)   
  implicit none

 
!Inputs
  integer ( kind = 4 ) EdgesTotal                                  
  integer ( kind = 4 ) TrianglesTotal 
  
  integer ( kind = 4 ), dimension (EdgesTotal) :: TrianglePlus
 
  integer ( kind = 4 ), dimension (EdgesTotal) :: TriangleMinus 
  
  real    ( kind = 8 ),dimension(EdgesTotal):: Edgelength
  real    (kind = 8), dimension(TrianglesTotal)::Area
  real    ( kind = 8 ), dimension(3,TrianglesTotal)::Center
  real    ( kind = 8 ),dimension (3,9,TrianglesTotal)::Center_  
  
  real     ( kind = 8 ), dimension(3,9,EdgesTotal):: RHO_P
  real     ( kind = 8 ), dimension(3,9,EdgesTotal):: RHO_M  
  
  real    ( kind = 8 ),dimension(3,9,EdgesTotal):: RHO__Plus
  real    ( kind = 8 ),dimension (3,9,EdgesTotal)::RHO__Minus
  
  complex  ( kind = 8)   Kc,Kc2
  complex  ( kind = 8 ), dimension(1,EdgesTotal):: FactorA 
  complex  ( kind = 8 ), dimension(1,EdgesTotal):: FactorFi 
  complex  ( kind = 8 ), dimension (1,EdgesTotal):: Factord
  complex  ( kind = 8 ), dimension (1,EdgesTotal):: FactorA2
  complex  ( kind = 8 ), dimension (1,EdgesTotal):: FactorA2_2
  complex  ( kind = 8 ), dimension (1,EdgesTotal):: FactorFi2
  complex  ( kind = 8 ), dimension (1,EdgesTotal):: FactorFi_2
  complex  ( kind = 8 ), dimension (1,EdgesTotal):: FactorFi2_2

!Output
  complex (kind = 8), dimension(EdgesTotal,EdgesTotal):: Matimp
  complex (kind = 8), dimension(2*EdgesTotal,2*EdgesTotal):: dMatimp
  

!Local variables

  integer (kind = 4) i,j,k,p,n,l
  integer (kind = 4), dimension (5) :: scratch
  integer (kind = 4), dimension (EdgesTotal) :: scratch_pos
  integer ( kind = 4 ),allocatable, dimension (:) ::Plus,Minus

  real (kind = 8),dimension(EdgesTotal,1):: temp_Edgelength
  real (kind = 8),dimension(TrianglesTotal,1)::temp_Area
  real (kind = 8), dimension(1,9,TrianglesTotal) ::R
  real (kind = 8), dimension(3,9,TrianglesTotal) ::D
  real (kind = 8), dimension(3,9,EdgesTotal)::DP,DM                    !plus and minus D arrays
  real (kind = 8), dimension(3,9,EdgesTotal) ::RHOcD__P,RHOcD__M       !cross product of RHO and D arrays
  real (kind = 8), dimension(3,9,EdgesTotal) ::RP,RdP
  real (kind = 8), dimension(1,9,EdgesTotal) ::sumRPpRHO_P,sumRPpRHO_M
  real (kind = 8), dimension(1,9,EdgesTotal)::sumRdPpRHOcD__P,sumRdPpRHOcD__M 


  complex (kind = 8), dimension(1,9,TrianglesTotal) ::g,g2,gd,gd2
  complex (kind = 8), dimension(1,9,EdgesTotal) ::gP,gP2,gM,gM2
  complex (kind = 8), dimension(1,9,EdgesTotal) ::gdP,gdP2,gdM,gdM2              !plus and minus green's functions for 'dec' case
  complex (kind = 8), dimension(1,EdgesTotal) ::Fi,Fi2
  complex (kind = 8), dimension(1,EdgesTotal) ::sumagp,sumagm
  complex (kind = 8), dimension(1,EdgesTotal) ::sumagp2,sumagm2
  complex (kind = 8), dimension(EdgesTotal,1)::ZF,Z1
  complex (kind = 8), dimension(EdgesTotal,1)::ZF_2,ZF2_2,Z1_2,Z12_2        !transposes for 'dec' case
  complex (kind = 8), dimension(EdgesTotal,1)::Zd1                              

  complex (kind = 8), dimension(1,EdgesTotal)::A,A2
  complex (kind = 8), dimension(1,EdgesTotal)::Dmat,Dmat2 
  complex (kind = 8), dimension(EdgesTotal,EdgesTotal)::Matimp_2,Matimp2_2               !same as Matimp but different prefactor for 'dec' case
  complex (kind = 8), dimension(EdgesTotal,EdgesTotal)::DMatsub1,DMatsub2    !D submatrices for 'dec' case                                      

!Memory allocation
  do i=1,EdgesTotal
     do j=1,EdgesTotal
        Matimp(i,j)=(0.0d+00,0.0d+00)
     enddo
  enddo

  do i=1,2*EdgesTotal
     do j=1,2*EdgesTotal
        dMatimp(i,j)=(0.0d+00,0.0d+00)
     enddo
  enddo
!Loop over integration triangles


  do p=1,TrianglesTotal
    
     scratch(:)=0
     j=0
     do i=1, EdgesTotal
        if(TrianglePlus(i)-p.eq.0) then
           j=j+1
           scratch(j)=i
        end if
     enddo
    
     i=1
     allocate(Plus(j))
     do while(scratch(i).ne.0)
        Plus(i)=scratch(i)

        i=i+1
     enddo
     
     scratch(:)=0
     j=0
     do i=1, EdgesTotal
        if(TriangleMinus(i)-p.eq.0) then
           j=j+1
           scratch(j)=i
    
           scratch_pos(i)=i
       
        end if
     enddo
   
     i=1
     allocate(Minus(j))
     do while(scratch(i).ne.0)
        Minus(i)=scratch(i)
        i=i+1
     enddo
  
 !********************************************************************************************
 !variables used in both 'pec' and 'dec' cases

     D=spread(spread(Center(:,p),DIM=2,NCOPIES=9),DIM=3,NCOPIES=TrianglesTotal)-Center_
     
     R(1,:,:)=sqrt(sum(D**2,DIM=1))         !%[1 9 TrianglesTotal]
    
     g=exp(-Kc*R)/R                               !%[1 9 TrianglesTotal]
     g2=exp(-Kc2*R)/R 

!Plus and Minus Green's functions for both dielectric constants for two different media in the dielectric case.
     gP=g(:,:,TrianglePlus)                         !%[1 9 EdgesTotal]
     gP2=g2(:,:,TrianglePlus)

     gM=g(:,:,TriangleMinus)                        !%[1 9 EdgesTotal]
     gM2=g2(:,:,TriangleMinus)
 
     gd=-1*(1+Kc*R)/(R)**2*g
     gd2=-1*(1+Kc2*R)/(R)**2*g2 
   
     gdP=gd(:,:,TrianglePlus)
     gdP2=gd2(:,:,TrianglePlus)

     gdM=gd(:,:,TriangleMinus)
     gdM2=gd2(:,:,TriangleMinus)

     temp_Area(:,1)=Area    
     temp_Edgelength(:,1)=Edgelength
 !********************************************************************************************
 !Calculation of Impedance matrix for 'pec' case or (A-mn - B-mn) submatrices for 'dec' case refer to Sarkar's PMCHW formulation.
  
     sumagp=sum(gP,dim=2)
     sumagm=sum(gM,dim=2)
     sumagp2=sum(gP2,dim=2)
     sumagm2=sum(gM2,dim=2)
   
     Fi=sumagp-sumagm                             !%[1 1 EdgesTotal]
     Fi2=sumagp2-sumagm2  

     ZF= transpose(FactorFi*Fi)                    !%[EdgesTotal 1]
     ZF_2= transpose(FactorFi*Fi+FactorFi_2*Fi2)
     ZF2_2= transpose(FactorFi2*Fi+FactorFi2_2*Fi2)
   
     do k=1,size(Plus)
        n=Plus(k)
        
        RP=spread(RHO__Plus(:,:,n),DIM=3,NCOPIES=EdgesTotal)  !%[3 9 EdgesTotal]
        sumRPpRHO_P(1,:,:)=sum(RP*RHO_P,DIM=1)!*Area(n)
        sumRPpRHO_M(1,:,:)=sum(RP*RHO_M,DIM=1)!*Area(n)
        
        A=sum(gP*sumRPpRHO_P,DIM=2)+sum(gM*sumRPpRHO_M,DIM=2)
        A2=sum(gP2*sumRPpRHO_P,DIM=2)+sum(gM2*sumRPpRHO_M,DIM=2)
        Z1=transpose(FactorA*A)
        Z1_2=transpose(FactorA*A+FactorA*A2)
        Z12_2=transpose(FactorA2*A+FactorA2_2*A2)
        

        
        Matimp(:,n)=Matimp(:,n)+temp_EdgeLength(n,1)*(Z1(:,1)+ZF(:,1))
        Matimp_2(:,n)=Matimp_2(:,n)+temp_EdgeLength(n,1)*(Z1_2(:,1)+ZF_2(:,1))
        Matimp2_2(:,n)=Matimp_2(:,n)+temp_EdgeLength(n,1)*(Z12_2(:,1)+ZF2_2(:,1)) 
     
     enddo

     do k=1,size(Minus)
        
        n=Minus(k)
        if(n.ge.1) then
           RP=spread(RHO__Minus(:,:,n),DIM=3,NCOPIES=EdgesTotal)        !%[3 9 EdgesTotal]
           sumRPpRHO_P(1,:,:)=sum(RP*RHO_P,DIM=1)!*Area(n)
           sumRPpRHO_M(1,:,:)=sum(RP*RHO_M,DIM=1)!*Area(n)

           A=sum(gP*sumRPpRHO_P,DIM=2)+sum(gM*sumRPpRHO_M,DIM=2)
           A2=sum(gP2*sumRPpRHO_P,DIM=2)+sum(gM2*sumRPpRHO_M,DIM=2)
           Z1=transpose(FactorA*A)
           Z1_2=transpose(FactorA*A+FactorA*A2)
           Z12_2=transpose(FactorA2*A+FactorA2_2*A2)

          
           Matimp(:,n)=Matimp(:,n)+temp_EdgeLength(n,1)*(Z1(:,1)-ZF(:,1))
           Matimp_2(:,n)=Matimp_2(:,n)+temp_EdgeLength(n,1)*(Z1_2(:,1)-ZF_2(:,1))
           Matimp2_2(:,n)=Matimp_2(:,n)+temp_EdgeLength(n,1)*(Z12_2(:,1)-ZF2_2(:,1)) 
        endif
           
     enddo

 !***********************************************************************************************
 !Calculation of Impedance submatrix D-mn for 'dec' case

     DP=D(:,:,TrianglePlus)                     !%[3 9 EdgesTotal]
     DM=D(:,:,TriangleMinus)                    !%[3 9 EdgesTotal]

     do l=1,size(Plus)
        n=Plus(l)
        
         RP=spread(RHO__Plus(:,:,n),DIM=3,NCOPIES=EdgesTotal)           !%[3 9 EdgesTotal]

         do i=1,EdgesTotal
            do j=1,9
               do k=1,3
               if(k.eq.1) then
               RHOcD__P(k,j,i)=RP(2,j,i)*DP(3,j,i)-RP(3,j,i)*DP(2,j,i)
               end if
               if(k.eq.2) then
               RHOcD__P(k,j,i)=RP(3,j,i)*DP(1,j,i)-RP(1,j,i)*DP(3,j,i)          !plus cross product
               endif
               if(k.eq.3) then
               RHOcD__P(k,j,i)=RP(1,j,i)*DP(2,j,i)-RP(2,j,i)*DP(1,j,i)
               endif     
               enddo
            enddo
         enddo

         do i=1,EdgesTotal
            do j=1,9
               do k=1,3
               if(k.eq.1) then
               RHOcD__M(k,j,i)=RP(2,j,i)*DM(3,j,i)-RP(3,j,i)*DM(2,j,i)
               end if
               if(k.eq.2) then
               RHOcD__M(k,j,i)=RP(3,j,i)*DM(1,j,i)-RP(1,j,i)*DM(3,j,i)          !minus cross product
               end if
               if(k.eq.3) then
               RHOcD__M(k,j,i)=RP(1,j,i)*DM(2,j,i)-RP(2,j,i)*DM(1,j,i)
               end if  
               enddo
            enddo
         enddo



             sumRdPpRHOcD__P(1,:,:)=sum(RHOcD__P*RHO_P,DIM=1)!*Area(n)
             sumRdPpRHOcD__M(1,:,:)=sum(RHOcD__M*RHO_M,DIM=1)!*Area(n)
        
 
             Dmat=sum(gdP*sumRdPpRHOcD__P,DIM=2)+sum(gdM*sumRdPpRHOcD__M,DIM=2)
             Dmat2=sum(gdP2*sumRdPpRHOcD__P,DIM=2)+sum(gdM2*sumRdPpRHOcD__M,DIM=2)
             Zd1=transpose(Factord*Dmat+Factord*Dmat2)

             
             DMatsub1(:,n)=DMatsub1(:,n)+temp_EdgeLength(n,1)*Zd1(:,1)

     
     enddo

     do l=1,size(Minus)
        
        n=Minus(l)
        if(n.ge.1) then


         RP=spread(RHO__Minus(:,:,n),DIM=3,NCOPIES=EdgesTotal)         !%[3 9 EdgesTotal]

         do i=1,EdgesTotal
            do j=1,9
               do k=1,3
               if(k.eq.1) then
               RHOcD__P(k,j,i)=RP(2,j,i)*DP(3,j,i)-RP(3,j,i)*DP(2,j,i)
               end if
               if(k.eq.2) then
               RHOcD__P(k,j,i)=RP(3,j,i)*DP(1,j,i)-RP(1,j,i)*DP(3,j,i)          !plus cross product
               end if
               if(k.eq.3) then  
               RHOcD__P(k,j,i)=RP(1,j,i)*DP(2,j,i)-RP(2,j,i)*DP(1,j,i)
               endif     
               enddo
            enddo
         enddo

         do i=1,EdgesTotal
            do j=1,9
               do k=1,3
               if(k.eq.1) then
               RHOcD__M(k,j,i)=RP(2,j,i)*DM(3,j,i)-RP(3,j,i)*DM(2,j,i)
               end if
               if(k.eq.2) then
               RHOcD__M(k,j,i)=RP(3,j,i)*DM(1,j,i)-RP(1,j,i)*DM(3,j,i)          !minus cross product
               end if
               if(k.eq.3) then
               RHOcD__M(k,j,i)=RP(1,j,i)*DM(2,j,i)-RP(2,j,i)*DM(1,j,i)
               end if  
               enddo
            enddo
         enddo



             sumRdPpRHOcD__P(1,:,:)=sum(RHOcD__P*RHO_P,DIM=1)!*Area(n)
             sumRdPpRHOcD__M(1,:,:)=sum(RHOcD__M*RHO_M,DIM=1)!*Area(n)


             Dmat=sum(gdP*sumRdPpRHOcD__P,DIM=2)+sum(gdM*sumRdPpRHOcD__M,DIM=2)
             Dmat2=sum(gdP2*sumRdPpRHOcD__P,DIM=2)+sum(gdM2*sumRdPpRHOcD__M,DIM=2)
             Zd1=transpose(Factord*Dmat+Factord*Dmat2)

                
             DMatsub1(:,n)=DMatsub1(:,n)+temp_EdgeLength(n,1)*Zd1(:,1)
                
        endif
           
     enddo

     DMatsub2=-1.0d+00*DMatsub1

   

     deallocate(Plus)
     deallocate(Minus)

  end do



 !Impedance matrix fill loops for 'dec' case


     do i=1,EdgesTotal
        do j=1,EdgesTotal
        dMatimp(i,j)=Matimp_2(i,j)
        end do
     end do

     do i=1,EdgesTotal
        do j=EdgesTotal+1,2*EdgesTotal    
        dMatimp(i,j)=DMatsub1(i,j-EdgesTotal)
        end do
     end do

     do i=EdgesTotal+1,2*EdgesTotal 
        do j=1,EdgesTotal
        dMatimp(i,j)=DMatsub2(i-EdgesTotal,j)
        end do
     end do

     do i=EdgesTotal+1,2*EdgesTotal
        do j=EdgesTotal+1,2*EdgesTotal    
        dMatimp(i,j)=Matimp2_2(i-EdgesTotal,j-EdgesTotal)
        end do
     end do  


end subroutine impmet


!subroutine computes scattered nearfield in a window of area (N*lamda/10)^2
!where N is the number of points sampled along one of the axes, so that approximately N^2 points are sampled in total
!Refer to "Surface integral formulation for 3D simulations of plasmonic and high permittivity nanostructures eq.35
subroutine nearfield( nField,EdgesTotal,TrianglesTotal,EdgeLength,Area,Kc,Center_,TrianglePlus,TriangleMinus,&
     & RHO_P,RHO_M,RHO__Plus,RHO__Minus,Vd,t,element_order,f,epsilon_,mu_,c_,omega,k,pi,N_pts)   
  implicit none


!*******************************************************************************
!Inputs
  integer ( kind = 4 ) N_pts 
  integer ( kind = 4 ) EdgesTotal                                  
  integer ( kind = 4 ) TrianglesTotal   
  integer ( kind = 4 ), dimension (EdgesTotal) :: TrianglePlus
  integer ( kind = 4 ), dimension (EdgesTotal) :: TriangleMinus 
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ), dimension (element_order, TrianglesTotal) :: t
  integer ( kind = 4) i,j,p,n
  
  real    ( kind = 8 ),dimension(EdgesTotal):: Edgelength
  real    (kind = 8), dimension(TrianglesTotal)::Area
  real    ( kind = 8 ),dimension (3,9,TrianglesTotal)::Center_  
  
  real     ( kind = 8 ), dimension(3,9,EdgesTotal):: RHO_P
  real     ( kind = 8 ), dimension(3,9,EdgesTotal):: RHO_M  
  
  real    ( kind = 8 ),dimension(3,9,EdgesTotal):: RHO__Plus
  real    ( kind = 8 ),dimension (3,9,EdgesTotal)::RHO__Minus

  real    ( kind = 8 )f,epsilon_,mu_,c_,omega,k,pi
  
  complex  ( kind = 8) Kc
  complex  ( kind = 8 ), dimension(2*EdgesTotal)::Vd

!Local
  integer ( kind = 4 ) ij_upper,ij_lower
  
  real ( kind = 8 ),allocatable, dimension(:,:)::fieldpt_xyz 
  real ( kind = 8 ),allocatable, dimension(:,:)::singularpts_xyz
  real ( kind = 8 ),allocatable, dimension(:,:)::allowedpts_xyz

  real ( kind = 8 ),dimension(3,9,TrianglesTotal)::D
  real ( kind = 8 ),dimension(3,9,EdgesTotal)::D_P,D_M 
  real ( kind = 8 ),dimension(1,9,TrianglesTotal)::R 

  real ( kind = 8 ),dimension(3,9,TrianglesTotal)::g
  real ( kind = 8 ),dimension(3,9,EdgesTotal)::g_P,g_M
  real ( kind = 8 ),dimension(3,9,TrianglesTotal)::grad_g
  real ( kind = 8 ),dimension(3,9,EdgesTotal)::grad_g_P,grad_g_M  
  real ( kind = 8 ),dimension(3,9,EdgesTotal)::RHOcD__P,RHOcD__M

  real ( kind = 8 ),dimension(3,1,EdgesTotal)::Efield_1P,Efield_1M
  real ( kind = 8 ),dimension(3,1,EdgesTotal)::Efield_2P,Efield_2M
  real ( kind = 8 ),dimension(3,1,EdgesTotal)::Efield_3P,Efield_3M
  real ( kind = 8 ) lamda

  complex ( kind = 8 ),parameter::im1=(0.0d+00,1.0d+00)
  complex ( kind = 8 ),dimension(3,1,EdgesTotal)::nField_int

!Output
  complex ( kind = 8 ), dimension(3,1,N_pts**2):: nField 

!*********************************************************************************

  !ij_upper=CEILING(5.0d+00*((N_pts-1.0d+00)/N_pts)*(N_pts/10.0d+00+1.0d+00))
  !ij_lower=FLOOR(5.0d+00*((N_pts-1.0d+00)/N_pts)*(N_pts/10.0d+00-1.0d+00))

allocate(fieldpt_xyz(3,N_pts**2))
allocate(allowedpts_xyz(3,N_pts**2))


lamda=c_/f
!loops fill field point array with sampling plane which is here arbitrarily chosen to be the xz-plane

 p=1
   do i=0,N_pts-1
      do j=0,N_pts-1
         fieldpt_xyz(1,p) = -(N_pts/2.0d+00)*(lamda/10.0d+00)+ 2.0d+00*i*(N_pts/(N_pts-1))*(1/2.0d+00)*(lamda/10.0d+00)
         fieldpt_xyz(2,p) = 0.0d+00
         fieldpt_xyz(3,p) = (N_pts/2.0d+00)*(lamda/10.0d+00)- 2.0d+00*j*(N_pts/(N_pts-1))*(1/2.0d+00)*(lamda/10.0d+00)

         !eliminates points inside or on surface of sphere of radius lamda/2 (adapt for arbitrary geometry)
         if(sqrt((fieldpt_xyz(1,p))**2+(fieldpt_xyz(2,p))**2+(fieldpt_xyz(3,p))**2).gt.(lamda/2.0d+00)) then 
         allowedpts_xyz(1,p)=fieldpt_xyz(1,p)
         allowedpts_xyz(2,p)=fieldpt_xyz(2,p)
         allowedpts_xyz(3,p)=fieldpt_xyz(3,p)
         else 
         allowedpts_xyz(1,p)=0.0d+00
         allowedpts_xyz(2,p)=0.0d+00
         allowedpts_xyz(3,p)=0.0d+00         
         end if

         p=p+1
      end do
   end do



 !scattered E-field integration loop
 do p=1,N_pts**2

 if((allowedpts_xyz(1,p)==0.0d+00).and.(allowedpts_xyz(2,p)==0.0d+00).and.(allowedpts_xyz(3,p)==0.0d+00)) then
    nField(:,1,p)=0.0d+00
 else 
    D=spread(spread(allowedpts_xyz(:,p),DIM=2,NCOPIES=9),DIM=3,NCOPIES=TrianglesTotal)-Center_
    D_P=D(:,:,TrianglePlus)
    D_M=D(:,:,TriangleMinus)
    R(1,:,:)=sqrt(sum(D**2,DIM=1))

    g=spread(exp(-im1*k*R(1,:,:))/(4*pi*R(1,:,:)),DIM=1,NCOPIES=3)
    g_P=g(:,:,TrianglePlus)
    g_M=g(:,:,TriangleMinus)

    grad_g=spread((1+im1*k*R(1,:,:))*exp(-im1*k*R(1,:,:))/(4*pi*R(1,:,:)**3),DIM=1,NCOPIES=3)
    grad_g_P=grad_g(:,:,TrianglePlus)
    grad_g_M=grad_g(:,:,TriangleMinus)
    

    Efield_1P(3,:,:)=-1.0d+00/(9.0d+00*k**2)*sum(grad_g_P*D_P,DIM=2)
    Efield_1M(3,:,:)=1.0d+00/(9.0d+00*k**2)*sum(grad_g_M*D_M,DIM=2)
    Efield_2P(3,:,:)=1.0d+00/18.0d+00*sum(g_P*RHO__Plus,DIM=2)
    Efield_2M(3,:,:)=1.0d+00/18.0d+00*sum(g_M*RHO__Minus,DIM=2)


    
!*****************************************************************************************************************
         do i=1,EdgesTotal
            do j=1,9
               do n=1,3
               if(n.eq.1) then
               RHOcD__P(n,j,i)=RHO__Plus(2,j,i)*D_P(3,j,i)-RHO__Plus(3,j,i)*D_P(2,j,i)
               end if
               if(n.eq.2) then
               RHOcD__P(n,j,i)=RHO__Plus(3,j,i)*D_P(1,j,i)-RHO__Plus(1,j,i)*D_P(3,j,i)          !plus cross product
               end if
               if(n.eq.3) then  
               RHOcD__P(n,j,i)=RHO__Plus(1,j,i)*D_P(2,j,i)-RHO__Plus(2,j,i)*D_P(1,j,i)
               endif     
               enddo
            enddo
         enddo

         do i=1,EdgesTotal
            do j=1,9
               do n=1,3
               if(n.eq.1) then
               RHOcD__M(n,j,i)=RHO__Minus(2,j,i)*D_M(3,j,i)-RHO__Minus(3,j,i)*D_M(2,j,i)
               end if
               if(n.eq.2) then
               RHOcD__M(n,j,i)=RHO__Minus(3,j,i)*D_M(1,j,i)-RHO__Minus(1,j,i)*D_M(3,j,i)         !minus cross product
               end if
               if(n.eq.3) then
               RHOcD__M(n,j,i)=RHO__Minus(1,j,i)*D_M(2,j,i)-RHO__Minus(2,j,i)*D_M(1,j,i)
               end if  
               enddo
            enddo
         enddo

!*****************************************************************************************************************

    Efield_3P(3,:,:)=1.0d+00/18.0d+00*sum(grad_g_P*RHOcD__P,DIM=2)
    Efield_3M(3,:,:)=1.0d+00/18.0d+00*sum(grad_g_M*RHOcD__M,DIM=2)

      do i=1,EdgesTotal
         nField_int(:,1,i)=-1.0d+00*Vd(i)*Edgelength(i)*omega*mu_*(Efield_1P(:,1,i)+Efield_1M(:,1,i)&
      & +Efield_2P(:,1,i)+Efield_2M(:,1,i))+Vd(i+EdgesTotal)*Edgelength(i)*(Efield_3P(:,1,i)+Efield_3M(:,1,i))
         
      end do    
          
    
    nField(:,:,p)=sum(nField_int,DIM=3)  

 end if   
 end do  


deallocate(fieldpt_xyz)
deallocate(allowedpts_xyz)

end subroutine nearfield


!Subroutine computes the scattered farfield for given values of theta and phi. Refer to Gedney's notes pg 39
subroutine farfield(fField,EdgesTotal,TrianglesTotal,EdgeLength,Area,Kc,Center_,TrianglePlus,TriangleMinus,&
     & RHO_P,RHO_M,RHO__Plus,RHO__Minus,Vd,t,element_order,f,epsilon_,mu_,c_,omega,k,pi)

 implicit none 

!*******************************************************************************
!Inputs
  integer ( kind = 4 ) EdgesTotal                                  
  integer ( kind = 4 ) TrianglesTotal   
  integer ( kind = 4 ), dimension (EdgesTotal) :: TrianglePlus
  integer ( kind = 4 ), dimension (EdgesTotal) :: TriangleMinus 
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ), dimension (element_order, TrianglesTotal) :: t
  
  real ( kind = 8 ),dimension(EdgesTotal):: Edgelength
  real (kind = 8), dimension(TrianglesTotal)::Area
  real (kind = 8),dimension (3,9,TrianglesTotal)::Center_  
  
  real  (kind = 8), dimension(3,9,EdgesTotal):: RHO_P
  real  (kind = 8), dimension(3,9,EdgesTotal):: RHO_M  
  
  real    ( kind = 8 ),dimension(3,9,EdgesTotal):: RHO__Plus
  real    ( kind = 8 ),dimension (3,9,EdgesTotal)::RHO__Minus

  real    ( kind = 8 )f,epsilon_,mu_,c_,omega,k,pi
  
  complex  ( kind = 8) Kc
  complex  ( kind = 8 ), dimension(2*EdgesTotal)::Vd

!Output
 complex (kind = 8),dimension(3) :: fField

!Local
  integer ( kind = 4) i,j,p,n
  
  real (kind = 8),parameter :: theta = 1.570d+00
  real (kind = 8),parameter :: phi = 0.0d+00
  real (kind = 8), dimension(3,9,EdgesTotal) :: Center_Plus, Center_Minus
  real (kind = 8), dimension(EdgesTotal) :: N_theta, N_phi                  !refer to "Antenna Theory" by C.A.Balanis pg660 3rd ed.
  real (kind = 8),dimension(1,9,Edgestotal) :: N_thP,N_thM,N_phP,N_phM
  real (kind = 8),dimension(1,9,Edgestotal) :: exp_factP,exp_factM,rP,rM

  complex (kind = 8) Etheta, Ephi

!********************************************************************************

    Center_Plus(:,:,:) = Center_(:,:,TrianglePlus)
    Center_Minus(:,:,:) = Center_(:,:,TriangleMinus)

    rP(1,:,:) = Center_Plus(1,:,:)*sin(theta)*cos(phi) + Center_Plus(2,:,:)*sin(theta)*sin(phi) + Center_Plus(3,:,:)*cos(theta)
    rM(1,:,:) = Center_Minus(1,:,:)*sin(theta)*cos(phi) + Center_Minus(2,:,:)*sin(theta)*sin(phi) + Center_Minus(3,:,:)*cos(theta)

    exp_factP = exp(Kc*rP)
    exp_factM = exp(Kc*rM)
  
    N_thP(1,:,:) = RHO__Plus(1,:,:)*cos(theta)*cos(phi) + RHO__Plus(2,:,:)*cos(theta)*sin(phi) - RHO__Plus(3,:,:)*sin(phi)
    N_thM(1,:,:) = RHO__Minus(1,:,:)*cos(theta)*cos(phi) + RHO__Minus(2,:,:)*cos(theta)*sin(phi) - RHO__Minus(3,:,:)*sin(phi)
    N_phP(1,:,:) = -1*RHO__Plus(1,:,:)*sin(phi) + RHO__Plus(2,:,:)*cos(phi) 
    N_phM(1,:,:) = -1*RHO__Minus(1,:,:)*sin(phi) + RHO__Minus(2,:,:)*cos(phi) 

    N_theta = sum(N_thP(1,:,:)*exp_factP(1,:,:) - N_thM(1,:,:)*exp_factM(1,:,:),DIM=2)           
    N_phi = sum(N_phP(1,:,:)*exp_factP(1,:,:) - N_phM(1,:,:)*exp_factM(1,:,:),DIM=2)    

  Etheta = 0.0d+00
  Ephi = 0.0d+00
    do i = 1,EdgesTotal 
  
    Etheta = Etheta + (-1.00d+00/18.00d+00)*Kc*sqrt(mu_/epsilon_)*Edgelength(i)*Vd(i)*N_theta(i)
    Ephi = Ephi + (-1.00d+00/18.00d+00)*Kc*sqrt(mu_/epsilon_)*Edgelength(i)*Vd(i)*N_phi(i)
 
   end do


   fField(1) = Etheta*cos(theta)*cos(phi)-Ephi*sin(phi)https://ch1prd0610.outlook.com/owa/?hm=1&wa=wsignin1.0
   fField(2) = Etheta*cos(theta)*sin(phi)+Ephi*cos(phi)
   fField(3) = -1.0d+00*Etheta*sin(theta)


end subroutine farfield


