Module Pole2LFunctions

Use Control
Use Mathematics
Use LoopFunctions
Use Settings

Real(dp), parameter :: epsD = 1E-8_dp, epsM = 1E-10_dp


 Interface Smalldp
  Module Procedure SmalldpR, SmalldpC
 End Interface

Contains

!!!!!!!!!!! TADPOLE FUNCTIONS


real(dp) Function TfSS(x,y,z,Q2)
Implicit None
real(dp), intent(in):: x,y,z,q2
real(dp) :: res

if (smallc(z/q2)) then
   res = 0._dp
else
  if (NewGBC) Then
   if(allsmall(q2,x,y)) Then
      res = 0._dp
   else
      res = real(-J0(z,Q2)*BB(x,y,q2),dp)
   end if
  Else
       res = real(-J0(z,Q2)*BB(x,y,q2),dp)
  End if

   
end if

if(res .ne. res) then
      write(*,*) "Problem in TfSS"
   end if


TfSS=res
end function tfss


real(dp) Function TfSSS(x,y,z,Q2)
Implicit None
real(dp), intent(in):: x,y,z,q2
real(dp) :: res

res=-real(II(x,y,z,q2),dp)

if(res .ne. res) then
      write(*,*) "Problem in TfSSS"
   end if



TfSSS=res
end function tfsss

real(dp) Function TfSSSS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

If (NewGBC) Then
 if(allsmall(q2,x,y)) then
   if(allsmall(q2,z,u)) then
      ! GGGG -> 0
      res=0._dp
   else
      ! R_SS(0,0,z,u,q2)
      res=real(RSS(z,u,q2),dp)
   end if
 else
   ! Just U_0
      If(smalldp(x,y)) then
         res=real(-1._dp*IIp(x,z,u,q2),dp)
      else
         res=Real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
      end if
 end if
Else 
 If(smalldp(x,y)) then
      res=real(-1._dp*IIp(x,z,u,q2),dp)
 else
   res=Real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
 end if
End if

if(res .ne. res) then
      write(*,*) "Problem in TfSSSS"
   end if

TfSSSS=res
end function tfssss



real(dp) Function TfSSFF(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

If (NewGBC) Then
 if(allsmall(q2,x,y)) then
   res=real(-II(z,u,0._dp,q2)-(z+u)*RSS(z,u,q2),dp)
 else
   res=real(-(z+u-y)*UU(x,y,z,u,q2)-II(x,z,u,q2)+BB(x,y,q2)*(J0(z,q2)+J0(u,q2)),dp)
 end if
Else
 res=real(-(z+u-y)*UU(x,y,z,u,q2)-II(x,z,u,q2)+BB(x,y,q2)*(J0(z,q2)+J0(u,q2)),dp)
End If

if(res .ne. res) then
      write(*,*) "Problem in TfSSFF"
   end if

TfSSFF=res
end function tfssff

real(dp) Function TfSSFbFb(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

If (NewGBC) Then
 if(allsmall(q2,x,y)) then
   res=real(RSS(z,u,q2),dp)
 else
   RES=real(UU(x,y,z,u,q2),dp) !!!!!!! NB THIS IS opposite of definition in paper !! April 17
 end if
Else
 RES=real(UU(x,y,z,u,q2),dp)
End if

if(res .ne. res) then
      write(*,*) "Problem in TfSSFbFb"
   end if


TfSSFbFb=res
end function tfssfbfb

real(dp) Function TfFFFbS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

RES=real(II(x,z,u,q2)-(x+y)/2._dp*UU(x,y,z,u,q2),dp)

if(res .ne. res) then
      write(*,*) "Problem in TfFFFbS"
   end if



TfFFFbS=res
end function TfFFFbS

real(dp) Function TfFFbFS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

RES=Fffs100(x,y,z,u,q2)

if(res .ne. res) then
      write(*,*) "Problem in TfFFbFS"
   end if



TfFFbFS=res
end function TfFFbFS

real(dp) Function TfFbFbFbS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

RES=real(UU(x,y,z,u,q2),dp)!!!! opposite sign definition to paper!


if(res .ne. res) then
      write(*,*) "Problem in TfFbFbFbS"
   end if

TfFbFbFbS=res
end function tffbfbfbs

real(dp) Function TfSV(x,Q2)
Implicit None
real(dp), intent(in):: x,q2
real(dp) :: res

if(smallc(x/q2)) then
   res = 0._dp
else
   res=real(3._dp*x*(-5._dp+4._dp*log(x/q2)-(log(x/q2))**2)-J0(x,q2)+2._dp*x,dp)
end if


if(res .ne. res) then
      write(*,*) "Problem in TfSV"
   end if

TfSV=res
end function tfsv

real(dp) Function TfFV(x,Q2)
Implicit None
real(dp), intent(in):: x,q2
real(dp) :: res


if(smallc(x/q2)) then
res=0._dp
tffv=res
return
else
   res=real(4._dp*x*(6._dp-7._dp*log(x/q2)+3._dp*(log(x/q2))**2),dp)
end if

if(rMS.ne.0._dp) then
   res=res+real(4._dp*x*(2._dp*log(x/q2)-1._dp),dp)
end if

if(res .ne. res) then
      write(*,*) "Problem in TfFV"
   end if
   

TfFV=res
end function tffv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MASS FUNCTIONS


real(dp) Function WfSSSS (p2,x_in,y_in,z_in,u_in,Q2)
Implicit None
real(dp), intent(in):: p2,x_in,y_in,z_in,u_in,q2
real(dp) :: res,x,y,z,u

x=x_in
y=y_in
z=z_in
u=u_in

If (NewGBC) Then 
 if(allsmall(q2,x,y)) then
   res=real(RSS(z,u,q2),dp)
 else
   If(smalldp(x,y)) then
      res=real(-1._dp*IIp(x,z,u,q2),dp)
   else
      res=Real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
   end if
 end if
Else
 If(smalldp(x,y)) then
   res=real(-1._dp*IIp(x,z,u,q2),dp)
 else
   res=Real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
 end if
End if

if(res .ne. res) then
      write(*,*) "Problem in WfSSSS"
   end if
WfSSSS=res

end function WfSSSS


real(dp) Function XfSSS(p2,x,y,z,Q2)

Implicit None
real(dp), intent(in):: p2,x,y,z,q2
real(dp) :: res

if (smallc(z/q2)) then
   res = 0._dp
else
 If (NewGBC) Then
   if(allsmall(q2,x,y)) then
      res=0._dp
   else      
      res = real(-J0(z,Q2)*BB(x,y,q2),dp)
   end if
 Else
  res = real(-J0(z,Q2)*BB(x,y,q2),dp)
 End if
end if

if(res .ne. res) then
      write(*,*) "Problem in XfSSS"
   end if

XfSSS = res
end function XfSSS


real(dp) Function YfSSSS(p2,x,y,z,u,Q2)
real(dp), intent(in):: p2,x,y,z,u,q2

real(dp) :: res


if(smallc(u/q2)) then
   res = 0._dp
else
 If (NewGBC) Then
   if(allsmall(q2,y,z)) then
      res = 0._dp
   else
      res = real(J0(u,q2)*CC(x,y,z,q2),dp)
   end if
 Else
    res = real(J0(u,q2)*CC(x,y,z,q2),dp)
 End if
end if

if(res .ne. res) then
      write(*,*) "Problem in YfSSSS", p2,x,y,z,u,q2
   end if

YfSSSS=res


end Function YfSSSS


real(dp) Function ZfSSSS(p2,x,y,z,u,Q2)
real(dp), intent(in):: p2,x,y,z,u,q2

real(dp) :: res


If (NewGBC) Then
 res=real(BBs(p2,x,y,q2)*BBs(p2,z,u,q2),dp)
Else
 res = real(BB(x,y,q2)*BB(z,u,q2),dp)
End if
ZfSSSS=res
end function ZfSSSS

real(dp) Function SfSSS(p2,x,y,z,Q2)
real(dp), intent(in):: p2,x,y,z,q2

real(dp) :: res


res = -real(II(x,y,z,q2),dp)

if(res .ne. res) then
      write(*,*) "Problem in SfSSS"
   end if

SfSSS=res
end function SfSSS

real(dp) Function UfSSSS(p2,x,y,z,u,Q2)
real(dp), intent(in):: p2,x,y,z,u,q2

real(dp) :: res

If (NewGBC) Then 
 res = real(UUs(p2,x,y,z,u,q2),dp)
Else
 res = UU(x,y,z,u,q2)
End if

if(res .ne. res) then
      write(*,*) "Problem in UfSSSS"
   end if

UfSSSS=res
end function UfSSSS








real(dp) Function VfSSSSS(p2,x_in,y_in,z_in,u_in,v_in,Q2)
real(dp), intent(in):: p2,x_in,y_in,z_in,u_in,v_in,q2

real(dp) :: res,u,v,x,y,z

if(abs(v_in) .gt. abs(u_in)) then
   v = v_in
   u = u_in
else
   v = u_in
   u = v_in
end if

x=x_in
y=y_in
z=z_in

If (NewGBC) Then
if(allsmall(q2,y,z)) then
   !V(x,0,0,u,v)
   if(allsmall(q2,u,v)) then
      !V(x,0,0,0,0) -> no coupling to this
      res=0._dp
   else
      if(smallc(x/q2)) then
         !V(0,0,0,u,v)
!         res=real(-BB(u,v,q2)*(1._dp-2._dp*clog(-p2/q2))/p2 +RSS(u,v,q2)/p2,dp)
         res=0._dp
         !! Need case for u=0
         if(smallc(u/q2)) then
            !V(0,0,0,0,v)
            res=res+real( clog(-p2/q2)/(-2._dp*v) +(5._dp+2._dp*clog(v/q2))/(4._dp*v),dp)
         else
            if(smalldp(u,v)) then
               res=res+real(clog(-p2/q2)/(-6._dp*u) +(11._dp+3._dp*clog(u/q2))/(18._dp*u),dp)
            else             
               res=res+real( (-clog(-p2/q2)/(2._dp*(v-u)**3))*(v**2-u**2-2*u*v*clog(v/u)),dp)
               ! res=res+real((0.25_dp/((v-u)**3))*(5._dp*(v-u)**2 + 2._dp*u*clog(u/q2)*(2._dp*v-u+2._dp*v*clog(v/q2)) +2._dp*v*clog(v/q2)*(v-2._dp*u-2._dp*u*clog(v/q2))-8._dp*u*v*CLi2_mod(cmplx(1._dp-u/v,0._dp,dp))),dp)
               !res=res+real((0.25_dp/((v-u)**3))*(5._dp*(v**2-u**2) + 2._dp*u*clog(u/q2)*(2._dp*v-u+2._dp*v*clog(u/q2)) +2._dp*v*clog(v/q2)*(v-2._dp*u-2._dp*u*clog(v/q2))-4._dp*u*v*(CLi2_mod(cmplx(1._dp-u/v,0._dp,dp)) - CLi2_mod(cmplx(1._dp-v/u,0._dp,dp)) )),dp)

               res=res+real((0.25_dp/((v-u)**3))*(5._dp*(v**2-u**2) + 2._dp*u*clog(u/q2)*(2._dp*v-u+2._dp*v*clog(v/q2)) +2._dp*v*clog(v/q2)*(v-2._dp*u-2._dp*u*clog(v/q2))-8._dp*u*v*CLi2_mod(cmplx(1._dp-u/v,0._dp,dp))),dp)
               

            end if
         endif
               
      else
         !! no momentum dependence necessary
         !res=real(UU(x,0._dp,u,v,q2)/x**2+(RSS(u,v,q2)+BB(u,v,q2)*(clog(x/q2)-1._dp))/x,dp)
         !res=real(UU(x,0._dp,u,v,q2)/x+(RSS(u,v,q2)+BB(u,v,q2)*(clog(x/q2)-1._dp))/x,dp)
         res=real(UU(x,0._dp,u,v,q2)/x+(-RSS(u,v,q2)+BB(u,v,q2)*(clog(x/q2)-1._dp))/x,dp)
 
      end if
   end if

   
 if(res .ne. res) then
      write(*,*) "Problem in special VfSSSSS",x,y,z,u,v,p2
   end if
   
else
   if(smalldp(x,y)) then
      res=real(-VV0(z,x,u,v,q2),dp)
   else
      if(smalldp(y,z)) then
         res=real(-VV0(x,y,u,v,q2),dp)
      else
         res=real((UU(x,y,u,v,q2)-UU(x,z,u,v,q2))/(y-z),dp)
      end if
   end if

   if(res .ne. res) then
      write(*,*) "Problem in regular VfSSSSS"
   end if
   
end if
Else 

if(smalldp(x,y)) then
       res=real(-VV0(z,x,u,v,q2),dp)
else
   if(smalldp(y,z)) then
      res=real(-VV0(x,y,u,v,q2),dp)
   else
      res=real((UU(x,y,u,v,q2)-UU(x,z,u,v,q2))/(y-z),dp)
   end if
end if

End if

! if(res .ne. res) then
!       write(*,*) "Problem in VfSSSSS"
!    end if
   
VfSSSSS=res
end function VfSSSSS

real(dp) Function MfSSSSS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=0._dp

If (NewGBC) Then
if(allsmall(q2,z,u,v)) then
   res=0._dp
else
   if(allsmall(q2,x,y,v)) then
      res=0._dp
   else
      if(allsmall(q2,x,z)) then
         if(allsmall(q2,y,u)) then
            ! M(0,0,0,0,v)
            res=-real((clog(-p2/v)**2-2._dp*clog(-p2/v)+Pi**2/3._dp)/v,dp)
         else
            ! -M(0,y,0,u,v)
            res=-real(MM00(p2,y,u,v,q2),dp)
         end if
      else
         if(allsmall(q2,y,u)) then
            ! -M(x,0,z,0,v)
            res=-real(MM00(p2,x,z,v,q2),dp)
         else
            !! FINITE!
            res=-real(MM(x,y,z,u,v,q2),dp)
         end if
      end if
   end if
end if
Else
 res=-real(MM(x,y,z,u,v,q2),dp)
End if
if(res .ne. res) then
      write(*,*) "Problem in MfSSSSS"
   end if

MfSSSSS=res
end function MfSSSSS


complex(dp) function AM(yin,uin,vin,q2)
  
  real(dp), intent(in) :: yin,uin,vin,q2
  real(dp) :: y,u,v
  complex(dp) :: res
  real(dp) :: invals(3)
  
  ! if( abs(yin) .le. abs(uin)) then
  !    if( abs(uin) .le. abs(vin)) then
  !       call atvals(y,u,v,x,y,yin,uin,vin)

  !    else
  !       if(
  !    end if
  ! else
  ! end if
  invals = (/ yin, uin,vin /)
  call sortabsvals(invals)
  y=invals(1)
  u=invals(2)
  v=invals(3)


  
  if(smallc(y/q2)) then
     !!! Are assuming only one value can be zero
     if(smalldp(u,v)) then
        res=-1._dp/u
        
     else
        res=clog(v/u)/(u-v)
     end if
  else
     if(smalldp(y,u)) then
        if(smalldp(y,v)) then
           !AM(y,y,y)
           res=-0.5_dp/y
        else
           !AM(y,y,v)
           res=v*clog(y/v)/((y-v)**2)-1._dp/(y-v)
        end if
     else
        if(smalldp(u,v)) then
           !AM(y,u,u)
           res=y*clog(u/y)/((u-y)**2)-1._dp/(u-y)
        else
           !AM(y,u,v)
           res=u*clog(u/q2)/((y-u)*(u-v))-y*clog(y/q2)/((y-u)*(y-v))-v*clog(v/q2)/((y-v)*(u-v))
        end if
     end if
  end if
  
  AM=res
end function AM


complex(dp) function MM00(p2,yin,uin,v,q2)
real(dp), intent(in):: p2,yin,uin,v,q2

complex(dp) :: res
complex(dp) :: rAM, BM
real(dp) :: y,u

!!! 
res=(0._dp,0._dp)

if(abs(yin) .le. abs(uin)) then
   y=yin
   u=uin
else
   y=uin
   u=yin
endif
   
rAM=AM(yin,uin,v,q2)

if(smallc(v/q2)) then
   if(smalldp(y,u)) then
      BM=(2._dp+clog(y/q2))/y
   else
      BM=clog(u/y)*(4._dp+clog(u/q2)+clog(y/q2))/(2._dp*(u-y))
   end if
else
   if(smalldp(u,v)) then
      if(smalldp(y,u)) then
!!! BM(y,y,y)
        BM=0.5_dp/y-(2._dp+clog(y/q2))*rAM
      else
         !! BM(y,u,u)
         BM=-(2._dp+clog(u/q2))*rAM+2._dp/(u-y)-(u+y)*CLi2_mod(cmplx(1._dp-y/u,0.,dp))/((u-y)**2)
      end if
   else
      if(smalldp(y,v)) then
         ! BM(y,u,y) = BM(u,y,y)
         BM=-(2._dp+clog(y/q2))*rAM+2._dp/(y-u)-(u+y)*CLi2_mod(cmplx(1._dp-u/y,0.,dp))/((u-y)**2)
      else
         if(smalldp(y,u)) then
            ! BM(y,y,v)
            BM=-(2._dp+clog(v/q2))*rAM+((v+y)*clog(y/v)+2._dp*v*CLi2_mod(cmplx(1._dp-y/v,0.,dp)))/((y-v)**2)
         else
            ! BM(y,u,v)
            BM=-(2._dp+clog(v/q2))*rAM+(u+v)*CLi2_mod(cmplx(1._dp-u/v,0.,dp))/((y-u)*(u-v))-(v+y)*CLi2_mod(cmplx(1._dp-y/v,0.,dp))/((y-u)*(y-v))
         end if
      
      end if
   end if
end if

res=BM+rAM*clog(-p2/q2)

if(res.ne.res) then
   write(*,*) "Problem in MM00!",p2,q2,y,u,v,rAM,BM
   end if

MM00=res

end function MM00


!!!!!!!!!!!!!!!!! SCALARS AND VECTORS


real(dp) Function WfSSSV(p2,x,Q2)
real(dp), intent(in):: p2,x,q2

real(dp) :: res

!res=real(3*II(0._dp,x,x,q2)-J0(x,q2)+2._dp*x,dp)
if(smallc(x/q2)) then
   res = 0._dp
else
   res=real(3._dp*x*(-5._dp+4._dp*log(x/q2)-(log(x/q2))**2)-J0(x,q2)+2._dp*x,dp)
end if


if(res .ne. res) then
      write(*,*) "Problem in WfSSSV"
   end if

WfSSSV=res

end function WfSSSV

real(dp) Function MfSSSSV(p2,x_in,y_in,Q2)
real(dp), intent(in):: p2,x_in,y_in,q2

real(dp) :: res,x,y

if(abs(y_in) .gt. abs(x_in)) then
   y = y_in
   x= x_in
else
   y = x_in
   x =y_in
end if

if(smalldp(x,y)) then
   ! call saveme(q2,y) !!!! y > x so just maybe ..
   if(smallc(x/q2)) then
      res=0._dp
      write(Errcan,*) "MfSSSSV called with small input: ",p2,x_in,y_in,q2
   else
      res = real( -1 + 5*log(y/q2) -3 *log(y/q2)**2,dp)
   end if

else
   if(smallc(x/q2)) then
     !!! NB the case that they are both small is dealt with above ...
      res=real(12 + (11*(y*log(y/q2))-3*(y*(log(y/q2))**2))/(y),dp)
   else
      res=real(-12 + (11*(x*log(x/q2)-y*log(y/q2))-3*(x*(log(x/q2))**2-y*(log(y/q2))**2))/(x-y),dp)
   end if
end if

MfSSSSV=res

end function MfSSSSV


!!!!!!!!!!!!!!!!! FERMIONS AND SCALARS

real(dp) Function WfSSFF (p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: res

If (NewGBC) Then
if(allsmall(q2,x,y)) then
   
      res = real( -II(z,u,0._dp,q2) -(z+u)*RSS(z,u,q2),dp)
   else
      
      res = real( &
           & -(z+u-y)*UU(x,y,z,u,q2)-II(x,z,u,q2)+BB(x,y,q2)*(J0(z,q2)+J0(u,q2)) &
           & ,dp)
  end if
Else
res = real( &
& -(z+u-y)*UU(x,y,z,u,q2)-II(x,z,u,q2)+BB(x,y,q2)*(J0(z,q2)+J0(u,q2)) &
& ,dp)
End if

if(res.ne.res) then
   write(*,*) "Problem in WfSSFF!!"
end if
  
WfSSFF=res

end function WfSSFF


real(dp) Function WfSSFbFb (p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: res

res=real(-2._dp*WfSSSS(p2,x,y,z,u,Q2),dp)

if(res.ne.res) then
   write(*,*) "Problem in WfSSFbFb!!"
   end if

WfSSFbFb=res

end function WfSSFbFb

real(dp) Function MfFbFbFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res



res=2._dp*real(MM(x,y,z,u,v,q2),dp)

MfFbFbFbFbS=res
end function MfFbFbFbFbS


real(dp) Function MfFFbFbFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real( &
& (y+z-v)*MM(x,y,z,u,v,q2)-UU(x,z,u,v,q2) - UU(u,y,x,v,q2) +BB(x,z,q2)*BB(y,u,q2) &
& ,dp)

MfFFbFbFS=res
end function MfFFbFbFS

real(dp) Function MfFFbFFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real( &
& (x+z)*MM(x,y,z,u,v,q2)-UU(y,u,z,v,q2)-UU(u,y,x,v,q2) &
& ,dp)

MfFFbFFbS=res
end function MfFFbFFbS

real(dp) Function MfFFFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real( &
& (x+y-v)*MM(x,y,z,u,v,q2)-UU(x,z,u,v,q2)-UU(y,u,z,v,q2)+BB(x,z,q2)*BB(y,u,q2) &
& ,dp)

MfFFFbFbS=res
end function MfFFFbFbS

real(dp) Function MfFFFFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res



res=real( &
& (x*u+y*z)*MM(x,y,z,u,v,q2)-x*UU(z,x,y,v,q2)-z*UU(x,z,u,v,q2) &
& -u*UU(y,u,z,v,q2)-y*UU(u,y,x,v,q2)+II(x,u,v,q2)+II(y,z,v,q2) &
& ,dp)

MfFFFFS=res
end function MfFFFFS





real(dp) Function MfSFbSFbFb(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

If (NewGBC) Then
if(allsmall(q2,x,z)) then
   res=2._dp*real(MM00(p2,y,u,v,q2),dp)
   else
      res=2._dp*real(MM(x,y,z,u,v,Q2),dp)
   end if
Else
 res=2._dp*real(MM(x,y,z,u,v,Q2),dp)
End if

if(res.ne.res) then
   write(*,*) "Problem in MfSFbSFbFb!!",p2,x,y,z,u,v,Q2
end if
   
MfSFbSFbFb=res
end function MfSFbSFbFb


real(dp) Function MfSFSFbF(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

If (NewGBC) Then
if(allsmall(q2,x,z)) then
   res=real((v+y)*MM00(p2,y,u,v,q2)+UU(y,u,z,v,q2)-UUs(p2,x,z,u,v,q2)-(2._dp-clog(-p2/q2))*BB(y,u,q2),dp)
  else 
     res=real((v-x+y)*MM(x,y,z,u,v,q2)+UU(y,u,z,v,q2)-UU(x,z,u,v,q2)-BB(x,z,q2)*BB(y,u,q2),dp)
end if
Else
res=real((v-x+y)*MM(x,y,z,u,v,q2)+UU(y,u,z,v,q2)-UU(x,z,u,v,q2)-BB(x,z,q2)*BB(y,u,q2),dp)
End if

if(res.ne.res) then
   write(*,*) "Problem in MfSFSFbF!!",p2,x,y,z,u,v,Q2,MM00(p2,y,u,v,q2),UUs(p2,x,z,u,v,q2)
   end if

MfSFSFbF=res
end function MfSFSFbF

real(dp) Function MfSFSFFb(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

If (NewGBC) Then
if(allsmall(q2,x,z)) then
   res=real((y+u)*MM00(p2,y,u,v,q2)-UUs(p2,x,z,u,v,q2)-UUs(p2,z,x,y,v,q2),dp)
   else
res=real((y+u)*MM(x,y,z,u,v,q2)-UU(x,z,u,v,q2)-UU(z,x,y,v,q2),dp)
end if
Else
res=real((y+u)*MM(x,y,z,u,v,q2)-UU(x,z,u,v,q2)-UU(z,x,y,v,q2),dp)
End if

if(res.ne.res) then
   write(*,*) "Problem in MfSFSFFb!!",p2,x,y,z,u,v,Q2
   end if

MfSFSFFb=res
end function MfSFSFFb

real(dp) Function VfSSSFbFb(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = -2._dp*VfSSSSS(p2,x,y,z,u,v,q2)

VfSSSFbFb=res
end function VfSSSFbFb

real(dp) Function VfSSSFF(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

If (NewGBC) Then
if(allsmall(q2,y,z)) then
   res = real(UUs(p2,x,0._dp,u,v,q2) -(u+v)*VfSSSSS(p2,x,0._dp,0._dp,u,v,q2),dp)
else
   res = UU(x,y,u,v,q2) +(z-u-v)*VfSSSSS(p2,x,y,z,u,v,q2)-(J0(u,Q2) +J0(v,Q2))*CC(x,y,z,q2)
end if
Else
res = UU(x,y,u,v,q2) +(z-u-v)*VfSSSSS(p2,x,y,z,u,v,q2)-(J0(u,Q2) +J0(v,Q2))*CC(x,y,z,q2)
End if

VfSSSFF=res
end function VfSSSFF

!!!! For these functions don't want to interfere with the small cases -> but if the fermion masses are "small", so are the prefactors -> set to zero

real(dp) Function VfFbFbFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

!if(smallc(x/q2,y/q2,z/q2,u/q2)) then
if(smallc(y/q2,z/q2)) then !! more conservative, just don't activate GB routines
   res = 0._dp
else
   res = real(-2._dp*VfSSSSS(p2,x,y,z,u,v,q2),dp)
end if


VfFbFbFbFbS=res
end function VfFbFbFbFbS

real(dp) Function VfFbFFbFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

if(smallc(x/q2,z/q2)) then
   res = 0._dp
else
res = real(-UU(x,y,u,v,q2) +(v-z-u)*VfSSSSS(p2,x,y,z,u,v,q2)-(J0(v,Q2) -J0(u,Q2))*CC(x,y,z,q2),dp)
end if 

VfFbFFbFS=res
end function VfFbFFbFS

real(dp) Function VfFbFFFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

!res = real(-2._dp*UU(x,y,u,v,q2) -2._dp*z*VfSSSSS(p2,x,y,z,u,v,q2),dp)
!if(smallc(x/q2,u/q2)) then
if(smallc(x/q2)) then !!! more conservative, just make sure we don't activate the GB routines
   res = 0._dp
else
   if(smallc(z/q2)) then
      res=real(-2._dp*UU(x,y,u,v,q2),dp)
   else
      res = real(-2._dp*UU(x,y,u,v,q2) -2._dp*z*VfSSSSS(p2,x,y,z,u,v,q2),dp)
   end if
end if


VfFbFFFbS=res
end function VfFbFFFbS

real(dp) Function VfFFbFbFS(p2,x_in,y_in,z_in,u_in,v_in,Q2)
real(dp), intent(in):: p2,x_in,y_in,z_in,u_in,v_in,q2

real(dp) :: res,x,y,z,v,u

call atvals(x,y,z,u,v,x_in,y_in,z_in,u_in,v_in)

!call saveme(q2,x,y,z,u)



if(smallc(y/q2,z/q2)) then
   res = 0._dp
else
   res=Fffs200(x,y,z,u,v,q2)
end if


VfFFbFbFS=res
end function VfFFbFbFS

real(dp) Function VfFFFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

!if(smallc(z/q2,u/q2)) then
if(smallc(z/q2)) then 
   res = 0._dp
else
   res = real( -UU(x,y,u,v,q2) - UU(y,z,u,v,q2) -(x+z)*VfSSSSS(p2,x,y,z,u,v,q2),dp)
end if


VfFFFbFbS=res
end function VfFFFbFbS

real(dp) Function VfFFFFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

if(smallc(x/q2)) then
   res = Fffs100(y,z,u,v,q2)
else
!!!! NB Fffs200 is symmetric in y<->x !!!!!!!!!!
!!!!! Nb that since x !=0 we will never activate the GB regulator
res = Fffs100(y,z,u,v,q2)+x*Fffs200(y,x,z,u,v,q2)   
end if



VfFFFFS=res
end function VfFFFFS


!!!!!!!!!!!!!!!!! FERMIONS AND VECTORS


real(dp) Function GfFFV(p2,x_in,y_in,Q2)
real(dp), intent(in) :: p2,x_in,y_in,q2

real(dp) :: res,x,y,z,u,v
call atvals(x,y,z,u,v,x_in,y_in)

if(smalldp(x,y)) then
   if(smallc(x/q2,y/q2)) then
      res = 0._dp 
      gfffv=res
      return 
   else
      res = real(x*(20 - 32*Log(x/Q2) + 24*Log(x/Q2)**2),dp)
   end if
else

res=Real( 2._dp*(  (x+y)*(3._dp*UU(x,y,x,0._dp,q2)+3._dp*UU(x,y,y,0._dp,q2) - 5._dp*BB(x,y,q2))-3._dp*II(x,x,0._dp,q2)-3._dp*II(y,y,0._dp,q2)+5._dp*J0(x,q2)+5._dp*J0(y,q2)-8*(x+y)),dp)


end if

if(rMS.ne.0._dp) then
res=res+real(4._dp*(x+y+J0(x,q2)+J0(y,q2)-(x+y)*BB(x,y,q2)),dp)
end if



GfFFV = res

end function GfFFV

real(dp) Function GfFbFbV(p2,x_in,y_in,Q2)
real(dp), intent(in) :: p2,x_in,y_in,q2

real(dp) :: res,x,y,z,u,v
call atvals(x,y,z,u,v,x_in,y_in)

if(smalldp(x,y)) then
   if(smallc(x/q2,y/q2)) then
      res = 0._dp !!! since the function will be multiplied by zero anyway
      gffbfbv=res
      return
   else
      res = real(-4 - 4*Log(x/Q2) + 12*Log(x/Q2)**2,dp)
   end if
else
   res=real(12._dp*UU(x,y,x,0._dp,q2)+12._dp*UU(x,y,y,0._dp,q2)-20._dp*BB(x,y,q2)-16._dp,dp) 

end if


if(rMS.ne.0._dp) then
res=res-real(4._dp*(-1.0_dp+2._dp*BB(x,y,q2)),dp)
end if


GfFbFbV = res


end function GfFbFbV




!!!!!!!!!!!!!!!!! FUNCTIONS FOR FIRST DERIVATIVES

real(dp) Function Fffs100(x,y,u,v,Q2)
real(dp), intent(in):: x,y,u,v,q2

real(dp) :: res

res = real(BB(x,y,q2)*(J0(v,q2)-J0(u,q2))+II(x,u,v,q2)-(y+u-v)*UU(x,y,u,v,q2),dp)

Fffs100 =res
end function Fffs100

real(dp) Function Fffs001(x,y,u,v,Q2)
real(dp), intent(in):: x,y,u,v,q2

real(dp) :: res

res = real(-(x+y-u)*UU(x,y,u,v,q2)-II(x,y,v,q2)+BB(u,v,q2)*(J0(x,q2)+J0(y,q2)),dp)

Fffs001 =res
end function Fffs001



!!!!!!!!!!!!!!!!! FUNCTIONS FOR SECOND DERIVATIVES *FROM* FIRST DERIVATIVES

real(dp) Function Fffs200(x,y,z,u,v,Q2)
real(dp), intent(in):: x,y,z,u,v,q2

real(dp) :: res

res = real(CC(x,y,z,q2)*(J0(u,q2)-J0(v,q2))-UU(x,z,u,v,q2)-(y+u-v)*VfSSSSS(0._dp,x,y,z,u,v,q2),dp)

Fffs200 =res
end function Fffs200




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HELPER FUNCTIONS NOW





real(dp) function RSS(xin,yin,Q2)
implicit none
real(dp), intent(in) :: xin,yin,Q2
real(dp) :: res,x,y

if(abs(xin) .le. abs(yin)) then
   x=xin
   y=yin
else
   x=yin
   y=xin 
end if


if(smallc(x/q2)) then
   if(smallc(y/q2)) then
      !RSS(0,0) -> return zero but really shouldn't be here
      res = 0._dp
   else
      !RSS(0,y)
      res=0.5_dp-Pi**2/6._dp-0.5_dp*clog(y/q2)**2
   end if
else
   
   if(smalldp(x,y)) then
      ! RSS(x,x)
         !      res = real(-3._dp-2._dp*J0(x,q2)/x-((J0(x,q2))**2)/(2._dp*x**2),dp)
         res = real(-1._dp -0.5_dp*(clog(x/q2)+1._dp)**2,dp)  
      else
         ! RSS(x,y)
      !   res = real(((x+y)**2+2._dp*J0(x,q2)*J0(y,q2)-2._dp*x*J0(x,q2)-2._dp*y*J0(y,q2)+(x+y)*II(x,y,0._dp,q2))/((y-x)**2),dp)
      res=real(0.5_dp +(0.5_dp/(x-y))*(2._dp*y*clog(x/q2)*clog(y/q2)-(x+y)*(clog(x/q2)**2+2._dp*CLi2_mod(cmplx(1._dp-y/x,0.,dp)))),dp)
   end if

end if

if(res .ne. res) then
   write(*,*) "Problem in RSS"
   end if

RSS=res
end function RSS


complex(dp) function BBs(p2,m1_in,m2_in,q2)
implicit none
real(dp), intent(in) :: p2,m1_in,m2_in,q2
complex(dp) ::res
real(dp) :: m1,m2

m1=m1_in
m2=m2_in

if(allsmall(q2,m1,m2)) then
   res=2._dp-clog(-p2/q2)
else
   If(smalldp(m1,m2)) then
!   call saveme(q2,m1) !!! Should not need that now
      res=-clog(m1/q2)
   else
      res=(J0(m1,q2) - J0(m2,q2))/(m2-m1) 
   end if
end if


BBs = res

end function BBs


complex(dp) function BB(m1_in,m2_in,q2)
implicit none
real(dp), intent(in) :: m1_in,m2_in,q2
complex(dp) ::res
real(dp) :: m1,m2

m1=m1_in
m2=m2_in

If (newGBC) Then
if(allsmall(q2,m1,m2)) then
      res=0._dp  !!! Need this to catch the goldstone boson cases in 1L routines
else
   If(smalldp(m1,m2)) then
            
      res=-clog(m1/q2)
   else
      res=(J0(m1,q2) - J0(m2,q2))/(m2-m1) 
   end if
end if
Else 
If(smalldp(m1,m2)) then
   call saveme(q2,m1)
   res=-clog(m1/q2)
else
   res=(J0(m1,q2) - J0(m2,q2))/(m2-m1) 
end if
End if

BB = res

end function BB


complex(dp) function CC(x_in,y_in,z_in,q2)        !!!!! 1/(k+x)(k+y)(k+z)
implicit none
real(dp), intent(in) :: x_in,y_in,z_in,q2
complex(dp) :: res
real(dp) :: x,y,z

real(dp) :: invals(3)
invals = (/ x_in, y_in,z_in /)
call sortabsvals(invals)
  x=invals(1)
  y=invals(2)
  z=invals(3)


! x=x_in
! y=y_in
! z=z_in

If (NewGBC) Then
If(smalldp(y,z)) then
   If(smalldp(x,y)) then
      !      call saveme(q2,x)
      if(allsmall(q2,x)) then
         res=0._dp  !!! Need this to catch the goldstone boson cases in 1L routines
      else
         res  =0.5_dp/x
      end if
   else
      if(allsmall(q2,x)) then
         res = 1.0_dp/y
      else
         res = (-x + y + x*clog(x/y))/(x - y)**2
      end if
   end if
else
   res=(BB(x,y,q2) - BB(x,z,q2))/(z-y)
end if
Else
If(smalldp(y,z)) then
   If(smalldp(x,y)) then
      call saveme(q2,x)
      res  =1/(2._dp*x)
   else
      res = (-x + y + x*clog(x/y))/(x - y)**2
   end if
else
   res=(BB(x,y,q2) - BB(x,z,q2))/(z-y)
end if
End if

if(res.ne.res) then
   write(*,*) "Problem in CC: ",x,y,z,q2
   end if
CC = res

end function CC

complex(dp) function CCtilde(x_in,y_in,z_in,q2)        !!!!! 1/(k+x)(k+y)(k+z) but zero if y=z=0.
implicit none
real(dp), intent(in) :: x_in,y_in,z_in,q2
complex(dp) :: res

if(allsmall(q2,y_in,z_in)) then
   res=0._dp
else
   res=CC(x_in,y_in,z_in,q2)
end if

CCtilde=res

end function CCtilde



complex(dp) function BBp(x,y_in,q2)
implicit none
!!! Provides B(x,y')
real(dp), intent(in) :: x,y_in,q2
complex(dp) :: res
real(dp) ::y

y=y_in


If(smalldp(x,y)) then
   !   call saveme(q2,y)
   if(smallc(y/q2)) then
      write(ErrCan,*) "BBp called with small inputs: ",x,y_in,q2 
      res = 0._dp !!!!!! Something wrong in this case but set to zero
   else

      res = -1._dp/(2._dp*y)
   end if

else
res = (x-y-x*clog(x/y))/((x-y)**2)
end if

BBp=res
end function BBp


complex(dp) function clog(x)
 implicit none
 real(dp), intent(in) :: x

clog=log(dcmplx(x,0._dp))


end function clog





complex(dp) function J0(m1,q2) 
  implicit none
  real(dp), intent(in) ::  m1,q2
  
  if(smallc(m1/q2)) then
     J0=0._dp
  else
     J0 = m1*(clog(m1/q2)-1)
  end if

end function J0


real(dp) function II0(x,q2) 
      !!!! THIS IS I(x,x,0,q2) !!!!!!
      implicit none
      real(dp), intent(in) :: x,q2
      real(dp) :: res
      
      if(abs(x/q2).lt.epsM) then
         res = 0._dp
      else
         res = real(-x* ( (clog(x/Q2))**2 - 4._dp  * cLog(x/q2) + 5._dp ),dp)
      end if

      II0 = res

end function II0



real(dp) Function IIp(x_in,y_in,z_in,Q2)
Implicit None
!!! Calculates the D[II(....),x] 
Real(dp),Intent(in) :: x_in,y_in,z_in,Q2
real(dp) :: x,y,z
Complex(dp) :: res


x=x_in
!call saveme(q2,x)

if(smallc(x/q2)) then
   write(ErrCan,*) "IIp called with small input: ",x,y_in,z_in,q2
   IIp=0._dp
   return
end if


if(abs(z_in) .lt. abs(y_in)) then
   y = z_in
   z = y_in
else
   y = y_in
   z = z_in
end if 

if(smallc(y/q2)) then
   if(smallc(z/q2)) then
      res=real(-0.5_dp - Pi**2/6._dp + cLog(x/Q2) - 0.5_dp*cLog(x/Q2)**2,dp)
   else
      res=real(-CLi2_mod(dcmplx(1 - z/x,0._dp))  - 0.5_dp*(1._dp-cLog(x/Q2))**2,dp)
      
   end if
else
if(smalldp(x,y) .and. smalldp(x,z)) then
res = 0.5_dp*real(-(1._dp- cLog(x/Q2))**2 + phimod(x, x, x),dp)
else
res = real(0.5_dp*(-1._dp + clog(y/Q2)* clog(z/Q2) - cLog(x/Q2) *(-2._dp + cLog(y/Q2) + cLog(z/Q2))) -0.5_dp*((x-y-z)/z)*phimod(x,y,z),dp)
end if

end if

IIp=res

end function IIp






real(dp) function UU(x,y,z,u,Q2)
implicit none
real(dp), intent(in) :: x,y,z,u,Q2
real(dp) :: res



if(smalldp(x,y)) then
res = -real(IIp(x,z,u,q2),dp)
else
res = real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
if(res.ne.res) then
   write(*,*) "Problem in UU: ",x,y,z,u,II(x,z,u,q2),II(y,z,u,q2)
   end if
end if

UU=res
end function UU

complex(dp) function UUs(p2,x,y,zin,uin,Q2)
implicit none
real(dp), intent(in) :: p2,x,y,zin,uin,Q2
complex(dp) :: res
real(dp) :: z,u

if(abs(zin) .gt. abs(uin)) then !!! order the arguments so that z <= u
   z=uin
   u=zin
else
   z=zin
   u=uin
end if


if(allsmall(q2,x,y)) then
   ! U(0,0,z,u) needs momentum
   if(allsmall(q2,z,u)) then
      res=1._dp+0.25_dp*((clog(-p2/q2)-3._dp)**2)

      if(res.ne.res) then
         write(*,*) "Problem in UUs(0,0,0,0): ",x,y,z,u,q2,p2
      end if
      
   else
      
      
      if(smalldp(z,u)) then
         res=clog(z/q2)*clog(-p2/q2)-1.5_dp-3._dp*clog(z/q2)-0.25_dp*(clog(z/q2)**2)
      else
         
         if(smallc(z/q2)) then
            res=(clog(u/q2)-1._dp)*clog(-p2/q2)-Pi**2/6._dp+2.5_dp-clog(u/q2)*(2._dp+0.5_dp*clog(u/q2))
         else
            res=-BB(q2,z,u)*clog(-p2/q2)+2.5_dp+(0.5_dp/(u-z))*(-(z+u)*((clog(u/q2))**2) +4._dp*z*clog(z/q2)-4._dp*u*clog(u/q2)+2._dp*z*clog(z/q2)*clog(u/q2)-2._dp*(z+u)*Cli2_mod((1._dp,0._dp)-z/u))
         end if
      end if
   end if
   if(res.ne.res) then
         write(*,*) "Problem in special UUs: ",x,y,z,u,q2,p2
      end if

   
else
   if(smalldp(x,y)) then
      res = -real(IIp(x,z,u,q2),dp)
   else
      res = real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
      ! if(res.ne.res) then
      !    write(*,*) "Problem in UU: ",x,y,z,u,II(x,z,u,q2),II(y,z,u,q2)
      ! end if
   end if
end if
if(res.ne.res) then
         write(*,*) "Problem in UUs: ",x,y,z,u,q2,p2
      end if
!write(*,*) "uus: ",res
UUs=res
end function UUs





real(dp) function VV0(x_in,y_in,u_in,v_in,Q2)
implicit none
real(dp), intent(in) :: x_in,y_in,u_in,v_in,Q2
real(dp) :: res,x,y,z,u,v

integer :: section

call atvals(x,y,u,v,z,x_in,y_in,u_in,v_in)    !!!!! note called with x,y,u,v ...

if(abs(u) .gt. abs(v)) then
   u = v_in
   v = u_in
end if


if(smalldp(x,y)) then
   
!!!!! 1/2d^2/dx^2 I(x,u,v)
   if(smallc(x/q2)) then ! always singular if x=0
      write(Errcan,*) "VV0 called with small inputs: ",x_in,y_in,u_in,v_in,Q2
      res=0._dp
   else
      if(smallc(u/q2)) then
         if(smallc(v/q2)) then 
            !! I(x'',0,0)/2
            !         call saveme(q2,x) !!!! always singular if x=0
            res = real((1._dp - cLog(x/Q2))/(2._dp*x),dp)      
         else
            !! I(x'',0,v)
            if(smalldp(x,v)) then
!!! I(x'',0,x)/2 
               res = real(-cLog(x/Q2)/(2.*x),dp)
            else
               res = real(-(v*(-v + x)*cLog(v/Q2) + v**2*cLog(v/x) + (v - x)*(-x + (v + x)*cLog(x/Q2)))/(2.*(v - x)*x**2),dp)
            end if
         end if
      else
!!! u,v are not small 
          res = real((-(deltanm(x,u,v)*(-2 + cLog(u/Q2) + cLog(v/Q2))) + (u + v - x)*((u - v)*cLog(v/u) + x*cLog(x**2/(u*v))))/(4.*x*deltanm(x,u,v)),dp)
         res=res+real(u*phimod(x,u,v)/deltanm(x,u,v),dp)


      end if
   end if
else
   res = real((IIp(y,u,v,q2) + UU(x,y,u,v,q2))/(y-x),dp)
   if(res.ne.res) then
         write(*,*) "problem VV0: ",x,y,u,v,IIp(y,u,v,q2),UU(x,y,u,v,q2)
         end if

end if

VV0=res
end function VV0






real(dp) function MM(x_in,y_in,z_in,u_in,v_in,q2)
implicit none
real(dp), intent(in) :: x_in,y_in,z_in,u_in,v_in,Q2
real(dp) :: res,x,y,z,u,v



if(smalldp(x_in,z_in)) then
   ! swap x,z and y,u
  call atvals(y,x,u,z,v,x_in,y_in,z_in,u_in,v_in)
else
   call atvals(x,y,z,u,v,x_in,y_in,z_in,u_in,v_in)
end if

if(smalldp(u,y)) then
   !   call saveme(q2,x,y) -> we shouldn't get here, the new routines should take over
   if(smallc(y/q2)) then
      write(Errcan,*) "MM somehow called with small inputs: ",x_in,y_in,z_in,u_in,v_in,q2
      MM=0._dp
      return
   end if

   if(smalldp(x,z)) then
      ! result is I(x',y',v)
      if(smallc(x/q2)) then
         write(Errcan,*) "MM somehow called with small inputs: ",x_in,y_in,z_in,u_in,v_in,q2
         MM=0._dp
         return
      end if

      if(smalldp(x,y)) then      
         !Result is I(x',x',v)
         if(smallc(v/q2)) then
            !Result is I(x',x',0)
            res = real(1/x,dp)
         else
            res=real((v - 2._dp*x)*phimod(x,x,v)/(v*(v-4._dp*x)),dp)
            res=res+real(2._dp*clog(x/v)/(v-4._dp*x),dp)
         end if
      else
         If(smallc(v/q2)) then
            ! result is I(x',y',0) with x =!= y
            res = real(clog(x/y)/(x-y),dp)
         else
            res = real((deltanm(x,y,v)*(v*(cLog(v/x)) + y*phimod(x,y,v)) - (v - x + y)*(v*((v - x)*cLog(v/x) - y*(cLog(y/v) + cLog(y/x))) + y*(-v - x + y)*phimod(x,y,v)))/(2.*v*y*deltanm(x,y,v)),dp)
         end if
      end if
      
      
   else
      !! Exploit symmetry of the function here ...
      res=real((IIp(y,x,v,q2)-IIp(y,z,v,q2))/(x-z),dp)
   end if
else

      res=real((UU(x,z,y,v,q2) - UU(x,z,u,v,q2))/(u-y),dp)


end if

MM=res
end function MM


 complex(dp) Function phimod(x,y,z)
 !--------------------------------------------------------------
 ! from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 ! version by Pietro Slavich
 ! 12.03.02: portation to f90
 !--------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y, z
  
  complex(dp) :: xd,yd,zd
  xd=x
  yd=y
  zd=z
  
  If(abs(x).Le.abs(z).And.abs(y).Le.abs(z)) Then
     phimod = myphimod(x,y,z)
  Elseif(abs(z).Le.abs(x).And.abs(y).Le.abs(x)) Then
     phimod = z/x*myphimod(z,y,x)
  Elseif(abs(z).Le.abs(y).And.abs(x).Le.abs(y)) Then
     phimod = z/y*myphimod(z,x,y)
  Endif


  if(phimod.ne.phimod) then
     write(*,*) "wtf is up with phimod? ",x,y,z
     end if

 Contains

  Real(dp) Function myphimod(x,y,z)
   Implicit None
   Real(dp), Intent(in) :: x,y,z
   complex(dp) :: u,v
   Complex(dp) :: clam,cxp,cxm,ccphi

   !     auxiliary variables
   u = x/z
   v = y/z
   If(abs(u)<=epsm) Then
     
     If(abs(v)/=1._dp) Then
       myphimod = real((Log(u)*Log(v)+2._dp*CLi2_mod(1._dp-v))/(1._dp-v),dp)
     Else
       myphimod = real(2._dp-Log(u),dp)
     Endif

     Elseif(abs(v)<=epsM) Then
     If(abs(u)/=1._dp) Then
        myphimod = real((Log(v)*Log(u)+2._dp*CLi2_mod(1._dp-u))/(1._dp-u),dp)
     Else
       myphimod = real(2._dp-Log(v),dp)
     Endif
   Else
     
       clam = Sqrt((1._dp-u-v)**2 - 4._dp*u*v)

    cxp = (1._dp+(u-v)-clam)/2._dp
    cxm = (1._dp-(u-v)-clam)/2._dp
       
     ccphi = (2._dp*Log(cxp)*Log(cxm) - Log(u)*Log(v)  &
& - 2._dp*(CLI2_mod(cxp) + CLI2_mod(cxm)) + Pi2o3)/clam
    myphimod = Real(ccphi,dp)
       
   Endif
   Return

  End Function myphimod
  
 End Function phimod  







real(dp) function deltaNM(x,y,z)
real(dp), intent(in)::  x,y,z

      deltaNM = x**2+y**2+z**2 -2._dp*(x*y+x*z+y*z)
end function deltanm


Complex(dp) Function RFunction(x,y,z)
Implicit None
Complex(dp), Intent(in) :: x,y,z
Complex(dp) :: temp
temp = x**2+y**2+z**2-2._dp*x*y-2._dp*x*z-2._dp*y*z
RFunction = sqrt(temp)

End Function RFunction




Complex(dp) Function II(x_in,y_in,z_in,Q2)
Implicit None
! x_in,y_in,z_in are supposed to be squared masses
! symmetric under interchange of any two of these
Real(dp),Intent(in) :: x_in,y_in,z_in,Q2
Complex(dp) :: res,x,y,z
x = x_in
y = y_in
z = z_in
If (Q2.le.0._dp) Then 
  write(*,*) 'invalid renormalization scale!'
  return
end if
Call SortArguments3C(x,y,z) ! sort x<=y<=z


if(smallc(real(x/q2))) then
   if(smallc(real(y/q2))) then
      if(smallc(real(z/q2))) then
         II = ZeroC ! I(0,0,0)=0
      else  !I(0,0,z)
         II = -0.5_dp*z*(log(z/Q2))**2 + 2*z*log(z/Q2)-2.5_dp*z-pi**2*z/6._dp
      end if
   else !I(0,y,z)
      if(smalldp(y,z)) then ! I(0,z,z)
         II = -z*(log(z/Q2))**2+4._dp*z*log(z/Q2)-5._dp*z
      else !I(0,y,z)
         II = (z-y)*(CLi2_mod(y/z)-log(z/y)*log((z-y)/Q2) &
   & +0.5_dp*(log(z/Q2))**2-pi**2/6._dp)  &
   & -2.5_dp*(z+y)+2*y*log(y/Q2)+2*z*log(z/Q2) & 
   & -z*log(z/Q2)*log(y/Q2)
      end if
   end if

else




res = Cmplx(0._dp,0._dp)
res=res + 0.5_dp*(x-y-z)*log(y/Q2)*log(z/Q2)
res=res + 0.5_dp*(y-x-z)*log(x/Q2)*log(z/Q2)
res=res + 0.5_dp*(z-x-y)*log(x/Q2)*log(y/Q2)
res=res + 2._dp*x*log(x/Q2)
res=res + 2._dp*y*log(y/Q2)
res=res + 2._dp*z*log(z/Q2)
res=res - 2.5_dp*(x+y+z)
res=res - 0.5*XiFunction(x,y,z)
II = res
if(res.ne. res) then
   write(*,*) "Problem in general II: ",x,y,z,xifunction(x,y,z)
end if

End If

End Function II








Complex(dp) Function XiFunction(x0,y0,z0)
Implicit None
Complex(dp), Intent(in) :: x0,y0,z0
Complex(dp) :: x,y,z
Complex(dp) :: res,R
x=x0
y=y0
z=z0
Call SortArguments3C(x,y,z)
R=Rfunction(x,y,z)
res = Cmplx(0.,0.,dp)
res = res + 2._dp*log((z+x-y-R)/(2._dp*z))*log((z+y-x-R)/(2._dp*z))
res = res - log(x/z)*log(y/z)
res = res - 2._dp*CLi2_mod((z+x-y-R)/(2._dp*z))
res = res - 2._dp*CLi2_mod((z+y-x-R)/(2._dp*z))
res = res + (pi**2)/3._dp
res = R*res
XiFunction = res
End Function XiFunction

!! sorts three real arguments a,b,c from small to large
!! and returns them
Subroutine SortArguments3(a,b,c)
Implicit None
Real(dp), Intent(inout) :: a,b,c
Real(dp) :: x,y,z
x=a ! assume a is smallest
y=b
z=c
If (b.le.a) Then
 x=b
 y=a
 z=c
End If
If (c.le.x) Then
 x=c
 y=a
 z=b
End If
! compare last two elements
If (z.le.y) Then
 a=x ! output
 b=z
 c=y
Else
 a=x !output
 b=y
 c=z
End If
End Subroutine SortArguments3




Subroutine SortArguments3C(a,b,c)
Implicit None
Complex(dp), Intent(inout) :: a,b,c
Complex(dp) :: x,y,z
x=a ! assume a is smallest
y=b
z=c
If (Abs(b).le.Abs(a)) Then
 x=b
 y=a
 z=c
End If
If (Abs(c).le.Abs(x)) Then
 x=c
 y=a
 z=b
End If
! compare last two elements
If (Abs(z).le.Abs(y)) Then
 a=x ! output
 b=z
 c=y
Else
 a=x !output
 b=y
 c=z
End If
End Subroutine SortArguments3C


 Complex(DP) Function CLI2_mod(Z)
 !--------------------------------------------------------------
 !     routine for the complex dilogarithm function
 !--------------------------------------------------------------
 Implicit None      
  !-------
  ! input
  !-------
  Complex(DP), Intent(in) :: Z

  !-----------------
  ! local variables
  !-----------------
  Complex(DP) :: sumC, cy, cz, cz2, CLI2
  Real(dp) :: rz, iz, az
  Integer :: i1, jsgn, ipi12
  !---------------------------------------------
  ! Bernulli numbers / (2 n + 1)!
  !---------------------------------------------
  Real(dp), Parameter :: bf(20) = (/ - 1._dp/4._dp, + 1._dp/36._dp        &
       & , - 1._dp/36.e2_dp, + 1._dp/21168.e1_dp,  - 1._dp/108864.e2_dp   &
       & , + 1._dp/52690176.e1_dp,  - 691._dp/16999766784.e3_dp           &
       & , + 1._dp/1120863744.e3_dp,  - 3617._dp/18140058832896.e4_dp     &
       & , + 43867._dp/97072790126247936.e3_dp                            &
       & , - 174611._dp/168600109166641152.e5_dp                          &
       & , + 77683._dp/32432530090601152512.e4_dp                         &
       & , - 236364091._dp/4234560341829359173632.e7_dp                   &
       & , + 657931._dp/5025632054039239458816.e6_dp                      &
       & , - 3392780147._dp/109890470493622010006470656.e7_dp             &
       & , +172.3168255201_dp/2355349904102724211909.3102313472e6_dp      &
       & , -770.9321041217_dp/4428491985594062112714.2791446528e8_dp      &
       & ,  4.1576356446138997196178996207752266734882541595116e-29_dp    &
       & , -9.9621484882846221031940067024558388498548600173945e-31_dp    &
       & ,  2.3940344248961653005211679878937495629342791569329e-32_dp    &
       & /)

  rz = Real(z,dp)
  iz = Aimag(z)
  az = Sqrt(rz**2 + iz**2)
  !------------------
  ! special cases
  !------------------
  If (iz.Eq.0._dp) Then
    If (rz.gt.1._dp) Then
    CLi2_mod = Cmplx( Li2(rz), -pi*log(rz),dp)
    Else
    CLi2_mod = Cmplx( Li2(rz), 0.,dp)
    End If
   Return
  Else If (az.Lt.Epsilon(1._dp)) Then
   CLi2_mod = z
   Return
  End If
  !-----------------------------------------------------
  ! transformation to |z|<1, Re(z)<=0.5
  !-----------------------------------------------------
  If (rz.Le.0.5_dp) Then
   If (az.Gt.1._dp) Then
    cy = -0.5_dp * Log(-z)**2
    cz = - Log(1._dp - 1._dp / z)
    jsgn = -1
    ipi12 = -2
   Else !  (az.le.1._dp)
    cy = 0
    cz = - Log(1._dp - z)
    jsgn = 1
    ipi12 = 0
   End If
  Else ! rz.gt.0.5_dp
   If (az.Le.Sqrt(2*rz)) Then
    cz = - Log(z)
    cy = cz * Log(1._dp - z)
    jsgn = -1
    ipi12 = 2
   Else ! (az.gt.Sqrt(2*rz))
    cy = -0.5_dp * Log(-z)**2
    cz = - Log(1._dp - 1._dp / z)
    jsgn = -1
    ipi12 = -2
   End If
  End If
  !--------------------------------------------
  ! the dilogarithm
  !--------------------------------------------
  cz2 = cz**2
  sumC = cz2 * bf(20)
  Do i1=3,19
   sumC = cz2 * ( sumC + bf(22-i1) )
  End Do
  ! watch the powers of z
  sumC = cz + cz2 * (bf(1) + cz * (bf(2) + sumC))
  If (jsgn.Eq.1) Then
   CLi2 = sumC + cy
  Else
   CLi2 = - sumC + cy
  End If
  If (ipi12.Ne.0) Cli2 = Cli2 + ipi12 * pi2o12
  
  CLI2_mod = CLI2
 End Function CLI2_mod




Logical Function SmalldpC(x,y) 
Implicit None
Complex(dp), Intent(in) :: x,y
If (Max(Abs(x),Abs(y)).lt.epsD) Then
 SmalldpC = .True.
Else If ((Abs(x-y)/Max(Abs(x),Abs(y))).le.epsD) Then
 SmalldpC = .True.
Else
 SmalldpC = .False.
End if
End Function SmalldpC


Logical Function SmalldpR(x,y) 
Implicit None
Real(dp), Intent(in) :: x,y
If (Max(Abs(x),Abs(y)).lt.epsD) Then
 SmalldpR = .True.
Else If ((Abs(x-y)/Max(Abs(x),Abs(y))).le.epsD) Then
 SmalldpR = .True.
Else
 SmalldpR = .False.
End if
End Function SmalldpR






Logical Function smallc(x,y,z,u,v)
Implicit None
real(dp), Intent(in) :: x
real(dp), Intent(in),optional :: y,z,u,v

Logical :: res
res = .False.
smallc=.False.

If(abs(x) .lt. epsM) then
   smallc = .True.
else 
   if(present(y)) then
      if(abs(y) .lt. epsM) then
         smallc = .True.
      end if
   end if
   if(present(z)) then
      if(abs(z) .lt. epsM) then
         smallc = .True.
      end if
   end if
   if(present(u)) then
      if(abs(u) .lt. epsM) then
         smallc = .True.
      end if
   end if
   if(present(v)) then
      if(abs(v) .lt. epsM) then
         smallc = .True.
      end if
   end if
end if





end function smallc


subroutine saveme(q2,x,y,z,u,v)
implicit none

real(dp), intent(in) ::q2
real(dp), Intent(inout) :: x
real(dp), intent(inout), optional :: y,z,u,v

real(dp) :: minm


minm=epsM*Q2

If(abs(x) .lt. minm) then
   x = minm
end if

if(present(y)) then
   if(abs(y) .lt. minM) then
!       write(*,*) "Saveme invoked!"
       y = minM
   end if
end if
if(present(z)) then
   if(abs(z) .lt. minM) then
!      write(*,*) "Saveme invoked!"
       z = minM
   end if
end if
if(present(u)) then
   if(abs(u) .lt. minM) then
!      write(*,*) "Saveme invoked!"
       u = minM
   end if
end if
if(present(v)) then
   if(abs(v) .lt. minM) then
!      write(*,*) "Saveme invoked!"
       v = minM
   end if
end if



end subroutine saveme

subroutine atvals(x,y,z,u,v,x_in,y_in,z_in,u_in,v_in)
implicit none

real(dp), intent(in), optional :: x_in,y_in,z_in,u_in,v_in
real(dp), intent(out) :: x,y,z,u,v

x=x_in

If(present(y_in)) then
   y=y_in
else
   y=0._dp
end if
If(present(z_in)) then
   z=z_in
else
   z=0._dp
end if
If(present(u_in)) then
   u=u_in
else
   u=0._dp
end if
If(present(v_in)) then
   v=v_in
else
   v=0._dp
end if

end subroutine atvals





Logical Function allsmall(q2,x,y,z)
Implicit None
real(dp), Intent(in) :: q2,x
real(dp), Intent(in),optional :: y,z
real(dp) :: gbm2

allsmall=.True.

gbm2=q2*epsM
If(abs(x) .gt. gbm2) then
   allsmall = .False.
end if
if(present(y)) then
      if(abs(y) .gt. gbm2) then
         allsmall= .False.
      end if
end if

if(present(z)) then
      if(abs(z) .gt. gbm2) then
         allsmall = .False.
      end if
end if



end function allsmall



SUBROUTINE sortabsvals(a)
    real(dp), INTENT(IN OUT) :: a(:)
    INTEGER :: i, minIndex
    real(dp) :: temp
    
    DO i = 1, SIZE(a)-1
       minIndex = MINLOC(abs(a(i:)), 1) + i - 1
       IF (abs(a(i)) > abs(a(minIndex))) THEN
          temp = a(i)
          a(i) = a(minIndex)
          a(minIndex) = temp
       END IF
    END DO
END SUBROUTINE sortabsvals






End module Pole2LFunctions
