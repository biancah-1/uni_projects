TITLE 'Bimorph Actuator'     
COORDINATES cartesian3 
VARIABLES      
u !x displacement
v !y displacement 
w !z displacement

SELECT      

DEFINITIONS    
mag = 0.1*globalmax(magnitude(x,y,z))/globalmax(magnitude(U,V,W))

!dimensions
Lx=2.5e-3
Ly=25e-3
Lz=1.5e-3

ratio=0.5
Lz_top=ratio*Lz !thickness of top layer 
Lz_bot=Lz-Lz_top !thicknes of bottom layer 

Vmax=12
Efieldx=0
Efieldy=0
Efieldz= if(z<0) then -Vmax/(Lz_bot) else Vmax/(Lz_top)


!stiffness matrix
C11 = 1.20346249598166*10^11, C21 = 7.51791312603156*10^10, C31 = 7.50900664786295*10^10, C41 = 0., C51 = 0., C61 = 0., C12 = 7.51791312603156*10^10, C22 = 1.20346249598166*10^11, C32 = 7.50900664786295*10^10, C42 = 0., C52 = 0., C62 = 0., C13 = 7.50900664786295*10^10, C23 = 7.50900664786295*10^10, C33 = 1.10867051061245*10^11, C43 = 0., C53 = 0., C63 = 0., C14 = 0., C24 = 0., C34 = 0., C44 = 2.10526315789474*10^10, C54 = 0., C64 = 0., C15 = 0., C25 = 0., C35 = 0., C45 = 0., C55 = 2.10526315789474*10^10, C65 = 0., C16 = -0., C26 = -0., C36 = -0., C46 = -0., C56 = -0., C66 = 2.25733634311512*10^10

!piezoelectric coupling coefficients
d11 = 0				d12 = 0		d13 = 0			d14 = 0		d15 = 584e-12	d16 = 0
d21 = 0				d22 = 0		d23 = 0			d24 = d15	d25 = 0			d26 = 0
d31 = -171e-12	d32 = d31	d33 = 374e-12	d34 = 0		d35 = 0	d36=0

!strain definition 
ex = dx(u)
ey = dy(v)
ez = dz(w)
gyz = dy(w) + dz(v)
gxz = dx(w) + dz(u)
gxy = dx(v) + dy(u)

!mechanical strain
exm  = ex - (d11*Efieldx+d21*Efieldy+d31*Efieldz)
eym  = ey   - (d12*Efieldx+d22*Efieldy+d32*Efieldz)
ezm  = ez  - (d13*Efieldx+d23*Efieldy+d33*Efieldz)
gyzm = gyz  - (d14*Efieldx+d24*Efieldy+d34*Efieldz)
gxzm = gxz  - (d15*Efieldx+d25*Efieldy+d35*Efieldz)
gxym = gxy  - (d16*Efieldx+d26*Efieldy+d36*Efieldz)

!hooke's law
sx  = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy  = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz  = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
syz = C41*exm + C42*eym + C43*ezm + C44*gyzm + C45*gxzm + C46*gxym
sxz = C51*exm + C52*eym + C53*ezm + C54*gyzm + C55*gxzm + C56*gxym
sxy = C61*exm + C62*eym + C63*ezm + C64*gyzm + C65*gxzm + C66*gxym


INITIAL VALUES
EQUATIONS      
u: dx(sx) + dy(sxy) + dz(sxz) = 0
v: dx(sxy) + dy(sy) + dz(syz) = 0
w: dx(sxz) + dy(syz) + dz(sz) = 0
EXTRUSION
surface 'bottom' z = -Lz_bot
surface 'mid' z = 0
surface 'top' z = Lz_top

BOUNDARIES    
  REGION 1       
    START(0,0) !y=0
		value(u) = 0 value(v) = 0 value(w) = 0 !fixes y=0 face, cantilever condition
    LINE TO (Lx,0) !x = Lx
		load(u)= 0 load(v) = 0 load(w) = 0 !to remove previous fixed conditions 
	LINE TO (Lx,Ly) !y = Ly
	LINE TO (0,Ly) !x = 0
	LINE TO CLOSE

PLOTS          
grid(x+mag*u, y+mag*v, z+mag*w)
elevation(w) from (Lx/2,0) to (Lx/2,Ly) 
contour(sz) on x=Lx/2 painted 
SUMMARY
report val(w,Lx/2,Ly,Lz/2)
END
