import matplotlib.pyplot as plt
import numpy as np
import subprocess
#flex code
BimorphFlexCode='''
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
ratio={RATIO}
Lz_top=ratio*Lz !thickness of top layer
Lz_bot=Lz-Lz_top !thicknes of bottom layer
Vmax=12
Efieldx=0
Efieldy=0
Efieldz= if(z<0) then -Vmax/(Lz_bot) else Vmax/(Lz_top)
!stiffness matrix
C11 = 1.20346249598166*10^11, C21 = 7.51791312603156*10^10, C31 =
7.50900664786295*10^10, C41 = 0., C51 = 0., C61 = 0., C12 =
7.51791312603156*10^10, C22 = 1.20346249598166*10^11, C32 =
7.50900664786295*10^10, C42 = 0., C52 = 0., C62 = 0., C13 =
7.50900664786295*10^10, C23 = 7.50900664786295*10^10, C33 =
1.10867051061245*10^11, C43 = 0., C53 = 0., C63 = 0., C14 = 0., C24 = 0., C34
= 0., C44 = 2.10526315789474*10^10, C54 = 0., C64 = 0., C15 = 0., C25 = 0.,
C35 = 0., C45 = 0., C55 = 2.10526315789474*10^10, C65 = 0., C16 = -0., C26 =
-0., C36 = -0., C46 = -0., C56 = -0., C66 = 2.25733634311512*10^10
!piezoelectric coupling coefficients
d11 = 0 d12 = 0 d13 = 0 d14 = 0 d15 = 584e12 d16 = 0
d21 = 0 d22 = 0 d23 = 0 d24 = d15 d25 = 0
d26 = 0
d31 = -171e-12 d32 = d31 d33 = 374e-12 d34 = 0 d35 = 0 d36=0
!strain definition
ex = dx(u)
ey = dy(v)
ez = dz(w)
gyz = dy(w) + dz(v)
gxz = dx(w) + dz(u)
gxy = dx(v) + dy(u)
!mechanical strain
exm = ex - (d11*Efieldx+d21*Efieldy+d31*Efieldz)
eym = ey - (d12*Efieldx+d22*Efieldy+d32*Efieldz)
ezm = ez - (d13*Efieldx+d23*Efieldy+d33*Efieldz)
gyzm = gyz - (d14*Efieldx+d24*Efieldy+d34*Efieldz)
gxzm = gxz - (d15*Efieldx+d25*Efieldy+d35*Efieldz)
gxym = gxy - (d16*Efieldx+d26*Efieldy+d36*Efieldz)
!hooke's law
sx = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
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
 value(u) = 0
 value(v) = 0
 value(w) = 0
 LINE TO (Lx,0) !x = Lx
 load(u)= 0
 load(v) = 0
 load(w) = 0

 LINE TO (Lx,Ly) !y = Ly
 LINE TO (0,Ly) !x = 0
 LINE TO CLOSE
PLOTS
grid(x+mag*u, y+mag*v, z+mag*w)
elevation(w) from (Lx/2,0) to (Lx/2,Ly) PrintOnly Export Format '#x#b#1' file
= 'BimorphData.txt'
SUMMARY
report val(w,Lx/2,Ly,Lz/2)
END
'''
UnimorphFlexCode = '''
TITLE 'Unimorph Actuator'
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
ratio={RATIO}
Lz_top=ratio*Lz !thickness of top layer
Lz_bot=Lz-Lz_top !thicknes of bottom layer
Vmax=12
Efieldx=0
Efieldy=0
Efieldz= if(z<0) then 0 else Vmax/(Lz_top)
!aluminum parameters
E= {EMOD}
nu= {NU}
!stiffness matrix
C11 =if(z<0) then E*(1-nu)/(1+nu)/(1-2*nu) else 1.20346249598166*10^11
C12 = if(z<0) then E*nu/(1+nu)/(1-2*nu) else 7.51791312603156*10^10
C13 = if(z<0) then E*nu/(1+nu)/(1-2*nu) else 7.50900664786295*10^10
C14 = 0 C15 = 0 C16 = 0
C21=C12 C22=C11 C23 = C13
C24 = 0 C25 = 0 C26 = 0
C31 = C13 C32 = C13
C33 = if(z<0) then E*(1-nu)/(1+nu)/(1-2*nu) else 1.10867051061245*10^11
C34 = 0 C35 = 0 C36 = 0
C41 = 0 C42 = 0 C43 = 0
C44 = if (z<0) then E/(2*(1+nu)) else 2.10526315789474*10^10
C45 = 0 C46 = 0
C51 = 0 C52 = 0 C53 = 0
C54 = 0 C55 = C44 C56 = 0
C61 = 0 C62 = 0 C63 = 0
C64 = 0 C65 = 0
C66 = if (z<0) then E/(2*(1+nu)) else 2.25733634311512*10^10
!piezoelectric coupling coefficients
d11 = 0 d12 = 0 d13 = 0 d14 = 0 d15 =
if(z<0) then 0 else 584e-12
d16 = 0
d21 = 0 d22 = 0 d23 = 0 d24 = d15 d25 = 0
d26 = 0
d31 =if(z<0) then 0 else -171e-12 d32 = d31 d33 =if(z<0) then 0 else 374e12 d34 = 0 d35 = 0 d36=0
!strain definition
ex = dx(u)
ey = dy(v)
ez = dz(w)
gyz = dy(w) + dz(v)
gxz = dx(w) + dz(u)
gxy = dx(v) + dy(u)
!mechanical strain
exm = ex - (d11*Efieldx+d21*Efieldy+d31*Efieldz)
eym = ey - (d12*Efieldx+d22*Efieldy+d32*Efieldz)
ezm = ez - (d13*Efieldx+d23*Efieldy+d33*Efieldz)
gyzm = gyz - (d14*Efieldx+d24*Efieldy+d34*Efieldz)
gxzm = gxz - (d15*Efieldx+d25*Efieldy+d35*Efieldz)
gxym = gxy - (d16*Efieldx+d26*Efieldy+d36*Efieldz)
!hooke's law
sx = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
syz = C41*exm + C42*eym + C43*ezm + C44*gyzm + C45*gxzm + C46*gxym
sxz = C51*exm + C52*eym + C53*ezm + C54*gyzm + C55*gxzm + C56*gxym
sxy = C61*exm + C62*eym + C63*ezm + C64*gyzm + C65*gxzm + C66*gxym
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
 value(u) = 0 value(v) = 0 value(w) = 0 !fixes y=0 face, cantilever
condition
 LINE TO (Lx,0) !x = Lx
 load(u)= 0 load(v) = 0 load(w) = 0 !to remove previous fixed condition
 LINE TO (Lx,Ly) !y = Ly
 LINE TO (0,Ly) !x = 0
 LINE TO CLOSE
PLOTS
grid(x+mag*u, y+mag*v, z+mag*w)
elevation(w) from (Lx/2,0) to (Lx/2,Ly) PrintOnly Export Format '#x#b#1' file
= '{EMOD}UnimorphData.txt'
SUMMARY
report val(w,Lx/2,Ly,Lz/2)
END
'''
ratios=np.arange(0.1,1,0.1)
#bimorph
BimorphFlexFileName="bimorph.pde"
def bimorph_best_disp(ratios_list,FlexFileName,FlexCode):
   print("Bimorph \n")
    BestRatio=-1
    BestTip=0
    wtips=[]
    for ratio in ratios:
        with open(FlexFileName, 'w') as f:
            print(FlexCode.format(RATIO=ratio),file=f)
        subprocess.run(['FlexPDE6s', '-S', FlexFileName])

        with open(FlexFileName,'r') as f:
            flexoutput=np.loadtxt('BimorphData.txt', skiprows=8)
        wtip=flexoutput[-1,-1]
        print("Ratio:", round(ratio,1),"\t Tip z-displacement:",wtip, "m")
        
        if wtip>BestTip:
            BestRatio=ratio
            BestTip=wtip
        wtips.append(wtip)

    plt.plot(ratios,wtips)
    plt.ylabel("Tip Displacement (m)")
    plt.xlabel("Ratio")
    plt.title("Bimorph Displacements at Varying Ratios")
    plt.show()
    print("For the bimorph, the maximum displacement was",BestTip,"which
    occured at a ratio of",round(BestRatio,1))
    return(BestRatio,BestTip)

#unimorph
UnimorphFlexFileName="Unimorph.pde"
#elastic material parameters to be tested
Emod=[68.9e9,78e9,130e9]
nu=[0.33,0.42,0.28]
materials=["6061 Aluminum","Gold","Silicon"]

print("\nUnimorph")
def unimorph_best_material(emods, nu, materialnames, ratios, FlexFileName,FlexCode):
    counter = 0
    bestlist = []
    for emod_val in emods:
        wtips = []
        BestRatio = -1
        BestTip = 0
        nu_val = nu[counter]
        33
        material = materialnames[counter]
        print("\n", material, "\t Elastic Modulus:", emod_val, "Pa \t Poisson's Ratio:", nu_val)
        
        for ratio in ratios:
            with open(FlexFileName, 'w') as f:
                print(FlexCode.format(RATIO=ratio, EMOD=emod_val, NU=nu_val), file=f)
            subprocess.run(['FlexPDE6s', '-S', FlexFileName])

            with open(FlexFileName, 'r') as f:
                flexoutput = np.loadtxt('{EMOD}UnimorphData.txt'.format(EMOD=emod_val), skiprows=8)
            wtip = flexoutput[-1, -1]
            print("Ratio:", round(ratio,1), "\t Tip z-displacement:", wtip,"m")

            if wtip > BestTip:
                BestRatio = ratio
                BestTip = wtip
            wtips.append(wtip)
            print("For", material, "a ratio of", round(BestRatio,1), "gave the max displacement of", BestTip, "m")
        
            #plot of the tip displacement as a function of the ratio
            plt.plot(ratios, wtips, label=material)
            counter += 1
            #keep the best tip displacement and ratio used for each material
            materialbest = [material, BestTip, BestRatio]
            bestlist.append(materialbest)

    #show and label plot of the tip displacements vs ratio for all materials
    plt.legend()
    plt.ylabel("Tip Displacement (m)")
    plt.xlabel("Ratio")
    plt.title("Unimorph Displacements at Varying Ratios")
    plt.show()

    #split bestlist into components to identify material that maximizes displacement
    material_name = [i[0] for i in bestlist]
    best_tip_disps = [i[1] for i in bestlist]
    best_ratios = [i[2] for i in bestlist]
    BestMaterial = None
    BestOverallTip = -1
    BestOverallRatio = 0
    counter2 = 0

    #comparing best tip displacement for each material to find the maximum
    displacement for all materials and ratios
    for tips in best_tip_disps:
        if tips > BestOverallTip:
            BestOverallTip = tips
            BestMaterial = material_name[counter2]
            BestOverallRatio = best_ratios[counter2]
        counter2 += 1
        print("\nFor the unimorph actuator, the best tip displacement was", BestOverallTip, "m. The material was", BestMaterial, "and the ratio used was", BestOverallRatio)
    return (BestMaterial, BestOverallRatio, BestOverallTip)


#compare unimorph and bimorph
best_bimorph_list=bimorph_best_disp(ratios,BimorphFlexFileName,BimorphFlexCode)
best_unimorph_list=unimorph_best_material(Emod,nu,materials,ratios,UnimorphFlexFileName,UnimorphFlexCode)

best_bimorph_ratio=best_bimorph_list[0]
best_bimorph_tip=best_bimorph_list[1]

best_unimorph_material=best_unimorph_list[0]
best_unimorph_ratio=best_unimorph_list[1]
best_unimorph_tip=best_unimorph_list[2]

if best_unimorph_tip>best_bimorph_tip:
    print("\nThe largest displacement was seen in the unimorph with",best_unimorph_material,"with a ratio of",best_unimorph_ratio,". This resulted in a tip displacement of", best_unimorph_tip,"m")
else:
    print("\nThe largest displacement was seen in the bimorph with a ratio of", best_bimorph_ratio, ". This resulted in a tip displacement of",best_bimorph_tip, "m")