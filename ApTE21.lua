application = cf.Application.GetInstance()
project = application:NewProject()

-- Create a magnetic dipoles on   
-- the aperture of a TE21 circular waveguide 
-- a is the raduis of the aperture
-- d is the spacing between any two adjacent dipoles.
-- Amx is current strength of the x-directed dipole
-- Amy is current strength of the y-directed dipole
-- the point (x,y) is the location of the dipole
-- the wavelength is assumed = 1.0

a=5
d=5/40
for q = -40, 40 do
for p = -40, 40 do
x=q*d
y=p*d

  rho=math.sqrt(x^2+y^2)
 if rho<=a then
 xo=x
 yo=y
 magneticDipolex =
 project.Contents.SolutionConfigurations.GlobalSources:AddMagneticDipole(cf.Point(xo,yo,0),90,0)
 magneticDipoley =
 project.Contents.SolutionConfigurations.GlobalSources:AddMagneticDipole(cf.Point(xo,yo,0),90,90)
 sum1=0
 sum3=0
 for k=0, 40 do
 Ak=1
 for n = 1, k do
 Ak = Ak*n
 end
 Bk = Ak*(k+1)*(k+2)*(k+3)
 sum3=sum3+(-1)^k*(3.0542*rho/(2*a))^(2*k+3)/(Ak*Bk)
 sum1=sum1+(-1)^k*(3.0542*rho/(2*a))^(2*k+1)/(Ak*Ak*(k+1))
 end
 if rho==0 then
 Amy=0
 Amx=0
 else
 Amy=-(3.0542/a)*(yo/rho)*(sum1+(3-4*(yo/rho)^2)*sum3)
 Amx=(3.0542/a)*(xo/rho)*(sum1-(4*(xo/rho)^2-3)*sum3)
 end
 magneticDipolex.Magnitude=Amx
 magneticDipoley.Magnitude=Amy
 
 end
end
end
