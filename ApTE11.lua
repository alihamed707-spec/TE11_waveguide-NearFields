application = cf.Application.GetInstance()
project = application:NewProject()

-- Create a magnetic dipoles on   
-- the aperture of a TE11 circular waveguide 
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
 magneticDipolex =
 project.Contents.SolutionConfigurations.GlobalSources:AddMagneticDipole(cf.Point(x,y,0),90,0)
 magneticDipoley =
 project.Contents.SolutionConfigurations.GlobalSources:AddMagneticDipole(cf.Point(x,y,0),90,90)
  sum0=0
 sum2=0
 for k=0, 40 do
 Ak=1
 for n = 1, k do
 Ak = Ak*n
 end
 Bk = Ak*(k+1)*(k+2)
 sum2=sum2+(-1)^k*(1.841*rho/(2*a))^(2*k+2)/(Ak*Bk)
 sum0=sum0+(-1)^k*(1.841*rho/(2*a))^(2*k)/(Ak*Ak)
 end
 if rho==0 then
 Amy=0
 Amx=(1.841/a)*sum0
 else
 Amy=-2*(1.841/a)*sum2*x*y/(rho^2)
 Amx=(1.841/a)*(sum0-(1-2*y^2/(rho^2))*sum2)
 end
 magneticDipolex.Magnitude=Amx
 magneticDipoley.Magnitude=Amy
 
 end
end
end
