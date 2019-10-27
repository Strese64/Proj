# v1.0-5 Beta at 2019-10-27
import Proj as pj 

# create a Proj-instance
t = pj.Proj(fromEPSG=4258, toEPSG=3035)

# set the longitude l and latitude b
l = 9.54152
b = 52.41955

# transform the coordinate to the Lambert LAEA projection
(x1, y1, z1) = t.transform(l,b)

# Output
print("Project ...")
print("EPSG:4258   l = {:.5f}°  b = {:.5f}°  to".format(l,b))
print("EPSG:3035   x = {:.1f}   y = {:.1f}".format(x1,y1))
