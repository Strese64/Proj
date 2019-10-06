import Proj as pj 

t = pj.Proj(von=4258, nach=3035)
l = 9.54152
b = 52.41955
(x1,y1) = t.transform(l,b)

print("Projeziere ...")
print("EPSG:4258   l = {:.5f}  b = {:.5f}".format(l,b))
print("EPSG:3035   x = {:.1f}  y = {:.1f}".format(x1,y1))
