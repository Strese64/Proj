#   -------------------------------------------------------------------
#   Proj    coordinate transformation class
#           Purpose Unittest for all 4 projection functions
#           Version 1.0-5 Beta, from 2019-10-27
#           Author  Dipl.-Ing. Michael Dreesmann
#           License MIT
#   -------------------------------------------------------------------
import unittest
from Proj import Proj
		
class ProjText(unittest.TestCase):
    def testCalculation(self):
        # check the "LAEA"-way: 4258 -> 3035 -> 4258
        l = 12.123450
        b = 51.234560
        t11 = Proj(fromEPSG=4258, toEPSG=3035)
        (x1, y1, z1) = t11.transform(l, b)
        t12 = Proj(fromEPSG=3035, toEPSG=4258)
        (x2, y2, z2) = t12.transform(x1, y1)
        self.assertEqual( int(abs(x2-l)*100000), 0 )
        self.assertEqual( int(abs(y2-b)*100000), 0 )

        # check the "LCC"-way: 4258 -> 3034 -> 4258
        l = 12.123450
        b = 51.234560
        t11 = Proj(fromEPSG=4258, toEPSG=3034)
        (x1, y1, z1) = t11.transform(l, b)
        t12 = Proj(fromEPSG=3034, toEPSG=4258)
        (x2, y2, z2) = t12.transform(x1, y1)
        self.assertEqual( int(abs(x2-l)*100000), 0 )
        self.assertEqual( int(abs(y2-b)*100000), 0 )

        # check the "ETRS89/32"-way: 4258 -> 3044 -> 4258
        l = 12.123450
        b = 51.234560
        t11 = Proj(fromEPSG=4258, toEPSG=3044)
        (x1, y1, z1) = t11.transform(l, b)
        t12 = Proj(fromEPSG=3044, toEPSG=4258)
        (x2, y2, z2) = t12.transform(x1, y1)
        self.assertEqual( int(abs(x2-l)*100000), 0 )
        self.assertEqual( int(abs(y2-b)*100000), 0 )

        # check the "ETRS89/33"-way: 4258 -> 5650 -> 4258
        l = 12.123450
        b = 51.234560
        t11 = Proj(fromEPSG=4258, toEPSG=5650)
        (x1, y1, z1) = t11.transform(l, b)
        t12 = Proj(fromEPSG=5650, toEPSG=4258)
        (x2, y2, z2) = t12.transform(x1, y1)
        self.assertEqual( int(abs(x2-l)*100000), 0 )
        self.assertEqual( int(abs(y2-b)*100000), 0 )
				
if __name__ == "__main__":
    unittest.main()
