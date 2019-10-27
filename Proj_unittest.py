#   -------------------------------------------------------------------
#   Proj    coordinate transformation class
#           Purpose Unittest for all 4 projection functions
#           Version 1.0-2 Beta, from 06.10.2019
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
        t11 = Proj(von=4258, nach=3035)
        (x1, y1) = t11.transform(l, b)
        t12 = Proj(von=3035, nach=4258)
        (x2, y2) = t12.transform(x1, y1)
        self.assertEqual( int(abs(x2-l)*100000), 0 )
        self.assertEqual( int(abs(y2-b)*100000), 0 )

        # check the "ETRS89/32"-way: 4258 -> 3044 -> 4258
        l = 12.123450
        b = 51.234560
        t11 = Proj(von=4258, nach=3044)
        (x1, y1) = t11.transform(l, b)
        t12 = Proj(von=3044, nach=4258)
        (x2, y2) = t12.transform(x1, y1)
        self.assertEqual( int(abs(x2-l)*100000), 0 )
        self.assertEqual( int(abs(y2-b)*100000), 0 )

        # check the "ETRS89/33"-way: 4258 -> 5650 -> 4258
        l = 12.123450
        b = 51.234560
        t11 = Proj(von=4258, nach=5650)
        (x1, y1) = t11.transform(l, b)
        t12 = Proj(von=5650, nach=4258)
        (x2, y2) = t12.transform(x1, y1)
        self.assertEqual( int(abs(x2-l)*100000), 0 )
        self.assertEqual( int(abs(y2-b)*100000), 0 )
				
if __name__ == "__main__":
    unittest.main()
