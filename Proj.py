""" Koordinatentransformation auf den mathematischen Grundlagen von Maarten Hooijberg (NL).
    Abgeleitet aus meiner C#-Bibliothek. Michael Dreesmann (DE), (C)opyright 2008-2019
    Lizenz: Open Source für private, freie und kommerzielle Anwendungen.
    Version 1.0 Beta, vom 05.10.2019
    ----------------------------------------------------------------------------------------
    Anwendung:
		t = Proj(von=4258, nach=3035)		# von EPSG:4258 nach EPSG:3035 (Lambert LAEA)
		(x, y) = t.transform(12.12, 51.23)
		(x, y) = t.transform(12.00, 51.00)
"""

import math as Math

class Proj:
    EPSGs = {   3035:("LAEA",0), 
                3044:("ETRS",32),  3045:("ETRS",33), 
                4647:("ETRS",32),  5650:("ETRS",32), 
                4258:("GEOD", 0),  4326:("GEOD", 0), 
               25832:("ETRS",32), 25833:("ETRS",33) }

    def __init__(self, von=4258, nach=3035):
        """ Konstruktor: von EPSG:<n> nach EPSG:<n> """
        chk_von  = -1
        chk_nach  = -1
        if int(von)  in self.EPSGs: chk_von  = int(von)
        if int(nach) in self.EPSGs: chk_nach = int(nach)
        self.von  = chk_von
        self.nach = chk_nach

    def transform(self, ost, nord):
        """ Transformation der Koordinaten """
        dle_ost = float(ost)
        dle_nrd = float(nord)

        (sys, zn) = self.EPSGs[self.von]
        if   sys == 'LAEA': 
            (x,y) = self.trnLAEA_GRS80(dle_ost, dle_nrd)
        elif sys == 'ETRS':
            if self.von == 4647: dle_ost -= 32000000.0
            if self.von == 5650: dle_ost -= 33000000.0
            (x,y) = self.trnETRS_GRS80(dle_ost, dle_nrd, zn)
        else:    #  Geodätische Koordinaten auf dem GRS80-Ellipsoid
            x = dle_ost
            y = dle_nrd
        
        (sys, zn) = self.EPSGs[self.nach]
        if   sys == 'LAEA': 
            (dle_ost, dle_nrd) = self.trnGRS80_LAEA(x, y)
        elif sys == 'ETRS':
            (dle_ost, dle_nrd) = self.trnGRS80_ETRS(x, y, zn)
            if self.nach == 4647: dle_ost += 32000000.0
            if self.nach == 5650: dle_ost += 33000000.0
        else:    #  Geodätische Koordinaten auf dem GRS80-Ellipsoid
            dle_ost = x
            dle_nrd = y

        return (dle_ost, dle_nrd)


    def trnLAEA_GRS80(self, xi, yi):
        """ Transformation: Lambert LAEA (Europa) nach GRS80 Lamda, Phi (Erde) """
        rho = 180.0 / 3.1415926535
        a = 6378137.0
        f = 1.0 / 298.2572221008827
        e = e2 = e4 = e6 = 0.0
        l0 = b0 = 0.0
        x0 = y0 = 0.0
        sinb0 =  q0 = qp = cosb0 = bt0 = rq = dd = 0.0
        cosbt0 = sinbt0 = 0.0
        p1 = p2 = p3 = pp = cc = b1 = t1 = t2 = t3 = 0.0
        l = 0.0
        b = 0.0

        # Ist Koordinate in Europa ?
        if (xi < 100000.0 or xi > 8000000.0): return (0.0,0.0)
        if (yi < 100000.0 or yi > 8000000.0): return (0.0,0.0)

        e2 = 2.0 * f - (f * f)
        e = Math.sqrt(e2)
        e4 = Math.pow(e, 4.0)
        e6 = Math.pow(e, 6.0)
        l0 = 10.0
        b0 = 52.0
        x0 = 4321000.0
        y0 = 3210000.0
        sinb0 = Math.sin(b0 / rho)
        cosb0 = Math.cos(b0 / rho)
        q0 = (1.0 - e2) * ((sinb0 / (1.0 - (e2 * sinb0 * sinb0))) - (1.0 / (2.0 * e)) * Math.log((1.0 - (e * sinb0)) / (1.0 + (e * sinb0))))
        qp = (1.0 - e2) * ((1.0 / (1.0 - e2)) - (1.0 / (2 * e) * Math.log((1.0 - e) / (1.0 + e))))
        bt0 = Math.asin(q0 / qp)
        cosbt0 = Math.cos(bt0)
        sinbt0 = Math.sin(bt0)
        rq = a * Math.sqrt(qp / 2.0)
        dd = (a * cosb0) / (Math.sqrt(1.0 - (e2 * sinb0 * sinb0)) * (rq * cosbt0))

        p1 = (xi - x0) / dd
        p2 = dd * (yi - y0)
        pp = Math.sqrt((p1 * p1) + (p2 * p2))
        cc = 2.0 * Math.asin(pp / (2.0 * rq))
        p3 = (dd * (yi - y0) * Math.sin(cc) * Math.cos(bt0)) / pp
        b1 = Math.asin(Math.cos(cc) * Math.sin(bt0) + p3)

        t1 = (e2 / 3.0) + (e4 * 31.0 / 180.0) + (e6 * 517.0 / 5040.0)
        t2 = (e4 * 23.0 / 360.0) + (e6 * 251.0 / 3780.0)
        t3 = (e6 * 761.0 / 45360.0)

        b = (b1 + (t1 * Math.sin(2.0 * b1)) + (t2 * Math.sin(4.0 * b1)) + (t3 * Math.sin(6.0 * b1))) * rho
        l = l0 + Math.atan((Math.sin(cc) * (xi - x0)) / ((dd * pp * Math.cos(bt0) * Math.cos(cc)) - (dd * dd * (yi - y0) * Math.sin(bt0) * Math.sin(cc)))) * rho
        return (l, b)


    def trnETRS_GRS80(self, xi, yi, z):
        """ Transformation: UTM/ETRS89 Zone <z> (Europa) nach GRS80 Lamda, Phi (Erde) """
        rho = 180.0 / 3.1415926535
        a = 6378137.0
        f = 1.0 / 298.2572221008827
        k0 = 0.9996
        e = e2 = ei2 = ei4 = ei6 = ei8 = 0.0
        l0 = b0 = 0.0
        x0 = y0 = 0.0
        sinb = cosb = sinb2 = cosb2 = 0.0
        n = n2 = c = r = u0 = u2 = u4 = u6 = v0 = v2 = v4 = v6 = w0 = s0 = rrf = t2 = n3 = a1 = a2 = a3 = a4 = a5 = a6 = a7 = gl = gl2 = 0.0
        bf = q = q2 = q3 = b2 = b3 = b4 = b5 = b6 = b7 = sinw = cosw = cosw2 = 0.0

        l = 0.0
        b = 0.0

        # Ist Koordinate in Europa?
        if (z < 26 or z > 39): return (0.0, 0.0)
        if (xi < 0.0 or xi > 1000000.0): return (0.0, 0.0)
        if (yi < 3000000.0 or yi > 7000000.0): return (0.0, 0.0)

        e2 = 2.0 * f - (f * f)
        e = Math.sqrt(e2)
        ei2 = e2 / (1.0 - e2)
        ei4 = ei2 * ei2
        ei6 = ei4 * ei2
        ei8 = ei6 * ei2
        b0 = 0.0
        x0 = 500000.0
        y0 = 0.0
        c = a / (Math.sqrt(1.0 - e2))
        n = f / (2.0 - f)
        n2 = n * n
        r = a * (1.0 + n2 / 4) / (1.0 + n)
        u0 = c * (((((11025.0 - (86625.0 * ei2 / 8.0)) * ei2 / 64.0 - 175.0) * ei2 / 4.0 + 45.0) * ei2 / 16.0 - 3.0) * ei2 / 4.0)
        u2 = c * ((((3675.0 - (17325.0 * ei2 / 4.0)) * ei2 / 256.0 - (175.0 / 12.0)) * ei2 + 15.0) * ei4 / 32.0)
        u4 = c * (735.0 * ei2 - 1493.0 / 2.0) * ei6 / 2048.0
        u6 = c * ((315.0 - 3465.0 * ei2 / 4.0) * ei8 / 1024.0)
        v0 = ((((16384.0 * ei2 - 11025.0) * ei2 / 64.0 + 175.0) * ei2 / 4.0 - 45.0) * ei2 / 16.0 + 3.0) * ei2 / 4.0
        v2 = (((19413.0 - 20464721.0 * ei2 / 120.0) * ei2 / 8.0 - 1477.0) * ei2 / 32.0 + 21.0) * ei4 / 32.0
        v4 = ((4737141.0 * ei2 / 28.0 - 17121.0) * ei2 / 32.0 + 151.0) * ei6 / 192.0
        v6 = (1097.0 - 427277.0 * ei2 / 35.0) * ei8 / 1024.0

        sinb = Math.sin(b0)
        cosb = Math.cos(b0)
        cosb2 = cosb * cosb
        w0 = b0 * r + sinb * cosb * (u0 + cosb2 * (u2 + cosb2 * (u4 + u6 * cosb2)))
        s0 = k0 * w0

        w0 = (yi - y0 + s0) / (k0 * r)
        sinw = Math.sin(w0)
        cosw = Math.cos(w0)
        cosw2 = cosw * cosw
        bf = w0 + (sinw * cosw) * (v0 + v2 * cosw2 + v4 * cosw2 * cosw2 + v6 * cosw2 * cosw2 * cosw2)
        sinb = Math.sin(bf)
        sinb2 = sinb * sinb
        cosb = Math.cos(bf)
        cosb2 = cosb * cosb
        rrf = (k0 * a) / Math.sqrt(1.0 - e2 * sinb2)
        q = (xi - x0) / rrf
        q2 = q * q
        q3 = q * q * q
        t = Math.tan(bf)
        t2 = t * t
        n3 = ei2 * cosb2
        b2 = t * (1.0 + n3) / -2.0
        b4 = (5.0 + 3.0 * t2 + n3 * (1.0 - 9.0 * t2) - 4.0 * n3 * n3) / -12.0
        b6 = (61.0 + 90.0 * t2 + 45.0 * t2 * t2 + n3 * (46.0 - 252.0 * t2 - 90.0 * t2 * t2)) / 360.0
        b3 = (1.0 + 2.0 * t2 + n3) / -6.0
        b5 = (5.0 + 28.0 * t2 + 24.0 * t2 * t2 + n3 * (6.0 + 8.0 * t2)) / 120.0
        b7 = (61.0 + 662.0 * t2 + 1320.0 * t2 * t2 + 720.0 * t2 * t2 * t2) / -5040.0

        b = bf + b2 * q2 * (1.0 + q2 * (b4 + b6 * q2))
        b *= rho
        gl = q * (1.0 + q2 * (b3 + q2 * (b5 + b6 * q3)))
        l0 = (z * 6.0 - 183.0)
        l = l0 + (gl / cosb)*rho
        return (l, b)


    def trnGRS80_LAEA(self, l, b):
        """ Transformation: GRS80 Lamda, Phi (Erde) nach Lambert LAEA (Europa) """
        rho = 180.0 / 3.1415926535
        a = 6378137.0
        f = 1.0 / 298.2572221008827
        e = e2 = 0.0
        l0 = b0 = 0.0
        x0 = y0 = 0.0
        sinb = sinb0 = q = q0 = qp = cosb0 = sinll = cosll = bt = bt0 = rq = bb = dd = 0.0
        cosbt = sinbt = cosbt0 = sinbt0 = 0.0

        x = 0.0
        y = 0.0

        # Ist Koordinate in Europa?
        if (l < -32.0 or l > 41.0): return false
        if (b < 27.0 or b > 82.0): return false

        e2 = 2.0 * f - (f * f)
        e = Math.sqrt(e2)
        l0 = 10.0
        b0 = 52.0
        x0 = 4321000.0
        y0 = 3210000.0
        sinb = Math.sin(b / rho)
        sinb0 = Math.sin(b0 / rho)
        cosb0 = Math.cos(b0 / rho)
        sinll = Math.sin((l - l0) / rho)
        cosll = Math.cos((l - l0) / rho)
        x = x0
        y = y0

        q = (1.0 - e2) * ((sinb / (1.0 - (e2 * sinb * sinb))) - (1.0 / (2.0 * e)) * Math.log((1.0 - (e * sinb)) / (1.0 + (e * sinb))))
        q0 = (1.0 - e2) * ((sinb0 / (1.0 - (e2 * sinb0 * sinb0))) - (1.0 / (2.0 * e)) * Math.log((1.0 - (e * sinb0)) / (1.0 + (e * sinb0))))
        qp = (1.0 - e2) * ((1.0 / (1.0 - e2)) - (1.0 / (2 * e) * Math.log((1.0 - e) / (1.0 + e))))
        bt = Math.asin(q / qp)
        bt0 = Math.asin(q0 / qp)
        cosbt = Math.cos(bt)
        cosbt0 = Math.cos(bt0)
        sinbt = Math.sin(bt)
        sinbt0 = Math.sin(bt0)
        rq = a * Math.sqrt(qp / 2.0)
        dd = (a * cosb0) / (Math.sqrt(1.0 - (e2 * sinb0 * sinb0)) * (rq * cosbt0))
        bb = rq * Math.sqrt(2.0 / (1.0 + (sinbt0 * sinbt) + (cosbt0 * cosbt * cosll)))

        y += (bb / dd * ((cosbt0 * sinbt) - (sinbt0 * cosbt * cosll)))
        x += (bb * dd * cosbt * sinll)
        return (x, y)


    def trnGRS80_ETRS(self, l, b, z):
        """ Transformation: UTM/ETRS89 Zone <z> (Europa) nach GRS80 Lamda, Phi (Erde) """
        rho = 180.0 / 3.1415926535
        a = 6378137.0
        f = 1.0 / 298.2572221008827
        k0 = 0.9996
        e = e2 = ei2 = ei4 = ei6 = ei8 = 0.0
        ll = bb = l0 = b0 = 0.0
        x0 = y0 = 0.0
        sinb = cosb = sinb0 = cosb0 = 0.0
        n = n2 = c = r = u0 = u2 = u4 = u6 = w0 = s0 = s = rr = t2 = n3 = a1 = a2 = a3 = a4 = a5 = a6 = a7 = gl = gl2 = 0.0

        x = 0.0
        y = 0.0

        # Ist Koordinate in Europa?
        if (z < 26 or z > 39): return (0.0, 0.0)
        if ((l + 186.0) / 6.0 < (z - 0.8)): return (0.0, 0.0)
        if ((l + 186.0) / 6.0 > (z + 1.8)): return (0.0, 0.0)
        if (l < -32.0 or l > 41.0): return (0.0, 0.0)
        if (b < 27.0 or b > 82.0): return (0.0, 0.0)

        ll = l / rho
        bb = b / rho

        e2 = 2.0 * f - (f * f)
        e = Math.sqrt(e2)
        ei2 = e2 / (1.0 - e2)
        ei4 = ei2 * ei2
        ei6 = ei4 * ei2
        ei8 = ei6 * ei2
        l0 = (z * 6.0 - 183.0) / rho
        b0 = 0.0
        x0 = 500000.0
        y0 = 0.0

        # Südhalbkugel ?
        if (bb < 0.0): y0 = 10000000.0

        c = a / (Math.sqrt(1.0 - e2))
        n = f / (2.0 - f)
        n2 = n * n
        r = a * (1.0 + n2 / 4) / (1.0 + n)

        sinb = Math.sin(b0)
        cosb = Math.cos(b0)
        cosb2 = cosb * cosb
        u0 = c * (((((11025.0 - (86625.0 * ei2 / 8.0)) * ei2 / 64.0 - 175.0) * ei2 / 4.0 + 45.0) * ei2 / 16.0 - 3.0) * ei2 / 4.0)
        u2 = c * ((((3675.0 - (17325.0 * ei2 / 4.0)) * ei2 / 256.0 - (175.0 / 12.0)) * ei2 + 15.0) * ei4 / 32.0)
        u4 = c * (735.0 * ei2 - 1493.0 / 2.0) * ei6 / 2048.0
        u6 = c * ((315.0 - 3465.0 * ei2 / 4.0) * ei8 / 1024.0)

        w0 = b0 * r + sinb * cosb * (u0 + cosb2 * (u2 + cosb2 * (u4 + u6 * cosb2)))
        s0 = k0 * w0

        sinb = Math.sin(bb)
        sinb2 = sinb * sinb
        cosb = Math.cos(bb)
        cosb2 = cosb * cosb

        gl = (ll - l0) * cosb
        gl2 = gl * gl
        t2 = Math.tan(bb) * Math.tan(bb)
        n3 = ei2 * cosb2
        w0 = bb * r + sinb * cosb * (u0 + cosb2 * (u2 + cosb2 * (u4 + u6 * cosb2)))
        s = k0 * w0
        rr = k0 * a / Math.sqrt(1.0 - e2 * sinb2)
        a1 = rr
        a3 = (1.0 - t2 + n3) / 6.0
        a5 = (5.0 - 18 * t2 + t2 * t2 + n3 * (14.0 - 58.0 * t2)) / 120.0
        a7 = (61.0 - 479 * t2 + 179.0 * t2 * t2 - t2 * t2 * t2) / 5040.0
        a2 = rr * Math.tan(bb) / 2.0
        a4 = (5.0 - t2 + n3 * (9.0 + 4.0 * n3)) / 12.0
        a6 = (61.0 - 58.0 * t2 + t2 * t2 + n3 * (270.0 - 330.0 * t2)) / 360.0

        x = x0 + a1 * gl * (1.0 + gl2 * (a3 + gl2 * (a5 + a7 * gl2)))
        y = s - s0 + y0 + a2 * gl2 * (1.0 + gl2 * (a4 + a6 * gl2))
        return (x, y)


if __name__ == "__main__":
    # Irgendwo an der A38 bei Sössen, Nahe der Grenze Zone 32 und Zone 33
    xi = 12.12345
    yi = 51.23456
    # Die Rechenunschärfe bitte in Meter am Äquator ausrechnen
    fk = 40000000.0 / 360.0
    print("Berechnung einer Rechenunschärfe ...")

    # Test1: GRS80 -> LAEA -> GRS80
    t11 = Proj(von=4258, nach=3035)
    (x1, y1) = t11.transform(xi, yi)
    t12 = Proj(von=3035, nach=4326)
    (x2, y2) = t12.transform(x1, y1)
    print("Test1 :: dLambda < {:.5f} m".format( abs(x2-xi)*fk ))
    print("            dPhi < {:.5f} m".format( abs(y2-yi)*fk ))

    # Test2: GRS80 -> ETRS89/32N -> GRS80
    t21 = Proj(von=4258, nach=3044)
    (x1, y1) = t21.transform(xi, yi)
    t22 = Proj(von=3044, nach=4326)
    (x2, y2) = t22.transform(x1, y1)
    print("Test2 :: dLambda < {:.5f} m".format( abs(x2-xi)*fk ))
    print("            dPhi < {:.5f} m".format( abs(y2-yi)*fk ))

    # Test3: GRS80 -> ETRS89/33N -> GRS80
    t31 = Proj(von=4258, nach=3045)
    (x1, y1) = t31.transform(xi, yi)
    t32 = Proj(von=3045, nach=4326)
    (x2, y2) = t32.transform(x1, y1)
    print("Test3 :: dLambda < {:.5f} m".format( abs(x2-xi)*fk ))
    print("            dPhi < {:.5f} m".format( abs(y2-yi)*fk ))
