""" Coordinate transformation based on mathematic formulars of Maarten Hooijberg (NL).
    Derived from my C# - API's for geodetic coordinates and functions.
    Dipl.-Ing. (Geodesy) Michael Dreesmann (DE), (C)opyright 2008-2019
    License: MIT-License
    ----------------------------------------------------------------------------------------
    Version 1.0-4 Beta, at 27.10.2019
    1.0-4   MD  Updated the module test
    1.0-3   MD  Erweiterung um die Systeme LCC, GOOG, XYZ
    1.0-2   MD  Fehlerkorrektur    
    ----------------------------------------------------------------------------------------
    Usage:
		t = Proj(von=4258, nach=3035)
		(x, y) = t.transform(12.12, 51.23)
		(x, y) = t.transform(12.00, 51.00)
"""
import math as Math

class Proj:
    EPSGs = { 3035:("LAEA", 0),  3034:( "LCC", 0),
              3857:("GOOG", 0),  7912:( "XYZ", 0),
              4258:("GEOD", 0),  # the official basis-world geodetic system for the earth
              3044:("ETRS",32),  3045:("ETRS",33), 
              4647:("ETRS",32),  5650:("ETRS",33), 
             25832:("ETRS",32), 25833:("ETRS",33) }

    # official basis-parameter of the Earth (see: ITF), named GRS80
    rho = 180.0 / Math.pi
    a   = 6378137.0
    f   = 1.0 / 298.2572221008827
    e   = Math.sqrt(2.0*f - f*f)
    e2  = (2.0*f - f*f)

#   ---------------------------------------------------------
#   Definition der Instanzen - Jede Klasse eine Transformation
#   ---------------------------------------------------------
    def __init__(self, von=4258, nach=3035):
        """ Konstruktor: von EPSG:<n> nach EPSG:<n> """
        chk_von   = int(von)
        chk_nach  = int(nach)
        if chk_von  == 4326:
            print("Official Basis-CRS of the Earth = EPSG:4258")
        if chk_nach == 4326:
            print("Official Basis-CRS of the Earth = EPSG:4258")
        if int(chk_von)  in self.EPSGs: self.von  = int(chk_von)
        if int(chk_nach) in self.EPSGs: self.nach = int(chk_nach)

#   ---------------------------------------------------------
#   Transformiere eine Koordinate
#   ---------------------------------------------------------
    def transform(self, ost, nord, Z=0.0):
        """ Transformation der Koordinaten """
        dle_ost = float(ost)
        dle_nrd = float(nord)
        dle_z   = float(Z)

        (sys, zn) = self.EPSGs[self.von]
        if   sys == 'LAEA': 
            (x,y) = self.trnLAEA_GRS80(dle_ost, dle_nrd)
        elif sys == 'XYZ': 
            (x,y, z) = self.trnXYZ_GRS80(dle_ost, dle_nrd, dle_z)
        elif sys == 'LCC': 
            (x,y) = self.trnLCC_GRS80(dle_ost, dle_nrd)
        elif sys == 'GOOG': 
            (x,y) = self.trnGOOG_GRS80(dle_ost, dle_nrd)
        elif sys == 'ETRS':
            if self.von == 4647: dle_ost -= 32000000.0
            if self.von == 5650: dle_ost -= 33000000.0
            (x,y) = self.trnETRS_GRS80(dle_ost, dle_nrd, zn)
        else:    #  Geodätische Koordinaten auf dem GRS80-Ellipsoid
            x = dle_ost
            y = dle_nrd
        
        dle_z = 0.0
        (sys, zn) = self.EPSGs[self.nach]
        if   sys == 'LAEA': 
            (dle_ost, dle_nrd) = self.trnGRS80_LAEA(x, y)
        elif sys == 'LCC': 
            (dle_ost, dle_nrd) = self.trnGRS80_LCC(x, y)
        elif sys == 'XYZ': 
            (dle_ost, dle_nrd, dle_z) = self.trnGRS80_XYZ(x, y, 0.0)
        elif sys == 'GOOG': 
            (dle_ost, dle_nrd) = self.trnGRS80_GOOG(x, y)
        elif sys == 'ETRS':
            (dle_ost, dle_nrd) = self.trnGRS80_ETRS(x, y, zn)
            if self.nach == 4647: dle_ost += 32000000.0
            if self.nach == 5650: dle_ost += 33000000.0
        else:    #  Geodätische Koordinaten auf dem GRS80-Ellipsoid
            dle_ost = x
            dle_nrd = y

        return (dle_ost, dle_nrd, dle_z)

#   ---------------------------------------------------------
#   3D - Transformation auf dem GRS80-Rotationsellipsoid
#   ---------------------------------------------------------
    def trnXYZ_GRS80(self, X, Y, Z):
        l = Math.atan(Y/X)
        b0 = Math.atan(Z/(1 - self.e*self.e)/Math.sqrt(X*X + Y*Y))
        N0 = self.a / (Math.sqrt(1 - self.e*self.e*Math.sin(b0)*Math.sin(b0)))
        if (b0*self.rho)<45.0:
            h0 = Math.sqrt(X*X + Y*Y)/Math.cos(b0) - N0
        else:
            h0 = Z / Math.sin(b0) - (1 - self.e*self.e)*N0
        b = Math.atan(Z/Math.sqrt(X*X + Y*Y)/(1 - (self.e*self.e*N0)/(N0+h0)))
        l *= self.rho
        b *= self.rho
        return (l, b, h0)

    def trnGRS80_XYZ(self, l, b, h):
        N = self.a / (Math.sqrt(1 - self.e*self.e*Math.sin(b/self.rho)*Math.sin(b/self.rho)))
        X = (N+h) * Math.cos(b/self.rho) * Math.cos(l/self.rho)
        Y = (N+h) * Math.cos(b/self.rho) * Math.sin(l/self.rho)
        Z = (N* (1 - self.e*self.e)+h) * Math.sin(b/self.rho)    
        return (X, Y, Z)

#   ---------------------------------------------------------
#   LCC - Transformation auf dem GRS80-Rotationsellipsoid
#   ---------------------------------------------------------
    def trnLCC_GRS80(self, ee, nn):
        ltl = 35.0 / self.rho
        ltu = 65.0 / self.rho
        ltc = 52.0 / self.rho
        lgc = 10.0 / self.rho
        e0 = 4000000.0
        n0 = 2800000.0
        e1 = e2 = 0.0

        e2 = self.e
        e1 = Math.sqrt(e2)

        f0 = f2 = f4 = f6 = f8 = 0.0

        f0 = e2 * (1.0 + e2 / 6.0 * (-1.0 + e2 * (0.2 + e2 / 84.0 * (31.0 / 5.0 + 5.0 * e2))))
        f2 = e2 * e2 * (7.0 / 6.0 + e2 / 5.0 * (-4.5 + e2 / 7.0 * (13.0 - 101.0 / 36.0 * e2)))
        f4 = e2 * e2 * e2 * (28.0 / 15.0 + 1.0 / 56.0 * e2 * (-467.0 / 3.0 + 117.0 * e2))
        f6 = e2 * e2 * e2 * e2 / 45.0 * (4279.0 / 28.0 - 344.0 * e2)
        f8 = e2 * e2 * e2 * e2 * e2 * 2087.0 / 315.0

        sl = qu = wu = ql = wl = qo = wo = q0 = w0 = 0.0

        sl = Math.sin(ltl)
        ql = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        wl = Math.sqrt(1.0 - (e2 * sl * sl))

        sl = Math.sin(ltu)
        qu = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        wu = Math.sqrt(1.0 - (e2 * sl * sl))

        sl = Math.sin(ltc)
        qo = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        wo = Math.sqrt(1.0 - (e2 * sl * sl))

        ss = m0 = k = r0 = 0.0

        ss = (Math.log(wu * Math.cos(ltl) / (wl * Math.cos(ltu)))) / (qu - ql)
        m0 = Math.atan(ss / Math.sqrt(1.0 - (ss * ss)))

        sl = Math.sin(m0)
        q0 = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        w0 = Math.sqrt(1.0 - (e2 * sl * sl))

        k  = self.a * Math.cos(ltl) * Math.exp(ql * ss) / (wl * ss)
        r0 = k / Math.exp(q0 * ss)

        rb = rc = 0.0

        rb = k / Math.exp(qo * ss)
        rc = n0 + rb
        # Ende: allgemeiner Teil

        nd = ed = lda1 = lon1 = qd = q = x = cx = cc = lat1 = 0.0

        nd = rc - nn
        ed = ee - e0
        lda1 = Math.atan(ed / nd)
        lon1 = lgc + lda1 / ss

        qd = Math.sqrt(nd * nd + ed * ed)
        q  = Math.log(k / qd) / ss
        x  = 2 * Math.atan((Math.exp(q) - 1.0) / (Math.exp(q) + 1.0))
        cx = Math.cos(x)
        cc = cx * cx
        lat1 = x + Math.sqrt(1.0 - cc) * cx * (f0 + cc * (f2 + cc * (f4 + cc * (f6 + f8 * cc))))

        wk = k1 = 0.0

        sl = Math.sin(lat1)
        wk = Math.sqrt(1.0 - (e2 * sl * sl))
        k1 = wk * qd * abs(ss) / (self.a * Math.cos(lat1))

        lat = lat1 * self.rho
        lon = lon1 * self.rho
        return (lon, lat)

    def trnGRS80_LCC(self, l, b):
        ltl = 35.0 / self.rho
        ltu = 65.0 / self.rho
        ltc = 52.0 / self.rho
        lgc = 10.0 / self.rho
        e0 = 4000000.0
        n0 = 2800000.0
        e1 = e2 = 0.0

        e2 = self.e
        e1 = Math.sqrt(e2)

        sl = qu = wu = ql = wl = qo = wo = q0 = w0 = 0.0

        sl = Math.sin(ltl)
        ql = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        wl = Math.sqrt(1.0 - (e2 * sl * sl))

        sl = Math.sin(ltu)
        qu = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        wu = Math.sqrt(1.0 - (e2 * sl * sl))

        sl = Math.sin(ltc)
        qo = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        wo = Math.sqrt(1.0 - (e2 * sl * sl))

        ss = m0 = k = r0 = 0.0

        ss = (Math.log(wu * Math.cos(ltl) / (wl * Math.cos(ltu)))) / (qu - ql)
        m0 = Math.atan(ss / Math.sqrt(1.0 - (ss * ss)))

        sl = Math.sin(m0)
        q0 = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        w0 = Math.sqrt(1.0 - (e2 * sl * sl))

        k  = self.a * Math.cos(ltl) * Math.exp(ql * ss) / (wl * ss)
        r0 = k / Math.exp(q0 * ss)

        rb = rc = 0.0

        rb = k / Math.exp(qo * ss)
        rc = n0 + rb
        # Ende: allgemeiner Teil

        qi = wi = ri = lda = 0.0

        sl = Math.sin(b/self.rho)
        qi = (Math.log((1.0 + sl) / (1.0 - sl)) - e1 * Math.log((1.0 + (e1 * sl)) / (1.0 - (e1 * sl)))) / 2.0
        wi = Math.sqrt(1.0 - (e2 * sl * sl))

        ri = k / Math.exp(qi * ss)
        lda = (lgc - (l/self.rho)) * ss

        y = rc - (ri * Math.cos(lda))
        x = e0 - (ri * Math.sin(lda))
        return (x, y)

#   ---------------------------------------------------------
#   Mercator-Transformation auf dem GRS80-Rotationsellipsoid
#   ---------------------------------------------------------
    def trnGOOG_GRS80(self, x, y):
        o = 2.0 * Math.pi * self.a / 2.0
        l = (x/o) * 180.0
        b = (y/o) * 180.0
        b = self.rho *(2.0 * Math.atan( Math.exp(b/self.rho)) - Math.pi/2.0)
        return (l, b)

    def trnGRS80_GOOG(self, l, b):
        o = 2.0 * Math.pi * self.a / 2.0
        x = l * o / 180.0
        y = Math.log( Math.tan((90.0+b)*Math.pi/360.0)) * self.rho
        y *= o / 180.0
        return (x, y)

#   ---------------------------------------------------------
#   LAEA - Transformation auf dem GRS80-Rotationsellipsoid
#   ---------------------------------------------------------
    def trnLAEA_GRS80(self, xi, yi):
        """ Transformation: Lambert LAEA (Europa) nach GRS80 Lamda, Phi (Erde) """
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

        e2 = self.e2
        e = Math.sqrt(e2)
        e4 = Math.pow(e, 4.0)
        e6 = Math.pow(e, 6.0)
        l0 = 10.0
        b0 = 52.0
        x0 = 4321000.0
        y0 = 3210000.0
        sinb0 = Math.sin(b0 / self.rho)
        cosb0 = Math.cos(b0 / self.rho)
        q0 = (1.0 - e2) * ((sinb0 / (1.0 - (e2 * sinb0 * sinb0))) - (1.0 / (2.0 * e)) * Math.log((1.0 - (e * sinb0)) / (1.0 + (e * sinb0))))
        qp = (1.0 - e2) * ((1.0 / (1.0 - e2)) - (1.0 / (2 * e) * Math.log((1.0 - e) / (1.0 + e))))
        bt0 = Math.asin(q0 / qp)
        cosbt0 = Math.cos(bt0)
        sinbt0 = Math.sin(bt0)
        rq = self.a * Math.sqrt(qp / 2.0)
        dd = (self.a * cosb0) / (Math.sqrt(1.0 - (e2 * sinb0 * sinb0)) * (rq * cosbt0))

        p1 = (xi - x0) / dd
        p2 = dd * (yi - y0)
        pp = Math.sqrt((p1 * p1) + (p2 * p2))
        cc = 2.0 * Math.asin(pp / (2.0 * rq))
        p3 = (dd * (yi - y0) * Math.sin(cc) * Math.cos(bt0)) / pp
        b1 = Math.asin(Math.cos(cc) * Math.sin(bt0) + p3)

        t1 = (e2 / 3.0) + (e4 * 31.0 / 180.0) + (e6 * 517.0 / 5040.0)
        t2 = (e4 * 23.0 / 360.0) + (e6 * 251.0 / 3780.0)
        t3 = (e6 * 761.0 / 45360.0)

        b = (b1 + (t1 * Math.sin(2.0 * b1)) + (t2 * Math.sin(4.0 * b1)) + (t3 * Math.sin(6.0 * b1))) * self.rho
        l = l0 + Math.atan((Math.sin(cc) * (xi - x0)) / ((dd * pp * Math.cos(bt0) * Math.cos(cc)) - (dd * dd * (yi - y0) * Math.sin(bt0) * Math.sin(cc)))) * self.rho
        return (l, b)

    def trnGRS80_LAEA(self, l, b):
        """ Transformation: GRS80 Lamda, Phi (Erde) nach Lambert LAEA (Europa) """
        e = e2 = 0.0
        l0 = b0 = 0.0
        x0 = y0 = 0.0
        sinb = sinb0 = q = q0 = qp = cosb0 = sinll = cosll = bt = bt0 = rq = bb = dd = 0.0
        cosbt = sinbt = cosbt0 = sinbt0 = 0.0

        x = 0.0
        y = 0.0

        # Ist Koordinate in Europa?
        if (l < -32.0 or l > 41.0): return (0.0, 0.0)
        if (b < 27.0 or b > 82.0): return (0.0, 0.0)

        e2 = self.e2
        e = Math.sqrt(e2)
        l0 = 10.0
        b0 = 52.0
        x0 = 4321000.0
        y0 = 3210000.0
        sinb = Math.sin(b / self.rho)
        sinb0 = Math.sin(b0 / self.rho)
        cosb0 = Math.cos(b0 / self.rho)
        sinll = Math.sin((l - l0) / self.rho)
        cosll = Math.cos((l - l0) / self.rho)
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
        rq = self.a * Math.sqrt(qp / 2.0)
        dd = (self.a * cosb0) / (Math.sqrt(1.0 - (e2 * sinb0 * sinb0)) * (rq * cosbt0))
        bb = rq * Math.sqrt(2.0 / (1.0 + (sinbt0 * sinbt) + (cosbt0 * cosbt * cosll)))

        y += (bb / dd * ((cosbt0 * sinbt) - (sinbt0 * cosbt * cosll)))
        x += (bb * dd * cosbt * sinll)
        return (x, y)

#   ---------------------------------------------------------
#   UTM - Transformation auf dem GRS80-Rotationsellipsoid
#   ---------------------------------------------------------
    def trnETRS_GRS80(self, xi, yi, z):
        """ Transformation: UTM/ETRS89 Zone <z> (Europa) nach GRS80 Lamda, Phi (Erde) """
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

        e2 = self.e2
        e = Math.sqrt(e2)
        ei2 = e2 / (1.0 - e2)
        ei4 = ei2 * ei2
        ei6 = ei4 * ei2
        ei8 = ei6 * ei2
        b0 = 0.0
        x0 = 500000.0
        y0 = 0.0
        c  = self.a / (Math.sqrt(1.0 - e2))
        n  = self.f / (2.0 - self.f)
        n2 = n * n
        r  = self.a * (1.0 + n2 / 4) / (1.0 + n)
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
        rrf = (k0 * self.a) / Math.sqrt(1.0 - e2 * sinb2)
        q  = (xi - x0) / rrf
        q2 = q * q
        q3 = q * q * q
        t  = Math.tan(bf)
        t2 = t * t
        n3 = ei2 * cosb2
        b2 = t * (1.0 + n3) / -2.0
        b4 = (5.0 + 3.0 * t2 + n3 * (1.0 - 9.0 * t2) - 4.0 * n3 * n3) / -12.0
        b6 = (61.0 + 90.0 * t2 + 45.0 * t2 * t2 + n3 * (46.0 - 252.0 * t2 - 90.0 * t2 * t2)) / 360.0
        b3 = (1.0 + 2.0 * t2 + n3) / -6.0
        b5 = (5.0 + 28.0 * t2 + 24.0 * t2 * t2 + n3 * (6.0 + 8.0 * t2)) / 120.0
        b7 = (61.0 + 662.0 * t2 + 1320.0 * t2 * t2 + 720.0 * t2 * t2 * t2) / -5040.0

        b = bf + b2 * q2 * (1.0 + q2 * (b4 + b6 * q2))
        b *= self.rho
        gl = q * (1.0 + q2 * (b3 + q2 * (b5 + b6 * q3)))
        l0 = (z * 6.0 - 183.0)
        l = l0 + (gl / cosb) * self.rho
        return (l, b)


    def trnGRS80_ETRS(self, l, b, z):
        """ Transformation: UTM/ETRS89 Zone <z> (Europa) nach GRS80 Lamda, Phi (Erde) """
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

        ll = l / self.rho
        bb = b / self.rho

        e2 = self.e2
        e = Math.sqrt(e2)
        ei2 = e2 / (1.0 - e2)
        ei4 = ei2 * ei2
        ei6 = ei4 * ei2
        ei8 = ei6 * ei2
        l0 = (z * 6.0 - 183.0) / self.rho
        b0 = 0.0
        x0 = 500000.0
        y0 = 0.0

        # Südhalbkugel ?
        if (bb < 0.0): y0 = 10000000.0

        c = self.a / (Math.sqrt(1.0 - e2))
        n = self.f / (2.0 - self.f)
        n2 = n * n
        r = self.a * (1.0 + n2 / 4) / (1.0 + n)

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
        rr = k0 * self.a / Math.sqrt(1.0 - e2 * sinb2)
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
    xi = 12.1
    yi = 51.2
    fk = 111111.111
    print("Calculate the compute blurring in [m]...")

    print("EPSG:{} > {:.3f} {:.3f} {:.3f}".format(4258,xi,yi,0.0))

    # Test1: GRS80 -> LAEA -> GRS80
    t11 = Proj(von=4258, nach=3035)
    (x1, y1, z1) = t11.transform(xi, yi)
    t12 = Proj(von=3035, nach=4258)
    (x2, y2, z2) = t12.transform(x1, y1)
    print("EPSG:{} > {:.3f} {:.3f} {:.3f}".format(3035,x1,y1,z1))
    print("           dLambda ~ {:.4f} m".format( abs(x2-xi)*fk ))
    print("              dPhi ~ {:.4f} m".format( abs(y2-yi)*fk ))

    # Test2: GRS80 -> ETRS89/32N -> GRS80
    t21 = Proj(von=4258, nach=3044)
    (x1, y1, z1) = t21.transform(xi, yi)
    t22 = Proj(von=3044, nach=4258)
    (x2, y2, z2) = t22.transform(x1, y1)
    print("EPSG:{} > {:.3f} {:.3f} {:.3f}".format(3044,x1,y1,z1))
    print("           dLambda ~ {:.4f} m".format( abs(x2-xi)*fk ))
    print("              dPhi ~ {:.4f} m".format( abs(y2-yi)*fk ))

    # Test3: GRS80 -> ETRS89/33N -> GRS80
    t31 = Proj(von=4258, nach=3045)
    (x1, y1, z1) = t31.transform(xi, yi)
    t32 = Proj(von=3045, nach=4258)
    (x2, y2, z2) = t32.transform(x1, y1)
    print("EPSG:{} > {:.3f} {:.3f} {:.3f}".format(3045,x1,y1,z1))
    print("           dLambda ~ {:.4f} m".format( abs(x2-xi)*fk ))
    print("              dPhi ~ {:.4f} m".format( abs(y2-yi)*fk ))

    # Test4: GRS80 -> LCC -> GRS80
    t41 = Proj(von=4258, nach=3034)
    (x1, y1, z1) = t41.transform(xi, yi)
    t42 = Proj(von=3034, nach=4258)
    (x2, y2, z2) = t42.transform(x1, y1)
    print("EPSG:{} > {:.3f} {:.3f} {:.3f}".format(3034,x1,y1,z1))
    print("           dLambda ~ {:.4f} m".format( abs(x2-xi)*fk ))
    print("              dPhi ~ {:.4f} m".format( abs(y2-yi)*fk ))

    # Test5: GRS80 -> GOOG -> GRS80
    t51 = Proj(von=4258, nach=3857)
    (x1, y1, z1) = t51.transform(xi, yi)
    t52 = Proj(von=3857, nach=4258)
    (x2, y2, z2) = t52.transform(x1, y1)
    print("EPSG:{} > {:.3f} {:.3f} {:.3f}".format(3857,x1,y1,z1))
    print("           dLambda ~ {:.4f} m".format( abs(x2-xi)*fk ))
    print("              dPhi ~ {:.4f} m".format( abs(y2-yi)*fk ))

    # Test6: GRS80 -> XYZ -> GRS80
    t61 = Proj(von=4258, nach=7912)
    (x1, y1, z1) = t61.transform(xi, yi)
    t62 = Proj(von=7912, nach=4258)
    (x2, y2, z2) = t62.transform(x1, y1, z1)
    print("EPSG:{} > {:.3f} {:.3f} {:.3f}".format(7912,x1,y1,z1))
    print("           dLambda ~ {:.4f} m".format( abs(x2-xi)*fk ))
    print("              dPhi ~ {:.4f} m".format( abs(y2-yi)*fk ))
