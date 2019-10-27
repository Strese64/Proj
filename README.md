# Proj
Map Projections for Germany is a Python Class.

## Fundamentals
The rotation ellipsoid GRS80 (Geodetic Reference System, at date 1980) is the basis for all projections in the world, except the U.S.A., they use WGS84, which is identical. The geodesist and surveyor in the world use the EPSG (European Petroleum Surveying Group) codes to identify the different CRS (Coordinate Reference Systems).

* EPSG:4258 (and EPSG:4326) for the geodetic coordinates (the angles longitude and latitude), world wide
* EPSG:3034 as the LCC Lambert-projection, for Europe
* EPSG:3035 as the LAEA Lambert-projection, for Europe
* EPSG:3044 (as well as EPSG:4647 and EPSG:25832) for ETRS89 coordinates of the 32nd UTM zone (here: Germany)
* EPSG:3045 (as well as EPSG:5650 and EPSG:25833) for ETRS89 coordinates of the 33rd UTM zone (here: Germany)
* EPSG:3038 - EPSG:3051 for whole Europe

Here we use only the number (EPSG-code). ETRS89 is the official coordinate system of all nations in the European Union.

## Usage

Just use the following example to transform from EPSG:4258 (GRS80) to EPSG:3035 (LAEA):

```
t = Proj(fromEPSG=4258, toEPSG=3035)
(x, y, z) = t.transform(12.12, 51.23)
(x, y, z) = t.transform(12.00, 51.00)
```

We use the parameter order east value and north value.

## Test automation

There is a module test and a unit test stored in this code.

## Licence

The base is the MIT-licence. Please add the following words in your code: "Python-Class Proj, developed by Michael Dreesmann (DE), based on the formula of Hooijberg/Annoni/Luzet/Gubler/Ihde/et al, v1.0 in 2019"

## Source

The mathematics based on:

* Maarten Hooijberg, Practical Geodesy (1997, Springer Verlag)
* Annoni Luzet Gubler & Ihde, Map Projections for Europe (2003, European Commission, Joint Research Center)

## Author

Dipl.-Ing. Michael Dreesmann made a surveying/geodesy-degree in 1991 at the university of applied science of Bochum. Since than he works as geomatics engineer as an software developer, at 2nd/3rd Level-Support, as a consultant and an innovation manager. He is one of the coauthors of the book "Geodateninfrastruktur, Grundlagen und Anwendungen" from 2004. From 2004 to 2007 he was the technical head of the SDI (Spatial Data Infrastructure) of the German state Brandenburg. From 2008 to 2011 he was a member of the Open Geospatial Consortium (OGC).
