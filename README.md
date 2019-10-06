# Proj
Map Projections for Germany is a Python Class.

## Grundlagen
Der Rotationsellipsoid GRS80 (als auch das spätere identische WGS84) ist die Basis für alle Projektionen bzw. Abbildungen. Zur Kodierung der verschiedenen geodätischen Koordinatensysteme wird auf die EPSG-Mimik zurückgegriffen. Hierbei steht EPSG für die European Petroleum Surveying Group, eine Vereinigung europäischer Geodäten (Vermessungsingenieure).

* EPSG:4258 (sowie EPSG:4326) als Reinform der geodätischen Koordinaten (die Winkel Lambda=Länge und Phi=Breite), für die weltweite Nutzung
* EPSG:3035 als Azimuthale Lambert-Projektion, ideal für eine Karten-Darstellung im europaweiten Kontext
* EPSG:3044 (sowie EPSG:4647 und EPSG:25832) für ETRS89-Koordinaten der 32. Zone als UTM-Abbildung mit variablen Maßstab
* EPSG:3045 (sowie EPSG:5650 und EPSG:25833) für ETRS89-Koordinaten der 33. Zone als UTM-Abbildung mit variablen Maßstab

Wir verwenden nur die Zahl und lassen das Kürzel "EPSG:" weg. Das ETRS89 in der UTM-Abbildung ist das offiziellen Koordinaten-Referenz-System (engl. CRS) der Europäischen Union. Es findet zentral im Liegenschaftskataster (den Geobasisdaten von Europa) Verwendung.

## Anwendung

Folgendes kleines Beispiel projeziert geodätische Koordinaten ins Azimuthale Lambert-System:

```
t = Proj(von=4258, nach=3035)
(x, y) = t.transform(12.12, 51.23)
(x, y) = t.transform(12.00, 51.00)
```

Wir verwenden in der Schnittstelle immer den Ostwert/Lambda als ersten Parameter und den Nordwert/Phi als zweiten Parameter.

## Testautomation

Wie in Python gewohnt testet sich die Schnittstelle selbst, wenn man die Klasse direkt (bspw. mit der Software VS-Code) ausführt.

## Lizenz

Jeder darf diese Python-Klasse für private, freie oder kommerzielle Anwendungen nutzen, wenn er den Inhalt unverändert belässt und mich wie folgt als Autor der Klasse benennt: "Python-Klasse Proj, entwickelt von Michael Dreesmann (DE), auf Basis von Hooijberg/Annoni/Luzet/Gubler/Ihde/et al, v1.0 in 2019"

## Quellen

Die Mathematik aus dieser Klasse wurde folgenden Quellen entnommen:

* Maarten Hooijberg, Practical Geodesy (1997, Springer Verlag)
* Annoni Luzet Gubler & Ihde, Map Projections for Europe (2003, European Commission, Joint Research Center)

## Zum Autor

Dipl.-Ing. Michael Dreesmann hat Ende der 80'er Jahre in Bochum Vermessungswesen (Geodäsie) studiert. Seitdem arbeitet er als Geoinformatiker in der Softwareentwicklung, im 2nd/3rd Level-Support oder als Consultant und Innovationsmanager in Wirtschaftsverbänden. Er ist ein Mitautor des Buches "Geodateninfrastruktur, Grundlagen und Anwendungen" von 2004 und war 1. Geschäftsstellenleiter der Geodateninfrastruktur des Landes Brandenburg.
