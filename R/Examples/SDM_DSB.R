###############################################################################################
## Radionica - Uvod u R i ekološku analizu | Mirza Cengic (mirzaceng@gmail.com) | 27.12.2016 ##
## Skripta sadrži osnove obrade prostornih podataka i osnove modeliranja distribucije vrsta  ##
## ############################################################################################
# --------------------------------
# Kod u RStudiu možemo pokrenuti ako pritisnemo Ctr + Enter, i linija
# na kojoj je kursor trenutno će biti pokrenuta, ili dio koda koje je trenutno
# selektovan!

# Ova funkcija će ukloniti sve objekte iz memorije u R-u. Dobra je praksa
# koristiti ovu funkciju na početku skripti. Ukoliko niste spasili pojedine
# objekte koje ne mozete reprodukovati skriptom, budite pažljivi sa ovom
# komandom.

rm(list = ls())

# Funkciju install.packages() je potrebno pokrenuti samo jednom za svaki paket nakon instaliranja ili update-ovanja R-a.
# Nakon što jednom instaliramo pojedine pakete, svaki slijedeći put potrebno je samo ih učitati u memoriju funkcijom library()
install.packages(c("rgdal", "sp", "raster"))

# Ova tri paketa će pretvoriti R u jako sposoban i moćan GIS software. Ova tri
# paketa su samo mali dio R GIS univerzuma. Za više informacija o dostupnim GIS
# paketima za R možete posjetiti
# https://cran.r-project.org/web/views/Spatial.html Postoje i brojni tutoriali
# dostupni online. Jedan od njih koji prikazuje neke osnovne mogućnosti R-a je
# https://pakillo.github.io/R-GIS-tutorial/
library(rgdal)
library(sp)
library(raster)

# Bitan korak je podesiti "working directory" (radni direktorij). Napravite
# folder u kojem će biti svi file-ovi potrebni za ovu skriptu. Put do
# direktorija treba da bude između navodnika (" "), i ukoliko koristite Windows
# kao operativni sistem, direktorije je potrebno odvojiti sa "/" ili sa "\\".
# Ukoliko koristite Linux ili Mac računar, separator za foldere je "/". Trenutni
# working directory možemo provjeriti komandom getwd(). 

# Put ispod je dat za primjer.
setwd("C:/Users/MirzaC/Desktop/DSB_workshop/R_directory/")
# -----------------------------------------------------------------------------------------------------------
# Da dobijemo lokaciju file-a kojeg želimo učitati u R, korisna je funkcija
# file.choose(). Ova funkcija će otvoriti novi prozor, i nakon što izaberemo
# file, njegova lokacija koju možemo prekopirati će biti ispisana u konzoli.
# Nađite file Salamandra_atra_loc_BiH.csv i prekopirajte put do file-a između zagrada u funkciji read.csv() ispod.
# HINT: put mora biti između znakova navodnika (" "). Ukoliko koristite RStudio, i file se nalazi u working direktoriju,
# dovoljno je da započnemo pisati ime file-a i pritisnemo dugme Tab za auto-fill (npr. read.csv("Sal + Tab))

file.choose()

# Funkcijom read.csv() učitavamo csv file-ove (comma separated value). CSV
# dokumenti sadrže podatke u vidu tabela koje su odvojene specijalnim karakterom
# - separatorom. Separator je najčešće zarez - "," (eng. comma), ali se često
# koriste i drugi karakteri.
atra_loc <- read.csv()


# Funkcija head() nam prikazuje prvih 6 linija objekta. Možemo također koristiti i funkciju tail() za ispis posljednjih 6 linija.
head(atra_loc)


# Korisne funkcije za pregled sadržaja objekta: head(), tail(), str(), names(), unique(), table(), summary()

str(atra_loc)
table(atra_loc$Location)

# Ovaj data.frame sadrži 10 varijabli. Koje sve tipove podataka možemo pronaći
# ovdje? Da li tip podatka prikladan svakoj varijabli? Treba li nešto izmjeniti?
# HINT: Pogledajte help dokument za funkciju read.csv(), i pročitajte opis za
# argument "stringsAsFactors".

?read.csv()

# -------------------------------------------------------------------------
# Prvo ćemo definisati da su kolone X i Y koordinate u geografskom prostoru.
# HINT: čest problem u ovoj fazi (ili prilikom unosa podataka) je zamjena X i Y
# (ili lat/long) koordinata. Ako znate da imate podatke za BiH, a na karti tačke
# se nalaze na području Saudijske Arabije, rješenje je zamijeniti X i Y
# koordinate.
browseURL("http://www.gps-coordinates.net/")

# Ukoliko dobijete grešku "Error in coordinates(atra_loc) <- ~Y + X : could not
# find function "coordinates<-", najvjerovatnije niste učitali paket sp
# (funkcija library(sp) iznad).

coordinates(atra_loc)<-~Y+X

# Definisanje koordinatnog sistema (u ovom slučaju u GPS-u sa kojim su
# zabilježene tačke izabran je WGS84 referentni sistem, pa ćemo definisati isti
# referentni sistem za naše tačke).
browseURL("https://www.google.com/search?q=wgs84+epsg&oq=wgs84+es&aqs=chrome.1.69i57j0l5.3564j0j4&sourceid=chrome&ie=UTF-8")

proj4string(atra_loc) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

plot(atra_loc)
# -----------------------------------------------------

# Postoji veliki broj paketa i fukcija u R-u koje nam mogu koristiti za skidanje
# raznih podataka. Jedan primjer je funkcija getData() iz paketa raster. Pomoću
# ove funkcije možemo dobiti podatke o administrativnim regijama, o klimi,
# nadmorskoj visini itd.
?getData()

# Funkcija getData() skida podatke sa interneta, i potrebno je da jednom skinemo
# file-ove nakon čega ih sa funcijom writeOGR možemo spasiti na disk. 
# Mijenjanjem argumenta "level" možemo skinuti granice "nižih" administrativnih
# regija (za BiH entiteti, kantoni i sl.). Dalje u skripti ćete naći komande
# gdje spašavamo file i gdje ga dalje opet učitavamo. Nije potrebno učitavati
# file-ove nakon što su spašeni jer će ostati u memoriji sve dok ne ugasimo R.
# Ukoliko restartujemo R potrebno je učitati spašene file-ove, ili je potrebno
# spasiti trenutnu sesiju (u RStudio Session -> Save workspace as).

BiH_boundary <- getData(name = "GADM", country = "BA", level = 0)
writeOGR(obj = BiH_boundary, dsn = ".", layer = "BiH_boundary", driver = "ESRI Shapefile")

# Pomoću funkcije readOGR() iz paketa rgdal možemo učitati prostorne podatke u
# različitim formatima kao što su kml, gpx ili shapefile formatima. Mali detalji
# su bitni da ova funkcija radi kako treba. Npr. ukoliko postoje separatori za
# foldere na kraju puta direktorija ("\\DSB_workshop_R\\"), funkcija će izbaciti
# grešku. Za prirodnije učitavanje shapefile-ove koje je sličnije funkciji
# read.csv(), možemo koristiti funkciju shapefile() iz paketa raster Napomena:
# ispod je dat samo primjer. Ukoliko želite učitati svoj shapefile koji ste
# snimili, potrebno je staviti pravilan put do file-a (HINT: file.choose())

BiH_boundary <- readOGR("C:\\Users\\MirzaC\\Desktop\\DSB_workshop\\R_directory", "BiH_boundary")


plot(BiH_boundary)

# Pomoću argumenata xlim i ylim možemo zumirati na određeni dio karte
plot(BiH_boundary, xlim = c(17,19), ylim = c(43, 44))

# Pomoću argumenta add = TRUE možemo dodavati dodatne elemente na kartu
plot(atra_loc, add = TRUE)

# Učitavanje podataka o nadmorskoj visini - DEM (Digital Elevation Model)
# SRTM (Shuttle Radar Topography Mission) je naziv projekta u kojem su korišteni podaci sa senzora Space shuttle-a Endeavour 
# iz kojeg je nastao DEM skoro cijele planete Zemlje u rezoluciji od 90m. Ovi podaci dolaze u pojedinim "pločama" (tile), 
#i ako želimo da imamo podatke za određeni dio Zemlje, obično je potrebno obraditi podatke. 

# Skidanje originalnih SRTM podataka. Raspored tile-ova možete pogledati na
# http://dwtkns.com/srtm/. Za BiH je potrebno da spojimo dva tile-a.

# Napomena: prilikom projere ove skripte serveri na kojima se nalaze SRTM i
# bioclim (vidi ispod) podaci nisu pravilno funkcionisali. Ukoliko budete imali
# problema sa skidanjem DEM sloja, preskočite slijedećih nekoliko linija i
# učitajte BiH_DEM.tif file na disku.

BiH_DEM_south <- getData(name = "SRTM", lon = 15, lat = 45) # 30mb
BiH_DEM_north <- getData(name = "SRTM", lon = 15, lat = 50)

# Koristeći funkciju mosaic() iz paketa raster možemo "zalijepiti" dva ili više rastera
# Napomena: može potrajati na sporijim računarima. Znamo da je komanda završila kada vidimo "prompt" u konzoli (simbol ">")
BiH_DEM <- mosaic(BiH_DEM_north, BiH_DEM_south, fun = mean)

# Uz pomoć funkcije crop(), možemo iskoristiti granicu BiH koju smo prethodno
# učitali da "izrežemo" mozaik dva rastera DEM-a. Funkcija crop() će izrezati
# raster u oblik pravougaonika (vidi funkciju extent()). Ako želimo da dobijemo
# raster koji neće izlaziti van granica, možemo korisiti funkciju mask().

BiH_DEM <- crop(BiH_DEM, BiH_boundary)

plot(BiH_DEM)
# Funkcija writeRaster() spavašava file na disk. Ukoliko ne naglasimo drugu
# lokaciju, funkcija u ovom obliku će spasiti file na lokaciju koju možemo
# dobiti uz pomoć funkcije getwd().

writeRaster(BiH_DEM, "BiH_DEM.tif", format = "GTiff")

# Za učitavanje DEM sloja možemo koristiti funkciju ispod (napomena: put do
# file-a ispod ne mora odgovarati putu gdje je spašen file na vašem PC-u). 
# Simbol "." (".\\Data\\") je skraćenica za working directory
BiH_DEM <- raster(".\\BiH_DEM.tif")

plot(BiH_DEM)
plot(BiH_boundary, add = TRUE)

# Da saznate kako spasiti plot u R-u koristeći komande, posjetite https://www.stat.berkeley.edu/~s133/saving.html

# Učitavanje podataka o klimi ----

# Funkcija getData() može biti korisna da skinemo worldclim podatke (Hijmans et al. 2005). 
# Ovi podaci su generisani interpolacijom podataka sa meteoroloških stanica, koristeći pri tome nadmorsku visinu, 
# geografsku širinu i dužinu kao nezavisne varijable. 
# Zbog ovoga potrebno je biti oprezan sa interpretacijom ovih podataka u područjima gdje je mali broj meteoroloških stanica. 
# Više informacija na: http://worldclim.org/.

# Definisanjem argumenta var kao "bio" skinut ćemo 19 bioclim varijabli (vidi
# link). Format u kojem R spašava ovakve podatke je raster stack. HINT: vidi
# argument "bylayer" u help file-u funkcije writeRaster() ako želiš spasiti ove
# podatke na disk.

# Napomena: ukoliko budete imali problema prilikom skidanja ovih file-ova, učitajte ove file-ove sa diska (vidi ispod funkciju stack())
BiH_bioclim <- getData("worldclim", var = "bio", res = 0.5, lon = 15, lat = 50) # 75 mb

# Ukoliko naš objekt ima više kolona ili "slotova" (vidi funkciju str()), možemo
# im pristupiti uz pomoć znaka "$" U našem slučaju imamo 19 pojedinačnih slojeva
# u objektu BiH_bioclim, a svakom pojedinačnom možemo pristupiti sa znakom "$", 
# ili u ovom slučaju sloju Bio1 (Annual mean temperature, vidi
# http://worldclim.org/bioclim).

# Učitavanje svih 19 spašenih file-ova istovremeno. Ukoliko dobivate grešku, provjerite put do file-ova.
BiH_bioclim <- stack(list.files(path= ".\\Bioclim\\",  pattern="*.tif$", full.names = TRUE))

plot(BiH_bioclim$Bio1_16)

BiH_bioclim <- crop(BiH_bioclim, BiH_boundary)

crs(BiH_bioclim) <- BiH_DEM@crs

# Dat ćemo imena bioclim slojevima (vidi http://worldclim.org/bioclim).
names(BiH_bioclim) <- c("BiH_annual_mean_temp",
                        "BiH_mean_diurnal_range",
                        "BiH_isothermality",
                        "BiH_temp_seasonality",
                        "BiH_max_temp_warmest_month",
                        "BiH_max_temp_coldest_month",
                        "BiH_temp_annual_range",
                        "BiH_mean_temp_wettest_quarter",
                        "BiH_mean_temp_driest_quarter",
                        "BiH_mean_temp_warmest_quarter",
                        "BiH_mean_temp_coldest_quarter",
                        "BiH_annual_precip",
                        "BiH_precip_wettest_month",
                        "BiH_precip_driest_month",
                        "BiH_precip_seasonality",
                        "BiH_precip_wettest_quarter",
                        "BiH_precip_driest_quarter",
                        "BiH_precip_warmest_quarter",
                        "BiH_precip_coldest_quarter")

# Sa argumentom bylayer podešenim na TRUE, svaki od 19 slojeva će biti spašen u
# poseban file. Grešku možete dobiti ukoliko nemate foldere Data i Bioclim u
# vašem radnom direktoriju. Možete ili promjeniti put gdje će file-ovi biti
# spašeni, ili napraviti iste direktorije.

writeRaster(BiH_bioclim, ".\\Bioclim\\", bylayer = TRUE,
            names(BiH_bioclim), format = "GTiff")



# Resampluj DEM (rezolucija 90m = veličina "pixela" je 90m) na rezoluciju bioclim podataka (rezolucija ~1km).
# Ukoliko povečavamo rezoluciju (upsampling), gubimo na kvaliteti podataka. Bolje rješenje je uraditi "downsampling", 
# ali u tom slučaju će trebati mnogo više vremena za računanje modela. U ovom slučaju nam je cilj da prođemo osnove bez puno
# zadržavanja, tako da ćemo prihvatiti lošiju soluciju koja će rezultirati "lošijim" modelima.

BiH_DEM_1km <- resample(BiH_DEM, BiH_bioclim, method = "bilinear")


crs(BiH_DEM_1km) <- BiH_DEM@crs
names(BiH_DEM_1km) <- "BiH_DEM"
# Iz DEM podataka moguće je dobiti razne podatke kao što su nagib terena
# (slope), ekspozicija (aspect) i sl. U ovom slučaju da dobijemo ove slojeve
# koristit ćemo funkciju terrain().

BiH_slope <- terrain(BiH_DEM_1km, opt = "slope", unit = "degrees", neighbors = 4)
BiH_aspect <- terrain(BiH_DEM_1km, opt = "aspect", unit = "degrees", neighbors = 4)

plot(BiH_slope)
plot(BiH_aspect)

# Budući da je aspect "kružna" varijabla (0-360°), gdje su vrijednosti od npr.
# 359° i 1° vrlo bliske u stvarnosti. Iz ovog razloga potrebno je preračunati
# ove varijable, i možemo koristiti funkcije sin() i/ili cos() da dobijemo mjeru
# ekspozicije prema istoku/sjeveru između -1 i 1
BiH_northness <- cos(BiH_aspect)
names(BiH_northness) <- "BiH_northness"
# Za slijedeći korak želimo vidjeti koje varijable su "kolinearne", tj. da li
# sadrže iste informacije što nije poželjno u statističkim modelima. U ovom
# slučaju izračunat ćemo VIF (Variance Inflation Factor) uz pomoć funkcije
# vifstep() iz paketa usdm. 
# http://www.statsmakemecry.com/smmctheblog/confusing-stats-terms-explained-multicollinearity.html

install.packages("usdm")
library(usdm)

# Budući da funkcija vifstep() može raditi na rasterstack formatu, kreirat ćemo novi stack od slojeva koje smo napravili iznad.
BiH_rasters <- stack(BiH_bioclim, BiH_DEM_1km, BiH_northness, BiH_slope)

# Ova funkcija će uklanjati visoko kolinearne varijable i rekalkulisati
# vrijednosti svaki korak sve dok sve varijable nemaju VIF vrijednost ispod 10. 
# VIF = 10 je odokativna vrijednost koja se uzima kao granica za kolinearnost.
# Restriktivniji modeli će imati VIF granice od 2 - 5.

# Napomena: može potrajati na sporijim računarima
vifstep(BiH_rasters, th = 10)

# Pregledajte tekst u konzoli. 13 varijabli su kolinearne, tako da ćemo
# nastaviti samo sa varijablama/slojevima koje imaju VIF<10
BiH_rasters_subset <- subset(BiH_rasters, c("BiH_mean_temp_warmest_quarter", "BiH_mean_temp_coldest_quarter", 
                                            "BiH_precip_wettest_month", "BiH_precip_driest_month",
                                            "BiH_precip_warmest_quarter", "BiH_precip_coldest_quarter",
                                            "BiH_DEM", "BiH_northness", "slope"))

# Za R je napravljeno više paketa za SDM (Species Distribution Modelling). 
# Jedan od paketa je biomod2 kojeg ćemo danas koristiti. Biomod2 ima 10 različitih statističih metoda za modelovanje. 
#

install.packages("biomod2", dependencies = TRUE)
library(biomod2)

# --------------------------------------------------------------------------------------------------
# Objasniti koncepte SDM-a i objasniti pojedine postavke u funkcijama ispod, i razloge zašto su tako postavljene je izvan dometa
# komentara u R skripti. Za više informacija o pojedinim argumentima možete pogledati help dokumentaciju za svaku funkciju, 
# a materijali vezani za SDM će biti također dostupni svima.

# Koristit ćemo default postavke za sve dostupne modele
BiomodOptions <- BIOMOD_ModelingOptions()

# Napomena: svaka od slijedeće tri komande će trajati duže vremena, posebno na sporijim računarima.

# Da proces bude brži, broj pojedinih parametara je namjerno smanjen (npr. PA.nb.rep i NbRunEval).
# U ovom slučaju ćemo napraviti samo 4 modela, dok bi za uobičajene analize bilo potrebno napraviti barem 100 modela.

# Pripremanje podataka za modelovanje
Atra_modeling_data <- BIOMOD_FormatingData(resp.var = rep(1, length(atra_loc$Subspecies)), expl.var = BiH_rasters_subset,
                                           resp.xy = atra_loc@coords, eval.resp.xy = NULL,
                                           PA.nb.rep = 2, 
                                           PA.nb.absences = 1000,
                                           PA.strategy = 'random', na.rm = TRUE, resp.name = "Sal_atra")



# Modelovanje uz pomoć boosted regression trees algoritma (GBM = Generalized Boosted Models)
Atra_model_out <- BIOMOD_Modeling(Atra_modeling_data, models = "GBM", models.options = BiomodOptions, NbRunEval=2, 
                                  DataSplit=80, VarImport=10, models.eval.meth = c("TSS","ROC"),
                                  SaveObj = FALSE, rescal.all.models = TRUE, do.full.models = FALSE,
                                  modeling.id = "Sal_atra_GBM")


# Projektovanje modela u geografski prostor
Atra_model_projection <- BIOMOD_Projection(modeling.output = Atra_model_out, new.env = BiH_rasters_subset, proj.name = "current", 
                                           selected.models = "all", binary.meth = "TSS", compress = "xz",
                                           build.clamping.mask = FALSE, output.format = ".img")

# ----------------------------------------------------------------------------------------------------------------------


#

# Sada ćemo pogledati koje su varijable (ekološki faktori) najznačajniji za
# model. Rezultati se mogu interpretirati (ali oprezno) i kao najznačajniji
# ekološki faktori za distribuciju. Vrijednosti su raspona od 0 - 1. Što je 
# vrijednost bliža 1, varijabla je "bitnija" za model. Za objašnjenje kako se 
# ova vrijednost računa, pogledajte help dokument za funkciju variables_importance()

var_importance <- as.data.frame(get_variables_importance(Atra_model_out))

var_importance

boxplot(t(var_importance), main = "Značaj varijabli")

# Izračunavanje metrika kvalitete modela (AUC i TSS)
Atra_model_accuracy <- as.data.frame(get_evaluations(Atra_model_out))

Atra_model_accuracy <- t(Atra_model_accuracy[, seq(1,16,4)])

# TSS i ROC (Receiving Operator Curve) su metrike kvalitete modela. Što je vrijednost bliža 1, možemo reći (ali jako oprezno) da je model bolji.
# U našem slučaju ove vrijednosti treba uzeti sa velikom pažnjom, jer smo
# resamplovali naše varijable na 1km, i samo imamo 4 modela da izvučemo
# zaključke, što je prenizak broj za bilo kakav statistički validan zaključak.

Atra_model_accuracy


# ----------------------------------------------------------------------------------------------------------------------
# Prikaz modeliranih projekcija
plot(Atra_model_projection)

# Uzet ćemo srednje vrijednosti naša 4 modela
Atra_model_mean <- mean(Atra_model_projection@proj@val)

plot(Atra_model_mean)
# Naš finalni sloj ima raspon vrijednosti od 0 - 1000. Iako modeli poput ovih
# imaju vrijednosti od 0 - 1, iz razloga što "float" brojevi (sa decimalnom
# tačkom) zauzimaju više memorije, u ovom slučaju vrijednosti su od 0 - 1000. Da
# vratimo vrijednosti od 0 - 1 možemo jednostavno podijeliti naš sloj sa 1000.

Atra_model_mean <- Atra_model_mean/1000
plot(Atra_model_mean)


# Za ljepši prikaz na karti, pretvorit ćemo sve vrijednosti niže od 0.1 u NA (Not Available) 
Atra_model_mean[Atra_model_mean < 0.1] <- NA

# Prikaz granice i originalnih tačaka
plot(BiH_boundary, add = TRUE)
plot(atra_loc, add = TRUE)

# Spašavanje mape kao kml file kojeg je moguće otvoriti u Google Earthu
KML(Atra_model_mean, file = "Atra_niche.kml", overwrite = TRUE)

# Instaliranje paketa leaflet za pravljenje web-mapa
install.packages("leaflet")
library(leaflet)

# Kreiranje palete boja za prikaz podataka 
pal <- colorNumeric(c("#FFFFCC", "#e9a3c9", "#af8dc3"), values(Atra_model_mean),
                    na.color = "transparent")

# Napravit ćemo leaflet mapu koja je pogodna za korištenje na web-u (mapa je interaktivna)
atra_leaflet_map <- leaflet() %>% addTiles() %>%
  addRasterImage(Atra_model_mean, colors = pal, opacity = 0.7) %>%
  addLegend(pal = pal, values = values(Atra_model_mean),
            title = "Vrijednost")

atra_leaflet_map

install.packages("htmlwidgets")
library(htmlwidgets)
saveWidget(atra_leaflet_map, file="Atra_leaflet_map.html")

# Otvorite spašeni html file

# Za sva pitanja: mirzaceng@gmail.com
