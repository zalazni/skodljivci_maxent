# skodljivci_maxent
Modeliranje porazdelitev vrst (koda in izračuni) s pomočjo metode maksimalne entropije, primer Ameriške koruzne sovke (Spodoptera frugiperda) in njenih naravnih sovražnikov Trichogramma pretiosum, Cotesia marginiventris, Telenomus remus in Eiphosoma laphygmae. 

Podatki, potrebni za izračun: 
- količina padavin,
- minimalna temperatura,
- maksimalna temperatura in
- evapotranspiracija

V tem primeru so podatki pridobljeni s portala Odprti podatki Slovenije, uporabimo podatke za referenčno obdobje 1981-2010 in podnebne projekcije za 2011-2040, 2041-2070 ter 2071-2100 (RCP4.5 in RCP8.5):
Podnebne spremembe: https://podatki.gov.si/data/search?q=podnebne+spremembe&page=1 

Podatki o skupni površini polj s koruzo:
AKTRP. 2022. Data on Crops Grown on Slovenian Arable Fields in 2022. https://rkg.gov.si/razno/portal_analysis/KMRS_2022.rar

Podatki o razširjenosti vrst pridobljeni s strani GBIF:
Global Biodiversity Information Facility (https://www.gbif.org/)

Model najprej konstruiramo na podatkih WorldClim, nato izračunamo ustreznost podnebja Slovenije (referenčno obdobje 1981-2010) in v tretjem delu skripte izračunamo ustreznost za pojav v podnebju podnebnih projekcij.

DOI: https://doi.org/10.5281/zenodo.20342534
