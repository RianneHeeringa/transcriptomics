# Genexpressie bij Reumatoïde artritis
![image](https://github.com/user-attachments/assets/6da2d9cb-315a-4ea1-a778-6e86430bd666)

Rianne Heeringa (5144086)/
LBM2A/ 
J2P4/
20-06-2025

## 📁 Inhoud/structuur

- `data/raw/` – fictionele datasets voor de analyse van spreuk effectiviteit, gevaar en welke spreuken het beste samengaan met verschillende types staf.  
- `data/processed` - verwerkte datasets gegenereerd met scripts 
- `scripts/` – scripts om prachtige onzin te genereren
- `resultaten/` - grafieken en tabellen
- `bronnen/` - gebruikte bronnen 
- `README.md` - het document om de tekst hier te genereren
- `assets/` - overige documenten voor de opmaak van deze pagina
- `data_stewardship/` - Voor de competentie beheren ga je aantonen dat je projectgegevens kunt beheren met behulp van GitHub. In deze folder kan je hulpvragen terugvinden om je op gang te helpen met de uitleg van data stewardship. 


## Inleiding
Reumatoïde artritis (RA) is een chronische, systematische auto-immuunziekte. De ziekte wordt gekenmerkt door ontstekingen in het gewrichtsslijmvlies, wat uiteindelijk kan leiden tot gewrichtsschade en verlies van functie.(Firestein & McInnes, 2017) Hoewel de precieze oorzaak nog niet bekend is, wordt wel aangenomen dat een combinatie van genetische aanleg, omgevingsfactoren en een ontspoord immuunsysteem een belangrijke rol speelt. De diagnose wordt doorgaans gesteld op basis van klinische symptomen, lichamelijk onderzoek en het aantonen van specifieke autoantistoffen, zoals reumafactor en anti-CCP.RA is een progressieve aandoening die, zonder behandeling, steeds ernstiger kan worden. Er is momenteel geen genezing mogelijk, maar met medicatie en vroegtijdige behandeling kunnen de symptomen onder controle worden gehouden en verdere gewrichtsschade worden vertraagd(Díaz-González & Hernández-Hernández, 2023). In dit onderzoek wordt uitgezocht welke genen hoger of lager tot expressie komen in personen met reuma. daarnaast wordt onderzocht welke metabole routes anders functioneren bij reuma. De deelvragen van dit onderzoek zijn dan ook; Welke genen zijn verhoogd en welke verlaagd bij personen met reuma? Welke metabole routes functioneren anders bij reuma?  [bronnenlijst](1.inleiding/bronnenlijst.pdf).

## Methode
Sequencing methode
In dit onderzoek werd gebruikgemaakt van RNA-sequencingdata afkomstig van synoviumbiopten van 4 gezonde personen en 4 personen met RA. Personen met RA waren positief getest op ACPA, personen zonder RA negatief. ACPA, ook wel anti-CCP genoemd, meet autoantistoffen tegen het CCP (cyclische gecitrullineerde peptiden) eiwit. De monsters waren geprepareerd met een paired-end sequencingmethode, waarbij zowel de voorwaartse als achterwaartse streng van het RNA-fragment werd uitgelezen. Deze methode bood een hogere nauwkeurigheid en betere mapping, met name in complexe gen regio’s.
Mappen
De verkregen humane data werden gemapt met behulp van het R-pakket “Rsubread”. Eerst werd een index gegenereerd met “buildindex”, waarna afzonderlijke reads uit de 8 monsters werden uitgelijnd met “align”. De gemapte reads werden opgeslagen in BAM-bestanden.
Differentiële genexpressie-analyse
Met de functie “featureCounts” werd een count-matrix gegenereerd, waarin per monster het aantal reads per gen werd geteld. Deze matrix werd vervolgens ingelezen in “DESeq2”, waarmee normalisatie, dispersieschatting en statistische toetsing werden uitgevoerd. Hierbij werd getest welke genen significant verschillend tot expressie kwamen tussen de personen met RA en gezonde personen (met p-adjust < 0,05 en |log2FoldChange| > 1). De resultaten werden gevisualiseerd in een volcano plot met behulp van het pakket “EnhancedVolcano”.
KEGG-pathwayanalyse
Na de differentiële genexpressie-analyse werd een KEGG-pathwayanalyse uitgevoerd met behulp van de R-pakketten “goseq”, “clusterProfiler” en “pathview”. Hierbij werd gecorrigeerd voor genlengte en gebruikgemaakt van humane genannotaties uit “org.Hs.eg.db”. De KEGG-pathways werden bepaald op basis van genen met een aangepaste p-waarde kleiner dan 0,05. Belangrijke biologische routes werden vervolgens gevisualiseerd met “pathview”, waarin te zien was welke genen meer of minder actief waren binnen die routes.
Gene Ontology-analyse
De Gene Ontology-analyse werd uitgevoerd met het R-pakket “goseq”. Omdat langere genen vaker worden gevonden in RNA-seq-data, werd eerst gecorrigeerd voor genlengte om vertekening te voorkomen. Met behulp van geninformatie uit “org.Hs.eg.db” werd vervolgens gekeken welke biologische processen vaker voorkwamen in de lijst met significante genen. GO-termen met een p-waarde kleiner dan 0,05 werden als statistisch significant beschouwd. De top 10 meest opvallende biologische processen werden uiteindelijk gevisualiseerd met ggplot2, waarbij de termen gesorteerd werden op relevantie. [script.R](2.methode/transcriptomics.script.R)

## Resultaten
In dit onderzoek is RNA-sequencingdata van gezonde mensen en mensen met reumatoïde artritis (RA) geanalyseerd met behulp van R. Dit is gedaan met de pakketten “Rsubread”, “DESeq2” en “goseq”. Na kwaliteitscontrole zijn de Reads uitgelijnd tegen het humane referentiegenoom (GRCh38.p14).  
Vervolgens is er differentiële genexpressieanalyse uitgevoerd met “DESeq2”, waarbij de monsters zijn gegroepeerd in normaal of reuma. Uit deze analyse kwamen meerdere significante gereguleerde genen naar voren. In totaal werden 166 genen significant differentieel tot expressie gebracht, waarvan 103 genen up gereguleerd waren bij RA. Belangrijke ontsteking gerelateerde genen zoals IL1B, IL6 en CXCL8 vertoonden verhoogde expressie in de RA-groep, zoals te zien in figuur 2. 
Om te bepalen welke biologische processen betrokken zijn, is een GO-analyse uitgevoerd met het “goseq” pakket. De top 10 meest GO- termen zijn weergeven in figuur 1. Hieruit blijkt dat genen betrokken bij onder andere “immune system process” en “immune respons” significant zijn bij RA. Dit komt overeen met de verwachting dat immuun activiteit een belangrijke rol speelt bij RA. 
Ook is er nog een KEGG-pathway geanalyseerd. Figuur 4 toont aan dat IB significant verhoogd is wat een rol speelt bij ontstekingscascades wat dus bij het beeld van RA past. Ook is CCL2 verhoogd, welke belangrijk is voor het aantrekken van immuun cellen naar gewrichten. In figuur 5 is de KEGG pathway van de cell cyclus te zien, bij gezonde personen verwacht je een normale en evenwichtige cell cyclus. In figuur 5 is een sterke up regulatie van de genen PP2A en APC/C te zien, dit kan wijzen op een verhoogde en minder gecontroleerde celdeling in immuun cellen.  In figuur 6 is de KEGG pathway complement and coagulation cascades te zien. Hier is te zien dat onder andere C4 verhoogd is wat een sterke ontstekingsmediator is, ook is er een verhoogd immuun activiteit te zien aan de verhoogde genen CR1, CR2 en CR4 en een verhoogde gevoeligheid voor inflammatoire mediatoren te zien door de verhoging van B1/B2.

## conclusie
In dit onderzoek is een duidelijk significant verschil in genexpressie tussen gezonde en mensen met RA gevonden. Er zijn 166 genen gevonden die verschillend tot expressie kwamen waarbij de meerderheid up gereguleerd is bij RA. Belangrijke ontstekingsgenen zoals IL1B, IL6 en CXCL8 hadden een verhoogde expressie, wat betekent dat ontstekingen een rol speelt bij RA. De GO-analyse toont aan dat biologische processen gerelateerd aan het immuunsysteem en immuunrespons sterk betrokken zijn. Daarnaast wijst de KEGG-analyses op een verhoogde immuun activiteit en een mogelijk ontregelde celdeling in immuun cellen, wat kan bijdragen aan de ziekteverloop. De complement en coagulase cascade lijken ook actiever te zijn wat de ontstekingsreactie versterkt. Voor vervolgonderzoek wordt aanbevolen om te kijken naar de verandering in genexpressie in samenhang met de ernst van de ziekte van de patiënten dit kan helpen bij betere behandeling. Verder kan het ook interessant zijn om te kijken naar de genen die betrokken zijn bij celdeling omdat deze mogelijk een rol spelen bij het ontstaan van RA.


