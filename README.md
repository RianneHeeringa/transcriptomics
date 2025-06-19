# transcriptomics
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
De Gene Ontology-analyse werd uitgevoerd met het R-pakket “goseq”. Omdat langere genen vaker worden gevonden in RNA-seq-data, werd eerst gecorrigeerd voor genlengte om vertekening te voorkomen. Met behulp van geninformatie uit “org.Hs.eg.db” werd vervolgens gekeken welke biologische processen vaker voorkwamen in de lijst met significante genen. GO-termen met een p-waarde kleiner dan 0,05 werden als statistisch significant beschouwd. De top 10 meest opvallende biologische processen werden uiteindelijk gevisualiseerd met ggplot2, waarbij de termen gesorteerd werden op relevantie. [flowschema]()




