Repository for OR2 project
TODO
#ifndef e #define dello stesso file serve per gestire problemi di confitto
#ifndef significa se questa variabile non � definita viene definita tutto quello fino all' endifndef altrimenti essendo gia definita non la ridifenisce
meglio fare la struttura per i punti.
FARE MATRICE CON I COSTI.
fare tsp_free per liberare l'istanza.
fare tsp_init per inizializzare istanza
fare sempre un assert quando inizializzo i vettori.
chiamare 100 volte rand a vuoto perche l'inizio della sequenza di generazione � polarizzata.

HW
EURISTICO vuol dire posso trovaRe una buona soluzione ma non posso garantire che sia buona.
GREEDY (NEAREST NEIGHBORH)
definiamo una citta iniziale scegliendo la citta 0, man mano creiamo il cammino andando sul punto piu vicino.
alla iterazione piu generica abbiamo un cammino (v0,v1,..,vK) dove tutti questi nodi sono visitati e vK � il last node visitato.
INDUZIONE-> il cammino si estende scegliendo vK+1 piu vicino a vK (quindi minimizza il costo cK,k+1)
alla fine devo chiudere il cammino facendo un arco con il nodo iniziale 0.
IMPLEMENTARE STRUTTURA GREEDY.
serve un vettore seq_nodes in cui tra gli indici 0 e len-1 ci sono i nodi visitati e nelle posizioni successive ci sono i nodi unvisited (last � il nodo SEQ[len - 1])_
Ad ogni iterazione in O(n) trovo quello a distanza minima, complessivamente la complessit� � O(n^2).
inizio nel vettore 0 c'� un nodo qualsiasi quindi posso fare oltre al metodo che funziona partendo quel nodo un ciclo for che prova tutte le varie soluzioni in cui per ogni soluzione cambia solo
il primo nodo in cui inizi.
ABBIAMO SEMPRE BEST_SOL (vettore con path ottimo) e BEST_VAL e BEST_TIME(tempo per trovarla)_
CE UN UNICA FUNZIONE DELEGATA AD AGGIORNARE BEST_SOL e BEST_VAL. (tsp_update_bestsol).
tsp_check_sol per checkare se la soluzione � un hamiltonian cycle (usare un flag di debug).
nel check della soluzione uso un array di contatori che salva il conteggio di visita dei nodi TUTTI DEVONO ESSERE VISITATI UNA VOLTA. verifico anche che il valore/costo
della soluzione sia uguale a best_val.(dare errore)
ATTENZIONE BISOGNA FARE UN METODO PER CONFRONTARE I DUOUBLE!! USARE UNA PRECISIONE COME COSTANTE.
ULTERIORE OTTIMIZZAZIONE se nel plot ci sono degli intrecci nei path. se ci sono interecci allora probabilmente (dimostrabile dice il prof) la soluzione non � cosi tanto buona.
si puo fare un otttimizzazione con la funzione 2_opt per eliminare gli intrecci.(non obbligatorio, cercare in rete come fare la funzione 2_opt il prof non aveva tempo).

LEZ 12/3
verbose da leggere da linea di comando, verbose = 0 non voglio vedere niente, piu aumenta piu roba vedo.
con verbosit� alta devo controllare che ogni funzione worki.
quando stampo le soluzioni controllare che la soluzione ci sia.
il prof dice di dichiarare le variaibili inline quando possibile. 
in opt2 si dovrebbe controllare tutti i possibili improvement e poi prendere il migliore. altrimenti potresti fare un improvement infinitesimo (e peggio potresti farne tantissimi), mentre devi fare gli improvement migliori!.
trovi coppia i,j migliore per fare lo scambio.
Alla fine di ogni greedy il costo della soluzione.
basta mettere .gitkeep per pushare cartelle senza il conten
LEZ 13/3
Il greedy parte da 0 e lui prosegue.
Ma ci sono vari euristici diversi, ma il greedy è uno dei più efficaci. Un altro euristico è extra milage che considera i pti nel pinao e considera un ciclo ragionevole che non copre tutti i punti.
dobbiamo estendere il ciclo fino a raggiungere un ciclo che copre tutti i vertici.
Iterativamente si segue una logica greedy.
Chiamiamo conn C i lati del ciclo per ogni vertice che non appartiene a V(C) considero tutti i lati i,j nel ciclo e conisidero tutti gli extramilage EM(h, (i,j)) = c_ih + c_hj - c_ij.
Quindi per ogni lato nel ciclo e per ogni vertice non coperto calolo tutti gli extramilage in O(n^2) e aggiungo l'h con l'EM minore.
ad ogni iteraz faccio la scelta piu logica, pero obvvimanete questo non porta alal soluzione ottima. UN punto delicato è l'inizialiizzazione. scegliete due vertici casuali e create un ciclo tra questi due vertici, la scelta di questi vertici è completamente arbitraria, di oslito viene suggerito di prenderli il piu lontano possibile. un altra possibilita è inizializzare il ciclo con un traingolo tra 3 vertci abbastanza distanti.
Un altra possibilita è inizializzare con tutti i vertici esterni. QUEST' ANNO NON VA IMPLEMENTATO EM!!!.
Quindi ce lo sta spiegando tanto per!
Questi algoritmi sono deterministici. è utile a volte randomizzare l'algoritmo in modo da trovare soluzioni diverse.
Come posso randomizzare un algoritmo deterministico?
TECNICA GRASP.
Fa parte del baggaglio di chi fa euristici.
Supponiamo di avere un path che va da 0 a last, proseguiamo scegliendo il nodo piu vicino a last. Questa cosa veiene fatta il 99%, ma con l'1% decidiamo di andare al secondo minimo oppure sceglie uno tra la posizione 2 e 6 dei minimim (a scelta nostra di progettaszione). Attenzione qeul 99% può essere modificato ma non tipo a 50% perche fa schifo, bisogna tenere una probabilita di deviare dalla scelta greedy con una probabilita bassa. L''idead di grasp puo essserer usata anche per altri algoritmi euristici. CURIOSITA Dentro chat gpt i modelli linguistici di intelligenza artificial usano una tecnica simile per generare il testo per generare non sempre lo stesso testo. QUando chiediamo ad un AI di darci un testo piu fantasioso semplicemente diminuisce il valore di p con cui prende la prima parola più spontanea come succesiiva.
Ha senso usare grasp_greedy con l'idea di fare piu esecuzioni.
Il prof dice: o provi tutti i vertici iniziali con greedy normali, oppure provi con grasp facendo all'inizio greedy normale con un nodo iniziale arbitrario e le successivele fai con grasp. ma cmq mi sembra molto a scelta nostra.
IL PROF DICE CHE NON SERVE IMPLEMENTARE GRASP.


HW DI OGGI
Tecnica che va sotto il nome di METAHEURISTICS.
Noi ci troviamo nella sitauzione in cui gia conosciamo 2_opt.
abbiamo due soluzioni TOur T1 e T2 io posso calcolare una distanza tra T1 e T2, la distanza tra due soluzioni si calcola vedendo quali lati ci sono diversi tra t1 e t2. Se t1 e t2 fossero implementati come sequenza di flag che indicano se i lati ci sono nel tour questa distanza sarebbe la distanza ddi hamming tra i due tour. DISTANZA TRA TOUR NUMERO DI LATI DIVERSI.
Data una soluzione T0 considero tutte le soluzioni cambiando solo due lati (al piu due lati) questo definisce un intorno in cui dentro ce anche la soluzione T0. Questo intorno viene chiamato intorno di 2-opt. Quindi noi abbiamo i nostri tour che avranno un loro costo specifico. noi partiamo da un tour e facciamo una neighb. search partendo da T0 e cerchiamo nell'intorno cerchiamo quella di costo minimo T1. A questo punto ci spostiamo da T1 con lo stesso principio. Ad un certo punto ci fermiamo perche siamo in un ottimo locale (che si fa??).
Un idea è il multistart cioe partire da una soluzione diversa da T0, molto probabilmente cadro in un ottimo locale diverso e potrei cadere nel ottimo globale.
Negli anni 70 FRED GROVER si e inentato una tecnica che si chiama TABU SEARCH per sfuggire dagli ottimi locali. Lui dice ok, voi siete partiti da una soluzione T0 e siete riusciti a miglioraRlo localmente. Arriviamo al local optimum di T0 ora cosa fare?
Ho capito che sono in un ottimo locale, avendo capito cio provo a fare uno scambio sconveniente (peggiorativo) e quindi la funzione obiettivo (costo) aumenta, ma il probelma che alla prossima iterazione tornerei indietro, ma FRED GROVER ha un idea, dopo avere fatot lo scambio peggiorativo mi creo una llista di scambi/mosse tabu, una lista di mosse vietate che non posso fare. quindi qunado arrivo nell ottimo faccio lo scambio peggiorativo e faccio un altra mossa peggiorativa perche ho dichiaraTo quella di tornare all'ottimo come tabu e iterativamente dichiaro anche la mossa di tornare nel punto precedente un tabu. L'idea è quella che peggiorando sempre prima o poi risalendo la curva della funzione costo trovo un punto in cui cambia la convessita della curva e ricomincio a scendere fino a che non trovo un nuovo minimo locale. se tale minimmmo locale ha un costo minore di quello precedente locale aggiorno anche l'incumbent che era fermo al vecchio minimo locale, con questa idea mi posso spostare lungo la funzione obiettivo definita su tutte le possibili soluzioni ordinate per intorni di 2-opt (due soluzioni sono adiacenti se èpossono essere ottenute da un 2-opt ottimo NO NON CHIARO, DAL SUO DISEGNO SEMBRA COSI MA IL SUO DISEGNO E' DUE DIMENSIONALE MENTRE REALISTICAMENTE SAREBBE MULTIDIMENSIONALE). Grover definisce una cardinalita massima di questa lista tabu ceh si chaima tenure che pou essre ad esempio di 1000 mosse. quando arrivo a  1000 cancello il tabu piu lontano. La lista è una fifo. Bisogna tarare il nostro algoritmo con il parametro giusto di tenure. ha senso avere una tenure che oscilla -> reactive tabular search. scegliere uan regola in cui si intenisifica la ricerca (tenure grande) e fase in cui si diversifica la ricerca (tenure piccolo).
Bisonga definire cosa è la mossa tabu, nell'it precedente ho fatto questo scambio, ho cambiato il lato i-j e il lato h-k con l'intreccio, devo evitare al colpo dopo di tornare indietro cioe tornare al lato i-h e lato j-k, oppure posso dire evito di tornare indietro dicendo all'algoritmo prendendo uno dei 2 vertici coinvolti nello scambio i o h e congelare tale vertice facendo in modo da proseguire solo con scambi che non coinvolgono questi vertci. Il prof decide uno dei due vertici lo dichiaro tabo e congelo i due lati che lo toccano, cioe non cambiero quei due lati che toccano il vertice finche il vertice è nella tabu list. La tenure con questa scekta non puo essere grande (circa n/10) il bello di questa scelta è che vi basta tenere un contatore di iterazione, poi avete un array in cui in ogni vertice i avete l'iterazione in cui lui è diventato tabo. vi basta confronare iter con tabu[h] se iter < tabu[h] + tenure, allora il vertice è ancora bloccato.
Per l'incumbent checkare sempre anche se scambi perggiorativi, per comodita.

18/03
come suggerimento mettere dentro la funzione di update il controllo che updata solo se è migliorativo il costo, FARE CONTROLLO CHE SIA EFFETTIBAMENTE UN TOUR.
Cercare di imparare il debug. tenure (1*sin(counter *1/50)+1)