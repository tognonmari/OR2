Repository for OR2 project
TODO
#ifndef e #define dello stesso file serve per gestire problemi di confitto
#ifndef significa se questa variabile non è definita viene definita tutto quello fino all' endifndef altrimenti essendo gia definita non la ridifenisce
meglio fare la struttura per i punti.
FARE MATRICE CON I COSTI.
fare tsp_free per liberare l'istanza.
fare tsp_init per inizializzare istanza
fare sempre un assert quando inizializzo i vettori.
chiamare 100 volte rand a vuoto perche l'inizio della sequenza di generazione è polarizzata.

HW
EURISTICO vuol dire posso trovaRe una buona soluzione ma non posso garantire che sia buona.
GREEDY (NEAREST NEIGHBORH)
definiamo una citta iniziale scegliendo la citta 0, man mano creiamo il cammino andando sul punto piu vicino.
alla iterazione piu generica abbiamo un cammino (v0,v1,..,vK) dove tutti questi nodi sono visitati e vK è il last node visitato.
INDUZIONE-> il cammino si estende scegliendo vK+1 piu vicino a vK (quindi minimizza il costo cK,k+1)
alla fine devo chiudere il cammino facendo un arco con il nodo iniziale 0.
IMPLEMENTARE STRUTTURA GREEDY.
serve un vettore seq_nodes in cui tra gli indici 0 e len-1 ci sono i nodi visitati e nelle posizioni successive ci sono i nodi unvisited (last è il nodo SEQ[len - 1])_
Ad ogni iterazione in O(n) trovo quello a distanza minima, complessivamente la complessità è O(n^2).
inizio nel vettore 0 c'è un nodo qualsiasi quindi posso fare oltre al metodo che funziona partendo quel nodo un ciclo for che prova tutte le varie soluzioni in cui per ogni soluzione cambia solo
il primo nodo in cui inizi.
ABBIAMO SEMPRE BEST_SOL (vettore con path ottimo) e BEST_VAL e BEST_TIME(tempo per trovarla)_
CE UN UNICA FUNZIONE DELEGATA AD AGGIORNARE BEST_SOL e BEST_VAL. (tsp_update_bestsol).
tsp_check_sol per checkare se la soluzione è un hamiltonian cycle (usare un flag di debug).
nel check della soluzione uso un array di contatori che salva il conteggio di visita dei nodi TUTTI DEVONO ESSERE VISITATI UNA VOLTA. verifico anche che il valore/costo
della soluzione sia uguale a best_val.(dare errore)
ATTENZIONE BISOGNA FARE UN METODO PER CONFRONTARE I DUOUBLE!! USARE UNA PRECISIONE COME COSTANTE.
ULTERIORE OTTIMIZZAZIONE se nel plot ci sono degli intrecci nei path. se ci sono interecci allora probabilmente (dimostrabile dice il prof) la soluzione non è cosi tanto buona.
si puo fare un otttimizzazione con la funzione 2_opt per eliminare gli intrecci.(non obbligatorio, cercare in rete come fare la funzione 2_opt il prof non aveva tempo).