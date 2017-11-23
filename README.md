# ISTRUZIONI PER COMPILAZIONE ED ESECUZIONE PROGRAMMA

Il programma gm.c deve essere eseguito su sistemi Unix (su windows manca una metodo delle librerie standard di C).


Per compilare il programma, è necessaria la libreria OpenMP.

`gcc gm.c -fopenmp linalg.c matrix.c scan.c -l gomp -lgmp -o gm -O3`

Per eseguire il programma

`./gm --options`


**options**:

--**input input_file**: viene utilizzato il file **input_file** come input, in assenza di tale opzione viene usato stdin come input

--**output file_name**: viene stampato il risultato sul file **file_name**, in assenza di tale opzione viene stampato il risultato su stdout

--**execute n**: per scegliere quale tecnica risolutiva utilizzare, **n** scelta nel seguente elenco (default 0)

0 -> standard (si moltiplicano tutti i polinomi)

1 -> confronto (si moltiplicano solo i polinomi il cui monomio di grado maggiore è cambiato)

2 -> eliminazione (non si moltiplicano i polinomi iniziali rimasti immutati)

--**reduced-expand**: moltiplica solamente i monomi formabili con grado massimo calcolato in base alle righe mancanti per la soluzione (vecchia moltiplicazione_ridotta)

--**set-expand n**: moltiplica solamente i monomi formabili con grado massimo pari a **n** (non funziona bene col metodo confronto)

--**manual-expand**: ad ogni fase di espansione viene chiesto di inserire il grado massimo dei monomi che verrano moltiplicati (non funziona bene col metodo confronto)

--**partial-gauss**: ad ogni passo chiede se effettuare l'eliminazione gaussiana su una porzione della matrice

--**verbose**: stampa la matrice completa ad ogni passo

--**loops n**: cambia la condizione di terminazione (quella standard è data dalla matrice rimasta immutato dopo un passo). n indica il numero massimo di passi consecutivi nella computazione in cui i gradi mancanti restano immutati.

--**verify x y**: per eseguire un confronto tra i risultati di due tecniche, x e y scelte dall'elenco precedente (generalmente usata per il testing). 

## ESEMPI DI ESECUZIONE

`./gm --execute 0 --input input.txt --output output.txt`

esegue il metodo standard sull'input contenuto nel file input.txt

`./gm --execute 1 --reduced-expand --input input.txt --output output.txt`

esegue il metodo confronto, con fase di moltiplicazione ridotta, sull'input contenuto nel file input.txt

`./gm --execute 2 --set-expand 3 --input input.txt --output output.txt`

esegue il metodo eliminazione, con grado massimo dei monomi in fase di moltiplicazione pari a 3, sull'input contenuto nel file input.txt

`./gm --execute 2 --manual-expand --input input.txt --output output.txt`

esegue il metodo eliminazione, chiedendo ad ogni passo il grado massimo dei monomi in fase di moltiplicazione, sull'input contenuto nel file input.txt


## FORMATO DEI FILES DI INPUT

Per eseguire correttamente il programma occorre fornire i dati di input nel corretto ordine e formato.
- modulo dei coefficienti
- grado massimo raggiungibile
- numero polinomi di partenza
- tipo di ordinamento
- elenco delle variabili utilizzate nei polinomi (ogni variabile è un singolo carattere, le variabili sono su una sola riga in ordine alfabetico non separate da spazi)
- elenco di polinomi (uno per riga, con il formato utilizzato da Magma, esempio nel file)

esempio di file di input https://github.com/FRSYH/Tesi/blob/master/input.txt

In caso di malfunzionamento del programma probabilmente l'input ha un formato scorretto.

## ORDINAMENTO POLINOMI

Ordinamenti implementati e relativi codici di input:
- grevlex_comparison -> 0


















