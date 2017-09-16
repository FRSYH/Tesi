# ISTRUZIONI PER COMPILAZIONE ED ESECUZIONE PROGRAMMA

Il programma gm.c deve essere eseguito su sistemi Unix (su windows manca una metodo delle librerie standard di C).


Per compilare il programma, necessaria la libreria OpenMP.

`gcc gm.c -fopenmp linalg.c matrix.c scan.c -l gomp -lgmp -o gm`

Per eseguire il programma

`./gm --options < input_file > output_file`


**options**:

--**verbose**: stampa delle matrici nei passi intermedi

--**test**: stampa di informazioni aggiuntive per testing

--**execute n**: per scegliere quale tecnica risolutiva utilizzare, n scelta nel seguente elenco

0 -> standard (normale moltiplicazione)

1 -> confronto (si moltiplicano le sole righe che sono state ridotte nel Gauss precedente)

2 -> eliminazione 

3 -> moltiplicazione_ridotta (si moltiplica la matrice fino al grado massimo mancante)

4 -> confronto_ridotto (tecnica 1+2)

5 -> eliminazione_ridotta (tecnica 1+3)


--**verify x y**: per eseguire un confronto tra i risultati di due tecniche, x e y scelte dall'elenco precedente. 


## ESEMPIO DI ESECUZIONE

`./gm --execute 0 < input.txt > output.txt`

esegue il metodo standard sull'input contenuto nel file input.txt



## FORMATO DEI FILES DI INPUT

Per eseguire correttamente il programma occorre fornire i dati di input nel seguente ordine e formato.
- modulo dei coefficienti
- grado massimo raggiungibile
- numero polinomi di partenza
- tipo di ordinamento
- elenco delle variabili utilizzate nei polinomi (su una sola riga in ordine alfabetico)
- elenco di polinomi (uno per riga, con il formato utilizzato da Magma, esempio nel file)

esempio di file di input https://github.com/FRSYH/Tesi/blob/master/input.txt


## ORDINAMENTO POLINOMI

Ordinamenti implementati e relativi codici di input:
- grevlex_comparison -> 0


















