# Tesi
Repo per sviluppo tesi crittografia

Librerie per Algebra Lineare:
- GSL
- BLAS
- LAPACK
- librsb (per matrici sparse, architetture parallele)
- ATLAS (basto su BLAS)

Link :
- https://en.wikipedia.org/wiki/Comparison_of_linear_algebra_libraries
- http://stackoverflow.com/questions/6977677/blas-lapack-routine-for-doing-gaussian-elimination  LU = Gauss
- https://en.wikipedia.org/wiki/LU_decomposition
- http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=104 Nessuna procedura per eseguire riduzione di Gauss diretta su LAPACK
- http://www.netlib.no/netlib/lapack/double/dgetrf.f LAPACK esegue LU decomposition su matrici di qualsiasi dimensione http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html


Papers:
- http://www.netlib.org/lapack/lawnspdf/lawn127.pdf  LU
- http://www.netlib.org/utk/people/JackDongarra/PAPERS/beautiful-code.pdf  LU on LAPACK
- http://www.math.iit.edu/~fass/477577_Chapter_7.pdf  perchè LU è meglio di Gauss
- https://www.math.ucdavis.edu/~linear/old/notes11.pdf LU su matrici non quadrate


Problema della imprecisione delle operazioni Floating point 
- https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
- http://www2.lawrence.edu/fast/GREGGJ/Math420/Section_6_2.pdf --> Bisogna usare partial pivoting.


GSL,LAPACK e BLAS non implementano aritmetica modulare

Librerie per aritmetica modulare:
- FLINT                                -> http://www.flintlib.org/
- Fast Galois Field Arithmetic Library -> http://web.eecs.utk.edu/~plank/plank/papers/CS-07-593/




ISTRUZIONI PER COMPILAZIONE ED ESECUZIONE PROGRAMMA

Tutto il programma è contenuto nel file gm.c 
Il programma deve essere eseguito su sistemi Unix (su windows manca una metodo delle librerie standard di C).

Per compilare il programma. NEW!!! rimosso fattoriale gsl, non servono più le librerie
>gcc gm.c -o gm

Il programma necessita di dati di input, per comodità si consiglia di utilizzare un file di input.

Per eseguire il programma occorre fornire i dati di input (esempio di file di input https://github.com/FRSYH/Tesi/blob/master/input.txt ) nel seguente ordine e formato.
- modulo dei coefficienti
- grado massimo raggiungibile
- numero polinomi di partenza
- tipo di ordinamento  <-- NEW!!!
- elenco delle variabili utilizzate nei polinomi (su una sola riga in ordine alfabetico)
- elenco di polinomi (uno per riga, con il formato utilizzato da Magma, esempio nel file)
>./gm < input.txt > output.txt


N.B. se l'input è in un formato scorretto il programma abortirà la computazione.

N.B.2 il parser che si occupa della lettura dell'input non è ancora testato in modo ottimale, si consiglia sempre di confrontare il risultato fornito con quello di Magma.

Il parser risulta comodo per inserire i dati, tuttavia il tempo impiegato per questa operazione non è da sottovalutare.

ORDINAMENTO POLINOMI
L'ordinamento dei polinomi è ora parametrico, ed è possibile scegliere in fase di input quale utilizzare inserendo l'apposito codice.

ordinamenti implementati e relativi codici di input:
- grevlex_comparison -> 0