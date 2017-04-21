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


ISTRUZIONI PER COMPILAZIONE ED ESECUZIONE gm.c

Per compilare è necessario aver installato le librerie GSL
>gcc gm.c -lgsl -lgslcblas -lm -o gm


Per eseguire il programma occorre creare un file di testo contenente l'input nel formato indicato https://github.com/FRSYH/Tesi/blob/master/input.txt
- modulo
- grado massimo
- numero polinomi di partenza
- elenco delle variabili utilizzate
- elenco di polinomi (uno per riga, con il formato utilizzato)
>./gm < input.txt > output.txt
