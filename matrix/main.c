#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main()
{
    double** proba2=beolvas("proba2.txt");
    double** proba3=beolvas("proba3.txt");
    double** negyzetes=beolvas("negyzetes.txt");
    double** pici=beolvas("pici.txt");
    double** mindenes=mult(pici,inverse(pici));
    double** ujrabeolvas=beolvas("kiir.txt");
    cout(pici);
    printf("\npici det:%g\n",det(pici));
    printf("\npici inverze:\n");
    cout(inverse(pici));
    printf("\npici * inverze:\n");
    cout(mindenes);
    letisztaz(mindenes);
    printf("\npici * inverze letisztazva:\n");
    cout(mindenes);
    printf("\npici + pici:\n");
    cout(sum(pici,pici));
    printf("\npici - pici:\n");
    cout(sub(pici,pici));
    printf("\nGauss Jordan:\n");
    cout(GJE(pici));
    printf("\nFajlba kiirt:\n");
    fkiir("kiir.txt",proba2);
    cout(ujrabeolvas);
    freem(ujrabeolvas);
    freem(proba2);
    freem(proba3);
    freem(negyzetes);
    freem(pici);
    freem(mindenes);
    return 0;
}
