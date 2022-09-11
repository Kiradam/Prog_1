#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"


void fkiir(char* hely, double** mit){
    /*
    Param�terek:
    char* hely: a kimeneti f�jl helye, neve (pl: c:\\users\\proba.txt)
    double** mit: az adott m�trixra mutat� pointer
    V�ltoz�k:
    int oszlopszam: oszlopok sz�m�t t�rolja
    int sorszam: sorok sz�m�t t�rolja
    int i,j: indexv�ltoz�k
    M�k�d�se:
    Ciklusok �ltal az adott f�jlba ki�rja a m�trixot olyan form�ban, ahogy a t�bbi f�ggv�ny azt haszn�lni tudja.
    */
    int oszlopszam=mit[0][1];
    int sorszam=mit[0][0];
    int i=1;
    int j=0;
    FILE *ptr;
    ptr = fopen(hely,"w");
    fprintf(ptr,"%d %d\n",sorszam,oszlopszam);
    while (i!=sorszam+1){
        j=0;
        while(j!=oszlopszam){
            if(j==oszlopszam-1){
                fprintf(ptr,"%g",mit[i][j]);
            }
            else{
            fprintf(ptr,"%g ",mit[i][j]);}
            j++;
        }
        fprintf(ptr,"\n");
        i++;
    }
    fclose(ptr);
}
void Mmalloc(double** M,int i, int j){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a leend� m�trix
    int i, j: sor �s oszlopsz�m
    V�ltoz�k:
    int k: ciklusv�ltoz�
    M�k�d�se:
    Lefoglal a mem�ri�ban egy megfelel� form�tum� �s m�ret� blokkot a m�trixnak
    */
    int k=1;
    M=(double**)malloc((i+1)*sizeof(double*));
    M[0]=(double*)malloc(2*sizeof(double));
    while (k!=i+1){
        M[k]=(double*)malloc(j*sizeof(double));
        k++;
    }
}
double** GJE(double** M){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a Gauss elimin�land� m�trix
    V�ltoz�k:
    double** eredmeny: eredm�ny m�trix t�rol�s�ra szolg�l, ezt adja vissza a f�ggv�ny
    double** sorv: f�ggv�nyek sz�m�ra kezelhet� form�tum�, val�j�ban egy sort tartalmaz� v�ltoz�, mely v�g�l az eredm�ny m�trix sora lesz
    double** kisebbitendo: f�ggv�nyek sz�m�ra kezelhet� form�tum�, val�j�ban egy sort tartalmaz� v�ltoz�, az adott v�ltoztatand� sort tartalmazza
    double** kivonando: f�ggv�nyek sz�m�ra kezelhet� form�tum�, val�j�ban egy sort tartalmaz� v�ltoz�, ennek alkalmas konstans-szorosai ker�lnek kivon�sra 'kisebbitendo'-b�l
    int i,j,b: ciklusv�ltoz�k
    double a: adott sor konstans szorz�ja az elimin�l�s sor�n
    M�k�d�se:
    Teljesen lem�solja a param�terk�nt megadott m�trixot, �gy az v�ltozatlan marad a folyamat v�g�re.
    El�sz�r az els� sorban megkeresi az els� nem 0 elemet, majd a t�bbi sorb�l kivonja ezen kiemelt sor megfelel� konstansszoros�t, hogy
    a kiemelt sor els� nem 0 eleme "alatt"(�s felett) csak 0-k legyenek. Az algoritmust ism�tli a 2., 3., stb. sorokra is. Amint az algoritmus befejez�d�tt, az eredm�ny m�trixb�l m�r
    leolvashat� a sort�r egy b�zisa, �m ezeket m�g soronk�nt az els� nem 0 elem szerint norm�lja. A v�g�n visszat�r a v�geredm�ny m�trixra mutat� pointerrel.
    */
    double** eredmeny=deepcpy(M);
    double** sorv;
    double** kisebbitendo;
    double** kivonando;
    int i=1;
    int j=0;
    int b=1;
    double a;
    while (i!=eredmeny[0][0]+1){
        while(b!=eredmeny[0][0]+1){
            j=0;
            while(eredmeny[i][j]==0&&j!=eredmeny[0][1]){
                j++;
            }
            if(b==i){
                //semmi
            }
            else{
                if(eredmeny[i][j]!=0){
                    a=eredmeny[b][j]/eredmeny[i][j];}
                else{
                    a=1;
                }
                kivonando=kulonsorra(eredmeny,i);
                kisebbitendo=kulonsorra(eredmeny,b);
                mulc(kivonando,a);
                free(eredmeny[b]);
                sorv=sub(kisebbitendo,kivonando);
                eredmeny[b]=(double*)malloc(sizeof(double)*eredmeny[0][1]);
                eredmeny[b]=sorv[1];
                freem(kisebbitendo);
                freem(kivonando);
                free(sorv[0]);
            }
            b++;
        }
        b=1;
        i++;
    }
    i=1;
    j=0;
    while(i!=eredmeny[0][0]+1){
        j=0;
        while(eredmeny[i][j]==0&&j!=eredmeny[0][1]){
            j++;
        }
        if(j==eredmeny[0][1]){
            break;
        }
        a=1/eredmeny[i][j];
        kisebbitendo=kulonsorra(eredmeny,i);
        mulc(kisebbitendo,a);
        free(eredmeny[i]);
        eredmeny[i]=kisebbitendo[1];
        free(kisebbitendo[0]);
        i++;
    }
    return eredmeny;
}
double** deepcpy(double** M){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a m�soland� m�trix
    V�ltoz�k:
    int i,j: ciklusv�ltoz�k
    double** deep: v�geredm�ny m�trixra mutat� pointer
    M�k�d�se:
    Lefoglal a mem�ri�ban egy ugyan olyan form�tum� �s m�ret� blokkot a m�trixnak, mint a param�terben megadott, majd a m�soland� m�trix elem�t egyenk�nt be�rja a megfelel� c�mbe. A v�g�n visszat�r a lem�solt m�trixra mutat� pointerrel.
    */
    double** deep=(double**)malloc(sizeof(double*)*(M[0][0]+1));
    int i=1;
    int j=0;
    deep[0]=(double*)malloc(2*sizeof(double));
    deep[0][0]=M[0][0];
    deep[0][1]=M[0][1];
    while (i!=deep[0][0]+1){
        j=0;
        deep[i]=(double*)malloc(sizeof(double)*deep[0][1]);
        while(j!=deep[0][1]){
            deep[i][j]=M[i][j];
            j++;
        }
        i++;
    }
    return deep;
}
void letisztaz(double** M){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a tiszt�zand� m�trix
    V�ltoz�k:
    int i,j: ciklusv�ltoz�k
    M�k�d�se:
    K�nyelmi szerepe van, az elimin�ci�k �s m�veletek sor�n double �rt�k korl�tolts�g�b�l ad�d�an el�fordulhat, hogy 0 helyett egy igen kicsi sz�mot rak a m�trixba. Ezen f�ggv�ny 1*10^(-10)-n�l kisebb 0-t�l val� elt�r�s eset�n az adott �rt�ket 0-ra �ll�tja a m�trixban.
    */
    int i=1;
    int j=0;
    while(i!=M[0][0]+1){
        j=0;
        while(j!=M[0][1]){
            if((M[i][j]-round(M[i][j]))<1.0e-010||(M[i][j]-round(M[i][j]))>-1.0e-010){
                if(round(M[i][j])==-0){
                    M[i][j]=0;}
                else{
                    M[i][j]=round(M[i][j]);
                }
                }
            j++;
        }
        i++;
    }
}
double** inverse(double** M){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, az invert�land�
    V�ltoz�k:
    double a: M determin�ns�nak reciproka
    double** eredmeny: eredm�nym�trix pointere
    M�k�d�se:
    Kisz�m�tja a m�trix inverz�t determin�nsos m�dszerrel: determin�ns reciprok�val beszorozza az adjung�lt m�trix transzpon�ltj�t, majd visszat�r az eredm�nnyel.
    */
    double a=1/det(M);
    if (det(M)==0){exit(1);}
    double** eredmeny=transponalt(adj(M));
    mulc(eredmeny,a);
    //letisztaz(eredmeny);
    return eredmeny;
}
double** adj(double** M){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, m�trix melynek az aldetermin�nsaib�l �ll� m�trixra vagyunk k�v�ncsiak
    V�ltoz�k:
    int i,j: ciklusv�ltoz�k
    double** inverse: leend� adjung�lt m�trix
    double ertek: k�ztes v�ltoz�, az adjung�lt m�trix elemeinek be�r�s�ra szolg�l
    M�k�d�se:
    Nem n�gyzetes m�trix eset�n kil�p (1). Lefoglal egy megfelel� m�ret� t�mb�t a mem�ri�ban az inverz m�trixnak. Ezut�n felt�lti elemenk�nt az eredeti m�trix megfelel� el�jel� aldetermin�nsaival. A v�g�n visszat�r az inverz m�trixra mutat� pointerrel.
    */
    double** inverse=NULL;
    int i=1;
    int j=0;
    double ertek;
    if(M[0][0]!=M[0][1]){
        printf("Wrong dimensions.");
        exit(1);
    }
    else{
        inverse=(double**)malloc((M[0][0]+1)*sizeof(double*));
        inverse[0]=(double*)malloc(2*sizeof(double));
        inverse[0][0]=M[0][0];
        inverse[0][1]=M[0][1];
        while (i!=M[0][1]+1){
            j=0;
            inverse[i]=(double*)malloc(M[0][1]*sizeof(double));
            while(j!=M[0][1]){
                ertek=det(almatrix(M,j,i-1));
                inverse[i][j]=ertek*pow(-1,(i-1)+j);
                j++;
            }
            i++;
        }
        return inverse;
    }

}
void sorcsere(double** M, int i, int j){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a cser�lend� sorokat tartalmaz� pointer
    int i, j: sor �s oszlopsz�m
    V�ltoz�k:
    double* koztes: ideiglenes pointer, mely r�vid ideig a cser�lend� sor �rt�k�vel b�r
    M�k�d�se:
    i-edik sort ideiglenesen a koztes pointerbe teszi, majd fel�l�rja a j-edik sorral. V�g�l a j-edik sort fel�l�rja a koztes sorral, azaz az eredeti i-edik sorral. A v�g�n felszabad�tja a koztes mem�ri�j�t.
    */
    double* koztes;
    koztes=M[i];
    M[i]=M[j];
    M[j]=koztes;
    koztes=NULL;
    free(koztes);
}
double** transponalt(double** M){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a transzpon�land� m�trix
    V�ltoz�k:
    int i,j: ciklusv�ltoz�k
    M�k�d�se:
    Lefoglal a mem�ri�ban egy megfelel� form�tum� �s m�ret� blokkot a m�trixnak, majd ciklusok r�v�n transzpon�lja a m�trixot az elemek megfelel� helyre val� �r�s�val. V�g�n visszaadja az eredm�nym�trixra mutat� pointert.
    */
    double** transp=(double**)malloc(sizeof(double*)*(M[0][1]+1));
    int i=1;
    int j=0;
    transp[0]=(double*)malloc(2*sizeof(double));
    transp[0][0]=M[0][1];
    transp[0][1]=M[0][0];
    while (i!=transp[0][0]+1){
        j=0;
        transp[i]=(double*)malloc(sizeof(double)*transp[0][1]);
        while(j!=transp[0][1]){
            transp[i][j]=M[j+1][i-1];
            j++;
        }
        i++;
    }
    return transp;
}
double** kulonsorra(double** M, int sor){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a m�trix melyb�l ki akarunk szervezni egy adott sort
    int sor: kiszervezend� sor indexe
    V�ltoz�k:
    double** sorvektor: f�ggv�nyek �ltali kezel�sre alkalmas form�tum�, a m�trix kiszervezend� sor�t tartalmaz� 'm�trix'
    int i: ciklusv�ltoz�
    M�k�d�se:
    Lefoglal a mem�ri�ban egy megfelel� form�tum� �s m�ret� blokkot a m�trixnak, majd a kiszervezend� sor elemeit belepakolja, �s visszat�r az eredm�nyre mutat� pointerre.
    */
    double** sorvektor=(double**)malloc(2*sizeof(double*));
    int i=0;
    sorvektor[0]=(double*)malloc(2*sizeof(double));
    sorvektor[0][0]=1;
    sorvektor[0][1]=M[0][1];
    sorvektor[1]=(double*)malloc((sorvektor[0][1])*sizeof(double));
    while (i!=sorvektor[0][1]){
        sorvektor[1][i]=M[sor][i];
        i++;
    }
    return sorvektor;
}
double** almatrix(double** M, int oszlop, int sor){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a csonk�tand� m�trix
    int oszlop, sor: oszlopindex, sorindex, ezeket v�gja ki a m�trixb�l
    V�ltoz�k:
    double ujsorszam: �j m�trix sorsz�ma
    double ujoszlopszam: �j m�trix oszlopsz�ma
    int i,j,k: ciklusv�ltoz�k
    double** ujmatrix: eredm�ny m�trix pointere
    M�k�d�se:
    Lefoglal a mem�ri�ban egy megfelel� form�tum� �s m�ret� blokkot a m�trixnak, majd felt�lti az eredeti m�trix azon elemeivel megfelel� sorrendben, melyekre igaz, hogy oszlopindex�k nem egyenl� a param�terben megadott 'oszlop'-pal �s sorindex�k nem egyenl� az ugyanitt megadott 'sor'-ral.
    A v�g�n visszat�r a csonk�tott m�trixra mutat� pointerrel.
    */
    double ujsorszam=M[0][0]-1;
    double ujoszlopszam=M[0][1]-1;
    int i=1;
    int j=0;
    int k=0;
    double** ujmatrix;
    ujmatrix=(double**)malloc((ujsorszam+1)*sizeof(double*));
    ujmatrix[0]=(double* )malloc(sizeof(double)*2);
    ujmatrix[0][0]=ujsorszam;
    ujmatrix[0][1]=ujoszlopszam;
    while (i-1!=ujsorszam+1){
        if(i<sor+1){
            j=0;
            k=0;
            ujmatrix[i]=malloc((ujoszlopszam)*sizeof(double));
            while (k!=M[0][1]){
                if (k!=oszlop){
                    ujmatrix[i][j]=M[i][k];
                    j++;
                }
                k++;
            }
        }
        else if (i>sor+1){
            j=0;
            k=0;
            ujmatrix[i-1]=malloc((ujoszlopszam)*sizeof(double));
            while (k!=M[0][1]){
                if (k!=oszlop){
                    ujmatrix[i-1][j]=M[i][k];
                    j++;
                }
                k++;
            }
        }
        i++;
    }
    return ujmatrix;
}
double det(double** A){
    /*
    Param�terek:
    double** A: pointerre mutat� pointer, a m�trix melynek determin�ns�t szeretn�nk tudni
    V�ltoz�k:
    int i: ciklusv�ltoz�
    double determinant: leend� determin�ns
    double** atmeneti: �tmeneti m�trix mely csak a rekurz�v determin�nssz�m�t�s k�ztes alm�trixainak t�rol�s�ra val�
    M�k�d�se:
    Ha a m�trix nem n�gyzetes, hib�t dob, �s kil�p(1).Ellenkez� esetben ha a m�trix 1x1-es, visszat�r az adott elemmel, ha nxn-es (n>=2), a fels� sorra alkalmazott kifejt�si szab�ly szerint rekurz�van kisz�molja a determin�nst.
    K�zben az �tmeneti m�trixokat sorra felszabad�tja. A v�g�n visszat�r az adott val�s sz�mmal.
    */
    int i=0;
    double** atmeneti;
    double determinant=0;
    if (A[0][0]!=A[0][1]){
        printf("It has no determinant.");
        exit(1);
    }
    else{
        if(A[0][0]==1){
            return A[1][0];
        }
        else{
            while(i!=A[0][0]){
                if (i/2==0&&i!=0){
                    atmeneti=almatrix(A,i,0);
                    determinant-=A[1][i]*det(atmeneti);
                    freem(atmeneti);
                }
                else{
                    atmeneti=almatrix(A,i,0);
                    determinant+=(A[1][i])*(det(atmeneti));
                    freem(atmeneti);
                }
                i++;
            }
            return determinant;
        }
    }
}

void freem(double** M){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a m�trix aminek a hely�t fel szeretn�nk szabad�tani
    V�ltoz�k:
    int sor: m�trix sorsz�ma
    int i: ciklusv�ltoz�
    M�k�d�se:
    El�sz�r a val�di �rt�kkel b�r� sorokat szabad�tja fel, majd a m�trix m�ret�t tartalmaz� 0. sort is.
    */
    int sor=M[0][0];
    int i=1;
    while (i!=sor+1){
        free(M[i]);
        i++;
    }
    free(M[0]);
    free(M);
}
double** mult(double** M1, double** M2){
    /*
    Param�terek:
    double** M1: pointerre mutat� pointer, a szorz�sban bal oldali m�trix
    double** M2: pointerre mutat� pointer, a szorz�sban jobb oldali m�trix
    V�ltoz�k:
    int sorelem: sorban fut� index
    int oszlop: oszlopban fut� index
    double a: k�ztes v�ltoz�, k�t m�trix megfelel� elemeinek szorzata itt ad�dik �ssze, ez lesz a szorzatm�trix eleme
    int j: ciklusv�ltoz�
    double** result: eredm�ny m�trix
    M�k�d�se:
    Nem megfelel� "alak�" m�trixok eset�n hib�t dob. Egy�bk�nt lefoglal az eredm�nym�trixnak megfelel� m�ret� t�mb�t �s a m�trixszorz�s szab�ly�nak megfelel�en ciklusok r�v�n felt�lti az eredm�ny m�trixot, majd visszaadja a r� mutat� pointert.
    */
    int sorelem=0;
    int oszlop=1;
    double a=0;
    int j=0;
    double** result=(double**)malloc((M1[0][0]+1)*sizeof(double*));
    result[0]=(double*)malloc(2*sizeof(double));
    result[0][0]=M1[0][0];
    result[0][1]=M2[0][1];
    if (M1[0][1]!=M2[0][0]){
        printf("Cannot be multiplicated.");
        freem(result);
        exit(1);
    }
    else{
        while(oszlop!=M1[0][0]+1){
            sorelem=0;
            result[oszlop]=(double*)malloc(M2[0][1]*sizeof(double));
            while(sorelem!=M2[0][1]){
                while(j!=M1[0][1]){
                    a+=M1[oszlop][j]*M2[j+1][sorelem];
                    j++;
                }
                result[oszlop][sorelem]=a;
                sorelem++;
                j=0;
                a=0;
            }
            oszlop++;
        }
    }
    return result;
}
void mulc(double** M, double lambda){
    /*
    Param�terek:
    double** M: pointerre mutat� pointer, a konstanssal szorzand� m�trix
    double lambda: a konstans amivel szorzunk
    V�ltoz�k:
    int sori: sorban fut� index
    int oszlopi: oszlopban fut� index
    M�k�d�se:
    V�gigmegy a param�terben megadott m�trix osszes elem�n �s beszorozza az adott lambda konstanssal.
    */
    int sori=0;
    int oszlopi=1;
    while (oszlopi!=M[0][0]+1){
            sori=0;
            while(sori!=M[0][1]){
                M[oszlopi][sori]*=lambda;
                sori++;
            }
            oszlopi++;
    }
}
double** sub(double** M1, double** M2){
    /*
    Param�terek:
    double** M1: pointerre mutat� pointer, a kisebb�tend� m�trix
    double** M2: pointerre mutat� pointer, a kivonand� m�trix
    V�ltoz�k:
    int sori: sorban fut� index
    int oszlopi: oszlopban fut� index
    double a: k�ztes v�ltoz�, k�t m�trix megfelel� elemeinek k�l�nbs�ge, ez lesz a k�l�nbs�gm�trix eleme
    double** kul: eredm�ny m�trix
    M�k�d�se:
    Nem azonos "alak�" m�trixok eset�n hib�t dob(1). Ellenkez� esetben elemenk�nt felt�lti az eredem�ny m�trixot a k�t param�term�trix megfelel� elemeinek k�l�nbs�g�vel, majd visszat�rt az eredem�nym�trixra mutat� pointerrel.
    */
    double** kul=(double**)malloc((M1[0][0]+1)*sizeof(double*));
    int sori=0;
    int oszlopi=1;
    double c;
    if ((M1[0][0]!=M2[0][0])||M1[0][1]!=M1[0][1]){
        printf("wrong sizes");
        free(kul);
        exit(1);
    }
    else {
        kul[0]=(double*)malloc(2*sizeof(double));
        kul[0][0]=M1[0][0];
        kul[0][1]=M1[0][1];
        while (oszlopi!=M1[0][0]+1){
            sori=0;
            kul[oszlopi]=(double*)malloc(M1[0][1]*sizeof(double));
            while(sori!=M1[0][1]){
                c=M1[oszlopi][sori]-M2[oszlopi][sori];
                kul[oszlopi][sori]=c;
                sori++;
            }
            oszlopi++;
        }
    }
    return kul;
}
double** sum(double** M1, double** M2){
    /*
    Param�terek:
    double** M1: pointerre mutat� pointer, az egyik �sszeadand� m�trix
    double** M2: pointerre mutat� pointer, a m�sik �sszeadand� m�trix
    V�ltoz�k:
    int sori: sorban fut� index
    int oszlopi: oszlopban fut� index
    double c: k�ztes v�ltoz�, k�t m�trix megfelel� elemeinek �sszege, ez lesz a �sszegm�trix eleme
    double** osszeg: eredm�ny m�trix
    M�k�d�se:
    Nem azonos "alak�" m�trixok eset�n hib�t dob(1). Ellenkez� esetben elemenk�nt felt�lti az eredem�ny m�trixot a k�t param�term�trix megfelel� elemeinek �sszeg�vel, majd visszat�rt az eredem�nym�trixra mutat� pointerrel.
    */
    double** osszeg=(double**)malloc((M1[0][0]+1)*sizeof(double*));
    int sori=0;
    int oszlopi=1;
    double c;
    if ((M1[0][0]!=M2[0][0])||M1[0][1]!=M1[0][1]){
        printf("wrong sizes");
        free(osszeg);
        exit(1);
    }
    else {
        osszeg=(double**)malloc((1+M1[0][0])*sizeof(double*));
        osszeg[0]=(double*)malloc(2*sizeof(double));
        osszeg[0][0]=M1[0][0];
        osszeg[0][1]=M1[0][1];
        while (oszlopi!=M1[0][0]+1){
            sori=0;
            osszeg[oszlopi]=(double*)malloc(M1[0][1]*sizeof(double));
            while(sori!=M1[0][1]){
                c=M1[oszlopi][sori]+M2[oszlopi][sori];
                osszeg[oszlopi][sori]=c;
                sori++;
            }
            oszlopi++;
        }
    }
    return osszeg;
}
void cout(double** matrix){
    /*
    Param�terek:
    double** matrix: pointerre mutat� pointer, a konzolra ki�rand� m�trix
    V�ltoz�k:
    int sori: sorban fut� index
    int oszlopi: oszlopban fut� index
    int sor: param�term�trix sorainak hossza
    int oszlop: param�term�trix oszlopainak hossza
    M�k�d�se:
    A m�trix elemeit egyes�vel, 4 sz�k�zzel elv�lasztva, sorokat 2 sort�r�ssel tagolva ki�rja a konzolablakba. A m�ret�t (azaz a 0. sort) nem �rja ki.
    */
    int sori=0;
    int oszlopi=0;
    int sor=matrix[0][1];
    int oszlop=matrix[0][0];
    while (oszlopi!=oszlop){
        while (sori!=sor){
             printf("%g    ",matrix[oszlopi+1][sori]);
            sori++;
        }
        sori=0;
        printf("\n\n");
        oszlopi++;
    }
}
double** beolvas(char* hely){
    /*
    Param�terek:
    char* hely: f�jl el�r�si helye
    V�ltoz�k:
    char szoveg[1000]: 1000 karakter hossz�s�gig beolvassa a f�jlb�l az adott sort
    double* sorolv: 'szoveg'-b�l k�sz�tett double t�mb, mely a m�trix sora lesz
    double** matrix: v�geredm�ny m�trix, pointerre mutat� pointer
    int sor: sorsz�m
    int oszlop: oszlopsz�m
    int i: ciklusv�ltoz�
    M�k�d�se:
    Megnyitja az adott sz�vegf�jlt. Az els� sorb�l, mely a m�trix m�ret�t tartalmazza, �t�ll�tja ennek megfelel�en 'sor' �s 'oszlop' �rt�keit. Ezut�n ezek f�ggv�ny�ben ism�tli az algoritmust, csak m�r a val�di m�trixelemeket tartalmaz� sorokra, �s az �gy k�sz�tett double t�mb�ket a megfelel� 'matrix' c�mhez t�rs�tja. A v�g�n visszat�r a beolvasott m�trixra mutat� pointerrel.
    */
    char szoveg[1000];
    double* sorolv;
    double** matrix;
    int sor;
    int oszlop;
    int i=1;
    FILE *fptr;
    fptr = fopen(hely,"r");
    fgets(szoveg,1000,fptr);
    sorolv=(double*)malloc(2*sizeof(double));
    sorolv=soronkent(szoveg);
    sor=sorolv[0];
    oszlop=sorolv[1];
    matrix=(double**)malloc((sor+1)*sizeof(double*));
    matrix[0]=sorolv;
    while (i!=sor+1){
        sorolv=(double*)malloc(oszlop*sizeof(double));
        fgets(szoveg,1000,fptr);
        sorolv=soronkent(szoveg);
        matrix[i]=sorolv;
        i++;
    }
    fclose(fptr);
    return matrix;
}
double* soronkent(char* sor){
    /*
    Param�terek:
    char* sor: karaktert�mbre mutat� pointer
    V�ltoz�k:
    int space: sz�k�z�ket sz�mol� v�ltoz�
    int i: ciklusv�ltoz�
    int j: ciklusv�ltoz�
    int szamlalo: ciklusv�ltoz�
    int elozo: ciklusv�ltoz�
    char* szamma: karakterpointer, k�s�bb ezt alak�tja sz�mm�
    double* elsosor: eredm�ny pointer, a karaktert�mbb�l sz�mm� alak�tott �rt�keket tartalmazza
    M�k�d�se:
    Megsz�molja az adott stringben tal�lhat� sz�k�z�ket. Ennek f�ggv�ny�ben megfelel� m�ret� mem�ri�t foglal 'elsosor'-nak. El�sz�r az els� sz�k�zig beolvassa az adatot, azt sz�mm� konvert�lja �s be�rja 'elsosor' megfelel� c�m�be. Megjegyzi az el�z� sz�k�z helyzet�t, �gy az algoritmust el�r�l kezdve a k�vetkez� sz�mot �rja 'elsosor'-ba, eg�szen a string v�g�ig. V�g�l visszat�r az adott double t�mbbel.
    */
    int space=0;
    int i=0;
    int j=0;
    int szamlalo=0;
    int elozo=0;
    char* szamma;
    double* elsosor;
    while (i<strlen(sor)){
        if (sor[i]==' '){
            space++;
        }
        i++;
    }
    i=0;
    elsosor=(double*)malloc((strlen(sor)-space)*sizeof(double));
    while (szamlalo!=strlen(sor)+1){
        if (sor[szamlalo]==' ' || sor[szamlalo]=='\n' || sor[szamlalo]=='\0'){
            szamma=malloc(szamlalo-elozo+1);;
            i=0;
            while (elozo!=szamlalo){
                szamma[i]=sor[elozo];
                i++;
                elozo++;
            }
            szamma[i]='\0';
            elsosor[j]=atof(szamma);
            //printf("%g      %d\n",elsosor[j],j);
            j++;
            free(szamma);
            elozo=szamlalo+1;
        }
        szamlalo++;
    }
    return elsosor;
}
