#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"


void fkiir(char* hely, double** mit){
    /*
    Paraméterek:
    char* hely: a kimeneti fájl helye, neve (pl: c:\\users\\proba.txt)
    double** mit: az adott mátrixra mutató pointer
    Változók:
    int oszlopszam: oszlopok számát tárolja
    int sorszam: sorok számát tárolja
    int i,j: indexváltozók
    Mûködése:
    Ciklusok által az adott fájlba kiírja a mátrixot olyan formában, ahogy a többi függvény azt használni tudja.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a leendõ mátrix
    int i, j: sor és oszlopszám
    Változók:
    int k: ciklusváltozó
    Mûködése:
    Lefoglal a memóriában egy megfelelõ formátumú és méretû blokkot a mátrixnak
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
    Paraméterek:
    double** M: pointerre mutató pointer, a Gauss eliminálandó mátrix
    Változók:
    double** eredmeny: eredmény mátrix tárolására szolgál, ezt adja vissza a függvény
    double** sorv: függvények számára kezelhetõ formátumú, valójában egy sort tartalmazó változó, mely végül az eredmény mátrix sora lesz
    double** kisebbitendo: függvények számára kezelhetõ formátumú, valójában egy sort tartalmazó változó, az adott változtatandó sort tartalmazza
    double** kivonando: függvények számára kezelhetõ formátumú, valójában egy sort tartalmazó változó, ennek alkalmas konstans-szorosai kerülnek kivonásra 'kisebbitendo'-bõl
    int i,j,b: ciklusváltozók
    double a: adott sor konstans szorzója az eliminálás során
    Mûködése:
    Teljesen lemásolja a paraméterként megadott mátrixot, így az változatlan marad a folyamat végére.
    Elõször az elsõ sorban megkeresi az elsõ nem 0 elemet, majd a többi sorból kivonja ezen kiemelt sor megfelelõ konstansszorosát, hogy
    a kiemelt sor elsõ nem 0 eleme "alatt"(és felett) csak 0-k legyenek. Az algoritmust ismétli a 2., 3., stb. sorokra is. Amint az algoritmus befejezõdött, az eredmény mátrixból már
    leolvasható a sortér egy bázisa, ám ezeket még soronként az elsõ nem 0 elem szerint normálja. A végén visszatér a végeredmény mátrixra mutató pointerrel.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a másolandó mátrix
    Változók:
    int i,j: ciklusváltozók
    double** deep: végeredmény mátrixra mutató pointer
    Mûködése:
    Lefoglal a memóriában egy ugyan olyan formátumú és méretû blokkot a mátrixnak, mint a paraméterben megadott, majd a másolandó mátrix elemét egyenként beírja a megfelelõ címbe. A végén visszatér a lemásolt mátrixra mutató pointerrel.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a tisztázandó mátrix
    Változók:
    int i,j: ciklusváltozók
    Mûködése:
    Kényelmi szerepe van, az eliminációk és mûveletek során double érték korlátoltságából adódóan elõfordulhat, hogy 0 helyett egy igen kicsi számot rak a mátrixba. Ezen függvény 1*10^(-10)-nél kisebb 0-tól való eltérés esetén az adott értéket 0-ra állítja a mátrixban.
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
    Paraméterek:
    double** M: pointerre mutató pointer, az invertálandó
    Változók:
    double a: M determinánsának reciproka
    double** eredmeny: eredménymátrix pointere
    Mûködése:
    Kiszámítja a mátrix inverzét determinánsos módszerrel: determináns reciprokával beszorozza az adjungált mátrix transzponáltját, majd visszatér az eredménnyel.
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
    Paraméterek:
    double** M: pointerre mutató pointer, mátrix melynek az aldeterminánsaiból álló mátrixra vagyunk kíváncsiak
    Változók:
    int i,j: ciklusváltozók
    double** inverse: leendõ adjungált mátrix
    double ertek: köztes változó, az adjungált mátrix elemeinek beírására szolgál
    Mûködése:
    Nem négyzetes mátrix esetén kilép (1). Lefoglal egy megfelelõ méretû tömböt a memóriában az inverz mátrixnak. Ezután feltölti elemenként az eredeti mátrix megfelelõ elõjelû aldeterminánsaival. A végén visszatér az inverz mátrixra mutató pointerrel.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a cserélendõ sorokat tartalmazó pointer
    int i, j: sor és oszlopszám
    Változók:
    double* koztes: ideiglenes pointer, mely rövid ideig a cserélendõ sor értékével bír
    Mûködése:
    i-edik sort ideiglenesen a koztes pointerbe teszi, majd felülírja a j-edik sorral. Végül a j-edik sort felülírja a koztes sorral, azaz az eredeti i-edik sorral. A végén felszabadítja a koztes memóriáját.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a transzponálandó mátrix
    Változók:
    int i,j: ciklusváltozók
    Mûködése:
    Lefoglal a memóriában egy megfelelõ formátumú és méretû blokkot a mátrixnak, majd ciklusok révén transzponálja a mátrixot az elemek megfelelõ helyre való írásával. Végén visszaadja az eredménymátrixra mutató pointert.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a mátrix melybõl ki akarunk szervezni egy adott sort
    int sor: kiszervezendõ sor indexe
    Változók:
    double** sorvektor: függvények általi kezelésre alkalmas formátumú, a mátrix kiszervezendõ sorát tartalmazó 'mátrix'
    int i: ciklusváltozó
    Mûködése:
    Lefoglal a memóriában egy megfelelõ formátumú és méretû blokkot a mátrixnak, majd a kiszervezendõ sor elemeit belepakolja, és visszatér az eredményre mutató pointerre.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a csonkítandó mátrix
    int oszlop, sor: oszlopindex, sorindex, ezeket vágja ki a mátrixból
    Változók:
    double ujsorszam: új mátrix sorszáma
    double ujoszlopszam: új mátrix oszlopszáma
    int i,j,k: ciklusváltozók
    double** ujmatrix: eredmény mátrix pointere
    Mûködése:
    Lefoglal a memóriában egy megfelelõ formátumú és méretû blokkot a mátrixnak, majd feltölti az eredeti mátrix azon elemeivel megfelelõ sorrendben, melyekre igaz, hogy oszlopindexük nem egyenlõ a paraméterben megadott 'oszlop'-pal és sorindexük nem egyenlõ az ugyanitt megadott 'sor'-ral.
    A végén visszatér a csonkított mátrixra mutató pointerrel.
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
    Paraméterek:
    double** A: pointerre mutató pointer, a mátrix melynek determinánsát szeretnénk tudni
    Változók:
    int i: ciklusváltozó
    double determinant: leendõ determináns
    double** atmeneti: átmeneti mátrix mely csak a rekurzív determinánsszámítás köztes almátrixainak tárolására való
    Mûködése:
    Ha a mátrix nem négyzetes, hibát dob, és kilép(1).Ellenkezõ esetben ha a mátrix 1x1-es, visszatér az adott elemmel, ha nxn-es (n>=2), a felsõ sorra alkalmazott kifejtési szabály szerint rekurzívan kiszámolja a determinánst.
    Közben az átmeneti mátrixokat sorra felszabadítja. A végén visszatér az adott valós számmal.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a mátrix aminek a helyét fel szeretnénk szabadítani
    Változók:
    int sor: mátrix sorszáma
    int i: ciklusváltozó
    Mûködése:
    Elõször a valódi értékkel bíró sorokat szabadítja fel, majd a mátrix méretét tartalmazó 0. sort is.
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
    Paraméterek:
    double** M1: pointerre mutató pointer, a szorzásban bal oldali mátrix
    double** M2: pointerre mutató pointer, a szorzásban jobb oldali mátrix
    Változók:
    int sorelem: sorban futó index
    int oszlop: oszlopban futó index
    double a: köztes változó, két mátrix megfelelõ elemeinek szorzata itt adódik össze, ez lesz a szorzatmátrix eleme
    int j: ciklusváltozó
    double** result: eredmény mátrix
    Mûködése:
    Nem megfelelõ "alakú" mátrixok esetén hibát dob. Egyébként lefoglal az eredménymátrixnak megfelelõ méretû tömböt és a mátrixszorzás szabályának megfelelõen ciklusok révén feltölti az eredmény mátrixot, majd visszaadja a rá mutató pointert.
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
    Paraméterek:
    double** M: pointerre mutató pointer, a konstanssal szorzandó mátrix
    double lambda: a konstans amivel szorzunk
    Változók:
    int sori: sorban futó index
    int oszlopi: oszlopban futó index
    Mûködése:
    Végigmegy a paraméterben megadott mátrix osszes elemén és beszorozza az adott lambda konstanssal.
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
    Paraméterek:
    double** M1: pointerre mutató pointer, a kisebbítendõ mátrix
    double** M2: pointerre mutató pointer, a kivonandó mátrix
    Változók:
    int sori: sorban futó index
    int oszlopi: oszlopban futó index
    double a: köztes változó, két mátrix megfelelõ elemeinek különbsége, ez lesz a különbségmátrix eleme
    double** kul: eredmény mátrix
    Mûködése:
    Nem azonos "alakú" mátrixok esetén hibát dob(1). Ellenkezõ esetben elemenként feltölti az eredemény mátrixot a két paramétermátrix megfelelõ elemeinek különbségével, majd visszatért az eredeménymátrixra mutató pointerrel.
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
    Paraméterek:
    double** M1: pointerre mutató pointer, az egyik összeadandó mátrix
    double** M2: pointerre mutató pointer, a másik összeadandó mátrix
    Változók:
    int sori: sorban futó index
    int oszlopi: oszlopban futó index
    double c: köztes változó, két mátrix megfelelõ elemeinek összege, ez lesz a összegmátrix eleme
    double** osszeg: eredmény mátrix
    Mûködése:
    Nem azonos "alakú" mátrixok esetén hibát dob(1). Ellenkezõ esetben elemenként feltölti az eredemény mátrixot a két paramétermátrix megfelelõ elemeinek összegével, majd visszatért az eredeménymátrixra mutató pointerrel.
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
    Paraméterek:
    double** matrix: pointerre mutató pointer, a konzolra kiírandó mátrix
    Változók:
    int sori: sorban futó index
    int oszlopi: oszlopban futó index
    int sor: paramétermátrix sorainak hossza
    int oszlop: paramétermátrix oszlopainak hossza
    Mûködése:
    A mátrix elemeit egyesével, 4 szóközzel elválasztva, sorokat 2 sortöréssel tagolva kiírja a konzolablakba. A méretét (azaz a 0. sort) nem írja ki.
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
    Paraméterek:
    char* hely: fájl elérési helye
    Változók:
    char szoveg[1000]: 1000 karakter hosszúságig beolvassa a fájlból az adott sort
    double* sorolv: 'szoveg'-bõl készített double tömb, mely a mátrix sora lesz
    double** matrix: végeredmény mátrix, pointerre mutató pointer
    int sor: sorszám
    int oszlop: oszlopszám
    int i: ciklusváltozó
    Mûködése:
    Megnyitja az adott szövegfájlt. Az elsõ sorból, mely a mátrix méretét tartalmazza, átállítja ennek megfelelõen 'sor' és 'oszlop' értékeit. Ezután ezek függvényében ismétli az algoritmust, csak már a valódi mátrixelemeket tartalmazó sorokra, és az így készített double tömböket a megfelelõ 'matrix' címhez társítja. A végén visszatér a beolvasott mátrixra mutató pointerrel.
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
    Paraméterek:
    char* sor: karaktertömbre mutató pointer
    Változók:
    int space: szóközöket számoló változó
    int i: ciklusváltozó
    int j: ciklusváltozó
    int szamlalo: ciklusváltozó
    int elozo: ciklusváltozó
    char* szamma: karakterpointer, késõbb ezt alakítja számmá
    double* elsosor: eredmény pointer, a karaktertömbbõl számmá alakított értékeket tartalmazza
    Mûködése:
    Megszámolja az adott stringben található szóközöket. Ennek függvényében megfelelõ méretû memóriát foglal 'elsosor'-nak. Elõször az elsõ szóközig beolvassa az adatot, azt számmá konvertálja és beírja 'elsosor' megfelelõ címébe. Megjegyzi az elõzõ szóköz helyzetét, így az algoritmust elõrõl kezdve a következõ számot írja 'elsosor'-ba, egészen a string végéig. Végül visszatér az adott double tömbbel.
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
