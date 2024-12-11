// Support.c, Version 1.1.2
// (C) Winfried Bruns, 2011

/*
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version. This program is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public
License along with this program; if not, see
http://www.gnu.org/licenses/
*/


// Liest den im File argv[1] gespeicherten glatten Fächer
// und bestimmt zu ihm Support-Polytope.
// Diese werden an ChPoly zur weiteren Untersuchung übergeben.

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


#define LONGLONG

// typedef int Zahl;

#include "ToricExp.c"


void main(int argc, char **argv)
{

    int NOS, Dim,i,j,EDim,Rang,q, qmax,MaxDimDiagSplit,MaxPicRankHilbBase,Rang0;
    Matrix OS, OST,A0,A1,Sub,Sub1,Pi,Id,L,LT,C0Bas;
    Zahl Nenner,Index;
    Vektor Diag;
    VList OneSkel, Conv, Conv0, ConvB, FConv, DummySupp, HilbBasVA,VAC,
                InclMin, EssHP;
    VLE *VP, *VP1, *VPNext, *FromList;
    SimplexList Fan;
    Simplex *Smp, *Smp1;
    Indizes FF,Facet;
    FILE *Ein, *Prot, *Aus;
    boolean IsSplit,ChPoly;

    char database[1000], DiagTest[5], Uebergabe[1000];

#define NUM_OPTIONS 7

    char *DefOptions[NUM_OPTIONS]=
    {
        "Timer true",
        "MaxPicRankHilbBase 10", // 0 to switch off
        "MaxModDiagSplit 5", // 0 to switch off
        "MaxDimDiagSplit 8",
        "Verbose false",
        "ChPoly true",
        "end_options"
    };
    char **Options;

    GetAndPrintPID(argv[0],Uebergabe);

    InitFreeLists();
    RunTime();
    TransferDefOptions(&Options,DefOptions,NUM_OPTIONS);
    ReadConfigFile("./Support",Options);
    BoolOption(&VERBOSE,"Verbose",Options);
    BoolOption(&TIMER,"Timer",Options);
    BoolOption(&ChPoly,"ChPoly",Options);
    if(!ChPoly)
        VERBOSE=true;

    IntOption(&MaxPicRankHilbBase,"MaxPicRankHilbBase",Options); // Max Modul für DiagSplit
    IntOption(&MaxDimDiagSplit,"MaxDimDiagSplit",Options); // Max Dim für DiagSplit
    IntOption(&qmax,"MaxModDiagSplit",Options);


    strcat(Uebergabe,".supp");

    Ein=fopen(argv[1],"r");
    fscanf(Ein,"%d",&NOS);  // Anzahl Strahlen des Fächers
    fscanf(Ein,"%d",&Dim); // hier Polytop-Dimension

    EDim=NOS-Dim+1;
    MaxDim=NOS+1;

    Diag=(Vektor) calloc(Dim,sizeof(Zahl));


    LiesTriang(Ein,&OneSkel,NOS,Dim,&Fan);
    fclose(Ein);

    for(Smp=Fan.St;Smp;Smp=Smp->R)  // Simplizes aufsteigend sortieren
        qsort(Smp->Vert,Dim,sizeof(int),intvergl_asc);


    IsSplit=false;
   if(MaxDimDiagSplit>=Dim)
   {
    for(q=3;q<=qmax;q++)
    {
        if(q%2==0 || (q>3 && q%3==0))
            continue;
        if(VERBOSE)
            printf("Test DiagSplit for q = %d\n",q);
        if(DiagSplit(OneSkel,Dim,q))
            {
                printf("In Support diagonally split for q = %d\n",q);
                IsSplit=true;
                break;
            }
    }
    if(q>qmax && qmax>0)
        printf("In Support not diagonally split for 3 <= q <= %d\n",qmax);
  }
  else
    printf("No diag split test\n");


    MatAlloc(&A0,Dim,Dim);
    MatAlloc(&A1,Dim,Dim);

    FuelleMatrixVL(A0,Dim,Dim,OneSkel,(Fan.St)->Vert); // erste Ecke auf

    InvertZ(A1,A0,Dim,Diag,&Nenner);                 // positiven Orthanten

    for(VP=OneSkel.St;VP;VP=VP->R)
        MxV(VP->V,A1,Dim,Dim,VP->V);                 // transformieren


    FF=(Indizes) calloc(Dim,sizeof(int));  // FF steht fuer "First facet" (in diesem Fall
                                           // maximaler Kegel)
                                           // entspricht Ecke des Polytops
                                           // die in Nullpunkt geschoben wird

    Facet=(Indizes) calloc(Dim-1,sizeof(int)); // Hier wirklich Facette eines max. Kegels

    memcpy(FF,Fan.St->Vert,Dim*sizeof(int));

    MatAlloc(&OS,NOS,Dim);
    VL2Mat(OS,OneSkel,Dim);

    MatAlloc(&OST,Dim,NOS);
    TrMat(OST,OS,NOS,Dim);

    MatAlloc(&Pi,NOS,Dim);
    MatAlloc(&Sub1,NOS,Dim);
    MatAlloc(&Sub,NOS,NOS);
    MatAlloc(&Id,NOS,NOS);
    MacheEinhMat(Id,NOS);

    MatAlloc(&L,NOS,NOS);
    MatAlloc(&LT,NOS,NOS);

    StartVL(&Conv);


    for(Smp=Fan.St;Smp;Smp=Smp->R) // Mache Ungleichungen fuer Konvexitaet
    {
        FuelleMatrix(A0 ,Dim, Dim,OS, Smp->Vert);
        InvertZ(A1,A0,Dim,Diag,&Nenner); // Gibt Inverse der Transponierten zurueck
        if(Nenner!=1)
        {
            printf("Support: non unimodular fan\n");
            exit(1);
        }
        for(i=0;i<NOS;i++) // Mache Projektionsmatrix fuer Facette Smp
            for(j=0;j<Dim;j++)
                Pi[i][j]=0;
        for(j=0;j<Dim;j++)
            Pi[Smp->Vert[j]][j]=1; // gemacht

        MxM0(Sub1,Pi,NOS,Dim,A1,Dim); // Mache Subtraktionsmatrix fuer Facette Smp
        MxM0(Sub, Sub1,NOS,Dim,OST,NOS);

        for(i=0;i<NOS;i++) // Subtrahiere von Id
            for(j=0;j<NOS;j++)
                L[i][j]=Id[i][j]-Sub[i][j];

        TrMat(LT, L, NOS,NOS); // Verwandle in Linearformen
        // Mat2VL(&FConv,LT,NOS,NOS);
        // FSchreibeVL(stdout,FConv,NOS);

        for(i=0;i<Dim-1;i++)        // nur eventuell wesentliche Ungleichungen
            Facet[i]=Smp->Vert[i];  // auswaehlen. Sie gehoeren zu Elementen von
                                    // OneSkel, die mit einer Facete von Smp
                                    // einen Kegel aufspannen (!=Smp)
                                    // Facet auf Startwert gesetzt

        for(i=0;i<Dim;i++)          // Wir suchen zu jeder Facette den "Gegenkegel"
        {
            if(i>0)
                Facet[Dim-1-i]=Smp->Vert[Dim-i]; // Facette fortschreiben
            for(Smp1=Fan.St;Smp1;Smp1=Smp1->R)
            {
                if(Smp1==Smp)
                    continue;
                j=FacetOfSimplex(Smp1->Vert,Facet,Dim);
                if(j>=0)
                    break;
            }

            if(!Smp1)
            {
                printf("Support: fan not complete\n"); // Keinen "Gegenkegel" gefunden
                exit(1);
            }

            InsertV(&Conv,LT[j],NOS);
        }


        /* for(VP=FConv.St;VP;VP=VPNext)
        {
            VPNext=VP->R;
            if(NullVekt(VP->V,NOS) || VLEinVL(VP,Conv,NOS))
            {
                FreeVLE(&FConv,VP);
                continue;
            }
            MoveVLE(&Conv,&FConv,VP);
        } */
    }

    CSchreibeVL(stdout,
                "Defining hyperplanes of nef cone in PL (not pointed)",Conv,NOS);

    StartVL(&ConvB);

    VP=NewVLE(&ConvB);
    for(i=0;i<EDim;i++)
        VP->V[i]=0;
    VP->V[EDim-1]=1; // "Rees"-Variable hat positiven Exponenten

    for(VP=Conv.St;VP;VP=VP->R) // Beschraenke auf Unterraum(phi|FF =0)
    {
        VP1=NewVLE(&ConvB);
        j=0;
        for(i=0;i<NOS;i++)
            if(!InVertices(i,FF,Dim)) // Doubletten aussortieren
            {
                VP1->V[j]=VP->V[i];
                j++;
            }
        VP1->V[EDim-1]=-1; // Wir wollen das Innere bestimmen
    }

    FreeVL(&Conv); // ist jetzt nach ConvB uebertragen worden



    printf("Made hyperpl restr nef cone ext by int\n");
    CSchreibeVL(stdout,
            "Defining hyperplanes of restricted nef cone in PL, ext by int",ConvB,EDim);

    RunTime();
    fflush(stdout);

    SuppHypVL(&VAC,ConvB,EDim); // VAC = Randstrahlen zu "very ample cone"

    // FreeVL(&ConvB);

    printf("#Ext rays %d\n",NinVL(VAC));
    CSchreibeVL(stdout,
                  "Ext rays of restricted nef cone in PL, extended by int",VAC,EDim);
    RunTime();
    fflush(stdout);

    StartVL(&HilbBasVA);


    if(EDim<=MaxPicRankHilbBase+1)
    {
        /* SelExtremeRRank(&EssHP, ConvB, VAC,EDim);
        printf("#Ess Hyperpl %d\n",NinVL(EssHP));
        fflush(stdout);
        FreeVL(&ConvB);
        CSchreibeVL(stdout,"Essential hyperplanes",EssHP,EDim);
        HilbBasHsp(EssHP,EDim,&HilbBasVA, &Rang,&C0Bas,&Rang0); */

        HilbBasTri(VAC, EDim, &HilbBasVA, &DummySupp, &Rang, &Index);
        FromList=HilbBasVA.St;

        printf("#Hilb bas %d\n",NinVL(HilbBasVA));
        CSchreibeVL(stdout,
                  "Hilbert basis of resttricted nef cone in PL, ext by int",HilbBasVA,EDim);
        RunTime();
        fflush(stdout);
    }
    else
        FromList=VAC.St;



    StartVL(&InclMin);
    for(VP=FromList;VP;VP=VP->R)  // Betrachte nur inklusionsminimale Polytope
        if(VP->V[EDim-1]==1)
            RedInsertPosOrth(&InclMin,EDim,VP->V);

    FreeVL(&HilbBasVA);
    FreeVL(&VAC);


    printf("Min support pol given by\n");
    FSchreibeVL(stdout,InclMin,EDim);
    RunTime();
    fflush(stdout);

    if(!ChPoly)
        exit(0);

    if(!InclMin.St)
    {
        printf("No integral support polytopes from extreme rays\n");
        exit(0);
    }

    Aus=fopen(Uebergabe,"w");

    for(VP=InclMin.St;VP;VP=VP->R)
    {
        // printf("Durchlauf %lld\n", VP->V[EDim-1]);
        // SchreibeV(stdout,VP->V,EDim);
        if(VP->V[EDim-1]==1) // Jetzt ueberfluessige Abfrage !!!
        {
            j=0;            // Mache Stuetzhyperebenen fuer Support-Polytop
            i=0;            // Koennen dafuer OneSkel benutzen und letzte Komponente
                            // setzen
            for(VP1=OneSkel.St;VP1;VP1=VP1->R)
            {
                if(InVertices(j,FF,Dim))  // Setze unterdrueckte Komponenten wieder ein
                    VP1->V[Dim]=0;
                else
                {
                    VP1->V[Dim]=VP->V[i];
                    i++;
                    // printf("Drin\n");
                }
                j++;
            }

            FSchreibeVL(Aus,OneSkel,Dim+1);
            fflush(stdout);

        }
    }

    fclose(Aus);

    if(DiagTest)
        run_pgm("ChPoly",Uebergabe," ","");
    else
        run_pgm("ChPoly",Uebergabe," ","-n");

    printf("Back in Support from ChPoly\n");

    FileExistsRm(Uebergabe,"");

}
