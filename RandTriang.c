// RandTriang.c, Version 1.2
// (C) Winfried Bruns, 2011
// Berechnet "zufällige" unimodulare projektive Fächer

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

#define LONGLONG

#include "ToricExp.h"

void CompleteRandTriang(VList *OneSkel, SimplexList *Triang, int Dim, Zahl Bound, int ExtraV)
{
    Matrix M;
    int i,j,NVect,EDim;
    Simplex *Smp, *Smp1;
    Indizes Extr;
    VList Dummy;
    VLE *VP, *VPR;
    SimplexList PrelTriang;

    // printf("Start Compl VZ %d FVZ %d\n",VektorZaehler, NinVL(VLFree));

    NVect=Dim+1+rand() % (ExtraV+1);
    EDim=Dim+1;

    MatAlloc(&M,NVect+1,Dim+1);
    Extr=calloc(NVect+1,sizeof(int));
    StartVL(OneSkel);
    StartSimplexL(Triang);

    MacheEinhMat(M,Dim);
    for(i=0;i<Dim;i++)
        M[Dim][i]=-(1+rand() % Bound);  // Damit haben wir die vertikale Achse in Inneren
    MakePrime(M[Dim],Dim);

    // SchreibeVL(stdout,*OneSkel,Dim);

    for(i=Dim+1;i<NVect;i++)
    {
        NZRandVect(M[i],Dim,Bound);
        MakePrime(M[i],Dim);
    }

    for(i=0;i<NVect;i++)
        M[i][Dim]=1+rand()%(2*Bound+1); // letzte Komponente positiv

    for(j=0;j<Dim;j++)
        M[NVect][j]=0;
    M[NVect][Dim]=-1; // Wir sehen uns das von unten an

    Mat2VL(OneSkel,M,NVect+1,Dim+1);

    /* printf("--------------------\n");
    FSchreibeVL(stdout,*OneSkel,Dim+1);
    printf("--------------------\n");*/

    TriangSupp(&PrelTriang, &Dummy, M,NVect+1, EDim, Extr);

    // SchreibeTriang(stdout,*OneSkel,Dim+1,PrelTriang);
    FreeVL(&Dummy);

    for(Smp=PrelTriang.St;Smp;Smp=Smp->R)
    {
        for(i=0;i<Dim+1;i++)
            if(Smp->Vert[i]==NVect)
                break;
        if(i==Dim+1)  // Kein "unteres" Simplex
            continue;
        Smp1=NewSimplexLE(Triang);
        j=-1;
        for(i=0;i<Dim+1;i++) // Ecken übernehmen bis auf untere
        {
            if(Smp->Vert[i]==NVect)
                continue;
            j++;
            Smp1->Vert[j]=Smp->Vert[i];
        }
    }
    FreeSimplexL(&PrelTriang);

    for(Smp=Triang->St;Smp;Smp=Smp->R) // Wirklich verwendete Vektoren registrieren
        for(i=0;i<Dim;i++)
            Extr[Smp->Vert[i]]=1;
    j=1;
    for(i=0;i<NinVL(*OneSkel);i++) // Durchzählen der verwendeten Vektoren
        if(Extr[i])                // Index um 1 zu groß
        {
            Extr[i]=j;
            j++;
        }

    /* printf("==========================\n");
    for(i=0;i<NinVL(*OneSkel);i++)
        printf("%d ",Extr[i]);
    printf("==========================\n"); */

    for(Smp=Triang->St;Smp;Smp=Smp->R) // Transformation der Indizes
        for(i=0;i<Dim;i++)
            Smp->Vert[i]=Extr[Smp->Vert[i]]-1; // Index wieder verkleinert

    j=-1;
    for(VP=OneSkel->St;VP;VP=VPR)  // Entfernen überflüssiger Vektoren
    {
        j++;
        VPR=VP->R;
        if(Extr[j])
            continue;
        FreeVLE(OneSkel,VP);
    }

    /* printf("==========================\n");
    printf("Start 0\n");
    SchreibeTriang(stdout,*OneSkel,Dim,*Triang);
    printf("==========================\n"); */

    // printf("Vor Set VZ %d FVZ %d OneSkel %d\n",VektorZaehler, NinVL(VLFree),NinVL(*OneSkel));


    SetCornSupp(*OneSkel,*Triang,Dim);

    // printf("Nach Set VZ %d FVZ %d OneSkel %d CornSupp %d\n",VektorZaehler, NinVL(VLFree),NinVL(*OneSkel)
    //                          , 2*Dim*NinSimplL(*Triang));


    FreeMat(M,NVect+1);
    free(Extr);
}

void main(int argc, char **argv)
{

    int Dim, i, Last,ExtraV,MaxOneSkel,MaxTriang,NAtt;
    int Bound;
    SimplexList Triang, NewTriang;
    VList OneSkel, NewOneSkel;
    Simplex *Smp;
    FILE *Aus, *Aus1, *Ein;
    clock_t ticks;
    boolean WriteFanToDepot,SearchFanInDatabase,Verbose;

    char database[100], FanDepot[100], DummyStr[100];

    char Uebergabe[100];

#define NUM_OPTIONS 8

    char *DefOptions[NUM_OPTIONS]=
    {
        "Timer false",
        "NameDatabaseFan ../Fan_Smooth", // (empty) to disable
        "NameDepotFan (empty)", // (empty) to disable
        "MaxPicRank 20", //
        "MaxConesInFan 150",
        "Verbose true",
        "seed (empty)",
        "end_options"
    };

    char **Options;

    GetAndPrintPID(argv[0],Uebergabe);

    if(argc<4)
    {
        printf("RandTriang: not enough parameters\n");
        exit(1);
    }

    InitFreeLists();
    RunTime();
    TransferDefOptions(&Options,DefOptions,NUM_OPTIONS);
    ReadConfigFile("./RandTriang",Options);
    BoolOption(&Verbose,"Verbose",Options);
    BoolOption(&TIMER,"Timer",Options);

    StringOption(DummyStr,"seed",Options);
    if(strcmp(DummyStr,"")==0)
        srand(time(&ticks));
    else
    {
        sscanf(DummyStr,"%d",&ticks);
        srand(ticks);
    }


    FixedLVSize=1000;

    sscanf(argv[1],"%d",&Dim); // Hier Polytop-Dimension
    MaxDim=Dim+1;  // +1 wichtig in CompleteRandTriang

    StringOption(FanDepot,"NameDepotFan",Options);
    if(strcmp(FanDepot,"")==0)
        WriteFanToDepot=false;
    else
        WriteFanToDepot=true;

    StringOption(database,"NameDatabaseFan",Options);
    if(strcmp(database,"")==0)
        SearchFanInDatabase=false;
    else
        SearchFanInDatabase=true;

    sscanf(argv[2],"%d",&Bound);

    sscanf(argv[3],"%d",&ExtraV);

    IntOption(&MaxOneSkel,"MaxPicRank",Options);
    IntOption(&MaxTriang,"MaxConesInFan",Options);

    strcat(Uebergabe,".ufn");

    printf("Start RandTriang ticks %d Dim %d Bound %d ExtraV %d MaxPicRank %d MaxConesInfan %d\n",
                  ticks,Dim,Bound,ExtraV,MaxOneSkel,MaxTriang);

    fflush(stdout);

      NAtt=0;

    RunTime();

while(!FileExistsRm("RandTriang.stop","")) // eternal loop
{

   NAtt++; // zaehlt Versuche bis zum naechsten Erfolg

   if(NAtt>10000)
   {
       printf("No new fans found in 10000 attempts\n");
       exit(0);
    }

   CompleteRandTriang(&OneSkel, &Triang,Dim,Bound,ExtraV);
   RunTime();


   printf("Random triang made\n");fflush(stdout);


    /* printf("==========================\n");
    printf("Start\n");
    SchreibeTriang(stdout,OneSkel,Dim,Triang); */

    UnimodTriang(&NewOneSkel,&NewTriang,OneSkel,Triang,Dim);

     printf("UniModular refinement made\n");
     RunTime();
     fflush(stdout);


    FreeVL(&OneSkel);
    FreeSimplexL(&Triang);

    if(NinVL(NewOneSkel)-Dim>MaxOneSkel || NinSimplL(NewTriang)>MaxTriang)
    {
        printf("Fan too large: %d > %d or %d > %d\n",
            NinVL(NewOneSkel)-Dim,MaxOneSkel,NinSimplL(NewTriang),MaxTriang);
        goto Freimachen;
    }

    if(!SearchFanInDatabase || !FindInBaseFan(NewOneSkel,NewTriang,Dim,
                  database, true, false))
    {

        printf("+++++++++++++++++++++++++++++++++++++++++\n");
        printf("Investigating new unimodular fan\n");
        RunTime();
        if(Verbose)
        {
            printf("Unimdoular fan\n");
            SchreibeTriang(stdout,NewOneSkel,Dim,NewTriang);
        }

        NAtt=0;  // Versuchszaehler zuruecksetzen

        if(WriteFanToDepot)
            HaengeTriangAn(FanDepot,NewOneSkel,Dim,NewTriang,"Fan");

        Aus1=fopen(Uebergabe,"w");
        SchreibeTriang(Aus1,NewOneSkel,Dim,NewTriang);
        fclose(Aus1);

        run_pgm("Support",Uebergabe,"","");

        printf("Investigation of fan finished\n");
        remove(Uebergabe);
        printf("+++++++++++++++++++++++++++++++++++++++++\n");
    }

    Freimachen:

        FreeVL(&NewOneSkel);

        FreeSimplexL(&NewTriang);

        if(VektorZaehler!=NinVL(VLFree))
        {
            printf("VZ %d FVZ %d\n",VektorZaehler, NinVL(VLFree));
            exit(0);
        }
}

    printf("RandTriang stopped\n");

}
