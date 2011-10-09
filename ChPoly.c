// ChPoly.c, Version 1.2
// (C) Winfried Bruns, 2011
// Überprüft glatte Polytope auf vermutete Eigenschaften

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


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>



#define LONGLONG

#include "ToricExp.c"

main(int argc, char **argv)
{

    int Dim,i,j,kk,Rang,NSupp,MaxDegCheck,MaxHBVol,MaxTotVol,q,qmax,
                    MaxDimDiagSplit,MaxMult,MaxHBConnCheck,NVert,NPol;
    Zahl Nenner,Index,Volu,dfact,SuppVolu;
    VList PolytV, HilbBasPolyt, PolytSupp,
             CutPolytV, CutPolytSupp, CheckPolytV,CheckPolytSupp;
    SimplexList Tri;
    VLE *VP, *VP1, *VPNext;
    FILE *Ein, *Prot, *Aus;
    boolean Normal, SuperConn,DiagSplitTest, SendMail,CheckCutPoly,
           CheckSuppPoly,CurrentFanDatabase,WritePolyToBase,Vertices,NoIsolatedPoint,IsDiagSplit;

    char database[100], database1[100], Address[100],Loeschen[100],
                Inputfile[100],Params[20];
    int HVNmz[4];

#define NUM_OPTIONS 13

    char *DefOptions[NUM_OPTIONS]=
    {
        "Timer true",
        "NameDepotPoly ../VertSmPol", // (empty) to disable
        "Address wbruns@uos.de", // (empty) to disable e-mail
        "CheckSuppPoly true",
        "CheckCutPoly true",
        "MaxModDiagSplit 5", // 0 to switch off
        "MaxDimDiagSplit 7",
        "AbsoluteMaxMult 100000000",
        "MaxVolCheck 100000",
        "MaxVolNormCheck 10000", // 0 to switch off
        "MaxHBDeg2Check 2000", // 0 to switch off
        "MaxHBConnCheck 10000", // 0 to switch off
        "end_options"
    };

    char **Options;

    GetAndPrintPID(argv[0],database);
    strcat(database,"ChPoly");

    InitFreeLists();
    RunTime();
    TransferDefOptions(&Options,DefOptions,NUM_OPTIONS);
    ReadConfigFile("./ChPoly",Options);


    IntOption(&MaxMult,"AbsoluteMaxMult",Options); // Absolute Obergrenze
                    // für Multiplizität des Polytops. Falls größer,
                    // nicht einmal CutFaces
    IntOption(&MaxTotVol,"MaxVolCheck",Options); // Obergrenze fuer alle Checks
    IntOption(&MaxHBVol,"MaxVolNormCheck",Options); // Nur bis zu diesem eukl. Vol.
                                                    // wird auf Norm. getestet
    IntOption(&MaxDegCheck,"MaxHBDeg2Check",Options); // Obergrenze für Deg2Check
    IntOption(&MaxHBConnCheck,"MaxHBConnCheck",Options); // Sool Zshggstest gemacht werden ?

    IntOption(&qmax,"MaxModDiagSplit",Options); // Max Modul für DiagSplit
    IntOption(&MaxDimDiagSplit,"MaxDimDiagSplit",Options); // Max Dim für DiagSplit

    // BoolOption(&VERBOSE,"Verbose",Options);
    BoolOption(&TIMER,"Timer",Options);
    BoolOption(&CheckSuppPoly,"CheckSuppPoly",Options);
    BoolOption(&CheckCutPoly,"CheckCutPoly",Options);
    // BoolOption(&CurrentFanDatabase,"CurrentFanDatabase",Options);
    CurrentFanDatabase=true;


    StringOption(Address,"Address",Options);
    if(strcmp(Address,"")==0)
        SendMail=false;
    else
        SendMail=true;
    StringOption(database1,"NameDepotPoly",Options);
    if(strcmp(database1,"")==0)
        WritePolyToBase=false;
    else
        WritePolyToBase=true;

    Vertices=false;
    DiagSplitTest=true;

    strcpy(Params,"");

    for(i=1;i<argc;i++)
    {
        if(strncmp(argv[i],"-",1))
            strcpy(Inputfile,argv[i]);
        else
            strcat(Params,argv[i]);
    }
    if(strchr(Params,'v'))
        Vertices=true;
    if(strchr(Params,'n'))
        DiagSplitTest=false;

    Dim=GetDim(Inputfile); // Hier Kegel-Dimension
    MaxDim=Dim;

    dfact=1;
    for(i=2;i<Dim;i++)
        dfact*=i;

    // printf("Dim = %d\n",Dim);

    Ein=fopen(Inputfile,"r");

NPol=0;
while(1)
{
    if(Vertices)
    {
        if(!ReadNextVL(Ein,&PolytV,&NVert,Dim))
            break;
        SuppHypVL(&PolytSupp,PolytV,Dim);  // Stützhyperebenen berechnen
    }
    else
    {
        if(!ReadNextVL(Ein,&PolytSupp,&NSupp,Dim))
            break;
        SuppHypVL(&PolytV,PolytSupp,Dim);  // Ecken berechnen
    }

    NPol++;
    printf("******** Investigating polytope %d\n",NPol);


    RunTime();
    printf("#vert supp pol %d\n",NinVL(PolytV));
    fflush(stdout);
    FreeVL(&PolytSupp); // wird neu berechnet

    Volu=VolPolyt(PolytV,Dim,&PolytSupp,&Tri);
    SuppVolu=Volu;
    FreeSimplexL(&Tri);
    printf("Mult supp pol %lld\n",Volu);
    RunTime();
    fflush(stdout);

    if(Volu>MaxMult)
    {
        printf("Support polytope too large for any action\n");
        break;
    }

for(kk=0;kk<=1;kk++)
{
    if(kk==0)
    {
        if(!CheckSuppPoly)
            continue;
        CheckPolytV=PolytV;
        CheckPolytSupp=PolytSupp;
        printf("Checking support polytope\n");
    }


    if(kk==1)
    {
        printf("Computing cut polytope\n");
        CutFaces(&CutPolytV, &CutPolytSupp, PolytV, PolytSupp,Dim);
        RunTime();
        fflush(stdout);
        FreeVL(&CutPolytSupp); // wird neu berechnet
        printf("Checking cut polytope\n");
        Volu=VolPolyt(CutPolytV,Dim,&CutPolytSupp,&Tri); // Volu = Multiplizität
        if(Volu==SuppVolu)
        {
            printf("Cut poly = supp poly\n");
            if(CheckSuppPoly)
                break;
        }
        FreeVL(&PolytV); // nicht mehr benötigt
        FreeVL(&PolytSupp);
        FreeSimplexL(&Tri);
        CheckPolytV=CutPolytV;
        CheckPolytSupp=CutPolytSupp;
        printf("Mult cut pol %lld\n",Volu);
        RunTime();
        fflush(stdout);
    }

    if(Volu>dfact*MaxTotVol)
    {
        printf("Polytope too large\n");
        break;
    }

    if(CurrentFanDatabase && FindInBase(CheckPolytV, CheckPolytSupp,Dim, database,true,true))
    {                          // Erst den laufenden Faecher ueberpruefen
            continue;
    }

    printf("#Vert check pol %d\n",NinVL(CheckPolytV));
    FSchreibeVL(stdout,CheckPolytV,Dim);
    fflush(stdout);

    if(WritePolyToBase)
        HaengeVLan(database1,CheckPolytV,Dim,"Vertices");
        
  IsDiagSplit=false;

  if(DiagSplitTest && MaxDimDiagSplit>=Dim-1)
   {
    for(q=3;q<=qmax;q++)
    {
        if(q%2==0 || (q>3 && q%3==0))
            continue;
        // printf("Test DiagSplit for %d\n",q);
        if(DiagSplit(CheckPolytSupp,Dim-1,q))
            {
                printf("In ChPoly diagonally split for q = %d\n",q);
                IsDiagSplit=true;
                break;
            }
    }
    if(q>qmax && qmax>0)
        printf("In ChPoly not diagonally split for 3 <= q <= %d\n",qmax);
  }
  else
    printf("No diag split test\n");

    FreeVL(&CheckPolytSupp); // wird neu berechnet
    

    if(Volu<MaxHBVol*dfact)
    {
        printf("Making Hilbert basis\n");
        fflush(stdout);

        HilbBasTri(CheckPolytV,Dim, &HilbBasPolyt,
            &CheckPolytSupp, &Rang, &Index);

        printf("#Hilbert basis of cone over polytope %d\n",
                     NinVL(HilbBasPolyt));
        RunTime();
        fflush(stdout);

        Normal=true;
        for(VP1=HilbBasPolyt.St;VP1;VP1=VPNext)
        {
            VPNext=VP1->R;
            if(VP1->V[Dim-1]!=1)
            {
                Normal=false;
                // FreeVLE(&HilbBasPolyt,VP1);
            }
        }
        NmzHV(HVNmz, HilbBasPolyt,Dim, Inputfile); // get h-vector, check Ehrhart
        if(!Normal)
        {
            printf("Counterexample: not Normal\n");
            if(SendMail)
            {
                TextMail(Address,"Nonnormal counterexample");
                VLMail(Address,CheckPolytV,Dim);
            }

        }
    }
    else
    {
        printf("Making only lattice points\n");
        fflush(stdout);
        FindLatPoints(CheckPolytV, Dim, &HilbBasPolyt,&CheckPolytSupp);
        printf("#Lattice points of polytope %d\n",
                     NinVL(HilbBasPolyt));
        RunTime();
        fflush(stdout);

    }

    if(NinVL(HilbBasPolyt)<MaxDegCheck && !IsDiagSplit)
    {
        if(!Deg2GenDeg3(HilbBasPolyt,Dim,HVNmz))
        {
            printf("Counterexample Deg 2\n");
            FSchreibeVL(stdout,HilbBasPolyt,Dim);
            if(SendMail)
            {
                TextMail(Address,"Degree > 2");
                VLMail(Address,CheckPolytV,Dim);
            }
        }
        printf("After deg2/deg3 check\n");
    }
    else
    {
        if(NinVL(HilbBasPolyt)<=MaxHBConnCheck && !IsDiagSplit)
        {
            if(!SyzBors(HilbBasPolyt, CheckPolytSupp,
                                 Dim, Dim*NinVL(HilbBasPolyt), 100))
            {
                if(SendMail)
                {
                    TextMail(Address,"Disconnected counterexample");
                    VLMail(Address,CheckPolytV,Dim);
                }
                printf("Disconnected counterexample\n");
            }
        printf("After deg2 random conn check\n");
        }
    }

    RunTime();
    fflush(stdout);


    if(!NewDeg2AffCover(CheckPolytV,HilbBasPolyt,CheckPolytSupp, Dim,&SuperConn,&NoIsolatedPoint))
    {
        printf("Counterexample: affine cover not of not degree 2\n");
        TextMail(Address,"Non deg 2 affine cover counterexample");
        VLMail(Address,CheckPolytV,Dim);
    }
    if(!NoIsolatedPoint)
    {
        if(SendMail)
        {
            TextMail(Address,"Example with isolated point");
            VLMail(Address,CheckPolytV,Dim);
        } 
        printf("Example with isolated point\n");   
    }
    
    if(!SuperConn)
        printf("Not superconnected\n");
        
    printf("After deg 2 affine cover\n");
    RunTime();
    fflush(stdout);
    
    if(!StronglyConnected(CheckPolytV, HilbBasPolyt, CheckPolytSupp, Dim))
    {
        if(SendMail)
        {
            TextMail(Address,"Counterexample to strong connectivity");
            VLMail(Address,CheckPolytV,Dim);
        }
        printf("Not strongly connected\n");
    }

    printf("After check for strong conn\n");
    RunTime();
    fflush(stdout);
    
    if(NinVL(HilbBasPolyt)<=MaxHBConnCheck)
    {
        if(!AbundantDeg2Rel(CheckPolytV,HilbBasPolyt,CheckPolytSupp,Dim))
        {
            TextMail(Address,"Empty Symm Pol");
            VLMail(Address,CheckPolytV,Dim);
        }
        printf("After abundance check\n");
        RunTime();
        fflush(stdout);
        
    }

    FreeVL(&HilbBasPolyt);

    printf("------------- Done with polytope\n");

    if(!CheckCutPoly)
        break;

} // kk

} // Ende ewige Schleife

    sprintf(Loeschen,"%s_%d.base",database,Dim);
    FileExistsRm(Loeschen,"");

}
