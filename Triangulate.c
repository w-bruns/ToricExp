// Triangulate.c, Version 1.2
// (C) Winfried Bruns, 2011
// Berechnet unimdodulare Triangulierungen

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


#define LONGLONG

#include "ToricExp.c"

void main(int argc, char **argv)
{

    FILE *Ein, *Aus;
    int i,N,Dim,Rank;
    Zahl Index;
    VList Vert, OneSkel,Supp,NewOneSkel;
    SimplexList Tri,NewTriang;
    VLE *VP;
    boolean ValidPar;
    char FileName[1000];

    PrintPID(argv[0]);

    ValidPar=false;
    if(argc<3)
    {
        printf("Invalid parameters\n");
        exit(1);
    }

    if(!strncmp(argv[2],"C",1))
    {
        ValidPar=true;
        Dim=GetDim(argv[1]);
        MaxDim=Dim;
        NLiesVL(argv[1],&OneSkel,&N,&Dim);
        RangIndexVL(OneSkel,Dim,&Rank,&Index);
        if(Rank<Dim)
        {
            printf("Cone not of maximal dimension\n");
            exit(1);
        }
        TriangSuppVL(&Tri, &Supp, OneSkel,Dim);
        SetCornSupp(OneSkel, Tri,Dim);
        UnimodTriang(&NewOneSkel,&NewTriang,OneSkel,Tri,Dim);
    }

    if(!strncmp(argv[2],"P",1))
    {

        ValidPar=true;
        Dim=GetDim(argv[1]);
        MaxDim=Dim;
        NLiesVL(argv[1],&Vert,&i,&N);
        if(!IsFlat(Vert,Dim))
        {
            printf("Polytopal data non-flat\n");
            exit(1);
        }
        RangIndexVL(Vert,Dim,&Rank,&Index);
        if(Rank<Dim)
        {
            printf("Polytope not of maximal dimension\n");
            exit(1);
        }
        SuppHypVL(&Supp,Vert,Dim);
        ResolveNormalFan(&NewOneSkel,&NewTriang, Vert,Supp,Dim);
        Dim--; // for output
    }

    if(!strncmp(argv[2],"F",1))
    {
        ValidPar=true;
        Dim=GetDim(argv[1]);
        MaxDim=Dim;
        Ein=fopen(argv[1],"r");
        FLiesTriang(Ein,&OneSkel,&N,&Dim,&Tri);
        RangIndexVL(OneSkel,Dim,&Rank,&Index);
        if(Rank<Dim)
        {
            printf("Fan not of maximal dimension\n");
            exit(1);
        }
        SetCornSupp(OneSkel,Tri,Dim);
        UnimodTriang(&NewOneSkel,&NewTriang,OneSkel,Tri,Dim);
    }

    if(!ValidPar)
    {
        printf("Invalid parameters\n");
        exit(1);
    }

    strcpy(FileName,argv[1]);
    strcat(FileName,".ufn");

    Aus=fopen(FileName,"w");
    SchreibeTriang(Aus,NewOneSkel,Dim,NewTriang);

}
