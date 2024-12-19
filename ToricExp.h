// ToricExp.h, Version 1.2
// Bibliothek mit Algorithmen für Kegel und affine monoide
// (C) Winfried Bruns 2011

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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Daten-Typen

#ifdef LONGLONG

typedef long long Zahl;

#else

typedef int Zahl;

#endif

typedef Zahl *Vektor;
typedef Vektor *Matrix;
typedef int *Indizes;  // fuer Index-Listen in Triangulierungen etc.

typedef int boolean;

const boolean true=1, false=0;

// Listen fuer "lange" Vektoren

typedef struct LongV
{
  Vektor V;  // Laenge beliebig, gegeben durch Size
  int Size;
  Zahl Mark; // Hier kann man sich was merken
} LongV;

typedef struct LongVLE // Element fuer Liste von Langvektoren
{
  LongV *LV;
  struct LongVLE *L, *R;
} LongVLE;

typedef struct LongVList // die Liste
{
  LongVLE *St, *E;
} LongVList;

// Vektorlisten

typedef struct VLE // analog
{
  Vektor V;   // Der eigentliche Vektor, wird mit Laende MaxDim allokiert
  Vektor LV;  // wird nicht automatisch allokiert, kopiert etc.
              // fuer Zuordnung von werten unter Support-Formen etc.
  LongVLE *LVP; // bei FixedLVSize>0 werden die LVs in einem Pool gehalten
                // vom typ LongVList. LVP verweist auf dieses Listenelement
                // und LV auf den entsprechenden Vektor
  Zahl Mark,Mark1,Mark2,Mark3;  // Hier kann man sich etwas merken, nicht automatisch kopiert
  struct VLE *L, *R;
} VLE;

typedef struct VList
{
  VLE *St, *E;
} VList;

// Simplex-Listen

typedef struct Simplex
{
  Indizes Vert;
  VList Corn, Supp;
  Zahl Vol;
  struct Simplex *L, *R;
} Simplex;

typedef struct SimplexList
{
  Simplex *St, *E;
} SimplexList;

// Globale Variable

int VektorZaehler=0, VektorZaehler1=0, MaxDim=-1; // VZs fuer Statistik, siehe NewVLE
            //
            // MaxDim bestimmt die Laenge der in VLEs allokierten Vektoren
            //
            // Achtung: MUSS UNBEDINGT GESETZT WERDEN
            //
VList VLFree; // Fuer das Recycling allokierter VLEs
LongVList LongVLFree; // dito LVLEs
LongVList LVPool; // Liste fuer in VLists gebundene LVs
SimplexList SimplexFree; // dito

boolean CheckLists=0; // Debugging-Hilfe, siehe FreeVL

boolean VERBOSE=0; // Damit kann CSchreibeVL aktiviert werden
boolean SENDMAIL=1; // Damit wird das Versenden von E-Mails aktiviert
boolean TIMER=0;    // Damit kann Timer aktiviert werden

int FixedLVSize=0, AnzLV=0, AnzLV1=0, MaxLVSize=0;
int AnzSpx=0, AnzSpx1=0;

// int TFV; // fuer Debugging von TransformValues benutzt

void InitFreeLists()
{
  VLFree.St=NULL;
  VLFree.E=NULL;
  LongVLFree.St=NULL;
  LongVLFree.E=NULL;
  LVPool.St=NULL;
  LVPool.E=NULL;
  SimplexFree.St=NULL;
  SimplexFree.E=NULL;
}


// Ein- und Ausgabe von Zahlen

#ifdef LONGLONG

void Schreibe(FILE *aus, Zahl Z)
{
    fprintf(aus,"%lld ",Z);
}

void Lies(FILE *ein, Zahl *Z)
{
    fscanf(ein,"%lld",Z);
}

#else

void Schreibe(FILE *aus, Zahl Z)
{
    fprintf(aus,"%d ",Z);
}

void Lies(FILE *ein, Zahl *Z)
{
    fscanf(ein,"%d",Z);
}

#endif

// Arithmetik fuer Zahlen

int intvergl(const void *m,const void *n)
{
    return(*(int *)n-*(int *)m);
}

int intvergl_asc(const void *m,const void *n)
{
    return(*(int *)m-*(int *)n);
}


int Zahlvergl(const void *m,const void *n)
{
    return(*(Zahl *)n-*(Zahl *)m);
}


inline void Vertausche(Zahl *u, Zahl *v)
{
  register Zahl Dummy;

  Dummy = *u;
  *u = *v;
  *v = Dummy;
}

int sign(Zahl u)
{
  int Result;

  Result = 0;
  if (u < 0)
    return -1;
  if (u != 0)
    return 1;
  return Result;
}

Zahl Zabs(Zahl u)
{
    if(u>=0) return(u);
    return(-u);
}


Zahl ggt(Zahl z1, Zahl z2)
{
  register Zahl a, b, r;

  if (z1 == 0)
    return Zabs(z2);
  if (z2 == 0)
    return Zabs(z1);

  a = Zabs(z1);
  r = Zabs(z2);
  do {
    b = r;
    r = a % b;
    a = b;
  } while (r != 0);
  return b;

}


Zahl kgv(Zahl u, Zahl v)
{
  if (u == 0 || v == 0)
    return 0;
  else
    return Zabs(u * v / ggt(u, v));
}

// Funktionen fuer Vektoren

void AddVekt(Vektor u, Vektor v, Vektor w, int Dim)
{
    register int i;
    for(i=0;i<Dim;i++)
        u[i]=v[i]+w[i];
}

void SubtVekt(Vektor u, Vektor v, Vektor w, int Dim)
{
    register int i;
    for(i=0;i<Dim;i++)
        u[i]=v[i]-w[i];
}

boolean NonNegDiff(Vektor u, Vektor v, Vektor w, int Dim)
{
    register int i;
    for(i=0;i<Dim;i++)
    {
        u[i]=v[i]-w[i];
        if(u[i]<0)
            return false;
    }
    return true;
}

void MacheNullVekt(Vektor u, int Dim)
{
    int i;
    for(i=0;i<Dim;i++)
        u[i]=0;
}

boolean NullVekt(Vektor V, int Dim)
// testet auf NullVektor
{
    int i;
    for(i=0;i<Dim;i++)
        if(V[i])
            return false;
    return true;
}

inline Zahl SkalProd(Vektor V1, Vektor V2, int Dim)
{
  register int i;
  register Zahl S;

  S = 0;
  for (i = 0; i < Dim; i++)
    S += V1[i] * V2[i];
  return S;
}

void SxV(Vektor Res, int Dim, Zahl F, Vektor V)
// Res = Skalar F mal Vektor V
{
    int i;
    for(i=0;i<Dim;i++)
        Res[i]=F*V[i];
}

Zahl MakePrime(Vektor V, int Dim)
// Macht die Komponenten eines Vektors teilerfremd
// Gibt ggT zurueck
{
  int i;
  Zahl G;

  G = V[0];
  if (G == 1 || G == -1)
    return 1;
  i = 1;
  while (i <Dim)
  {
    G = ggt(G, V[i]);
    if (G == 1)
      return 1;
    i++;
  }
  for (i = 0; i < Dim; i++)
    V[i] /= G;
  return G;
}

void NZRandVect(Vektor V, int Dim, int Bound)
{
    boolean Zero;
    int i;

    Zero=true;
    do
    {
        for(i=0;i<Dim;i++)
        {
            V[i]=rand()%(2*Bound+1)-Bound;
            if(V[i])
                Zero=false;
        }

    }
    while(Zero);
}


void SchreibeV(FILE *aus, Vektor V, int Dim)
{
    int i;
    for(i=0;i<Dim;i++)
        Schreibe(aus,V[i]);
    fprintf(aus, "\n");
}

void LiesV(FILE *ein, Vektor V, int Dim)
{
    int i;
    for(i=0;i<Dim;i++)
        Lies(ein,&(V[i]));
}

int VCompare(Vektor V1, Vektor V2, int Dim)
// vergleicht Vektoren komponentenweise
// -1  <, 0 ==, 1 >, 2 unvergleichbar
{
    int i;
    for(i=0;i<Dim;i++)
        if(V1[i]!=V2[i]) break;
    if(i==Dim) return(0); // Gleichheit
    if(V1[i]<V2[i])
    {
        for(i++;i<Dim;i++)
            if(V1[i]>V2[i]) return(2); //unvergleichbar
        return(-1);
    }
    // jetzt nur noch V1[i]>V2[i] moeglich
    for(i++;i<Dim;i++)
        if(V1[i]<V2[i]) return(2);
    return(1);
}

// Matrizen und lineare Gleichungssysteme

void MatAlloc(Matrix *mat, int ZZ, int SZ)
// Allokiert matrix mat mit ZZ Zeilen und SZ Spalten
{
    register int i;

    *mat= (Matrix) calloc(ZZ,sizeof(Zahl*));
    for(i=0;i<ZZ;i++)
        (*mat)[i]= (Vektor) calloc(SZ,sizeof(Zahl));
}

void FreeMat(Matrix mat, int ZZ)
// Macht mat wieder frei
{
    int i;

    for(i=ZZ-1;i>=0;i--)
    {
        free(mat[i]);
    }
    free(mat);
}

void MacheEinhMat(Matrix E, int Dim)
{
    int i,j;
    for(i=0;i<Dim;i++)
    {
        for(j=0;j<Dim;j++) E[i][j]=0;
        E[i][i]=1;
    }
}

void KopiereMat(Matrix Ziel, int ZZ, int SZ, Matrix Ori)
{
    int i;
    for(i=0;i<ZZ;i++)
        memcpy(Ziel[i],Ori[i],SZ*sizeof(Zahl));
}

void TrMat(Matrix Ziel, Matrix Ori, int ZZ, int SZ)
{
    int i,j;
    for(i=0;i<ZZ;i++)
        for(j=0;j<SZ;j++)
            Ziel[j][i]=Ori[i][j];
}

void ZeilenTausch(Matrix Mat, int SZ, int P, int k)
{
  int i;
  for (i=0;i<SZ; i++)
    Vertausche(&Mat[P][i], &Mat[k][i]);
}

void ZeilenTauschP(Matrix Mat, int SZ, int P, int k)
{
  int i;
  for (i=P;i<SZ; i++)
    Vertausche(&Mat[P][i], &Mat[k][i]);
}

void SpaltenTausch(Matrix Mat, int ZZ, int P, int k)
{
  int i;
  for (i=0;i<ZZ; i++)
    Vertausche(&Mat[i][P], &Mat[i][k]);
}

void SpaltenTauschP(Matrix Mat, int ZZ, int P, int k)
{
  int i;
  for (i=P;i<ZZ; i++)
    Vertausche(&Mat[i][P], &Mat[i][k]);
}


void FormeZeilenUm(Matrix Mat, int ZZ, int SZ, int P)
{
  register int i, j;
  register Zahl F;

  for (i=P+1;i<ZZ;i++)
  {
    if (Mat[i][P] != 0)
    {
        F = Mat[i][P] / Mat[P][P];
        for (j=P; j<SZ;j++)
            Mat[i][j] -= F*Mat[P][j];
    }
  }
}


void DreieckMat(Matrix Mat, int ZZ, int SZ)
// Destruktiv
// Formt in Dreicksmatrix um mittels elementarer Zeilen- und Spaltenumformungen
{
  Zahl Test, Min;
  boolean Ausgeraeumt;
  int j, k, MinZ, mn, PS, PZ;

  mn=SZ;
  if(ZZ<mn) mn=ZZ;

  for (j=0; j<mn; j++)
  {
      for(PZ=j;PZ<ZZ;PZ++)
      {
        for(PS=j;PS<SZ;PS++)
            if(Mat[PZ][PS])break;
        if(PS<SZ) break;
      }
      if(PZ==ZZ) break;

      if(PZ!=j)
        ZeilenTausch(Mat,SZ,PZ,j);
      if(PS!=j)
        SpaltenTausch(Mat,ZZ,PS,j);

      do
      {
         Min = Zabs(Mat[j][j]);
         MinZ = j;
         Ausgeraeumt = true;
         for (k=j+1; k<ZZ; k++)
         {
           Test = Zabs(Mat[k][j]);
           Ausgeraeumt = (Ausgeraeumt && !Test);
           if (Test != 0 && (Test < Min))
           {
               Min = Test;
               MinZ = k;
           }
         }
         if (MinZ != j)
             ZeilenTausch(Mat,SZ,j,MinZ);

         FormeZeilenUm(Mat,ZZ,SZ,j);
      } while(!Ausgeraeumt);
  }
}


void RankIndexMat0(Matrix Mat, int ZZ,int SZ, int *rang, Zahl *index)
// Destruktiv
// bestimmt Rang und Index=|Det| (nach Umformung in quad Matrix)
// Falls Rang < Spaltenzahl, Index = -1
{

    int i,j,mn;

    DreieckMat(Mat,ZZ,SZ);

    mn=SZ;
    if(ZZ<mn) mn=ZZ;
    for(i=mn-1;i>=0;i--)
    {
        if(Mat[i][i]) break;
    }
    *rang=++i;

    if(*rang<SZ)
    {
        *index=-1;
        return;
    }

    *index=1;
    for(i=0;i<*rang;i++)
        *index*=Mat[i][i];
    if(*index<0) *index*=-1;
    return;
}

int RankMat0(Matrix Mat, int ZZ,int SZ)
// Destruktiv
// wie oben, jedoch nur Rang
{

    int i,j,mn, rang;

    DreieckMat(Mat,ZZ,SZ);

    mn=SZ;
    if(ZZ<mn) mn=ZZ;
    for(i=mn-1;i>=0;i--)
    {
        if(Mat[i][i]) break;
    }
    rang=++i;
    return(rang);
}

void RankIndex(Matrix Mat, int ZZ,int SZ, int *rang, Zahl *index)
// Nicht destruktiv
{
    Matrix Mat1;

    MatAlloc(&Mat1,ZZ,SZ);
    KopiereMat(Mat1,ZZ,SZ,Mat);
    RankIndexMat0(Mat1, ZZ,SZ, rang, index);
    FreeMat(Mat1,ZZ);
}

int RankMat(Matrix Mat, int ZZ,int SZ)
// Nicht destruktiv
{
    Matrix Mat1;
    int r;

    MatAlloc(&Mat1,ZZ,SZ);
    KopiereMat(Mat1,ZZ,SZ,Mat);
    r=RankMat0(Mat1, ZZ,SZ);
    FreeMat(Mat1,ZZ);
    return r;
}
void FuelleMatrix(Matrix mat, int ZZ, int Dim, Matrix Erz,  Indizes V)
// Fuellt Matrix mat durch Auswahl von Zeilen aus Erz
// auszuwaehlende Zeilen in V
{
  int i, j;
  for (i = 0; i < ZZ; i++)
  {
     for (j = 0; j < Dim; j++)
       mat[i][j] = Erz[V[i]][j];
  }
}


boolean VollRang(Matrix Erz, int Dim, Indizes C, int ZZ)
// Untersucht, ob Rang(Erz)=ZZ
{
    Matrix M;
    int r;

    MatAlloc(&M,ZZ,Dim);

    FuelleMatrix(M,ZZ,Dim,Erz,C);
    r=RankMat0(M,ZZ,Dim);
    FreeMat(M,ZZ);

    if(r==ZZ) return true;
    return false;
}

void FindeInitSimpl0(Indizes C, Matrix Erz, int AnzErz, int Dim, int Stufe)
// Bestimmt in Erz lexikograpisch kleinste Zeilenauswahl mit Rang=dim
// fuer rekursiven Aufruf
{
   if(Stufe==Dim)
      return;

   if(Stufe==0) C[Stufe]=0;
   else C[Stufe]=C[Stufe-1]+1;

   for(;C[Stufe]<AnzErz;C[Stufe]++)
   {
      if(VollRang(Erz,Dim,C,Stufe+1))
      {
          FindeInitSimpl0(C, Erz, AnzErz,Dim, Stufe+1);
          return ;
      }
   }
}


void FindeInitSimpl(Matrix SimpErz, Matrix Erz, int AnzErz, int Dim)
// Die von aussen aufzurufende Routine
// fuer InitSimpl0
{
    Indizes C;

    C=(Indizes) calloc(Dim,sizeof(int));
    FindeInitSimpl0(C,Erz,AnzErz,Dim,0);
    FuelleMatrix(SimpErz,Dim,Dim,Erz,C);
}

Zahl Vol0(Matrix Mat, int Dim)
// bestimmt Vol=|Det|
// Destruktiv
{
    int rang;
    Zahl index;

    RankIndexMat0(Mat,Dim,Dim, &rang, &index);
    if(index<0)
        index=0;
    return(index);
}

Zahl Vol(Matrix Mat, int Dim)
// Nicht destruktiv
{
    int rang;
    Zahl index;

    RankIndex(Mat,Dim,Dim, &rang, &index);
    if(index<0)
        index=0;
    return(index);
}

void SchreibeMatrix(FILE *aus, Matrix M, int ZZ, int SZ)
{
    int i;
    for(i=0;i<ZZ;i++)
        SchreibeV(aus,M[i],SZ);
}

void FSchreibeMatrix(FILE *aus, Matrix M, int ZZ, int SZ)
{

    fprintf(aus,"%d\n%d\n",ZZ,SZ);
    SchreibeMatrix(aus, M,ZZ,SZ);
}

void LiesMatrix(FILE *ein, Matrix M, int ZZ, int SZ)
{
    int i,j;
    for(i=0;i<ZZ;i++)
        for(j=0;j<SZ;j++) Lies(ein,&M[i][j]);
}

void VxM(Vektor R, Vektor F1, int Dim, Matrix F2, int SZ)
// R = Vektor x Matrix
// Dabei darf R die gleiche Adresse wie F1 haben
{
    int k,j;
    Vektor F1C;
    Zahl s;

    F1C=(Zahl*) calloc(Dim,sizeof(Zahl));
    memcpy(F1C,F1,Dim*sizeof(Zahl));

    for(j=0;j<SZ;j++)
    {
        s=0;
        for(k=0;k<Dim;k++) s+=F1C[k]*F2[k][j];
        R[j]=s;
    }
    free(F1C);
}

void MxV(Vektor R, Matrix F2, int ZZ, int Dim, Vektor F1)
// analog, nur M x V
{
    int k,j;
    Zahl s;

    Vektor F1C;

    F1C=(Zahl*) calloc(Dim,sizeof(Zahl));
    memcpy(F1C,F1,Dim*sizeof(Zahl));

    for(j=0;j<ZZ;j++)
    {
        s=0;
        for(k=0;k<Dim;k++) s+=F2[j][k]*F1C[k];
        R[j]=s;
    }
    free(F1C);
}

void MxM0(Matrix R, Matrix F1, int NZ, int NS1, Matrix F2, int NS2)
// multipliziert Matrizen: R=F1*F2
// Achtung: R muss von F1 und F2 verschieden sein
{
    int i,j,k;
    Zahl S;

    for(i=0;i<NZ;i++)
        for(j=0;j<NS2;j++)
        {
            S=0;
            for(k=0;k<NS1;k++)
                S+=F1[i][k]*F2[k][j];
            R[i][j]=S;
        }
}



void FormeZeilenUmMM(Matrix Mat, Matrix Mat1, int ZZ, int SZ, int P)
// Zeilenumformung simultan mit Mat und Mat1, P Pivotzeile
{
  register int i, j;
  register Zahl F;

  for (i=P+1;i<ZZ;i++)
  {
    if (Mat[i][P] != 0)
    {
        F = Mat[i][P] / Mat[P][P];
        for (j=0; j<SZ;j++)
        {
            Mat[i][j] -= F*Mat[P][j];
            Mat1[i][j] -= F*Mat1[P][j];
        }
    }
  }
}

void MacheLoesung( Vektor L, Matrix Mat, int Dim, Vektor RS, Zahl Nenner)
// Bestimmt Loesung zu linearem Gleichungssystem in Dreiecksform
// Mat ist quadratisch, L Loesung, RS rechte Seite
{
   register int i, j;
   register Zahl S;

   for (i=Dim-1;i>= 0;i--)
   {
        S =Nenner*RS[i];
        for (j =i+1; j<Dim; j++)
            S -= Mat[i][j] * L[j];
        L[i] = S / Mat[i][i];
   }
}

void InvertZ(Matrix L, Matrix U, int Dim, Vektor Diag, Zahl *Nenner)
// invertiert Matrix U, Inverse der Transponierten zurueckgegeben in L
// Diag ist die Diagonale nach Umformung von U auf Dreiecksform
// Nenner ist kgV ueber Diag
// Nicht destruktiv
{

  register int i, j, k,  MinZ, rank;
  boolean Ausgeraeumt;
  register Zahl Test,  Min;
  Matrix Mat, RS;
  Vektor LV, RSV;
  Zahl index;

    MatAlloc(&Mat, Dim,Dim);
    MatAlloc(&RS,Dim,Dim);
    LV= (Zahl*) calloc(Dim,sizeof(Zahl));
    RSV= (Zahl*) calloc(Dim,sizeof(Zahl));
    KopiereMat(Mat,Dim,Dim,U);
    MacheEinhMat(RS,Dim);

  for (j=0; j<Dim; j++)
  {
    do {
      Min = Zabs(Mat[j][j]);
      MinZ = j;
      Ausgeraeumt = true;
      for (k =j+1; k<Dim; k++)
      {
            Test = Zabs(Mat[k][j]);
            Ausgeraeumt = (Ausgeraeumt && Test == 0);
            if (Test != 0 && (Test < Min || Min == 0))
            {
                Min = Test;
                MinZ = k;
            }
      }
      if (!Ausgeraeumt && MinZ != j)
      {
            ZeilenTausch(Mat, Dim, j, MinZ);
            ZeilenTausch(RS, Dim, j, MinZ);
      }
      if (Mat[j][j] != 0)
            FormeZeilenUmMM(Mat, RS, Dim, Dim, j);

    } while (!(Ausgeraeumt || Mat[j][j] == 0));
  }

  *Nenner=Diag[0]=Zabs(Mat[0][0]);

  for(i=1;i<Dim;i++)
  {
    *Nenner *=Mat[i][i];
    Diag[i]=Zabs(Mat[i][i]);
  }
  *Nenner=Zabs(*Nenner);


  for(i=0;i<Dim;i++)
  {
     for(j=0;j<Dim;j++)
        RSV[j]=RS[j][i];
     MacheLoesung(LV, Mat, Dim, RSV,*Nenner);
     memcpy(L[i],LV,Dim*sizeof(Zahl));

  }


  for(i=0;i<Dim;i++)
    if(SkalProd(U[i],L[i],Dim)!=*Nenner)
        {
            printf("Arithmetic overflow InvertZ\n");
            SchreibeMatrix(stdout,Mat,Dim,Dim);
            printf("-----------------\n");
            SchreibeMatrix(stdout,RS,Dim,Dim);
            printf("-----------------\n");
            SchreibeMatrix(stdout,U,Dim,Dim);
            for(j=0;j<Dim;j++)
                Schreibe(stdout,SkalProd(U[j],L[j],Dim));
            printf("\n-----------------\n");
            printf(" Nenner = "); Schreibe(stdout,*Nenner);printf("\n");
            printf("-----------------\n");
            SchreibeMatrix(stdout,L,Dim,Dim);
            exit(1);
        }


  FreeMat(Mat,Dim);
  FreeMat(RS,Dim);
  free(LV);
  free(RSV);

}

void FormeZeilenUmMitRS(Matrix Mat, int ZZ, int SZ, Vektor RS, int P)
// Zeilenumformung inklusive rechte Seite
// P ist Pivot-zeile
{
  register int i, j;
  register Zahl F;

  for (i=P+1;i<ZZ;i++)
  {
    if (Mat[i][P] != 0)
    {
        F = Mat[i][P]/Mat[P][P];
        for (j=P; j<SZ;j++)
            Mat[i][j] -= F*Mat[P][j];
            RS[i] -= F*RS[P];
    }
  }
}


void Loese0(Vektor L, Matrix Mat, int ZZ, int Dim, Vektor RS, Zahl *Nenner, boolean *ExistUniq)
// Loest lin GLS "Mat x L =RS" mit ganzzahligen Koeff ueber den rationalen Zahlen
// L ganzzahlig mit Nenner in *Nenner
// Voraussetzung: eindeutig loesbar, Anzeige ueber *ExistsUniq
// Destruktiv
{

  int j, k, MinZ;
  Zahl Test, Min;
  boolean Ausgeraeumt;

  for (j=0;j<Dim;j++)
  {
    do
    {
      Min = Zabs(Mat[j][j]);
      MinZ = j;
      Ausgeraeumt = true;
      for (k=j+1; k<ZZ; k++)
      {
            Test = Zabs(Mat[k][j]);
            Ausgeraeumt = (Ausgeraeumt && Test == 0);
            if (Test != 0 && (Test < Min || Min == 0))
            {
                Min = Test;
                 MinZ = k;
            }
      }
      if (!Ausgeraeumt && MinZ != j)
      {
            ZeilenTausch(Mat, Dim, j, MinZ);
            Vertausche(&RS[j], &RS[MinZ]);
      }
      if (Mat[j][j] != 0)
            FormeZeilenUmMitRS(Mat, ZZ, Dim, RS, j);

    } while (!(Ausgeraeumt || Mat[j][j] == 0));
  }
  *Nenner=1;
  for(j=0;j<Dim;j++)
    if(Mat[j][j]==0)
    {
        *ExistUniq=false;
        return;
    }
    else
        *Nenner*=Mat[j][j];

  for(j=Dim;j<ZZ;j++)
    if(RS[j])
    {
        *ExistUniq=false;
        return;
    }

  *ExistUniq=true;
  *Nenner=Zabs(*Nenner);
  MacheLoesung(L, Mat, Dim, RS, *Nenner);
}

void Loese(Vektor L, Matrix Mat, int ZZ, int Dim, Vektor RS, Zahl *Nenner, boolean *ExistUniq)
// Wie Loese0, aber
// nicht destuktiv
{
    Matrix M;
    Vektor R;

    MatAlloc(&M,ZZ,Dim);
    R=(Vektor) calloc(ZZ,sizeof(Zahl));

    KopiereMat(M,ZZ,Dim,Mat);
    memcpy(R,RS,ZZ*sizeof(Zahl));

    Loese0(L,M,ZZ,Dim,R,Nenner,ExistUniq);

    FreeMat(M,ZZ);
    free(R);
}

Zahl SimpSupp( Matrix S, Matrix E, int Dim)
// Bestimmt Stuetzhyperebenen zu simp Mon erzeugt von E
// Rueckgabe in S
// Rueckgabewert Vol=|Det|
{
    Zahl Nenner, Vol, *Diag;
    int i;

    Diag= (Zahl*) calloc(Dim,sizeof(Zahl));

    InvertZ(S, E, Dim, Diag, &Nenner);

    Vol=1;
    for(i=0;i<Dim;i++)
        Vol*=Diag[i];
    if(Vol<0)
        Vol=-Vol;
    for(i=0;i<Dim;i++)
        MakePrime(S[i],Dim);
    free(Diag);
    return(Vol);
}




// Weiterhin lineare Algebra ueber Z
// jetz aber die Dinge, die auf Diagonalisieren beruhen

void ReduziereSpalte(Matrix mat, int ZZ, int SZ, int P)
// Hilfsfunktion fuer Diagonalisiere, Elementare Zeilentransformation
// Pivotzeile P
{
    int i,j,KleinstZ,StartZ;
    Zahl Min,F;

    StartZ=P+1;
    while((StartZ<ZZ) && !mat[StartZ][P]) StartZ++;
    if(StartZ==ZZ) return;

    if(mat[P][P])
    {
        KleinstZ=P;
        Min=Zabs(mat[P][P]);
    }
    else
    {
        KleinstZ=StartZ;
        Min=Zabs(mat[KleinstZ][P]);
    }
    for(i=P;i<ZZ;i++)
        if(mat[i][P] && Zabs(mat[i][P])<Min)
        {
            KleinstZ=i;
            Min=Zabs(mat[KleinstZ][P]);
        }

    if(KleinstZ!=P)
        for(i=P;i<SZ;i++)
            Vertausche(&mat[P][i],&mat[KleinstZ][i]);

    for(i=P+1;i<ZZ;i++)
    {
        F=mat[i][P]/mat[P][P];
        for(j=P;j<SZ;j++)
            mat[i][j]=mat[i][j]-F*mat[P][j];
    }
}

void ReduziereZeile(Matrix mat, int ZZ, int SZ, int P, Matrix Tr, Matrix InvTr)
// Hilfsfunktion fuer Diagonalisiere, Elementare Spaltentransformation
// Pivotzeile P
// Tr und InvTr speichern die Transformation
{
    int i,j,KleinstS,StartS;
    Zahl Min,F;

    StartS=P+1;
    while((StartS<SZ) && !mat[P][StartS]) StartS++;
    if(StartS==SZ) return;

    if(mat[P][P])
    {
        KleinstS=P;
        Min=Zabs(mat[P][P]);
    }
    else
    {
        KleinstS=StartS;
        Min=Zabs(mat[P][KleinstS]);
    }
    for(i=P;i<SZ;i++)
        if(mat[P][i] && Zabs(mat[P][i])<Min)
        {
            KleinstS=i;
            Min=Zabs(mat[P][KleinstS]);
        }

    if(KleinstS!=P)
    {
        for(i=P;i<ZZ;i++)
            Vertausche(&mat[i][P],&mat[i][KleinstS]);
        for(i=0;i<SZ;i++)
        {
            Vertausche(&Tr[i][P],&Tr[i][KleinstS]);
            Vertausche(&InvTr[P][i],&InvTr[KleinstS][i]);
        }
    }

    for(j=P+1;j<SZ;j++)
    {
        F=mat[P][j]/mat[P][P];
        for(i=P;i<ZZ;i++)
            mat[i][j]=mat[i][j]-F*mat[i][P];
        for(i=0;i<SZ;i++)
        {
            Tr[i][j]=Tr[i][j]-F*Tr[i][P];
            InvTr[P][i]=InvTr[P][i]+F*InvTr[j][i];
        }
    }
}

boolean Ausgeraeumt(Matrix mat, int ZZ, int SZ, int P)
{
    int i;


    for(i=P+1;i<ZZ;i++)
        if(mat[i][P]) break;
    if(i<ZZ)
        return(false);
    for(i=P+1;i<SZ;i++)
            if(mat[P][i]) break;
    if(i<SZ) return(false);
    return(true);
}

void Diagonalisiere0(Matrix mat, int ZZ, int SZ, Matrix Tr, Matrix InvTr, Vektor ElT,
            int *Rang, int *Index, int *Exp)
// Destruktiv
{

    int i, P, Pmax;


        MacheEinhMat(Tr,SZ);
        MacheEinhMat(InvTr,SZ);

        Pmax=ZZ;
        if(Pmax>SZ) Pmax=SZ;
        for(P=0;P<Pmax;P++)
            do
            {
                ReduziereSpalte(mat,ZZ,SZ,P);
                ReduziereZeile(mat,ZZ,SZ,P,Tr,InvTr);
            }
            while(!Ausgeraeumt(mat,ZZ,SZ,P));

        *Rang=0;
        *Index=1;
        *Exp=1;

        for(i=0; i<Pmax;i++)
        {
            if(!mat[i][i]) break;
            (*Rang)++;
            ElT[i]=Zabs(mat[i][i]);
            *Index*=ElT[i];
            *Exp=kgv(*Exp,ElT[i]);
        }
        for(;i<SZ;i++)
            ElT[i]=0;
}

void Diagonalisiere(Matrix mat, int ZZ, int SZ, Matrix Tr, Matrix InvTr, Vektor ElT,
            int *Rang, int *Index, int *exp)
// Nicht destruktiv
// Diagonalisiert mat durch Zeilen- und Spaltentransformationen. Der Basiswechsel
// in ZZ^n, n=SZ, wird in Tr und InvTr gespeichert, ElT liefert die Diagonale,
// Index das Produkt, Exp das kgV der Elemente in ElT (Laenge=Rang).
// Die neue Basis f_1,...,f_n hat die Eigenschaft, dass g_i=ElT[i]*f_i, i=1,...,r,
// eine Basis von U bilden, wobei U der von den Zeilen von mat aufgespannte Unter-
// modul ist.
// x*Tr stellt x in der neuen Basis dar, y*invTr mach dies rückgängig.
// Um die Elemente von U in der Basis g_1,...,g_r darzustellen, muss man die
// Darstellung in f_1,..,f_r komponentenweise durch die ElT[i] teilen
// Bei Rücktransformation entsprechend multiplizieren
// In den Transformationsfunktionen entspricht emb=0 der Basis f_1,..,f_r des
// kleinsten direkten Summanden, der U enthält, emb=1 der Basis g_1,..,g_r
// von U.
//
{
    Matrix mat0;

    MatAlloc(&mat0,ZZ,SZ);
    KopiereMat(mat0,ZZ,SZ,mat);
    Diagonalisiere0(mat0, ZZ, SZ, Tr, InvTr, ElT, Rang, Index, exp);
    FreeMat(mat0,ZZ);
}

void Kern(Matrix U, int ZZ, int SZ, Matrix *K, int *rangK,Matrix *Quot)
// Bestimmt Kern der linearen Abbildung ZZ^n --> ZZ^m gegeben durch
// mxn Matrix U, (m=ZZ, n=SZ) also den Durchschnitt der Kerne der durch die Zeilen
// von U gegebenen Linearformen
// Kern enthält Basis des Kerns, dessen Rang durch rangK angegeben ist
// Quot enthält eine geliftete Basis von ZZ^n/K
{
    Matrix Tr, InvTr;
    Vektor ElT;
    int i,j,Rang, Index, Exp;

    MatAlloc(&Tr,SZ,SZ);
    MatAlloc(&InvTr,SZ,SZ);
    ElT=calloc(SZ,sizeof(Zahl));

    Diagonalisiere(U,ZZ,SZ,Tr, InvTr,ElT,&Rang,&Index,&Exp); // ElT, Index, Exp werden
                                            // hier nicht gebraucht, ebenso InvTr
    *rangK=SZ-Rang;
    MatAlloc(Quot,Rang,SZ);
    MatAlloc(K,*rangK,SZ);

    for(i=0;i<Rang;i++)
        for(j=0;j<SZ;j++)
            (*Quot)[i][j]=Tr[j][i];
    for(i=Rang;i<SZ;i++)
        for(j=0;j<SZ;j++)
            (*K)[i-Rang][j]=Tr[j][i];

    for(i=0;i<ZZ;i++)
        for(j=0;j<*rangK;j++)
            if(SkalProd(U[i],(*K)[j],SZ))
            {
                printf("Arithmetic overflow Kern\n");
                SchreibeMatrix(stdout,U,ZZ,SZ);
                printf("---------------\n");
                SchreibeMatrix(stdout,*K,*rangK,SZ);
            }

    FreeMat(Tr,SZ);
    FreeMat(InvTr,SZ);
}


// Vektorlisten

void StartVL(VList *VL)
{
    VL->St=NULL;
    VL->E=NULL;
}

int NinVL(VList vectors)
// Zaehlt Elemente in Vectors
{
    VLE *VP;
    int i;

    i=0;
    for(VP=vectors.St;VP;VP=VP->R) i++;
    return(i);
}

int RNinVL(VList vectors)
// dito, nur rueckwaerts
{
    VLE *VP;
    int i;

    i=0;
    for(VP=vectors.E;VP;VP=VP->L) i++;
    return(i);
}

void MakeFirst(VList *VL,VLE *VP)
// Macht das Eelement, auf das VP zeigt, zum ersten
// Element der Liste
{
    if(VP==VL->St)
        return;

    (VP->L)->R=VP->R;
    if(VP==VL->E)
        VL->E=VP->L;
    else
        (VP->R)->L=VP->L;

    (VL->St)->L=VP;
    VP->R=VL->St;
    VL->St=VP;
    VP->L=NULL;
}


void MoveVLE(VList *Target, VList *Source, VLE *El)
// Bewegt VLE von Source nach Target
// Es wird an das Ende von Target angehaengt
//
{
  if (Source->St == El)
    Source->St = El->R;
  else
    El->L->R = El->R;
  if (Source->E == El)
    Source->E = El->L;
  else
    El->R->L = El->L;
  El->R = NULL;
  if (Target->E != NULL) {
    Target->E->R = El;
    El->L = Target->E;
  } else
    El->L = NULL;
  Target->E = El;
  if (Target->St == NULL)
    Target->St = El;
}


VLE *NewVLE(VList *VL)
// Erzeugt ein neues Element in VL am Ende
// und gibt Zeiger darauf zurueck
{
  VLE *NewV;

    VektorZaehler1++;
    if (VLFree.St == NULL)
    {
       VektorZaehler++;
       // if (!(VektorZaehler % 10000)) printf("VZ %d %d\n", VektorZaehler, VektorZaehler1);
       VLFree.St = calloc(1,sizeof(VLE));
       VLFree.St->V= (Vektor) calloc(MaxDim,sizeof(Zahl));
       VLFree.St->LVP=NULL; // Zeigt an, dass kein LV zugeordnet
       VLFree.E = VLFree.St;
       VLFree.St->L = NULL;
       VLFree.St->R = NULL;
    }
    NewV = VLFree.St;
    MoveVLE(VL, &VLFree,NewV);
    return NewV;
}

void MoveLongVLE(LongVList *Target, LongVList *Source, LongVLE *El)
// Hier, weil es gleich gebraucht wird
{
  if (Source->St == El)
    Source->St = El->R;
  else
    El->L->R = El->R;
  if (Source->E == El)
    Source->E = El->L;
  else
    El->R->L = El->L;
  El->R = NULL;
  if (Target->E != NULL) {
    Target->E->R = El;
    El->L = Target->E;
  } else
    El->L = NULL;
  Target->E = El;
  if (Target->St == NULL)
    Target->St = El;
}

void StripLVE(VLE *VP)
{
    if(VP->LVP)
        {
            if(FixedLVSize)
            {
                MoveLongVLE(&LongVLFree,&LVPool,VP->LVP);
                VP->LVP=NULL;
                VP->LV=NULL;
                return;

            }

            free(VP->LV);
            VP->LVP=NULL;
            VP->LV=NULL;
        }
}

void StripLV(VList VL)
// Entfernt die LVs aus der Liste VL
{
    VLE *VP;

    for(VP=VL.St;VP;VP=VP->R)
        StripLVE(VP);
}

void FreeVL(VList *VL)
// Macht VL frei, nachdem die LVs entfernt worden sind
{

    int fv,fr,lv,lr;

    StripLV(*VL);

    if(CheckLists) // Checken von VLFree und *VL durch Vorwaerts/Rueckwaertszaehlen
    {
        fv=NinVL(VLFree);
        fr=RNinVL(VLFree);
        lv=NinVL(*VL);;
        lr=RNinVL(*VL);
        if( fr!=fv || lr!=lv)
        {
          printf("Vor FreeVL %d %d\n",NinVL(VLFree), NinVL(*VL));
          printf("Vor FreeVL Rueckw %d %d\n",RNinVL(VLFree), RNinVL(*VL));
          exit(1);
        }
    }
    if (!VL->St) return;
    if (VLFree.St == NULL)
        VLFree = *VL;
    else
    {
        VLFree.St->L = VL->E;
        VL->E->R = VLFree.St;
        VLFree.St = VL->St;
    }
    VL->St = NULL;
    VL->E = NULL;
}

void FreeVLE(VList *VL, VLE *VP)
// macht einzelnes Element frei
{
    StripLVE(VP);
    MoveVLE(&VLFree,VL,VP);
}

void CopyConcatVL(VList *Target, VList Source1, VList Source2,int Dim)
// Bildet Target durch Verketten von S1 und S2. Letztere unveraendert
{
    VLE *S;

    Target->St = NULL;
    Target->E = NULL;
    S = Source1.St;
    while (S != NULL)
    {
        memcpy(NewVLE(Target)->V, S->V, Dim*sizeof(Zahl));
        S = S->R;
    }
    S = Source2.St;
    while (S != NULL)
    {
        memcpy(NewVLE(Target)->V, S->V, Dim*sizeof(Zahl));
        S = S->R;
    }
}

void ConcatVLFree(VList *Full, VList *Add)
// Add wird an Full angehaengt und selbst geleert
{
    if(!(Add->St)) return; // nichts hinzuzufuegen

    if(!(Full->St)) Full->St=Add->St;  // Full leer
    else  // beide nichtleer
    {
        Full->E->R=Add->St;
        Add->St->L=Full->E;
    }
    Full->E=Add->E;
    Add->St=NULL;
    Add->E=NULL;
}

void CopyVL(VList *New, VList Old,int Dim)
// LVs und Marks nucht beruecksichtigt
{
    VLE *VP,*VP1;

    StartVL(New);

    for(VP=Old.St;VP;VP=VP->R)
    {
       VP1=NewVLE(New);
       memcpy(VP1->V,VP->V,Dim*sizeof(Zahl));
    }
}

VLE* CopyVLE0(VList *Target, VLE *VP, int Dim)
// Kopiert ein einzelnes VLE
// Marks und LVs nicht beruecksichtigt
{
    VLE *VP1;

    VP1=NewVLE(Target);
    memcpy(VP1->V,VP->V,Dim*sizeof(Zahl));
    return(VP1);
}


void LiesVL(FILE *ein, VList *VL, int N, int dim)
// liest VL aus ein, wobei N Vektoren der Laenge dim gelesen werden
{
    int i,j;
    VLE *VP;

    StartVL(VL);

    for(i=0;i<N;i++)
    {
        VP=NewVLE(VL);
        for(j=0;j<dim;j++)
        {
            Lies(ein,&(VP->V[j]));
            // Schreibe(stdout,VP->V[j]);printf("\n");
        }
    }
}

int GetDim(char *FName)
// Liest Spaltenzahl der Matrix in File FName
{
    FILE *Ein;
    int Dim;
    char Dummy[1000];

    Ein=fopen(FName,"r");
    if(!Ein)
    {
        printf("File %s does not exist\n",FName);
        exit(1);
    }
    while(1)
    {
        fscanf(Ein,"%s",&Dummy); // Lesen bis Zahl !=0 gefunden ist
        if(atoi(Dummy))
            break;
    }
    fscanf(Ein,"%d",&Dim); // jetzt wirklich gelesen
    fclose(Ein);
    return(Dim);
}

boolean FLiesVL(FILE *Ein,VList *VL, int *N, int *Dim)
// Wie LiesVL, aber bestimmt Anzahl N der Vektoren und Dimension
{
    if(fscanf(Ein,"%d",N)<=0) // Fileende erreicht
        return false;
    fscanf(Ein,"%d",Dim); // Dim
    if(*Dim>MaxDim)
    {
        printf("Dim %d > MaxDim %d\n",*Dim,MaxDim);
        exit(1);
    }
    LiesVL(Ein,VL,*N,*Dim);
    return true;
}

void NLiesVL(char *FName,VList *VL, int *N, int *Dim)
// Wie FLiesVL, aber File durch Name gegeben
// File wrid hinterher geschlossen
{
    FILE *Ein;

    Ein=fopen(FName,"r");
    if(!Ein)
    {
        printf("File %s does not exist\n",FName);
        exit(1);
    }
    FLiesVL(Ein,VL,N,Dim);
    fclose(Ein);
}

boolean ReadNextVL(FILE *Ein, VList *VL,int *NVL, int Dim)
{
        int Dummy;
        char Line[100],*Test,*PureLine;

        // printf("In ReadVL\n");

        Test=fgets(Line,100,Ein);
        if(!Test)
        {
            // printf("Raus1");
            return false;
        }
        PureLine=strtok(Line," \r\t\n");
        // printf("Pureline%sPP\n",PureLine);
        if(!PureLine || !strcmp(PureLine,""))
        {
            // printf("Raus2");
            return false;
        }
        sscanf(Line,"%d",NVL);
        // printf("NVL = %d\n",*NVL);
        fscanf(Ein,"%d",&Dummy);
        // printf("Dummy = %d\n",Dummy);
        if(Dummy!=Dim)
        {
            printf("Input has wrong dimension\n");
            exit(1);
        }
        LiesVL(Ein,VL,*NVL,Dim);
        fgets(Line,100,Ein); // Auf Zeilenanfang vorgehen
        return true;
}


void SchreibeVL(FILE *aus, VList VL, int dim)
// die andere Richtung
{
    VLE *VP;

    for(VP=VL.St;VP;VP=VP->R)
        SchreibeV(aus,VP->V,dim);
}

void SchreibeVLVal(FILE *aus, VList VL, int Anz)
// schreibt die LVs der Vektoren
{
    VLE *VP;

    for(VP=VL.St;VP;VP=VP->R)
        SchreibeV(aus,VP->LV,Anz);
}

void FSchreibeVL(FILE *aus, VList VL, int Dim)
{
    fprintf(aus,"%d\n",NinVL(VL));
    fprintf(aus,"%d\n",Dim);
    SchreibeVL(aus,VL,Dim);
}

void HaengeVLan(char *database, VList VL, int Dim, char *caption)
// haengt VL and Datenbasis an.
// Wartet bis diese nicht mehr blockiert ist,
// angezeigt durch Existenz von File mit ".lock" am Ende
{
    char filename[1000], filenamelock[1000];
    int i,dummy;
    FILE *aus,*lock;

    sprintf(filename,"%s_%d.base",database,Dim);
    strcpy(filenamelock,filename);
    strcat(filenamelock,".lock");

    while(access(filenamelock, F_OK) == 0)
    {
        printf("File %s locked\n",filename);
        for(i=0;i<1000000;i++)
        {
            dummy=i;
        }
    }

    lock=fopen(filenamelock,"w");
    aus=fopen(filename,"a");
    fprintf(aus,"%s\n",caption);
    FSchreibeVL(aus,VL,Dim);
    fclose(aus);
    fclose(lock);
    remove(filenamelock);
}

void CSchreibeVL(FILE *Aus, char *Cap, VList VL, int Dim)
// Schreibt VL mit "Caption", wenn GlobalTest aktiviert
{
    if(!VERBOSE) return;
    fprintf(Aus,"%s\n",Cap);
    FSchreibeVL(Aus,VL,Dim);
    fprintf(Aus,"--------- %s\n",Cap);
}



VLE* VLEinVL(VLE *VP1, VList VL, int Dim)
// Prueft ob der Vektor von VP1 in VL vorhanden ist
{
    VLE *VP;

    for(VP=VL.St;VP;VP=VP->R)
     if(!memcmp(VP1->V,VP->V, Dim*sizeof(Zahl)))
        break;
    return(VP);
}

VLE* VinVL(Vektor V, VList VL, int Dim)
// Prueft ob der Vektor V in VL vorhanden ist
{
    VLE *VP;

    for(VP=VL.St;VP;VP=VP->R)
     if(!memcmp(V,VP->V, Dim*sizeof(Zahl)))
        break;
    return(VP);
}

boolean InsertV(VList *VL, Vektor V, int Dim)
// Fuegt V in VL ein, falls noch nicht vorhanden
{
    if(VinVL(V,*VL,Dim))
        return false;
    memcpy(NewVLE(VL)->V,V,Dim*sizeof(Zahl));
    return true;
}

VLE* VLEind(VList VL, int Choice)
// Gibt den Zeiger auf das Element mit Index Choice
// in VL zurueck. Zaehlung beginnt bei 0
{
    int k;
    VLE *VP;

    k=0;
    for(VP=VL.St;VP;VP=VP->R)
    {
        if(k==Choice)
            return(VP);
        k++;
    }
    printf("VLEind out of bounds\n");
    exit(1);
}

int IndVLEinVL(VLE *VP,VList VL)
// Sucht den Index zum VLE VP in VL
// gibt -1 zurueck, falls nicht gefunden
{
    int i;
    VLE *VP1;

    i=0;
    for(VP1=VL.St;VP1;VP1=VP1->R)
    {
        if(VP==VP1)
            return i;
    i++;
    }
    return -1;
}

void Registriere(VList Source, VLE ***Reg, int *NSource)
{
    VLE *VP;
    int i;

    *NSource=NinVL(Source);
    *Reg = (VLE**) calloc(*NSource,sizeof(VLE*));  // Vektoren registrieren
    VP=Source.St;                               // damit Autoreduktion mit Umordnung
    for(i=0;i<*NSource;i++)                       // moeglich ist
    {
        (*Reg)[i]=VP;
        VP=VP->R;
    }
}


void Diag2ResCl(VList *ResCl, Vektor Diag, int Dim, boolean MitNull)
// Bestimmt Repraesentantetnsystem ResCl von Restklassen, wenn ZZ^n/U,
// n= Dim, Filtrierung mit zyklischen Quotienten bez. Standard-Basis
// wobei die Ordnungen durch die Eintraege von Diag gegeben sind.
// (U ist unwesentlich)
// MitNull entscheidet, ob die Nullklasse zurueckgegeben wird
{
    int i,Last;
    VLE *VP;
    Vektor C;

    C=(Zahl*) calloc(Dim,sizeof(Zahl));
    StartVL(ResCl);

    if(MitNull)
    {
        VP=NewVLE(ResCl);
        memcpy(VP->V,C,Dim*sizeof(Zahl));
    }

    while(1)
    {
       Last = -1;
       for (i=0;i<Dim;i++) if (C[i]<Diag[i]-1) Last = i;
       if (Last == -1) break;
       C[Last]++;
       for(i=Last+1;i<Dim;i++) C[i]=0;
       VP=NewVLE(ResCl);
       memcpy(VP->V,C,Dim*sizeof(Zahl));
    }

    free(C);
}

void VL2Mat(Matrix mat, VList VL,int dim)
// Macht aus VList eine Matrix
// mat muss schon allokiert sein
{
    int i;
    VLE *VP;

    i=0;
    for(VP=VL.St;VP;VP=VP->R)
    {
        memcpy(mat[i],VP->V,dim*sizeof(Zahl));
        i++;
    }
}

void FuelleMatrixVL(Matrix mat, int ZZ, int Dim, VList Erz,  Indizes V)
// Fuellt Matrix mat durch Auswahl von Zeilen aus Kektorliste Erz
// Indizes der auszuwaehlenden Zeilen in V
{
  int i;
  for (i = 0; i < ZZ; i++)
    memcpy(mat[i],VLEind(Erz,V[i])->V,Dim*sizeof(Zahl));
}

void Mat2VL(VList *VL, Matrix mat, int ZZ, int SZ)
// die andere Richtung
{
    int i;

    StartVL(VL);
    for(i=0;i<ZZ;i++)
        memcpy(NewVLE(VL)->V,mat[i],SZ*sizeof(Zahl));
}

void RangIndexVL(VList VL, int dim, int *rang, Zahl *index)
// Wie oben, nur von VList ausgehend
{
    int n;
    Matrix mat;

    n=NinVL(VL);
    MatAlloc(&mat,n,dim);
    VL2Mat(mat,VL,dim);
    RankIndexMat0(mat,n,dim,rang,index);
    FreeMat(mat,n);
}



void FindeInitSimplVL(Matrix SimpErz, VList Gen, int Dim)
// dito fuer Vektorlisten, Rueckgabe als Matrix
{
    Indizes C;
    Matrix Erz;
    int k;

    k=NinVL(Gen);
    MatAlloc(&Erz,k,Dim);


    VL2Mat(Erz,Gen,Dim);

    FindeInitSimpl(SimpErz,Erz,k,Dim);

    FreeMat(Erz,NinVL(Gen));
}




void Transf(VList *New, VList Old, Matrix Tr, Vektor ElT, int EmbDim, int Dim, boolean emb)
// Siehe Diagonalisiere, EmbDim= n, wenn Old in ZZ^n lebt, Dim = Rang des Untermoduls
// U bzw. V oben.
//
{
    VLE *VP, *VP1;
    int i;

    StartVL(New);
    for(VP=Old.St;VP;VP=VP->R)
    {
        VP1=NewVLE(New);
        VxM(VP1->V,VP->V,EmbDim, Tr,Dim);
        if(emb)
            for(i=0;i<Dim;i++)
                VP1->V[i]/=ElT[i];
    }
}

void InvTrans(VList *New, VList Old, Matrix InvTr, Vektor ElT, int EmbDim, int Dim, boolean emb)
// Rücktransformation der Vektoren aus V bzw. U
{
    VLE *VP, *VP1;
    int i;
    Vektor VV;

    VV=(Vektor) calloc(Dim,sizeof(Zahl));

    StartVL(New);
    for(VP=Old.St;VP;VP=VP->R)
    {
        VP1=NewVLE(New);
        if(emb)
        {
            for(i=0;i<Dim;i++)
                VV[i]=VP->V[i]*ElT[i];
            VxM(VP1->V,VV,Dim, InvTr,EmbDim);
        }
        else
            VxM(VP1->V,VP->V,Dim, InvTr,EmbDim);
    }
    free(VV);
}

void InvTransLF(VList *NewF, VList OldF, Matrix Tr, Vektor ElT, int EmbDim, int Dim, boolean emb)
// liftet Linearformen von V bzw. U zurueck nach ZZ^n, n=EmbDim.
// Dass geht im Fall von U möglicherweise erst nach Multiplikation mi´t
// dem Produkt der ElT[i].
// Danach wird aber die Linearform auf ZZ^n primitiv gemacht
{
    VLE *VP, *VP1;
    int i;
    Vektor VV;
    Zahl E;

    VV=(Vektor) calloc(Dim,sizeof(Zahl));

    E=1;
    for(i=0;i<Dim;i++)
        E*=ElT[i];

    StartVL(NewF);
    for(VP=OldF.St;VP;VP=VP->R)
    {
        VP1=NewVLE(NewF);
        if(emb)
        {
            for(i=0;i<Dim;i++)
                VV[i]=VP->V[i]*E/ElT[i];
            MxV(VP1->V,Tr,EmbDim,Dim,VV);
        }
        else
            MxV(VP1->V,Tr,EmbDim,Dim,VP->V);
        MakePrime(VP1->V,EmbDim);
    }
    free(VV);
}


// Simplex-Listen

void StartSimplexL(SimplexList *SimplexL)
{
    SimplexL->St=NULL;
    SimplexL->E=NULL;
}


void MoveSimplex(SimplexList *Target, SimplexList *Source, Simplex *El)
// Bewegt El von Source an das Ende von Target
{

  if (Source->St == El)
    Source->St = El->R;
  else
    El->L->R = El->R;
  if (Source->E == El)
    Source->E = El->L;
  else
    El->R->L = El->L;
  El->R = NULL;
  if (Target->E != NULL) {
    Target->E->R = El;
    El->L = Target->E;
  } else
    El->L = NULL;
  Target->E = El;
  if (Target->St == NULL)
    Target->St = El;
}

Simplex *NewSimplexLE(SimplexList *SimplexL)
// Erzeugt neuen Eintrag am Ende von SimplexL
// und gibt Zeiger darauf zurueck
{
  Simplex *NewS;

  AnzSpx1++;
  if (SimplexFree.St == NULL)
  {
    AnzSpx++;
    SimplexFree.St = calloc(1,sizeof(Simplex));
    SimplexFree.St->Vert = (Indizes) calloc(MaxDim,sizeof(int));
    SimplexFree.E = SimplexFree.St;
    SimplexFree.St->L = NULL;
    SimplexFree.St->R = NULL;
    StartVL(&(SimplexFree.St->Corn));
    StartVL(&(SimplexFree.St->Supp));
  }
  NewS = SimplexFree.St;
  MoveSimplex(SimplexL,&SimplexFree,NewS);
  return NewS;
}


int NinSimplL(SimplexList Triang)
// Zaehlt Elemente in Triang
{
    Simplex *Smp;
    int i;

    i=0;
    for(Smp=Triang.St;Smp;Smp=Smp->R) i++;
    return(i);
}

boolean InVertices(int k, Indizes V, int Dim)
{
    int i;
    for(i=0;i<Dim;i++)
        if(V[i]==k)
            return true;
    return false;
}

void StripCornSupp(SimplexList Triang)
// Entfernt Erzeuger und Stuetzhyperebenen
//
{
    Simplex *Smp;
    for(Smp=Triang.St;Smp;Smp=Smp->R)
    {
        // printf("In Strip %d %d\n",Smp->Corn.St,Smp->Supp.St);
        FreeVL(&(Smp->Corn));
        FreeVL(&(Smp->Supp));
        // printf("In Strip VZ %d FVZ %d\n",VektorZaehler, NinVL(VLFree));
    }
}

void FreeSimplexL(SimplexList *SimplexL)
// Leert die Simplex-Liste
{

  if(!SimplexL->St)
    return;
  StripCornSupp(*SimplexL);
  if (SimplexFree.St == NULL)
    SimplexFree = *SimplexL;
  else {
    SimplexFree.St->L = SimplexL->E;
    SimplexL->E->R = SimplexFree.St;
    SimplexFree.St = SimplexL->St;
  }
  SimplexL->St = NULL;
  SimplexL->E = NULL;
}


void FreeSimplexLE(SimplexList *Triang, Simplex *Smp)
// Verschiebt Simplex in SimplexFree, entfernt Erzeuger und
// Stuetzhyperbenen
{
    FreeVL(&(Smp->Corn));
    FreeVL(&(Smp->Supp));
    MoveSimplex(&SimplexFree,Triang,Smp);
}

void SetCornSuppE(VList OneSkel,Simplex *Smp, int Dim)
// Setzt Erzeuger und Stuetzhyperebenen fuer einzelnes Simplex
{
    Matrix M,S;
    int i;

    MatAlloc(&M,Dim,Dim);
    MatAlloc(&S,Dim,Dim);

    /* SchreibeVL(stdout,OneSkel,Dim);
    printf("++++++++++\n"); */

    FuelleMatrixVL(M,Dim,Dim,OneSkel,Smp->Vert);
    /* for(i=0;i<Dim;i++)
        printf("%d ",Smp->Vert[i]);printf("\n");
    SchreibeMatrix(stdout,M,Dim,Dim);*/
    Smp->Vol=SimpSupp(S,M,Dim);
    Mat2VL(&(Smp->Corn),M,Dim,Dim);
    Mat2VL(&(Smp->Supp),S,Dim,Dim);

    FreeMat(M,Dim);
    FreeMat(S,Dim);
}

void SetCornSupp(VList OneSkel, SimplexList Triang, int Dim)
// Setzt Erzeuger und Stuetzhyperebenen
{
    Simplex *Smp;

    if(Triang.St->Corn.St) // Bereits gesetzt
        return;
    for(Smp=Triang.St;Smp;Smp=Smp->R)
        SetCornSuppE(OneSkel,Smp,Dim);
}

int FacetOfSimplex(Indizes Simp, Indizes Facet, int Dim)
// Prueft of Facet in Simp enthalten.
// Wenn nein, Rueckgabe -1
// Wenn ja, Rueckgabe der zu Facet komplementaren Ecke
{
    int i,j;

    for(i=0;i<Dim-1;i++)
    {
        if(Facet[i]>Simp[i]) // Simp[i] nicht in Facet
            break;
        if(Facet[i]<Simp[i]) // Facet[i] nicht in Simp
            return -1;
    }
    j=i; // Falls durchgelaufen Simp[Dim-1] nicht in Facet
    for(;i<Dim-1;i++)
        if(Facet[i]!=Simp[i+1])
            return -1;
    return(Simp[j]);
}

void CopySimplexLE(SimplexList *Target, Simplex *Smp, int Dim)
// Kopiert einzelnes Simplex nach Target
{
    Simplex *NewSmp;

    NewSmp=NewSimplexLE(Target);
    memcpy(NewSmp->Vert,Smp->Vert,Dim*sizeof(int));
    NewSmp->Vol=Smp->Vol;
    CopyVL(&(NewSmp->Corn), Smp->Corn, Dim);
    CopyVL(&(NewSmp->Supp), Smp->Supp, Dim);
}

void CopySimplexList(SimplexList *Target,SimplexList Source, int Dim)
// Kopiert ganze Liste
{
    Simplex *Smp;

    StartSimplexL(Target);

    for(Smp=Source.St;Smp;Smp=Smp->R)
        CopySimplexLE(Target,Smp,Dim);
}



void SchreibeSimplexL(FILE *Aus,SimplexList Triang, int Dim)
// Analog SchreibeVL
{
    Simplex *Smp;
    int i,k;

    k=NinSimplL(Triang);
    for(Smp=Triang.St;Smp;Smp=Smp->R)
    {
        for(i=0;i<Dim;i++)
            fprintf(Aus,"%d ",Smp->Vert[i]+1);
        fprintf(Aus,"\n");
    }

}

void FSchreibeSimplexL(FILE *Aus,SimplexList Triang, int Dim)
// Analog FSchreibeVL
{
    fprintf(Aus,"%d\n",NinSimplL(Triang));
    fprintf(Aus,"%d\n",Dim);
    SchreibeSimplexL(Aus,Triang,Dim);
}

void SchreibeTriang(FILE *Aus,VList OneSkel, int Dim, SimplexList Triang)
{
    FSchreibeVL(Aus,OneSkel,Dim);
    FSchreibeSimplexL(Aus,Triang,Dim);

}


void HaengeTriangAn(char *database, VList OneSkel, int Dim, SimplexList Tri, char *caption)
// haengt Faecher an Datenbasis an.
// Wartet bis diese nicht mehr blockiert ist,
// angezeigt durch Existenz von File mit ".lock" am Ende
{
    char filename[1000], filenamelock[1000];
    int i,dummy;
    FILE *aus,*lock;

    sprintf(filename,"%s_%d.base",database,Dim);
    strcpy(filenamelock,filename);
    strcat(filenamelock,".lock");

    while(access(filenamelock, F_OK) == 0)
    {
        printf("File %s locked\n",filename);
        for(i=0;i<1000000;i++)
        {
            dummy=i;
        }
    }

    lock=fopen(filenamelock,"w");
    aus=fopen(filename,"a");
    fprintf(aus,"%s\n",caption);
    SchreibeTriang(aus,OneSkel,Dim,Tri);
    fclose(aus);
    fclose(lock);
    remove(filenamelock);
}


void LiesSimplexL(FILE *Ein,SimplexList *Triang, int Dim)
// Analog LiesVL, aber Anzahl Simplizes und Dim werden intern gelesen
{
    Simplex *Smp;
    int i,j,k,n;

    StartSimplexL(Triang);

    fscanf(Ein,"%d",&k);
    fscanf(Ein,"%d",&n); // Diemsnions-Angabe ueberlesen
    for(i=0;i<k;i++)
    {
        Smp=NewSimplexLE(Triang);
        for(j=0;j<Dim;j++)
        {
            fscanf(Ein,"%d ",&n);
            Smp->Vert[j]=n-1;
        }
    }
}

void LiesTriang(FILE *Ein,VList *OneSkel, int N, int Dim, SimplexList *Triang)
{
    LiesVL(Ein,OneSkel,N,Dim);
    LiesSimplexL(Ein,Triang,Dim);
}

void FLiesTriang(FILE *Ein,VList *OneSkel, int *N, int *Dim, SimplexList *Triang)
{
    FLiesVL(Ein,OneSkel,N,Dim);
    LiesSimplexL(Ein,Triang,*Dim);
}


// Spezielle Debugging-Hilfen nur fuer 32bit genacht

#ifndef LONGLONG

void SchreibeVLEDebug(FILE *Aus, VLE *VP, int Dim, int LVSize)
// schreibt den gesamten VLE auf Aus
{
    int i;
    for(i=0;i<Dim;i++)
        fprintf(Aus,"%d ", VP->V[i]);
    fprintf(Aus,"|| %d ",VP->Mark);
    fprintf(Aus,"? %d ",VP->Mark1);
    fprintf(Aus,"? %d ",VP->Mark2);
    fprintf(Aus,"? %d || ",VP->Mark3);
    for(i=0;i<LVSize;i++)
        fprintf(Aus,"%d ", VP->LV[i]);
    fprintf(Aus,"V %d LVP %d LV %d",VP->V,VP->LVP,VP->LV);
    fprintf(Aus,"\n");
}

void SchreibeVLDebug(FILE *Aus, VList VL, int Dim, int LVSize)
// macht das fuer die gesamte Liste
{
    VLE *VP;
    for(VP=VL.St;VP;VP=VP->R)
        SchreibeVLEDebug(Aus,VP,Dim,LVSize);
}

#endif


// Listen fuer "lange" vektoren


LongV* MakeLongV(int Size)
{
  LongV *LV;
  LV= (LongV*) calloc(1,sizeof(LongV));
  LV->Size = Size;
  if(Size>MaxLVSize)
      MaxLVSize=Size; // fuer ststistische Zwecke
  if(FixedLVSize)
  {
        if(Size>FixedLVSize)
        {
            printf("Size = %d > FixedLVSize = %d\n", Size, FixedLVSize);
            exit(1);
        }
        else
            Size=FixedLVSize;
  }
  LV->V= (Vektor) calloc(Size,sizeof(Zahl));
  return LV;
}

LongVLE *NewLongVLE(LongVList *LongVL, int Size)
{
  LongVLE *NewLongV;
  int i;

    AnzLV1++;
    if (LongVLFree.St == NULL)
    {
        AnzLV++;
        LongVLFree.St = calloc(1,sizeof(LongVLE));
        LongVLFree.E = LongVLFree.St;
        LongVLFree.St->L = NULL;
        LongVLFree.St->R = NULL;
        if(FixedLVSize)
            LongVLFree.St->LV=MakeLongV(Size);
    }
    NewLongV = LongVLFree.St;
    MoveLongVLE(LongVL,&LongVLFree,NewLongV);
    /*if(FixedLVSize)
        for(i=0;i<FixedLVSize;i++)
            (NewLongV->LV)->V[i]=i;*/

    if(!FixedLVSize)
        NewLongV->LV=MakeLongV(Size);
    NewLongV->LV->Size=Size;
    return NewLongV;
}

int NinLVL(LongVList vectors)
// Zaehlt Elemente in Vectors
{
    LongVLE *VP;
    int i;

    i=0;
    for(VP=vectors.St;VP;VP=VP->R) i++;
    return(i);
}

void FreeLongVL1(LongVList *LongVL)
{
  if(!LongVL->St) return;
  if (LongVLFree.St == NULL)
    LongVLFree = *LongVL;
  else {
    LongVLFree.St->L = LongVL->E;
    LongVL->E->R = LongVLFree.St;
    LongVLFree.St = LongVL->St;
  }
  LongVL->St = NULL;
  LongVL->E = NULL;
}

void FreeLongVLE(LongVList *LVL, LongVLE *LVLP)
{
    if(!FixedLVSize)
        free((LVLP->LV)->V);
    LVLP->LV->Size=0;
    MoveLongVLE(&LongVLFree, LVL,LVLP);
}


void FreeLongVL(LongVList *LongVL)
{
  LongVLE *LongV_;

  LongV_ = LongVL->St;
  while (LongV_ != NULL)
  {
    if(!FixedLVSize)
        free(LongV_->LV->V);
    LongV_->LV->Size=0;
    LongV_ = LongV_->R;
  }
  FreeLongVL1(LongVL);
}

void SchreibeLVL(FILE *aus, LongVList LVL, int SZ)
// analog SchreibeVL, Laenge der Vektoren hier SZ
{
    LongVLE *LVP;
    for(LVP=LVL.St;LVP;LVP=LVP->R)
        SchreibeV(aus,LVP->LV->V,SZ);
}

void LVL2VL(VList *VL,LongVList LVL, int SZ)
// Uebertraegt LongVList in VList
{
    LongVLE *LVLP;

    StartVL(VL);
    if(SZ>MaxDim)
        printf("Alarm SZ (%d) > MaxDim (%d)\n",SZ,MaxDim);
    for(LVLP=LVL.St;LVLP;LVLP=LVLP->R)
        memcpy(NewVLE(VL)->V,(LVLP->LV)->V,SZ*sizeof(Zahl));
}

// Langvektoren fuer Vektorlisten

void AllocLVE(VLE *VP,int Size)
// Allokiert LangVektor fuer VP der Laenge Size unter Beroecksichtigung
// von FixedLVSize
{

    if(Size>MaxLVSize)
      MaxLVSize=Size; // fuer ststistische Zwecke

    if(FixedLVSize)
    {
        if(Size>FixedLVSize)
        {
            printf("Size = %d > FixedLVSize = %d\n", Size, FixedLVSize);
            exit(1);
        }
        if(VP->LVP) // hat schon einen Langvektor
        {
            return;
        }

        VP->LVP=NewLongVLE(&LVPool,Size);
        VP->LV=VP->LVP->LV->V;
    }
    else
    {
        if(VP->LVP)
            free(VP->LV); // hier koennte nam auch realloc einsetzen
        VP->LV=(Vektor) calloc(Size,sizeof(Zahl));
        /* if(!VP->LV)
        {
            printf("Platzmangel\n");
        } */
        VP->LVP=(LongVLE *) 1; // um es ungleich Null zu machen
    }
}

VLE* CopyVLE(VList *Target, VLE *VP, int Dim, int LVSize)
// Kopiert ein einzelnes VLE
// LVs beruecksichtigt, aber Marks nicht
{
    VLE *VP1;

    VP1=NewVLE(Target);
    memcpy(VP1->V,VP->V,Dim*sizeof(Zahl));
    AllocLVE(VP1,LVSize);
    memcpy(VP1->LV,VP->LV,LVSize*sizeof(Zahl));
    return(VP1);
}


void AllocLV2VL(VList Gen, int Size)
// Allokiert LangVektoren der Laenge Size für die Vektoren in Gen
{
    VLE *VP;
    for(VP=Gen.St;VP;VP=VP->R)
        AllocLVE(VP,Size);
}

void AddValLV(VList Gen, VList Supp, int Dim)
// Fuegt den Vektoren in Gen "Wertvektoren" hinzu
// die sich als Skalarprodukte der Elemente in Gen
// mit denen in Supp ergeben
{
    VLE *VP, *VP1;
    int i,NSupp;

    NSupp=NinVL(Supp);
    AllocLV2VL(Gen, NSupp);
    for(VP=Gen.St;VP;VP=VP->R)
    {
        i=0;
        for(VP1=Supp.St;VP1;VP1=VP1->R)
        {
            VP->LV[i]=SkalProd(VP->V,VP1->V,Dim);
            i++;
        }
    }
}


// Reduktion

int RedInsert(VList *Bas, int Dim, Vektor New, VList Supp, int NSupp)
// Reduziert New gegen Bas und umgekehrt, Test ob Differenz im Kegel
// ueber Stuetzformen in Supp, NSupp nur zur Zeitersparnis mitgefuehrt
{
    VLE *VP,*VP1,*VPNext,*VP1Next;
    int i, test;
    Vektor Val;
    boolean inserted;

    for(i=0;i<Dim;i++)
        if(New[i]) break;
    if(i==Dim)
        return(0);  // Nullvektor

    Val=calloc(NSupp,sizeof(Zahl));

    i=0;
    for(VP=Supp.St;VP;VP=VP->R)
    {
        Val[i]=SkalProd(VP->V,New,Dim);
        i++;
    }

    inserted=false;
    for(VP=Bas->St;VP;VP=VPNext)
    {
        VPNext=VP->R;

        test=VCompare(VP->LV,Val,NSupp);
        if(test<=0) // Neuer Vektor ist >= und erledigt
        {
            free(Val);
            return(0);
        }
        if(test==1) // Neuer Vektor <
        {
            if(!inserted) // muss noch eingesetzt werden
            {
                memcpy(VP->V,New,Dim*sizeof(Zahl));
                memcpy(VP->LV,Val,NSupp*sizeof(Zahl));
                inserted=true;
                continue;
            }
            FreeVLE(Bas,VP);
        }
    }

    // falls !inserted, mit allen unvergleichbar

    if(!inserted)
    {
        VP=NewVLE(Bas); // mit keinem vergleichbar, also einsetzen
        AllocLVE(VP,NSupp);
        memcpy(VP->V,New,Dim*sizeof(Zahl));
        memcpy(VP->LV,Val,NSupp*sizeof(Zahl));
    }
    free(Val);
    return(1);
}

boolean RedInsertPosOrth(VList *Bas, int Dim, Vektor New)
// Reduziert New gegen Bas und umgekehrt, Test ob Differenz im positiven Orthanten
// Gibt "true" zuruck, falls New NICHT reduziert wird
{
    VLE *VP,*VPNext;
    int i, test;
    boolean inserted;

    for(i=0;i<Dim;i++)
        if(New[i]) break;
    if(i==Dim)
        return false;  // Nullvektor

    inserted=false;
    for(VP=Bas->St;VP;VP=VPNext)
    {
        VPNext=VP->R;

        test=VCompare(VP->V,New,Dim);
        if(test<=0) // Neuer Vektor ist >= und erledigt
        {
            return false;
        }
        if(test==1) // Neuer Vektor <
        {
            if(!inserted) // muss noch eingesetzt werden
            {
                memcpy(VP->V,New,Dim*sizeof(Zahl));
                inserted=true;
                continue;  // Rest überspringen
            }
            FreeVLE(Bas,VP); // altes Element loeschen, falls schon "inserted"
        }
    }
    if(!inserted)
    {
        VP=NewVLE(Bas); //
        memcpy(VP->V,New,Dim*sizeof(Zahl));
    }
    return true;
}


boolean Irred(VList VL, Vektor V, int Bis)
// Reduziert WERTVEKTOR V gegen VL durch Vergleich mit Wertvektoren in VL
// Gibt true zurueck, falls irreduzibel
{
    VLE *VP;
    int test,i;

    for (VP=VL.St;VP;VP=VP->R)
    {
        test=VCompare(VP->LV,V,Bis);

        if(test<=0)
            return(false); // Neuer vektor ist reduzibel
        if(test==1)
            VP->Mark2=1; // Markiere VP als reduzibel
    }
    return(true);
}



// Simpliziale Monoide


void PrepPar(VList *InPar, int Dim, Matrix Inv, Zahl Nenner)
// Bereitet die Darstellung der Elemente in Par(U) vor, wobei
// U indirekt ueber die Transponierte der Inversen gegeben ist.
// Restklassen von InPar repraesentiert. Diese Liste wird transformiert.
// Rueckgabe von Vektoren, die mit Nenner multipliziert und in der Basis
// von U dargestellt
{
    VLE *VP;
    int i;

    for(VP=InPar->St;VP;VP=VP->R)
    {
        MxV(VP->V,Inv,Dim,Dim,VP->V);
        for(i=0;i<Dim;i++)
        {
            VP->V[i] %= Nenner;
            if(VP->V[i]<0) VP->V[i]+=Nenner;
        }
    }
}

void MakePar(VList *InPar, int Dim, Matrix U, Zahl Nenner)
// Setzt das Werk von PepPar oder PrepParHilbBas fort und bestimmt DSarstellung in
// kanonischer Basis und dividiert durch Nenner
{
    VLE *VP;
    int i,j;
    Zahl APrioriBd;

    APrioriBd=0;

    for(i=0;i<Dim;i++)
        for(j=0;j<Dim;j++)
            if(APrioriBd < Zabs(U[i][j])) APrioriBd=Zabs(U[i][j]);
    APrioriBd*=Dim;

    for(VP=InPar->St;VP;VP=VP->R)
    {
        VxM(VP->V,VP->V,Dim,U,Dim);
        for(i=0;i<Dim;i++)
        {
            VP->V[i]/=Nenner;
            if(Zabs(VP->V[i])>APrioriBd)
            {
                printf("Arithmetic overflow apriori bound\n");
                exit(1);
            }
        }
    }
}

void Par(VList *InPar, Matrix U, int Dim)
// Bestimmt Repräsentantensystem von ZZ^n/U in Par(U)
{
    Zahl Nenner;
    Matrix Inv;
    Vektor Diag;
    VList ResCl;

    MatAlloc(&Inv,Dim,Dim);
    Diag= calloc(Dim,sizeof(Zahl));

    InvertZ(Inv, U, Dim, Diag, &Nenner);

    Diag2ResCl(InPar, Diag, Dim, 1);
    PrepPar(InPar, Dim, Inv, Nenner);
    MakePar(InPar, Dim,U,Nenner);

    FreeMat(Inv,Dim);
    free(Diag);
}

void PrepParHilbBasOld(VList *InPar, Vektor Diag, Matrix Inv, int Dim, Zahl Nenner)
// Analogon zu PrepPar, aber enthält Diag2ResCl und fuer Hilbert-Basis-Bestimmung
// Elemente werden direkt reduziert
// alte Version
{
    int i,Last;
    VLE *VP;
    Vektor C,V;

    C=(Vektor) calloc(Dim,sizeof(Zahl));
    V=(Vektor) calloc(Dim,sizeof(Zahl));

    StartVL(InPar);

    while(1)
    {
        Last = -1;
        for (i=0;i<Dim;i++) if (C[i]<Diag[i]-1) Last = i;
        if (Last == -1) break;
        C[Last]++;
        for(i=Last+1;i<Dim;i++) C[i]=0;

        MxV(V,Inv,Dim,Dim,C);
        for(i=0;i<Dim;i++)
        {
            V[i] %= Nenner;
            if(V[i]<0) V[i]+=Nenner;
        }

        RedInsertPosOrth(InPar,Dim,V);
    }

    free(C);
    free(V);
}

void MakeFacesSimplex(SimplexList *Faces, int Dim)
// Macht alle Seiten des Simplex mit den Ecken 0,...,Dim
{
    int i,j,k,n;
    Simplex *Simpl,*Start,*End,*NewS;

    StartSimplexL(Faces);

    NewS=NewSimplexLE(Faces);
    for(n=0;n<Dim;n++)
        NewS->Vert[n]=-1;

    for(i=0;i<Dim;i++)
    {
        if(i)
            Start=End->R;
        else
            Start=Faces->St;
        End=Faces->E;
        for(Simpl=Start;;Simpl=Simpl->R)
        {
            if(!i)
                j=0;
            else
                j=Simpl->Vert[i-1]+1;
            for(k=j;k<Dim;k++)
            {
                NewS=NewSimplexLE(Faces);
                memcpy(NewS->Vert,Simpl->Vert,i*sizeof(int));
                NewS->Vert[i]=k;
                for(n=i+1;n<Dim;n++)
                    NewS->Vert[n]=-1;
            }
            if(Simpl==End)
                break;
        }
    }
}



void HilbParOld(VList *InPar, Matrix U, int Dim)
// Bestimmt Hilbert-Basis fuer ganzen Abschluss des durch U erzeugten
// simplizialen Monoids, aber ohne die Erzeuger von U !!
// alte Version
{
    Zahl Nenner;
    Matrix Inv;
    Vektor Diag;

    MatAlloc(&Inv,Dim,Dim);
    Diag= calloc(Dim,sizeof(Zahl));

    InvertZ(Inv, U, Dim, Diag, &Nenner);
    PrepParHilbBasOld(InPar, Diag, Inv, Dim, Nenner);
    MakePar(InPar, Dim,U,Nenner);

    FreeMat(Inv,Dim);
    free(Diag);
}

typedef struct KandV{
        Vektor V;
        Zahl TotGrad;
    } KandV;

int VEKTLAENGE;

int VecVergl(const void *m,const void *n)
{
    int i;
    Vektor Vm,Vn;
    Zahl Tm,Tn;

    Vm=(*(KandV *)m).V;
    Vn=(*(KandV *)n).V;
    Tm=(*(KandV *)m).TotGrad;
    Tn=(*(KandV *)n).TotGrad;

    if(Tm>Tn)
        return 1;
    if(Tn>Tm)
        return -1;

    for(i=0;i<VEKTLAENGE;i++)
    {
        if(Vm[i]==Vn[i])
            continue;
        if(Vm[i]>Vn[i])
            return 1;
        return -1;
    }
    return 0;
}

void PrepParHilbBas(VList *InPar, Vektor Diag, Matrix Inv, int Dim, Zahl Nenner)
// Analogon zu PrepPar, aber enthält Diag2ResCl und Reduktion
{
    int i,j,Last,NKand,Gr;
    VLE *VP, *VP1;
    Vektor C,V;
    KandV *Kand;
    VList FirstList;
    boolean Reduziert;
    Zahl TGrad;
    Zahl S, Nvgl;

    C=(Vektor) calloc(Dim,sizeof(Zahl));
    V=(Vektor) calloc(Dim,sizeof(Zahl));

    StartVL(InPar);
    StartVL(&FirstList);

    while(1)
    {
        Last = -1;
        for (i=0;i<Dim;i++) if (C[i]<Diag[i]-1) Last = i;
        if (Last == -1) break;
        C[Last]++;
        for(i=Last+1;i<Dim;i++) C[i]=0;

        MxV(V,Inv,Dim,Dim,C);
        for(i=0;i<Dim;i++)
        {
            V[i] %= Nenner;
            if(V[i]<0) V[i]+=Nenner;
        }
        VP=NewVLE(&FirstList);
        memcpy(VP->V,V,Dim*sizeof(Zahl));
    }

    NKand=NinVL(FirstList);
    Kand = (KandV*) calloc(NKand,sizeof(KandV));

    VP=FirstList.St;
    for(i=0;i<NKand;i++)
    {
        Kand[i].V=VP->V;
        S=0;
        for(j=0;j<Dim;j++)
            S+=VP->V[j];
        Kand[i].TotGrad=S;  // Totalgrad ist 0-te Komponente
        VP=VP->R;
    }


    VEKTLAENGE=Dim;
    qsort(Kand,NKand,sizeof(KandV),VecVergl);

    // Nvgl=0;

    for(i=0;i<NKand;i++)
    {
        TGrad=Kand[i].TotGrad; // TotalGrad übernehmen

        Gr=0;
        Reduziert=false;


        for(VP=InPar->St;VP;VP=VP->R)
        {
            // Nvgl++;

            if(2*VP->Mark>TGrad)
                break;

            if(Kand[i].V[Gr]<VP->V[Gr])
               continue;

            for(j=0;j<Dim;j++)
            {
                if(Kand[i].V[j]<VP->V[j])
                    break;
            }
            Gr=j;
            if(j==Dim) // reduziert
            {
                Reduziert=true;
                MakeFirst(InPar,VP);
                break;
            }
        }
        if(Reduziert)
            continue;


    InHilbBas:
        VP1=NewVLE(InPar); // jetzt irreduzibel
        memcpy(VP1->V,Kand[i].V,Dim*sizeof(Zahl));
        VP1->Mark=TGrad;
    }

    // printf("%d %d %lld\n",NKand,NinVL(*InPar),Nvgl);
    fflush(stdout);

    FreeVL(&FirstList);
    free(Kand);
    free(C);
    free(V);
}

void HilbPar(VList *InPar, Matrix U, int Dim)
// Bestimmt Hilbert-Basis fuer ganzen Abschluss des durch U erzeugten
// simplizialen Monoids, aber ohne die Erzeuger von U !!
{
    Zahl Nenner;
    Matrix Inv;
    Vektor Diag;

    MatAlloc(&Inv,Dim,Dim);
    Diag= calloc(Dim,sizeof(Zahl));

    InvertZ(Inv, U, Dim, Diag, &Nenner);
    PrepParHilbBas(InPar, Diag, Inv, Dim, Nenner);
    MakePar(InPar, Dim,U,Nenner);

    FreeMat(Inv,Dim);
    free(Diag);
}




// Duale Kegel und Triangulierungen

void MakeNewSupp(LongVList *Val, int size, int Tot, LongVLE *ValP, LongVLE *ValN)
{

        Zahl t1,t2,g;
        LongVLE *ValTest;
        int i;

        ValTest=NewLongVLE(Val, size);
        t1=(ValP->LV)->Mark;
        t2=(ValN->LV)->Mark;
        g=ggt(t1,t2);
        if(g!=1)
        {
            t1/=g;
            t2/=g;
        }
        for (i=0;i<=Tot;i++)
              (ValTest->LV)->V[i]=t1*(ValN->LV)->V[i]-t2*(ValP->LV)->V[i];
        (ValTest->LV)->Mark=0;
        MakePrime((ValTest->LV)->V,Tot);
}

int DIM2; // globale Variable

int RegVerglDim2(const void *m,const void *n)
{
    return(memcmp(m,n,DIM2*sizeof(int)));
}

void TransformValues(LongVList *Val, int Dim, int Tot, Matrix Erz, Indizes E)
//
// Fuehrt den wesentliche Transformationsschritt fuer Bestimmung
// des dualen Kegels durch
// dabei wird aus der Kiste der Stuetzhyperebenen von RR_+(x_1,...,x_n)
// entsprechende Liste fuer RR_+(x_1,...,x:n,x_{n+1}) erzeugt
// die ersten Dim Elemente eines jeden LV von Val sind die Koordinaten
// der Stuetzhyperebenen, danach folgen die werte auf den x_i
// Tot gibt jeweils die nach der Transformation erreichte Laenge
// der Vektoren in Val an
// Die Werte auf x_[n+1} sind werden in den "Marks" transportiert
// Erz sind die Erzeuger des Kegels und E ist das "Scharnier" zwischen
// der Reihenfolge in Erz und der Reihenfolge, in der die Vektoren
// verarbeitet werden: E[Dim+i-1] gibt die Zeile in erz an, die zum i-ten
// vektor in der Verarbeitungsreihenfolge gehört.
{
  LongVLE *ValP, *ValN, *ValTest;
  Indizes ZB, ZP;
  int i, NZB, NZP, size;
  boolean Einzig;


  int NSimpl,j,NVal,NPos,NNeg,P,N,NAll,NNonSimp,NSubF,NNeutr,Ne,
             k,NNegSimp, NegNonSimp,Start,TestMSize,NEff,NMono,NegNeutr,NegNeg;
  boolean *PosSimp, *NegSimp, NV, IsSimp,RangTest,PMono,Gefunden;

  Indizes TestMZ;
  Indizes PVis,NVis,DVis,Search;

  LongVLE **NonSimp, **Pos, **Neg, **Neutr;

  Matrix NBlReg,TestM;

  typedef struct RegE
  {
    int Z[Dim-2];
    LongVLE *VP;
  } RegE;

  RegE *NReg, *SubFac;

  // if(Tot-Dim>17) exit(0);


  // printf("-----------------\n");

  DIM2=Dim-2; // fuer Vergleich oben

  size=(Val->St->LV)->Size;

  NVal=NinLVL(*Val);
  // printf("Tot-Dim %d NSupp %d\n",Tot- Dim,NVal);

  PosSimp=calloc(NVal,sizeof(boolean));
  NegSimp=calloc(NVal,sizeof(boolean));
  NonSimp=calloc(NVal,sizeof(LongVLE*));
  Pos=calloc(NVal,sizeof(LongVLE*));
  Neg=calloc(NVal,sizeof(LongVLE*));
  Neutr=calloc(NVal,sizeof(LongVLE*));
  ZP=(Indizes) calloc(size, sizeof(int));
  ZB=(Indizes) calloc(size, sizeof(int));
  Search=(Indizes) calloc(Dim-2, sizeof(int));

  PVis=calloc(Tot,sizeof(boolean));
  NVis=calloc(Tot,sizeof(boolean));
  DVis=calloc(Tot,sizeof(boolean));

  TestMSize=0;

  NNeg=-1;
  NPos=-1;
  NNeutr=-1;
  NSimpl=0;
  NNonSimp=-1;
  NNegSimp=0;

  for(ValP=Val->St;ValP;ValP=ValP->R)
  {
     (ValP->LV)->V[Tot]=(ValP->LV)->Mark;
     IsSimp=false;
     NZP=0;
     for(i=Dim;i<Tot;i++)
       if(!(ValP->LV->V[i]))
       {
           NZP++;
           if((ValP->LV)->Mark>0)
           {
                PVis[i]=true;
                continue;
           }
           if((ValP->LV)->Mark<0)
                NVis[i]=true;
       }
    if(NZP==Dim-1)
    {
       IsSimp=true;
       NSimpl++;
    }
    else
    {
        NNonSimp++;
        NonSimp[NNonSimp]=ValP;
    }
    if((ValP->LV)->Mark>0)
    {
        NPos++;
        Pos[NPos]=ValP;
        // printf("POs %d\n",ValP);
        if(IsSimp)
            PosSimp[NPos]=true;
        continue;
    }
    if((ValP->LV)->Mark<0)
    {
        NNeg++;
        Neg[NNeg]=ValP;
        if(IsSimp)
        {
            NegSimp[NNeg]=true;
            NNegSimp++;
        }
        continue;
    }

    NNeutr++;
    Neutr[NNeutr]=ValP;
  }
  NPos++;
  NNeg++;
  NNeutr++;
  NNonSimp++;

  // printf("Pos %d Neg %d Neutr %d Simpl %d\n",NPos,NNeg,NNeutr,NSimpl);

  j=0;
  for(i=Dim;i<Tot;i++)
   if(PVis[i]&&NVis[i])
   {
        j++;
        DVis[i]=true;
   }
  // printf("D %d\n",j);

  NEff=0;
  NMono=0;



  for(ValP=Val->St;ValP;ValP=ValP->R)  // Nur Statistik
  {
    if(!(ValP->LV)->Mark)
        continue;
    NZB=0;
    for(i=Dim;i<Tot;i++)
        if(DVis[i]&& !(ValP->LV)->V[i])
        {
            ZB[NZB]=i;
            NZB++;
        }
    if(NZB>=Dim-2)
    {

            NEff++;
            if(NZB==Dim-2)
                NMono++;
    }
  }

  // printf("NEff %d NMono %d\n",NEff,NMono); // Statistik



  RangTest=false;

  if(NNonSimp>Dim*Dim*Dim/6)  // Dim*Dim*Dim
  {
        RangTest=true;
        // printf("RangTest\n");
        TestMZ=(Indizes) calloc(Tot-Dim,sizeof(int));
  }

  NReg=calloc(NNegSimp*(Dim-1), sizeof(RegE));

  // NNegSimp=0;
  NegNonSimp=0;
  NSubF=0;
  for(N=0;N<NNeg;N++)
  {

    j=0;
    for(i=Dim;i<Tot;i++)
        if(DVis[i] && !Neg[N]->LV->V[i])
        {
            ZB[j]=i; // Registrieren der Erzeuger in Facette
            j++;
        }
    if(j<Dim-2)
        continue;  // Eliminieren ineffektiver Facetten

    if(!NegSimp[N])  // Nicht simplizial
    {
        Neg[NegNonSimp]=Neg[N];
        NegNonSimp++;
        continue;
    }

    ValP=Neg[N];  // Simplizial
    // printf("Neg %d\n",ValP);

    if(j==Dim-2)  // Nur eine Subfacette ist Kandidat
    {
        NReg[NSubF].VP=ValP;
        memcpy(NReg[NSubF].Z,ZB,j*sizeof(int));
        NSubF++;
        continue;
     }

     for(i=0;i<Dim-1;i++)  // Alle Subfacetten sind Kandidaten
     {
        NReg[NSubF].VP=ValP;
        memcpy(NReg[NSubF].Z,ZB,i*sizeof(int));
        memcpy(&(NReg[NSubF].Z[i]),&(ZB[i+1]),(Dim-2-i)*sizeof(int));
        NSubF++;
     }
  }

  qsort(NReg,NSubF,sizeof(RegE),RegVerglDim2);
  // printf("Sortiert\n");

  j=0;
  NegNeg=0;
  NegNeutr=0;
  for(N=0;N<NSubF;N++)
  {
    if( N<NSubF-1 && !memcmp(NReg[N].Z,NReg[N+1].Z,(Dim-2)*sizeof(int)))
    {
        N++; // Auslassen, weil doppelt
        continue;
    }

    // jetzt nur einfach vorhanden

    Gefunden=false;
    for(Ne=0;Ne<NNeutr;Ne++) // pruefen, ob benachbart zu neutraler Facette
    {
        ValP=Neutr[Ne];
        for(i=0;i<Dim-2;i++)
            if(ValP->LV->V[NReg[N].Z[i]])
                break;
        if(i<Dim-2)
            continue;
        Gefunden=true;
        NegNeutr++;
        break;
    }
    if(Gefunden) continue;

    for(Ne=0;Ne<NegNonSimp;Ne++) // pruefen, ob benachbart zu nonsimpl neg Facette
    {
        ValP=Neg[Ne];
        for(i=0;i<Dim-2;i++)
            if(ValP->LV->V[NReg[N].Z[i]])
                break;
        if(i<Dim-2)
            continue;
        NegNeg++;
        Gefunden=true;
        break;
    }
    if(Gefunden) continue;

    memcpy(&(NReg[j]),&(NReg[N]),sizeof(RegE)); // zu pos Fac benachbart
    j++;
  }
  //printf("%d SubFac %d doppelt %d NegNonSimp %d Neutr %d Pos\n",
  //             NSubF,NSubF-NegNeutr-NegNeg-j, NegNeg,NegNeutr,j);

  NSubF=j; // neu setzen


  for(P=0;P<NPos;P++) // Simpliziale Pos vs. simpliziale Neg
  {
    if(!PosSimp[P])
        continue;
    ValP=Pos[P];
    NZP=0;
    for (i=Dim;i<Tot; i++)
    {
      if (!(ValP->LV)->V[i] && DVis[i])
      {
            ZP[NZP]=i; // Enthaelt die Indizes, wo ValP = 0 ist
            NZP++;
      }
    }

    if(NZP<Dim-2)  // ineffektiv
        continue;

    memcpy(Search,ZP,(Dim-2)*sizeof(int));
    SubFac=(RegE*) bsearch(Search,NReg,NSubF,sizeof(RegE),RegVerglDim2);
    if(SubFac)
    {
        ValN=SubFac->VP;
        MakeNewSupp(Val,size,Tot,ValP, ValN);
        SubFac->VP=NULL;  // stillegen, weil verbraucht
    }

    if(NZP==Dim-2)  // Nur eine Subfacette Kandidat
        continue;

    for(i=Dim-2;i>0;i--)
    {
        Search[i-1]=ZP[i];
        SubFac=(RegE*) bsearch(Search,NReg,NSubF,sizeof(RegE),RegVerglDim2);
        if(SubFac)
        {
            ValN=SubFac->VP;
            MakeNewSupp(Val,size,Tot,ValP, ValN);
            SubFac->VP=NULL; // stillegen, weil verbraucht
        }
    }
  } // PosSimpl vs NegSimpl

  j=0;
  for(i=0;i<NSubF;i++)
    if(NReg[i].VP)
        j++;
  // printf("Subfac PosNonSimp %d\n",j);

  for(N=0;N<NSubF;N++) //NegSimpl vs PosNonSimpl
  {
        if(!NReg[N].VP)  // Subfacette schon erledigt
            continue;

        for(P=0;P<NPos;P++)
        {
            if(PosSimp[P])
                continue; // simplizial
            ValP=Pos[P];
            for(i=0;i<Dim-2;i++)
                if(ValP->LV->V[NReg[N].Z[i]])
                    break;
            if(i<Dim-2)
                continue;
            ValN=NReg[N].VP;
            MakeNewSupp(Val,size,Tot,ValP, ValN);
            break; // Subfacette verbraucht
        }


  } //NegSimpl vs PosNonSimpl

  // printf("NegSimpl abgearbeitet\n");

  for(P=0;P<NPos;P++) // Pos vs NegNonSimpl
  {
    ValP=Pos[P];
    NZP=0;
    for (i=Dim;i<Tot; i++)
    {
      if (!(ValP->LV)->V[i] && DVis[i])
      {
            ZP[NZP]=i; // Enthaelt die Indizes, wo ValP = 0 ist
            NZP++;
      }
    }

    if(NZP<Dim-2)  // ineffektiv
        continue;

    PMono=false;
    if(NZP==Dim-2)
        PMono=true;

    for(N=0;N<NegNonSimp;N++)
    {
        ValN=Neg[N];
        NZB=0;
        NV=false;
        for (i=0;i<NZP; i++)
        {
              if(!(ValN->LV)->V[ZP[i]])
              {
                  ZB[NZB]=ZP[i]; // Enthaelt die Indizes, wo ValP und ValN beide 0 sind
                  NZB++;
              }
              else
              {
                if(PosSimp[P]&&NV)
                    break;
                NV=true;
              }

        }
        if (NZB<Dim-2) continue; // Durchschnitt keine Subfacette


        if(PosSimp[P])  // Falls Pos simplizial, ist es Subfacette
        {
            MakeNewSupp(Val,size,Tot,ValP, ValN);
            if(PMono)
                goto PDone;
            continue;
        }

        if(RangTest) // Falls Pos nicht simplizial weiter testen
        {
            if(NZB>TestMSize)
            {
                if(TestMSize)
                    FreeMat(TestM,TestMSize);
                MatAlloc(&TestM,NZB,Dim);
                TestMSize=NZB;
            }
            for(i=0;i<NZB;i++)
                TestMZ[i]=E[ZB[i]];
            FuelleMatrix(TestM,NZB,Dim, Erz, TestMZ);
            if(RankMat0(TestM,NZB,Dim)>=Dim-2)
            {
                MakeNewSupp(Val,size,Tot,ValP, ValN);
                if(PMono)
                    goto PDone;
            }
            continue;
        }

        Einzig=true;
        for(NAll=0;NAll<NNonSimp;NAll++)
        {
              ValTest=NonSimp[NAll];
              if (ValTest==ValP || ValTest==ValN) continue; //
              i=0;
              while (i< NZB && !(ValTest->LV)->V[ZB[i]])
                  i++;
              if (i==NZB)
              {
                  Einzig=false;
                  break;
              }
        }
        if (Einzig)
        {
            MakeNewSupp(Val,size,Tot,ValP, ValN);
            if(PMono)
                goto PDone;
        }

      }  // for ValN

  PDone:   i=i;


  } // for ValP

  free(ZB);
  free(ZP);
  free(PosSimp);
  free(NegSimp);
  free(NonSimp);
  free(Pos);
  free(Neg);
  free(NReg);
  free(Neutr);
  free(Search);
  free(PVis);
  free(NVis);
  free(DVis);
  if(TestMSize)
    FreeMat(TestM, TestMSize);
  if(RangTest)
    free(TestMZ);

  ValP=Val->St;
  while (ValP)
  {
    ValN=ValP->R;
    if ((ValP->LV)->Mark < 0)
    {
        FreeLongVLE(Val, ValP);
    }
    ValP=ValN;
  }
}






Zahl SimplexSuppVal(LongVList *Values, Matrix Erz, int AnzErz, int Dim,
                       Indizes C)
// Bestimmt Stuetzhyperebenen und Werte vektoren fuer Gebrauch in TransformVal
// fuer simplizialen Kegel, dessen Erzeuger durch C aus Erz ausgewaehlt werden
// Rueckgabewert Vol=|Det|
{
  Matrix E, Supp;
  Zahl Vol;
  int i, j,ii;
  VLE *VP;
  LongVLE *LVP, *ValP;

  MatAlloc(&E,Dim,Dim);
  MatAlloc(&Supp,Dim,Dim);

  FuelleMatrix(E, Dim, Dim, Erz, C);
  Vol=SimpSupp(Supp, E, Dim);

  for (i=0;i<Dim;i++)
  {
      LVP = NewLongVLE(Values,AnzErz+Dim);
      for (j=0;j<Dim; j++)
      {
          (LVP->LV)->V[j]=Supp[i][j];
          (LVP->LV)->V[j+Dim]=SkalProd(E[j], Supp[i],Dim);
      }
  }
  FreeMat(E,Dim);
  FreeMat(Supp,Dim);
  return(Vol);
  // SchreibeLVL(stdout,*Values,2*Dim);
}


void AddSimpl(SimplexList *Tri, Matrix Erz, int Dim, LongVList Val,
int NewVert,
                int Tot, Indizes Col, Indizes E)
//
// Erweitert die Trriangulierung von C=RR_+(x_1,...,x_n) zu Traingulierung von RR_+(x_1,...,x_{n+1})
// wobei x_{n+1} aus Erz durch NewVert ausgewählt wird.
// Stuetzhyperebenen von C repraesentiert durch die Wertvektoren in Val (Werte auf den x_i)
// j=Col[i] gibt an, welche Komponente Val->LV->V[j] dem Vektor Erz[i] entspricht
{
   Simplex *Ende, *S, *Rette, *Next;
   int i, AnzP, NZ,j,AnzZ;
   LongVLE *ValP;
   Simplex *NewS;
   Zahl test;
   boolean Weiter;

   Ende=Tri->E;

   for(ValP=Val.St;ValP;ValP=ValP->R)
   {
      if ((ValP->LV)->Mark >= 0)
      {
        // printf("Pos\n");
        continue;
      }

      AnzZ=0;
      for(i=Dim;i<Tot;i++)
      {
         if(!ValP->LV->V[i])
            AnzZ++;
      }
      if(AnzZ==Dim-1) // simplizial
      {
          NewS=NewSimplexLE(Tri);
          j=0;
          for(i=Dim;i<Tot;i++)
            if(!ValP->LV->V[i])
            {
                NewS->Vert[j]=E[i];
                j++;
            }
          NewS->Vert[Dim-1]=NewVert;
          // printf("Simp %d\n",NinSimplL(*Tri));
          continue;
      }


      Weiter=true;
      for(S=Tri->St;Weiter;S=S->R)
      {
         Weiter=(S!=Ende);
          AnzP=0;
          for(i=0;i<Dim;i++)
          {
              test=(ValP->LV)->V[Col[S->Vert[i]]];
              if(test)
              {
                    AnzP++;
                    NZ=i;
              }
              if(AnzP>1) break;
          }
          if (AnzP>1)
              continue;
          else
          {
              NewS=NewSimplexLE(Tri);
              memcpy(NewS->Vert,S->Vert,Dim*sizeof(int));
              NewS->Vert[NZ]=NewVert;
          }
      }
      // printf("NonSimp %d\n",NinSimplL(*Tri));
   }
}

void TriangSupp(SimplexList *Tri, VList *Supp, Matrix Erz, int AnzErz, int Dim, Indizes Extr)
//
// Bestimmt Triangulierung Tri und Stuetzhyperebenen Supp zu dem von Erz erzeugten Kegel
// Falls Extr[i] != 0, ist bekannt, dass Erz[i] ein Extremstrahl ist. Nach Moeglichkeit
// werden Extremstrahlen zuerst eingebaut, damit moeglichst wenige ecken in der Triangulierung
// verwendet werden
// Zu Col siehe oben, ebenso E in transformVakues
{
   int i,size,kk;
   Indizes C,E;
   LongVLE *ValP;
   Indizes Col;
   boolean ReallyNew, *InTriang;
   Zahl L;
   Simplex *NewS;
   LongVList Val;

   Matrix Sch;

   C= (Indizes) calloc(Dim,sizeof(int));
   E= (Indizes) calloc(Dim+AnzErz,sizeof(int));
   InTriang=(boolean*) calloc(AnzErz,sizeof(boolean));
   Col=(Indizes) calloc(AnzErz,sizeof(int));

   // In Col[i] steht derjenige Index j bzgl. Val, der den wert
   // an der Stelle i liefert.

   for (i=0;i<AnzErz;i++)
       InTriang[i]=false;

   FindeInitSimpl0(C,Erz,AnzErz,Dim,0);
   for(i=0;i<Dim;i++)
   {
        InTriang[C[i]]=true;
        Col[C[i]]=i+Dim;
        E[Dim+i]=C[i];
   }

   Tri->St=Tri->E=NULL;
   Val.St=NULL;
   Val.E=NULL;

   NewS=NewSimplexLE(Tri);
   memcpy(NewS->Vert,C,Dim*sizeof(int));
   NewS->Vol=SimplexSuppVal(&Val,Erz,AnzErz,Dim,C);

   size=2*Dim; // Bisherige effektive Laenge der Val-Vektoren


for(kk=0;kk<=1;kk++) // 2 Durchgaenge, im ersten moeglichst nur Extremstrahlen
   for(i=0;i<AnzErz;i++)
   {

      if(InTriang[i])
            continue;
      if(!kk && !Extr[i])
        continue;

      ReallyNew=false;


      for(ValP=Val.St;ValP;ValP=ValP->R)
      {
          L=SkalProd((ValP->LV)->V, Erz[i],Dim);
          if (L<0) ReallyNew=true;
          (ValP->LV)->Mark=L;
      }

      if(ReallyNew)
      {
        InTriang[i]=true;
        AddSimpl(Tri,Erz,Dim,Val,i,size,Col,E);
        // AddSimpl(Tri,Erz,Dim,Val,i,Col);
        E[size]=i;
        TransformValues(&Val,Dim,size,Erz,E);
        Col[i]=size;
        size++; // Val-Vektoren sind um ein Feld laenger geworden

      }

  }


  LVL2VL(Supp,Val,Dim);  // Supp-Formen aus Val herausziehen

  FreeLongVL(&Val);
  free(Col);
  free(C);
  free(E);
  free(InTriang);
}

void TriangSuppVL(SimplexList *Tri, VList *Supp, VList Gen, int Dim)
// Das Gleiche, nur von VL ausgehend
{
    Matrix Erz;
    int AnzErz;
    Indizes Extr;

    AnzErz=NinVL(Gen);
    MatAlloc(&Erz,AnzErz,Dim);
    VL2Mat(Erz,Gen,Dim);
    Extr=calloc(AnzErz,sizeof(int));
    TriangSupp(Tri,Supp,Erz,AnzErz,Dim,Extr);
    FreeMat(Erz,AnzErz);
    free(Extr);
}


void SuppHyp(VList *Supp, Matrix Erz, int AnzErz, int Dim, Indizes Extr)
//
// Bestimmt  Stuetzhyperebenen Supp zu dem von Erz erzeugten Kegel
// Falls Extr[i] != 0, ist bekannt, dass Erz[i] ein Extremstrahl ist. Nach Moeglichkeit
// werden Extremstrahlen zuerst eingebaut
// Zu C und E siehe oben
{
   int i,size,kk;
   Indizes C, E;
   LongVLE *ValP;
   boolean ReallyNew, *InTriang;;
   Zahl L;
   Simplex *NewS;
   LongVList Val;

   // Matrix Sch;

   InTriang=(boolean*) calloc(AnzErz,sizeof(boolean));
   C= (Indizes) calloc(Dim,sizeof(int));
   E= (Indizes) calloc(Dim+AnzErz,sizeof(int));

   for (i=0;i<AnzErz;i++)  // Hiermerken wir uns, ob Vektor i schon verarbeitet
       InTriang[i]=false;

   FindeInitSimpl0(C,Erz,AnzErz,Dim,0);
   for(i=0;i<Dim;i++)
   {
        InTriang[C[i]]=true;
        E[Dim+i]=C[i];
   }

   Val.St=NULL;
   Val.E=NULL;

   SimplexSuppVal(&Val,Erz,AnzErz,Dim,C); // Stuetzhyperebenen zu bisherigem Simplex

   size=2*Dim; // Bisherige effektive Laenge der Val-Vektoren


for(kk=0;kk<=1;kk++) // 2 Durchgaenge, im ersten moeglichst nur Extremstrahlen
   for(i=0;i<AnzErz;i++)
   {

      // printf("Starting generator %d\n",i);
      if(InTriang[i])
            continue;
      if(!kk && !Extr[i])
        continue;

      ReallyNew=false;


      for(ValP=Val.St;ValP;ValP=ValP->R)
      {
          L=SkalProd((ValP->LV)->V, Erz[i],Dim);
          if (L<0) ReallyNew=true;
          (ValP->LV)->Mark=L;
      }

      if(ReallyNew)
      {
        InTriang[i]=true;
        E[size]=i;
        TransformValues(&Val,Dim,size,Erz,E);
        size++; // Val-Vektoren sind um ein Feld laenger geworden
      }
  }

  LVL2VL(Supp,Val,Dim);  // Supp-Formen aus Val herausziehen

  FreeLongVL(&Val);
  free(C);
  free(E);
  free(InTriang);
}

void SuppHypVL(VList *Supp, VList Gen, int Dim)
// Das Gleiche, nur von VL ausgehend
{
    Matrix Erz;
    int AnzErz;
    Indizes Extr;

    AnzErz=NinVL(Gen);
    MatAlloc(&Erz,AnzErz,Dim);
    VL2Mat(Erz,Gen,Dim);
    Extr=calloc(AnzErz,sizeof(int));
    SuppHyp(Supp,Erz,AnzErz,Dim,Extr);
    FreeMat(Erz,AnzErz);
    free(Extr);
}



// Hilbert-Basen und Gitterpunkte

void HilbBasTriOld(VList ConeGen, int EmbDim, VList *HilbBas,
                  VList *Supp, int *Rang, Zahl *Index)
// Bestimmt Hilbert-Basis und Stuetzhyperebenen zu dem von
// ConeGen erzeugten Kegel
// Tri: normaliz-Algorithmus
// Setzt voraus, dass ConeGen volldimensionalen spitzen Kegel erzeugen
// Sonst nur Rueckgabe von Rang des von ConeGen erzeugten Z-Moduls
// Index (im volldimenionalen Fall) ist der Index dieses Untermoduls
// im vollen Gitter
// alte Version
{
    int NConeGen,i,Dim,NSupp;
    Matrix Erz,SpxErz;
    SimplexList Tri;
    Simplex *Spx;
    VList HBPar;
    VLE *VP;
    Indizes Extr;

    StartVL(HilbBas);
    StartVL(Supp);


    RangIndexVL(ConeGen,EmbDim,Rang,Index);
    Dim=*Rang;
    if(Dim<EmbDim)
        return;

    NConeGen=NinVL(ConeGen);
    MatAlloc(&Erz,NConeGen,Dim);
    VL2Mat(Erz,ConeGen,Dim);
    MatAlloc(&SpxErz,Dim,Dim);
    Extr=(Indizes) calloc(NConeGen,sizeof(int));

    TriangSupp(&Tri,Supp,Erz,NConeGen,Dim,Extr);
    NSupp=NinVL(*Supp);

    for(i=0;i<NConeGen;i++)
        RedInsert(HilbBas,Dim,Erz[i],*Supp,NSupp);

    for(Spx=Tri.St;Spx;Spx=Spx->R)
    {
            // if(Spx->Vol==1) continue;

            FuelleMatrix(SpxErz,Dim,Dim,Erz,Spx->Vert);
            HilbParOld(&HBPar,SpxErz,Dim);

            for(VP=HBPar.St;VP;VP=VP->R)
                RedInsert(HilbBas,Dim,VP->V,*Supp,NSupp);

            FreeVL(&HBPar);
    }

    FreeMat(Erz,NConeGen);
    FreeMat(SpxErz,Dim);
    FreeSimplexL(&Tri);
    free(Extr);
}

void HilbBasTri(VList ConeGen, int EmbDim, VList *HilbBas,
                  VList *Supp, int *Rang, Zahl *Index)
// Bestimmt Hilbert-Basis und Stuetzhyperebenen zu dem von
// ConeGen erzeugten Kegel
// Tri: normaliz-Algorithmus
// Setzt voraus, dass ConeGen volldimensionalen spitzen Kegel erzeugen
// Sonst nur Rueckgabe von Rang des von ConeGen erzeugten Z-Moduls
// Index (im volldimenionalen Fall) ist der Index dieses Untermoduls
// im vollen Gitter
{
    int NConeGen,i,Dim,NSupp,NKand,j,Gr;
    Matrix Erz,SpxErz,SuppMat;
    SimplexList Tri;
    Simplex *Spx;
    VList HBPar;
    VLE *VP,*VP1;
    Indizes Extr;
    Vektor ValV;
    VList FirstList;
    Zahl S,Tiefe,ZHilb, TGrad,Comp0;
    boolean Reduziert;



    KandV *Kand;

    StartVL(HilbBas);
    StartVL(Supp);


    RangIndexVL(ConeGen,EmbDim,Rang,Index);
    Dim=*Rang;
    if(Dim<EmbDim)
        return;

    NConeGen=NinVL(ConeGen);
    MatAlloc(&Erz,NConeGen,Dim);
    VL2Mat(Erz,ConeGen,Dim);
    MatAlloc(&SpxErz,Dim,Dim);
    Extr=(Indizes) calloc(NConeGen,sizeof(int));

    TriangSupp(&Tri,Supp,Erz,NConeGen,Dim,Extr);
    NSupp=NinVL(*Supp);

    // printf("%d Simpl\n",NinSimplL(Tri));

    ValV= (Vektor) calloc(NSupp,sizeof(Zahl));

    MatAlloc(&SuppMat,NSupp,Dim);
    VL2Mat(SuppMat,*Supp,Dim);

    CopyVL(&FirstList,ConeGen,Dim);

    for(Spx=Tri.St;Spx;Spx=Spx->R)
    {
            // if(Spx->Vol==1) continue;

            // for(i=0;i<Dim;i++) printf("%d ",Spx->Vert[i]); printf("\n");

            FuelleMatrix(SpxErz,Dim,Dim,Erz,Spx->Vert);
            HilbPar(&HBPar,SpxErz,Dim);
            ConcatVLFree(&FirstList,&HBPar);
    }
    NKand=NinVL(FirstList);

    Kand = (KandV*) calloc(NKand,sizeof(KandV));
    VP=FirstList.St;
    for(i=0;i<NKand;i++)
    {
        Kand[i].V=VP->V;
        MxV(ValV,SuppMat,NSupp,Dim,Kand[i].V);
        S=0;
        for(j=0;j<NSupp;j++)
            S+=ValV[j];
        Kand[i].TotGrad=S;  // Totalgrad ist 0-te Komponente
        VP=VP->R;
    }

    // printf("%d Kandidaten unsortiert\n",NKand);

    VEKTLAENGE=Dim;
    qsort(Kand,NKand,sizeof(KandV),VecVergl);

    // printf("%d Kandidaten sortiert\n",NKand);

    Tiefe=0;
    ZHilb=0;

    Comp0=0;


    for(i=0;i<NKand;i++)
    {
        if( i<NKand-1 && !memcmp(Kand[i].V,Kand[i+1].V,Dim*sizeof(Zahl)))
            continue;

        MxV(ValV,SuppMat,NSupp,Dim,Kand[i].V);
        TGrad=Kand[i].TotGrad; // TotalGrad übernehmen
        // printf("TG %lld ",TGrad);

        Reduziert=false;
        Gr=0;


        for(VP=HilbBas->St;VP;VP=VP->R)
        {
            ZHilb++;

            if(2*VP->Mark>TGrad)
                break;

            if(ValV[Gr]<VP->LV[Gr])
               continue;

            Comp0++;

            for(j=0;j<NSupp;j++)
            {
                if(ValV[j]<VP->LV[j])
                    break;
            }
            Gr=j;
            Tiefe+=j+1;
            if(j==NSupp) // reduziert
            {
                Reduziert=true;
                // printf("Red\n");
                VP->Mark1++;
                MakeFirst(HilbBas,VP);
                break;
            }
        }
        if(Reduziert)
            continue;


    InHilbBas:
        VP1=NewVLE(HilbBas); // jetzt irreduzibel
        // printf("HB\n");

        memcpy(VP1->V,Kand[i].V,Dim*sizeof(Zahl));
        AllocLVE(VP1,NSupp);
        memcpy(VP1->LV,ValV,NSupp*sizeof(Zahl));
        VP1->Mark=TGrad;
        VP1->Mark1=0;
    }

    /* for(VP=HilbBas->St;VP;VP=VP->R)
    {
        printf("%lld\n",VP->Mark1);
    } */

    // printf("%lld total Vergl, %lld Vergl Start bei 0\n",ZHilb,Comp0);

    // printf("Im Schnitt Red nach %f VTiefe %f\n",(float) ZHilb/NKand,(float) Tiefe/ZHilb);

    free(Kand);
    FreeVL(&FirstList);
    FreeMat(SuppMat,NSupp);
    FreeMat(Erz,NConeGen);
    FreeMat(SpxErz,Dim);
    FreeSimplexL(&Tri);
    free(Extr);
    free(ValV);
}

void HilbBasDualTri(VList Supp,int Dim,VList *HilbBas)
// Finde Hilbertbasis ausgehend von Supp
// setzt implizit voraus, dass Kegel volldimensional und spitz
{
    Zahl Index;
    int Rang;
    VList Supp1, ConeGen;

    SuppHypVL(&ConeGen,Supp,Dim); // Finde Erz von Kegel
    HilbBasTri(ConeGen, Dim, HilbBas,&Supp1,&Rang,&Index);
    FreeVL(&ConeGen);
    FreeVL(&Supp1);
}

void HilbSimp(VList *HilbBas, Matrix SimpGen, int Dim)
// Macht Hilbert_basis des simplizialen Kegels erzeugt von SimpGen
// setzt volle Dimension voraus und nimmt an, dass SimpGen zur Hilbert-basis
// gehoeren
{
    VList SGen;

    HilbPar(HilbBas,SimpGen,Dim);
    Mat2VL(&SGen,SimpGen,Dim,Dim);
    ConcatVLFree(HilbBas,&SGen);
}

void HilbSimpDual(VList *HilbBas, Matrix Supp, int Dim)
// Das gleiche, nur ausgehend von Stuetzhyperebenen
{
    Matrix SimpGen;

    MatAlloc(&SimpGen,Dim,Dim);

    SimpSupp(SimpGen,Supp,Dim); // find generators of simplicial overcone
    HilbSimp(HilbBas,SimpGen,Dim);

    FreeMat(SimpGen,Dim);
}

void LatPinSimpl(VList *LatPInt, VList *LatPBd, Matrix U, int Dim)
// Bestimmt Gitterpunkte in Simplex gegeben durch U (auf Hoehe 1)
// Innere Punkte in LatPInt, Randpunkte in LatPBd zurueckgegeben
{
    int i,Last,Dim1;
    Vektor C,V,Diag;
    Zahl Nenner;
    Matrix Inv;
    boolean Bd;

    Dim1=Dim-1;
    MatAlloc(&Inv,Dim,Dim);
    Diag= calloc(Dim,sizeof(Zahl));
    C=(Vektor) calloc(Dim,sizeof(Zahl));
    V=(Vektor) calloc(Dim,sizeof(Zahl));

    InvertZ(Inv, U, Dim, Diag, &Nenner);

    StartVL(LatPInt);
    StartVL(LatPBd);

    while(1)
    {
        Last = -1;
        for (i=0;i<Dim;i++) if (C[i]<Diag[i]-1) Last = i;
        if (Last == -1) break;
        C[Last]++;
        for(i=Last+1;i<Dim;i++) C[i]=0;

        MxV(V,Inv,Dim,Dim,C);
        Bd=false;
        for(i=0;i<Dim;i++)
        {
            V[i] %= Nenner;
            if(!V[i])
            {
                Bd=true;
                continue;
            }
            if(V[i]<0)
                V[i]+=Nenner;
        }

        VxM(V,V,Dim,U,Dim);
        for(i=0;i<Dim;i++)
            V[i]/=Nenner;
        if(V[Dim1]!=1)
            continue;
        if(Bd)
            memcpy(NewVLE(LatPBd)->V,V,Dim*sizeof(Zahl));
        else
            memcpy(NewVLE(LatPInt)->V,V,Dim*sizeof(Zahl));
    }

    FreeMat(Inv,Dim);
    free(Diag);
    free(C);
    free(V);

}

void FindLatPoints(VList ConeGen, int EmbDim, VList *LatPoints, VList *Supp)
//  Bestimmt Gitterpunkte und Stuetzhyperebenen zu Polytop gegeben durch ConeGen (auf Hoehe 1)
// setzt voraus, dass Kegel volldimensional und spitz ist
{
    int NConeGen,i,Dim,Dim1;
    Matrix Erz,SpxErz;
    SimplexList Tri;
    Simplex *Spx;
    VList LatParInt, LatParBd, LatBd;
    VLE *VP;
    Indizes Extr;

    StartVL(LatPoints);
    StartVL(&LatBd);  // nimmt zunaechst die Punkte in Raendern von Simplizes auf
    StartVL(Supp);
    Dim=EmbDim;
    Dim1=Dim-1;

    NConeGen=NinVL(ConeGen);
    MatAlloc(&Erz,NConeGen,Dim);
    VL2Mat(Erz,ConeGen,Dim);
    MatAlloc(&SpxErz,Dim,Dim);
    Extr=(Indizes) calloc(NConeGen,sizeof(int));

    TriangSupp(&Tri,Supp,Erz,NConeGen,Dim,Extr);

    for(i=0;i<NConeGen;i++)
    {
        if(Erz[i][Dim1]==1)
            memcpy(NewVLE(LatPoints)->V,Erz[i],Dim*sizeof(Zahl));
    }

    for(Spx=Tri.St;Spx;Spx=Spx->R)
    {
            FuelleMatrix(SpxErz,Dim,Dim,Erz,Spx->Vert);
            LatPinSimpl(&LatParInt,&LatParBd,SpxErz,Dim);
            /* printf("Simplex\n");
            FSchreibeVL(stdout,LatParInt,Dim);
            FSchreibeVL(stdout,LatParBd,Dim); */
            ConcatVLFree(LatPoints,&LatParInt);
            for(VP=LatParBd.St;VP;VP=VP->R)
                InsertV(&LatBd,VP->V,Dim); // um Doppeleintraege zu vermiden
            FreeVL(&LatParBd);
    }
    ConcatVLFree(LatPoints,&LatBd); // jetzt zusammenfuegen
    FreeMat(Erz,NConeGen);
    FreeMat(SpxErz,Dim);
    FreeSimplexL(&Tri);
    free(Extr);
}


// Extremstrahlen, innere Punkte und Homogenitaet

void ExtremeR(VList Gen,int Dim, int NSupp)
//
// Bestimmt Extremstrahlem des von Gen erzeugten Kegels
// Setzt voraus, dass die Erzeuger Gen mit Wert-Vektoren ausgestattet sind
// Deren Laenge ist NSupp
// Anzeige ob extrem, ueber Mark
{
    int i,k,l;
    VLE *VP, *VP1;
    Indizes ZP;

    ZP=(Indizes) calloc(NSupp,sizeof(int));

    for(VP=Gen.St;VP;VP=VP->R)
        VP->Mark=1;

    for(VP=Gen.St;VP;VP=VP->R)
    {
        k=0;
        for(i=0;i<NSupp;i++)
            if(!VP->LV[i])
            {
                ZP[k]=i;
                k++;
            }
        if(k<Dim-1) // in nicht genuegend Facetten
        {
            VP->Mark=0;
            continue;
        }

        for(VP1=Gen.St;VP1;VP1=VP1->R)
        {
            if(VP1==VP) continue; // nicht mit sich selbst vergleichen
            if(!VP1->Mark) continue; // nicht mit bekannt nicht-extremen vergleichen
            l=0;
            for(i=0;i<k;i++)
                if(!VP1->LV[ZP[i]]) l++;
            if(l>=k)
            {
                VP->Mark=0;
                break;
            }
        }
    }

    free(ZP);
}

void SelExtremeR0(VList *Extr, VList Gen, VList Supp, int Dim)
// Waehlt Extremstrahlen aus, macht Werte in Gen, setzt Marks in Gen!
{
    VLE *VP;

    AddValLV(Gen,Supp,Dim);
    ExtremeR(Gen,Dim,NinVL(Supp));
    StartVL(Extr);
    for(VP=Gen.St;VP;VP=VP->R)
        if(VP->Mark)
            CopyVLE0(Extr,VP,Dim);
}

void SelExtremeR(VList *Extr, VList Gen, VList Supp, int Dim)
// Waehlt Extremstrahlen aus, veraendert Gen nicht
{
    VLE *VP, *VPNext;

    CopyVL(Extr,Gen,Dim);
    AddValLV(*Extr,Supp,Dim);
    ExtremeR(*Extr,Dim,NinVL(Supp));
    StripLV(*Extr);
    for(VP=Extr->St;VP;VP=VPNext)
    {
        VPNext=VP->R;
        if(!VP->Mark)
            FreeVLE(Extr,VP);
    }
}

boolean TesteHom(Vektor Hom, VList Gen, int Dim)
// Ermittelt, ob die Vektoren Gen in einer Hyperebene mit Abstand 1 vom Nullpunkt liegen
// Diese wird gegebenenfalls ueber Hom zurueckgegeben
// Setzt voraus, dass Gen volldimensional ist
{
    int AnzGen, i;
    Zahl Nenner;
    boolean ExistUniq;
    Vektor RS;
    Matrix GLS;

    AnzGen=NinVL(Gen);
    MatAlloc(&GLS,AnzGen,Dim);
    VL2Mat(GLS,Gen,Dim);
    RS=(Vektor) calloc(AnzGen,sizeof(Zahl));
    for(i=0;i<AnzGen;i++)
        RS[i]=1;
    Loese(Hom,GLS,AnzGen,Dim,RS,&Nenner,&ExistUniq);
    FreeMat(GLS,AnzGen);
    free(RS);
    if(!ExistUniq)
        return(false);
    for(i=0;i<Dim;i++)
        if(Hom[i] % Nenner)
            break;
    if(i<Dim)
        return(false);
    for(i=0;i<Dim;i++)
        Hom[i]/=Nenner;
    return(true);
}

boolean InInt(Vektor Vek, VList Supp, int Dim)
// Prueft ob Vek im Innern des Kegels mit Stuetzhyp Supp liegt
{
    VLE *VP;
    for(VP=Supp.St;VP;VP=VP->R)
        if(SkalProd(VP->V,Vek,Dim)<=0)
            return false;
    return true;
}

boolean InCone(Vektor Vek, VList Supp, int Dim)
// Prueft ob Vek im Kegel mit Stuetzhyp Supp liegt
{
    VLE *VP;
    for(VP=Supp.St;VP;VP=VP->R)
        if(SkalProd(VP->V,Vek,Dim)<0)
            return false;
    return true;
}

// Gitter-Polytope


boolean IsFlat(VList Gen, int Dim)
// Prueft, ob alle Elemente von Gen auf Hoehe 1
{
    VLE *VP;
    for(VP=Gen.St;VP;VP=VP->R)
        if(VP->V[Dim-1]!=1)
            return false;
    return true;
}

void CornerCone(VList *Corner,VList LatPoints,VList *CornerSupp, VList Supp,
                            Vektor V, int EmbDim)
// Mach den Eckenkegel in der Ecke V, inklusive Stuetzhyperebenen
// Falls nur CornerSupp bestimmt werden soll, CornerHyp benutzen
{

    int i;
    VLE *VP,*VP1;

    StartVL(Corner);

    if(V[EmbDim-1]!=1)
    {
        printf("Corner not flat 1\n");
        exit(1);
    }
    for(VP=LatPoints.St;VP;VP=VP->R)
    {
        if(VP->V[EmbDim-1]!=1)
        {
            printf("Corner not flat 2\n");
            exit(1);
        }
        VP1=NewVLE(Corner);
        for(i=0;i<EmbDim-1;i++)
            VP1->V[i]=VP->V[i]-V[i];
        if(NullVekt(VP1->V,EmbDim-1))
            FreeVLE(Corner,VP1);
    }

    StartVL(CornerSupp);
    for(VP=Supp.St;VP;VP=VP->R)
    {
        if(!SkalProd(VP->V,V,EmbDim))
            CopyVLE0(CornerSupp,VP,EmbDim-1);
    }
}

void CornerHyp(VList *CornerSupp, VList Supp, Vektor V, int EmbDim)
// Macht Stuetzhyperebenen des Eckenkegels in der Ecke V
{

    VLE *VP;

    StartVL(CornerSupp);
    for(VP=Supp.St;VP;VP=VP->R)
    {
        if(!SkalProd(VP->V,V,EmbDim))
            CopyVLE0(CornerSupp,VP,EmbDim-1);
    }
}

boolean NormalCorner(Vektor C, VList Corner, VList CornerSupp, VList HilbBas, int Dim)
// Prueft ob die Gitterpunkte im Polytop Hilbert-Basis an der Ecke C definieren
// HilbBas ist dabei globale Hilbert-Basis
{
    VLE *VP,*VP1;
    int i,NCS;
    Zahl Deg;
    Vektor V,VC,ValVC;


    AddValLV(Corner,CornerSupp,Dim-1);
    NCS=NinVL(CornerSupp);
    VC=calloc(Dim-1,sizeof(Zahl));
    ValVC=calloc(NCS,sizeof(Zahl));

    for(VP=HilbBas.St;VP;VP=VP->R)
    {
        V=VP->V; //Abkuerzung
        Deg=V[Dim-1];
        if(Deg==1)
            continue;

        for(i=0;i<Dim-1;i++)
            VC[i]=V[i]-Deg*C[i];

        i=0;
        for(VP1=CornerSupp.St;VP1;VP1=VP1->R)
        {
            ValVC[i]=SkalProd(VC,VP1->V,Dim-1);
            i++;
        }
        if(Irred(Corner,ValVC,NCS))
        {
            StripLV(Corner);
            free(VC);
            free(ValVC);
            return false;
        }
    }

    StripLV(Corner);
    free(VC);
    free(ValVC);
    return true;
}

void VASS(VList LatPoints, VList Supp, VList HilbBas, int EmbDim, boolean Normal,
                    boolean *VA, boolean *Simple,boolean *Smooth)
// Prueft auf Very Ample, Simple, Smooth
// HilbBas globale Hilbert-Basis, Supp globale Stuetzhyperebenen
{
    VList Corner,CornerSupp,Extr;
    Zahl Index;
    // boolean Destr;
    int Rang, Check;

    VLE *VP;

    // Check=VektorZaehler-NinVL(VLFree);

    // printf("Starte VeryAmple\n");

    *VA=*Simple=*Smooth=true;

    StartVL(&Extr);
    StartVL(&Corner);
    StartVL(&CornerSupp);


    SelExtremeR(&Extr, LatPoints, Supp, EmbDim);
    for(VP=Extr.St;VP;VP=VP->R)
    {
        // printf("%d %d %d\n",*VA,*Simple,*Smooth);
        if(!(*VA) && !(*Simple))
            break;

        CornerCone(&Corner,LatPoints,&CornerSupp,Supp,VP->V,EmbDim);
        if(*VA && !Normal)
        {
            *VA=NormalCorner(VP->V,Corner,CornerSupp, HilbBas, EmbDim);
            if(!*VA)
                *Smooth=false;
        }

        if(!*Simple)
        {
            FreeVL(&Corner);
            FreeVL(&CornerSupp);
            continue;
        }

        if(NinVL(CornerSupp)!=EmbDim-1)
        {
            *Simple=false;
            *Smooth=false;
        }

        if(!*Smooth)
        {
            FreeVL(&CornerSupp);
            FreeVL(&Corner);
            continue;
        }


        RangIndexVL(CornerSupp,EmbDim-1, &Rang, &Index);
        if(Index!=1)
            *Smooth=false;

        FreeVL(&CornerSupp);
        FreeVL(&Corner);
    }

    FreeVL(&Extr);

    // printf("Heraus\n");

    // if(VektorZaehler-NinVL(VLFree)!=Check)
    //               printf("Alarm !!!\n!");

}

void HomothIm(VList *Image, VList Gen, int Dim, Vektor Center, Zahl fact)
// Macht Bild von Gen unter Streckung mit Faktor fact und Zentrum Center
{
    VLE *VP,*VP1;
    Zahl fact1;
    int i;

    StartVL(Image);

    fact1=fact-1;
    for(VP=Gen.St;VP;VP=VP->R)
    {
        VP1=NewVLE(Image);
        for(i=0;i<Dim;i++)
            VP1->V[i]=fact*VP->V[i]-fact1*Center[i];
    }
}


void MakePolPos(VList Gens,int Dim)
// Mache alle Ecken Gens nichtnegativ durch Parallelverschiebung
{
    VLE *VP;
    int i;
    Zahl Min;

    for(VP=Gens.St;VP;VP=VP->R)
        if(VP->V[Dim-1]!=1)
        {
            printf("MakePolPos not flat\n");
            exit(1);
        }

    for(i=0;i<Dim-1;i++)
    {
        VP=Gens.St;
        Min=VP->V[i];
        for(VP=VP->R;VP;VP=VP->R)
            if(VP->V[i]<Min)
                Min=VP->V[i];
        for(VP=Gens.St;VP;VP=VP->R)
            VP->V[i]-=Min;
    }
}

Zahl VolPolyt(VList ConeGen, int EmbDim,VList *Supp, SimplexList *Tri)
// Bestimmt Volumen und Stuetzhyperebenen und Triangulierung
// des von ConeGen erzeugten Polytops
// Setzt voraus, dass Polytope volldimensional, sonst Rueckgabewert 0
{
    int NConeGen,i,Dim;
    Matrix Erz,SpxErz;
    Simplex *Spx;
    Indizes Extr;
    Zahl Index, Volu;


    StartVL(Supp);


    RangIndexVL(ConeGen,EmbDim,&Dim,&Index);
    if(Dim<EmbDim)
        return 0;

    NConeGen=NinVL(ConeGen);
    MatAlloc(&Erz,NConeGen,Dim);
    VL2Mat(Erz,ConeGen,Dim);
    MatAlloc(&SpxErz,Dim,Dim);
    Extr=(Indizes) calloc(NConeGen,sizeof(int));

    TriangSupp(Tri,Supp,Erz,NConeGen,Dim,Extr);

    Volu=0;

    for(Spx=Tri->St;Spx;Spx=Spx->R)
    {
            FuelleMatrix(SpxErz,Dim,Dim,Erz,Spx->Vert);
            Volu+=Vol(SpxErz,Dim);
    }

    FreeMat(Erz,NConeGen);
    FreeMat(SpxErz,Dim);
    free(Extr);

    return Volu;
}



// Stellare Unterteilung und unimodulare Triangulierung

void StellarSubDiv(VList *OneSkel, SimplexList *Triang, Vektor V, int Dim)
// Macht stellare Unterteilung und gibt OneSkel und Triang veraendert zurueck
// Es wird davon ausgegangen, dass Triang mit Corn und Supp versehen ist.
{
    Simplex *Smp,*ThisSmp,*NextSmp,*EndTriang, *NewSmp;
    VLE *VP;
    int i,j,NewVert;
    Indizes Exchange;
    boolean InCone;
    Zahl Test;

    if(!Triang->St) return;

    Exchange=calloc(Dim,sizeof(int));

    EndTriang=Triang->E;
    VP=NewVLE(OneSkel);
    memcpy(VP->V,V,Dim*sizeof(Zahl));
    NewVert=NinVL(*OneSkel)-1;

    // printf("Hurra %d\n",NewVert);

    Smp=Triang->St;
    do
    {
        ThisSmp=Smp;
        NextSmp=Smp->R;
        InCone=true;
        i=-1;
        j=-1;
        for(VP=Smp->Supp.St;VP;VP=VP->R)
        {
            i++;
            Test=SkalProd(VP->V,V,Dim);
            if(Test<0)
            {
                InCone=false;
                break;
            }
            if(Test>0)
            {
                j++;
                Exchange[j]=i;
            }
        }
        if(InCone)
        {
            for(i=0;i<=j;i++)
            {
                NewSmp=NewSimplexLE(Triang);
                memcpy(NewSmp->Vert,Smp->Vert,Dim*sizeof(int));
                NewSmp->Vert[Exchange[i]]=NewVert;
                qsort(NewSmp->Vert,Dim,sizeof(int),intvergl);
                SetCornSuppE(*OneSkel,NewSmp,Dim);
            }
            FreeSimplexLE(Triang,Smp);
        }
        Smp=NextSmp;
    }
    while(ThisSmp!=EndTriang && j!=Dim-1);

    free(Exchange);
}

void RandomOrder(VList *VL)
{
    int N,i;
    VLE *VP;

    N=NinVL(*VL);
    for(i=0;i<2*N;i++)
    {
        VP=VLEind(*VL,rand() % N);
        MakeFirst(VL,VP);
    }
}

void UnimodTriang(VList *NewOneSkel, SimplexList *NewTriang,
                        VList OneSkel, SimplexList Triang, int Dim)
// Verfeinert Triangulierung (OneSkel,Triang) zu einer unimodularen
// (NewOneSkel,NewTriang)
// ACHTUNG: Setzt voraus, dass Corn und Supp gesetzt sind
{
    VList AccuHilb, InPar;
    VLE *VP,*VPNext;
    Simplex *Smp;
    Matrix U;

    MatAlloc(&U,Dim,Dim);

    CopySimplexList(NewTriang,Triang,Dim);
    CopyVL(NewOneSkel,OneSkel,Dim);
    StartVL(&AccuHilb);

    do
    {
        FreeVL(&AccuHilb);

        for(Smp=NewTriang->St;Smp;Smp=Smp->R)
        {
            if(Smp->Vol==1)
                continue;
            VL2Mat(U,Smp->Corn,Dim);
            HilbPar(&InPar,U,Dim);
            // ConcatVLFree(&AccuHilb,&InPar);
            for(VP=InPar.St;VP;VP=VPNext)
            {
                VPNext=VP->R;
                if(!VLEinVL(VP,AccuHilb,Dim))
                    MoveVLE(&AccuHilb,&InPar,VP);
            }
            FreeVL(&InPar);
        }

        /* printf("Hier %d\n",NinVL(AccuHilb));
        SchreibeVL(stdout,AccuHilb,Dim);
        printf("================\n");*/

        RandomOrder(&AccuHilb);

        for(VP=AccuHilb.St;VP;VP=VP->R)
            StellarSubDiv(NewOneSkel, NewTriang, VP->V,Dim);
    }
    while(AccuHilb.St);

    FreeMat(U,Dim);
}

void ResolveNormalFan(VList *NewSupp,SimplexList *NormalFan, VList Extr, VList Supp,int Dim)
// ACHTUNG: Dim ist Kegeldimension von Extr
{
    SimplexList Triang, TriE;
    Simplex *Simp,*SimpE;
    VList SuppE,Dummy;
    VLE *VP,*VP1,*VP2;

    int Dim1,i,j;
    Indizes Reg;

    // SchreibeVL(stdout,Supp,Dim);
    //    printf("==========\n");

    StartSimplexL(&Triang);
    StartVL(&SuppE);
    Dim1=Dim-1;
    Reg=calloc(NinVL(Supp),sizeof(int));

    for(VP=Extr.St;VP;VP=VP->R)
    {
        j=-1;
        i=-1;
        for(VP1=Supp.St;VP1;VP1=VP1->R) // waehle Supp von Eckenkegel aus
        {
            j++;
            if(!SkalProd(VP->V,VP1->V,Dim))
            {
                i++;
                Reg[i]=j;
                VP2=NewVLE(&SuppE);
                memcpy(VP2->V,VP1->V,Dim1*sizeof(Zahl));
            }
        }

        //SchreibeVL(stdout,SuppE,Dim1);
        // printf("----------\n");
        TriangSuppVL(&TriE, &Dummy, SuppE, Dim1); // Triangulierung des Duals von Eckenkegel
        FreeVL(&Dummy);
        FreeVL(&SuppE);

        for(SimpE=TriE.St;SimpE;SimpE=SimpE->R)
        {
            Simp=NewSimplexLE(&Triang);
            for(i=0;i<Dim1;i++)               // Rueckuebersetzung
            {
                Simp->Vert[i]=Reg[SimpE->Vert[i]];
            }
        }
        FreeSimplexL(&TriE);

    }  // jetzt haben wir simplizialen Faecher

    // SchreibeVL(stdout,Extr,Dim);

    // SchreibeTriang(stdout,Supp,Dim1,Triang);

    SetCornSupp(Supp,Triang,Dim1);

    UnimodTriang(NewSupp, NormalFan, Supp,Triang,Dim1); // unimodular unterteilen
    FreeSimplexL(&Triang);
    free(Reg);

}


// System-Interface

// int CallZ=0;

void run_pgm(char *pgm, char *par1, char *par2, char *par3)
// Baut einen Befehl auf, der dann extern ausgefuehrt wird
// Setzt sich zusammen aus den 4 Parametern, wobei hinter pgm
// ein Blank eingefuegt wird
{
    char *FullCommand;

    fflush(stdout);
    FullCommand= (char*) calloc (strlen(pgm)+1+strlen(par1)+
          strlen(par2)+strlen(par3)+3,sizeof(char));
    strcpy(FullCommand,pgm);
    strcat(FullCommand," ");
    strcat(FullCommand,par1);
    strcat(FullCommand,par2);
    strcat(FullCommand,par3);
    // CallZ++;
    // printf("Rufe %s %d\n",FullCommand,CallZ);
    printf("Calling %s\n",FullCommand);
    fflush(stdout);
    if(system(NULL))
        system(FullCommand);
    else
    {
        printf("Nix geht mehr\n");
        exit(1);
    }
    fflush(stdout);
}

boolean FileExists(char *name, char *suff)
// true wenn das File name.suff existiert
{
   char *FullFileName;
   FILE *Ein;

   FullFileName=(char*) calloc(strlen(name)+strlen(suff)+1,sizeof(char));

  strcpy (FullFileName, name);
  strcat (FullFileName, suff);
  if(Ein=fopen(FullFileName,"r"))
  {
    fclose(Ein);
    return(true);
  }
  return(false);
}

boolean FileExistsRm(char *name, char *suff)
// wie oben, loescht aber File
{
  FILE *Ein;
  char *FullFileName;

   FullFileName=(char*) calloc(strlen(name)+strlen(suff)+1,sizeof(char));

  strcpy (FullFileName, name);
  strcat (FullFileName, suff);
  if(Ein=fopen(FullFileName,"r"))
  {
    fclose(Ein);
    remove(FullFileName);
    return(true);
  }
  return(false);
}

void FileMail(char *Addr, char *FileName)
{
    if(!SENDMAIL) return;
    run_pgm("mail", Addr," < ",FileName);
}

void TextMail(char *Addr, char *Text)
{
    FILE *Aus;

    if(!SENDMAIL) return;

    Aus=fopen("TextMail.Mail","w");
    fprintf(Aus,"%s\n", Text);
    fclose(Aus);
    FileMail(Addr, "TextMail.Mail");
}

void VLMail(char *Addr, VList VL, int Dim)
{
    FILE *Aus;

    if(!SENDMAIL) return;

    Aus=fopen("TextMail.Mail","w");
    FSchreibeVL(Aus,VL,Dim);
    fclose(Aus);
    FileMail(Addr, "TextMail.Mail");
}

clock_t OldTicks=0;

void RunTime()
{
    clock_t ticks;
    if(!TIMER)
        return;
    ticks=clock();
    printf("Time %0.2f\n",(float) (ticks-OldTicks)/CLOCKS_PER_SEC);
    OldTicks=ticks;
}

void GetAndPrintPID(char *progname,char *PIDStr)
{
    int PID;

    PID=getpid();
    printf("*******************\n");
    printf("PID %s %d\n",progname, PID);
    printf("*******************\n");
    sprintf(PIDStr,"%d",PID);
}

void PrintPID(char *progname)
{
    int PID;

    PID=getpid();
    printf("*******************\n");
    printf("PID %s %d\n",progname, PID);
    printf("*******************\n");
}

// Suchen in Datenbasis

// Achtung: im Folgenden Felder ab 1 indiziert
// Deshalb einige Routinen verdoppelt


void liesmatrix1(Matrix mat, int m, int n, FILE *ein)
{
    int i,j;
    for(i=1;i<=m;i++)
        for(j=1;j<=n;j++) Lies(ein,&mat[i][j]);
}

void schreibematrix1(Matrix mat, int m, int n, FILE *aus, int start)
{
    int i,j;
    fprintf(aus,"%d\n",m);
    fprintf(aus,"%d\n",n);
    for(i=start;i<=m;i++)
    {
        for(j=start;j<=n;j++) Schreibe(aus,mat[i][j]);
        fprintf(aus,"\n");
    }
}

void kopieremat_tr1(Matrix mat, Matrix kopie, Matrix tr, int anzerz,int dim)
{
    int i,j;
    for(i=0;i<=anzerz;i++)
    {
        kopie[i][0]=0;
        for(j=1;j<=dim;j++)
        {
            kopie[i][j]=mat[i][j];
            tr[j][i]=mat[i][j];
        }
    }
}

void findepartner(Matrix mat1, Matrix mat2, int m, int n, Matrix zuord)
{
    int i,j;
    for(i=1;i<=m;i++)
    {
        for(j=0;j<=m;j++) zuord[i][j]=0; // Muss wegen mgl. Wiederverwendung
                            // vor Bestimmung zurueckgesetzt werden
        for(j=1;j<=m;j++)
            if (!memcmp(mat1[i], mat2[j], (n+1)*sizeof(Zahl)))
            {
                zuord[i][0]++;
                zuord[i][zuord[i][0]]=j;
            }
    }
}


int matpermvgl_vorb_alloc(Matrix mat1, Matrix mat1_zeis, Matrix mat1_spas, Matrix  mat1_zeiz, Matrix mat1_spaz,
                         Matrix mat2, Matrix mat2_zeis, Matrix mat2_spas, Matrix  mat2_zeiz, Matrix mat2_spaz,
                                int anzerz, int dim, int ***passpal)
// Version mit vorbereiteten Daten fuer mat1 und allokierten Feldern
{
    int i,j,k,zuzuz, zz, zzz, testsp, Rueckgabe;
    int *belegt, *belegt1, *zugez;
    FILE *ein, *aus;

    belegt = (int*) calloc(anzerz+1, sizeof(int));
    belegt1 = (int*) calloc(dim+1, sizeof(int));
    zugez = (int*) calloc(anzerz+1, sizeof(int));


    kopieremat_tr1(mat2,mat2_zeis,mat2_spas,anzerz,dim);
    for(i=0;i<=anzerz;i++) qsort(mat2_zeis[i],dim+1,sizeof(Zahl),Zahlvergl);
    for(j=0;j<=dim;j++) qsort(mat2_spas[j],anzerz+1,sizeof(Zahl),Zahlvergl);

    /* printf("mat1 zeiz\n");
    schreibematrix1(mat1_zeiz,anzerz,anzerz,stdout,0); */

    findepartner(mat1_zeis,mat2_zeis,anzerz,dim,mat2_zeiz);
    /* printf("mat2 zeiz\n");
    schreibematrix1(mat2_zeiz,anzerz,anzerz,stdout,0); */
    findepartner(mat1_spas,mat2_spas,dim,anzerz,mat2_spaz);

    // schreibematrix1(mat2_zeiz, anzerz, anzerz,stdout, 0);
    // schreibematrix1(mat2_spaz, dim, dim, stdout, 0);

    /* printf("Nochmal mat1 zeiz\n");
    schreibematrix1(mat1_zeiz,anzerz,anzerz,stdout,0);
    printf("Nochmal mat2 zeiz\n");
    schreibematrix1(mat2_zeiz,anzerz,anzerz,stdout,0); */


    for(i=1;i<=anzerz;i++) if(mat1_zeiz[i][0] != mat2_zeiz[i][0])
    {
        Rueckgabe=0;
        // printf("Ruas 1111\n");
        goto Ende;
    }
    for(i=1;i<=dim;i++) if(mat1_spaz[i][0] != mat2_spaz[i][0])
    {
        Rueckgabe=0;
        // printf("Ruas 2\n");
        goto Ende;
    }

    // jetzt stimmen Zeilentypen und Haeufigkeiten ueberein

    zuzuz=0;  // Vorbereiten Zeilenzordnung
    for(i=1;i<=anzerz;i++) belegt[i]=0; // wegen Mehrfachverwendung

    for(i=0;i<=dim;i++) // Vorbereiten Tabelle mit passenden Spalten
        for(j=0;j<=dim;j++)
            passpal[0][i][j]=mat2_spaz[i][j];


    while(zuzuz<anzerz)
    {
        zuzuz++; // Versuche nächste Zeile zuzuordnen
        zugez[zuzuz]=0; // Index innerhalb des Feldes mit zuordnungbaren Zeilen
        while(1)
        {
            zugez[zuzuz]++;
            zz=zugez[zuzuz]; // Abkuerzung
            if(zz>mat2_zeiz[zuzuz][0])  // keine Zeile mehr vhd?
            {
                 if(zuzuz==1) // erste Zeile kann nicht mehr zugeordnet werden
                 {
                    Rueckgabe=0;
                    // printf("Ruas 3\n");
                    goto Ende;
                 }
                 else
                 {
                    zuzuz--; // springe Zeile zurueck
                    belegt[mat2_zeiz[zuzuz][zugez[zuzuz]]]=0; // gib deren bisherige Belegung frei
                    continue;
                 }
            }

            zzz=mat2_zeiz[zuzuz][zz]; // Versuche Zuordnung zuzuz --> zzz
            if(belegt[zzz]) continue;  // schon belegt

            // jetzt gleiche Spalten in mat2 finden, Zeilen 1 bis zuzuz in mat 1, zugeodnete in mat2
            // dabei Info in passpal ausnutzen und nur neue Zeile vergleichen

            for(i=1;i<=dim;i++)
            {
                passpal[zuzuz][i][0]=0;
                for(j=1;j<=passpal[zuzuz-1][i][0];j++)
                {
                    testsp=passpal[zuzuz-1][i][j]; // Abkuerzung
                    if(mat1[zuzuz][i]==mat2[zzz][testsp])
                    {
                        passpal[zuzuz][i][0]++;
                        passpal[zuzuz][i][passpal[zuzuz][i][0]]=testsp;
                    }
                }
            }

            // Gleiche Spalten gefunden (oder auch nicht). Versuche Zuordnung.

            for(j=1;j<=dim;j++) belegt1[j]=0;

            for(j=1;j<=dim;j++)
            {
                for(k=1;k<=passpal[zuzuz][j][0];k++) // belege die erste Spalte,
                                                     //  die noch frei ist
                {
                    testsp=passpal[zuzuz][j][k]; // Abkuerzung
                    if(!belegt1[testsp])
                    {
                        belegt1[testsp]=1;
                        break;
                    }
                }
            }
            for(j=1;j<=dim;j++) if(!belegt1[j]) break; // Zuordnung surjektiv ??
            if(j<=dim) continue; // Nein

            belegt[zzz]=1; // Erfolgreich bisher. Belege Zeile zzz
            break; // Versuche naechste Zeile zuzuordnen

        }  // Ende while(1)


    } // Ende while(zuzuz<anzerz)

    Rueckgabe=1; // Erfolgreich bis zur letzten Zeile

    Ende:

    free(belegt);
    free(belegt1);
    free(zugez);

    return Rueckgabe;

}

int sucheinfile(Matrix mat1, int anzerz, int dim, char *filename, int schreiben)
{
    Matrix mat1_zeis, mat1_spas, mat1_zeiz, mat1_spaz;
    Matrix mat2, mat2_zeis, mat2_spas, mat2_zeiz, mat2_spaz;
    int ***passpal;

    int i,j, anzerz2, dim2, dummy, Rueckgabe;
    int iso;
    char zeile[1000], filenamelock[100];

    FILE *ein, *aus,*lock;

    if (!(ein=fopen(filename,"r")))
    {                                           // Datenbasis leer
        if(!schreiben) return(0);

        aus=fopen(filename,"w");
        fprintf(aus,"Values\n");
        schreibematrix1(mat1,anzerz,dim, aus,1);
        fclose(aus);
        return(-1);
    }


    MatAlloc(&mat2,anzerz+1,dim+1);
    MatAlloc(&mat1_zeis,anzerz+1,dim+1); // fuer sortierte Zeilen
    MatAlloc(&mat2_zeis,anzerz+1,dim+1);
    MatAlloc(&mat1_spas,dim+1,anzerz+1);  // fuer sortierte Spalten
    MatAlloc(&mat2_spas,dim+1,anzerz+1);
    MatAlloc(&mat1_zeiz,anzerz+1,anzerz+1); // fuer Indizes der nach Sortierung
    MatAlloc(&mat2_zeiz,anzerz+1,anzerz+1); // gleichen Zeilen bzw.
    MatAlloc(&mat1_spaz,dim+1,dim+1);   // Spalten. Dabei mat1 fuer Vergleich
                                        // mat1 <--> mat1, mat2 fuer mat1 <--> mat2
    MatAlloc(&mat2_spaz,dim+1, dim+1);

    passpal = (int***) calloc(anzerz+1, sizeof(int**)); // Indizes mit gleichen Spalten nach
                                                        // partieller Zuordnung
    for (i=0; i<=anzerz; i++)
    {
        passpal[i] = (int **) calloc(dim+1, sizeof (int*));
        for(j=0; j<= dim; j++) passpal[i][j] = (int*) calloc(dim+1, sizeof(int));
    }

    kopieremat_tr1(mat1,mat1_zeis,mat1_spas, anzerz, dim);
    for(i=0;i<=anzerz;i++) qsort(mat1_zeis[i],dim+1,sizeof(Zahl),Zahlvergl);
    for(j=0;j<=dim;j++) qsort(mat1_spas[j],anzerz+1,sizeof(Zahl),Zahlvergl);
    findepartner(mat1_zeis,mat1_zeis,anzerz,dim,mat1_zeiz); // Haeufigkeiten der sortierten
    findepartner(mat1_spas,mat1_spas,dim,anzerz,mat1_spaz); // Zeilen und Spalten finden

    iso=0;

    while(fscanf(ein,"%s", zeile)>=1)
    {
        if(strncmp(zeile, "Values",6)) continue; // Noch nicht Anfang gefunden
        iso++;

        fscanf(ein,"%d", &anzerz2);
        fscanf(ein,"%d", &dim2);
        if(dim!=dim2 || anzerz!=anzerz2) continue; // Format passt nicht

        liesmatrix1(mat2, anzerz, dim, ein); // Jetzt passt es

        if(matpermvgl_vorb_alloc(mat1, mat1_zeis, mat1_spas, mat1_zeiz, mat1_spaz,
                         mat2, mat2_zeis, mat2_spas,  mat2_zeiz, mat2_spaz,
                           anzerz, dim, passpal))
        {
            fclose(ein);
            Rueckgabe=iso;
            goto Ende;
        }

    }

    fclose(ein);

    if(schreiben)
    {
        strcpy(filenamelock,filename);
        strcat(filenamelock,".lock");
        while(FileExists(filenamelock,""))
        {
            printf("File %s locked\n",filename);
            for(i=0;i<1000000;i++)
            {
                dummy=i;
            }
        }
        lock=fopen(filenamelock,"w");
        aus=fopen(filename,"a");
        fprintf(aus,"Values\n");
        schreibematrix1(mat1,anzerz,dim, aus,1);
        fclose(aus);
        fclose(lock);
        remove(filenamelock);
        Rueckgabe=-(iso+1);
        goto Ende;
    }

    Rueckgabe=0;

    Ende:

    FreeMat(mat2,anzerz+1);
    FreeMat(mat1_zeis,anzerz+1);
    FreeMat(mat2_zeis,anzerz+1);
    FreeMat(mat1_spas,dim+1);
    FreeMat(mat2_spas,dim+1);
    FreeMat(mat1_zeiz,anzerz+1);
    FreeMat(mat2_zeiz,anzerz+1);
    FreeMat(mat1_spaz,dim+1);
    FreeMat(mat2_spaz,dim+1);
    for (i=0; i<=anzerz; i++)
    {
        for(j=0; j<= dim; j++) free(passpal[i][j]);
        free(passpal[i]);
    }
    free(passpal);

    return Rueckgabe;
}


int matinfile(Matrix Gegeben, int anzerz, int dim, char *database, int BasicDim, boolean schreiben, boolean verbose)
// Gibt <0 zurueck, falls nicht vorhanden, sonst Nummer der schon vhd Matrix oder 0,
// falls neu, aber nicht geschrieben
// Schreibt Matrix auf von 0 auf 1-Konvention um. Fuer sucheinfile usw. benoetigt.
{

    Matrix mat1;

    int i,j, iso;

    char filename2[100];

    FILE *ein, *aus;

    // printf("In matvergleich\n");

    sprintf(filename2,"%s_%d.base",database,BasicDim);

    // printf("*** %s \n",filename2);

    MatAlloc(&mat1, anzerz+1,dim+1);

    for(i=0;i<anzerz;i++)
        for(j=0;j<dim;j++)
            mat1[i+1][j+1]=Gegeben[i][j];

    /* printf("Mat 1\n");
    schreibematrix1(mat1,anzerz,dim,stdout,0); */


    iso=sucheinfile(mat1, anzerz, dim, filename2,schreiben);

    // printf("iso %d\n",iso);

    FreeMat(mat1,anzerz+1);

    if(iso>0)
    {
        if(verbose)
            printf("Already in  %s Typ %d\n", filename2, iso);
        return(iso);
    }
    if(iso<0)
    {
        printf("New in  %s Typ %d\n", filename2, -iso);
        return(iso);
    }

    printf("Not in  %s\n", filename2);
    return 0;

}

boolean FindInBase(VList Gen, VList Supp, int Dim, char *database,
           boolean Write, boolean Verbose)
// Suche in Datenbasis. Falls neu und Schreibe==true, eintragen
// Falls Verbose==false, Meldung "Already ..." nicht ausgeben
// Erzeugt "val"-Matrix durch Auswerten von Supp auf Gen
{
    VLE *VP, *VP1;
    Zahl Wert;
    int i,j,k;
    Matrix mat;
    boolean vhd;

    MatAlloc(&mat,NinVL(Gen),NinVL(Supp));

   i=0;
   for(VP=Gen.St;VP;VP=VP->R)
   {
     j=0;
     for(VP1=Supp.St;VP1;VP1=VP1->R)
     {
          Wert=0;
          for(k=0;k<Dim;k++)
            Wert+=VP->V[k]*VP1->V[k];
          mat[i][j]=Wert;
          j++;
     }
     i++;
   }

   vhd=false;

   if(matinfile(mat,i,j,database, Dim,Write,Verbose)>0)
      vhd=true;

   FreeMat(mat,i);

   return vhd;
}

boolean FindInBaseMat(Matrix Total,int ZZ, int SZ ,int Dim,char *database,
                    boolean Write, boolean Verbose)
// Analog zu FindIn Base, nur dass die "val"-Matrix schon in Total gegeben ist
{
    boolean vhd;

    vhd=false;

   if(matinfile(Total,ZZ,SZ,database, Dim,Write,Verbose)>0)
      vhd=true;

   return vhd;

}

// Spezielle Routinen fuer "glatte" Polytope

void MakeTransitive(Matrix Con, int N)
// makes a relation matrix Con of size NxN transitive
{
    boolean NewCon;
    int i,j,k;

    do
    {
        NewCon=false;
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                for(k=0;k<N;k++)
                    if(Con[i][j] && Con[j][k] && !Con[i][k])
                    {
                        NewCon=true;
                        Con[i][k]=1;
                    }

    }
    while(NewCon);
}

void Connect(VList LattPoints, VList Steps, int Dim, Matrix Con)
// "Verbindet" Gitterpunkte, die sich durch ein Element aus Steps
// unterscheiden, Con gibt die Zusammenhangsmatrix zurueck
{
    int i,j,k,NLattPoints;
    Vektor Diff;
    VLE *VP,*VP1;

    Diff=calloc(Dim,sizeof(Zahl));
    NLattPoints=NinVL(LattPoints);

    for(i=0;i<NLattPoints;i++)
        for(j=0;j<NLattPoints;j++)
            Con[i][j]=0;

    i=-1;
    for(VP=LattPoints.St;VP;VP=VP->R)
    {
        i++;
        j=-1;
        for(VP1=LattPoints.St;VP1;VP1=VP1->R)
        {
            j++;
            if(i==j)
            {
                Con[i][j]=1;
                continue;
            }
            for(k=0;k<Dim;k++)
                Diff[k]=VP1->V[k]-VP->V[k];
            if(VinVL(Diff,Steps,Dim))
                Con[i][j]=1;
        }
    }

    MakeTransitive(Con,NLattPoints);

    free(Diff);
}

boolean HilbConnected(VList PolytGen, int Dim, boolean IsVeryAmple,
                        boolean *HConn, boolean *VConn, boolean *SConn)
// Prueft auf Hilbert-Zusammenhang (was auch immer das ist)
// PolytGen spannt ein Giterpolytop der Dimension Dim-1 auf
// IsVeryAmple zeigt an, ob schon bekannt, dass PolytGen die Gitterpunkte
// eines very ample Polytops enthaelt
// HConn: Hilb connected, VConn: sogar strongly, SConn: sogar super
// Rueckgabewert ist false, wenn HilbConn nicht angewendet werden konnte
{
    VList LattPoints,Supp, Extr,Corner, CornerBas,CornerSupp,AccuHilb,HilbBas;
    int Rang, i,j,k,NLattPoints;
    Zahl Index;
    VLE *VP, *VP1;
    Matrix Con;
    Vektor VeryCon;
    boolean HilbConnected, Simple,Smooth,VA;

    if(!IsVeryAmple)
    {
        HilbBasTri(PolytGen,Dim,&HilbBas,&Supp,&Rang, &Index);
        if(Rang<Dim)
        {
            printf("Polytope for HilbConnected not of full dimension\n");
            FreeVL(&Supp);
            FreeVL(&LattPoints);
            FreeVL(&HilbBas);
            return false;
        }
        StartVL(&LattPoints);
        for(VP=HilbBas.St;VP;VP=VP->R)
            if(VP->V[Dim-1]==1)
                CopyVLE0(&LattPoints,VP,Dim);
        VASS(LattPoints,Supp,HilbBas,Dim,false,&VA,&Simple,&Smooth);
        FreeVL(&HilbBas);
        if(!VA)
        {
            printf("Polytope for HilbConnected not very ample\n");
            FreeVL(&Supp);
            FreeVL(&LattPoints);
            return false;
        }

    }
    else
    {
        LattPoints=PolytGen;
        SuppHypVL(&Supp,LattPoints,Dim);
    }
    NLattPoints=NinVL(LattPoints);

    MatAlloc(&Con,NLattPoints,NLattPoints);
    VeryCon=calloc(NLattPoints, sizeof(int));
    StartVL(&AccuHilb);

    *SConn=true;

    SelExtremeR(&Extr,LattPoints,Supp,Dim);

    for(VP=Extr.St;VP;VP=VP->R)
    {
        CornerCone(&Corner,LattPoints,&CornerSupp,Supp,VP->V, Dim);
        FreeVL(&CornerSupp);
        HilbBasTri(Corner,Dim-1,&CornerBas,&CornerSupp,&Rang,&Index);
        FreeVL(&CornerSupp);
        FreeVL(&Corner);

        Connect(LattPoints,CornerBas,Dim-1,Con);

        VP1=VLEinVL(VP,LattPoints,Dim-1);
        i=IndVLEinVL(VP1,LattPoints);

        for(j=0;j<NLattPoints;j++)
        {
            VeryCon[j]=VeryCon[j] || Con[i][j];
            if(!Con[i][j])
            {
                *SConn=false;
                // SchreibeV(stdout,VP1->V,Dim-1);
                // printf("nicht verbunden mit %d\n",j);
            }
        }
        ConcatVLFree(&AccuHilb,&CornerBas);

    }

    *VConn=true;
    for(j=0;j<NLattPoints;j++)
        if(!VeryCon[j])
        {
             // printf("Not very connected\n");
             *VConn=false;
        }

    FreeVL(&Extr);
    FreeVL(&Supp);

    // SchreibeVL(stdout,AccuHilb,Dim-1);

    *HConn=true;

    if(!*VConn)
    {
        Connect(LattPoints,AccuHilb,Dim-1,Con);
        for(i=0;i<NLattPoints;i++)
        {
            for(j=0;j<NLattPoints;j++)
                if(!Con[i][j])
                {
                    *HConn=false;
                    goto Ende;
                }
        }
    }

    Ende:

    // SchreibeMatrix(stdout,Con,NLattPoints,NLattPoints);

    FreeVL(&AccuHilb);
    free(VeryCon);
    FreeMat(Con,NLattPoints);
    if(!IsVeryAmple)
        FreeVL(&LattPoints);
    return true;
}

boolean Deg2AffCover(VList PolytGen, int Dim, VList PolytSupp, boolean IsVeryAmple,
        boolean HasSupp, boolean *SuperConn)
//
// OBSOLET Ersetzt durch NewDeg2AffCover
//
// Prueft, ob affine Ueberdeckung glattet torischer Varietaet
// im Grad 2 definiert ist. 
// außerdem Check auf "Super Hilb connected" (was eine gewisse
// Verschaerfung davon darstellt)
{
    VList LattPoints,Supp, Extr,Corner, CornerBas,CornerSupp,HilbBas;
    int Rang, i,j,k,NLattPoints, Dim1;
    Zahl Index;
    VLE *VP, *VP1, *VP2;
    Vektor Diff,Hom;
    boolean Reduced, SConn;

    boolean NotDeg2, Simple,Smooth,VA;

    Dim1=Dim-1;

    if(!IsVeryAmple)
    {
        HilbBasTri(PolytGen,Dim,&HilbBas,&Supp,&Rang, &Index);
        if(Rang<Dim)
        {
            printf("Polytope for Deg2FinL not of full dimension\n");
            FreeVL(&Supp);
            FreeVL(&LattPoints);
            FreeVL(&HilbBas);
            return false;
        }
        StartVL(&LattPoints);
        for(VP=HilbBas.St;VP;VP=VP->R)
            if(VP->V[Dim-1]==1)
                CopyVLE0(&LattPoints,VP,Dim);
        VASS(LattPoints,Supp,HilbBas,Dim,false,&VA,&Simple,&Smooth);
        FreeVL(&HilbBas);
        if(!Smooth)
        {
            printf("Polytope for Deg2FinL not smooth\n");
            FreeVL(&Supp);
            FreeVL(&LattPoints);
            return false;
        }

    }
    else
    {
        LattPoints=PolytGen;
        if(!HasSupp)
            SuppHypVL(&Supp,LattPoints,Dim);
        else
            Supp=PolytSupp;
    }

    SelExtremeR(&Extr,LattPoints,Supp,Dim);
    Diff=calloc(Dim,sizeof(Zahl));
    Hom=calloc(Dim1,sizeof(Zahl));

    Reduced=true;
    *SuperConn=true;

    for(VP=Extr.St;VP;VP=VP->R)
    {
        // printf("Hurra\n");

        // printf("Ecke ");SchreibeV(stdout,VP->V,Dim1);
        CornerCone(&Corner,LattPoints,&CornerSupp,Supp,VP->V, Dim);

        // HilbBasTri(Corner,Dim1,&CornerBas,&CornerSupp,&Rang,&Index);

        SuppHypVL(&CornerBas,CornerSupp,Dim1);

        /* CSchreibeVL(stdout,"Corner",Corner,Dim1);
        CSchreibeVL(stdout,"CornerSupp",CornerSupp,Dim1);
        CSchreibeVL(stdout,"CornerBas",CornerBas,Dim1);*/


        TesteHom(Hom,CornerBas,Dim1);
        //SchreibeV(stdout,Hom,Dim1);

        for(VP1=Corner.St;VP1;VP1=VP1->R)
            VP1->Mark1=SkalProd(Hom,VP1->V,Dim1);

        for(VP1=Corner.St;VP1;VP1=VP1->R)
        {
            //printf("Hurra2\n");
            // SchreibeV(stdout,VP1->V,Dim1);

            if(VP1->Mark1==1)
                continue;

            // printf("Untersuche "); SchreibeV(stdout,VP1->V,Dim1);

            Reduced=false;
            SConn=false;

            for(VP2=CornerBas.St;VP2;VP2=VP2->R) // Erst die Nachbarn abziehen
            {
                for(k=0;k<Dim1;k++)
                    Diff[k]=VP1->V[k]-VP2->V[k];
                if(!InCone(Diff,CornerSupp,Dim1))
                    continue;

                for(k=0;k<Dim1;k++)             // Ecke wieder addieren
                    Diff[k]=Diff[k]+VP->V[k];
                Diff[Dim1]=1;                   // fertig damit

                if(InCone(Diff,Supp,Dim))
                {
                    // printf("Abgezogen "); SchreibeV(stdout,VP2->V,Dim1);
                    Reduced=true;
                    SConn=true;
                    break;
                }
            }

            if(!SConn)
                *SuperConn=false;
            if(Reduced) continue;

            // printf("Kein Abzug\n");

            for(VP2=Corner.St;VP2;VP2=VP2->R)
            {
                if(VP2->Mark1 >= VP1->Mark1)
                    continue;
                for(k=0;k<Dim1;k++)
                    Diff[k]=VP1->V[k]-VP2->V[k];
                if(!InCone(Diff,CornerSupp,Dim1))
                    continue;

                for(k=0;k<Dim1;k++)             // Ecke wieder addieren
                    Diff[k]=Diff[k]+VP->V[k];
                Diff[Dim1]=1;                  // fertig damit

                if(InCone(Diff,Supp,Dim))
                {
                    Reduced=true;
                    break;
                }
            }
            if(!Reduced)
                break;
        }

        FreeVL(&CornerSupp);
        FreeVL(&Corner);
        FreeVL(&CornerBas);
        if(!Reduced)
            break;
    }

    if(!IsVeryAmple)
    {
        FreeVL(&LattPoints);
        FreeVL(&Supp);
    }
    else
        if(!HasSupp)
            FreeVL(&Supp);

    FreeVL(&Extr);

    free(Diff);
    free(Hom);

    return Reduced;
}

boolean NewDeg2AffCover(VList Extr, VList LattPoints, VList Supp,int Dim,boolean *SuperConn, boolean *NoIsolatedPoint)
// Neue Version, spart Speicherplatz
// Extr: Ecken des Polytops
// LattPoints:  alle Gitterpunkte
// Prueft, ob affine Ueberdeckung glattet torischer Varietaet
// im Grad 2 definiert ist. 
// außerdem Check auf "Super Hilb connected" (was eine gewisse
// Verschaerfung davon darstellt)
// Mark1 wird fuer "Eckengrad" benutzt
// Mark3 fuer isolierten Punkt
{
    VList CornerBas,CornerSupp;
    int Rang, i,j,k,NLattPoints, Dim1;
    Zahl ValExtr;
    VLE *VP, *VP1, *VP2;
    Vektor Diff,Hom,Ecke,Punkt;
    boolean Reduced, SConn;

    boolean NotDeg2, Simple,Smooth,VA;

    Dim1=Dim-1;
    
    Diff=calloc(Dim,sizeof(Zahl));
    Hom=calloc(Dim1,sizeof(Zahl));

    Reduced=true;
    *SuperConn=true;
    *NoIsolatedPoint=true;
    
    for(VP1=LattPoints.St;VP1;VP1=VP1->R) // Vorbereiten Test fuer "isolated point"
        VP1->Mark3=0;                     // Erst einmal alle Punkte isoliert
        
    for(VP=Extr.St;VP;VP=VP->R)
    {
        // printf("Hurra Ecke\n");

        Ecke=VP->V;

        CornerHyp(&CornerSupp,Supp,Ecke, Dim);
        
        SuppHypVL(&CornerBas,CornerSupp,Dim1);

        /* //CSchreibeVL(stdout,"Corner",Corner,Dim1);
        CSchreibeVL(stdout,"CornerSupp",CornerSupp,Dim1);
        CSchreibeVL(stdout,"CornerBas",CornerBas,Dim1); */
        
        TesteHom(Hom,CornerBas,Dim1);
        //SchreibeV(stdout,Hom,Dim1);

        ValExtr=SkalProd(Hom,Ecke,Dim1); // Graduierung im Eckenkegel machen
        for(VP1=LattPoints.St;VP1;VP1=VP1->R)
            VP1->Mark1=SkalProd(Hom,VP1->V,Dim1)-ValExtr; // gemacht

        for(VP1=CornerBas.St;VP1;VP1=VP1->R) // Ecke wieder aufaddieren
                                             // so dass Nachbarpunkte berechnet werden
            for(k=0;k<Dim1;k++)
                VP1->V[k]+=Ecke[k];

        // CSchreibeVL(stdout,"CornerBas nach Add",CornerBas,Dim1);

        for(VP1=LattPoints.St;VP1;VP1=VP1->R)
        {
            Punkt=VP1->V;

            if(VP1->Mark1<=1)    //Ecke selbst oder Nachbar
            {
                VP1->Mark3=1;    // ist an Ecke angebunden
                continue;
            }
            
            // printf("Untersuche "); SchreibeV(stdout,Punkt,Dim1);

            Reduced=false;
            SConn=false;

            for(VP2=CornerBas.St;VP2;VP2=VP2->R) // Erst die Nachbarn abziehen
            {
                for(k=0;k<Dim1;k++)
                    Diff[k]=Punkt[k]-VP2->V[k];
                if(!InCone(Diff,CornerSupp,Dim1)) // Differenz nicht einmal im Eckenkegel
                    continue;

                for(k=0;k<Dim1;k++)             // Ecke addieren
                    Diff[k]+=Ecke[k];
                Diff[Dim1]=1;                   // fertig damit

                if(InCone(Diff,Supp,Dim))
                {
                    // printf("Abgezogen "); SchreibeV(stdout,VP2->V,Dim1);
                    Reduced=true;
                    SConn=true;
                    VP1->Mark3=1;
                }
            }

            if(!SConn)
                *SuperConn=false;
            if(Reduced) continue;

            // printf("Kein Abzug\n");

            for(VP2=LattPoints.St;VP2;VP2=VP2->R)  // Jetzt mit allen Punkten arbeiten
            {
                if(VP2->Mark1 >= VP1->Mark1) // Nur solche mit kleinerem Grad testen
                    continue;
                for(k=0;k<Dim1;k++)
                    Diff[k]=Punkt[k]-VP2->V[k];
                if(!InCone(Diff,CornerSupp,Dim1))
                    continue;

                for(k=0;k<Dim1;k++)             // Ecke wieder addieren
                    Diff[k]=Diff[k]+Ecke[k];
                Diff[Dim1]=1;                  // fertig damit

                if(InCone(Diff,Supp,Dim))
                {
                    Reduced=true;
                    break;
                }
            }  
            if(!Reduced) 
                break;
                
        } // Schleife ueber Gitterpunkte

        FreeVL(&CornerSupp); // Ecke ist abgearbeitet
        FreeVL(&CornerBas);
        if(!Reduced)
            break;
    }  // Schleife ueber Ecken
    
    for(VP1=LattPoints.St;VP1;VP1=VP1->R)
        if(!VP1->Mark3)
        {
            *NoIsolatedPoint=false;
            printf("********* Isolated point ");
            SchreibeV(stdout,VP1->V,Dim);
            fflush(stdout);
        }

    free(Diff);
    free(Hom);

    return Reduced;
}

boolean AbundantDeg2Rel(VList Extr, VList LattPoints, VList Supp, int Dim)
{
    VLE *VP,*VP1,*VPTest;
    int i,NSupp; 
    Vektor Summe; 

    NSupp=NinVL(Supp);
    Summe=calloc(NSupp,sizeof(Zahl));
        
    for(VP=LattPoints.St;VP;VP=VP->R)
    {
        VP->Mark3=0;
        if(VLEinVL(VP,Extr,Dim))
            VP->Mark3=1;    
    }
    AddValLV(LattPoints,Supp,Dim);
    
    for(VP=LattPoints.St;VP;VP=VP->R)
    {
        if(VP->Mark3)                 // Ecken gegen alle anderen Punkte in Deg2AffCover gecheckt
            continue;
        for(VP1=VP;VP1;VP1=VP1->R)     // Testen auch x^2=yz
        {
            if(VP1->Mark3)  // dito
                continue;
            AddVekt(Summe,VP->LV,VP1->LV,NSupp);
            
            for(VPTest=LattPoints.St;VPTest;VPTest=VPTest->R)
            {
                if(VPTest==VP || VPTest==VP1 || !VPTest->Mark3) // x und y ausschließen, zunaechst nur Ecken
                    continue;
                for(i=0;i<NSupp;i++)
                    if(Summe[i]-VPTest->LV[i]<0)
                        break;
                if(i==NSupp)  // Gitterpunkt in P \cap (x+y-P) gefunden
                    break;
            }
            if(VPTest) 
                continue;// Gefunden und naechste Schleife ueberspringen            
                
            for(VPTest=LattPoints.St;VPTest;VPTest=VPTest->R)
            {
                if(VPTest==VP || VPTest==VP1 || VPTest->Mark3) // x und y ausschließen, ebenso Ecken
                    continue;
                for(i=0;i<NSupp;i++)
                    if(Summe[i]-VPTest->LV[i]<0)
                        break;
                if(i==NSupp)  // Gitterpunkt in P \cap (x+y-P) gefunden
                    break;
            }
            if(!VPTest)  // Keiner gefunden
            {
                free(Summe);
                return(false);
            }
                        
        }  
    }
    free(Summe);
    return(true); 
}

boolean StronglyConnected(VList Extr, VList LattPoints, VList Supp, int Dim)
{
    VLE *VP, *VP1, *VPEinLP, **Conn;  // Conn nimmt die schon mit Ecke verbundenen Punkte auf
    
    VList CornerBas,CornerSupp;
    int Rang, i,j,k, Dim1,NConn;
    Zahl ValExtr,Grad;
    Vektor Diff,Hom,Ecke;
    boolean Weiter;
    
    Dim1=Dim-1;
    
    // printf("Punkte %d\n",NinVL(LattPoints));
    
    Diff=calloc(Dim,sizeof(Zahl));
    Hom=calloc(Dim1,sizeof(Zahl));
    Conn=calloc(NinVL(LattPoints),sizeof(VLE*));
    
    for(VP1=LattPoints.St;VP1;VP1=VP1->R) // Vorbereiten Test fuer "nicht verbunden mit irgendeiner Ecke"
        VP1->Mark3=0;                     // Erst einmal kein Punkt verbunden

    for(VP=Extr.St;VP;VP=VP->R)
    {
        // printf("Hurra Ecke "); SchreibeV(stdout,VP->V,Dim);

        Ecke=VP->V;  // Abkuerzung
        
        VPEinLP=VinVL(Ecke,LattPoints,Dim);  // Ecke in Gitterpunkten finden

        CornerHyp(&CornerSupp,Supp,Ecke, Dim); // Machen erst CorneSupp        
        SuppHypVL(&CornerBas,CornerSupp,Dim1); // CornerBas mittels dualisieren 
        FreeVL(&CornerSupp); // Nicht mehr gebraucht
        
        // SchreibeVL(stdout,CornerBas,Dim1);
        // printf("***********\n");
        
        TesteHom(Hom,CornerBas,Dim1);  // Ist naruerlich homogen,aber Hom muss berechnet werden

        ValExtr=SkalProd(Hom,Ecke,Dim1); // Graduierung im Eckenkegel machen
        VPEinLP->Mark1=ValExtr;          // In LatPoints Grad der Ecke richtig setzen
        for(VP1=LattPoints.St;VP1;VP1=VP1->R)
            VP1->Mark1=SkalProd(Hom,VP1->V,Dim1)-ValExtr; // gemacht  
            
        for(VP1=LattPoints.St;VP1;VP1=VP1->R) // Vorbereiten Test fuer "nicht verbunden mit irgendeiner Ecke"
            VP1->Mark2=0;                     // Erst einmal kein Punkt verbunden
            
        VPEinLP->Mark2=1; // Ecke mit sich selbst bverbunden
        Conn[0]=VPEinLP;
        NConn=1;
        
        for(i=0;i<NConn;i++)
        {
            
            Grad=Conn[i]->Mark1+1;
                        
            for(VP1=LattPoints.St;VP1;VP1=VP1->R)
            {

                // printf("Schleife %d ** ",VP1->Mark1); SchreibeV(stdout,VP1->V,Dim);
                
                if(VP1->Mark2 || VP1->Mark1!=Grad)  // falscher Grad oder schon verbunden
                    continue;
                SubtVekt(Diff,VP1->V,Conn[i]->V,Dim1);
                
                // printf("Diff ");SchreibeV(stdout,Diff,Dim1);
                
                if(VinVL(Diff,CornerBas,Dim1))
                {
                    VP1->Mark2=1;  // jetzt verbunden
                    Conn[NConn]=VP1; // als verbunden registrieren
                    NConn++;
                    // printf("Verbunden %d %lld\n",NConn,VP1->Mark1);
                }
            }
            
            // if(i==NinVL(LattPoints)) exit(1);
        
        }
        
        // exit(1);
        
        FreeVL(&CornerBas); // Ecke ist abgearbeitet
        
        Weiter=false;
        for(VP1=LattPoints.St;VP1;VP1=VP1->R)  // Fuer jeden Punkt registrieren, ob er verbunden ist
        {
            if(VP1->Mark2)
            {
                VP1->Mark3=1;
                // printf("+++ ");SchreibeV(stdout,VP1->V,Dim);
            }
            
            if(!VP1->Mark3)
                Weiter=true;
            
        }
        
        if(!Weiter)  // Schon alle Punkte verbunden
            break;
            
        
    }  // Schleife ueber Ecken
        
    free(Conn);
    free(Diff);
    free(Hom);
    
    if(!Weiter)  // Bei !Weiter alle Punkte verbunden
        return true;
    else
    {        
        for(VP1=LattPoints.St;VP1;VP1=VP1->R)  // Testen ob jeder Punkt mit irgendeiner Ecke verbunden
            if(!VP1->Mark3)
            {
                printf("***** Not connected to a vertex ");
                SchreibeV(stdout,VP1->V,Dim);
                fflush(stdout);
            }
    
        return false; 
    }   
}

void MakeNormalFan(SimplexList *NormalFan,VList Extr,VList Supp,int Dim)
// Macht Normalenfaecher von einfachem Polytop mit Ecken Extr und Hyperebenen Supp
// Dim=Polytop-Dimension+1
// Corn und Supp der simplizialen Kegel werden NICHT gesetzt
{
    int i,j;
    VLE *VP,*VP1;
    Simplex *Simpl;

    StartSimplexL(NormalFan);

    for(VP=Extr.St;VP;VP=VP->R)
    {
        Simpl=NewSimplexLE(NormalFan);
        i=0;
        j=-1;
        for(VP1=Supp.St;VP1;VP1=VP1->R)
        {
            j++;
            if(!SkalProd(VP->V,VP1->V,Dim))
            {
                if(i>=Dim-1)
                {
                    printf("MakeNormalFan non-simple polytope\n");
                    exit(0);
                }
                Simpl->Vert[i]=j;
                i++;
            }
        }
    }
}

void CutCorn(VList *NewExtr, VList *NewSupp, VList Extr, VList Supp, int Dim)
// Versucht Polytop durch Abschneiden von Ecken zu verkleinern
// Spezialfall von DoCut (s.u.), kann aber ohne Vorbereitung eingesetzt werden
{

    VList Corner,CornerSupp,CornerBas;
    VLE *VP,*VP1,*VP2,*VPMax;
    Zahl Test,Min,Max,Dist;
    Vektor Diff;
    int i,j,Dim1;

    CopyVL(NewExtr,Extr,Dim);
    CopyVL(NewSupp,Supp,Dim);
    Dim1=Dim-1;
    Diff=calloc(Dim1,sizeof(Zahl));

    while(1)
    {
        Max=0;

        /* printf("-------\nStart\n-----------\n");
        FSchreibeVL(stdout,*NewExtr,Dim);
        FSchreibeVL(stdout,*NewSupp,Dim);
        printf("--------\n"); */


        for(VP=NewExtr->St;VP;VP=VP->R)
        {
            /* printf("Teste ");
            SchreibeV(stdout,VP->V, Dim); */
            CornerCone(&Corner,*NewExtr,&CornerSupp,*NewSupp,VP->V,Dim);
            Min=0;

            for(VP1=NewExtr->St;VP1;VP1=VP1->R) // Ermittle Minimaldistanz von VP zu Nachbarn
            {
                if(VP1==VP) continue; // Ecke nicht mit sich selbst vergleichen
                for(i=0;i<Dim1;i++)
                    Diff[i]=VP1->V[i]-VP->V[i];
                j=0;
                for(VP2=CornerSupp.St;VP2;VP2=VP2->R)
                {
                    Test=SkalProd(VP2->V,Diff,Dim1);
                    if(Test)
                    {
                        if(j) break; // Keine Nachbarecke
                        Dist=Test;
                        j++;
                    }
                }

                /* if(!VP2)
                {
                    printf("Nachbar ");
                    SchreibeV(stdout,VP1->V,Dim);
                    printf("Abstand %d\n",Dist);
                }*/
                if(!VP2) // Nachbarecke
                    if(Min==0 || Dist <Min)
                        Min=Dist;
            }
                    // Minimaldistanz fuer diese Ecke gefunden

            if(Max<Min)  // neuer maximaler Minimalabstand
            {
                Max=Min;
                VPMax=VP; // neue Maximalecke
            }

            FreeVL(&Corner);
            FreeVL(&CornerSupp);
        }

        if(Max==1) break; // nix zum Abschneiden gefunden

        // jetzt Maximalecke abschneiden, wenn Abstand >1

        /* printf("Max %d\n",Max);
        SchreibeV(stdout,VPMax->V,Dim1); */

        CornerCone(&Corner,*NewExtr,&CornerSupp,*NewSupp,VPMax->V,Dim);
        SuppHypVL(&CornerBas,CornerSupp,Dim1);
        /*printf("CornerSupp\n");
        SchreibeVL(stdout,CornerSupp,Dim1);
        printf("CornerBas\n");
        SchreibeVL(stdout,CornerBas,Dim1);
        printf("---------\n"); */


        for(VP=CornerBas.St;VP;VP=VP->R) //Neue Ecken machen
        {
            VP1=NewVLE(NewExtr);
            for(i=0;i<Dim1;i++)
                VP1->V[i]=VPMax->V[i]+(Max-1)*VP->V[i];
            VP1->V[Dim1]=1;
        }

        FreeVLE(NewExtr,VPMax); // alte Ecke entfernen
        /* printf("Neue Ecken\n");
        SchreibeVL(stdout,*NewExtr,Dim1); */
        FreeVL(&CornerBas);
        FreeVL(&Corner);
        FreeVL(&CornerSupp);
        FreeVL(NewSupp);
        SuppHypVL(NewSupp,*NewExtr,Dim);
    }

    free(Diff);
}


boolean DoCut(VList *Extr, VList *Supp, int Dim, int CutDim, VList Parallels)
{
// Schneidet von (glattem) Polytop ein Stück ab, das parallel zu einer Seite
// der Dimension (Dim-1)-CutDim ist. CutDim beziegt sich auf eine Seite des Normalen-
// faechers. Also: Wenn parallel zu den Facetten geschnitten werden soll ist
// CutDim=0 usw. Polytop-Dimension= Dim-1.
// Aus den Parallelen wird diejenige ausgesucht, die am tiefsten in das
// Polytop hineinschneidet, ohne vorhandene Facetten zu vernichten.
// Die Parallele wird an Supp angehängt (wenn CutDim>0).
// Dann werden die Ecken neu berechnet.

    Zahl Min,Max,Test;
    VLE *VP,*VP1,*VPMax;


    Max=0; // Maximalabstand der Hyperebenen zu jeweils den Ecken, die nicht drin liegen
    /* printf("Parallelen\n");
    FSchreibeVL(stdout,Parallels,Dim); */

    for(VP=Parallels.St;VP;VP=VP->R)
    {
        Min=0; // Minimalwert der Linearform auf den Ecken, auf denen sie NICHT verschwindet
        for(VP1=Extr->St;VP1;VP1=VP1->R) // Der wird jetzt ermittelt
        {
            Test=SkalProd(VP->V,VP1->V,Dim);
            if(!Test) continue;
            if(Min==0 || Test<Min)
                Min=Test;
        }
        if(Max==0 || Min>Max)
        {
            Max=Min;  // der neue Maximalabstand
            VPMax=VP;
        }
    }

    /* printf("CutDim %d Max %lld\n",CutDim,Max);
    printf("Gewaehlt ");SchreibeV(stdout,VPMax->V,Dim); */

    if(Max<=1)
        return false;

    if(!CutDim)
        VP=VLEinVL(VPMax,*Supp,Dim-1); // Vorhandene Hyperebene wird verschoben
    else
    {
        VP=NewVLE(Supp);
        memcpy(VP->V,VPMax->V,Dim*sizeof(Zahl)); // Neue Hyperebene wird hinzugefügt
    }
    VP->V[Dim-1]-=Max-1;  // Wir können unsere Nyperebene um Maximalabstand-1 verschieben

    FreeVL(Extr);
    SuppHypVL(Extr,*Supp,Dim);
    /* CSchreibeVL(stdout,"Supp",*Supp,Dim);
    CSchreibeVL(stdout,"Extr",*Extr,Dim); */
    return true;
}

boolean DoCutAvoidForbidden(VList *Extr, VList *Supp, VList Forbidden, int Dim, int CutDim, VList Parallels)
{
// Schneidet von (glattem) Polytop ein Stück ab, das parallel zu einer Seite
// der Dimension (Dim-1)-CutDim ist. CutDim beziegt sich auf eine Seite des Normalen-
// faechers. Also: Wenn parallel zu den Facetten geschnitten werden soll ist
// CutDim=0 usw. Polytop-Dimension= Dim-1.
// Aus den Parallelen wird diejenige ausgesucht, die am tiefsten in das
// Polytop hineinschneidet, ohne vorhandene Facetten zu vernichten
// oder VERBOTENE PUNKTE ZU ENTFERNEN.
// IM GEGENSATZ ZU DEN ECKEN DARF DURCH SIE GESCHNITTEN WERDEN
// OHNE "SMOOTH" ZU VERLETZEN (Wenn sie keine Ecken sind)
// Die Parallele wird an Supp angehängt (wenn CutDim>0).
// Dann werden die Ecken neu berechnet.

    Zahl Max,Test,SchnittAbstand,ForbiddenMin,EckenMin;
    VLE *VP,*VP1,*VPMax;


    Max=0; // IN DIESER ROUTINE WIRD Max ETWAS ANDERS BEHANDELT ALS IN DoCut
    /* printf("---------------\n");
    printf("Parallelen\n");
    FSchreibeVL(stdout,Parallels,Dim);
    printf("---------------\n"); */

    for(VP=Parallels.St;VP;VP=VP->R)
    {
        // printf("Parallele ");SchreibeV(stdout,VP->V,Dim);
        EckenMin=0; 
        for(VP1=Extr->St;VP1;VP1=VP1->R) // Ermitteln wie weit die Parallele verschoben werden kann
        {
            Test=SkalProd(VP->V,VP1->V,Dim);
            if(!Test) continue;
            if(EckenMin==0 || Test <EckenMin)
                EckenMin=Test;
        }
        
        ForbiddenMin=-1;
        for(VP1=Forbidden.St;VP1;VP1=VP1->R) // Hier werden die verbotenen Punkte berücksichtigt
        {
            Test=SkalProd(VP->V,VP1->V,Dim);
            // printf("Wert %lld auf ",Test);SchreibeV(stdout,VP1->V,Dim);
            if(ForbiddenMin==-1 || Test<ForbiddenMin)
                ForbiddenMin=Test;
        }
        // printf("EckenMin %lld ForbiddenMin %lld\n",EckenMin, ForbiddenMin);
        SchnittAbstand=EckenMin-1;
        if(ForbiddenMin<SchnittAbstand)
            SchnittAbstand=ForbiddenMin;
        if(Max==0 || SchnittAbstand>Max)
        {
            Max=SchnittAbstand;
            VPMax=VP;
        }
        // printf("DchnittAbstand %lld Max %lld\n",SchnittAbstand,Max);
    }

    /* printf("CutDim %d Max %lld\n",CutDim,Max);
    printf("Gewaehlt ");SchreibeV(stdout,VPMax->V,Dim); */

    if(Max<=0)
        return false;

    if(!CutDim)  // Vorhandene Hyperebene wird verschoben
        VP=VLEinVL(VPMax,*Supp,Dim-1);
    else
    {
        VP=NewVLE(Supp); // Neue Hyperebene wird hinzugefügt
        memcpy(VP->V,VPMax->V,Dim*sizeof(Zahl));
    }
    VP->V[Dim-1]-=Max; // Hier -1 entfernt, weil schon oben berfücksichtigt

    FreeVL(Extr);
    SuppHypVL(Extr,*Supp,Dim);
    /* CSchreibeVL(stdout,"Supp",*Supp,Dim);
    CSchreibeVL(stdout,"Extr",*Extr,Dim); */
    return true;
}

void MakeFacesParallels(VList *Parallels, VList Supp,SimplexList NormalFan,int Dim,int CutDim,
                        SimplexList FacesSimp)
// Macht Seitenparallelen zu einer ausgewählten Dimension CutDim im dualen
// Fächer. Wird abgelegt in Parallels. Normaler Fächer durch Supp und NormalFan repräsentiert
// FasesSimp enthaelt Seitenverband des Einheitssimplex (muss vorher gemacht werden)
{

    Simplex *Simpl,*Simpl1;
    Vektor NewPar;
    int i,j,k;
    VLE *VP;

    StartVL(Parallels);
    NewPar=calloc(Dim,sizeof(Zahl));

    for(Simpl=FacesSimp.St;Simpl;Simpl=Simpl->R)
    {
        if(Simpl->Vert[CutDim]<0) // zu kleine Dimension ueberspringen
            continue;
        if(CutDim < Dim-2 && Simpl->Vert[CutDim+1]>=0) // Dimension schon zu gross
            break;

        for(Simpl1=NormalFan.St;Simpl1;Simpl1=Simpl1->R)
        {
            for(i=0;i<Dim;i++)
                NewPar[i]=0;
            for(i=0;i<=CutDim;i++) // Linearformen in der Auswahl summieren
            {
                j=Simpl1->Vert[Simpl->Vert[i]]; // von Standardsimplex auf Simplex umsetzen
                VP=VLEind(Supp,j);
                for(k=0;k<Dim;k++)
                    NewPar[k]+=VP->V[k];
            }
            if(VinVL(NewPar,*Parallels,Dim))  // Schon vorhanden ?
                continue;
            VP=NewVLE(Parallels);
            memcpy(VP->V,NewPar,Dim*sizeof(Zahl));
        }
    }

    free(NewPar);

}



void CutFaces(VList *NewExtr, VList *NewSupp,VList Extr, VList Supp, int Dim)
// Versucht, glattes Polytop durchg Schneiden mit Hyperebenen zu verkleinern.
// Hyperebenen bestimmt durch Hoehe 1 Punkte ueber der Seite, m.a.W. durch Summe
// der Linearformen, die die Seite ausschneiden. Im Normalenfaecher entspricht
// das zentraler stellarer Unterteilung der dualen Seite.
// Voraussetzung ist, dass alle Ecken ausserhalb der gewaehlten Seite mindestens
// Hoehe 2 darueber haben. Es wird solange geschnitten, bis ein "irreduzibles" Polytop
// erreicht ist.
// Beim Schneideversuch werden die Seiten in absteigender Dimension getestet
// (aufsteigende Dimension im Normalenfaecher).
{

    SimplexList FacesSimp, NormalFan;
    VList Parallels;
    int i;
    Vektor NewFacet;

    NewFacet=calloc(Dim,sizeof(Zahl));

    MakeFacesSimplex(&FacesSimp,Dim-1);
    CopyVL(NewExtr,Extr,Dim);
    CopyVL(NewSupp,Supp,Dim);
    MakeNormalFan(&NormalFan,*NewExtr,*NewSupp,Dim);

    while(1)
    {
        for(i=0;i<Dim-1;i++) // aufsteigende Dimension im Normalenfaecher
        {
            MakeFacesParallels(&Parallels,*NewSupp,NormalFan,Dim,i,FacesSimp);
            if(DoCut(NewExtr,NewSupp,Dim,i,Parallels))
            {
                FreeSimplexL(&NormalFan);
                MakeNormalFan(&NormalFan,*NewExtr,*NewSupp,Dim);
                FreeVL(&Parallels);
                break;
            }
            FreeVL(&Parallels);
        }

        if(i==Dim-1)
            break;
    }

    FreeSimplexL(&FacesSimp);
    FreeSimplexL(&NormalFan);
    free(NewFacet);
}

void CutFacesAvoidForbidden(VList *NewExtr, VList *NewSupp,VList Extr, VList Supp, VList Forbidden,int Dim)
// Versucht, glattes Polytop durchg Schneiden mit Hyperebenen zu verkleinern.
// Hyperebenen bestimmt durch Hoehe 1 Punkte ueber der Seite, m.a.W. durch Summe
// der Linearformen, die die Seite ausschneiden. Im Normalenfaecher entspricht
// das zentraler stellarer Unterteilung der dualen Seite.
// Voraussetzung ist, dass alle Ecken ausserhalb der gewaehlten Seite mindestens
// Hoehe 2 darueber haben. Es wird solange geschnitten, bis ein "irreduzibles" Polytop
// erreicht ist.
// Beim Schneideversuch werden die Seiten in absteigender Dimension getestet
// (aufsteigende Dimension im Normalenfaecher).
// PUNKTE IN FGORBIDDEN WERDEN NICHT BGESCHNITTEN
{

    SimplexList FacesSimp, NormalFan;
    VList Parallels;
    int i;
    Vektor NewFacet;

    NewFacet=calloc(Dim,sizeof(Zahl));

    MakeFacesSimplex(&FacesSimp,Dim-1);
    CopyVL(NewExtr,Extr,Dim);
    CopyVL(NewSupp,Supp,Dim);
    MakeNormalFan(&NormalFan,*NewExtr,*NewSupp,Dim);

    while(1)
    {
        for(i=0;i<Dim-1;i++) // aufsteigende Dimension im Normalenfaecher
        {
            MakeFacesParallels(&Parallels,*NewSupp,NormalFan,Dim,i,FacesSimp);
            if(DoCutAvoidForbidden(NewExtr,NewSupp,Forbidden,Dim,i,Parallels))
            {
                FreeSimplexL(&NormalFan);
                MakeNormalFan(&NormalFan,*NewExtr,*NewSupp,Dim);
                FreeVL(&Parallels);
                break;
            }
            FreeVL(&Parallels);
        }

        if(i==Dim-1)
            break;
    }

    FreeSimplexL(&FacesSimp);
    FreeSimplexL(&NormalFan);
    free(NewFacet);
}



boolean DiagSplitOld(VList Supp, int Dim, int q)
// ERSETZT DURCH neue Funktion DiagSplit, sihe unten
// Testet "diagonal splitting" nach Payne
// Eingabe sind Stuetzhyperebenen des POlytops (ohne rechte Seiten)
// Dim = Poytopdimension, q = zu testender Modul
{
    VList FSupp,Dummy,LatPoints,ConeGen;
    int i,j,k,NResCl;
    boolean *ResCl,IsSplit;
    VLE *VP,*VPNext;

    CopyVL(&Dummy,Supp,Dim);  // Stuetzhyperebenen des Testpolytops
    for(VP=Dummy.St;VP;VP=VP->R)  // zusammenstellen
    {                         // Wir nehmen direkt nur Punkte im Inneren
        VP->V[Dim]=q-1;       // Deshalb q-1 statt q
    }
    CopyVL(&FSupp,Supp,Dim);
    for(VP=FSupp.St;VP;VP=VP->R)
    {
        SxV(VP->V, Dim,-1,VP->V);
        VP->V[Dim]=q-1;
    }
    ConcatVLFree(&FSupp,&Dummy); // Jetzt alle in FSupp gesammelt
    // FSchreibeVL(stdout,FSupp,Dim+1);

    SuppHypVL(&ConeGen,FSupp,Dim+1); // Kegel-Erzeuger bestimmen

    // FSchreibeVL(stdout,ConeGen,Dim+1);

    FindLatPoints(ConeGen,Dim+1, &LatPoints, &Dummy); // Gitterpunkte im (Innern des)
    FreeVL(&Dummy);                                  // Testpolytop
    FreeVL(&ConeGen);

    NResCl=1;                         // Anzahl der Restklassen mod q
    for(i=0;i<Dim;i++)
        NResCl*=q;

    IsSplit=false;
    // printf("NResCl %d IntLP %d\n",NResCl,NinVL(LatPoints));

    if(NResCl>NinVL(LatPoints))       // Nicht genug Gitterpunkte
        goto Ende;

    IsSplit=true;
    ResCl=calloc(NResCl,sizeof(boolean));
    for(VP=LatPoints.St;VP;VP=VP->R)  // Getroffene Restklassen ermitteln
    {
        // SchreibeV(stdout,VP->V,Dim+1);
        k=0;
        for(i=0;i<Dim;i++)      // Restklasse durch q-adische Rntwicklung
        {                       // repräsentiert
            j=VP->V[i] % q;
            if(j<0)
                j+=q;
            k=k*q+j;
        }
        ResCl[k]=true;
    }

    for(k=0;k<NResCl;k++)
        if(!ResCl[k])
        {
            IsSplit=false;
            break;
        }

    free(ResCl);

    Ende:

    FreeVL(&LatPoints);
    FreeVL(&FSupp);
    return IsSplit;
}

boolean NextVinCube(Vektor V, int Dim, int q)
{
    int i,j;

    for(i=Dim-1;i>=0;i--)
        if(V[i]<q)
            break;
    if(i<0)
        return false;

    V[i]++;
    for(j=i+1;j<Dim;j++)
        V[j]=0;
    return true;
}

boolean DiagSplit(VList Supp, int Dim, int q)
// Testet "diagonal splitting" nach Payne
// Eingabe sind Stuetzhyperebenen des POlytops (ohne rechte Seiten)
// WIR SETZEN VORAUS: Die Koordinatenhyperebenen gehören dazu (wird geprüft)
// Dim = Poytopdimension, q = zu testender Modul
// Prinzip: alle Restkallsen mod q bestimmen  und im Polytop suchen
{
    VList FSupp,Dummy;
    int i,j,k,l,NResCl,NInd;
    boolean *ResCl,IsSplit;
    VLE *VP,*VPNext;
    Vektor Test, Test1,Ind;
    Matrix IndV;

    if(MaxDim<Dim+1)
    {
        printf("IsDiagSplit: MaxDim = %d < Dim+1 = %d\n",MaxDim,Dim+1);
        exit(1);
    }

    CopyVL(&Dummy,Supp,Dim);  // Stuetzhyperebenen des Testpolytops
    for(VP=Dummy.St;VP;VP=VP->R)  // zusammenstellen
    {
        VP->V[Dim]=q;
    }
    CopyVL(&FSupp,Supp,Dim);
    for(VP=FSupp.St;VP;VP=VP->R)
    {
        SxV(VP->V, Dim,-1,VP->V);
        VP->V[Dim]=q;
    }
    ConcatVLFree(&FSupp,&Dummy); // Jetzt alle in FSupp gesammelt

    // FSchreibeVL(stdout,FSupp,Dim+1);

    Ind=calloc(Dim+1,sizeof(Zahl)); // Indikator-Vektoren über alle teilmangen
    NInd=1;                            // bilden
    for(i=0;i<Dim;i++) // Mache 2^Dim
        NInd*=2;
    MatAlloc(&IndV,NInd,Dim);
    k=0;
    while(NextVinCube(Ind,Dim,1))
    {
        k++;
        // printf("k %d ",k); SchreibeV(stdout,Ind,Dim);
        memcpy(IndV[k],Ind,Dim*sizeof(Zahl));
    }


    NResCl=1;                         // Anzahl der Restklassen mod q
    for(i=0;i<Dim;i++)
        NResCl*=q;
    // printf("NResCl %d \n",NResCl);
    ResCl=calloc(NResCl,sizeof(boolean));

    Test1=calloc(Dim,sizeof(Zahl)); // repräsentiert Restklasse
    Test=calloc(Dim+1,sizeof(Zahl)); // der eigentliche Testvektor
    Test[Dim]=1; // Kegel-Konvention

    IsSplit=true;
    k=0;

    while(NextVinCube(Test1,Dim,q-1)) // q-adisch hochzählen, Nullvektor braucht
    {                                 // nicht getestet zu werden
        k++;
        // printf(" %d Basis ",k);SchreibeV(stdout,Test1,Dim);
        if(ResCl[k]) // Schon getestet
        {
            // printf("Schon getestet\n");
            continue;
        }

        for(i=0;i<NInd;i++) // gehe über alle Repräsentanten der Restklasse
                            // im ürfel zwischen -(q-1) und q-1
        {
            memcpy(Test,Test1,Dim*sizeof(Zahl));
            for(j=0;j<Dim;j++)
                if(IndV[i][j] && Test[j])
                    Test[j]-=q;

            // printf("Mod   ");SchreibeV(stdout,Test,Dim+1);
            if(InInt(Test,FSupp,Dim+1)) // Restklasse vertreten
            {
                l=0;   // wegen Symmetrie auch -(diese Klasse) vertreten
                       // Müssen q-adischen Index bestimmen
                for(j=0;j<Dim;j++)
                {
                    l*=q;
                    if(Test[j])
                        l+=(q-Test1[j]);
                }

                /* printf("Setze %d\n",l);
                if(l<0 || l>= NResCl)
                    printf("Alarm\n");*/

                ResCl[l]=true;

                break;
            }
        }

        if(i==NInd) // Nicht gefunden
        {
            IsSplit=false;
            break;
        }
    }

    free(Test);
    free(ResCl);
    free(Test1);
    free(Ind);
    FreeMat(IndV,NInd);
    FreeVL(&FSupp);
    return IsSplit;
}


// Grad 2 Checker (fuer beliebeige Gitterpolytope)


typedef struct Deg2Entry{
    int f1;
    int f2;
    VList HigherDeg;
} Deg2Entry;

void MakeDeg2Rel(Deg2Entry*** Deg2M, VList LattPoints, int Dim)
// Sammelt Grad 2 Relationen in Deg2M. Ist Deg2M[i][j].f1 >=0,
// so gibt es die relation x(i)x(j)-x(f1)x(f2), wobei f1 und f2
// die entsprechenden Einträge von Deg2M[i][j] sind. Dabei ist
// x(i)x(j) das Leitmonom in der betrachteten monomialen Ordnung.
// Spaeter kann in HigherDeg eine Vektorliste angehaengt werden
// die Relationen hoeheren Grades repraesentiert, deren Leitmonom
// mit x(i)x(j) anfängt.
// Dim ist hier Kegeldimension.
{
    VList Standard;
    int i,  m1,m2,Dim1;
    VLE *VP1,*VP2,*VPS,*VPN;
    Vektor Sum;

    Dim1=Dim-1;

    Sum=calloc(Dim1,sizeof(Zahl));
    StartVL(&Standard);

    m2=NinVL(LattPoints);

    *Deg2M=calloc(m2,sizeof(Deg2Entry*));
    for(i=0;i<m2;i++)
    {
        (*Deg2M)[i]=calloc(i+1,sizeof(Deg2Entry));
    }
    for(VP2=LattPoints.E;VP2;VP2=VP2->L)
    {
        m2--;
        // printf("m2 %d\n",m2);
        m1=m2+1;
        for(VP1=VP2;VP1;VP1=VP1->L)
        {
            m1--;
            // printf("m1 %d\n",m1);
            for(i=0;i<Dim1;i++)
                Sum[i]=VP1->V[i]+VP2->V[i];
            VPS=VinVL(Sum,Standard,Dim1);
            if(!VPS)
            {
                VPN=NewVLE(&Standard);
                memcpy(VPN->V,Sum,Dim1*sizeof(Zahl));
                VPN->Mark1=m1;
                VPN->Mark2=m2;
                (*Deg2M)[m2][m1].f1=-1;
                StartVL(&((*Deg2M)[m2][m1].HigherDeg));
            }
            else
            {
                // printf("y(%d)*y(%d) - y(%d)*y(%d),\n",m1+1, m2+1, VPS->Mark1+1, VPS->Mark2+1);
                (*Deg2M)[m2][m1].f1=VPS->Mark1;
                (*Deg2M)[m2][m1].f2=VPS->Mark2;
                StartVL(&((*Deg2M)[m2][m1].HigherDeg));
            }
        }
    }
    FreeVL(&Standard);
}

void Mult3(Indizes u, int i, int j, int k)
// Ordnet i,j,k aufsteigend in u. ANNAHME: j <=k
// Interpretation: Das Monom x(j)x(k) wird mit x(i)
// multipliziert
{
    if(i<=j)
    {
        u[0]=i;u[1]=j;u[2]=k;
        return;
    }
    if(i<=k)
    {
        u[0]=j;u[1]=i;u[2]=k;
        return;
    }
    u[0]=j;u[1]=k;u[2]=i;
    return;
}

inline void VertauscheInt(int *u, int *v)
{
  register int Dummy;

  Dummy = *u;
  *u = *v;
  *v = Dummy;
}


void Swap3(Indizes u,Indizes v)
// Vertauscht inhaltlich 2 int-Vektoren der Laenge 3
{
    VertauscheInt(&(u[0]),&(v[0]));
    VertauscheInt(&(u[1]),&(v[1]));
    VertauscheInt(&(u[2]),&(v[2]));
}

void Order3(Indizes u, Indizes v)
// Ordnet das Binom "u-v" um, so dass Leitterm vorn steht
{
    if(u[2]<v[2])
        return;
    if(u[2]>v[2])
    {
        Swap3(u,v);
        return;
    }
    if(u[1]<v[1])
        return;
    if(u[1]>v[1])
    {
        Swap3(u,v);
        return;
    }
    if(u[0]<v[0])
        return;
    if(u[0]>v[0])
    {
        Swap3(u,v);
        return;
    }
}

void SPol3(Indizes u, Indizes v, Deg2Entry **Deg2M, int i, int j, int k, int l, int Pos)
//Bildet S-Polynom u-v zu den Grad 2 Relationen, die durch die Eintraaege bei
// (j,i) und (l,k) in Deg2M definiert sind.
// Pos==1 bedeutet i==k, Pos==2 bedeutet j==l, Pos 3 bedeutet j==k
{
    int g,h;

    //printf("SPol von %d %d - %d %d und %d %d - %d %d\n",i,j, Deg2M[j][i].f1,Deg2M[j][i].f2,
    //                          k,l,Deg2M[l][k].f1,Deg2M[l][k].f2);

    if(Pos==1)
    {
        g=j;
        h=l;
    }
    if(Pos==2)
    {
        g=i;
        h=k;
    }
    if(Pos==3)
    {
        g=i;
        h=l;
    }

    Mult3(u,h,Deg2M[j][i].f1,Deg2M[j][i].f2);
    Mult3(v,g,Deg2M[l][k].f1,Deg2M[l][k].f2);
    Order3(u,v);
}

void Reduce3(Deg2Entry **Deg2M, Indizes u,Indizes v, int *Ini3)
// reduziert Binome des Grades 3 und fuegt irreduzible in Deg2M ein.
{

    int g[3];
    Deg2Entry Cand, *CandP;
    boolean DreiTest;
    VLE *VP;

    // printf("Blabla %d\n",Deg2M[5][0].HigherDeg.St);

    // printf("red %d %d %d - %d %d %d\n",u[0],u[1],u[2],v[0],v[1],v[2]);
    // printf("Ini3 %d\n",*Ini3);

    if(!memcmp(u,v,3*sizeof(int))) // u==v, also u-v==0
    {
        // printf("Gleich\n");
        return;
    }

    Cand=Deg2M[u[1]][u[0]];
    DreiTest=false;
    if(Cand.f1<0)
    {
        /*if(u[1]==5 && u[0]==0)
        {
            printf("Drin\n");
            FSchreibeVL(stdout,Cand.HigherDeg,3);
            printf("$$$$$$$$$$$$$$$\n");
        }*/
        if(Cand.HigherDeg.St)
        {
            DreiTest=true;
            // printf("Setze Dreitest %d %d\n",u[0],u[1]);
        }
    }
    else
    {
        Mult3(g,u[2],Cand.f1,Cand.f2);
        Order3(g,v);
        Reduce3(Deg2M,g,v,Ini3);
        return;
    }

    Cand=Deg2M[u[2]][u[0]];
    if(Cand.f1>=0)
    {
        Mult3(g,u[1],Cand.f1,Cand.f2);
        Order3(g,v);
        Reduce3(Deg2M,g,v,Ini3);
        return;
    }

    Cand=Deg2M[u[2]][u[1]];
    if(Cand.f1>=0)
    {
        Mult3(g,u[0],Cand.f1,Cand.f2);
        Order3(g,v);
        Reduce3(Deg2M,g,v,Ini3);
        return;
    }

    if(DreiTest)
    {
        // printf("In DreiTest\n");
        Cand=Deg2M[u[1]][u[0]];
        for(VP=Cand.HigherDeg.St;VP;VP=VP->R)
        {
            if(memcmp(u,VP->V,3*sizeof(int)))
                continue;
            memcpy(g,VP->LV,3*sizeof(int));
            Order3(g,v);
            Reduce3(Deg2M,g,v,Ini3);
            return;
        }
    }

    // printf("Irreduzibel %d %d\n",u[1],u[0]);

    CandP=&(Deg2M[u[1]][u[0]]); // irreduzibel
    VP=NewVLE(&(CandP->HigherDeg));
    AllocLVE(VP,3);
    memcpy(VP->V,u,3*sizeof(int));
    memcpy(VP->LV,v,3*sizeof(int));
    (*Ini3)++;
    /* printf("Blabla %d\n",*Ini3);
    FSchreibeVL(stdout,CandP->HigherDeg,3);
    printf("Blabla %d %d\n",&(CandP->HigherDeg.St),Cand.HigherDeg.St);
    printf("%%%%%%%%%%%%%%\n"); */

}

void Deg3GB(Deg2Entry **Deg2M, int NLP, long long *HF2, long long *HF3)
// Berechnet Gröbner-Basis bis Grad 3, ausgehend
// von der in Deg2M kodierten GB im Grad 2.
// HF2 bestimmt dim_K (R/I)_2
// HF3 bestimmt dim_K (R/I)_3
//
{
    int i,j,k,l, Ini3;
    int u[3],v[3];

    *HF2=NLP*(NLP+1)/2;
    *HF3=*HF2*(NLP+2)/3;
    Ini3=0;  // zaehlt Initialmon im Grad 3

    for(j=0;j<NLP;j++)
        for(i=0;i<=j;i++)
    {
        // printf("----------------------\n");
        // printf("**** j %d i %d\n",j,i);
        if(Deg2M[j][i].f1<0)
            continue;

        (*HF2)--;  // Neues Initialmonom im Grad 2
        *HF3-=NLP; // Max Beitrag im Grad 3

        for(k=i+1;k<=j;k++)
        {
            // printf("+++ k %d\n",k);
            if(Deg2M[j][k].f1<0)
                continue;
            SPol3(u,v,Deg2M,i,j,k,j,2);
            Reduce3(Deg2M,u,v,&Ini3);
            if(Deg2M[k][i].f1<0 || k==j)  // Kompensation für doppeltes Erzeugen von ijk.
                (*HF3)++;        // Wenn drittes ijk teilendes Monom ebenfalls
                                  // Initialmonom ist, keine Kompensation hier.
            // printf("Zurück 1\n");
        }
        for(l=j+1;l<NLP;l++)
        {
            // printf("/// l %d\n",l);
            if(Deg2M[l][i].f1<0)
                continue;
            SPol3(u,v,Deg2M,i,j,i,l,1);
            Reduce3(Deg2M,u,v,&Ini3);
            (*HF3)++;             // Kompensation für doppeltes Erzeugen von ijk
            // printf("Zurück 2\n");
        }

        if(i==j) continue;

        for(l=j+1;l<NLP;l++)
        {
            // printf("--- l %d\n",l);
            if(Deg2M[l][j].f1<0)
                continue;
            SPol3(u,v,Deg2M,i,j,j,l,3);
            Reduce3(Deg2M,u,v,&Ini3);
            (*HF3)++;            // Kompensation für doppeltes Erzeugen von ijk
            // printf("Zurück 3\n");
        }
    }

    *HF3-=Ini3;
}

void HilbF3(int *HF, Deg2Entry **Deg2M, int Dim, int NLP, int Ini3)
// Ermittelt h-Vektor bis t^3 ausgehend von dem in Deg2M repraesentierten
// Leitmonomial-Ideal
// WIRD NICHT MEHR GEBRAUCHT, DA ZAEHLUNG SCHON IN Deg3GB STATTFINDET
{
    int i,j,k,m,n;
    int Ak, v[3],r,s;
    long long H2, H3;
    VList Vhd;

    StartVL(&Vhd);

    H2=0;
    H3=0;

    for(i=0;i<NLP;i++)
    {
        k=0;
        for(j=i;j<NLP;j++)
            if(Deg2M[j][i].f1>=0)
            {
                H2++;

                n=0;
                for(m=0;m<i;m++)
                    if(Deg2M[j][m].f1>=0 || Deg2M[i][m].f1>=0)
                        n++;

                H3=H3+NLP-k-n;

                /*s=0;
                for(r=0;r<NLP;r++)
                {
                    Mult3(v,r,i,j);
                    if(InsertV(&Vhd,v,3))
                        s++;
                }

                if(NLP-k-n!=s)
                {
                    printf("\nNLP-k %d s %d\n",NLP-k-n,s);
                    printf("j %d i %d k %d\n",j,i,k);
                }*/

                k++;

            }
    }

    HF[0]=1;
    HF[1]=NLP;
    HF[2]=NLP*(NLP+1)/2-H2;
    HF[3]=NLP*(NLP+1)*(NLP+2)/6-H3-Ini3;
}

void NmzHV(int *HV, VList LattP, int Dim, char *FileName)
// Besorgt den von Norm64 bestimmten h-Vektor bis t^3.
// Er wird in HV zurückgegeben
{
    FILE *Aus, *Ein;
    char key[1000];
    char FileNameNmz[1000],FullFileName[100];
    int i,Test;

    strncpy(FileNameNmz,FileName,950);
    strcat(FileNameNmz,"Nmz");

    strcpy(FullFileName,FileNameNmz);
    strcat(FullFileName,".in");

    Aus=fopen(FullFileName,"w");
    FSchreibeVL(Aus,LattP,Dim);
    fprintf(Aus,"cone\n");
    fclose(Aus);

    run_pgm("normaliz -pf -x=6",FileNameNmz," > /dev/null","");

    strcpy(FullFileName,FileNameNmz);
    strcat(FullFileName,".inv");

    Ein=fopen(FullFileName,"r");

    int num_components;

    while(fscanf(Ein,"%s", key)>=1)
    {
        if(strncmp(key, "vector",6))
            continue;
        fscanf(Ein, "%d", num_components);
        fscanf(Ein,"%s", key);
        if(strncmp(key, "hilbert_series_num", 18 ))
            continue;
        fscanf(Ein,"%s", key); // = ueberlesen
        fscanf(Ein,"%d",&HV[0]);
        fscanf(Ein,"%d",&HV[1]);
        if(num_components >= 3)
            fscanf(Ein,"%d",&HV[2]);
        else
            HV[2] = 0;

        if(Dim > 3){
            if(num_components >= 4)
                fscanf(Ein,"%d",&HV[3]);
            else HV[3] = 0;
        }
        // printf("\nHNmz0 %d HNmz1 %d HNmz2 %d HNmz3 %d\n",HV[0],HV[1],HV[2],HV[3]);
        break;
    }
    
    while(fscanf(Ein,"%s", key)>=1)
    {
        if(strncmp(key, "hilbert_polynomial",18)) continue;
         fscanf(Ein,"%s", key); // = ueberlesen
        // printf("Found Ehrhart\n");
        for(i=0;i<Dim;i++)
        {
            fscanf(Ein,"%d",&Test);
            // printf("%d ",Test);
            if(Test<0)
            {
                printf("Non-positive Ehrhart poynomial\n");
                TextMail("wbruns@uos.de","Non-positive Ehrhart");
                VLMail("wbruns@uos.de",LattP,Dim);
                break;
            }
        }
        // printf("\n");
        break;
    }
            
    fclose(Ein);

    strcpy(FullFileName,FileNameNmz);
    strcat(FullFileName,"*");
    run_pgm("rm",FullFileName,"","");

}

boolean Deg2GenDeg3(VList LattPoints, int Dim, int* HVNmz)
{
    Deg2Entry **Deg2M;
    int  NLattPoints,i,j;
    long long HF2,HF3, HF[4],HV[4];

    NLattPoints=NinVL(LattPoints);

    MakeDeg2Rel(&Deg2M,LattPoints,Dim);

    // printf("Deg 2 gefundeden\n");

    Deg3GB(Deg2M,NLattPoints,&HF2,&HF3);

    // printf("GB gefundeden\n");


    for(i=0;i<NLattPoints;i++)
    {
        for(j=0;j<=i;j++)
            if(!Deg2M[i][j].HigherDeg.St)
                FreeVL(&(Deg2M[i][j].HigherDeg));
        free(Deg2M[i]);
    }
    free(Deg2M);

    HF[0]=1;
    HF[1]=NLattPoints;
    HF[2]=HF2;
    HF[3]=HF3;

    memcpy(HV,HF,4*sizeof(long long)); //h-Vektor aus Hilbert-Funkt.
    for(i=0;i<Dim;i++)
    {
        for(j=3;j>=1;j--)
            HV[j]=HV[j]-HV[j-1];
    }

    printf("h-vector GB3 %lld %lld %lld %lld\n",HV[0],HV[1],HV[2],HV[3]);

    // NmzHV(HVNmz, LattPoints,Dim,FileName); Now in ChPoly

    printf("h-vector Nmz %d %d %d %d\n",HVNmz[0],HVNmz[1],HVNmz[2],HVNmz[3]);

    for(j=0;j<=3;j++)
        if(HV[j]!=HVNmz[j])
        {
            if(j<3)
            {
                printf("ALARM\n");
                TextMail("wbruns@uos.de","HILB ALARM");
                exit(1);
            }
            else
                return(false);
        }
    return(true);
}

boolean ConnectedDiv1(VList LP, VList Supp, int Dim, Vektor z, int NBel_z)
// Prueft Zusammenhang des Teilergraphen (im Sinne von "squarefree divisor complex")
// von z.
// Beachte: z ist Wertvektor
// x und y unter z sind durch Kante verbunden, wenn x+y <= z
// Welche x unter z liegen, ist in SyzBors schon festgestellt und durch Mark==1
// gekennzeichnet.
// Ihre Zahl ist NBel_z
// ACHTUNG: Funktion setzt Vorbereitung wie in SyzBors voraus.
// Wenn nicht von SyzBors aufgerufen, muss diese erst erfolgen.
{

    int NConn_z,i,k,l,MaxDist,NSupp;
    Vektor w,w1;
    boolean Connected;
    VLE *VP, **Conn_z; // In Conn_z die registriert, die schon mit Sartpunkt verbunden

    if(NBel_z==1)
        return true;

    NSupp=NinVL(Supp);

    w= (Vektor) calloc(NSupp,sizeof(Zahl));
    w1= (Vektor) calloc(NSupp,sizeof(Zahl));
    Conn_z=calloc(NBel_z,sizeof(VLE*));


    for(VP=LP.St;VP;VP=VP->R) // suche Startpunkt
        if(VP->Mark)
        {
            Conn_z[0]=VP;
            VP->Mark1=0;
            VP->Mark2=0; // fuer Statistik
            break;
        }

    NConn_z=1; // zaehlt die mit Startpunkt verbundenen

    Connected=false;

    k=0;l=0;MaxDist=0;  // Statistik

    for(i=0;i<NConn_z;i++) // Wir versuchen and die schon verbundenen neue Punte
    {                      // anzuhaengen. Dabei waechst die Liste der verbundenen
        k++;
        SubtVekt(w,z,Conn_z[i]->LV,NSupp);  // w=z-x

        for(VP=Conn_z[0];VP;VP=VP->R) // versuche noch nicht verbundene Punkte
            if(VP->Mark && VP->Mark1) // an den schon verbundenen Conn_z[i] anzuhaengen
            {
                l++;

                if(NonNegDiff(w1,w,VP->LV,NSupp)) // bilde w1=w-y=z-(x+y),falls>=0 x und y verbunden
                {
                    VP->Mark1=0;                  // als verbunden kennzeichnen
                    VP->Mark2=Conn_z[i]->Mark2+1; // Nur fuer Ermittlung Maximaldistanz
                    if(VP->Mark2>MaxDist)
                        MaxDist=VP->Mark2;
                    Conn_z[NConn_z]=VP; // an die Liste der schon verbundenen anhaengen
                    NConn_z++;
                    if(NConn_z==NBel_z) // alle verbunden
                    {
                        Connected=true;
                        goto Freimachen;
                    }
                }
            }
    }

    Freimachen:

    free(w);
    free(w1);
    free(Conn_z);

    // printf("%d Runden, %d Tests, %d MaxDist\n",k,l,MaxDist);

    return Connected;
}

boolean SyzBors(VList LP, VList Supp, int Dim, int NLP3, int NTS)
// Fuehrt "SyzBor-Test" aus
// LP enthaelt Gitterpunkte des Polytops in C(P)
// Supp enthaelt Stuetzhyberebenen von C(P)
// Dim ist Kegeldimension
// NLP3 sagt, wie viele Punkte in 3P getestet werden sollen
// NTS sagt, wie viele Paare auf SyzBor-Eigenschaft
// fuer jedes z in 3P getestet werden sollen.
{
    int i,j,k,l,m,n,NLP1,NSBz,NBel_zTot,NSyzBorTot,NNotSyzBor,NSupp;

    float AverageBel_z, AverageSyzBor;

    Vektor z,w,x,y,w1,zV;
    VLE *VP, **LPReg;
    boolean Conn;

    Conn=true;

    NSupp=NinVL(Supp);
    z= (Vektor) calloc(NSupp,sizeof(Zahl));
    w= (Vektor) calloc(NSupp,sizeof(Zahl));
    w1= (Vektor) calloc(NSupp,sizeof(Zahl));

    zV= (Vektor) calloc(Dim,sizeof(Zahl));

    Registriere(LP,&LPReg,&NLP1); // fuer wahlfreien Zugriff


    AddValLV(LP,Supp,Dim); // Fuegen Wert-Vektoren hinzu

    NBel_zTot=0;
    NSyzBorTot=0;

    for(i=0;i<NLP3;i++)
    {
        // printf("i %d\n",i);
        MacheNullVekt(z,NSupp);    // wird Punkt in 3P
        m=rand()%NLP1;
        n=rand()%NLP1;
        l=rand()%NLP1;
        AddVekt(z,LPReg[m]->LV,LPReg[n]->LV,NSupp);
        AddVekt(z,z,LPReg[l]->LV,NSupp);
        AddVekt(zV,LPReg[m]->V,LPReg[n]->V,Dim); // der Vektor zum Wertvektor z
        AddVekt(zV,zV,LPReg[l]->V,Dim);
        // printf("zV ");SchreibeV(stdout,zV,Dim);

        k=0;

        NNotSyzBor=0;

        for(VP=LP.St;VP;VP=VP->R)
        {

            if(NonNegDiff(w,z,VP->LV,NSupp))
            {
                VP->Mark=1;
                VP->Mark1=1; // fuer ConnectedDiv1 setzen
                k++;
            }
            else
                VP->Mark=0;
        }

        NBel_zTot+=k;
        // printf("NBel_zTot %d\n",NBel_zTot);
        NSBz=0;

        if(k==1) // Nur 1 Punkt unter z
        {
            NSyzBorTot+=NTS; // wird voll gezaehlt
            continue;
        }

        for(l=0;l<NTS;l++)
        {
            // printf("llll %d\n",l);

           while(1)  //Zufallspunkt finden, der unter z liegt
           {
               m=rand()%NLP1;
               if(LPReg[m]->Mark)
                  break;
           }
           while(1) // den zweiten davon verschieden
           {
                n=rand()%NLP1;
                if(m!=n && LPReg[n]->Mark)
                    break;
           }
           x=LPReg[m]->LV;
           y=LPReg[n]->LV;

           // printf("mmm %d\n",l);

           AddVekt(w,x,y,NSupp);  // Ausprobieren ob sogar x+y unter z
           if(NonNegDiff(w1,z,w,NSupp))
           {
               NSBz++;
               continue;
           }

           for(VP=LP.St;VP;VP=VP->R) // jetzt gemeinsamen oberen Nachbarn von
                                     // x und y in 2P suchen, der unter z liegt
           {
                if(!VP->Mark) // liegt nicht unter z
                    continue;
               SubtVekt(w,z,VP->LV,NSupp); // w=(z-dieser Vektor)

               if(NonNegDiff(w1,w,x,NSupp))  // liegt x unter w ?
               {
                   if(NonNegDiff(w1,w,y,NSupp))  // liegt auch noch y darunter ??
                       break;  // ja, sie sind SyzBors
               }
           }
           if(VP)
               NSBz++;
           else
           {
                // printf("Not SyzBors\n");
                NNotSyzBor++;
           }
        }

        NSyzBorTot+=NSBz;

        if((float) NNotSyzBor > 0.15*NTS)
        {
            // printf("%d suspicious\n",i);
            // SchreibeV(stdout,zV,Dim);
            if(!ConnectedDiv1(LP,Supp,Dim,z,k))
            {
                SchreibeV(stdout,zV,Dim);
                printf("Counterexample: a syz graph not connected\n");
                Conn=false;
            }
        }
    }

    /* AverageBel_z=(float) NBel_zTot / (float) NLP3;

    AverageSyzBor=(float) NSyzBorTot / ((float) NLP3*NTS);

    printf("%d points z in 3P tested, %d pairs (x,y) each\n",NLP3,NTS);
    printf("%f average points below z of %d points in P, %f of total pairs are syzbors\n",
                  AverageBel_z,NLP1,AverageSyzBor); */

    free(z);
    free(w);
    free(w1);
    free(zV);
    free(LPReg);
    StripLV(LP);
    return Conn;
}


// Optionsverwaltung

void ExtractNameValue(char **Name, char **Value, char *Line)
//extracts the first (at most) two tokens, delimited by whitespace, from Line
{
        *Name=strtok(Line," \r\t\n");
        *Value=strtok(NULL," \r\t\n");
}

int SearchOption(char *OptName, char **Options)
// searches option OptName in Options and returns its index i
{
    int i;
    char *Name,*Val,Line[100];
    for(i=0;;i++)
    {
        strncpy(Line,Options[i],100);
        ExtractNameValue(&Name,&Val,Line);
        if(!strcmp(Name,"end_options"))
        {
            printf("SearchOption: invalid option name %s\n",OptName);
            exit(1);
        }
        if(Name==NULL || strcmp(Name,OptName))
            continue;
        return i;
    }
}

void SetOption(char *Value, char *OptName, char **Options)
// sets Value of option OptName in Options
{
    int i;
    if(OptName==NULL || Value==NULL)
    {
        printf("SetOption: option with empty name or value\n");
        exit(1);
    }
    i=SearchOption(OptName,Options);
    sprintf(Options[i],"%s %s",OptName,Value);
}

void GetOption(char *Value, char *OptName, char **Options)
// retrieves Value of option OptName in Options
{
    char *Name,*Val,Line[100];
    int i;

    i=SearchOption(OptName,Options);
    strncpy(Line,Options[i],100);
    ExtractNameValue(&Name,&Val,Line);
    if(Val==NULL)
    {
        printf("GetOptValue: option %s has no value\n",Name);
        exit(1);
    }
    strcpy(Value,Val);
}

void ReadConfigFile(char *ProgName, char **Options)
// reads file "ProgName.cfg" and interprets the line
// as names and values of options in Options
// empty lines and comment lines allowed
// Warning: the Progname contained in argv[0] may depend on the system
{
    FILE *cfg;
    char FullName[100],Line[100],*OptName,*OptValue,*Test;

    sprintf(FullName,"%s.%s",ProgName,"cfg");
    cfg=fopen(FullName,"r");
    if(!cfg)                    // no configuration file
        return;

    while(1)
    {
        Test=fgets(Line,100,cfg);
        if(!Test)  // end of file reached
            break;
        ExtractNameValue(&OptName,&OptValue,Line);
        if(OptName==NULL)                                         // empty line
            continue;
        if(strlen(OptName)>=2 && OptName[0]=='/' && OptName[1]=='/') // comment line
            continue;
        if(OptValue==NULL)
        {
            printf("ReadConfigFile: option %s has no value\n",OptName);
            exit(1);
        }
        SetOption(OptValue,OptName,Options);
    }
}


void StringOption(char *Value, char *OptName, char **Options)
// searches the option Name in Options and returns value as a string in Value
{
    strcpy(Value,"");
    GetOption(Value,OptName,Options);
    if(strcmp(Value,"(empty)")==0)
    {
        strcpy(Value,"");
    }
}

void IntOption(int *Value, char *Name, char **Options)
// ditto, but returns value as an int
{
    char StringVal[100];
    StringOption(StringVal,Name,Options);
    if(strcmp(StringVal,"")==0)
    {
        printf("IntOption: empty string\n");
        exit(1);
    }
    sscanf(StringVal,"%d",Value);
}

void BoolOption(boolean *Value, char *Name, char **Options)
// ditto, but returns value as a boolean
{
    char StringVal[100];

    StringOption(StringVal,Name,Options);
    // printf("Val %d %s\n",strlen(StringVal),StringVal);
    if(!strcmp(StringVal,"true"))
    {
        *Value=true;
        return;
    }
    if(!strcmp(StringVal,"false"))
    {
        *Value=false;
        return;
    }
    printf("BoolOption: illegal value %s of option %s\n",StringVal,Name);
        exit(1);
}

void TransferDefOptions(char ***Options, char **DefOptions, int NumOptions)
// transfers default options from DefOptions to Options
// space for options is allocated
// NumOptions is he number of options = #(lines on DefOptions)
{
    int i;
    char Line[100],*Name,*Val;

    *Options=calloc(NumOptions,sizeof(char*));

    for(i=0;i<NumOptions;i++)
    {
        // printf("Reserviere %d %d\n",i,strlen(DefOptions[i]));
        (*Options)[i]=calloc(100,sizeof(char));
        strncpy((*Options)[i],DefOptions[i],100);
        strcpy(Line,(*Options)[i]);
        // printf("%s \n",Line);
        ExtractNameValue(&Name,&Val,Line);
        if(Name==NULL || (i!=NumOptions-1 && Val==NULL))
        {
            printf("TransferDefOptions: empty name or value in deafault options\n");
            exit(1);
        }
    }

    strcpy(Line,(*Options)[NumOptions-1]);
    ExtractNameValue(&Name,&Val,Line);
    if(strcmp(Name,"end_options"))
    {
        printf("TransferDefOptions: default options do not end with end_options\n");
        exit(1);
    }
}

// ShrinkHsp.c hier eingebaut

void MoveChain(VList *Target, VList *Source, VLE *From, VLE* To)
// Bewegt Kette von VLEs von Source nach Target
// Es wird an das Ende von Target angehaengt
//
{
  if (Source->St == From)
    Source->St = To->R;
  else
    From->L->R = To->R;
  if (Source->E == To)
    Source->E = From->L;
  else
    To->R->L = From->L;
  To->R = NULL;
  if (Target->E != NULL) {
    Target->E->R = From;
    From->L = Target->E;
  } else
    From->L = NULL;
  Target->E = To;
  if (Target->St == NULL)
    Target->St = From;
}


int VLEVerglLV(const void *m,const void *n)
// VLEs werden verglichen mittels Totalgrad (in Mark) und LVs
//
{
    int i;
    Vektor Vm,Vn;
    Zahl Tm,Tn;

    Vm=(*(VLE **)m)->LV;
    // printf("m %x Vm %x\n",m,Vm);
    //printf("*** ");SchreibeV(stdout,Vm,VEKTLAENGE);
    Vn=(*(VLE **)n)->LV;
    //printf("/// ");SchreibeV(stdout,Vn,VEKTLAENGE);
    Tm=(*(VLE **)m)->Mark;
    Tn=(*(VLE **)n)->Mark;

    if(Tm>Tn)
        return 1;
    if(Tn>Tm)
        return -1;

    for(i=0;i<VEKTLAENGE;i++)
    {
        if(Vm[i]==Vn[i])
            continue;
        if(Vm[i]>Vn[i])
            return 1;
        return -1;
    }
    return 0;
}


void RemDoubleLV(VList *Source, int NSupp)
{
    VLE *VP, *VPNext;

    if(!Source->St) return; // Liste leer

    for(VP=Source->St;;VP=VPNext)
    {
        VPNext=VP->R;
        if(!VPNext) break;
        if(!memcmp(VP->LV,VPNext->LV,NSupp*sizeof(Zahl)))
            FreeVLE(Source,VP);
    }
}

void OrderVListLV(VList *Source, int NSupp, VLE ***Kand, int *NKand)
// im Falle NSupp<0 keine Doubletten entfernen
{

    VLE *VP,*From,*To;
    VList Sorted;
    int i,NKand0,NKand1;
    VLE **Order;


    Registriere(*Source, &Order,&NKand0);

    if(!Source->St) // Liste leer
    {
        *NKand=0;
        *Kand=Order;
        return;
    }

    NKand1=NKand0-1;

    qsort(Order,NKand0,sizeof(VLE*),VLEVerglLV);
    StartVL(&Sorted);

    for(i=0;i<NKand0;i++)
    {
        From=Order[i];
        To=Order[i];
        while(i<NKand1 && Order[i+1]==To->R)
        {
            i++;
            To=Order[i];
        }
        MoveChain(&Sorted,Source,From,To);
    }
    free(Order);
    *Source=Sorted;
    if (NSupp>=0)
        RemDoubleLV(Source,NSupp);
    Registriere(*Source, Kand, NKand);
}

boolean Reduzibel(Zahl TGrad,Vektor Val, VList *Red, int NSupp,
                        boolean Umordnen, boolean Sortiert)
{
    int j,Gr;
    VLE *VP;
    Gr=0;


    for(VP=Red->St;VP;VP=VP->R)
    {
        if(Sortiert && 2*VP->Mark>TGrad)
           break;

        if(VP->Mark >= TGrad)
           continue;

        if(Val[Gr]<VP->LV[Gr])
           continue;

        for(j=0;j<NSupp;j++)
        {
            if(Val[j]<VP->LV[j])
                break;
        }
        Gr=j;
        if(j==NSupp) // reduziert
        {
             if(Umordnen)
               MakeFirst(Red,VP);
            return true;
        }
    }
    return false;
}

void RedVLvsVL(VList *Source, VList *Red, int NSupp)
{
    VLE *VP, *VPNext,**Reg,**Dummy;
    int NSource,i,NDummy;

    OrderVListLV(Source,NSupp,&Reg,&NSource);
    if(Source!=Red)  // Keine Autoreduktion
        OrderVListLV(Red,NSupp,&Dummy,&NDummy);

    for(i=0;i<NSource;i++)
    {
        VP=Reg[i];
        if(!VP->Mark2) continue;  // absolut irreduzibel
        if(Reduzibel(VP->Mark,VP->LV,Red,NSupp,true,true))
            FreeVLE(Source,VP);
    }
    free(Reg);
    if(Source!=Red)
        free(Dummy);
}

void AutoReduceVL(VList *Source, int NSupp)
{
    RedVLvsVL(Source,Source,NSupp);
}

boolean ReduzibelExtr(Zahl TGrad,Vektor Val, VList *Red, int NSupp)
{
    int j,Gr;
    VLE *VP;
    Gr=0;


    for(VP=Red->St;VP;VP=VP->R)
    {

        if(VP->Mark >= TGrad)
           break;

        if(Val[Gr]<VP->LV[Gr])
           continue;

        for(j=0;j<NSupp;j++)
        {
            if(Val[j]<VP->LV[j])
                break;
        }
        Gr=j;
        if(j==NSupp) // reduziert
        {
            MakeFirst(Red,VP);
            return true;
        }
    }
    return false;
}

boolean ReduzibelExtr1(Zahl TGrad,Vektor Val, VList *Red, int NSupp)
{
    int j,k;
    VLE *VP;
    Indizes Z;

    Z=calloc(NSupp,sizeof(int));
    k=0;
    for(j=0;j<NSupp;j++)
        if(!Val[j])
        {
            Z[k]=j;
            k++;
        }

    for(VP=Red->St;VP;VP=VP->R)
    {

        if(VP->Mark >= TGrad)
           break;

        for(j=0;j<k;j++)
        {
            if(VP->LV[Z[j]])
                break;
        }
        if(j==k) // reduziert
        {
            MakeFirst(Red,VP);
            free(Z);
            return true;
        }
    }
    free(Z);
    return false;
}

void RedVLvsVLExtr(VList *Extr, int NSupp)
// Fuer Bestimmung von Extremstrahlen
{
    VLE *VP, *VPNext,**Reg;
    int NSource,i;

    OrderVListLV(Extr,NSupp, &Reg,&NSource);

    for(i=0;i<NSource;i++)
    {
        VP=Reg[i];
        if(Reduzibel(VP->Mark,VP->LV,Extr,NSupp,true,false))

        {
            FreeVLE(Extr,VP);
        }
    }
    free(Reg);
}

void RedExtr(VList *Kand, int NSupp)
// Fuer Bestimmung von Extremstrahlen
// mittels Reduktion von 0-1-Vektoren
{
    VLE *VP, *VPNext,**Reg;
    int NSource,i;
    VList Extr;

    OrderVListLV(Kand,NSupp, &Reg,&NSource);
    StartVL(&Extr);

    for(i=0;i<NSource;i++)
    {
        VP=Reg[i];
        if(ReduzibelExtr1(VP->Mark,VP->LV,&Extr,NSupp))
            FreeVLE(Kand,VP);
        else
            MoveVLE(&Extr,Kand,VP);
    }
    free(Reg);
    *Kand=Extr;
}



void SelExtremeR1(VList *Extr, VList Gen, VList Supp, int Dim)
// Waehlt Extremstrahlen aus, veraendert Gen nicht
{
    VLE *VP, *VPNext;
    int i,k,NSupp, MaxOne;

    CopyVL(Extr,Gen,Dim);

    AddValLV(*Extr,Supp,Dim);
    NSupp=NinVL(Supp);
    MaxOne=NSupp-Dim+1;
    for(VP=Extr->St;VP;VP=VPNext)
    {
        VPNext=VP->R;
        k=0;
        for(i=0;i<NSupp;i++)
            if(VP->LV[i])
            {
                // VP->LV[i]=1;
                k++;
            }
        if(k>MaxOne)
        {
            FreeVLE(Extr,VP);
            continue;
        }
        VP->Mark=k;
    }
    VEKTLAENGE=NSupp;
    // OrderVListLV(Extr,NSupp);
    // RedVLvsVLExtr(Extr,NSupp);
    printf("Kand %d\n",NinVL(*Extr));
    RedExtr(Extr,NSupp);
}

void SelExtremeRRank(VList *Extr, VList Gen, VList Supp, int Dim)
// Waehlt Extremstrahlen aus, veraendert Gen nicht
{
    VLE *VP, *VPNext;
    int i,k,NSupp, MaxOne;
    Matrix SuppMat, TestMat;
    Indizes Z;

    NSupp=NinVL(Supp);
    MatAlloc(&SuppMat,NSupp,Dim);
    VL2Mat(SuppMat,Supp,Dim);
    MatAlloc(&TestMat,NSupp,Dim);
    Z=calloc(NSupp,sizeof(int));

    CopyVL(Extr,Gen,Dim);

    AddValLV(*Extr,Supp,Dim);
    MaxOne=NSupp-Dim+1;
    for(VP=Extr->St;VP;VP=VPNext)
    {
        VPNext=VP->R;
        k=0;
        for(i=0;i<NSupp;i++)
            if(!VP->LV[i])
            {
                Z[k]=i;
                k++;
            }
        if(k<Dim-1)
        {
            FreeVLE(Extr,VP);
            continue;
        }
        FuelleMatrix(TestMat,k,Dim,SuppMat,Z);
        if(RankMat0(TestMat,k,Dim)<Dim-1)
            FreeVLE(Extr,VP);
    }
    FreeMat(SuppMat,NSupp);
    FreeMat(TestMat,NSupp);
    free(Z);
}


void Markiere(VLE *VP,Zahl M,Zahl M1, Zahl M2, Zahl M3)
{
    VP->Mark=M;
    VP->Mark1=M1; // Generation (0: Grosseltern (und aelter), 1: Eltern, 2: Kinder
    VP->Mark2=M2; // 0: sicher irreduzibel, 1: vielleicht reduzibel
    VP->Mark3=M3; // Wert des Vaters
}

void TrenneGen(VList Source, VLE ***OldGen,int *AnzOldGen,
                            VLE ***NewGen, int *AnzNewGen)
{
    int i,j;
    VLE *VP;

    *AnzOldGen=*AnzNewGen=0;
    for(VP=Source.St;VP;VP=VP->R)
    {
        if(!VP->Mark1)
            (*AnzOldGen)++;
        else
            (*AnzNewGen)++;
    }

    *OldGen=(VLE**) calloc(*AnzOldGen,sizeof(VLE*));
    *NewGen=(VLE**) calloc(*AnzNewGen,sizeof(VLE*));

    i=j=0;
    for(VP=Source.St;VP;VP=VP->R)
    {
        if(!VP->Mark1)
        {
            (*OldGen)[i]=VP;
            i++;
        }
        else
        {
            (*NewGen)[j]=VP;
            j++;
        }
    }
}

void DivByHypIt0(VList *Hilb, Vektor HspGen, int EmbDim, Vektor Hyp,int IndHyp, int Size,
                VList *NewHilb, boolean C0HSP)
// Bestiimt Hilbertbasis NewHilb fuer den Durchschnitt des von Hilberbasis Hilb
// erzeugten Kegels C mit positivem Halbraum bezueglich Hyp.
// Im Fall dass C0HSP gesetzt ist, verkleinert sich der max. Unterraum C_0
// Die Erzeuger eines der Halbraeume von C_0 bez. Hyp wird ueber HspGen eingebracht.
// IndHyp ist der Index, der Hyp entspricht in den Wertvektoren
// Size ist die Laenge, die die wertvektoren maximal erreichen koennen
{
    VLE *VP,*VP1,*VP2,*VP3;
    VList Pos,Neg,Neutr,PosIrred,NegIrred,NeutrIrred,DummyVL;
    int i,j,T1,pp,nn,NPos,NNeg,AnzNeutrIrred,Round,kk,
          AnzNewGenPos,AnzOldGenPos,AnzNewGenNeg,AnzOldGenNeg,AnzPosIrred,AnzNegIrred;
    Vektor SumLV,SumV, V1,V2;
    Zahl Diff, PosW, NegW;
    boolean NewNeutr, NewPos,NewNeg;
    Zahl Quot,Hoehe,Wert,WN,S;

    VLE *Dummy, **PosRun,**NegRun, **PosIrredReg,**NegIrredReg, **NeutrIrredReg,
                **OldGenPos,**NewGenPos,**OldGenNeg,**NewGenNeg;

    StartVL(&Pos);
    StartVL(&Neg);
    StartVL(&Neutr);
    StartVL(&PosIrred);
    StartVL(&NegIrred);
    StartVL(&NeutrIrred);
    T1=IndHyp+1;
    // SumLV=calloc(T1,sizeof(Zahl));
    // SumV=calloc(EmbDim,sizeof(Zahl));
    StartVL(&DummyVL);
    Dummy=NewVLE(&DummyVL); // fuer das Suchen
    AllocLVE(Dummy,Size);
    SumLV=Dummy->LV;
    SumV=Dummy->V;


    if(C0HSP)
    {
        Hoehe=SkalProd(HspGen,Hyp,EmbDim);
        if(Hoehe<0)
        {
            Hoehe*=-1;
            SxV(HspGen,EmbDim,-1,HspGen); // mache ihn zum Erzeuger des pos. Halbraums
        }

        for(VP=Hilb->St;VP;VP=VP->R)  // modulo Halbraumerzeuger reduzieren
        {
            Wert=SkalProd(VP->V,Hyp,EmbDim);
            WN=1;
            if(Wert<0)
            {
                WN=-1;
                Wert*=-1;
            }
            Quot=WN*(Wert/Hoehe);
            for(i=0;i<EmbDim;i++)
                    VP->V[i]=VP->V[i]-Quot*HspGen[i];
        }
    }

    for(VP=Hilb->St;VP;VP=VP1) // aufteilen in Pos, Neg, Neutr
    {
        VP1=VP->R;
        S=SkalProd(VP->V,Hyp,EmbDim);
        if(S==0)
        {
            MoveVLE(&NeutrIrred,Hilb,VP);
            VP->LV[IndHyp]=0;
        }
        if(S>0)
        {
            MoveVLE(&PosIrred,Hilb,VP);
            VP->LV[IndHyp]=S;
        }
        if(S<0)
        {
            MoveVLE(&NegIrred,Hilb,VP);
            VP->LV[IndHyp]=-S;
        }
        S=0; // Totalgrad bestimmen
        for(i=0;i<T1;i++)
            S+=VP->LV[i];
        Markiere(VP,S,1,0,0);
    }

    if(C0HSP) // Halbraum-Erzeuger in Pos und Neg einsetzen
    {
        VP=NewVLE(&PosIrred);
        memcpy(VP->V,HspGen,EmbDim*sizeof(Zahl));
        AllocLVE(VP,Size);
        for(i=0;i<IndHyp;i++)
            VP->LV[i]=0;
        VP->LV[IndHyp]=Hoehe;
        Markiere(VP,Hoehe,0,0,0);

        VP=NewVLE(&NegIrred);
        SxV(VP->V,EmbDim,-1,HspGen);
        AllocLVE(VP,Size);
        for(i=0;i<IndHyp;i++)
            VP->LV[i]=0;
        VP->LV[IndHyp]=Hoehe;
        Markiere(VP,Hoehe,0,0,0);
    }

    VEKTLAENGE=T1;  // ABSOLUT WICHTIG FUER ALLE VERGLEICHE

    OrderVListLV(&PosIrred,T1,&PosIrredReg,&AnzPosIrred);   // Aufsteigend anordnen
    OrderVListLV(&NegIrred,T1,&NegIrredReg,&AnzNegIrred);
    OrderVListLV(&NeutrIrred,T1,&NeutrIrredReg,&AnzNeutrIrred);

    /* printf("Pos %d Neg %d Null %d\n",
              AnzPosIrred,AnzNegIrred,AnzNeutrIrred); */

    fflush(stdout);

    Round=0;

    while(1)
    {
        Round++;

        TrenneGen(PosIrred,&OldGenPos,&AnzOldGenPos,&NewGenPos,&AnzNewGenPos);
        TrenneGen(NegIrred,&OldGenNeg,&AnzOldGenNeg,&NewGenNeg,&AnzNewGenNeg);

    for(kk=0;kk<3;kk++)
    {
        if(kk==0) // Alt Pos gegen Neu Neg
        {
            PosRun=OldGenPos;
            NPos=AnzOldGenPos;
            NegRun=NewGenNeg;
            NNeg=AnzNewGenNeg;
        }
        if(kk==1) // Neu Pos gegen Neu Neg
        {
            PosRun=NewGenPos;
            NPos=AnzNewGenPos;
        }
        if(kk==2) // Neu Pos gegen Alt Neg
        {
            NegRun=OldGenNeg;
            NNeg=AnzOldGenNeg;
        }

        if(NPos==0 || NNeg==0)
            continue;

        for(pp=0;pp<NPos;pp++)
        {
          VP=PosRun[pp];
          for(nn=0;nn<NNeg;nn++)
          {
              VP1=NegRun[nn];

              V1= VP->LV;
              V2=VP1->LV;
              PosW=V1[IndHyp];
              NegW=V2[IndHyp];

              if( (VP->Mark3 && (NegW>=VP->Mark3)) ||
                               (VP1->Mark3 && (PosW>=VP1->Mark3)))
                  continue; // Addition haette schon Vater auf die
                            // andere Seite gebracht. Deshalb reduzibel

              for(i=0;i<IndHyp;i++)
                  SumLV[i]=V1[i]+V2[i];
              Diff=PosW-NegW;

              if(Diff==0)
              {
                  SumLV[IndHyp]=0;
                  // memcpy(Dummy->LV,SumLV,T1*sizeof(Zahl));
                  S=VP->Mark+VP1->Mark-PosW-NegW; // neuer Totalgrad
                  Dummy->Mark=S;
                  // printf("Dummy %x %x\n",Dummy, Dummy->LV);

                  if(Reduzibel(S,SumLV,&NeutrIrred,IndHyp,true,false))
                    continue; // reduzibel

                  // if(bsearch(&Dummy,NeutrIrredReg,AnzNeutrIrred,sizeof(VLE*),VLEVerglLV))
                  //  continue; // Doublette

                  VP2=NewVLE(&Neutr);
                  Wert=0;
              }
              if(Diff>0)
              {
                  SumLV[IndHyp]=Diff;
                  // memcpy(Dummy->LV,SumLV,T1*sizeof(Zahl));
                  S=VP->Mark+VP1->Mark+Diff-PosW-NegW;
                  Dummy->Mark=S;
                  // printf("Dummy %x %x\n",Dummy, Dummy->LV);

                  if(Reduzibel(S,SumLV,&NeutrIrred,IndHyp,true,false))
                     continue; // reduzibel
                  if(Reduzibel(S,SumLV,&PosIrred,T1,true,false))
                    continue; // reduzibel

                  // if(bsearch(&Dummy,PosIrredReg,AnzPosIrred,sizeof(VLE*),VLEVerglLV))
                  //   continue; // Doublette

                  VP2=NewVLE(&Pos);
                  Wert=PosW;
              }
              if(Diff<0)
              {
                  SumLV[IndHyp]=-Diff;
                  // memcpy(Dummy->LV,SumLV,T1*sizeof(Zahl));
                  S=VP->Mark+VP1->Mark-Diff-PosW-NegW;
                  Dummy->Mark=S;
                  // printf("Dummy %x %x\n",Dummy, Dummy->LV);

                  if(Reduzibel(S,SumLV,&NeutrIrred,IndHyp,true,false))
                       continue; // reduzibel
                  if(Reduzibel(S,SumLV,&NegIrred,T1,true,false))
                    continue; // reduzibel

                  // if(bsearch(&Dummy,NegIrredReg,AnzNegIrred,sizeof(VLE*),VLEVerglLV))
                  //    continue; // Doublette

                  VP2=NewVLE(&Neg);
                  Wert=NegW;
              }

              for(i=0;i<EmbDim;i++)  // Summenvektor einsetzen
                SumV[i]=VP->V[i]+VP1->V[i];
              memcpy(VP2->V,SumV,EmbDim*sizeof(Zahl));
              AllocLVE(VP2,Size);
              memcpy(VP2->LV,SumLV,T1*sizeof(Zahl));
              Markiere(VP2,S,2,1,Wert);



          } // neg Schleife
        } // pos Schleife
    }     // kk Schleife

        if(AnzPosIrred) free(PosIrredReg);
        if(AnzNegIrred) free(NegIrredReg);
        if(AnzNeutrIrred) free(NeutrIrredReg);

        if(AnzOldGenPos) free(OldGenPos);
        if(AnzNewGenPos) free(NewGenPos);
        if(AnzNewGenNeg) free(NewGenNeg);
        if(AnzOldGenNeg) free(OldGenNeg);

        /* printf("Vor Red Pos %d Neg %d Null %d\n",
              NinVL(Pos),NinVL(Neg),NinVL(Neutr));
        printf("Für Schleife ");RunTime(); */

        // Jetzt alle Summen gebildet und in Pos, Neg, Neutr einsortiert

        NewNeutr=(Neutr.St!=NULL);
        NewPos=(Pos.St!=NULL);
        NewNeg=(Neg.St!=NULL);


        if(NewNeutr) AutoReduceVL(&Neutr,IndHyp);  // Neue Vektoren gegen sich red                                                // selbst reduzieren
        if(NewPos) AutoReduceVL(&Pos,T1);
        if(NewNeg) AutoReduceVL(&Neg,T1); //WICHTIG: jetzt alle ohne Doubletten
        if(NewNeutr && NewPos) RedVLvsVL(&Pos,&Neutr,T1);
        if(NewNeutr && NewNeg) RedVLvsVL(&Neg,&Neutr,T1);

        /* printf("NewPos %d NewNeg %d NewNull %d\n",
              NinVL(Pos),NinVL(Neg),NinVL(Neutr));
        printf("Für Red Neu gegen Neu ");RunTime();
        fflush(stdout); */

        NewPos=(Pos.St!=NULL);
        NewNeg=(Neg.St!=NULL);

        if(NewNeutr) // alte gegen neue reduzieren und Listen vereinigen und ordnen
        {
            if(Round>1) RedVLvsVL(&NeutrIrred,&Neutr,IndHyp);
            if(Round>1) RedVLvsVL(&PosIrred,&Neutr,T1);
            if(Round>1) RedVLvsVL(&NegIrred,&Neutr,T1);
            ConcatVLFree(&NeutrIrred,&Neutr);
        }
        OrderVListLV(&NeutrIrred,IndHyp,&NeutrIrredReg,&AnzNeutrIrred);

        if(NewPos)
        {
            if(Round>1) RedVLvsVL(&PosIrred,&Pos,T1);
            ConcatVLFree(&PosIrred,&Pos);
        }
        OrderVListLV(&PosIrred,T1,&PosIrredReg,&AnzPosIrred);

        if(NewNeg)
        {
            if(Round>1) RedVLvsVL(&NegIrred,&Neg,T1);
            ConcatVLFree(&NegIrred,&Neg);
        }
        OrderVListLV(&NegIrred,T1,&NegIrredReg,&AnzNegIrred);

        /* printf("Für Red Alt gegen Neu ");RunTime();
        printf("--------------------\n"); */

        if(!NewPos && !NewNeg)
            break;

        // printf("Vor Generationen angepasst\n");

        for(VP=PosIrred.St;VP;VP=VP->R) // Generation anpassen
            if(VP->Mark1)
                VP->Mark1--;
        for(VP=NegIrred.St;VP;VP=VP->R)
            if(VP->Mark1)
                VP->Mark1--;

        // printf("Generationen angepasst\n");

   } // while(1)


    *NewHilb=PosIrred;
    ConcatVLFree(NewHilb,&NeutrIrred);
    FreeVL(&NegIrred);
    FreeVL(&DummyVL);
    // free(SumV);
    // free(SumLV);
}


void DivByHypIt(VList *Hilb, Matrix C0,int Rang0, int EmbDim,Vektor Hyp,int IndHyp,int Size,
                VList *NewHilb, Matrix *NewC0,int *NewRang0)
// Unterteilt Kegel C erzeugt von Hilb entlang Hyperebene Hyp
// Hilb ist Hilbert-Basis der Projektion auf C/max. Unterraum. Letzterer erzeugt von
// C0, seine Dimension ist Rang0.
// IndHyp zeigt an, mit der wievielten Hyperebene jetzt geschnitten wird (erste=0)
// Size ist die Antahl aller Hyperebenen, mit denen geschnitten wird
// Neue Daten fuer Schnitt von C mit positivem Halbraum bez. Hyp werden in "New"
// zurueckgeliefert
// It im Namen zeigt an, dass ueber Supp (implizit) iteriert wird
// In dieser Funktion nur Vorbereitung durch Bestimmung des "neuen" C_0
// und eines Halbraumerzeugers (falls notwendig)
{

    Vektor LF,HspGen;
    Matrix LFDummy; // um aus Vektor einzeilge Matrix zu machen
    Matrix K,Quot;
    boolean C0HSP;

    *NewRang0=Rang0;
    C0HSP=false;

    if(Rang0)
    {
        /* SchreibeV(stdout,Hyp,EmbDim);
        printf("Auf \n");
        SchreibeMatrix(stdout,C0,Rang0,EmbDim); */

        LF=calloc(Rang0,sizeof(Zahl));
        MxV(LF,C0,Rang0,EmbDim,Hyp); // Schraenke Linearform auf bisherigen
                                  // max. Unterraum ein

        /* printf("Kern eingeschraenkt von ");
        SchreibeV(stdout,LF,Rang0); */

        if(!NullVekt(LF,Rang0)) // Falls Linearform auf bisherigem C_0 verschwindet
                                // direkt zur Halbraum-Unterteilung von C
        {
           C0HSP=true;
           MatAlloc(&LFDummy,1,Rang0);    // Mache aus Linearform Matrix
           memcpy(LFDummy[0],LF,Rang0*sizeof(Zahl));

           Kern(LFDummy,1,Rang0,&K,NewRang0,&Quot); // Bestimme neuen max. Unterraum
                                                        //  als Unterraum des alten
           FreeMat(LFDummy,1);

           // printf("Quot \n"); SchreibeMatrix(stdout,Quot,Rang0-*NewRang0,Rang0);

           MatAlloc(NewC0,*NewRang0,EmbDim);
           MxM0(*NewC0,K,*NewRang0,Rang0,C0,EmbDim); // lifte Kern zurueck
           FreeMat(K,*NewRang0);

           // printf("Kern geliftet\n"); SchreibeMatrix(stdout,*NewC0,*NewRang0,EmbDim);

           HspGen=calloc(EmbDim,sizeof(Zahl)); // Lifte Erzeuger eines der Halbraeume
           VxM(HspGen,Quot[0],Rang0,C0,EmbDim);

           // printf("Quot geliftet "); SchreibeV(stdout,HspGen,EmbDim);

           FreeMat(Quot,1);

           C0HSP=true;
        }

        free(LF);
    }

    DivByHypIt0(Hilb, HspGen, EmbDim,Hyp,IndHyp, Size, NewHilb,C0HSP);

    // printf("Raus aus Null\n");

    if(C0HSP)
        free(HspGen);
}

void HilbBasHsp(VList Supp, int EmbDim, VList *Hilb, int *Rang, Matrix *C0Bas, int *Rang0)
// Supp definiert Kegel C durch Stuetzhyperebenen im Eaum der Dimension EmbDim
// Hilb gibt die geliftete Hilberbasis von C/C_0 zurueck (C_0 max Unterraum)
// Die zurueckgegebenen Vektoren sind mit Hoehen (Values) ueber Supp versehen
// Rang ist die Dimension von C
// C0Bas eine Basis von C_0, wird nur allokiert wenn Rang0>0 (dann Rang0 x EmbDim)
// Rang0 Dimension von C_0
{

    VList OldHilb;
    Matrix OldC0Bas, RangTest0;
    int OldRang0, IndHyp, Size,i,k,NPos,NNeg, NHyp,hh,NNeutr,ss;
    VLE *VP, *VP1,*VP2, **SuppOrd;
    Zahl DummyIndex,S;
    float Min;

    MatAlloc(&OldC0Bas,EmbDim,EmbDim);
    MacheEinhMat(OldC0Bas,EmbDim);
    OldRang0=EmbDim;
    StartVL(Hilb);
    StartVL(&OldHilb);

    Registriere(Supp, &SuppOrd,&Size);

    k=0;
    for(VP=Supp.St;VP;VP=VP->R)
    {
        S=0;
        for(i=0;i<EmbDim;i++)
            S+=Zabs(VP->V[i]);
        VP->Mark=S*Size+k;
        k++;
    }
    VEKTLAENGE=0;
    qsort(SuppOrd,Size,sizeof(VLE*),VLEVerglLV);
    /* printf("Geornet \n");
    for(ss=0;ss<Size;ss++)
    {
        Schreibe(stdout,SuppOrd[ss]->Mark);printf(" -- ");
        SchreibeV(stdout,SuppOrd[ss]->V,EmbDim);
    }
    printf("*******************\n"); */


    IndHyp=-1;
    for(ss=0;ss<Size;ss++)
    {
        IndHyp++;
        VP=SuppOrd[ss];

        DivByHypIt(&OldHilb, OldC0Bas,OldRang0, EmbDim,VP->V,IndHyp,Size,
                Hilb, C0Bas,Rang0);

        /* printf("Fertig mit Hyperebene %d\n", IndHyp +1);
        printf("================\n"); */
        // if(IndHyp>2) exit(0);
        /* SchreibeVL(stdout,*Hilb,EmbDim);
        printf("Rang0 = %d\n",*Rang0);
        if(Rang0)
                SchreibeMatrix(stdout,*C0Bas,*Rang0,EmbDim);
        printf("=============\n"); */

        // if(hh==NHyp)

        // if(!VP->R)
        if(ss==Size-1)
            break; // fertig
        OldHilb=*Hilb;
        OldRang0=*Rang0;
        KopiereMat(OldC0Bas,OldRang0,EmbDim,*C0Bas);
        /* printf("NeuerKern\n");
        SchreibeMatrix(stdout,OldC0Bas,*Rang0,EmbDim);*/
    }

    FreeMat(OldC0Bas,EmbDim);
    free(SuppOrd);

    if(!NinVL(*Hilb))
    {
        *Rang=*Rang0;
        return;
    }

    if(!*Rang0)
    {
        RangIndexVL(*Hilb,EmbDim,Rang,&DummyIndex);
        return;
    }

    k=NinVL(*Hilb);
    MatAlloc(&RangTest0,k+*Rang0,EmbDim);
    VL2Mat(RangTest0,*Hilb,EmbDim);
    for(i=0;i<*Rang0;i++)
        memcpy(RangTest0[i+k],*C0Bas[i],EmbDim*sizeof(Zahl));

    RankIndexMat0(RangTest0,k+*Rang0,EmbDim,Rang,&DummyIndex);
    FreeMat(RangTest0,k+*Rang0);

}

void LatPointsHsp(VList Supp, int EmbDim, VList *LattPoints)
// Supp definiert Kegel C durch Stuetzhyperebenen im Eaum der Dimension EmbDim
// LattPoints gibt die geliftete Hilberbasis von C/C_0 zurueck (C_0 max Unterraum)
// Die zurueckgegebenen Vektoren sind mit Hoehen (Values) ueber Supp versehen
// Rang ist die Dimension von C
// C0Bas eine Basis von C_0, wird nur allokiert wenn Rang0>0 (dann Rang0 x EmbDim)
// Rang0 Dimension von C_0
{

    VList OldHilb,Hilb;
    Matrix OldC0Bas, C0Bas, RangTest0;
    int OldRang0, Rang0, IndHyp, Size,i,k;
    VLE *VP,*VP1,*VPNext;
    Zahl DummyIndex;

    MatAlloc(&OldC0Bas,EmbDim,EmbDim);
    MacheEinhMat(OldC0Bas,EmbDim);
    OldRang0=EmbDim;
    StartVL(&OldHilb);
    StartVL(&Hilb);
    IndHyp=-1;
    Size=NinVL(Supp);

    for(VP=Supp.St;VP;VP=VP->R)
    {
        IndHyp++; printf("Hyp %d\n",IndHyp);
        VPNext=VP->R;  // fuer Pruefung ob im Kegel

        // printf("Drin %d\n", IndHyp+1);

        DivByHypIt(&OldHilb, OldC0Bas,OldRang0, EmbDim,VP->V,IndHyp,Size,
                &Hilb, &C0Bas,&Rang0);

        /* printf("Raus %d\n", IndHyp +1);
        SchreibeVL(stdout,Hilb,EmbDim);
        printf("Rang0 = %d\n",Rang0);
        if(Rang0)
                SchreibeMatrix(stdout,C0Bas,Rang0,EmbDim);
        printf("=============\n"); */

        FreeVL(&OldHilb);
        if(!VP->R)
            break; // fertig
        OldHilb=Hilb;
        OldRang0=Rang0;
        KopiereMat(OldC0Bas,OldRang0,EmbDim,C0Bas);
        /* printf("NeuerKern\n");
        SchreibeMatrix(stdout,OldC0Bas,Rang0,EmbDim); */

        if(OldRang0)  // noch kein positiver Kegel
            continue;

        for(VP1=Hilb.St;VP1;VP1=VP1->R) // Schon im oberen Halbraum?
            if(VP1->V[EmbDim-1]<1)
            {
                // SchreibeV(stdout,VP1->V,EmbDim);
                break;
            }
        if(!VP1)                       // ja
            break;
    }

    // printf("Raus mit %d\n",IndHyp+1);


    FreeMat(OldC0Bas,EmbDim);

    StartVL(LattPoints);
    for(VP=Hilb.St;VP;VP=VP->R)
    {
        if(VP->V[EmbDim-1]!=1) // Auf Hoehe 1??
            continue;
        for(VP1=VPNext;VP1;VP1=VP1->R)
            if(SkalProd(VP1->V,VP->V,EmbDim)<0) // Im Kegel ??
                break;
        if(VP1)                // Nein
            continue;
        CopyVLE0(LattPoints,VP,EmbDim);
    }

    FreeVL(&Hilb);
}

void MakeFanMatrix(Matrix *Total, int *NOneSkel, int *NTri, int *NCols, VList OneSkel,
                      SimplexList Tri, int Dim)
{
    int i,j,k;
    VLE *VP,*VP1;
    Simplex *Simpl;
    boolean New, Freimachen;

    *NOneSkel=NinVL(OneSkel);
    *NTri=NinSimplL(Tri);
    *NCols=(Dim+1)*(*NTri);

    MatAlloc(Total,*NOneSkel,*NCols);

    j=-1;
    k=*NTri-1;

    Freimachen=!Tri.St->Corn.St;

    if(Freimachen)
        SetCornSupp(OneSkel,Tri,Dim);

    for(Simpl=Tri.St;Simpl;Simpl=Simpl->R)
    {
        j++;
        for(i=0;i<*NOneSkel;i++)
            (*Total)[i][j]=0;
        for(i=0;i<Dim;i++)
        {

            // printf("%d %d %d\n",j,i,Simpl->Vert[i]);
            (*Total)[Simpl->Vert[i]][j]=-1;
        }

        for(VP=Simpl->Supp.St;VP;VP=VP->R)
        {
            k=k+1;
            i=-1;

            for(VP1=OneSkel.St;VP1;VP1=VP1->R)
            {
                i++;
                (*Total)[i][k]=SkalProd(VP->V,VP1->V,Dim);
            }
        }
    }
    if(Freimachen)
        StripCornSupp(Tri);
}

int SucheFanInFile(Matrix mat1, VList OneSkel, SimplexList Tri, int NOneSkel,
             int NTri,int NCols, int BasicDim,char *filename, int schreiben)
{
    Matrix mat1_zeis, mat1_spas, mat1_zeiz, mat1_spaz;
    Matrix mat2, mat2_zeis, mat2_spas, mat2_zeiz, mat2_spaz;
    int ***passpal;

    VList TestOneSkel;
    SimplexList TestTriang;
    Matrix Total;

    int anzerz,dim,Dummy1,Dummy2,Dummy3;
    int i,j, anzerz2, dim2, dummy, Rueckgabe;
    int iso;
    char zeile[1000], filenamelock[100];

    FILE *ein, *aus,*lock;

    anzerz=NOneSkel; // Anpassung
    dim=NCols;

    if (!(ein=fopen(filename,"r")))
    {                                           // Datenbasis leer
        if(!schreiben) return(0);

        aus=fopen(filename,"w");
        fprintf(aus,"Fan\n");
        SchreibeTriang(aus,OneSkel,BasicDim,Tri);
        fclose(aus);
        return(-1);
    }


    MatAlloc(&mat2,anzerz+1,dim+1);
    MatAlloc(&mat1_zeis,anzerz+1,dim+1); // fuer sortierte Zeilen
    MatAlloc(&mat2_zeis,anzerz+1,dim+1);
    MatAlloc(&mat1_spas,dim+1,anzerz+1);  // fuer sortierte Spalten
    MatAlloc(&mat2_spas,dim+1,anzerz+1);
    MatAlloc(&mat1_zeiz,anzerz+1,anzerz+1); // fuer Indizes der nach Sortierung
    MatAlloc(&mat2_zeiz,anzerz+1,anzerz+1); // gleichen Zeilen bzw.
    MatAlloc(&mat1_spaz,dim+1,dim+1);   // Spalten. Dabei mat1 fuer Vergleich
                                        // mat1 <--> mat1, mat2 fuer mat1 <--> mat2
    MatAlloc(&mat2_spaz,dim+1, dim+1);

    passpal = (int***) calloc(anzerz+1, sizeof(int**)); // Indizes mit gleichen Spalten nach
                                                        // partieller Zuordnung
    for (i=0; i<=anzerz; i++)
    {
        passpal[i] = (int **) calloc(dim+1, sizeof (int*));
        for(j=0; j<= dim; j++) passpal[i][j] = (int*) calloc(dim+1, sizeof(int));
    }

    kopieremat_tr1(mat1,mat1_zeis,mat1_spas, anzerz, dim);
    for(i=0;i<=anzerz;i++) qsort(mat1_zeis[i],dim+1,sizeof(Zahl),Zahlvergl);
    for(j=0;j<=dim;j++) qsort(mat1_spas[j],anzerz+1,sizeof(Zahl),Zahlvergl);
    findepartner(mat1_zeis,mat1_zeis,anzerz,dim,mat1_zeiz); // Haeufigkeiten der sortierten
    findepartner(mat1_spas,mat1_spas,dim,anzerz,mat1_spaz); // Zeilen und Spalten finden

    iso=0;

    while(fscanf(ein,"%s", zeile)>=1)
    {
        if(strncmp(zeile, "Fan",3)) continue; // Noch nicht Anfang gefunden
        iso++;

        fscanf(ein,"%d", &anzerz2);
        fscanf(ein,"%d", &dim2);
        if(NOneSkel!=anzerz2) continue; // Format passt nicht

        LiesTriang(ein,&TestOneSkel, NOneSkel, BasicDim, &TestTriang); // Jetzt passt es

        if(NinSimplL(TestTriang)!=NTri) // Format passt nicht
        {
            FreeVL(&TestOneSkel);
            FreeSimplexL(&TestTriang);
            continue;
        }

        MakeFanMatrix(&Total, &Dummy1,&Dummy2,&Dummy3,TestOneSkel,
                      TestTriang,BasicDim);
        for(i=0;i<NOneSkel;i++)
            for(j=0;j<NCols;j++)
                mat2[i+1][j+1]=Total[i][j];
        FreeMat(Total,NOneSkel);
        FreeVL(&TestOneSkel);
        FreeSimplexL(&TestTriang);



        /* printf("Mat 2\n");
        schreibematrix2(mat1,anzerz,dim,stdout,0); */

        if(matpermvgl_vorb_alloc(mat1, mat1_zeis, mat1_spas, mat1_zeiz, mat1_spaz,
                         mat2, mat2_zeis, mat2_spas,  mat2_zeiz, mat2_spaz,
                           anzerz, dim, passpal))
        {
            fclose(ein);
            Rueckgabe=iso;
            goto Ende;
        }

    }

    fclose(ein);

    if(schreiben)
    {
        strcpy(filenamelock,filename);
        strcat(filenamelock,".lock");
        while(FileExists(filenamelock,""))
        {
            printf("File %s locked\n",filename);
            for(i=0;i<1000000;i++)
            {
                dummy=i;
            }
        }
        lock=fopen(filenamelock,"w");
        aus=fopen(filename,"a");
        fprintf(aus,"Fan\n");
        SchreibeTriang(aus,OneSkel,BasicDim,Tri);
        fclose(aus);
        fclose(lock);
        remove(filenamelock);
        Rueckgabe=-(iso+1);
        goto Ende;
    }

    Rueckgabe=0;

    Ende:

    FreeMat(mat2,anzerz+1);
    FreeMat(mat1_zeis,anzerz+1);
    FreeMat(mat2_zeis,anzerz+1);
    FreeMat(mat1_spas,dim+1);
    FreeMat(mat2_spas,dim+1);
    FreeMat(mat1_zeiz,anzerz+1);
    FreeMat(mat2_zeiz,anzerz+1);
    FreeMat(mat1_spaz,dim+1);
    FreeMat(mat2_spaz,dim+1);
    for (i=0; i<=anzerz; i++)
    {
        for(j=0; j<= dim; j++) free(passpal[i][j]);
        free(passpal[i]);
    }
    free(passpal);

    return Rueckgabe;
}


int FanInFile(Matrix Total,VList OneSkel, SimplexList Tri,int NOneSkel,int NTri,int NCols,int BasicDim,
                   char *database,boolean schreiben, boolean verbose)
// Gibt <0 zurueck, falls nicht vorhanden, sonst Nummer der schon vhd Matrix oder 0,
// falls neu, aber nicht geschrieben
// Analog matinfile
{
    int i,j, iso;
    Matrix mat1;
    char filename2[100];
    FILE *ein, *aus;

    sprintf(filename2,"%s_%d.base",database,BasicDim);

    MatAlloc(&mat1, NOneSkel+1,NCols+1+1);
    for(i=0;i<NOneSkel;i++)
        for(j=0;j<NCols;j++)
            mat1[i+1][j+1]=Total[i][j];

    iso=SucheFanInFile(mat1,OneSkel,Tri,NOneSkel,NTri,NCols,BasicDim,filename2,schreiben);

    FreeMat(mat1,NOneSkel+1);

    if(iso>0)
    {
        if(verbose)
            printf("Already in  %s Typ %d\n", filename2, iso);
        return(iso);
    }
    if(iso<0)
    {
        printf("New in  %s Typ %d\n", filename2, -iso);
        return(iso);
    }

    printf("Not in  %s\n", filename2);
    return 0;

}



boolean FindInBaseFan(VList OneSkel, SimplexList Tri, int Dim, char *database, boolean write, boolean verbose)
{
    Matrix Total;
    int NOneSkel,NTri,NCols;
    boolean New;

    New=false;

    MakeFanMatrix(&Total,&NOneSkel,&NTri,&NCols,OneSkel,Tri,Dim);
    if(FanInFile(Total,OneSkel,Tri,NOneSkel,NTri,NCols,Dim,database,write,verbose)>0)
        New=true;
    FreeMat(Total,NOneSkel);

    return New;
}
