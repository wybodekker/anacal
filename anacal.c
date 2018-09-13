#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <nr.h>
#include <nrutil.h>
#define ERROR 2
#define EOL 1
#define OK 0
#define MAX_ELEMENTS 20       /* max no of elements in a bruto formula */
#define MAX_FORMULAE 26       /* max no of components / formulas */

FILE *in,*out;
typedef char str[4];
typedef struct {
   int nel;			/* number of elements in formula */
   float mw;			/* molecular weight */
   int elnr[MAX_ELEMENTS];	/* element nr (C=6, H=1 etc) */
   str elsym[MAX_ELEMENTS];	/* element symbol */
   int elcnt[MAX_ELEMENTS];	/* element counts */
   float elfrac[MAX_ELEMENTS];	/* element fractions (analyzed elements only!) */
   char title[80];
}
formula;
formula form[MAX_FORMULAE];
formula comp;		/* contains all data for the current composition */
str symbol;
float d;
int nf,kort=0,weight=0;
int indx[MAX_FORMULAE+1];
char temp[80];

#include <elements.h>

int comp_nr=0; /* counts compositions (just for printing references in the tables) */

void helpsrt() {
   system("anacal -H|sed -ns '/^= Syno/,/^= Desc/p'|sed '/^=/d'");
   exit(0);
}

void helpall() {
   fprintf(stdout,"\n= anacal - calculation of chemical composition from elemental analysis results\n\
\n\
= Synopsis\n\
anacal [options] file[.ea] [max_dev]	\n\
\n\
== Options\n\
-s	use the short format\n\
-w	in the short format, report weight instead of mole fractions\n\
-h	prints this help\n\
-H	prints full documentation\n\
\n\
= Description\n\
Given, for a chemical composition:\n\
- the elemental percentages for (some of the) elements and\n\
- a set of bruto formulae, one for each possible component \n\
anacal calculates compositions with one or more of the formula's\n\
resulting in a good (that is: better than a certain standard deviation) fit\n\
between calculated and experimental element percentages. \n\
\n\
Input is taken from |file.ea|, output is written to |stdout|.\n\
The extension .ea is automatically appended to the input file.\n\
\n\
Run without parameters for help.\n\
\n\
= Input\n\
The input file should have the following structure (numbers are line\n\
numbers):\n\
1	a title line\n\
2	elementsymbol followed by percentage for each element determined\n\
3+	bruto formula's of the possible components, one per line;\n\
	comments may be added with the #-character;\n\
	remaining characters up to end-of-line are skipped.\n\
\n\
Sample input file (case sensitive):\n\
\n\
   6718 DEBO 145\n\
   C54.53 H6.31 N8.10 S14.48\n\
   C15 H20 N2 O3 S\n\
   C11 H16 N2 O S\n\
   C4 H12 N2 S2\n\
   C2 H4 O S\n\
   S	#sulphur\n\
   H Cl	#hydrochloric acid\n\
   H2 O	#water\n\
\n\
= Author\n\
[Wybo Dekker](wybo@dekkerdocumenten.nl)\n\
\n\
= Copyright\n\
Released under the [GNU General Public License](www.gnu.org/copyleft/gpl.html)\n");
   exit(0);
}

int detsol(float **a, int n, float *b) {
   float d,mx;
   int i,j;
   /* check for singularity */
   for (i=1;i<=n;i++) {
      mx=0.0;
      for (j=1;j<=n;j++)
         if (fabs(a[i][j]) > 0.0)
            mx=1;
      if (mx == 0.0)
         return ERROR;
   }
   ludcmp(a,n,indx,&d);
   lubksb(a,n,indx,b);
   return OK;
}

int leeschar() {
   char s;
nextchar:
   s=fgetc(in);
   if (isalnum(s))	return s;
   if (s=='.')		return s;
   if (s=='\n')		return s;
   if (s==EOF)		return s;
   if (s=='#')		fscanf(in," %[^\n]",form[nf].title);
   goto nextchar;
}
int get_elnr_count(int *elnr,float *count,int *last) {
   char s,sym[3];
   int flag;

   /* read element symbol: capital, eventually followed by lower case letter */
   if (!isupper(sym[0]=leeschar()))
      return sym[0]==EOF ? EOF : ERROR;
   if (islower(sym[1]=leeschar()))
      sym[2]=0;
   else {
      ungetc(sym[1],in);
      sym[1]=0;
   }

   *elnr=0;
   flag=ERROR;
   while(elements[++(*elnr)].symbol[0])	/* last symbol in elements is empty */
      if (!strcmp(sym,elements[*elnr].symbol)) {
         flag=OK;
         break;
      }
   if (flag==ERROR) {
      fprintf(stderr,"unknown element: %s\n",sym);
      return ERROR;
   }

   /* read number */
   s=leeschar();
   ungetc(s,in);
   if (s=='.' || isdigit(s))
      fscanf(in,"%g",count);
   else
      *count=1;
   if (!(*last=(s=leeschar())=='\n'))
      ungetc(s,in);
   return 0;
}

int main(int argc, char *argv[]) {
   int i,j,k,l,m,c,max_used,wf,nf_used,eol;
   int form_nrs_used[MAX_ELEMENTS+1];
   int argn=1;
   float elfrac_calcd[MAX_ELEMENTS];	/* parallells comp.elfrac (=fractions found) */
   float frac_of_form[MAX_FORMULAE];
   float molfrac_of_form[MAX_FORMULAE];
   float **bs,p,q,som,limit=1;
   char title[200],label[20],infile[20],outfile[20],command[20]="mk ";
   char *optptr,*argptr=argv[1]; //pointer to first argument
   char coldefs[20];
   int wid_of_formula_col;
   char *hyphens="---------------------------------------------------------------------------------";
   if (argc<2)
      helpsrt();
   while (argptr[0]=='-') 	// if first argument is an option, like -ws
   { optptr=argptr+1;		// point to the 'w', then the 's'
      while(*optptr) {
         switch(*optptr) {
            case 's': kort=1; break;
            case 'w': weight=1; break;
            case 'V': printf("1.03\n");exit(0);
            case 'h': helpsrt(); break;
            case 'H': helpall(); break;
         }
         optptr++;
      }
      argc--;
      argptr=argv[++argn];
      ;
   }

   strcat(command,argptr);
   strcpy(infile,argptr);
   strcpy(label,argptr);
   if (!strstr(infile,".ea")) {
      strcat(infile,".ea");
   }
   if (!(in=fopen(infile,"rt"))) {
      fprintf(stderr,"%s: inputfile %s not found\n",argv[0],infile);
      exit(1);
   }
   strcpy(outfile,argptr);
   strcat(outfile,".tex");

   out=stdout;
   if (argc>3)
      helpsrt(); // there must be a one and only one filename plus perhaps a limit
   if (argc==3)
      limit=atof(argv[argn+1]); //sscanf(argptr,"%f",&limit); else limit=1;
   limit/=100;

   /* copy title line from input to output */
   fgets(title,80,in);
   title[strlen(title)-1]=0;

   /* read bruto formula´s;  the first is the composition found */

   nf=-1;        /* first formula read is for comp */
nextformula:
   if (nf>=MAX_FORMULAE) {
      fprintf(stderr,"Too many formulae\n");
      exit(1);
   }
   k=0;
nextelement:
   if (k>=MAX_ELEMENTS) {
      fprintf(stderr,"Too many elements\n");
      exit(1);
   }
   i=get_elnr_count(&j,&d,&eol);
   if (i==ERROR) {
      fprintf(stderr,"Error reading input file line %d\n",nf+3);
      exit(1);
   }
   if (i!=EOF) {
      if (nf<0) {
         comp.elnr[k]=j;
         strcpy(comp.elsym[k],elements[j].symbol);
         comp.elfrac[k]=d/100;
      } else {
         form[nf].elnr[k]=j;
         strcpy(form[nf].elsym[k],elements[j].symbol);
         form[nf].elcnt[k]=d+.5;
         d*=elements[j].chemweight;  /* this elements' contribution to mw */
         for (m=0;m<comp.nel;m++)
            if (comp.elnr[m]==j)     /* if element is an analyzed one */
            { form[nf].elfrac[m]=d;  /* update calcd value */
               break;
            }
         form[nf].mw+=d;             /* and update molecular weight */
      }
      k++;
      if (eol) {
         if (nf<0)
            comp.nel=k;
         else {
            form[nf].nel=k;
            for (i=0;i<comp.nel;i++)
               form[nf].elfrac[i]/=form[nf].mw;
         }
         nf++;
         goto nextformula;
      } else
         goto nextelement;
   }
   if (kort) {
      sprintf(coldefs,"r*{%d}{r}c*{%d}{r}",comp.nel+1,nf);
   } else {
      strcpy(coldefs,"rlrrrrrrrl");
   }
   wid_of_formula_col=-2*comp.nel-2;

   fprintf(out,
           "ANACAL - error limit %g%% - %s\n\n"
           "Candidate components:\n"
           "nr molwgt %*s name\n"
           ,limit*100
           ,title
           ,wid_of_formula_col
           ,"formula"
          );

   for (i=0;i<nf;i++) {
      fprintf(out,"%2d %6.2f",i+1,form[i].mw);
      *temp='\0';
      for (j=0;j<form[i].nel;j++) { // collect formula in temp
         strcat(temp,form[i].elsym[j]);
         if (form[i].elcnt[j]>1)
            sprintf(temp+strlen(temp),"%d",form[i].elcnt[j]);
      }
      fprintf(out," %*s %s\n",wid_of_formula_col,temp,form[i].title);
   }
   fprintf(out,"\n");

   if (kort) {
         char *h1=strndup(hyphens,(3+comp.nel*6-19)/2);
         char *h2=strndup(hyphens,(nf*5-6)/2);
         fprintf(out,"     %sfound-calcd (%%x100)%s  %s%s‰%s",h1,h1,h2,weight ? "wght" : "mole",h2);
         fprintf(out,"\nnr   std");
         for (i=0;i<comp.nel;i++)
            fprintf(out,"%6s",comp.elsym[i]);
         for (i=0;i<nf;i++)
            fprintf(out," %4d",i+1);
         fprintf(out,"\n");
   } else {
      fprintf(out,
              "                             --percent  by--\n"
              "       found  calcd st.dev   weight     mole  nr   molwgt   formula\n"
             );
   }

   // the following generates all possible permutations of the formulae
   // for example, if the maximum number of formulae is 4 then the following
   // combinations are generated (and placed in form_nrs_used and evaluated): 
   // 0
   // 1
   // 2
   // 3
   // 01
   // 02
   // 03
   // 12
   // 13
   // 23
   // 012
   // 013
   // 023
   // 123
   // 0123

   max_used=nf>comp.nel ? comp.nel+1 : nf;
   nf_used=1;
   wf=0;

next_combination:
   i=nf_used;
   if (wf==nf) {
L6:
      if (i==1) {
         if (nf_used==max_used) {
            exit(0);
         }
         nf_used++;
         wf=0;
      } else {
         j=wf;
         i--;
         wf=form_nrs_used[i-1]+1;
         if (wf+1==j)
            goto L6;
      }
   }
   for (;i<=nf_used;i++)
      form_nrs_used[i-1]=wf++;

   // now one combination is in form_nrs_used and needs evaluation

   if (nf_used>1) {
      bs=matrix(1,nf_used-1,1,nf_used-1);
      // detsol needs fortran dimensioned arrays, so we declare bs with indices atarting at 1
      // and fill it with indices incremented with 1
      // and we transfer the frac_of_form array as frac_of_form-1
      for (i=0;i<nf_used-1;i++) {
         j=form_nrs_used[i];
         for (k=0;k<nf_used-1;k++)
            bs[i+1][k+1]=0;
         frac_of_form[i]=0;
         for (c=0;c<comp.nel;c++) {
            q=form[wf-1].elfrac[c];
            p=form[j].elfrac[c]-q;
            for (k=0;k<nf_used-1;k++)
               bs[i+1][k+1]+=(form[form_nrs_used[k]].elfrac[c]-q)*p;
            frac_of_form[i]+=(comp.elfrac[c]-q)*p;
         }
      }
      if (detsol(bs,nf_used-1,frac_of_form-1)==ERROR)
         goto next_combination;
      free_matrix(bs,1,nf_used-1,1,nf_used-1);
   }
   som=0;
   q=0;
   for (i=0;i<nf_used;i++) {
      if (i==nf_used-1)
         frac_of_form[i]=p=1-q;
      else
         p=frac_of_form[i];
      if (p<0)
         goto next_combination;     /* skip negative concentrations */
      q+=p;
      p/=form[form_nrs_used[i]].mw;
      molfrac_of_form[i]=p;
      som+=p;
   }
   for (i=0;i<nf_used;i++)
      molfrac_of_form[i]/=som;
   for (q=c=0;c<comp.nel;c++) {
      for (p=i=0;i<nf_used;i++)
         p+=form[form_nrs_used[i]].elfrac[c]*frac_of_form[i];
      q+=(comp.elfrac[c]-p)*(comp.elfrac[c]-p);
      elfrac_calcd[c]=p;
   }
   q=sqrt(q/comp.nel);
   if (q>limit)
      goto next_combination; /* print only if better than limit */
   m=(nf_used>comp.nel)? nf_used : comp.nel;

   if (kort) {
      fprintf(out,"%2d %5.0f",++comp_nr,q*10000);
      for (c=0;c<comp.nel;c++)
         fprintf(out," %5.0f",comp.elfrac[c]*10000-elfrac_calcd[c]*10000);
      for (c=k=0;k<nf;k++) {
         if (k==form_nrs_used[c]) {
            fprintf(out," %4.0f",(weight ? frac_of_form[c] : molfrac_of_form[c])*1000);
            c++;
         } else {
            fprintf(out,"     ");
         }
      }
      fprintf(out,"\n");
   } else {
      fprintf(out,"\n");
      for (c=0;c<m;c++) {
         if (c<comp.nel) {
            if (c==0)
               fprintf(out,"%2d ",++comp_nr); // composition number
            else
               fprintf(out,"   ");
            fprintf(out,"%-2s%7.2f%7.2f",
                    comp.elsym[c],comp.elfrac[c]*100,elfrac_calcd[c]*100);
            	    // symbol        found              calcd */
            if (c==0)
               fprintf(out,"%7.3f ",q*100); // misfit
            else
               fprintf(out,"        ");
         } else
            fprintf(out,"                           ");

         if (c<nf_used) {
            k=form_nrs_used[c];
            fprintf(out,"%8.2f %8.2f %3d  %7.2f   ",
                    frac_of_form[c]*100,
                    molfrac_of_form[c]*100,
                    k+1,
                    form[k].mw);
            /* pct by weight       pct by mole                nr  molwght */
            for (l=0;l<form[k].nel;l++) {
               i=form[k].elcnt[l];
               if (i>1)
                  fprintf(out,"%s%d ",form[k].elsym[l],i);
               /* bruto formula */
               else
                  fprintf(out,"%s ",form[k].elsym[l]);
            }
         }
         fprintf(out,"\n");
      }
   }
   goto next_combination;
   exit(0);
}
