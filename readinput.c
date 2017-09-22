#include "head.h"

// read the input file
int readinput(FILE *fp, struct MOLE *M1, struct MOLE *M2, struct MOLE *Ms, struct PARS *pars, int *nodes, int *snaps, int *mode, int *nsolvent, int *nsolute)
{
	int i,j,num1,num2;
	char temp[FLN];
	char *pp,*tok;

	// initialize
	M1->idi = M1->idf = M1->n = 0;
	M2->idi = M2->idf = M2->n = 0;
	Ms->idi = Ms->idf = Ms->n = 0;
	for(i=0;i<MAXNATOM;i++)
		for(j=0;j<4;j++)
			M1->ATOMS[i][j]=M2->ATOMS[i][j]=Ms->ATOMS[i][j]='Z';
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			pars->lifetime[i][j]=pars->persist[i][j]=pars->tolerance[i][j]=0;
	
	printf("\nThe input file is:\n");
	// read the input file
	rewind(fp);
	while(fgets(temp,sizeof(temp),fp)!=NULL)
	{
		if((pp=strstr(temp,"[END]"))!=NULL)
			break;

		if((pp=strstr(temp,"[NUMBER OF NODES]"))!=NULL)
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			*nodes=atoi(tok);
		}

		if((pp=strstr(temp,"[NUMBER OF SNAPS]"))!=NULL)
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			*snaps=atoi(tok);
		}

		if((pp=strstr(temp,"[MODE]"))!=NULL)
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			*mode=atoi(tok);
		}

		if((pp=strstr(temp,"[NUMBER OF SOLVENT]"))!=NULL)
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			*nsolvent=atoi(tok);
			num1=*nsolvent * 2;
		}

		if((pp=strstr(temp,"[NUMBER OF SOLUTE]"))!=NULL)
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			*nsolute=atoi(tok);
			num2=*nsolute * 2;
		}

		if((pp=strstr(temp,"[IDS OF SOLVENT1]"))!=NULL && num1>0) //read the ids for solvent1
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			M1->idi=atoi(tok);
			tok=strtok(NULL," \n");
			M1->idf=atoi(tok);
			num1=num1-1;
		}

		if((pp=strstr(temp,"[ATOMS OF SOLVENT1]"))!=NULL && num1>0) //read the atoms for solvent1
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			M1->n=atoi(tok);
			for(i=1;i<=M1->n;i++)   // remember atoms from 1 to MAXNATOM
			{
				tok=strtok(NULL," \n");
				strcpy(M1->ATOMS[i],tok);
			}
			num1=num1-1;
		}

		if((pp=strstr(temp,"[IDS OF SOLVENT2]"))!=NULL && num1>0) //read the ids for solvent2
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			M2->idi=atoi(tok);
			tok=strtok(NULL," \n");
			M2->idf=atoi(tok);
			num1=num1-1;
		}

		if((pp=strstr(temp,"[ATOMS OF SOLVENT2]"))!=NULL && num1>0) //read the atoms for solvent2
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			M2->n=atoi(tok);
			for(i=1;i<=M2->n;i++)
			{
				tok=strtok(NULL," \n");
				strcpy(M2->ATOMS[i],tok);
			}
			num1=num1-1;
		}

		if((pp=strstr(temp,"[IDS OF SOLUTE]"))!=NULL && num2>0) //read the ids for solute
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			Ms->idi=atoi(tok);
			tok=strtok(NULL," \n");
			Ms->idf=atoi(tok);
			num2=num2-1;
		}

		if((pp=strstr(temp,"[ATOMS OF SOLUTE]"))!=NULL && num2>0) //read the atoms for solute
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			Ms->n=atoi(tok);
			for(i=1;i<=M2->n;i++)
			{
				tok=strtok(NULL," \n");
				strcpy(Ms->ATOMS[i],tok);
			}
			num2=num2-1;
		}

		if((pp=strstr(temp,"[LIFETIME OF SOLVENT1 SOLVENT1]"))!=NULL) //read parameters for solvent1
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			pars->lifetime[1][1]=atoi(tok);
			if(pars->lifetime[1][1]==1)
			{
				tok=strtok(NULL," \n");
				pars->persist[1][1]=atoi(tok);
				tok=strtok(NULL," \n");
				pars->tolerance[1][1]=atoi(tok);
			}
		}

		if((pp=strstr(temp,"[LIFETIME OF SOLVENT2 SOLVENT2]"))!=NULL) //read parameters for solvent1
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			pars->lifetime[2][2]=atoi(tok);
			if(pars->lifetime[2][2]==1)
			{
				tok=strtok(NULL," \n");
				pars->persist[2][2]=atoi(tok);
				tok=strtok(NULL," \n");
				pars->tolerance[2][2]=atoi(tok);
			}
		}

		if((pp=strstr(temp,"[LIFETIME OF SOLVENT1 SOLVENT2]"))!=NULL) //read parameters for solvent1
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			pars->lifetime[1][2]=atoi(tok);
			if(pars->lifetime[1][2]==1)
			{
				tok=strtok(NULL," \n");
				pars->persist[1][2]=atoi(tok);
				tok=strtok(NULL," \n");
				pars->tolerance[1][2]=atoi(tok);
			}
		}

		if((pp=strstr(temp,"[LIFETIME OF SOLVENT2 SOLVENT1]"))!=NULL) //read parameters for solvent1
		{
			printf("%s",pp);
			tok=strtok(pp," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			tok=strtok(NULL," \n");
			pars->lifetime[2][1]=atoi(tok);
			if(pars->lifetime[2][1]==1)
			{
				tok=strtok(NULL," \n");
				pars->persist[2][1]=atoi(tok);
				tok=strtok(NULL," \n");
				pars->tolerance[2][1]=atoi(tok);
			}
		}

	}


	return(0);
}

