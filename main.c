/*
 * this is the main function of the code
 * updated 2015.09.21
 * only check HB breakage in terms of not-changing-partner or changing-partners
 * there's no need inspect wether a HB switch changes via NHB or BHB
 */

#include "head.h"

int main(int argc, char *argv[])
{
// pre-run parameters set

	char finput[FLN];      // intput file
	int mode,nnodes,nsnaps,nsolv,nsolu,*pmode,*pnodes,*psnaps,*pnsolv,*pnsolu;
	FILE *fip;
	struct MOLE mole1,mole2,moles,*pm1,*pm2,*pms;
	struct PARS pars,*ppars;

	pmode=&mode;pnodes=&nnodes;psnaps=&nsnaps;pnsolv=&nsolv;pnsolu=&nsolu;
	pm1=&mole1;pm2=&mole2;pms=&moles,ppars=&pars;

// Start reading the input file
 	if(argc != 2)
	{
		printf("Usuage: %s input-corr\n", argv[0]);
		exit(-1);
	}
	sprintf(finput,"%s",argv[1]);
	if((fip=fopen(finput,"r"))==NULL)
	{
		printf("In main: can not find the inputfile\n");
		exit(-1);
	}

	if(readinput(fip,pm1,pm2,pms,ppars,pnodes,psnaps,pmode,pnsolv,pnsolu)!=0)
	{
		printf("In main: error when reading the inputfile\n");
		exit(-1);
	}

	fclose(fip);

	printf("\n------Output file format:------\n");
	printf("  Degrees-i-j: the degree destribution of solvent i with solvent j\n");
	printf("    atom 1 degree 2 observation 110 , means, there are 110 configurations when atom 1 is interacting with 2 molecules\n");
	printf("  Lifetimes-i-j: the lifetime of the interaction of solvent i with solvent j\n");
	printf("    1 4 3 2 11 50 R N , means, molecule 1 (atom 3) is interacting with molecule 4 (atom 2) from snap 11 to snap 50;\n");
	printf("      this is a Real HB and No switch is involved\n");
	printf("      the two letters: R is short for Real, or T means Transient; and N means No switch, or Y means switch, changing partner\n");
	printf("  Switches-i-j: the switch events of the interaction of solvent i with solvent j\n");
	printf("    R Y 1 146 191 2 27 36 37 50 , means, molecule 1 (atom 2) changes its partner from intial molecule 146 to final molecule 191;\n");
	printf("      the initial interaction lasts from snap 27 to snap 36, and the final interaction is from snap 37 to snap 50\n");
	printf("      the two letters: R is short for Real, or T means Transient; and N means No switch, or Y means switch, changing partner\n");
	printf("    these results are obtained by inspection the dynamic history of solvent i\n\n");

	printf("\n------Start to run %s------\n",argv[0]);

	if(mode==1)
		printf(" Running this code with correcting transient breaks\n\n");
	else if(mode==2)
		printf(" Running this code with correcting transient bonds\n\n");
	else if(mode==3)
		printf(" Running this code with correcting both transient breaks and transient bonds\n\n");
	else if(mode==0)
		printf(" Running this code without correction\n\n");

// End of the pre-run parameters set

	int sym=0;

	if(*pnsolv==1)
	{
		if(sym == 0 && ppars->lifetime[1][1] == 1)
		{
			sym=corrcode(1,1,pm1,pm1,pnsolu,pms,ppars,pnodes,psnaps,pmode); //HB lifetime for solvent1 solvent1
			if(sym != 0)
				printf("In main: error for the HB lifetime for solvent1 solvent1\n");
		}	
	}
	else if(*pnsolv==2)
	{
		if(sym == 0 && ppars->lifetime[1][1] == 1)
		{
			printf("Lifetime calculation of solvent%d solvent%d\n",1,1);
			sym=corrcode(1,1,pm1,pm2,pnsolu,pms,ppars,pnodes,psnaps,pmode); //HB lifetime for solvent1 solvent1
			if(sym != 0)
				printf("In main: error for the HB lifetime for solvent1 solvent1\n");
		}
	
		if(sym == 0 && ppars->lifetime[2][2] == 1)
		{
			printf("Lifetime calculation of solvent%d solvent%d\n",2,2);
			sym=corrcode(2,2,pm1,pm2,pnsolu,pms,ppars,pnodes,psnaps,pmode);
			if(sym != 0)
				printf("In main: error for the HB lifetime for solvent2 solvent2\n");
		}
	
		if(sym == 0 && ppars->lifetime[1][2] == 1)
		{
			printf("Lifetime calculation of solvent%d solvent%d\n",1,2);
			sym=corrcode(1,2,pm1,pm2,pnsolu,pms,ppars,pnodes,psnaps,pmode);
			if(sym != 0)
				printf("In main: error for the HB lifetime for solvent1 solvent2\n");
		}
	
		if(sym == 0 && ppars->lifetime[2][1] == 1)
		{
			printf("Lifetime calculation of solvent%d solvent%d\n",2,1);
			sym=corrcode(2,1,pm1,pm2,pnsolu,pms,ppars,pnodes,psnaps,pmode);
			if(sym != 0)
				printf("In main: error for the HB lifetime for solvent2 solvent1\n");
		}
	}

	return(0);
}



