
/* 
 * get the water-water H-bond lifetime
 * Tiecheng Zhou 
 * tiecheng.zhou@email.wsu.edu
 * 2015.0922
 *
 * needs all the graphgeod files, with the format as '12 24 0 0 1 1 3 1.53 120.3', two molecule ids, 3 pbc label, 2 atom ids, HB distance and angle
 *
 * this is the head file
 */


#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define PERCENT 0.5 // default parameter in defining transient bond 
#define Maxneighb 150 // maximum possiblenumber of neighbor for the a molecule
#define FLN 10000 // maximum string length for each line
#define MAXNATOM 20 // maximum number of atoms in a molecule
#define Ncarbon 3

#define Xsize 36.0956137      // 2017.04.17, for propane box, the size of the simulation box, for calculating water:solute distance
#define Ysize 36.0956137
#define Zsize 36.0956137

struct MOLE  // define the struct of a molecule
{
	int idi;
	int idf; // initial and final ids for the molecule
	int n; // number of atoms in the molecule
	char ATOMS[MAXNATOM+1][4];
};

struct PARS // define the struct of parameters for lifetime
{
	int lifetime[3][3];
	int persist[3][3];
	int tolerance[3][3];
};


//read the H-bonding history at snapshot(snap) for target water(watt), and remember the neighbor ids, atoms from target and neighbor
int loadhistory(int **nei, int **att, int **atn, int snap, int watt);


// calculate the total degree, the donating degree and the accepting degree at a certain snapshot(snap)
int caldegrees(int d[MAXNATOM+1], int snap, int dor, int acc, struct MOLE *mt, struct MOLE *pt, int **nei, int **att, int **atn);


//check if the molecule belongs to solvent1 or solvent2
int chkblg(int id, struct MOLE *s1, struct MOLE *s2);


// to find if there is a specific H-bond target(atmt) with water(atmn) at certain snapshot(snap)
int findHbond(int water, int atmt, int atmn, int snap, int **nei, int **att, int **atn, int *nsnaps);


// find the continuous bonded segment starting from snap, get the period ti to tf in array seg[]
int findseg(int seg[2], int water, int atm1, int atm2, int snap, int **nei, int **att, int **atn, int *nsnaps);


// to check if a particular H-bond from snap1 to snap2 is transient bond
int itransbond(int target, int water, int atmt, int atmn, int snap1, int snap2, int dor, int acc,  int **nei, int **att, int **atn, int *nsnaps, struct PARS *ppars);


// find the new H-bonded partner and snapshot with target(atmt)
int findatom(int b[Maxneighb+1], int atm1, int tsnap1, int tsnap2, int **nei, int **att, int **atn, int *nsnaps, int dor, int acc, struct PARS *pars);


// copy a specific H-bond information from refsnap to snap
int myadd(int target, int water, int atom1, int atom2, int snap, int **nei, int **att, int **atn);


// remove a specific H-bond from the original graphgeod
int myremove(int target, int water, int atom1, int atom2, int snap, int **nei, int **att, int **atn);


// read the input files
int readinput(FILE *fp, struct MOLE *M1, struct MOLE *M2, struct MOLE *Ms, struct PARS *pars, int *nodes, int *snaps, int *mode, int *nsolvent, int *nsolute);


// the lifetime correction code in 'corrcode.c'
int corrcode(int dor, int acc, struct MOLE *pm1, struct MOLE *pm2, int *nsolu, struct MOLE *pms, struct PARS *pars, int *nnodes, int *nsnaps, int *ind);


//update the graphgeod files with the new array information
int myupdate(int target, int **nei, int **att, int **atn, int *nsnaps);


//update graphgeod files with one snapshot
int myupdateone(int target, int **nei, int **att, int **atn, int snap);


//calculating the distance from the target molecule to a solute, 2016.03.08
float caldist(float **result, int target, int t, int dor, int acc, struct MOLE *m1, struct MOLE *m2, struct MOLE *mt, int *nnodes);

