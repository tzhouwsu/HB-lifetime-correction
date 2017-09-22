/*
 * this file has all the functions that are used in the code
 */


#include "head.h"

int loadhistory(int **nei, int **att, int **atn, int snap, int watt)
{
	FILE *fp;
	char file[100];
	char buffer[100];
	int number,node1,node2,pbcx,pbcy,pbcz,atm1,atm2;
	int dup,k,m;
	sprintf(file, "Input.water%d.xyz.solB%d.xyz.GraphGeod",snap,snap); // the name of the graphgeod file here, two "%d" will be filled with value of integer number "snap"
	if((fp=fopen(file,"r"))==NULL)
	{
		printf("Error: can't read the graphgeod file, snap %d, water %d\n",snap,watt);
		fclose(fp);
		return(-1);
	}

	rewind(fp);
	number=0;
	while(fscanf(fp, "%d %d %d %d %d %d %d", &node1,&node2,&pbcx,&pbcy,&pbcz,&atm1,&atm2)!=EOF) // read two molecule ids
	{
		fgets(buffer,sizeof(buffer),fp);
		
		if(node1==watt)
		{
			dup=0;
			for(k=1;k<=number;k++)         // check duplicate partners
				if(nei[snap][k]==node2 && att[snap][k]==atm1 && atn[snap][k]==atm2)
					dup=1;
			if(dup==0)                     // it's a new neighbor
			{
				number += 1;
				nei[snap][number]=node2; // water id that is H-bonded with solute
				att[snap][number]=atm1;
				atn[snap][number]=atm2;

			}

		}
		else if(node2==watt)
		{
			dup=0;
			for(k=1;k<=number;k++)
				if(nei[snap][k]==node1 && att[snap][k]==atm2 && atn[snap][k]==atm1)   // bug found 2015-12-10, atm1,atm2 changed
					dup=1;
			if(dup==0)
			{
				number += 1;
				nei[snap][number]=node1; // water id that is H-bonded with solute
				att[snap][number]=atm2;
				atn[snap][number]=atm1;

			}
		}

	}
	nei[snap][0]=number; // the number of neighbors is kept
	att[snap][0]=number;
	atn[snap][0]=number;

	fclose(fp);
	return(0);
}

// calculate the total degree, the donating degree and the accepting degree at a certain snapshot(snap)
int caldegrees(int d[MAXNATOM+1], int snap, int dor, int acc, struct MOLE *mt, struct MOLE *pt, int **nei, int **att, int **atn)
{
	int j,atomt,tot;
	int partner,group;

	for(j=0;j<=MAXNATOM;j++)
		d[j]=0;

	// the 0th field of array d[] keeps the total degree, others keep degree with respect to each atom
	for(j=1;j<=nei[snap][0];j++)
	{
		atomt=att[snap][j];
		partner=nei[snap][j];
		group=chkblg(partner,mt,pt); // check whether the partner belongs to solvent1 or solvent2
		if(group==acc)
		{
			d[0] += 1;
			d[atomt] += 1;
		}
	}

	return(0);
}

//check if the molecule belongs to solvent1 or solvent2
int chkblg(int id, struct MOLE *s1, struct MOLE *s2)
{
	int result,result1,result2;
	int k;

	result=result1=result2=0;

	for(k=s1->idi;k<=s1->idf;k++)
		if(k==id)
		{
			result1=1;
			break;
		}

	for(k=s2->idi;k<=s2->idf;k++)
		if(k==id)
		{
			result2=2;
			break;
		}

	if(result1==1 && result2==0)
		result=1;
	else if(result1==0 && result2==2)
		result=2;
	else if(result1==1 && result2==2)
		result=1;

	return(result);

}

// to find if there is a specific H-bond target(atmt) with water(atmn) at certain snapshot(snap)
int findHbond(int water, int atmt, int atmn, int snap, int **nei, int **att, int **atn, int *nsnaps) 
{
	int result;
	int j;

	if(snap<1 || snap>*nsnaps)
	{
		printf("Error: snapshot %d is out of simulation range\n",snap);
		return(-1);
	}
	
	result=0;
	for(j=1;j<=nei[snap][0];j++)
		if(nei[snap][j]==water && atn[snap][j]==atmn && att[snap][j]==atmt)
			result=1;
	
	return(result);
}

// find the continuous bonded segment starting from snap, get the period ti to tf in array seg[]
int findseg(int seg[2], int water, int atm1, int atm2, int snap, int **nei, int **att, int **atn, int *nsnaps)
{
	int i,j;
	int t;
	int find;

	if(snap>*nsnaps || snap<1)
	{
		printf("ERROR: in finding the time period of a H-bond starting from %d\n\n",snap);
		return(-1);
	}

	i=findHbond(water,atm1,atm2,snap,nei,att,atn,nsnaps);
	if(i==1)
	{
		t=snap;
		find=findHbond(water,atm1,atm2,t,nei,att,atn,nsnaps);
		while(find==1 && t>=1)
		{
			t -= 1;
			if(t<1)
				break;
			else
				find=findHbond(water,atm1,atm2,t,nei,att,atn,nsnaps);
		}
		seg[0]=t+1;
	}
	else if(i==0) // to find the starting snapshot
	{
		t=snap;
		find=findHbond(water,atm1,atm2,t,nei,att,atn,nsnaps);
		while(find==0 && t<=*nsnaps)
		{
			t += 1;
			if(t > *nsnaps)
				break;
			else
				find=findHbond(water,atm1,atm2,t,nei,att,atn,nsnaps);
		}
		seg[0]=t;
	}  

	if(seg[0]>*nsnaps)
		seg[0]=seg[1]=*nsnaps+1;    // if we can not find the next period, then set to "*nsnaps+1"
	else
	{
		t=seg[0];
		find=findHbond(water,atm1,atm2,t,nei,att,atn,nsnaps);
		while(find==1 && t<=*nsnaps) // to find the end snapshot
		{
			t += 1;
			if(t > *nsnaps)
				break;
			else
				find=findHbond(water,atm1,atm2,t,nei,att,atn,nsnaps);
		}
		seg[1]=t-1;
	}

	return(0);
}

// to check if a particular H-bond from snap1 to snap2 is transient bond
int itransbond(int target, int water, int atmt, int atmn, int snap1, int snap2, int dor, int acc,  int **nei, int **att, int **atn, int *nsnaps, struct PARS *ppars)
{
	int j,output;
	int occur,total,left,right;
	float percent;
	int find,k;

	left=(int)(snap1+snap2-ppars->persist[dor][acc])/2;
	right=(int)(snap1+snap2+ppars->persist[dor][acc])/2;

	if(left<=0)
		left=1;

	if(right>*nsnaps)
		right=*nsnaps;

	find=0;occur=0;total=0;
	for(j=left;j<=right;j++)
	{
		find=findHbond(water,atmt,atmn,j,nei,att,atn,nsnaps);
		total += 1;
		if(find==1)
			occur += 1;
	}
	percent = (occur+0.0)/total;

	if(percent <= PERCENT)
		output=1;
	else
		output=0;

	return(output);
}

// find the new H-bonded partner and snapshot with target(atm1)
int findatom(int b[Maxneighb+1], int atm1, int tsnap1, int tsnap2, int **nei, int **att, int **atn, int *nsnaps, int dor, int acc, struct PARS *pars)
{
	int temp,find;
	int j,newp,newa,num,dup,k;
	int ipersist,checki,checkf,ison,sum,count;
	float prob;

	if(tsnap2 >= *nsnaps || tsnap1 < 1)
		return(*nsnaps+1);

	ipersist=pars->persist[dor][acc];   // the check time interval eqauls to persistence value, this is to determine if it is an old partner, 2015.09.25
	checkf=tsnap2;  
	checki=tsnap2-ipersist+1;
	if(checki<1)
		checki=1;

	for(temp=tsnap2+1;temp<=*nsnaps;temp++) // the old partner survives from tsnap1 to tsnap2, check new partner from tsnap2+1
	{
		num=0;find=0;
		for(j=1;j<=nei[temp][0];j++)
		{
			if(att[temp][j]==atm1)
			{
				newp=nei[temp][j];newa=atn[temp][j];

				/* this part from 2015.09.25 */
				count=sum=0;
				for(k=checki;k<=checkf;k++)
				{
					sum += 1;
					ison=findHbond(newp,atm1,newa,k,nei,att,atn,nsnaps); // check if this is an old partner that already occurs at k;
					if(ison==1)
						count += 1;
				}
				if(sum>0)
					prob=(count+0.0)/sum;
				else
					prob=0.0;

				if(prob<0.5)      // the new partner should not be there
				{
					find=1;
					dup=0;
					for(k=1;k<=num;k++)
						if(b[k]==newp)
							dup=1;
					if(dup==0)
					{
						num += 1;
						b[num]=newp;
					}
				}
				/* this part is from 2015.09.24
				isold1=findHbond(newp,atm1,newa,tsnap1,nei,att,atn,nsnaps); // check if this is an old partner that already occurs at tsnap1;
				isold2=findHbond(newp,atm1,newa,tsnap2,nei,att,atn,nsnaps); // check if this is an old partner that already occurs at tsnap2;
				center= (int) ((tsnap1+tsnap2)/2);
				isold3=findHbond(newp,atm1,newa,center,nei,att,atn,nsnaps); // check if this is an old partner that already occurs at center snapshot;

				// if this new partner does not occur at tsnap1, tsnap2, or center, then it is highly possible a new partner, 2015.09.24
				if((isold1==0 && isold2==0) || (isold1==0 && isold3==0) || (isold2==0 && isold3==0)) // at most one identity equals to 1
				{
					find=1;    // this means that we find a new partner
					dup=0;
					for(k=1;k<=num;k++)
						if(b[k]==newp)
							dup=1;
					if(dup==0)
					{
						num += 1;
						b[num]=newp; 
					}
				} */

			}  
		} // end of for(j=1;j<=..)

		if(find==1)
		{
			b[0]=num;
			break;
		}
	}

	if(temp > *nsnaps)
		temp = *nsnaps+1;

	return(temp);
}

//remove a specific H-bond from the original graphgeod
int myremove(int target, int water, int atom1, int atom2, int snap, int **nei, int **att, int **atn)
{
	int find,j,k,i;

	k=nei[snap][0];

	find=0;
	for(j=1;j<=k;j++)
		if(nei[snap][j]==water && att[snap][j]==atom1 && atn[snap][j]==atom2)
		{
			find=1;
			break;
		}

	nei[snap][j]=nei[snap][k];//overwrite the 'j' information (to be removed) with the last HB
	att[snap][j]=att[snap][k];
	atn[snap][j]=atn[snap][k];

	nei[snap][0]=att[snap][0]=atn[snap][0]=k-1;

	return(0);
}

//copy a specific H-bond information from refsnap to snap
int myadd(int target, int water, int atom1, int atom2, int snap, int **nei, int **att, int **atn)
{
	int j,k,result;

	k=nei[snap][0];

	result=0;
	for(j=1;j<=k;j++) // check if this HB information is already there
		if(nei[snap][j]==water && att[snap][j]==atom1 && atn[snap][j]==atom2)
			result=1;
	
	if(result==0)
	{
		k=k+1;
		nei[snap][k]=water;
		att[snap][k]=atom1;
		atn[snap][k]=atom2;
	}

	nei[snap][0]=att[snap][0]=atn[snap][0]=k;

	return(0);

}

//update graphgeod files with one snapshot
int myupdateone(int target, int **nei, int **att, int **atn, int snap)
{

	FILE *fupdate;
	char newname[FLN];
	int j,pbc;

	pbc=0;    // currently, the pbc information is not considered in correction
	sprintf(newname,"Input.water%d.xyz.solB%d.xyz.GraphGeod-updated",snap,snap);
	fupdate=fopen(newname,"a");

	for(j=1;j<=nei[snap][0];j++)
		fprintf(fupdate,"%d %d %d %d %d %d %d ...\n",target,nei[snap][j],pbc,pbc,pbc,att[snap][j],atn[snap][j]);   // it may print duplicated interactions
	
	fclose(fupdate);

	return(0);
}


//update the graphgeod files with the new array information
int myupdate(int target, int **nei, int **att, int **atn, int *nsnaps)
{
	int snap,index;

	for(snap=1;snap<=*nsnaps;snap++)
	{
		index=0;
		index=myupdateone(target,nei,att,atn,snap);
		if(index!=0)
			break;


	}

	return(0);
}

//calculating the distance from the target molecule to a solute, 2016.03.08
float caldist(float **result, int target, int t, int dor, int acc, struct MOLE *m1, struct MOLE *m2, struct MOLE *mt, int *nnodes) 
{
//	float result[12];   // for alkanes, read maximum of 11 Carbons
	float output=0.0;
	FILE *fw,*fs;
	char filename[FLN],buff[FLN],*tok,label,ind1,ind2;
	int index1,natoms,i,j,cnum;
	float watx,waty,watz,solx,soly,solz;  // water Oxygen coordinates and solute first atom coordinates
	float distx,disty,distz;

	for(i=0;i<Ncarbon;i++)
		result[t][i]=0.0;

	sprintf(filename,"water%d.xyz",t);
	if((fw=fopen(filename,"r"))==NULL)
	{
		printf("In function 'caldist': can not find the xyz file at snapshot %d\n",t);
		output=-1.0;
	}
	sprintf(filename,"solB%d.xyz",t);
	if((fs=fopen(filename,"r"))==NULL)
	{
		printf("In function 'caldist': can not find the xyz file at snapshot %d\n",t);
		output=-1.0;
	}

	if(output!=0.0)
	{
		fclose(fw);
		fclose(fs);
		return output;
	}
	else
	{
		if(mt->idi > *nnodes || mt->idf > *nnodes)
			output=-2.0;  // error
		else
		{
			index1=chkblg(target,m1,m2);
			if(index1==1)
				natoms=m1->n;     // I assume water is the first solvent
			else if(index1==2)
			{
				printf("Warning in function 'caldist': the target molecule is solute?\n");
				fclose(fw);
				fclose(fs);
				return(-3.0);
			}

			rewind(fw);
			fgets(buff,sizeof(buff),fw);
			fgets(buff,sizeof(buff),fw);  // skip the head line in xyz file
			for(i=1;i<target;i++)
				for(j=1;j<=natoms;j++)
					fgets(buff,sizeof(buff),fw); // skip the coordinates before target molecule

			// there will be problem with strtok refering to 'NULL' in omp parallel, so change to fscanf, 2016.03.09
//			fgets(buff,sizeof(buff),fw);  // the first atom of target molecule;
//			tok = strtok(buff," \n");   // the atom label
//			tok = strtok(NULL," \n");    // the x coordinate
//			watx = atof(tok);
//			tok = strtok(NULL," \n");    // the y coordinate
//			waty = atof(tok);
//			tok = strtok(NULL," \n");    // the z coordinate
//			watz = atof(tok);
			if(fscanf(fw,"%c%c %f %f %f",&ind1,&ind2,&watx,&waty,&watz)!=5)
			{
				fclose(fw);
				fclose(fs);
				return(-4.0);
			}

			rewind(fs);      // 2015.12.21, I only have one solute molecule
			fgets(buff,sizeof(buff),fs);
			fgets(buff,sizeof(buff),fs);  // skip the head line in xyz file
			cnum=0;
			while(fscanf(fs,"%c%c %f %f %f",&ind1,&ind2,&solx,&soly,&solz)==5)
			{
				fgets(buff,sizeof(buff),fs); // skip the '\n' character in a line
		//		tok = strtok(buff," \n");
				label = ind1;
				if(label == 'C')   // calculate the distance to each Carbon atoms
				{
					cnum += 1;

					// there will be problem with strtok refering to 'NULL' in omp parallel, so change to fscanf, 2016.03.09
//					tok = strtok(NULL," \n");
//					solx = atof(tok);
//					tok = strtok(NULL," \n");
//					soly = atof(tok);
//					tok = strtok(NULL," \n");
//					solz = atof(tok);
			
					distx=watx-solx;
					if(distx > Xsize/2)
						distx = Xsize - distx;
					else if(distx < -Xsize/2)
						distx = Xsize + distx;
					disty=waty-soly;
					if(disty > Ysize/2)
						disty = Ysize - disty;
					else if(disty < - Ysize/2)
						disty = Ysize + disty;
					distz=watz-solz;
					if(distz > Zsize/2)
						distz = Zsize - distz;
					else if(distz < -Zsize/2)
						distz = Zsize + distz;
		
					result[t][cnum] = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
				}
			}
			result[t][0]=cnum;
		}

		fclose(fw);
		fclose(fs);
		return(output);
	}

}

