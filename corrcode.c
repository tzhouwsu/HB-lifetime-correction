/* the source code for the correction algorithm */

#include "head.h"

int corrcode(int dor, int acc, struct MOLE *pm1, struct MOLE *pm2, int *nsolu, struct MOLE *pms, struct PARS *pars, int *nnodes, int *nsnaps, int *ind)
{
	FILE *fop1,*fop2,*fop3,*fop4;
	char fname[FLN];
		
	sprintf(fname,"Degrees-%d-%d",dor,acc);
	fop1=fopen(fname,"w");
	sprintf(fname,"Lifetimes-%d-%d",dor,acc);
	fop2=fopen(fname,"w");
	sprintf(fname,"Switches-%d-%d",dor,acc);
	fop3=fopen(fname,"w");

	sprintf(fname,"degree-data-%d-%d",dor,acc);    // 2015.12.21, add the distance to solute information
	fop4=fopen(fname,"w");

	int i,j,k,m,di;    // dummy labels
	int **neighb,**atmt,**atmn;
	neighb = (int **)malloc(sizeof(int *)*(*nsnaps+1));
	for(i=0;i<=(*nsnaps);i++)
		neighb[i] = (int *)calloc(sizeof(int),Maxneighb+1);
	atmt = (int **)malloc(sizeof(int *)*(*nsnaps+1));
	for(i=0;i<=(*nsnaps);i++)
		atmt[i] = (int *)calloc(sizeof(int),Maxneighb+1);
	atmn = (int **)malloc(sizeof(int *)*(*nsnaps+1));
	for(i=0;i<=(*nsnaps);i++)
		atmn[i] = (int *)calloc(sizeof(int),Maxneighb+1);

	int **label;
	label = (int **)malloc(sizeof(int *)*(MAXNATOM));
	for(i=0;i<=MAXNATOM;i++)
		label[i] = (int *)calloc(sizeof(int),(*nnodes+1));

	int deg[MAXNATOM+1],degrees[Maxneighb+1][MAXNATOM+1],statnum;
	int sym,iapp,symi,syma,symr,symu,tagf,ipt,symd;
	int tgt,number,atomt,atomn,nb,nbn[Maxneighb+1],nbnew,atomnew;
	int period[2],t,t1,t2,t3,t4,tn,ti,tj,tsnap;
	int nbbelong,tgtbelong,tgtidi,tgtidf;
	float **dist,symt,tempdist[Ncarbon];     // 2015.12.21, distance of the target molecule to the solute, at each snapshot
	dist = (float **)malloc(sizeof(float *)*(*nsnaps+1));
	for(i=0;i<=*nsnaps;i++)
		dist[i] = (float *)malloc(sizeof(float)*Ncarbon);

	if(dor==1) // either the first molecule or the second molecule is treated as the target molecule
	{
		tgtidi=pm1->idi;
		tgtidf=pm1->idf;
	}
	else if(dor==2)
	{
		tgtidi=pm2->idi;
		tgtidf=pm2->idf;
	}

	for(i=0;i<=Maxneighb;i++)  // initialize the degrees
		for(j=0;j<=MAXNATOM;j++)
			degrees[i][j]=0;

	for(tgt=tgtidi;tgt<=tgtidf;tgt++)
	{
		printf("  For target %d\n",tgt);

		// re-initialize for every target molecule
		for(t=0;t<=(*nsnaps);t++)
			for(j=0;j<=Maxneighb;j++)
				neighb[t][j]=atmt[t][j]=atmn[t][j]=0;
		for(i=0;i<MAXNATOM;i++)
			for(j=0;j<=*nnodes;j++)
				label[i][j]=0;
		for(t=0;t<=(*nsnaps);t++)
			for(i=0;i<Ncarbon;i++)
				dist[t][i]=0.0;

//		printf("   loading history for %d\n",tgt);
		// load the HB history of the target molecle, which acts as the donor
#pragma omp parallel for shared(neighb,atmt,atmn,tgt)
		for(t=1;t<=(*nsnaps);t++)
		{
			sym+=loadhistory(neighb,atmt,atmn,t,tgt);  // regard solute individually
		//	if(sym!=0)
		//		break;
		}
		if(sym<0)
		{
			printf("In funct corrcode.c: error in loading history\n");
			break;
		}

		// make sure that this target water has HB history
		iapp=0;
		for(t=1;t<=*nsnaps;t++)
			if(neighb[t][0]!=0)
				iapp=1;
		if(iapp!=1) // if this target molecule does not have HB history, skip it
			continue;

		// load the distances to the solute
//		printf("   loading distance for %d\n",tgt);
		tgtbelong=chkblg(tgt,pm1,pm2);
		if(tgtbelong == 1)
		{	
#pragma omp parallel for shared(dist,tgt,dor,acc,pm1,pm2,pms,nnodes)
			for(t=1;t<=(*nsnaps);t++)
			{
				symt+=caldist(dist,tgt,t,dor,acc,pm1,pm2,pms,nnodes);
//				if(symt<0.0) // if there is error
//					break;
			}
			if(symt<0)
			{
				printf("In corrcode.c: error in loading distances\n");
				break;
			}
		}

//		printf("   calculate degrees for %d\n",tgt);
		// count the degrees of the target solvent
		for(t=1;t<=(*nsnaps);t++)
		{
			for(j=0;j<=MAXNATOM;j++)
				deg[j]=0;

			symd=caldegrees(deg,t,dor,acc,pm1,pm2,neighb,atmt,atmn);
			if(symd!=0)
				break;
			else
			{
				for(j=0;j<=MAXNATOM;j++)
					degrees[deg[j]][j] += 1;
			}
	//		if(deg[0]>7)    // I found a bug with duplicated degree counting, 2015-12-10
	//			for(j=0;j<=MAXNATOM;j++)
	//				printf("    %d %d %d %d\n",tgt,t,j,deg[j]);

			if(*ind == 0)
			{
				fprintf(fop4,"%d :",deg[0]);
				if(dor==1)
				{
					for(j=1;j<=pm1->n;j++)
						fprintf(fop4," %d",deg[j]);
				}
				else if(dor==2)
				{
					for(j=1;j<=pm2->n;j++)
						fprintf(fop4," %d",deg[j]);
				}
				fprintf(fop4," : ");
				for(j=1;j<=dist[t][0];j++)
					fprintf(fop4,"%f ",dist[t][j]);
				fprintf(fop4,"\n");
			}

		}

		if(symd!=0 || symt<0.0)   // if there is some error in calculating degrees and distances, the code stops
			break;

//		printf("   calculating lifetime %d\n",tgt);

		// calculate the lifetimes associated with target molecule, there are three modes,
		tagf=syma=symr=symi=symu=0; // initialize all the error indicators
		if(*ind==0)   // mode0, analyze the result, treat all bonds and breaks as true real events, decompose HB lifetime incormation with respect to different mechanisms
		{
			for(t=1;t<=(*nsnaps);t++)
			{
				number=neighb[t][0];
				for(j=1;j<=number;j++)
				{
					atomt=atmt[t][j];atomn=atmn[t][j];nb=neighb[t][j];
					nbbelong=chkblg(nb,pm1,pm2); // check wether the molecule 'nb' belongs to solvent1 or solvent2
					tgtbelong=chkblg(tgt,pm1,pm2);
					if(tgtbelong != dor || nbbelong != acc)  // only check for the specific dor:acc interaction
						continue;
	
					if(label[atomt][nb]<t) // this interaction has not been checked before, otherwise skip
					{
						period[0]=period[1]=0;
						tagf=findseg(period,nb,atomt,atomn,t,neighb,atmt,atmn,nsnaps); // pick out this interation interval
						if(tagf!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop j
							break;
						t1=period[0];t2=period[1];
						if(t2==*nsnaps)
							t3=t4=*nsnaps+100;   // set this number larger enough
						else
						{
							period[0]=period[1]=0;
							tagf=findseg(period,nb,atomt,atomn,t2+1,neighb,atmt,atmn,nsnaps); // pick out the next interation interval
							if(tagf!=0)
								break;
							t3=period[0];t4=period[1];

							tn=*nsnaps+1;
							for(k=0;k<=Maxneighb;k++)
								nbn[k]=0;
							tn=findatom(nbn,atomt,t1,t2,neighb,atmt,atmn,nsnaps,dor,acc,pars); // find the new H-bonded partnera with respect to target atom

						}

						symi=itransbond(tgt,nb,atomt,atomn,t1,t2,dor,acc,neighb,atmt,atmn,nsnaps,pars); // determine if it is transint bond or not

						if(t3>*nsnaps && tn>*nsnaps) // this interaction (t1-t2) does not reoccur, neither it switches partner, cannot determine its mechanism
						{
							if(symi==1) // this interaction (t1-t2) is transient bond, and no information after it
							{
								fprintf(fop2,"%d %d %d %d %d %d T U ",tgt,nb,atomt,atomn,t1,t2); // T means transient, U means undetermined
								fprintf(fop3,"T U %d %d %d %d %d %d %d %d ",tgt,nb,nb,atomt,t1,t2,*nsnaps+1,*nsnaps+1);
								for(di=1;di<=dist[t][0];di++)   // loop over all carbon site, i.e. distance from target water to every carbon
								{
									fprintf(fop2,"%f ",dist[t][di]);
									fprintf(fop3,"%f ",dist[t][di]);
								}
								fprintf(fop2," \n");
								fprintf(fop3," \n");
							}
							else if(symi==0)  // it is not transient bond, yet no information after it
							{
								fprintf(fop2,"%d %d %d %d %d %d R U ",tgt,nb,atomt,atomn,t1,t2); // R means real, U means undetermined
								fprintf(fop3,"R U %d %d %d %d %d %d %d %d ",tgt,nb,nb,atomt,t1,t2,*nsnaps+1,*nsnaps+1);
								for(di=1;di<=dist[t][0];di++)
								{
									fprintf(fop2,"%f ",dist[t][di]);
									fprintf(fop3,"%f ",dist[t][di]);
								}
								fprintf(fop2," \n");
								fprintf(fop3," \n");
							}
						}
						else    // there's some interaction after t1-t2
						{
							if(t3<=tn && t3<*nsnaps) // this interaction (t1-t2) reoccurs
							{
								if(symi==1)  // this interaction(t1-t2) is transient bond, and it reoccurs at t3
								{
									fprintf(fop2,"%d %d %d %d %d %d T N ",tgt,nb,atomt,atomn,t1,t2); // T means transient, N means it does not change partner
									fprintf(fop3,"T N %d %d %d %d %d %d %d %d ",tgt,nb,nb,atomt,t1,t2,t3,t4);
									for(di=1;di<=dist[t][0];di++)
									{
										fprintf(fop2,"%f ",dist[t][di]);
										fprintf(fop3,"%f ",dist[t][di]);
									}
									fprintf(fop2," \n");
									fprintf(fop3," \n");
								}
								else if(symi==0)
								{
									fprintf(fop2,"%d %d %d %d %d %d R N ",tgt,nb,atomt,atomn,t1,t2); // R means real, N means it does not change partner
									fprintf(fop3,"R N %d %d %d %d %d %d %d %d ",tgt,nb,nb,atomt,t1,t2,t3,t4);
									for(di=1;di<=dist[t][0];di++)
									{
										fprintf(fop2,"%f ",dist[t][di]);
										fprintf(fop3,"%f ",dist[t][di]);
									}
									fprintf(fop2," \n");
									fprintf(fop3," \n");
								}
							}

							if(tn<t3 && tn<*nsnaps) // new partner approaches
							{
								ipt=nbn[0];
								for(k=1;k<=nbn[0];k++)   // loop over all the new partners, it is highly possible that there's only one new partner
								{
									nbnew=nbn[k]; // the new partner id
									for(m=1;m<=neighb[tn][0];m++)
										if(nbnew==neighb[tn][m])
										{
											atomnew=atmn[tn][m];  // the atom from the new partner
											break;
										}

									period[0]=period[1]=0;
									tagf=findseg(period,nbnew,atomt,atomnew,tn,neighb,atmt,atmn,nsnaps); 
									if(tagf!=0)   // if there's some error, the code should stop, here it breaks from loop k
										break;

									ti=period[0];tj=period[1];

									if(symi==1)
									{
										if(ipt>0) // the indicator for printing the Lifetime information, since t1-t2 lifetime should only be printed once
										{
											fprintf(fop2,"%d %d %d %d %d %d T Y ",tgt,nb,atomt,atomn,t1,t2); // T means transient, Y means it does not change partner
											for(di=1;di<=dist[t][0];di++)
											{
												fprintf(fop2,"%f ",dist[t][di]);
											}
											fprintf(fop2," \n");
											ipt=0;
										}
										fprintf(fop3,"T Y %d %d %d %d %d %d %d %d ",tgt,nb,nbnew,atomt,t1,t2,ti,tj);
										for(di=1;di<=dist[t][0];di++)
										{
											fprintf(fop3,"%f ",dist[t][di]);
										}
										fprintf(fop3," \n");
									}
									else if(symi==0)
									{
										if(ipt>0) // the indicator for printing the Lifetime information, since t1-t2 lifetime should only be printed once
										{
											fprintf(fop2,"%d %d %d %d %d %d R Y ",tgt,nb,atomt,atomn,t1,t2); // R means transient, Y means it does not change partner
											for(di=1;di<=dist[t][0];di++)
											{
												fprintf(fop2,"%f ",dist[t][di]);
											}
											fprintf(fop2," \n");
											ipt=0;
										}
										fprintf(fop3,"R Y %d %d %d %d %d %d %d %d ",tgt,nb,nbnew,atomt,t1,t2,ti,tj);
										for(di=1;di<=dist[t][0];di++)
										{
											fprintf(fop3,"%f ",dist[t][di]);
										}
										fprintf(fop3," \n");
									}
								} // end of for(k=1;k<=nbn[0];k++)
							} // end of if(tn<t3 && tn<*nsnaps) // new partner approaches
						} // end of  else  // there's some interaction after t1-t2

						label[atomt][nb]=t2; // set the label to the lastest processed time-snap
						
					} // end of if(label[atomt][nb]<t)
					if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop j
						break;
				}  // end of for(j=1;j<=number;j++)
				if(tagf!=0 || syma!=0 || symr!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop t
					break;
			} // end of for(t=1;t<=(*nsnaps);t++)
			if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop tgt
				break;


		} // end of mode0
		else if(*ind==1) // mode1, correct transient breaks
		{
			for(t=1;t<=(*nsnaps);t++)
			{
				number=neighb[t][0];
				for(j=1;j<=number;j++)
				{
					atomt=atmt[t][j];atomn=atmn[t][j];nb=neighb[t][j];
					nbbelong=chkblg(nb,pm1,pm2); // check wether the molecule 'nb' belongs to solvent1 or solvent2
					tgtbelong=chkblg(tgt,pm1,pm2);
					if(tgtbelong != dor || nbbelong != acc)  // only check for the specific dor:acc interaction
						continue;
	
					if(label[atomt][nb]<t) // this interaction has not been checked before, otherwise skip
					{
						period[0]=period[1]=0;
						tagf=findseg(period,nb,atomt,atomn,t,neighb,atmt,atmn,nsnaps); // pick out this interation interval
						if(tagf!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop j
							break;
						t1=period[0];t2=period[1];
						if(t2==*nsnaps)
							t3=t4=*nsnaps+100; // set this number large enough
						else
						{
							period[0]=period[1]=0;
							tagf=findseg(period,nb,atomt,atomn,t2+1,neighb,atmt,atmn,nsnaps); // pick out the next interation interval
							if(tagf!=0)
								break;
							t3=period[0];t4=period[1];
						}

						if(t3<*nsnaps && t3-t2 < pars->tolerance[dor][acc]) // do correction if the break duration t3-t2 is less than the TOLERENCE
						{
							for(tsnap=t2+1;tsnap<t3;tsnap++) // modify the array here, then update the new GraphGeod later
							{
								syma=myadd(tgt,nb,atomt,atomn,tsnap,neighb,atmt,atmn);
								if(syma!=0) // if there's some error in calling this function, the code should stop
									break;
							}
						}

						label[atomt][nb]=t2; // set the label to the lastest processed time-snap
					} // end of if(label[atomt][nb]<t)

					if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop j
						break;
				}  // end of for(j=1;j<=number;j++)
				if(tagf!=0 || syma!=0 || symr!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop t
					break;
			} // end of for(t=1;t<=(*nsnaps);t++)
			if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop tgt
				break;

		} // end of mode1
		else if(*ind==2) // mode2, correct transient bonds
		{
			for(t=1;t<=(*nsnaps);t++)
			{
				number=neighb[t][0];
				for(j=1;j<=number;j++)
				{
					atomt=atmt[t][j];atomn=atmn[t][j];nb=neighb[t][j];
					nbbelong=chkblg(nb,pm1,pm2); // check wether the molecule 'nb' belongs to solvent1 or solvent2
					tgtbelong=chkblg(tgt,pm1,pm2);
					if(tgtbelong != dor || nbbelong != acc)  // only check for the specific dor:acc interaction
						continue;
	
					if(label[atomt][nb]<t) // this interaction has not been checked before, otherwise skip
					{
						period[0]=period[1]=0;
						tagf=findseg(period,nb,atomt,atomn,t,neighb,atmt,atmn,nsnaps); // pick out this interation interval
						if(tagf!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop j
							break;
						t1=period[0];t2=period[1];

						symi=itransbond(tgt,nb,atomt,atomn,t1,t2,dor,acc,neighb,atmt,atmn,nsnaps,pars); // determine if it is transint bond or not
						if(symi==1) // this interaction, t1-t2, is a transient bond
						{
							for(tsnap=t1;tsnap<=t2;tsnap++)
							{
								symr=myremove(tgt,nb,atomt,atomn,tsnap,neighb,atmt,atmn); // modify the array, update the GraphGeod file later
								if(symr!=0)  // if there's some error in calling this function, the code should stop
									break;
							}
						}

						label[atomt][nb]=t2; // set the label to the lastest processed time-snap
					} // end of if(label[atomt][nb]<t)

					if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop j
						break;
				}  // end of for(j=1;j<=number;j++)
				if(tagf!=0 || syma!=0 || symr!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop t
					break;
			} // end of for(t=1;t<=(*nsnaps);t++)
			if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop tgt
				break;

		} // end of mode2
		else if(*ind==3) // mode3, correct both transient breaks and transient bonds
		{
			for(t=1;t<=(*nsnaps);t++)
			{
				number=neighb[t][0];
				for(j=1;j<=number;j++)
				{
					atomt=atmt[t][j];atomn=atmn[t][j];nb=neighb[t][j];
					nbbelong=chkblg(nb,pm1,pm2); // check wether the molecule 'nb' belongs to solvent1 or solvent2
					tgtbelong=chkblg(tgt,pm1,pm2);
					if(tgtbelong != dor || nbbelong != acc)  // only check for the specific dor:acc interaction
						continue;
	
					if(label[atomt][nb]<t) // this interaction has not been checked before, otherwise skip
					{
						period[0]=period[1]=0;
						tagf=findseg(period,nb,atomt,atomn,t,neighb,atmt,atmn,nsnaps); // pick out this interation interval
						if(tagf!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop j
							break;
						t1=period[0];t2=period[1];
						if(t2==*nsnaps)
							t3=t4=*nsnaps+100; // set this number large enough
						else
						{
							period[0]=period[1]=0;
							tagf=findseg(period,nb,atomt,atomn,t2+1,neighb,atmt,atmn,nsnaps); // pick out the next interation interval
							if(tagf!=0)
								break;
							t3=period[0];t4=period[1];
						}

						if(t3<*nsnaps && t3-t2 < pars->tolerance[dor][acc]) // do correction if the break duration t3-t2 is less than the TOLERENCE
						{
							for(tsnap=t2+1;tsnap<t3;tsnap++) // modify the array here, then update the new GraphGeod later
							{
								syma=myadd(tgt,nb,atomt,atomn,tsnap,neighb,atmt,atmn);
								if(syma!=0) // if there's some error in calling this function, the code should stop
									break;
							}
						}

						// here I correct the transient bonds after the transient breaks, so by default: first transient breaks and then transient bonds
						symi=itransbond(tgt,nb,atomt,atomn,t1,t2,dor,acc,neighb,atmt,atmn,nsnaps,pars); // determine if it is transint bond or not
						if(symi==1) // this interaction, t1-t2, is a transient bond
						{
							for(tsnap=t1;tsnap<=t2;tsnap++)
							{
								symr=myremove(tgt,nb,atomt,atomn,tsnap,neighb,atmt,atmn); // modify the array, update the GraphGeod file later
								if(symr!=0)  // if there's some error in calling this function, the code should stop
									break;
							}
						}

						label[atomt][nb]=t2; // set the label to the lastest processed time-snap
						
					} // end of if(label[atomt][nb]<t)

					if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop j
						break;
				}  // end of for(j=1;j<=number;j++)
				if(tagf!=0 || syma!=0 || symr!=0)   // if there's some error in calling this function, the code should stop, here it breaks from loop t
					break;
			} // end of for(t=1;t<=(*nsnaps);t++)
			if(tagf!=0 || syma!=0 || symr!=0) // if there's some error in calling this function, the code should stop, here it breaks from loop tgt
				break;

		} // end of mode3

		if(*ind==1 || *ind==2 || *ind==3)  // update the GraphGeod files, for correction
		{
#pragma omp parallel for shared(neighb,atmt,atmn,tgt)	
			for(tsnap=1;tsnap<=*nsnaps;tsnap++)
				symu+=myupdateone(tgt,neighb,atmt,atmn,tsnap);  // update the GraphGeod file based on the new arrays
			if(symu<0)  // if there's some error in calling this function, it should stop, here it breaks from loop tgt
				break;
		}

	} // end of for(tgt...

	// print the degrees

	if(*ind==0)
	{
		fprintf(fop1,"\nDegrees distribution for solvent%d from solvent%d :\n\n\n",dor,acc);
		for(j=0;j<=MAXNATOM;j++)
		{
			for(i=0;i<Maxneighb;i++)
				fprintf(fop1,"atom %d degree %d observation %d\n",j,i,degrees[i][j]);
			fprintf(fop1,"\n\n");
		}
	}
	else if(*ind==1 || *ind==2 || *ind==3)
	{
		fprintf(fop1,"\nNo Degrees information printed with mode %d\n\n",*ind); // mode1 mode2 and mode3, do correction, and no output information
		fprintf(fop2,"\nNo Lifetimes information printed with mode %d\n\n",*ind);
		fprintf(fop3,"\nNo Switches information printed with mode %d\n\n",*ind);
	}

	fclose(fop1);
	fclose(fop2);
	fclose(fop3);
	fclose(fop4);

	for(i=0;i<=(*nsnaps);i++)
	{
		free(neighb[i]);
		free(atmt[i]);
		free(atmn[i]);
	}
	free(neighb);
	free(atmt);
	free(atmn);

	for(i=0;i<MAXNATOM;i++)
		free(label[i]);
	free(label);


	return(0);

}


