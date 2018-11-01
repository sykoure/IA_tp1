#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define PROFONDEUR_MAX 5
#define NB_LIGNES 10
#define NB_COLONNES 10
#define NB_PION_MAX 14
#define VALEUR_PION_MAX 22
#define ALPHA 0.4			//Nombre de pion  0.4
#define BRAVO 0.25		//Distance moyenne 0.25
#define CHARLIE 0.1		//Distance minimale 0.1
#define DELTA 0.25				//Somme valeur pion

//TODO CHANGER L'ORDRE DE LECTURE DOUBLE FOR
//#define DEBUG


typedef struct pion_s
{
	int couleur;
	int valeur;
}Pion;

Pion *plateauDeJeu;
int nb_appels = 0;
int nb_noeuds_visites = 0;
int nb_noeuds_totals = 0;
double total_time = 0.0;

void f_affiche_plateau(Pion *plateau);
int f_convert_char2int(char c);
char f_convert_int2char(int i);
Pion *f_init_plateau();
void data();
double f_max(Pion *plateau, int joueur, int profondeur, int *l1, int *c1,int *l2,int *c2, double* alpha, double* beta);
double f_min(Pion *plateau, int joueur, int profondeur, int *l1, int *c1,int *l2,int *c2, double* alpha, double* beta);

void data(Pion *plateau,double tabEuristique[8]){
	double distance_noire_min = FLT_MAX;
	double distance_blanche_min = FLT_MAX;
	double distance_noire_moyenne = 0;
	double distance_blanche_moyenne = 0;
	double nombre_noire = 0;
	double nombre_blanche = 0;
	int i,j;
	double somme_blanche = 0;
	double somme_noire = 0;
	double somme_valeur_blanche = 0;
	double somme_valeur_noire = 0;

	//f_affiche_plateau(plateau);

	for(i = 0;i < NB_LIGNES;i++){
		for(j = 0;j < NB_COLONNES;j++){
			if(plateau[i*NB_COLONNES+j].couleur == -1){
				nombre_blanche++;
				somme_valeur_blanche = somme_valeur_blanche + plateau[i*NB_COLONNES+j].valeur;
				somme_blanche = somme_blanche + (NB_LIGNES-1-i);
				if(distance_blanche_min >= NB_LIGNES-1-i){
					distance_blanche_min = NB_LIGNES-1-i;
				}
			}
			else if(plateau[i*NB_COLONNES+j].couleur == 1){
				nombre_noire++;
				somme_valeur_noire = somme_valeur_noire + plateau[i*NB_COLONNES+j].valeur;
				somme_noire = somme_noire + i;
				if(distance_noire_min >= i){
					distance_noire_min = i;
				}
			}
		}
	}
	distance_blanche_moyenne = somme_blanche / nombre_blanche;
	distance_noire_moyenne = somme_noire / nombre_noire;

	double tab[8] = {nombre_noire,distance_noire_moyenne,distance_noire_min,nombre_blanche,distance_blanche_moyenne,
	distance_blanche_min,somme_valeur_blanche,somme_valeur_noire};

	#ifdef DEBUG
		printf("Nombre de pion noirs restants : %f\n",tab[0]);
		printf("Distance moyenne des pions noirs vers l'arrivée : %f\n",tab[1]);
		printf("Distance minimale d'un pion noir vers l'arrivée : %f\n",tab[2]);
		printf("Nombre de pion blancs restants : %f\n",tab[3]);
		printf("Distance moyenne des pions blancs vers l'arrivée : %f\n",tab[4]);
		printf("Distance minimale d'un pion blanc vers l'arrivée : %f\n",tab[5]);
		printf("Somme valeur des pions blancs : %f\n",tab[6]);
		printf("Somme valeur des pions noirs : %f\n",tab[7]);
	#endif

	for(i = 0;i < 8;i++){
		tabEuristique[i] = tab[i];
	}
}

int f_convert_char2int(char c)
{

	if(c>='A' && c<='Z')
		return (int)(c-'A');
	if(c>='a' && c<='z')
		return (int)(c-'a');
	return -1;
}

char f_convert_int2char(int i)
{


	return (char)i+'A';
}

Pion *f_init_plateau0()
{
	int i, j;
	Pion *plateau=NULL;




	plateau = (Pion *)malloc(NB_LIGNES*NB_COLONNES*sizeof(Pion));
	if(plateau == NULL)
	{
		printf("error: unable to allocate memory\n");
		exit(EXIT_FAILURE);
	}

	for(i=0; i<NB_LIGNES; i++)
	{
		for(j=0; j<NB_COLONNES; j++)
		{
			plateau[i*NB_COLONNES+j].couleur = 0;
			plateau[i*NB_COLONNES+j].valeur = 0;
		}
	}
	return plateau;
}

Pion *f_init_plateau()
{
	int i, j;
	Pion *plateau=NULL;



	plateau = (Pion *)malloc(NB_LIGNES*NB_COLONNES*sizeof(Pion));
	if(plateau == NULL)
	{
		printf("error: unable to allocate memory\n");
		exit(EXIT_FAILURE);
	}

	for(i=0; i<NB_LIGNES; i++)
	{
		for(j=0; j<NB_COLONNES; j++)
		{
			plateau[i*NB_COLONNES+j].couleur = 0;
			plateau[i*NB_COLONNES+j].valeur = 0;
		}
	}

	plateau[9*NB_COLONNES+5].couleur = 1;
	plateau[9*NB_COLONNES+5].valeur = 1;

	plateau[9*NB_COLONNES+6].couleur = 1;
	plateau[9*NB_COLONNES+6].valeur = 2;

	plateau[9*NB_COLONNES+7].couleur = 1;
	plateau[9*NB_COLONNES+7].valeur = 3;

	plateau[9*NB_COLONNES+8].couleur = 1;
	plateau[9*NB_COLONNES+8].valeur = 2;

	plateau[9*NB_COLONNES+9].couleur = 1;
	plateau[9*NB_COLONNES+9].valeur = 1;

	plateau[8*NB_COLONNES+0].couleur = 1;
	plateau[8*NB_COLONNES+0].valeur = 1;

	plateau[8*NB_COLONNES+1].couleur = 1;
	plateau[8*NB_COLONNES+1].valeur = 3;

	plateau[8*NB_COLONNES+2].couleur = 1;
	plateau[8*NB_COLONNES+2].valeur = 3;

	plateau[8*NB_COLONNES+3].couleur = 1;
	plateau[8*NB_COLONNES+3].valeur = 1;

	plateau[8*NB_COLONNES+6].couleur = 1;
	plateau[8*NB_COLONNES+6].valeur = 1;

	plateau[8*NB_COLONNES+7].couleur = 1;
	plateau[8*NB_COLONNES+7].valeur = 1;

	plateau[8*NB_COLONNES+8].couleur = 1;
	plateau[8*NB_COLONNES+8].valeur = 1;

	plateau[7*NB_COLONNES+1].couleur = 1;
	plateau[7*NB_COLONNES+1].valeur = 1;

	plateau[7*NB_COLONNES+2].couleur = 1;
	plateau[7*NB_COLONNES+2].valeur = 1;

	plateau[2*NB_COLONNES+7].couleur = -1;
	plateau[2*NB_COLONNES+7].valeur = 1;

	plateau[2*NB_COLONNES+8].couleur = -1;
	plateau[2*NB_COLONNES+8].valeur = 1;

	plateau[1*NB_COLONNES+1].couleur = -1;
	plateau[1*NB_COLONNES+1].valeur = 1;

	plateau[1*NB_COLONNES+2].couleur = -1;
	plateau[1*NB_COLONNES+2].valeur = 1;

	plateau[1*NB_COLONNES+3].couleur = -1;
	plateau[1*NB_COLONNES+3].valeur = 1;

	plateau[1*NB_COLONNES+6].couleur = -1;
	plateau[1*NB_COLONNES+6].valeur = 1;

	plateau[1*NB_COLONNES+7].couleur = -1;
	plateau[1*NB_COLONNES+7].valeur = 3;

	plateau[1*NB_COLONNES+8].couleur = -1;
	plateau[1*NB_COLONNES+8].valeur = 3;

	plateau[1*NB_COLONNES+9].couleur = -1;
	plateau[1*NB_COLONNES+9].valeur = 1;

	plateau[0*NB_COLONNES+0].couleur = -1;
	plateau[0*NB_COLONNES+0].valeur = 1;

	plateau[0*NB_COLONNES+1].couleur = -1;
	plateau[0*NB_COLONNES+1].valeur = 2;

	plateau[0*NB_COLONNES+2].couleur = -1;
	plateau[0*NB_COLONNES+2].valeur = 3;

	plateau[0*NB_COLONNES+3].couleur = -1;
	plateau[0*NB_COLONNES+3].valeur = 2;

	plateau[0*NB_COLONNES+4].couleur = -1;
	plateau[0*NB_COLONNES+4].valeur = 1;


return plateau;
}

void f_affiche_plateau(Pion *plateau)
{
	int i,j,k;



	printf("\n    ");
	for(k=0; k<NB_COLONNES; k++)
		printf("%2c ",f_convert_int2char(k));
	printf("\n    ");
	for(k=0; k<NB_COLONNES; k++)
		printf("-- ");
	printf("\n");
	for(i=NB_LIGNES-1; i>=0; i--)
	{
		printf("%2d ", i);
		for(j=0; j<NB_COLONNES; j++)
		{
			printf("|");
			switch(plateau[i*NB_COLONNES+j].couleur)
			{
			case -1:
				printf("%do",plateau[i*NB_COLONNES+j].valeur);
				break;
			case 1:
				printf("%dx",plateau[i*NB_COLONNES+j].valeur);
				break;
			default:
				printf("  ");
			}
		}
		printf("|\n    ");
		for(k=0; k<NB_COLONNES; k++)
			printf("-- ");
		printf("\n");
	}
	printf("    ");

}


int f_gagnant()
{
	int i, j, somme1=0, somme2=0;



	//Quelqu'un est-il arrive sur la ligne de l'autre
	for(i=0; i<NB_COLONNES; i++)
	{
		if(plateauDeJeu[i].couleur == 1)
			return 1;
		if(plateauDeJeu[(NB_LIGNES-1)*NB_COLONNES+i].couleur == -1)
			return -1;
	}

	//taille des armees
	for(i=0; i<NB_LIGNES; i++)
	{
		for(j=0; j<NB_COLONNES; j++)
		{
			if(plateauDeJeu[i*NB_COLONNES+j].couleur == 1)
				somme1++;
			if(plateauDeJeu[i*NB_COLONNES+j].couleur == -1)
				somme2++;
		}
	}
	if(somme1==0)
		return -1;
	if(somme2==0)
		return 1;

	return 0;
}


/**
 * Prend comme argument la ligne et la colonne de la case
 * 	pour laquelle la bataille a lieu
 * Renvoie le couleur du gagnant
 * */
int f_bataille(int l, int c)
{
	int i, j, mini, maxi, minj, maxj;
	int somme=0;
	mini = l-1<0?0:l-1;
	maxi = l+1>NB_LIGNES-1?NB_LIGNES-1:l+1;
	minj = c-1<0?0:c-1;
	maxj = c+1>NB_COLONNES-1?NB_COLONNES-1:c+1;

	for(i=mini; i<=maxi; i++)
	{
		for(j=minj; j<=maxj; j++)
		{
			somme += plateauDeJeu[i*NB_COLONNES+j].couleur*plateauDeJeu[i*NB_COLONNES+j].valeur;
		}
	}
	somme -= plateauDeJeu[l*NB_COLONNES+c].couleur*plateauDeJeu[l*NB_COLONNES+c].valeur;

	if(somme < 0)
		return -1;
	if(somme > 0)
		return 1;

	return plateauDeJeu[l*NB_COLONNES+c].couleur;
}


/**
 * Prend la ligne et colonne de la case d'origine
 * 	et la ligne et colonne de la case de destination
 * Renvoie 1 en cas d'erreur
 * Renvoie 0 sinon
 * */
int f_test_mouvement(Pion *plateau, int l1, int c1, int l2, int c2, int couleur)
{
#ifdef DEBUG
	printf("de (%d,%d) vers (%d,%d)\n", l1, c1, l2, c2);
#endif
	/* Erreur, hors du plateau */
	if(l1 < 0 || l1 >= NB_LIGNES || l2 < 0 || l2 >= NB_LIGNES ||
			c1 < 0 || c1 >= NB_COLONNES || c2 < 0 || c2 >= NB_COLONNES)
		return 1;
	/* Erreur, il n'y a pas de pion a deplacer ou le pion n'appartient pas au joueur*/
	if(plateau[l1*NB_COLONNES+c1].valeur == 0 || plateau[l1*NB_COLONNES+c1].couleur != couleur)
		return 1;
	/* Erreur, tentative de tir fratricide */
	if(plateau[l2*NB_COLONNES+c2].couleur == plateau[l1*NB_COLONNES+c1].couleur)
		return 1;

	if(l1-l2 >1 || l2-l1 >1 || c1-c2 >1 || c2-c1 >1 || (l1==l2 && c1==c2))
		return 1;
	return 0;
}


/**
 * Prend la ligne et colonne de la case d'origine
 * 	et la ligne et colonne de la case de destination
 *  et effectue le trantement de l'operation demandée
 * Renvoie 1 en cas d'erreur
 * Renvoie 0 sinon
 * */
int f_bouge_piece(Pion *plateau, int l1, int c1, int l2, int c2, int couleur)
{
	int gagnant=0;


	if(f_test_mouvement(plateau, l1, c1, l2, c2, couleur) != 0)
		return 1;


	/* Cas ou il n'y a personne a l'arrivee */
	if(plateau[l2*NB_COLONNES+c2].valeur == 0)
	{
		plateau[l2*NB_COLONNES+c2].couleur = plateau[l1*NB_COLONNES+c1].couleur;
		plateau[l2*NB_COLONNES+c2].valeur = plateau[l1*NB_COLONNES+c1].valeur;
		plateau[l1*NB_COLONNES+c1].couleur = 0;
		plateau[l1*NB_COLONNES+c1].valeur = 0;
	}
	else
	{
		gagnant=f_bataille(l2, c2);
		/* victoire */
		if(gagnant == couleur)
		{
			plateau[l2*NB_COLONNES+c2].couleur = plateau[l1*NB_COLONNES+c1].couleur;
			plateau[l2*NB_COLONNES+c2].valeur = plateau[l1*NB_COLONNES+c1].valeur;
			plateau[l1*NB_COLONNES+c1].couleur = 0;
			plateau[l1*NB_COLONNES+c1].valeur = 0;
		}
		/* defaite */
		else if(gagnant != 0)
		{
			plateau[l1*NB_COLONNES+c1].couleur = 0;
			plateau[l1*NB_COLONNES+c1].valeur = 0;
		}
	}

	return 0;
}

//Calcul du nombre de pions sur le plateau du joueur
int f_nbPions(Pion* jeu, int joueur)
{
	int nbPion=0;
	int i, j;
	for (i = 0; i < NB_COLONNES; ++i)
	{
		for (j = 0; j < NB_LIGNES; ++j)
		{
			if (jeu[i * NB_COLONNES + j].couleur == joueur)
			{
				++nbPion;
			}
		}
	}
	return nbPion;
}

//Calcul de la valeur de tous les pions du joueur
int f_valeur(Pion* jeu, int joueur)
{
	int i, j;
	int valeur=0;
	for (i = 0; i < NB_COLONNES; ++i)
	{
		for (j = 0; j < NB_LIGNES; ++j)
		{
			if (jeu[i * NB_COLONNES + j].couleur == joueur)
			{
				valeur += jeu[i * NB_COLONNES + j].valeur;
			}
		}
	}
	return valeur;
}

//fonction d'évaluation
double f_eval(Pion * plateau,int joueur)
{
	double tabEuristique[8];
	data(plateau,tabEuristique);
	double valeur = 0;
	//Joueur pion noir
	if(joueur == -1){
		tabEuristique[0] = ALPHA * tabEuristique[0] / NB_PION_MAX;
		tabEuristique[7] = BRAVO * tabEuristique[7] / VALEUR_PION_MAX;
		tabEuristique[1] = CHARLIE * tabEuristique[1] / (NB_LIGNES-1);
		tabEuristique[2] = DELTA * tabEuristique[2] / (NB_LIGNES - 1);
		valeur = tabEuristique[0] + tabEuristique[7] -tabEuristique[1] -tabEuristique[2];
		#ifdef DEBUG
		printf("Valeur du pion noir : %f\n",valeur);
		#endif
	}
	//Joueur pion blanc
	else{
		tabEuristique[3] = ALPHA * tabEuristique[3] / NB_PION_MAX;
		tabEuristique[6] = BRAVO * tabEuristique[6] / VALEUR_PION_MAX;
		tabEuristique[4] = CHARLIE * tabEuristique[4] / (NB_LIGNES-1);
		tabEuristique[5] = DELTA * tabEuristique[5] / (NB_LIGNES - 1);
		valeur = tabEuristique[3] + tabEuristique[6] - tabEuristique[4] -tabEuristique[5];
		#ifdef DEBUG
			printf("Valeur du pion blanc : %f\n",valeur);
		#endif
	}
	return valeur;
}

//copie du plateau
void f_copie_plateau(Pion* source, Pion* destination)
{
	int i, j;
	for (i = 0; i < NB_LIGNES; i++)
	{
		for (j = 0; j < NB_COLONNES; j++)
		{
			destination[i * NB_COLONNES + j].couleur = source[i * NB_COLONNES + j].couleur;
			destination[i * NB_COLONNES + j].valeur = source[i * NB_COLONNES + j].valeur;
		}
	}
}

//mise a zero du plateau
Pion* f_raz_plateau()
{
	Pion* jeu = NULL;
	int i, j;
	jeu = (Pion *) malloc(NB_LIGNES * NB_COLONNES * sizeof (Pion));
	for (i = 0; i < NB_LIGNES; i++)
	{
		for (j = 0; j < NB_COLONNES; j++)
		{
			jeu[i * NB_COLONNES + j].couleur = 0;
			jeu[i * NB_COLONNES + j].valeur = 0;
		}
	}
	return jeu;
}

//Fonction min trouve le minimum des noeuds fils
double f_min(Pion *plateau, int joueur, int profondeur, int *l1, int *c1,int *l2,int *c2, double* alpha, double* beta)
{
	int ll1 = *l1;
	int lc1 = *c1;
	int ll2 = *l2;
	int lc2 = *c2;
	double lbeta = *beta;
	int i, j, k, l;

	#ifdef DEBUG
		printf("ll1 = %d; lc1 = %d; ll2 = %d; lc2 = %d ------ MIN \n",ll1,lc1,ll2,lc2);
		printf("profondeur = %d  ---- MIN\n",profondeur);
	#endif

	Pion* copie[NB_COLONNES*NB_LIGNES];
	double val_min = FLT_MAX;
	if(profondeur == PROFONDEUR_MAX){
		return f_eval(plateau,joueur);
	}

	else{
		// printf("(%d)\nFMIN : [%f,%f]\n", nb_appels, *alpha, *beta);
		for(i = 0;i < NB_COLONNES;i++){
			for(j = 0;j < NB_LIGNES;j++){
				if(plateau[i*NB_COLONNES+j].couleur == joueur){
					for(k = j-1;k <= j+1; k++){
						for(l = i-1;l <= i+1;l++){
							if((!f_test_mouvement(plateau, i, j, l, k, joueur))&&((l != i)&&(k !=j))){
								f_copie_plateau(plateau,copie);
								f_bouge_piece(copie,i,j,l,k,joueur);
								nb_noeuds_visites ++;
								int coupure = f_max(copie, -joueur, profondeur+1, l1, c1, l2, c2, alpha, &lbeta);
								*beta = lbeta < *beta ? lbeta : *beta;

								#ifdef DEBUG
									printf("valeur beta %f ------ fonction MIN\n",*beta);
								#endif

								if(*beta >= *alpha){
									*alpha = *beta;
									return 1;
								}
								else{
									// val_min = val;
									ll1 = i;
									lc1 = j;
									ll2 = l;
									lc2 = k;
								}
							}
						}
					}
				}
			}
		}
		nb_noeuds_totals += nb_noeuds_visites;
		nb_noeuds_visites = 0;

		*l1 = ll1;
		*c1 = lc1;
		*l2 = ll2;
		*c2 = lc2;

		return 0;
	}
}

//Fonction max trouve le maximum des noeuds fils
double f_max(Pion *plateau, int joueur, int profondeur, int* l1, int* c1, int* l2, int* c2, double* alpha, double* beta)
{
	int ll1 = *l1;
	int lc1 = *c1;
	int ll2 = *l2;
	int lc2 = *c2;
	double lalpha = *alpha;
	int i, j, k, l;
	#ifdef DEBUG
		printf("ll1 = %d; lc1 = %d; ll2 = %d; lc2 = %d ------ MAX  \n",ll1,lc1,ll2,lc2);
		printf("profondeur = %d  ---- MAX\n",profondeur);
	#endif
	//printf("")
	Pion* copie[NB_COLONNES*NB_LIGNES];

	double val_max = -FLT_MAX+1;
	if(profondeur == PROFONDEUR_MAX){
		return f_eval(plateau,joueur);
	}

	else{
		// printf("(%d)\tFMAX : [%f,%f]\n", nb_appels, *alpha, *beta);
		for(i = 0;i < NB_COLONNES;i++){
			for(j = 0;j < NB_LIGNES;j++){
				if(plateau[i*NB_COLONNES+j].couleur == joueur){
					for(k = j-1;k <= j+1; k++){
						for(l = i-1;l <= i+1;l++){
							if((!f_test_mouvement(plateau, i, j, l, k, joueur))&&((l != i)&&(k !=j))){

								#ifdef DEBUG
									printf("i = %d;j = %d;k = %d;l = %d\n",i,j,k,l);
								#endif

								f_copie_plateau(plateau,copie);
								f_bouge_piece(copie,i,j,l,k,joueur);
								nb_noeuds_visites ++;
								int coupure = f_min(copie, -joueur, profondeur+1, l1, c1, l2, c2, &lalpha, beta);
								*alpha = lalpha > *alpha ? lalpha : *alpha;

								#ifdef DEBUG
									printf("valeur alpha %f ------ fonction MAX\n", *alpha);
								#endif

								if(*alpha >= *beta){
									*beta = *alpha;
									return 1;
								}
								else{
									ll1 = i;
									lc1 = j;
									ll2 = l;
									lc2 = k;
								}
							}
						}
					}
				}
			}
		}
		nb_noeuds_totals += nb_noeuds_visites;
		nb_noeuds_visites = 0;

		*l1 = ll1;
		*c1 = lc1;
		*l2 = ll2;
		*c2 = lc2;

		#ifdef DEBUG
			printf("RETURN %f\n",val_max);
		#endif
		return 0;
	}
}

/**
 * Calcule et joue le meilleur cout
 * */
void f_IA(int joueur,Pion* plateau)
{
	int l1 = 0;
	int c1 = 0;
	int l2 = 0;
	int c2 = 0;
	double alpha = -FLT_MAX+1;
	double beta = FLT_MAX;

	struct timeval  begin, end;


	nb_appels ++;
	gettimeofday(&begin, NULL);
	alpha = f_max(plateau,joueur,0,&l1,&c1,&l2,&c2, &alpha, &beta);
	gettimeofday(&end, NULL);

	double tmp = (double) (end.tv_usec - begin.tv_usec) / 1000000 + (double) (end.tv_sec - begin.tv_sec);
	printf("alpha/beta : [%f,%f]\n", alpha, beta);
	printf ("Temps de calcul: %f sec\n", tmp);
	total_time += tmp;
	f_bouge_piece(plateau,l1,c1,l2,c2,joueur);

#ifdef DEBUG
	printf("dbg: exiting %s %d\n", __FUNCTION__, __LINE__);
#endif
}


/**
 * Demande le choix du joueur humain et calcule le coup demande
 * */
void f_humain(int joueur)
{
	char c1, c2;
	char buffer[32];
	int l1, l2;



	printf("joueur ");
	switch(joueur)
	{
	case -1:
		printf("o ");
		break;
	case 1:
		printf("x ");
		break;
	default:
		printf("inconnu ");
	}
	printf("joue:\n");
	while(1)
	{
		fgets(buffer, 32, stdin);
		if(sscanf(buffer, "%c%i%c%i\n", &c1, &l1, &c2, &l2) == 4)
		{
			if(f_bouge_piece(plateauDeJeu, l1, f_convert_char2int(c1), l2, f_convert_char2int(c2), joueur) == 0)
				break;
		}
		fflush(stdin);
		printf("mauvais choix\n");
	}

}

int main(int argv, char *argc[])
{
	int fin = 0,mode=0 , ret, joueur = -1;

	plateauDeJeu = f_init_plateau();
	f_eval(plateauDeJeu,joueur);
	printf("1 humain vs IA\n2 humain vs humain\n3 IA vs IA\n");
	scanf("%d",&mode);
	while (!fin)
	{
		f_affiche_plateau(plateauDeJeu);
		if(mode==1)
		{
			f_affiche_plateau(plateauDeJeu);
			if(joueur>0)
				f_humain(joueur);
			else
				f_IA(joueur,plateauDeJeu);
		}
		else if(mode==2)
		{
			f_humain(joueur);
		}
		else
		{
			f_IA(joueur,plateauDeJeu);
		}

		if ((ret = f_gagnant()) != 0)
		{
			switch (ret)
			{
			case 1:
				f_affiche_plateau(plateauDeJeu);
				printf("\njoueur x gagne!\n");
				fin = 1;
				break;
			case -1:
				f_affiche_plateau(plateauDeJeu);
				printf("\njoueur o gagne!\n");
				fin = 1;
				break;
			}
		}
		joueur = -joueur;
	}

	printf("\n###########################################################\n");
	printf("\t\t\tSTATISTIQUES\n");
	printf("Profondeur : %d\n", PROFONDEUR_MAX);
		printf("Nb coups : %d\n", nb_appels);
	printf("Temps total d'execution : %f sec\t Moyenne de %f ms par tour\n", total_time,
																																				total_time/nb_appels*1000);
	printf("%d noeuds explores au total\t Moyenne de %d par tour\n", nb_noeuds_totals,
																														 			 nb_noeuds_totals/nb_appels);
	printf("###########################################################\n");

	return 0;
}
