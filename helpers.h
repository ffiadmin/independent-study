/*
 * Chris Prosser
 * 
 * These functions provide a nice wrapper
 * to Dr. Falcetta's procedural-oriented
 * functions.  
 *
 * Future goal: 
 * Make these functions work indepentently 
 * in an object-oriented way using Psi4
 */

#include "functions.h"
#include "InputParser.h"

double calculateOverlap(Basis* first, Basis* second){
	return overlap(	first->alpha, first->l, first->m, first->n, first->x, first->y, first->z,
					second->alpha, second->l, second->m, second->n, second->x, second->y, second->z);
}

double calculateKineticEnergy(Basis* first, Basis* second){
	return kinetic(	first->alpha, first->l, first->m, first->n, first->x, first->y, first->z,
					second->alpha, second->l, second->m, second->n, second->x, second->y, second->z);
}

Complex calculateAttraction(Basis* first, Basis* second, Nuclei* nuc){
	return nuclear_attraction(first->x, first->y, first->z, first->norm, first->l, first->m, first->n, first->alpha,
						second->x, second->y, second->z, second->norm, second->l, second->m, second->n, second->alpha,
						nuc->nx, nuc->ny, nuc->nx);
}

