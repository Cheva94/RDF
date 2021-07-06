#include <stdio.h>
#include <stdlib.h>
#include <math.h>

# define M_PI
# define BOX

void init(void){
	ngr = 0;
	delg = box / (2*nhis);	// bin size
	for (int i=1; nhis; i++){ // nhis total number of bins
		g(i) = 0;
	}
}

void sample(void){	// podría llamarse histograma
	ngr += 1;
	for (int i=1; N-1; i++){
		for (int j=i+1; N; j++){
			xr = x(i) - x(j);
			xr -= box * round (xr/box);	// PBC
			r = sqrt(xr*xr);
			if (r<box/2){	// sólo dentro de media caja
				ig = int (r/delg);
				g(ig) += 2;	// contribución de las partículas
}

void results(void){	// podría llamarse rdf o algo así ya que acá es donde se hace
	for (int i=1; nhis; i++){
		r = delg*(i+0.5);	// distancia r
		vb = (pow(i+1, 3)-pow(i,3))*pow(delg,3);	// volumen entre bin i y i+1
		nid = (4/3) * M_PI * vb * rho;	// cantidad de partículas ideales en vb
		g(i) /= ngr*N*nid;	// normalizacion de la rdf
}

int main(void) {
    int N; // cantidad de partículas
    double xr;

    double ngr, box, r, ig, nhis, delg, vb, nid, rho, g;

	init();
	sample();
	results();

	return 0
}
