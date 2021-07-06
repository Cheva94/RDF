#include <omp.h> // omp_get_wtime()
#include <stdio.h> // printf()
#include <math.h>

#define N (1<<26)

struct point {
	float dx, dy, dz, dw;
} aos[N];

struct soa {
	float dx[N], dy[N], dz[N], dw[N];
} soa;

float max_dist_aos(void);
float max_dist_soa(void);

int main(int argc, char ** argv)
{
	float max_dist_0 = 0.0f, max_dist_1 = 0.0f;

	max_dist_0 += max_dist_aos(); // warm up
	double t = omp_get_wtime();
	max_dist_0 *= max_dist_aos();
	printf("max dist aos: %lf\n", omp_get_wtime()-t);

	max_dist_1 += max_dist_soa(); // warm up
	t = omp_get_wtime();
	max_dist_1 *= max_dist_soa();
	printf("max dist soa: %lf\n", omp_get_wtime()-t);

	return (int)(max_dist_0+max_dist_1);
}


float max_dist_aos(void) {
	float max_dist = 0.0f;
	for (size_t i=0; i<N; ++i) {
		float dist = sqrtf(aos[i].dx*aos[i].dx + aos[i].dy*aos[i].dy + aos[i].dz*aos[i].dz + aos[i].dw*aos[i].dw);
		if (max_dist<dist)
			max_dist = dist;
	}

	return max_dist;
}


float max_dist_soa(void) {
	float max_dist = 0.0f;
	for (size_t i=0; i<N; ++i) {
		float dist = sqrtf(soa.dx[i]*soa.dx[i] + soa.dy[i]*soa.dy[i] + soa.dz[i]*soa.dz[i] + soa.dw[i]*soa.dw[i]);
		if (max_dist<dist)
			max_dist = dist;
	}

	return max_dist;
}
