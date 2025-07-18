/* laplace_solver.c
 *
 * Reads:
 *   Nx Ny Nz                   (ints)
 *   hx hy hz                   (doubles)
 *   type_grid[Nx*Ny*Nz]        (0 interior, 1 Dirichlet, 2 Neumann)
 *   dirichlet_grid[Nx*Ny*Nz]   (values for Dirichlet)
 *
 * Writes to stdout:
 *   phi[Nx*Ny*Nz]              (computed potentials, row-major)
 *
 * Uses Gauss–Seidel until max change < tol or max_iter reached.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main() {
    int Nx, Ny, Nz;
    if (scanf("%d %d %d", &Nx, &Ny, &Nz) != 3) {
        fprintf(stderr, "Error reading grid dimensions.\n");
        return 1;
    }

    double hx, hy, hz;
    if (scanf("%lf %lf %lf", &hx, &hy, &hz) != 3) {
        fprintf(stderr, "Error reading grid spacings.\n");
        return 1;
    }

    long N = (long)Nx * Ny * Nz;
    int  *type = malloc(N * sizeof(int));
    double *dirich = malloc(N * sizeof(double));
    double *phi = malloc(N * sizeof(double));

    if (!type || !dirich || !phi) {
        fprintf(stderr, "Allocation failed.\n");
        return 1;
    }

    for (long i = 0; i < N; i++) {
        if (scanf("%d", &type[i]) != 1) {
            fprintf(stderr, "Error reading type_grid[%ld]\n", i);
            return 1;
        }
    }
    for (long i = 0; i < N; i++) {
        if (scanf("%lf", &dirich[i]) != 1) {
            fprintf(stderr, "Error reading dirichlet_grid[%ld]\n", i);
            return 1;
        }
    }

    /* Initialize phi: Dirichlet points to given value, others to 0 */
    for (long i = 0; i < N; i++) {
        phi[i] = (type[i] == 1 ? dirich[i] : 0.0);
    }

    /* Precompute coefficient weights for non-uniform grid spacing */
    double idx2 = 1.0/(hx*hx), idy2 = 1.0/(hy*hy), idz2 = 1.0/(hz*hz);
    double diag = 2*(idx2 + idy2 + idz2);

    const int max_iter = 20000;
    const double tol = 1e-6;
    double maxdiff;

    for (int iter = 0; iter < max_iter; iter++) {
        maxdiff = 0.0;

        for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nz; k++) {
            long idx = (long)i*Ny*Nz + j*Nz + k;

            if (type[idx] == 0) {
                /* interior: update via 7‑point FD */
                double sum = 0.0;
                double w = 0.0;

                /* +x neighbor */
                if (i+1 < Nx) {
                    long idxp = idx + Ny*Nz;
                    sum += idx2 * phi[idxp];
                    w   += idx2;
                } else {
                    /* Neumann: zero‐flux => use self */
                    sum += idx2 * phi[idx];
                    w   += idx2;
                }
                /* -x neighbor */
                if (i-1 >= 0) {
                    long idxm = idx - Ny*Nz;
                    sum += idx2 * phi[idxm];
                    w   += idx2;
                } else {
                    sum += idx2 * phi[idx];
                    w   += idx2;
                }
                /* +y */
                if (j+1 < Ny) {
                    long idyp = idx + Nz;
                    sum += idy2 * phi[idyp];
                    w   += idy2;
                } else {
                    sum += idy2 * phi[idx];
                    w   += idy2;
                }
                /* -y */
                if (j-1 >= 0) {
                    long idym = idx - Nz;
                    sum += idy2 * phi[idym];
                    w   += idy2;
                } else {
                    sum += idy2 * phi[idx];
                    w   += idy2;
                }
                /* +z */
                if (k+1 < Nz) {
                    long idzp = idx + 1;
                    sum += idz2 * phi[idzp];
                    w   += idz2;
                } else {
                    sum += idz2 * phi[idx];
                    w   += idz2;
                }
                /* -z */
                if (k-1 >= 0) {
                    long idzm = idx - 1;
                    sum += idz2 * phi[idzm];
                    w   += idz2;
                } else {
                    sum += idz2 * phi[idx];
                    w   += idz2;
                }

                double newphi = sum / w;
                double diff = fabs(newphi - phi[idx]);
                if (diff > maxdiff) maxdiff = diff;
                phi[idx] = newphi;
            }
            /* type==1 (Dirichlet): phi already set, skip */
            /* type==2 (Neumann at boundary): treat same as interior but no fixed value */
        }}}

        if (maxdiff < tol) {
            fprintf(stderr, "Converged in %d iters (Δ=%g)\n", iter, maxdiff);
            break;
        }
        if (iter == max_iter-1) {
            fprintf(stderr, "Warning: max_iter reached, Δ=%g\n", maxdiff);
        }
    }

    /* Output phi flat */
    /*
    for (long i = 0; i < N; i++) {
        printf("%.10e\n", phi[i]);
    }*/

    // Save solution
    FILE *fphi = fopen("solution.txt", "w");
    if (!fphi) {
        fprintf(stderr, "Error: could not open solution.txt for writing.\n");
        return 1;
    }
    for (long i = 0; i < N; i++) {
        fprintf(fphi, "%.10e\n", phi[i]);
    }
    fclose(fphi);

    // Save dimensions
    FILE *fdim = fopen("grid_Nx.txt", "w");
    if (!fdim) { fprintf(stderr, "Error writing Nx\n"); return 1; }
    fprintf(fdim, "%d\n", Nx); fclose(fdim);
    fprintf(stderr, "gridNx saved! (%d) \n", Nx);

    //print so that the output is received in the python code
    for (int i = 0; i < N; i++) {
        fprintf(stderr, "%f ", phi[i]);
        printf("%f ", phi[i]);
    }
    fprintf(stderr, "\n");


    fdim = fopen("grid_Ny.txt", "w");
    if (!fdim) { fprintf(stderr, "Error writing Ny\n"); return 1; }
    fprintf(fdim, "%d\n", Ny); fclose(fdim);
    fprintf(stderr, "gridNy saved! \n");

    fdim = fopen("grid_Nz.txt", "w");
    if (!fdim) { fprintf(stderr, "Error writing Nz\n"); return 1; }
    fprintf(fdim, "%d\n", Nz); fclose(fdim);
    fprintf(stderr, "gridNz saved! \n");

    free(type);
    free(dirich);
    free(phi);
    return 0;
}
