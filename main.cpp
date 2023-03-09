#include <cstdio>
#include <cmath>
#include <malloc.h>
#include <ctime>


#define  ll long long
const double PI = 3.141592653589793;
static int prec = 10000;
double ***table;
ll numSpacePointX, numTimePoint, numSpacePointY;

struct TimeParam {
    double t_in, t_sur, delta_t, t_lim;
};

struct SpaceParam {
    double length, diffusion, delta_x, width, delta_y;
};

void calNumTime(struct TimeParam t) {
    numTimePoint = (ll) (t.t_lim / t.delta_t);
}


void calNumSpace(struct SpaceParam s) {
    numSpacePointX = (ll) (s.length / s.delta_x);
    numSpacePointY = (ll) (s.width / s.delta_y);
}

double ***createTable(ll x, ll y, ll t) {
    double ***arr = (double ***) malloc(t * sizeof(double **));
    for (ll i = 0; i < t; ++i) {
        arr[i] = (double **) malloc(x * sizeof(double *));
        for (int j = 0; j < x; ++j) {
            arr[i][j] = (double *) malloc(y * sizeof(double));
        }
    }
    return arr;
}

double getValue(struct TimeParam time, struct SpaceParam space, double x, double t, double y, ll precision) {
    double sum = 0.0, exponential, spaceXTerm, spaceYTerm;
    for (ll m = 1; m < precision; ++m) {
        exponential = exp(-space.diffusion * t *
                          (pow(PI * (double) m / space.length, 2) + pow(PI * (double) m / space.width, 2)));
        spaceXTerm = (1 - pow(-1, (double) m) / ((double) m * PI)) * sin((double) m * PI * x / space.length);
        spaceYTerm = (1 - pow(-1, m) / ((double) m * PI)) * sin((double) m * PI * y / space.width);
        sum += exponential * spaceXTerm * spaceYTerm;
    }
    return time.t_sur + 2 * (time.t_in - time.t_sur) * sum;
}

double getNextStep(struct TimeParam time, struct SpaceParam space, ll x, ll t, ll y) {
    double newY = (double) y * space.delta_y;
    double newX = (double) x * space.delta_x;
    double newT = (double) t * time.delta_t;
    return getValue(time, space, newX, newT, newY, prec);
}

void solve(struct SpaceParam space, struct TimeParam time) {
    for (ll t = 0; t < numTimePoint; ++t) {
        for (ll x = 0; x < numSpacePointX; ++x) {
            for (int y = 0; y < numSpacePointY; ++y) {
                table[t][x][y] = getNextStep(time, space, x, t, y);
            }
        }
    }
}


int main() {
    clock_t t;
    t = clock();
    FILE *fptr;
    fptr = fopen("temps.txt", "w");
    struct TimeParam time = {100, 300, 0.01, 1};
    struct SpaceParam space = {1, 0.1, 0.05, 1, 0.05};
    calNumSpace(space);
    calNumTime(time);
    table = createTable(numSpacePointX, numSpacePointY, numTimePoint);
    solve(space, time);

    for (ll i = 0; i < numTimePoint; ++i) {
        for (ll j = 0; j < numSpacePointX; ++j) {
            for (int k = 0; k < numSpacePointY; ++k) {
                fprintf(fptr, "%f ", table[i][j][k]);
            }
            fprintf(fptr, "\n");
        }
        fprintf(fptr, "\n\n");
    }
    printf("%lld %lld %lld",numTimePoint,numSpacePointX,numSpacePointY);
    fclose(fptr);
    for (ll i = 0; i < numTimePoint; ++i) {
        for (ll j = 0; j < numSpacePointX; ++j) {
            free(table[i][j]);
        }
        free(table[i]);
    }
    t = clock() - t;
    double time_taken = ((double) t) / CLOCKS_PER_SEC; // in seconds
    printf("took %f seconds to execute \n", time_taken);

    return 0;
}

