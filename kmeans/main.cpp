#include <string>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <sys/time.h>

#include <omp.h>

#define DEFAULT_K 5
#define DEFAULT_NITERS 20

struct point_t {
    double x;
    double y;
};

double get_time() {
    struct timeval tv;
    gettimeofday(&tv, (struct timezone *) 0);
    return ((double) tv.tv_sec + (double) tv.tv_usec / 1000000.0);
}

void read_points(std::string filename, point_t *p, int n) {
    std::ifstream infile{filename};
    double x, y;
    int i = 0;
    while (infile >> x >> y) {
        if (i >= n) {
            printf("WARNING: more points in input file '%s' than read: stopping after %d lines\n", filename.c_str(), i);
            return;
        }
        p[i].x = x;
        p[i].y = y;
        i++;
    }
}

void write_memory(std::string filename, int niters, point_t *memory, int k) {
    std::ofstream outfile{filename};
    for (int iter = 0; iter < niters + 1; ++iter) {
        for (int i = 0; i < k; ++i) {
            outfile << iter << ' ' << memory[iter * k + i].x << ' ' << memory[iter * k + i].y << '\n';
        }
    }
}

void init_centroids(point_t *centroids, int k) {
    for (int i = 0; i < k; ++i) {
        centroids[i].x = rand() % 100;
        centroids[i].y = rand() % 100;
    }
}

void k_means(uint8_t niters, point_t *points, 
            point_t *centroids, uint32_t *assignment, 
            point_t *memory, uint32_t n, 
            uint16_t k, uint32_t* count, 
            double *sum_x, double *sum_y) {

    for (uint8_t iter = 0; iter < niters; ++iter) {
        // determine nearest centroids
        #pragma omp target teams distribute parallel for
        for (uint32_t i = 0; i < n; ++i) {
            double optimal_dist = DBL_MAX;
            for (uint16_t j = 0; j < k; ++j) {

                // using normal arithmetic operators without square root (logic is still correct)
                // to reduce calculation workload
                double dist = (points[i].x - centroids[j].x)*(points[i].x - centroids[j].x) + 
                              (points[i].y - centroids[j].y)*(points[i].y - centroids[j].y);
                if (dist < optimal_dist) {
                    optimal_dist = dist;
                    assignment[i] = j;
                }
            }
        }

        #pragma omp target update from(assignment[0:n])
        for (uint16_t j = 0; j < k; ++j) {
            count[j] = 0;
            sum_x[j] = 0;
            sum_y[j] = 0;
        }

        for (uint32_t i = 0; i < n; ++i) {
            count[assignment[i]]++;
            sum_x[assignment[i]]+= points[i].x;
            sum_y[assignment[i]] += points[i].y;
        }

        // update centroid positions
        for (uint16_t j = 0; j < k; ++j) {
            if (count[j] != 0) {
                centroids[j].x = sum_x[j] / count[j];
                centroids[j].y = sum_y[j] / count[j];
            }
            memory[(iter + 1) * k + j].x = centroids[j].x;
            memory[(iter + 1) * k + j].y = centroids[j].y;
        }

        #pragma omp target update to(centroids[0:k])
    }

}

int main(int argc, const char *argv[]) {
    srand(1234);
    if (argc < 3 || argc > 5) {
        printf("Usage: %s <input file> <num points> <num centroids> <num iters>\n", argv[0]);
        return EXIT_FAILURE;
    }
    const char *input_file = argv[1];
    const uint32_t n = atoi(argv[2]);
    const uint16_t k = (argc > 3 ? atoi(argv[3]) : DEFAULT_K);
    const uint8_t niters = (argc > 4 ? atoi(argv[4]) : DEFAULT_NITERS);

    point_t *points = (point_t *) malloc(n * sizeof(point_t));
    point_t *centroids = (point_t *) malloc(k * sizeof(point_t));
    uint32_t *assignment = (uint32_t *) malloc(n * sizeof(int));
    // reserve extra space to save the initial centroid placement
    point_t *memory = (point_t *) malloc((niters + 1) * k * sizeof(point_t));

    uint32_t *count = (uint32_t *) malloc(k * sizeof(point_t));
    double *sum_x = (double *) malloc(k * sizeof(point_t));
    double *sum_y = (double *) malloc(k * sizeof(point_t));

    read_points(input_file, points, n);
    init_centroids(centroids, k);

    for (uint16_t j = 0; j < k; ++j) {
        memory[j].x = centroids[j].x;
        memory[j].y = centroids[j].y;
    }


    #pragma omp target teams distribute parallel for map(tofrom: count[0:k], sum_x[0:k], sum_y[0:k])
    for (uint16_t j = 0; j < k; ++j) {
        count[j] = 0;
        sum_x[j] = 0;
        sum_y[j] = 0;
    }

    printf("Executing k-means ?|   %d iterations with %d points and %d centroids...\n", niters, n, k);
    double runtime = get_time();
    #pragma omp target enter data map(to:points[0:n], centroids[0:k], assignment[0:n], count[0:k], sum_x[0:k], sum_y[0:k])
    k_means(niters, points, centroids, assignment, memory, n, k,count, sum_x, sum_y);
    runtime = get_time() - runtime;

    printf("Time Elapsed: %f s\n", runtime);

    write_memory("memory.out", niters, memory, k);

    free(assignment);
    free(centroids);
    free(points);
    free(memory);

    return EXIT_SUCCESS;
}