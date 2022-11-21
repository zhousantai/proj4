#include <stdio.h>

#include "matrix.c"
#include "cblas.h"

int main() {
    matrix *a = malloc(sizeof(matrix));
    matrix *b = malloc(sizeof(matrix));
    matrix *c = malloc(sizeof(matrix));
    matrix *d = malloc(sizeof(matrix));
    for (int i = 1; i <10 ; ++i) {
        printf("--------- it is %d th time\n",i);

        createMatrix(a, i*1000, i*1000);
        createMatrix(b, i*1000, i*1000);
        createMatrix(c, i*1000, i*1000);
        createMatrix(d, i*1000, i*1000);

        for (int i = 0; i < a->row*a->col; ++i) {
            a->data[i]=rand()%10;
        }

        TIME_START
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, i*1000, i*1000, i*1000, 1.0f, a->data, i*1000, b->data, i*1000, 0,
                    c->data, i*1000);
        TIME_END("blas")

        //  c= mulMatirx_plain(a,b);

        d = mulMatirx_improved(a, b);

        bool testForA=true;
        for (int i = 0; i < a->row*b->col; ++i) {
            if (c->data[i]!=d->data[i])testForA=false;
        }
        printf("%d\n",testForA);
//    clock_t start1,finish1;
//    double duration1;
//    start1=clock();

        d=mulMatrix_devided(a, b);
        for (size_t  i = 0; i < a->row*b->col; ++i) {
            if (c->data[i]!=d->data[i]){ printf("i: %d\n",i);testForA=false;
                break;}
        }
        printf("%d\n",testForA);

    }


    return 0;
}
