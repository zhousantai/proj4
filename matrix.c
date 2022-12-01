
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <stdbool.h>
#include <string.h>
#include "matrix.h"
#include "omp.h"
#include "time.h"
#include <sys/time.h>


#include <arm_neon.h>


#define US_PER_SEC 1000000.f
#define TIME_START   struct timeval start3,end3;gettimeofday(&start3, NULL);
#define TIME_END(zst)  gettimeofday(&end3, NULL);printf("Time spent: %.5fs in %s mul\n",\
    end3.tv_sec - start3.tv_sec +(end3.tv_usec - start3.tv_usec) / US_PER_SEC,zst);

matrix *mulmatrix_plain(const matrix *A, const matrix *B) {
   // CHECK((&A),(&B))
    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);
//    clock_t start,finish;
//    double duration;
//    start=clock();
    TIME_START

    for (size_t i = 0; i < result->row; i++) {
        for (size_t j = 0; j < result->col; j++) {
            for (size_t k = 0; k < A->col; k++) {
                result->data[i * result->col + j] += A->data[i * A->col + k] * B->data[k * B->col + j];
            }
        }
    }
    TIME_END("plain")
//    finish=clock();
//    duration=(double)(finish-start)/CLOCKS_PER_SEC;
//    printf("duration for plain: %f\n",duration);
    return result;

}
matrix *mulmatrix_openMP(const matrix *A, const matrix *B) {

    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);
//    clock_t start,finish;
//    double duration;
//    start=clock();
    TIME_START
#pragma omp parallel for
    for (size_t i = 0; i < result->row; i++) {
        for (size_t j = 0; j < result->col; j++) {
            for (size_t k = 0; k < A->col; k++) {
                result->data[i * result->col + j] += A->data[i * A->col + k] * B->data[k * B->col + j];
            }
        }
    }
    TIME_END("openMP")
//    finish=clock();
//    duration=(double)(finish-start)/CLOCKS_PER_SEC;
//    printf("duration for plain: %f\n",duration);
    return result;

}
matrix *mulmatrix_openMP_ikj(const matrix *A, const matrix *B) {

    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);
    TIME_START
#pragma omp parallel for
    for (size_t i = 0; i < result->row; i++) {
        for (size_t k = 0; k < A->col; k++) {
            for (size_t j = 0; j < result->col; j++) {
                result->data[i * result->col + j] += A->data[i * A->col + k] * B->data[k * B->col + j];
            }
        }
    }
    TIME_END("openMP+ikj")
    return result;

}
matrix *mulmatrix_openMP_ikj_register(const matrix *A, const matrix *B) {

    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);

    TIME_START
#pragma omp parallel for
    for (register_t i = 0; i < result->row; i++) {
        for (register_t k = 0; k < A->col; k++) {
            for (register_t j = 0; j < result->col; j++) {
                result->data[i * result->col + j] += A->data[i * A->col + k] * B->data[k * B->col + j];
            }
        }
    }
    TIME_END("openMP+ikj+register")
    return result;

}
matrix *mulmatrix_ikj(const matrix *A, const matrix *B) {

    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);
//    clock_t start,finish;
//    double duration;
//    start=clock();
    TIME_START
    for (size_t i = 0; i < result->row; i++) {
        for (size_t k = 0; k < A->col; k++) {
            for (size_t j = 0; j < result->col; j++) {
                result->data[i * result->col + j] += A->data[i * A->col + k] * B->data[k * B->col + j];
            }
        }
    }
    TIME_END("ikj")
//    finish=clock();
//    duration=(double)(finish-start)/CLOCKS_PER_SEC;     clock（）失败的trial
//    printf("duration for plain: %f\n",duration);
    return result;

}

matrix *mulMatrix_devided_improved(const matrix*A,const matrix*B)
{TIME_START
    matrix* A1=(matrix*) malloc(sizeof(matrix));
    matrix* A2=(matrix*) malloc(sizeof(matrix));
    matrix* A3=(matrix*) malloc(sizeof(matrix));
    matrix* A4=(matrix*) malloc(sizeof(matrix));

    division(A,2,A1,A2,A3,A4);
    matrix* B1=(matrix*) malloc(sizeof(matrix));
    matrix* B2=(matrix*) malloc(sizeof(matrix));
    matrix* B3=(matrix*) malloc(sizeof(matrix));
    matrix* B4=(matrix*) malloc(sizeof(matrix));
    division(B,2,B1,B2,B3,B4);
    matrix * ans1 = addMatrix(mulmatrix_improved(A1,B1), mulmatrix_improved(A2,B3));
    matrix * ans2 = addMatrix(mulmatrix_improved(A1,B2), mulmatrix_improved(A2,B4));
    matrix * ans3 = addMatrix(mulmatrix_improved(A3,B1), mulmatrix_improved(A4,B3));
    matrix * ans4 = addMatrix(mulmatrix_improved(A3,B2), mulmatrix_improved(A4,B4));
    matrix* ans = (matrix*)(malloc(sizeof(matrix)));
    createMatrix(ans,A->row,A->col);
    merge(ans,ans1,ans2,ans3,ans4);
    free(ans1);free(A1);free(A2);free(A3);free(A4);free(ans2);free(ans3);free(ans4);free(B1);free(B2);free(B3);free(B4);
    TIME_END("devided_improved")
    return ans;
}matrix *mulMatrix_devided_plain(const matrix*A,const matrix*B)
{TIME_START
    matrix* A1=(matrix*) malloc(sizeof(matrix));
    matrix* A2=(matrix*) malloc(sizeof(matrix));
    matrix* A3=(matrix*) malloc(sizeof(matrix));
    matrix* A4=(matrix*) malloc(sizeof(matrix));

    division(A,2,A1,A2,A3,A4);
    matrix* B1=(matrix*) malloc(sizeof(matrix));
    matrix* B2=(matrix*) malloc(sizeof(matrix));
    matrix* B3=(matrix*) malloc(sizeof(matrix));
    matrix* B4=(matrix*) malloc(sizeof(matrix));
    division(B,2,B1,B2,B3,B4);
    matrix * ans1 = addMatrix(mulmatrix_plain(A1,B1), mulmatrix_plain(A2,B3));
    matrix * ans2 = addMatrix(mulmatrix_plain(A1,B2), mulmatrix_plain(A2,B4));
    matrix * ans3 = addMatrix(mulmatrix_plain(A3,B1), mulmatrix_plain(A4,B3));
    matrix * ans4 = addMatrix(mulmatrix_plain(A3,B2), mulmatrix_plain(A4,B4));
    matrix* ans = (matrix*)(malloc(sizeof(matrix)));
    createMatrix(ans,A->row,A->col);
    merge(ans,ans1,ans2,ans3,ans4);
    free(ans1);free(A1);free(A2);free(A3);free(A4);free(ans2);free(ans3);free(ans4);free(B1);free(B2);free(B3);free(B4);
    TIME_END("devided_plain")
    return ans;
}
void division(const matrix*A,int t,matrix* A1,matrix *A2,matrix*A3,matrix*A4)
{

TIME_START
    A1->col=A->col/2;
    A1->row=A->row/2;
    free(A1->data);
    A1->data=(float *) malloc(4*A1->row*A1->col);

#pragma omp parallel for
    for (register_t i = 0; i < A1->row; ++i) {
        memcpy(A1->data+(i*A1->col),A->data+(i*A->col),A1->col*sizeof(float ));
//        for (register_t j = 0; j <A1->col ; ++j) {
//            A1->data[i*A1->col+j]=A->data[i*A->col+j];
//            printf("A1->data[i*A1->col+j]:%f %ld\n",A1->data[i*A1->col+j],i*A1->col+j);
//            printf("A->data[i*A->col+j]:%f %ld\n",A->data[i*A->col+j],i*A->col+j);
//            printf("j:%ld\n",j);
//        }
    }

    A2->col=A->col/2;
    A2->row=A->row/2;
    free(A2->data);
    A2->data=(float *) malloc(4*A2->row*A2->col);
#pragma omp parallel for
    for (register_t  i = 0; i < A2->row; ++i) {
        memcpy(A2->data+(i*A2->col),A->data+(i*A->col+A->col/2),A2->col*sizeof(float ));
//        for (register_t  j = 0; j <A2->col ; ++j) {
//            A2->data[i*A2->col+j]=A->data[i*A->col+A->col/2+j];
//        }
    }
    A3->col=A->col/2;
    A3->row=A->row/2;
    free(A3->data);
    A3->data=(float *) malloc(4*A3->row*A3->col);
#pragma omp parallel for
    for (register_t  i = 0; i < A3->row; ++i) {
        memcpy(A3->data+(i*A3->col),A->data+(i+A->row/2)*A->col,A3->col*sizeof(float ));
//        for (register_t j = 0; j <A3->col ; ++j) {
//            A3->data[i*A3->col+j]=A->data[(i+A->row/2)*A->col+j];
//        }
    }
    A4->col=A->col/2;
    A4->row=A->row/2;
    free(A4->data);
    A4->data=(float *) malloc(4*A3->row*A3->col);
#pragma omp parallel for
    for (register_t i = 0; i < A4->row; ++i) {
        memcpy(A4->data+(i*A4->col),A->data+((i+A->row/2)*A->col+A->col/2),A4->col*sizeof(float ));
//        for (register_t  j = 0; j <A4->col ; ++j) {
//            A4->data[i*A4->col+j]=A->data[(i+A->row/2)*A->col+A->col/2+j];
//        }
    }
    TIME_END("for copy")
}
void merge(matrix*ans,matrix* ans1,matrix* ans2,matrix* ans3,matrix* ans4){
    TIME_START
#pragma omp parallel for
    for (register_t i = 0; i < ans1->row; ++i) {
        memcpy(ans->data+(i*ans->col),ans1->data+(i*ans1->col),ans->col*sizeof(float ));
//        for (register_t j = 0; j <ans1->col ; ++j) {
//            ans->data[i*ans->col+j]=ans1->data[i*ans1->col+j];
//        }
    }
#pragma omp parallel for
    for (register_t i = 0; i < ans2->row; ++i) {
        memcpy(ans->data+(i*ans->col),ans2->data+(i*ans2->col),ans->col*sizeof(float ));
//        for (register_t j = 0; j <ans2->col ; ++j) {
//            ans->data[i*ans->col+ans->col/2+j]=ans2->data[i*ans2->col+j];
//        }
    }
#pragma omp parallel for
    for (register_t i = 0; i < ans3->row; ++i) {
        memcpy(ans->data+(i*ans->col),ans3->data+(i*ans3->col),ans->col*sizeof(float ));
//        for (register_t j = 0; j <ans3->col ; ++j) {
//            ans->data[(i+ans->row/2)*ans->col+j]=ans3->data[i*ans3->col+j];
//        }
    }
#pragma omp parallel for
    for (register_t i = 0; i < ans4->row; ++i) {
        memcpy(ans->data+(i*ans->col),ans4->data+(i*ans4->col),ans->col*sizeof(float ));
//        for (register_t j = 0; j <ans4->col ; ++j) {
//            ans->data[(i+ans->row/2)*ans->col+ans->col/2+j] =ans4->data[i*ans4->col+j] ;
//        }
    }
    TIME_END("for merge")
}
matrix *mulmatrix_improved(const matrix *A, const matrix *B) {
//struct timeval start3,end3;

    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);
//    clock_t start1, finish1;
//    double duration1;
//    start1 = clock();
    transpose(B);
    int result_row = result->row;
    int result_col = result->col;
    int A_col = A->col;
    int B_col = B->col;
    size_t a = result->row;
    size_t b = result->col;
    size_t c = A->col;
    TIME_START

#pragma omp parallel for
    for (register_t i = 0; i < a; i++) {
        for (register_t j = 0; j < b; j++) {
          //  for (register_t k = 0; k < c; k += 4)
            {
                //float sum1[4] = {0};

                float* sum1=(float *) malloc(16);
                float32x4_t a1, b1;
                float32x4_t c1 = vdupq_n_f32(0);
                size_t l;
                for ( l= 0; l < c; l += 4) {
                    a1 = vld1q_f32(A->data + (i * A_col + l));
                    b1 = vld1q_f32(B->data + (j * B_col + l));
                    c1 = vaddq_f32(c1, vmulq_f32(a1, b1));
                }
                if (l-c>0){
                    for (int k = 0; k < l-c; ++k) {
                        result->data[i * result_col + j]+=A->data[j*A->col+l-4+k]*B->data[j*B->col+l-4+k];
                    }
                }//处理不是4的倍数的情况

                vst1q_f32(sum1, c1);
                result->data[i * result_col + j] += sum1[0] + sum1[1] + sum1[2] + sum1[3];
//                register_t temp = c - k;
//                if (temp >= 8) {
//
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 0)] * B->data[(k + 0) * B_col + j];
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 1)] * B->data[(k + 1) * B_col + j];
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 2)] * B->data[(k + 2) * B_col + j];
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 3)] * B->data[(k + 3) * B_col + j];
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 4)] * B->data[(k + 4) * B_col + j];
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 5)] * B->data[(k + 5) * B_col + j];
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 6)] * B->data[(k + 6) * B_col + j];
//                    result->data[i * result_col + j] += A->data[i * A_col + (k + 7)] * B->data[(k + 7) * B_col + j];
//                } else {
//                    for (int l = 0; l < temp; ++l) {
//                        result->data[i * result_col + j] += A->data[i * A_col + (k + l)] * B->data[(k + l) * B_col + j];
//                    }
//                }
                free(sum1);


            }
        }
    }
    TIME_END("improved")
//    finish1=clock();
//    duration1=(double)(finish1-start1)/CLOCKS_PER_SEC;
//    int k = omp_get_num_procs();
//    printf("%d\n",k);
//    duration1/=k;
//    printf("duration for improved: %f\n",duration1);

    transpose(B);//做完运算，transpose回来
    return result;
}



/////////////////////////////////后续为与proj3中重复的代码
/////////////////////////////////后续为与proj3中重复的代码
/////////////////////////////////后续为与proj3中重复的代码
bool createMatrix(matrix *A, size_t row, size_t col) {


    if (row <= 0 || col <= 0 || A == NULL)return false;
    A->row = row;
    A->col = col;
    A->data = (float *) aligned_alloc(128,sizeof(float) * A->row * A->col);
//    for (int i = 0; i < A->row; i++) {
//        A->data[i] = (float *) malloc(sizeof(float) * A->col);
//    }
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->data[i * A->col + j] = 0;
        }
    }
    return true;
}

bool createMatrixN(matrix *A, size_t row, size_t  col, float q) {//给用户设置一个所有元素均为same certain number的构造
//    Node *newNode;
//    newNode = (Node *) malloc(sizeof(Node));
//    newNode->next = head.next;
//    head.next = newNode;
//    newNode->val = (size_t) A;

    if (row <= 0 || col <= 0 || A == NULL)return false;
    A->row = row;
    A->col = col;
    A->data = (float *)  aligned_alloc(128,sizeof(float) * A->row * A->col);
//    for (int i = 0; i < A->row; i++) {
//        A->data[i] = (float *) malloc(sizeof(float) * A->col);
//    }
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->data[i * A->col + j] = q;
        }
    }
    return true;
}

bool deleteMatrix(matrix *A) {
    if (A==NULL)return false;
    else{
        free(A->data);
        free(A);
        return true;
    }
}


bool initMatrix(matrix *A) {
    // 初始化矩阵
    if (!check(A))return false;
    printf("please input element for matrix with %ld row and %ld column:\n", A->row, A->col);
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {

            char *str = (char *) malloc(sizeof(char) * 100); //创建长度为100的字符数组
            //读入6个字符，stdin声明在stdio.h中，表示标准输入，从键盘输入
            scanf("%s", str);
            char **endptr = NULL;
            float ans = strtod(str, endptr);
            if (endptr != NULL)return false;//假如输入的非法，则返回flase！
            A->data[i * A->col + j] = ans;
            free(str);
            free(endptr);
        }
    }
    return true;
}

matrix *addMatrix(const matrix *A, const matrix *B) {
    // 判断矩阵A和矩阵B是否为同型矩阵
    if (!check(A) || !check(B) || A->row != B->row || A->col != B->col)return NULL;

    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i * A->col + j] + B->data[i * B->col + j];
            double temp2 = result->data[i * A->col + j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i * A->col + j] += A->data[i * A->col + j] + B->data[i * B->col + j];
            //需要防止float溢出
        }
    }
    return result;
}

matrix *subMatrix(const matrix *A, const matrix *B) {
    // 判断矩阵A和矩阵B是否为同型矩阵

    if (!check(A) || !check(B) || A->row != B->row || A->col != B->col)return NULL;
    // 计算结果并返回
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i * A->col + j] - B->data[i * B->col + j];
            double temp2 = result->data[i * A->col + j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i * A->col + j] += A->data[i * A->col + j] - B->data[i * B->col + j];
        }
    }
    return result;
}

matrix *mulMatrix(const matrix *A, const matrix *B) {
    // 判断矩阵A的列与矩阵B的行是否相等
    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            for (int k = 0; k < A->col; k++) {
                double temp1 = A->data[i * A->col + j] * B->data[i * B->col + j];
                double temp2 = result->data[i * A->col + j] + temp1;
                if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
                result->data[i * result->col + j] += A->data[i * A->col + k] * B->data[k * B->col + j];
            }
        }
    }
    return result;
}


matrix *addScalar(const matrix *A, float b) {
    if (!check(A))return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i * A->col + j] + b;
            double temp2 = result->data[i * A->col + j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i * result->col + j] += A->data[i * A->col + j] + b;
            //需要防止float溢出
        }
    }
    return result;
}

matrix *subScalar(const matrix *A, float b) {
    if (!check(A))return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i * A->col + j] - b;
            double temp2 = result->data[i * A->col + j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i * A->col + j] += A->data[i * A->col + j] - b;
        }
    }
    return result;
}

matrix *mulScalar(const matrix *A, float b) {
    if (!check(A))return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i * A->col + j] * b;
            double temp2 = result->data[i * A->col + j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;
            result->data[i * A->col + j] += A->data[i * A->col + j] * b;
        }
    }
    return result;
}

bool findMax(const matrix *A, float *p) {
    if (!check(A))return false;
    float max = A->data[0];
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->col; ++j) {
            max = max > A->data[i * A->col + j] ? max : A->data[i * A->col + j];
        }
    }
    *p = max;
    return true;
}

bool findMin(const matrix *A, float *p) {
    if (!check(A))return false;
    float min = A->data[0];
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->col; ++j) {
            min = min < A->data[i * A->col + j] ? min : A->data[i * A->col + j];
        }
    }
    *p = min;
    return true;
}

bool copyMatrix(matrix *A, const matrix *B) {
    if (!check(B) || A == NULL || A->data == NULL)return false;
    A->row = B->row;
    A->col = B->col;

    free(A->data);

    A->data = (float **) malloc(sizeof(float *) * A->row);

    A->data = (float *) malloc(sizeof(float) * A->col * A->row);

    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->data[i * A->col + j] = B->data[i * B->col + j];
        }
    }
    return true;
}

bool printMatrix(matrix *A) {
    // 输出矩阵
    if (!check(A))return false;
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            printf("%f ", A->data[i * A->col + j]);
        }
        printf("\n");
    }
    return true;
}

bool transpose(matrix *A) {
    // 转置
    if (!check(A))return false;
    matrix temp = *A;
    createMatrix(A, A->col, A->row);
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->data[i * A->col + j] = temp.data[j * temp.col + i];
        }
    }
    return true;
}

bool check(const matrix *A) {
    //指针假如是一个坏的指针，那么就是NULL，假如user传的不是NULL，是另外的野指针，那不是我的lib需要考虑的问题，是user的事，我不需要大包大揽
    //这是c语言针对pointer的共识
    {
        if (A->row <= 0 || A->col <= 0 || A->data == NULL || A == NULL)return false;
        else return true;
    }
}

bool rowSwap(matrix *A, size_t  i, size_t  j) {
    if (!check(A) || i <= 0 || j <= 0 || i >= A->row || j >= A->row)return false;
    float temp;
    for (int k = 0; k < A->col; k++) {
        temp = A->data[i * A->col + k];
        A->data[i * A->col + k] = A->data[j * A->col + k];
        A->data[j * A->col + k] = temp;
    }
    return true;
}

bool rowSub(matrix *A, size_t  i, size_t  j, float multiple) {
    if (!check(A) || i < 0 || j < 0 || i >= A->row || j >= A->row)return false;
    for (int k = 0; k < A->col; k++) {
        A->data[i * A->col + k] -= A->data[j * A->col + k] * multiple;
        if (fabsf(A->data[i * A->col + k]) < __FLT_EPSILON__) {
            A->data[i * A->col + k] = 0;
        }
    }
    return true;
}

bool colSwap(matrix *A, size_t  i, size_t  j) {
    if (!check(A) || i < 0 || j < 0 || i >= A->col || j >= A->col)return false;
    float temp;
    for (int k = 0; k < A->col; k++) {
        temp = A->data[k * A->col + i];
        A->data[k * A->col + i] = A->data[k * A->col + j];
        A->data[k * A->col + j] = temp;
    }
    return true;
}

bool colSub(matrix *A, size_t  i, size_t  j, float multiple) {
    if (!check(A) || i <= 0 || j <= 0 || i >= A->col || j >= A->col)return false;
    for (int k = 0; k < A->col; k++) {
        A->data[k * A->col + i] -= A->data[k * A->col + j] * multiple;
        if (fabsf(A->data[k * A->col + i]) < __FLT_EPSILON__) {
            A->data[k * A->col + i] = 0;
        }
    }
    return true;
}

bool rowEchelon(matrix *A) {
    if (!check(A))return false;
    long i = 0;
    long j = 0;
    while (i < A->row && j < A->col) {
        if (A->data[i * A->col + j] != 0) {
            float divisor = A->data[i * A->col + j]; // 置1
            for (long k = j; k < A->col; k++) {
                A->data[i * A->col + k] /= divisor;
                if (fabsf(A->data[i * A->col + k]) < __FLT_EPSILON__) {
                    A->data[i * A->col + k] = 0;
                }
            }
            // 同一列下方元素置0
            for (long k = i + 1; k < A->row; k++) {
                rowSub(A, k, i, A->data[k * A->col + j]);
            }
            i++;
            j++;
        } else {
            // 向下寻找matrix[k][j]不为0的行
            int k;
            for (k = i + 1; k < A->row; k++) {
                if (A->data[k * A->col + j] != 0) {
                    break;
                }
            }
            if (k >= A->row) {
                j++;
                continue;
            }
            rowSwap(A, i, k);
            float divisor = A->data[i * A->col + j];
            for (long t = j; t < A->col; t++) {
                A->data[i * A->col + t] /= divisor;
                if (fabsf(A->data[i * A->col + t]) < __FLT_EPSILON__) {
                    A->data[i * A->col + t] = 0;
                }
            }
            for (long t = i + 1; t < A->row; t++) {
                rowSub(A, t, i, A->data[t * A->col + j]);
            }
            i++;
            j++;
        }
    }
    return true;

}


bool colEchelon(matrix *A) {
    if (!check(A))return false;
    int i = 0;
    int j = 0;
    while (i < A->row && j < A->col) {
        if (A->data[i * A->col + j] != 0) {
            // 首元素置1
            float divisor = A->data[i * A->col + j];
            for (int k = i; k < A->row; k++) {
                A->data[k * A->col + j] /= divisor;
                if (fabsf(A->data[k * A->col + j]) < __FLT_EPSILON__) {
                    A->data[k * A->col + j] = 0;
                }
            }
            // 同一行右方元素置0
            for (int k = j + 1; k < A->col; k++) {
                colSub(A, k, j, A->data[i * A->col + k]);
            }
            i++;
            j++;
        } else {
            int k;
            for (k = j + 1; k < A->col; k++) {
                if (A->data[i * A->col + k] != 0) {
                    break;
                }
            }
            if (k >= A->col) {
                i++;
                continue;
            }
            colSwap(A, j, k);
            float divisor = A->data[i * A->col + j];
            for (int t = i; t < A->row; t++) {
                A->data[t * A->col + j] /= divisor;
                if (fabsf(A->data[t * A->col + j]) < __FLT_EPSILON__) {
                    A->data[t * A->col + j] = 0;
                }
            }
            for (int t = j + 1; t < A->col; t++) {
                colSub(A, t, j, A->data[i * A->col + t]);
            }
            i++;
            j++;
        }
        if (A->data[i * A->col + j] == 0) {
            // 向右寻找data[i][k]不为0的列
            int k;
            for (k = j + 1; k < A->col; k++) {
                if (A->data[i * A->col + k] != 0) {
                    break;
                }
            }
            if (k >= A->col) {
                i++;
                continue;
            }
            colSwap(A, j, k);
        }
        // 情形1和2共有的步骤
        float divisor = A->data[i * A->col + j];
        for (int k = i; k < A->row; k++) {
            A->data[k * A->col + j] /= divisor;//置1
            if (fabsf(A->data[k * A->col + j]) < __FLT_EPSILON__) {
                A->data[k * A->col + j] = 0;
            }
        }
        // 同一行右方元素置0
        for (int k = j + 1; k < A->col; k++) {
            colSub(A, k, j, A->data[i * A->col + k]);
        }
        i++;
        j++;
    }
    return true;
}

long rank(matrix *A) {
    rowEchelon(A);// 计算矩阵的秩等价于寻找对角线有多少不为0的元素
    long i = 0;
    long j = 0;
    while (i < A->row && j < A->col) {
        if (A->data[i * A->col + j] == 0) {
            break;
        }
        i++;
        j++;
    }
    return i;
}

bool setElement(matrix *A, size_t  i, size_t  j, float value) {
    if (!check(A) || i >= A->row || j >= A->col)return false;
    else {
        A->data[i * A->col + j] = value;
        return true;
    }
}





