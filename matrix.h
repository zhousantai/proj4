//
// Created by winnie泰 on 2022/10/24.
//
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <stdbool.h>
#include <string.h>


#ifndef PROJ3TEST_TEST_H
#define PROJ3TEST_TEST_H


typedef struct {
    long row;            // 行
    long col;            // 列
    float *data;     // 元素
} matrix;

typedef struct {
    size_t val;
    struct Node *next;
} Node;

bool createMatrix(matrix *A, long row, long col);// 创建矩阵，不允许用户自己传入一个指针来当数据，必须按照我的要求一个一个输入
bool createMatrixN(matrix *A, long row, long col, float q);

bool deleteMatrix(matrix *A);                                           // 释放
bool initMatrix(matrix *A);                                    // 初始化
bool copyMatrix(matrix *A, const matrix *B);

bool findMax(const matrix *A, float *p);

bool findMin(const matrix *A, float *p);

bool transpose(matrix *A);

bool printMatrix(matrix *A);

bool check(const matrix *A);

bool rowSwap(matrix *A, int i, int j);                   // 交换两行
bool rowSub(matrix *A, int i, int j, float multiple);    // 两行相减
bool colSwap(matrix *A, int i, int j);                   // 交换两列
bool colSub(matrix *A, int i, int j, float multiple);    // 两列相减
bool rowEchelon(matrix *A);                                    // 将矩阵转换为行阶梯型
bool colEchelon(matrix *A);                                    // 将矩阵转换为列阶梯型
long rank(matrix *A);                                             // 求矩阵的秩
bool setElement(matrix *A, int i, int j, float value);//给用户修改矩阵单个数据的接口；

matrix *addMatrix(const matrix *A, const matrix *B);

matrix *subMatrix(const matrix *A, const matrix *B);

matrix *mulMatrix(const matrix *A, const matrix *B);

matrix *addScalar(const matrix *A, float b);

matrix *subScalar(const matrix *A, float b);

matrix *mulScalar(const matrix *A, float b);

matrix *mulMatirx_plain(const matrix *A, const matrix *B);

matrix *mulMatirx_improved(const matrix *A, const matrix *B);

void devision(const matrix*A,int t,matrix* A1,matrix *A2,matrix*A3,matrix*A4);
void merge(matrix*ans,matrix* ans1,matrix* ans2,matrix* ans3,matrix* ans4);

#endif //PROJ3TEST_TEST_H
