
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <stdbool.h>
#include <string.h>
#include "matrix.h"

Node head;

bool createMatrix(matrix *A, long row,long col){

    Node *newNode;
    newNode = (Node *) malloc(sizeof(Node));
    newNode->next = head.next;
    head.next = newNode;
    newNode->val = (size_t) A;

    if (row <= 0 || col <= 0 || A == NULL)return false;
    A->row = row;
    A->col = col;
    A->data = (float **) malloc(sizeof(float *) * A->row);
    for (int i = 0; i < A->row; i++) {
        A->data[i] = (float *) malloc(sizeof(float) * A->col);
    }
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->data[i][j] = 0;
        }
    }
    return true;
}

bool createMatrixN(matrix *A,long row,long col,float q) {//给用户设置一个所有元素均为same certain number的构造
    Node *newNode;
    newNode = (Node *) malloc(sizeof(Node));
    newNode->next = head.next;
    head.next = newNode;
    newNode->val = (size_t) A;

    if (row <= 0 || col <= 0 || A == NULL)return false;
    A->row = row;
    A->col = col;
    A->data = (float **) malloc(sizeof(float *) * A->row);
    for (int i = 0; i < A->row; i++) {
        A->data[i] = (float *) malloc(sizeof(float) * A->col);
    }
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->data[i][j] = q;
        }
    }return true;
}

bool deleteMatrix(matrix *A) {
    Node *p;
    p = &head;
    while (p->next != NULL) {
        Node *temp = p->next;
        if (temp->val == A) {
            p->next = temp->next;//跳过temp这个节点，从链表中删除
            free(temp);
            if (A == NULL)return false;
            for (int i = 0; i < A->row; i++) {
                free(A->data[i]);
                A->data[i] = NULL;
            }
            free(A->data);
            A->data = NULL;
            return true;
        } else p = p->next;
    }
    return false;//假如到最后还没有遍历到这个指针 说明有误

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
            A->data[i][j] = ans;
            free(str);
            free(endptr);
        }
    }
    return true;
}//只允许按我的要求一位一位读入，不允许直接传一个float**指针来当作data，以保证lib的安全性

matrix *addMatrix(const matrix *A, const matrix *B){
    // 判断矩阵A和矩阵B是否为同型矩阵
    if (!check(A) || !check(B) || A->row != B->row || A->col != B->col)return NULL;

    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i][j] + B->data[i][j];
            double temp2 = result->data[i][j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i][j] += A->data[i][j] + B->data[i][j];
            //需要防止float溢出
        }
    }
    return result;
}

matrix *subMatrix(const matrix *A, const matrix *B){
    // 判断矩阵A和矩阵B是否为同型矩阵

    if (!check(A) || !check(B) || A->row != B->row || A->col != B->col)return NULL;
    // 计算结果并返回
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i][j] - B->data[i][j];
            double temp2 = result->data[i][j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i][j] += A->data[i][j] - B->data[i][j];
        }
    }
    return result;
}

matrix *mulMatrix(const matrix *A, const matrix *B){
    // 判断矩阵A的行与矩阵B的列是否相等
    if (!check(A) || !check(B) || A->col != B->row)return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, B->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            for (int k = 0; k < A->col; k++) {
                double temp1 = A->data[i][j] * B->data[i][j];
                double temp2 = result->data[i][j] + temp1;
                if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
                result->data[i][j] += A->data[i][k] * B->data[k][j];
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
            double temp1 = A->data[i][j] + b;
            double temp2 = result->data[i][j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i][j] += A->data[i][j] + b;
            //需要防止float溢出
        }
    }
    return result;
}

matrix *subScalar(const matrix *A, float b){
    if (!check(A))return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i][j] - b;
            double temp2 = result->data[i][j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;//防止爆float的范围，如果爆了，就返回NULL！
            result->data[i][j] += A->data[i][j] - b;
        }
    }
    return result;
}

matrix *mulScalar(const matrix *A, float b){
    if (!check(A))return NULL;
    matrix *result = (matrix *) malloc(sizeof(matrix));
    createMatrix(result, A->row, A->col);
    for (int i = 0; i < result->row; i++) {
        for (int j = 0; j < result->col; j++) {
            double temp1 = A->data[i][j] * b;
            double temp2 = result->data[i][j] + temp1;
            if (fabsf(temp1) >= __FLT_MAX__ || fabsf(temp2) >= __FLT_MAX__)return NULL;
            result->data[i][j] += A->data[i][j] * b;
        }
    }
    return result;
}

bool findMax(const matrix *A, float *p) {
    if (!check(A))return false;
    float max = A->data[0][0];
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->col; ++j) {
            max = max > A->data[i][j] ? max : A->data[i][j];
        }
    }
    *p = max;
    return true;
}

bool findMin(const matrix *A,float *p){
    if (!check(A))return false;
    float min = A->data[0][0];
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->col; ++j) {
            min = min < A->data[i][j] ? min : A->data[i][j];
        }
    }
    *p = min;
    return true;
}

bool copyMatrix(matrix *A,const matrix *B) {
    if (!check(B) || A == NULL || A->data == NULL)return false;
    A->row = B->row;
    A->col = B->col;

    for (int i = 0; i < A->row; i++) {
        free(A->data[i]);
        A->data[i] = NULL;
    }
    free(A->data);//释放光原来的；

    A->data = (float **) malloc(sizeof(float *) * A->row);
    for (int i = 0; i < A->row; i++) {
        A->data[i] = (float *) malloc(sizeof(float) * A->col);
    }
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->data[i][j] = B->data[i][j];
        }
    }
    return true;
}

bool printMatrix(matrix *A) {
    // 输出矩阵
    if (!check(A))return false;
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            printf("%f ", A->data[i][j]);
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
            A->data[i][j] = temp.data[j][i];
        }
    }
    return true;
}

bool check(const matrix *A) {
    bool flag = false;
    Node *p;
    p = &head;
    while (p->next != NULL) {
        Node *temp = p->next;
        if (temp->val == A) {
            flag = true;
            break;
        } else p = p->next;
    }
    if (flag) {
        if (A->row <= 0 || A->col <= 0 || A->data == NULL || A == NULL)return false;
        else return true;
    }
    return true;
}

bool rowSwap(matrix *A, int i, int j) {
    if (!check(A) || i <= 0 || j <= 0 || i >= A->row || j >= A->row)return false;
    float temp;
    for (int k = 0; k < A->col; k++) {
        temp = A->data[i][k];
        A->data[i][k] = A->data[j][k];
        A->data[j][k] = temp;
    }
    return true;
}

bool rowSub(matrix *A, int i, int j, float multiple) {
    if (!check(A) || i < 0 || j < 0 || i >= A->row || j >= A->row)return false;
    for (int k = 0; k < A->col; k++) {
        A->data[i][k] -= A->data[j][k] * multiple;
        if (fabsf(A->data[i][k]) < __FLT_EPSILON__) {
            A->data[i][k] = 0;
        }
    }
    return true;
}

bool colSwap(matrix *A, int i, int j) {
    if (!check(A) || i < 0 || j < 0 || i >= A->col || j >= A->col)return false;
    float temp;
    for (int k = 0; k < A->col; k++) {
        temp = A->data[k][i];
        A->data[k][i] = A->data[k][j];
        A->data[k][j] = temp;
    }
    return true;
}

bool colSub(matrix *A, int i, int j, float multiple) {
    if (!check(A) || i <= 0 || j <= 0 || i >= A->col || j >= A->col)return false;
    for (int k = 0; k < A->col; k++) {
        A->data[k][i] -= A->data[k][j] * multiple;
        if (fabsf(A->data[k][i]) < __FLT_EPSILON__) {
            A->data[k][i] = 0;
        }
    }
    return true;
}

bool rowEchelon(matrix *A) {
    if (!check(A))return false;
    long i = 0;
    long j = 0;
    while (i < A->row && j < A->col) {
        if (A->data[i][j] != 0) {
            float divisor = A->data[i][j]; // 置1
            for (long k = j; k < A->col; k++) {
                A->data[i][k] /= divisor;
                if (fabsf(A->data[i][k]) < __FLT_EPSILON__) {
                    A->data[i][k] = 0;
                }
            }
            // 同一列下方元素置0
            for (long k = i + 1; k < A->row; k++) {
                rowSub(A, k, i, A->data[k][j]);
            }
            i++;
            j++;
        } else {
            // 向下寻找matrix[k][j]不为0的行
            int k;
            for (k = i + 1; k < A->row; k++) {
                if (A->data[k][j] != 0) {
                    break;
                }
            }
            if (k >= A->row) {
                j++;
                continue;
            }
            rowSwap(A, i, k);
            float divisor = A->data[i][j];
            for (long t = j; t < A->col; t++) {
                A->data[i][t] /= divisor;
                if (fabsf(A->data[i][t]) < __FLT_EPSILON__) {
                    A->data[i][t] = 0;
                }
            }
            for (long t = i + 1; t < A->row; t++) {
                rowSub(A, t, i, A->data[t][j]);
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
        if (A->data[i][j] != 0) {
            // 首元素置1
            float divisor = A->data[i][j];
            for (int k = i; k < A->row; k++) {
                A->data[k][j] /= divisor;
                if (fabsf(A->data[k][j]) < __FLT_EPSILON__) {
                    A->data[k][j] = 0;
                }
            }
            // 同一行右方元素置0
            for (int k = j + 1; k < A->col; k++) {
                colSub(A, k, j, A->data[i][k]);
            }
            i++;
            j++;
        } else {
            int k;
            for (k = j + 1; k < A->col; k++) {
                if (A->data[i][k] != 0) {
                    break;
                }
            }
            if (k >= A->col) {
                i++;
                continue;
            }
            colSwap(A, j, k);
            float divisor = A->data[i][j];
            for (int t = i; t < A->row; t++) {
                A->data[t][j] /= divisor;
                if (fabsf(A->data[t][j]) < __FLT_EPSILON__) {
                    A->data[t][j] = 0;
                }
            }
            for (int t = j + 1; t < A->col; t++) {
                colSub(A, t, j, A->data[i][t]);
            }
            i++;
            j++;
        }
        if (A->data[i][j] == 0) {
            // 向右寻找data[i][k]不为0的列
            int k;
            for (k = j + 1; k < A->col; k++) {
                if (A->data[i][k] != 0) {
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
        float divisor = A->data[i][j];
        for (int k = i; k < A->row; k++) {
            A->data[k][j] /= divisor;//置1
            if (fabsf(A->data[k][j]) < __FLT_EPSILON__) {
                A->data[k][j] = 0;
            }
        }
        // 同一行右方元素置0
        for (int k = j + 1; k < A->col; k++) {
            colSub(A, k, j, A->data[i][k]);
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
        if (A->data[i][j] == 0) {
            break;
        }
        i++;
        j++;
    }
    return i;
}

bool setElement(matrix *A, int i, int j, float value) {
    if (!check(A) || i >= A->row || j >= A->col)return false;
    else {
        A->data[i][j] = value;
        return true;
    }
}





