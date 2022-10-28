#include <stdio.h>

#include "matrix.c"

int main() {

    matrix a;
    createMatrix(&a, 1, 2);
    setElement(&a, 0, 0, 1);
    setElement(&a, 0, 1, 1);
    printf("matrix a after setting is :\n");
    printMatrix(&a);
    matrix b;
    createMatrix(&b, 2, 2);
    setElement(&b, 0, 0, 0);
    setElement(&b, 0, 1, 1);
    setElement(&b, 1, 0, 1);
    setElement(&b, 1, 1, 1);
    printf("matrix b after setting is :\n");
    printMatrix(&b);
    matrix ans;
    createMatrix(&ans, 2, 2);
    setElement(&ans, 0, 0, 1);
    setElement(&ans, 0, 1, 0);
    setElement(&ans, 1, 0, 0);
    setElement(&ans, 1, 1, 1);
    printf("matrix ans is :\n");
    printMatrix(&ans);

    int c;
    scanf("%d", &c);
    while (c) {
        if (c & 1) ans = *mulMatrix(&ans, &b);
        c >>= 1;
        b = *mulMatrix(&b, &b);
    }
    ans = *mulMatrix(&a, &ans);
    printMatrix(&ans);

    matrix test1;
    createMatrix(&test1, 3, 3);
    initMatrix(&test1);
    printMatrix(&test1);

    matrix test2;
    matrix *p2 = &test2;
    createMatrix(p2, 4, 4);
    setElement(p2, 0, 0, 4);
    setElement(p2, 0, 1, 3);
    setElement(p2, 0, 2, 2);
    setElement(p2, 0, 3, 1);
    setElement(p2, 1, 0, 1);
    setElement(p2, 1, 1, 1);
    setElement(p2, 1, 2, 1);
    setElement(p2, 1, 3, 1);
    setElement(p2, 2, 0, 0.5f);
    setElement(p2, 2, 1, 1.5f);
    setElement(p2, 2, 2, 0);
    setElement(p2, 2, 3, 0);
    setElement(p2, 3, 0, 0);
    setElement(p2, 3, 1, 0);
    setElement(p2, 3, 2, 0);
    setElement(p2, 3, 3, 1);

    printf("p2 as follows: \n");
    printMatrix(p2);

    rowEchelon(p2);
    printf("matrix after rowEchelon is :\n");
    printMatrix(p2);

    rowSub(p2, 1, 2, 0.1f);
    printf("matrix after rowSub is :\n");
    printMatrix(p2);

    rowSwap(p2, 2, 3);
    printf("matrix after rowSwap is :\n");
    printMatrix(p2);

    colSub(p2, 1, 2, 0.1f);
    printf("matrix after colSub is :\n");
    printMatrix(p2);

    colSwap(p2, 2, 3);
    printf("matrix after colSwap is :\n");
    printMatrix(p2);



    printf("rank of matrix is :%ld\n", rank(p2));
    float test3;
    float *p3 = &test3;
    findMax(p2, p3);
    printf("max of matrix is :%f\n", *p3);
    findMin(p2, p3);
    printf("min of matrix is :%f\n", *p3);

    bool flag = colSwap(p2, -1, -1);
    printf("Input invalid data into method cloSwap and return :%d\n",flag);

    bool flag2= deleteMatrix(p3);//假如传入一个野针 看会不会被检测出来
    if (!flag2) printf("error in delete matrix for input a wrong pointer\n");

    matrix test5;
    matrix test6;
    createMatrixN(&test5,10,5,2.0f);
    createMatrixN(&test6,9,5,3.0f);
    bool flag3= addMatrix(&test5,&test6);
    bool flag4= subMatrix(&test5,&test6);
    bool flag5= mulMatrix(&test5,&test6);
    printf("return value of addMatrix is: %d\n",flag3);
    printf("return value of subMatrix is: %d\n",flag4);
    printf("return value of mulMatrix is: %d\n",flag5);

    matrix test4;
    matrix *p4 = &test4;
    createMatrix(p4, 4, 4);
    setElement(p4, 0, 0, 4);
    setElement(p4, 0, 1, 3.2f);
    setElement(p4, 0, 2, 2);
    setElement(p4, 0, 3, 1);
    setElement(p4, 1, 0, 1);
    setElement(p4, 1, 1, 1.1f);
    setElement(p4, 1, 2, 1);
    setElement(p4, 1, 3, 1);
    setElement(p4, 2, 0, 0.5f);
    setElement(p4, 2, 1, 1.5f);
    setElement(p4, 2, 2, 3.3f);
    setElement(p4, 2, 3, 0);
    setElement(p4, 3, 0, 0);
    setElement(p4, 3, 1, 0);
    setElement(p4, 3, 2, 0);
    setElement(p4, 3, 3, 1);

    printf("p2 before copy is:\n");
    printMatrix(p2);
    copyMatrix(p2,p4);
    printf("p2 after copy is:\n");
    printMatrix(p2);

    transpose(p2);
    printMatrix(p2);

    createMatrixN(p2,4,4,2);
    printMatrix(p2);




    return 0;
}
