#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define _M 100000
#define _N 1000
#define _T 100
// #define DEBUG

struct formula
{
    double coeff[_N];
    char term[_N][_T];
    struct formula *next;
}*Mem[_M] = {0};
int MemFlag[_M] = {0};

struct formula *createFml();
struct formula *generateFml(int N, int SS, int MMs);
struct formula *addFml(struct formula *f1, struct formula *f2);
struct formula *coupleFml(struct formula *f1, int dSS, int dMMs, int SS, int MMs);
struct formula *connectFml(struct formula *f1, struct formula *f2);
void printFml(struct formula *fml);
void delFml(struct formula *fml);

int main(void)
{
    struct formula *fml;
    struct formula *first;
    int N, SS, MMs;
    double S, Ms;

    while (1)
    {
        printf("Enter: N S Ms (Quit if N = 0)\n");
        scanf("%d", &N);
        if (N <= 0)
            break;
        scanf("%lf %lf", &S, &Ms);
        SS = (int)2.0 * (S + 0.1);
        MMs = (int)2.0 * (Ms + 0.1);
        fml = generateFml(N, SS, MMs);
        printFml(fml);
        delFml(fml);
    }

    return 0;
}

struct formula *createFml()
{
    struct formula *fml;
    fml = (struct formula *)malloc(sizeof(struct formula));
    memset(fml, 0, sizeof(struct formula));
    fml->next = NULL;

    return fml;
}

struct formula *generateFml(int N, int SS, int MMs)
{
    struct formula *fml, *f1, *f2;
    struct formula *fb1, *fb2, *fa1, *fa2;

    if (SS < 0 || SS > N || MMs > SS || MMs < 0 - SS)
        return NULL;
    if (N % 2 != SS % 2 || SS % 2 != abs(MMs) % 2)
    {
        printf("Invalid Input!\n");
        return NULL;
    }
    if (N == 1)
    {
        struct formula *fml;
        fml = createFml();
        fml->coeff[0] = 1.0;
        fml->term[0][0] = (MMs == 1) ? 'a' : 'b';
        
        return fml;
    }

    f1 = generateFml(N - 1, SS + 1, MMs - 1);
    fb1 = coupleFml(f1, 1, -1, SS, MMs);
    f2 = generateFml(N - 1, SS + 1, MMs + 1);
    fb2 = coupleFml(f2, 1, 1, SS, MMs);
    fa1 = addFml(fb1, fb2);
    // #ifdef DEBUG
    //     printf("(fa1)N=%d, 2S=%d, 2Ms=%d\n", N, SS, MMs);
    //     printFml(fa1);
    //     printf("END\n");
    // #endif
    delFml(fb1);
    delFml(fb2);
    f1 = generateFml(N - 1, SS - 1, MMs - 1);
    fb1 = coupleFml(f1, -1, -1, SS, MMs);
    f2 = generateFml(N - 1, SS - 1, MMs + 1);
    fb2 = coupleFml(f2, -1, 1, SS, MMs);
    fa2 = addFml(fb1, fb2);
    // #ifdef DEBUG
    // printf("(fa2)N=%d, 2S=%d, 2Ms=%d\n", N, SS, MMs);
    //     printFml(fa2);
    //     printf("END\n");
    // #endif
    delFml(fb1);
    delFml(fb2);
    fml = connectFml(fa1, fa2);
    delFml(f1);
    delFml(f2);
#ifdef DEBUG
    printf("(fml)N=%d, 2S=%d, 2Ms=%d\n", N, SS, MMs);
    printFml(fml);
    printf("END\n\n");
#endif

    return fml;
}

struct formula *addFml(struct formula *q1, struct formula *q2)
{
    int i = 0, j = 0;
    int ii = 0, beginning = 1;
    struct formula *fml, *head, *p;

    if (q1 == NULL && q2 == NULL)
        return NULL;

    while (q1 != NULL || q2 != NULL)
    {
        fml = createFml();
        ii = 0;
        if (q1 != NULL)
        {
            for (; q1->term[ii][0] != 0; ii++)
            {
                fml->coeff[ii] = q1->coeff[ii];
                for (int jj = 0; q1->term[ii][jj] != 0; jj++)
                    fml->term[ii][jj] = q1->term[ii][jj];
            }
            q1 = q1->next;
        }
        if (q2 != NULL)
        {
            for (int m = 0; q2->term[m][0] != 0; m++)
            {
                fml->coeff[m + ii] = q2->coeff[m];
                for (int jj = 0; q2->term[m][jj] != 0; jj++)
                    fml->term[m + ii][jj] = q2->term[m][jj];
            }
            q2 = q2->next;
        }
        if (beginning)
        {
            head = fml;
            p = fml;
            beginning = 0;
        }
        else
        {
            p->next = fml;
            p = p->next;
        }
    }

    return head;
}

struct formula *coupleFml(struct formula *f1, int dSS, int dMMs, int SS, int MMs)
{
    double S, Ms, cff = 0;
    int len, beginning = 1;
    char spin;
    struct formula *fml, *head, *p;
    

    if (f1 == NULL)
        return NULL;

    S = SS / 2.0;
    Ms = MMs / 2.0;
    if (dSS == 1)
    {
        if (dMMs == 1)
        {
            cff = sqrt((S + Ms + 1.0) / ((S + S + 2.0)));
            spin = 'b';
        }
        else if (dMMs == -1)
        {
            cff = 0.0 - sqrt((S - Ms + 1.0) / (S + S + 2.0));
            spin = 'a';
        }
    }
    else if (dSS == -1)
    {
        if (dMMs == 1)
        {
            cff = sqrt((S - Ms) / (S + S));
            spin = 'b';
        }
        else if (dMMs == -1)
        {
            cff = sqrt((S + Ms) / (S + S));
            spin = 'a';
        }
    }
    // #ifdef DEBUG
    //     printf("(cff)2S=%d, 2Ms=%d (%d, %d)\n", SS + dSS, MMs + dMMs, dSS, dMMs);
    //     printf("%lf %lf %lf\n", S, Ms, cff);
    //     printf("END\n");
    // #endif
    while (f1 != NULL)
    {
        fml = createFml();        
        for (int i = 0; f1->term[i][0] != 0; i++)
        {
            fml->coeff[i] = cff * f1->coeff[i];
            for (len = 0; f1->term[i][len] != 0; len++)
                fml->term[i][len] = f1->term[i][len];
            fml->term[i][len] = spin;
        }
        f1 = f1->next;
        if (beginning)
        {
            head = fml;
            p = fml;
            beginning = 0;
        }
        else
        {
            p->next = fml;
            p = p->next;
        }
    }

    return head;
}

struct formula *connectFml(struct formula *fa1, struct formula *fa2)
{
    struct formula *head;

    if (fa1 == NULL)
        head = fa2;
    else
    {
        head = fa1;
        while (fa1->next != NULL)
            fa1 = fa1->next;
        fa1->next = fa2;
    }

    return head;
}

void printFml(struct formula *fml)
{
    int len, num = 0, beginning;
    double cff, count;
    char sgn;

    while (fml != NULL)
    {
        beginning = 1;
        count = 0;
        for (len = 0; fml->term[len][0] != 0; len++)
            ;
        for (int i = 0; i < len; i++)
        {
            cff = fml->coeff[i];
            if (cff > 1e-8)
                sgn = '+';
            if (cff < -1e-8)
            {
                sgn = '-';
                cff = 0.0 - cff;
            }
            count += cff * cff;
            if (beginning)
            {
                printf("%c:", 'A' + num);
                num++;
                beginning = 0;
            }
            printf(" %c %.8f|%s>", sgn, cff, fml->term[i]);
        }
        printf("\n    Sigma c^2 = %.2f\n", count);
        fml = fml->next;
    }
    printf("\n");

    return;
}

void delFml(struct formula *q)
{
    struct formula *p;
    while (q != NULL)
    {
        p = q->next;
        free(q);
        q = p;
    }

    return;
}
