#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>

#define _N 1000

using namespace std;

int main(void)
{
    char s[_N];
    s = "@\n";
    strcat(s, "hello,world!\n");

    printf("%s", s);

    return 0;
}