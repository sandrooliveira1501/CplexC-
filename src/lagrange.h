#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"

using namespace std;

class Lagrange {
public:

    float alpha[100];
    float beta[100][100][100];
    int ZUB;
    int LB;
    int ZLB;
    int l;
    int *N;
    vector<vector<Arc>> O;
    int n, o;

    IloArray<IloBoolVarArray> z;
    IloArray<IloArray<IloArray<IloIntVarArray>>> x;

    Lagrange(int l, int *N,
          vector<vector<Arc>> O, int n, int o, int LB);

    IloModel prepareModel(IloEnv env);

    void execute();
};

