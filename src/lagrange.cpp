#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"
#include "lagrange.h"
#include <Python.h>


using namespace std;
ILOSTLBEGIN

void llbp(int l, int N[], vector<vector<Arc>> O, int n, int o){

    IloEnv env;

    float alpha[l];

    float beta[l][n][n];

    for(int i = 0; i < l; i++){
        alpha[i] = 0;
        for(int a = 0; a < n; a++){
            for(int b = 0; b < n; b++){
                beta[i][a][b] = 0;
            }
        }
    }

    try{

        IloModel model(env);

        IloArray<IloBoolVarArray> z(env, l);

        for(int i = 0; i < l; i++){

            z[i] = IloBoolVarArray (env, o);

        }

        IloArray<IloArray<IloArray<IloIntVarArray>>> x(env, l);
        for (int i = 0; i < l; i++) {
            x[i] = IloArray<IloArray<IloIntVarArray> > (env, n);
            for (int j = 0; j < n; j++) {
                x[i][j] = IloArray<IloIntVarArray> (env, n);
                for (int k = 0; k < n; k++) {
                    x[i][j][k] = IloIntVarArray (env,n,0,INFINITY);
                }
            }
        }

        IloExpr obj(env);

        for(int k = 0; k < l;k++){
            for(int u = 0; u < o; u++){

                float c;

                c = 1 - alpha[k];

                for(int e = 0; e < n; e++){
                    int a = O[u][e].a;
                    int b = O[u][e].b;

                    c -= beta[k][a][b];

                }

                obj+= z[k][u]*c;

            }

            obj += -z[k][0];

            obj += alpha[k];

            for(int i = 0; i < n; i++){

                for(int a = 0; a < n; a++){

                    for(int b = 0; b < n; b++){

                        obj += beta[k][a][b]*x[k][i][a][b];

                    }

                }

            }

        }

        model.add(IloMinimize(env, obj));


        // outgoing
        for(int i = 0; i < n; i++){
           model.add(IloSum(x[0][i][i]) == 1);
        }

        // ingoing
        for(int i = 0; i < n; i++){
            IloExpr expr(env);
            for(int j = 0; j < n; j++){
                expr += x[l-1][i][j][N[i]];
            }
            model.add(expr == 1);
        }

        //conservation
        for(int i = 0; i < n; i++){
            for(int k = 1; k < l; k++){
                for(int a = 0; a < n; a++){
                    IloExpr expr(env);

                    for(int j = 0; j < n; j++){
                        expr +=  x[k][i][a][j];
                        expr -=  x[k-1][i][j][a];
                    }

                    model.add(expr == 0);
                }

            }
        }

        //capacite
        for(int k = 0; k < l; k++){
            for(int a = 0; a < n; a++){
                for(int b = 0; b < n; b++){
                    IloExpr expr(env);

                    for(int i = 0; i < n; i++){

                        expr += x[k][i][a][b];

                    }

                    model.add(expr <= 1);

                }
            }
        }

        IloCplex cplex(model);
        cplex.extract(model);
        cplex.exportModel("/home/alexsandro/llp.lp");

        if (cplex.solve()) {
            /*for (int k = 0; k < l; k++){
                for(int i = 0; i < n; i++){
                    cout << " k " << k << " i " << i << endl;
                    for(int a = 0; a < n; a++){

                        for(int b = 0; b < n; b++){
                            cout << "a " << a << " b " << b << " x " << cplex.getValue(x[k][i][a][b]) << " - ";
                        }

                    }
                    cout << endl;

                }

            }*/

            for(int k = 0; k < l; k++){
                for(int u = 0; u < o; u++){

                    //cout << "z_" << k << "_" << u << " = " << cplex.getValue(z[k][u]) << endl;

                }
            }

            cout << "Optimal value: " << cplex.getObjValue() << endl;
            cout << cplex.getTime() << endl;
        }else{
            cout << "timeout" << endl;
        }

        cplex.end();
        model.end();
        obj.end();


    }catch(IloException& e){
        cerr << "Concert exception caught: " << e << endl;
    }catch(...){
        cerr << "Unknow exception caught" << endl;
    }

    env.end();
}
