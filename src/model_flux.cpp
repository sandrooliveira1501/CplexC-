#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"

#include "model_flux.h"

using namespace std;
ILOSTLBEGIN

void modelFlux(int l, int N[],int  ord[], vector<vector<Arc>> O, int n, int o){

    IloEnv env;

    try{

        IloModel model(env);

        // var z_{k,mi}
        IloArray<IloBoolVarArray> z(env, l);

        for(int i = 0; i < l; i++){

            z[i] = IloBoolVarArray (env, o);

        }

        // var x_{kiab}
        IloArray<IloArray<IloArray<IloIntVarArray>>> x(env, l);
        for (int i = 0; i < l; i++) {
            x[i] = IloArray<IloArray<IloIntVarArray> > (env, n);
            for (int j = 0; j < n; j++) {
                x[i][j] = IloArray<IloIntVarArray> (env, n);
                for (int k = 0; k < n; k++) {
                    x[i][j][k] = IloIntVarArray (env,n,0,IloIntMax);
                }
            }
        }

        //objective function

        IloExpr obj(env);

        for(int k = 0; k < l;k++){
            for(int i = 1; i < o; i++){
                obj += z[k][i];
            }
        }

        model.add(IloMinimize(env, obj));

        //one operator by level
        for(int k = 0; k < l; k++){
            model.add(IloSum(z[k]) == 1);
        }

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
            expr.end();
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
                    expr.end();
                }

            }
        }

        // just the arcs that are in mi can be used


            for(int a = 0; a < n; a++){
                for(int b = 0; b < n; b++){
                    for(int k  = 0; k < l; k++){
                        IloExpr flowCapacityTerm1(env);
                        IloExpr flowCapacityTerm2(env);

                        for(int i = 0; i < n; i++){
                            flowCapacityTerm1 += x[k][i][a][b];
                        }

                        for(int y = 0; y < o; y++){
                            Arc arc(a,b);
                            if(contains(arc, O[y])){
                                flowCapacityTerm2 += z[k][y];
                            }
                        }

                        model.add(flowCapacityTerm1 <= flowCapacityTerm2);

                        flowCapacityTerm1.end();
                        flowCapacityTerm2.end();
                }
            }
        }


        // constraint to u0 be always in the end
        for(int k = 0; k < l-1; k++){
            model.add(z[k][0] <= z[k+1][0]);
        }

        //extra

        /*for (int k = 0; k < l-1; k++){

            //cout << "k = " << k << endl;
            for (int i = 0; i < (n-1); i++){

                for(int j = 1; (j+k) < l; j++){

                    //cout << "i = " << i << endl;
                    for(int a = 0; a < n; a++){

                        //cout << "a = " << a << endl;
                        for(int a2 = 0; a2 < n; a2++){

                            //cout << "a2 = " << a2 << endl;
                            for(int b = 0; b < (n-1); b++){

                                //cout << "b = " << b << endl;
                                for(int b2 = 0; b2 < (n-1); b2++){

                                    //cout << "b2 = " << b2 << endl;
                                    //cout << ord[i] << endl;
                                    //cout << ord[i+1] << endl;
                                    model.add(x[k][ord[i]][a][b] + x[k][ord[i+1]][a2][b+1] + x[k+j][ord[i]][b][b2] <= x[k+j][ord[i+1]][b+1][b2+1] + 2);

                                }
                            }
                        }
                    }
                }
            }
        }*/

        for(int k = 0; k < (l-1); k++){
            for(int a = 0; a < n; a++){
                model.add(x[k][ord[0]][a][1] <= x[k+1][ord[0]][1][1]);
            }
        }

        for(int k = 0; k < (l-1); k++){
            for(int a = 0; a < n; a++){
                model.add(x[k][ord[n-1]][a][n-1] <= x[k+1][ord[n-1]][n-1][n-1]);
            }
        }

        //Solving the problem
        //env.setOut(env.getNullStream());
        IloCplex cplex(model);
        //cplex.setOut(env.getNullStream());
        //cplex.setWarning(env.getNullStream());
        //cplex.setError(env.getNullStream());

        cplex.extract(model);



        //timeout
        cplex.setParam(IloCplex::Param::TimeLimit,7200);

        //cplex.exportModel("home/alexsandro/model1.lp");
        if (cplex.solve()) {
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
