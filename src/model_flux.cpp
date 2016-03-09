#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"
#include "model_flux.h"

using namespace std;
ILOSTLBEGIN

void modelFlux(int l, int N[], vector<vector<Arc>> O, int n, int o){

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
                }
            }
        }


        // constraint to u0 be always in the end
        for(int k = 0; k < l-1; k++){
            model.add(z[k][0] <= z[k+1][0]);
        }

        //Solving the problem
        env.setOut(env.getNullStream());
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setError(env.getNullStream());

        cplex.extract(model);



        //timeout
        cplex.setParam(IloCplex::Param::TimeLimit,7200);

        cplex.exportModel("/home/alexsandro/modelFlux.lp");
        if (cplex.solve()) {
            /*for (int k = 0; k < l; k++){
                for(int i = 0; i < n; i++){
                    cout << " k " << k << " i " << i << endl;
                    for(int a = 0; a < n; a++){

                        for(int b = 0; b < n; b++){
                            if(cplex.getValue(x[k][i][a][b]) == 1){
                                cout << "a " << a << " b " << b << " ";
                            }
                        }

                    }
                    cout << endl;

                }

            }*/
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
